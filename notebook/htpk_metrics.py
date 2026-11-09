"""Metrics comparing predicted and measured plasma concentration-time (Cp) curves.

Standalone module used by ``data_read.ipynb`` (requires numpy, pandas, matplotlib).
It reproduces the metric calculation of the original HTPK benchmark analysis:

1. **Measured curve of a case.** A single replicate is used as is. Several replicates
   are linearly interpolated onto the union of their time points, restricted to the
   interval covered by *all* replicates, and averaged. Their bounds are combined
   pessimistically (lowest lower bound, highest upper bound).
2. **Predicted curve on the measured grid.** The simulated curve is averaged in
   1-minute bins (labelled by the bin start), truncated at the last measured time and
   linearly interpolated onto the measured time points. Values beyond the last bin
   are held constant. ``resample_minutes=None`` skips the binning (exact interpolation).
3. **Error bounds (optional).** If a case has complete lower/upper bounds, a prediction
   outside ``[lower, upper]`` is compared with the violated bound instead of the
   measurement. With ``zero_error_within_bounds`` a prediction inside the bounds counts
   as exact. If any bound of the case is missing, the plain comparison is used.
4. **Exclusion of ~zero values.** Points where the prediction or the observation is at
   or below ``min_threshold`` are excluded (typically t = 0).
5. **Metrics.** Absolute, relative and log2-ratio errors; the main metric is the mean
   absolute log2 error (MAL2E): 1 = two-fold error on average.

Times are float hours, concentrations ng/mL.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Column names of data/Cptimecurves_machine_readable.csv and notebook/htpk_predictions.csv
CASE_COLUMNS = ["Compound", "Smiles", "Dose (mg)", "Route", "Species"]
TIME = "Time (h)"
CONC = "Concentration (ng/mL)"
REPLICATE = "Replicate"
LOWER = "Lower bound"
UPPER = "Upper bound"

METRIC_NAMES = [
    "n_points", "n_valid_points", "n_excluded_points", "bounds_used",
    "MSE", "RMSE", "MAE", "Mean_Error",
    "Mean_Relative_Error", "Mean_Absolute_Relative_Error", "RMSE_Relative",
    "Mean_Log2_Error", "Mean_Absolute_Log2_Error", "RMSE_Log2",
    "R_squared", "Correlation",
]


@dataclass(frozen=True)
class MetricSettings:
    use_error_bounds: bool = True
    zero_error_within_bounds: bool = True
    min_threshold: float = 1e-10
    resample_minutes: float | None = 1.0  # original analysis: 1-minute bin means


@dataclass
class Observed:
    """Measured curve on its time grid (aggregated over replicates if needed)."""

    time_h: np.ndarray
    conc: np.ndarray
    lower: np.ndarray
    upper: np.ndarray

    @property
    def has_bounds(self) -> bool:
        return bool(np.isfinite(self.lower).all() and np.isfinite(self.upper).all())


# --------------------------------------------------------------------------- data


def case_label(case: pd.Series | dict) -> str:
    return f"{case['Compound']} | {case['Species']} | {case['Route']} | {case['Dose (mg)']:g} mg"


def replicate_frames(case_data: pd.DataFrame) -> list[pd.DataFrame]:
    """Measured rows of one case -> one frame per replicate, sorted by time."""
    return [g.sort_values(TIME) for _, g in case_data.groupby(REPLICATE, sort=True)]


def case_category(case_data: pd.DataFrame) -> str:
    n_rep = case_data[REPLICATE].nunique()
    frac = case_data[LOWER].notna().mean()
    curves = "multiple curves" if n_rep > 1 else "single curve"
    bounds = "full bounds" if frac == 1 else "partial bounds" if frac > 0 else "no bounds"
    return f"{curves}, {bounds}"


# ------------------------------------------------------------------- measured side


def _interp_valid(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Linear interpolation through the finite points of y; NaN if there are none."""
    ok = np.isfinite(y)
    if not ok.any():
        return np.full(len(x_new), np.nan)
    return np.interp(x_new, x[ok], y[ok])


def aggregate_replicates(replicates: list[pd.DataFrame]) -> Observed:
    """Combine replicate curves into one observed curve (see module docstring, step 1)."""
    if not replicates:
        raise ValueError("No measured replicates")
    if len(replicates) == 1:
        r = replicates[0]
        return Observed(r[TIME].to_numpy(float), r[CONC].to_numpy(float),
                        r[LOWER].to_numpy(float), r[UPPER].to_numpy(float))

    t_start = max(r[TIME].min() for r in replicates)
    t_end = min(r[TIME].max() for r in replicates)
    grid = np.unique(np.concatenate([r[TIME].to_numpy(float) for r in replicates]))
    grid = grid[(grid >= t_start) & (grid <= t_end)]

    def stack(col: str) -> np.ndarray:
        return np.vstack([_interp_valid(grid, r[TIME].to_numpy(float), r[col].to_numpy(float)) for r in replicates])

    with warnings.catch_warnings():  # all-NaN bound columns are expected (no bounds)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return Observed(grid, np.nanmean(stack(CONC), axis=0),
                        np.nanmin(stack(LOWER), axis=0), np.nanmax(stack(UPPER), axis=0))


# ------------------------------------------------------------------ predicted side


def bin_mean(time_h: np.ndarray, conc: np.ndarray, minutes: float) -> tuple[np.ndarray, np.ndarray]:
    """Mean per left-labelled time bin, as pandas ``resample(f"{minutes}min").mean()``."""
    bins = np.floor(np.asarray(time_h) * 60.0 / minutes + 1e-9).astype(np.int64)
    uniq, inv = np.unique(bins, return_inverse=True)
    means = np.bincount(inv, weights=conc) / np.bincount(inv)
    return uniq * minutes / 60.0, means


def predicted_on_grid(sim_time_h: np.ndarray, sim_conc: np.ndarray, grid_h: np.ndarray, t_last_h: float,
                      resample_minutes: float | None) -> np.ndarray:
    """Simulated curve evaluated at the measured times (module docstring, step 2)."""
    if len(sim_time_h) == 0:
        raise ValueError("No simulated curve")
    order = np.argsort(sim_time_h)
    t, c = np.asarray(sim_time_h, float)[order], np.asarray(sim_conc, float)[order]
    grid_h = np.asarray(grid_h, float)
    if resample_minutes:
        # Regular bin grid: empty bins are linearly interpolated between non-empty ones,
        # the grid is truncated at the last bin label <= last measured time, and values
        # beyond that label are held constant.
        t, c = bin_mean(t, c, resample_minutes)
        t_cut = np.floor(t_last_h * 60.0 / resample_minutes + 1e-9) * resample_minutes / 60.0
        grid_h = np.minimum(grid_h, t_cut)
    return np.interp(grid_h, t, c)


# -------------------------------------------------------------------------- metrics


def effective_observed(pred: np.ndarray, obs: np.ndarray, lower: np.ndarray, upper: np.ndarray) -> np.ndarray:
    """Measured value, replaced by the violated bound when pred lies outside it."""
    eff = obs.copy()
    above, below = pred > upper, pred < lower
    eff[above] = upper[above]
    eff[below] = lower[below]
    return eff


def timeseries_metrics(pred: np.ndarray, observed: Observed, settings: MetricSettings) -> dict[str, float]:
    """Metrics of predicted values (on observed.time_h) against the observation."""
    pred = np.asarray(pred, dtype=float).copy()
    obs = observed.conc.astype(float).copy()

    bounds_used = False
    if settings.use_error_bounds and observed.has_bounds and not np.any(observed.lower > observed.upper):
        bounds_used = True
        lo, up = observed.lower, observed.upper
        obs_eff = effective_observed(pred, obs, lo, up)
        if settings.zero_error_within_bounds:
            inside = (pred >= lo) & (pred <= up)
            pred[inside] = obs[inside]
        obs = obs_eff

    finite = np.isfinite(pred) & np.isfinite(obs)
    pred, obs = pred[finite], obs[finite]
    n_points = len(pred)
    if n_points == 0:
        raise ValueError("No valid data points")

    thr = settings.min_threshold
    keep = (np.abs(pred) > thr) & (np.abs(obs) > thr)
    p, o = pred[keep], obs[keep]
    if len(p) == 0:
        raise ValueError("No valid data points after excluding zeros")

    err = p - o
    rel = err / o
    log2 = np.log2(np.maximum(np.abs(p), thr) / np.maximum(np.abs(o), thr))
    corr = np.nan
    if len(p) > 1 and np.std(p) > 0 and np.std(o) > 0:
        corr = float(np.corrcoef(p, o)[0, 1])

    return {
        "n_points": n_points,
        "n_valid_points": int(len(p)),
        "n_excluded_points": int(n_points - len(p)),
        "bounds_used": bounds_used,
        "MSE": float(np.mean(err**2)),
        "RMSE": float(np.sqrt(np.mean(err**2))),
        "MAE": float(np.mean(np.abs(err))),
        "Mean_Error": float(np.mean(err)),
        "Mean_Relative_Error": float(np.mean(rel)),
        "Mean_Absolute_Relative_Error": float(np.mean(np.abs(rel))),
        "RMSE_Relative": float(np.sqrt(np.mean(rel**2))),
        "Mean_Log2_Error": float(np.mean(log2)),
        "Mean_Absolute_Log2_Error": float(np.mean(np.abs(log2))),
        "RMSE_Log2": float(np.sqrt(np.mean(log2**2))),
        "R_squared": corr**2 if np.isfinite(corr) else np.nan,
        "Correlation": corr,
    }


def case_metrics(replicates: list[pd.DataFrame], prediction: pd.DataFrame, settings: MetricSettings) -> dict[str, float]:
    """All metrics of one case: measured replicates vs predicted curve (columns Time (h), Concentration (ng/mL))."""
    observed = aggregate_replicates(replicates)
    t_last = max(r[TIME].max() for r in replicates)
    pred = predicted_on_grid(prediction[TIME].to_numpy(), prediction[CONC].to_numpy(),
                             observed.time_h, t_last, settings.resample_minutes)
    return timeseries_metrics(pred, observed, settings)


# ------------------------------------------------------------------------- plotting

PREDICTED_COLOR = "#2a78d6"
MEASURED_COLOR = "#eb6834"
REPLICATE_MARKERS = ["o", "s", "^", "D", "v", "P"]


def plot_case(ax: plt.Axes, replicates: list[pd.DataFrame], prediction: pd.DataFrame,
              title: str = "", log_y: bool = False, predicted_label: str = "Predicted (HTPK)") -> None:
    """Measured replicates (+ bounds band) and the predicted curve up to the last measured time."""
    t_end = max(r[TIME].max() for r in replicates)
    for i, rep in enumerate(replicates):
        name = f"Measured (replicate {rep[REPLICATE].iloc[0]})" if len(replicates) > 1 else "Measured"
        if rep[LOWER].notna().any():  # NaN gaps are left unfilled
            ax.fill_between(rep[TIME], rep[LOWER], rep[UPPER], color=MEASURED_COLOR, alpha=0.12, linewidth=0)
        ax.plot(rep[TIME], rep[CONC], color=MEASURED_COLOR, linewidth=1, alpha=0.8,
                marker=REPLICATE_MARKERS[i % len(REPLICATE_MARKERS)], markersize=5,
                markeredgecolor="white", markeredgewidth=0.8, label=name)
    p = prediction[prediction[TIME] <= t_end]
    ax.plot(p[TIME], p[CONC], color=PREDICTED_COLOR, linewidth=2, label=predicted_label)
    if log_y:
        ax.set_yscale("log")
    if any(r[LOWER].notna().any() for r in replicates):
        ax.fill_between([], [], color=MEASURED_COLOR, alpha=0.12, label="Measured bounds")
    ax.set_xlabel(TIME)
    ax.set_ylabel(CONC)
    ax.set_title(title, loc="left", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize="small")

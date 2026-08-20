# pk-datasets

A curated, machine-readable benchmark database of plasma concentration–time
profiles for evaluating high-throughput mechanistic pharmacokinetic (PK)
prediction workflows.

This repository is the versioned, continually maintained release of the dataset
described in:

> Jones J.O., Bachorz R.A., Lawless M.S., Miller D.W., et al.
> *A curated benchmark database for high-throughput mechanistic pharmacokinetic prediction.*
> Research Square preprint, DOI: [10.21203/rs.3.rs-10671182/v1](https://doi.org/10.21203/rs.3.rs-10671182/v1)

## What's in the dataset

| | |
|---|---|
| Cases (compound × species × route × dose) | 180 |
| Unique compounds | 76 |
| Species | human, cynomolgus monkey (*Macaca fascicularis*), rhesus monkey (*Macaca mulatta*), dog, rat, mouse |
| Routes | IV bolus, immediate-release tablet (IRT), immediate-release solution (IRS) |

Each case carries a digitized mean concentration–time curve, upper/lower
uncertainty bounds where the source reported them (SEM or 95% CI), a
standardized chemical structure, and dose metadata. Fourteen cases are averaged
across multiple independent studies.

## Repository structure

```
.
├── README.md
└── data
    └── Cptimecurves_machine_readable.csv
```

## Data dictionary

Columns of `data/Cptimecurves_machine_readable.csv`:

| Column | Description |
|---|---|
| `compound` | Compound name |
| `SMILES` | Standardized chemical structure |
| `replicate` | Replicate / source-study index within a case |
| `species` | Species the profile was measured in |
| `route` | Route of administration (IV, IRT, IRS) |
| `dose` | Administered dose |
| `ECCS_class` | Predicted Extended Clearance Classification System class |
| `time` | Time after dose |
| `concentration` | Mean plasma concentration |
| `concentration_lower` / `concentration_upper` | Reported uncertainty bounds (SEM or 95% CI), where available |

Derived PK endpoints provided per case: **Cmax**, **Tmax**, **Thalf**,
**AUCinf**, **AUCt**, **CL**, **Vd**, and **Fa**.

> Check units in the file itself — they are carried alongside the values and
> should not be assumed uniform across cases.

## Provenance

Concentration–time profiles were extracted from published figures using
[PlotDigitizer](https://plotdigitizer.com), then manually curated curve by
curve. Bounds reflect what the original publication reported; absence of a
bound means the source did not provide one, not that uncertainty is zero.

## Intended use

The dataset is meant as a reproducible benchmark for PK prediction methods —
comparing predicted against observed *full profiles* rather than only summary
endpoints, which is the point of publishing the curves rather than just Cmax
and AUC. Typical uses: scoring PBPK or high-throughput mechanistic PK
workflows, cross-species comparison, and absorption-model evaluation.

## Citation

If you use this dataset, please cite the preprint above. A BibTeX entry:

```bibtex
@article{jones2025pkbenchmark,
  title  = {A curated benchmark database for high-throughput mechanistic
            pharmacokinetic prediction},
  author = {Jones, Jeremy O. and Bachorz, Rafa{\l} A. and Lawless, Michael S.
            and Miller, David W.},
  year   = {2025},
  doi    = {10.21203/rs.3.rs-10671182/v1},
  note   = {Research Square preprint}
}
```

## License

TBD
<!-- Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
matching the license of the associated preprint. -->

## Contact

Maintained by [Simulations Plus, Inc.](https://www.simulations-plus.com) Please
open an issue for corrections, additions, or questions about provenance of a
particular curve.

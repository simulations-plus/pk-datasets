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
| Unique compound identifiers | 75 |
| Species labels in the CSV | dog, human, monkey, mouse, rat |
| Routes | intravenous (IV) , immediate-release tablet (IRT), immediate-release solution (IRS) |

Each case contains digitized concentration–time observations, a SMILES
structure, dose metadata, and upper/lower bounds where reported by the source.

## Repository structure

```
.
├── README.md
├── LICENSE
└── data
  ├── Cptimecurves_machine_readable.csv
  └── REFERENCES.md
```

## Data dictionary

Columns of `data/Cptimecurves_machine_readable.csv`:

| Column | Description |
|---|---|
| `Compound` | Compound identifier |
| `Smiles` | Compound structure in SMILES format |
| `Dose (mg)` / `Dose (mg/kg)` | Administered dose in milligrams and milligrams per kilogram |
| `Route` | Route label (`IV`, `IRT`, or `IRS`) |
| `Species` | Species label (`dog`, `human`, `monkey`, `mouse`, `rat`); monkey is a combination of data from Macaca fascicularis and Macaca mulatta |
| `Source` | Source label mapped to the full reference in [REFERENCES.md](data/REFERENCES.md) |
| `Time (h)` | Time after dose in hours |
| `Concentration (ng/mL)` | Plasma concentration in nanograms per milliliter |
| `Replicate` | Replicate index within a case |
| `Lower bound` / `Upper bound` | Concentration bounds in nanograms per milliliter, where available |
| `ADMET_Risk` | ADMET risk score as calculated in ADMET Predictor v14|
| `ECCS_Class` | Extended Clearance Classification System class as calculated in ADMET Predictor v14 |
| `RuleOf5`, `CyclPept_filter`, `Macrocycle`, `PROTAC_filter` | Compound classification or filter fields as calculated in ADMET Predictor v14 |

This CSV contains concentration–time observations and metadata; derived PK
endpoints such as **Cmax**, **Tmax**, **AUC**, **CL**, and **Vd** are not included.

## Provenance

Plots were digitized using PlotDigitizer (plotdigitizer.com), followed by manual review of all curves against the source publications by the researcher who performed the extraction. Times and concentrations were adjusted when necessary. Bounds reflect what the original publication reported; absence of a bound means the source did
not provide one, not that uncertainty is zero.

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

This project is licensed under the Creative Commons Attribution 4.0 International License - see the [LICENSE](LICENSE) file for details.

[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: https://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://licensebuttons.net/l/by/4.0/88x31.png

## Contact

Maintained by [Simulations Plus, Inc.](https://www.simulations-plus.com) Please
open an issue for corrections, additions, or questions about provenance of a
particular curve.

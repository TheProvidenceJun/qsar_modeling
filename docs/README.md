# QSAR Modeling for Soil Sorption Coefficient (logKoc)

![Status](https://img.shields.io/badge/Status-Manuscript%20Hardened-brightgreen)
![Python](https://img.shields.io/badge/Python-3.x-blue)
![Machine Learning](https://img.shields.io/badge/Task-QSAR%20Modeling-success)

## Project Overview

The project develops, validates, and documents a quantitative structure-activity relationship pipeline for predicting the soil organic carbon-water partition coefficient, \(logK_{oc}\), for environmentally relevant organic chemicals. The finalized workflow uses QDB.177 from Gramatica et al. (2014) as a controlled historical benchmark and evaluates whether modern descriptor clustering and parsimonious feature selection can exceed the published external predictivity while preserving the identical eight-descriptor constraint.

- **Target Endpoint:** \(logK_{oc}\)
- **Reference Dataset:** QDB.177 associated with Gramatica et al. (2014)
- **Validated Dataset Size:** Six hundred forty-two unique curated compounds
- **Literature Benchmark:** Gramatica et al. (2014) Model Four external \(R^2_{ext} = 0.794\)
- **Champion Model:** `Hierarchical_MC_MLR`
- **Champion Algorithm:** Multiple Linear Regression trained on eight locked molecular descriptors
- **Modeling Standard:** Training-set-only filtering, cross-validation-driven model selection, external validation, MCCV, Y-scrambling, and applicability-domain analysis under the five OECD validation principles

## Workflow & Pipeline

The computational workflow follows a strict anti-leakage design in which all filtering and model-selection decisions were fitted within the training space before external test-set evaluation.

- **Data Preprocessing:** QDB.177 endpoint parsing, RDKit-based SMILES validation, salt removal, canonicalization, and duplicate exclusion secured 642 unique valid compounds.
- **Rational Data Splitting:** Y-ranking split preserved the response-domain distribution and generated 514 training compounds plus 128 external test compounds with zero canonical SMILES overlap.
- **Feature Extraction:** Mordred generated 1,613 raw two-dimensional molecular descriptors, and native PyQSAR3-compatible preprocessing reduced the matrix to 377 robust variables.
- **Feature Clustering and Selection:** Descriptor-wise Hierarchical, K-Means, and SOM feature clustering structured the molecular descriptor space before Monte Carlo feature selection identified the final eight-descriptor subset.
- **Model Development and Validation:** Linear and non-linear candidate models were screened under a strict \(Q^2_{cv}\)-driven selection policy, followed by MCCV, Y-scrambling, and Williams Plot applicability-domain verification.

## Environments & Installation

The project uses separated environments to isolate descriptor calculation, PyQSAR3 modeling, and server-side non-linear model execution.

```bash
# Mordred descriptor calculation environment
conda env create -f envs/env_mordred.yml
conda activate mordred

# PyQSAR3 modeling environment
conda env create -f envs/env_pyqsar3.yml
conda activate pq3
```

The server-side non-linear and validation scripts were designed for a dedicated `qsar_ml` environment containing stable scientific Python dependencies, including `pandas`, `numpy`, `scikit-learn`, `deap`, `joblib`, and `matplotlib`.

## Directory Structure

```text
qsar_modeling/
├── data/           # Curated datasets, descriptor matrices, metrics, figures, and validation outputs
├── envs/           # Conda environment definition files
├── notebooks/      # Phase-based Jupyter notebooks for data preparation, descriptors, modeling, and metrics
├── docs/           # Manuscript draft, project roadmap, handoff state, and technical README
├── examples/       # Mordred and PyQSAR3 reference examples
├── run_nonlinear_models.py
├── run_mccv.py
└── run_yscrambling.py
```

## Final Results & Conclusion

The finalized champion model, `Hierarchical_MC_MLR`, achieved stronger external predictivity than the historical Gramatica et al. (2014) benchmark while retaining the same eight-descriptor budget and avoiding additional mathematical complexity.

| Model | Descriptor Count | \(R^2_{train}\) | \(Q^2_{cv}\) | External \(Q^2_{ext\ F2}\) | MCCV Mean \(R^2_{ext}\) | Mean \(R^2_{y-sc}\) |
|---|---:|---:|---:|---:|---:|---:|
| Gramatica et al. (2014) Model Four | 8 | 0.790 | 0.780 | 0.794 | Not reported | Not reported |
| `Hierarchical_MC_MLR` | 8 | 0.820 | 0.812 | 0.814 | 0.806 | 0.016 |

The locked champion descriptor set contains the following molecular variables:

- `ABC`
- `BCUTs-1h`
- `C1SP2`
- `ETA_shape_y`
- `FilterItLogS`
- `NdS`
- `SlogP_VSA1`
- `SlogP_VSA2`

The validation portfolio supports the scientific conclusion that the champion model captures a stable and interpretable structure-property relationship for soil sorption. MCCV confirmed partition-independent generalizability, Y-scrambling rejected chance correlation, and the Williams Plot verified applicability-domain compliance for the external validation chemicals.

## Peer-Review Hardening Roadmap

The manuscript hardening plan records four additional validation initiatives designed to address foreseeable reviewer objections and strengthen regulatory defensibility.

| Vulnerability Vector | Planned Hardening Action | Scientific Purpose |
|---|---|---|
| Low Silhouette score for hierarchical feature clustering | Execute a clustered versus non-clustered ablation study | Quantify whether descriptor clustering improves model efficiency and predictive robustness beyond direct feature selection |
| Modest external \(R^2\) gain over the historical benchmark | Compile RMSE and MAE error reductions with statistical significance testing | Demonstrate whether the observed performance delta supports the additional workflow structure in regulatory comparisons |
| Potential post-hoc descriptor interpretation bias | Run cross-model descriptor consensus analysis across top-ranked models | Establish whether selected molecular features remain stable across adjacent high-performing modeling paths |
| Historical dataset size constraint | Frame QDB.177 as a controlled benchmarking vehicle | Clarify that the study evaluates algorithm topology under identical literature constraints rather than claiming large-scale data coverage |

## Replication Instructions

The verified production-grade Phase Four source code is embedded in Section Six of `docs/Paper_Draft.md`. The Appendix code retrains the locked `Hierarchical_MC_MLR` champion model, generates the predicted-versus-experimental figure, generates the Williams Plot, and compiles the final comprehensive metrics table.

Recommended replication sequence:

```bash
conda activate pq3
jupyter nbconvert --to notebook --execute notebooks/04_result_metrics.ipynb --inplace
```

The primary output artifacts are:

- `data/features/figure_predicted_vs_experimental.png`
- `data/features/figure_williams_plot.png`
- `data/features/final_comprehensive_metrics.csv`
- `data/features/benchmark_oecd_comparison.csv`

## References

- Gramatica, P.; Cassani, S.; Chirico, N. QSARINS-Chem: Insubria Datasets and New QSAR/QSPR Models for Environmental Pollutants in QSARINS. *J. Comput. Chem.* **2014**, *35* (13), 1036–1044.
- QSAR DataBank. QDB.177 Archive. DOI: 10.15152/QDB.177.
- Organisation for Economic Co-operation and Development. *Guidance Document on the Validation of (Quantitative) Structure-Activity Relationship [(Q)SAR] Models*; OECD Series on Testing and Assessment, Number 69; OECD Environment Directorate: Paris, 2007.

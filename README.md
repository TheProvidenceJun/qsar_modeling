# QSAR logKoc Modeling Project

## Project Overview

This repository contains a curated quantitative structure-activity relationship workflow for predicting the soil organic carbon-water partition coefficient, $logK_{oc}$. The project uses the QDB.177 dataset associated with Gramatica et al. (2014) as a controlled benchmarking vehicle for evaluating descriptor organization, feature selection, and model validation under fixed historical constraints.

The final validated model is `Hierarchical_MC_MLR`, a multiple linear regression model built from hierarchical descriptor clustering and Monte Carlo feature selection. The workflow preserves the identical eight-descriptor budget used by the Gramatica 2014 Model 4 benchmark while improving external predictive performance and maintaining full OECD-style validation evidence.

## Core Pipeline Architecture

The production modeling chain follows a strict leakage-controlled design:

`RDKit Data Curation` $\rightarrow$ `Y-Ranking Split` $\rightarrow$ `Mordred Descriptor Calculation` $\rightarrow$ `Training-Only Pre-filtering` $\rightarrow$ `Hierarchical Feature Clustering` $\rightarrow$ `Monte Carlo Feature Selection` $\rightarrow$ `Multiple Linear Regression`

Key pipeline dimensions:

| Stage | Output |
| --- | --- |
| RDKit structural curation | 642 unique valid compounds |
| Y-ranking split | 514 training compounds and 128 external test compounds |
| Mordred descriptor calculation | 1,613 raw two-dimensional descriptor variables |
| Strict train-set-only filtering | 377 robust descriptor variables |
| Feature-space organization | Hierarchical descriptor clustering |
| Final feature selection | Monte Carlo selection of eight descriptors |
| Final regression engine | Multiple Linear Regression |

## Champion Model Performance

The final champion model, `Hierarchical_MC_MLR`, was selected strictly by the highest cross-validation robustness score, $Q^2_{cv}$, across the evaluated model matrix. External validation metrics were used only for final verification to prevent data snooping.

| Model | $R^2_{train}$ | $Q^2_{cv}$ | External $Q^2_{ext\ F2}$ | MCCV Mean $R^2_{ext}$ | Mean $R^2_{y-sc}$ |
| --- | --- | --- | --- | --- | --- |
| Gramatica 2014 Model 4 | 0.790 | 0.780 | 0.794 | Not reported | Not reported |
| `Hierarchical_MC_MLR` | 0.820 | 0.812 | 0.814 | 0.806 | 0.016 |

The champion model exceeded the Gramatica 2014 external benchmark under the identical eight-descriptor constraint. The 100-iteration Monte Carlo Cross-Validation result showed stable predictive behavior, and the 100-iteration Y-scrambling result rejected chance correlation.

## Locked Champion Descriptors

The final eight selected descriptors are:

- `ABC`
- `BCUTs-1h`
- `C1SP2`
- `ETA_shape_y`
- `FilterItLogS`
- `NdS`
- `SlogP_VSA1`
- `SlogP_VSA2`

These descriptors encode molecular size, skeletal branching, conjugation, aqueous solubility, hydrophobic surface distribution, sulfur composition, and electronic interaction capacity. The final descriptor set supports mechanistic interpretation under OECD Principle 5.

## Repository Outputs

Primary manuscript and documentation artifacts:

- `docs/Paper_Draft.md`
- `docs/PROJECT.md`
- `docs/HANDOFF.md`

Primary modeling and validation artifacts:

- `data/features/final_comprehensive_metrics.csv`
- `data/features/benchmark_oecd_comparison.csv`
- `data/features/mccv_summary.json`
- `data/features/yscrambling_summary.json`
- `data/features/figure_predicted_vs_experimental.png`
- `data/features/figure_williams_plot.png`

## How to Run

The verified production-grade Python workflow for final metric integration and diagnostic visualization is embedded in Section 6 of `docs/Paper_Draft.md`. Users seeking to reproduce the final figures and the comprehensive metric table should follow the Appendix code block and execute the workflow from the repository root or from the `notebooks/` directory after activating the project environment.

The Appendix workflow retrains the locked `Hierarchical_MC_MLR` model, generates the predicted-versus-experimental plot, generates the Williams Plot, and compiles the final 16-model comprehensive metrics table.

## Peer-Review Hardening Roadmap

The repository now tracks four follow-up analyses intended to strengthen peer-review defense:

- Feature-clustering ablation against a non-clustered baseline.
- Statistical testing of the external performance gain over the Gramatica 2014 benchmark.
- Cross-model descriptor-consensus analysis across high-ranking models.
- Explicit controlled-benchmarking framing for the 642-compound historical dataset.

## Scientific Position

The validated `Hierarchical_MC_MLR` workflow demonstrates that careful descriptor curation, feature-space organization, and Monte Carlo selection can improve $logK_{oc}$ predictivity without increasing model complexity. The final model remains parsimonious, interpretable, externally predictive, stable under repeated resampling, resistant to chance correlation, and aligned with OECD validation principles.

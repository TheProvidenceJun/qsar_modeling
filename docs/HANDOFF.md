# Project Handoff

## Current Status
Phase 5 - Step 5.2 is complete, and the project is now in Phase 5 - Step 5.3 for Results and Discussion compilation.
- The QDB.177 `M2.logKoc` endpoint was parsed from `data/raw/properties/M2.logKoc/values`.
- Raw compound identifiers were mapped to SMILES files under `data/raw/compounds/<Compound Id>/daylight-smiles`.
- RDKit validation and canonicalization secured 642 unique valid compounds from 643 raw entries.
- A Y-Ranking 80/20 split was generated with zero SMILES overlap between Training and Test sets.
- Raw 2D Mordred descriptors were extracted independently for the Training and Test sets.
- The raw descriptor feature space was preserved without missing-value filtering, zero-variance filtering, or correlation filtering.
- Native PyQSAR3 preprocessing reduced the modeling tables to 379 columns while preserving strict train-set-only filtering.
- The 377 filtered descriptors were clustered variable-wise using native PyQSAR3 `FeatureCluster` modules, producing 20-25 descriptor clusters for the upcoming GA/MC feature-selection engines.
- The 12-model linear baseline matrix was regenerated under the strict Gramatica 2014 constraint of exactly 8 descriptors per model.
- Comprehensive internal and external metrics were calculated for all 12 linear baselines.
- The server-side non-linear modeling script evaluated four Hierarchical-track non-linear candidates: GA-SVR, MC-SVR, GA-RF, and MC-RF.
- All 16 candidate models have now been generated and evaluated: 12 linear baselines plus 4 non-linear models.
- The rigorous `Q²_cv` comparison selected `Hierarchical_MC_MLR` as the Champion Model.
- The 100-iteration MCCV robustness test on the locked `Hierarchical_MC_MLR` Champion Model yielded an average `R²_ext` of 0.806.
- The 100-iteration Y-scrambling test yielded an average scrambled `R²_y-sc` of 0.016, confirming that the Champion Model is not based on chance correlation.
- Phase 4 generated high-resolution diagnostic plots for the Champion Model:
  - `data/features/figure_predicted_vs_experimental.png`
  - `data/features/figure_williams_plot.png`
- The Williams Plot established the warning leverage threshold at `h* = 0.0525` and confirmed broad applicability-domain compliance under OECD Principle 3.
- All 16 evaluated model records were integrated with MCCV and Y-scrambling summaries in `data/features/final_comprehensive_metrics.csv` with shape `(16, 33)`.
- Phase 5 manuscript rules have been sealed in `docs/PROJECT.md`, including formal academic tone, no raw numerals at sentence starts, tense consistency, explicit pronouns, disciplined logical transitions, and the requirement to build every paper expansion from the current Phase 1-4 baseline.
- Phase 5 - Step 5.1 generated `data/features/benchmark_oecd_comparison.csv`, confirming direct comparison against Gramatica 2014 Model 4.
- Phase 5 - Step 5.2 completed the manuscript Introduction and Materials and Methods sections in `docs/Paper_Draft.md`.

## Active Project
- **Project:** QSAR logKoc Modeling Project
- **Dataset:** QDB.177, Gramatica et al. 2014
- **Target property:** logKoc
- **Main objective:** Build and validate an ensemble/consensus QSAR model.
- **Original benchmark:** Test R² = 0.794
- **Benchmark descriptor count:** 8 descriptors

## Current Phase
**Phase 5 - Step 5.3: Results and Discussion Compiling**
This is the next action.

## Next Action
Proceed to **Phase 5 - Step 5.3: Results and Discussion Compiling**.

Required Step 5.3 task:
- Compile Section 3, `Results and Discussion`, in `docs/Paper_Draft.md`.
- Discuss descriptor clustering performance, including Silhouette scores for Hierarchical, SOM, and K-Means tracks.
- Discuss the 16-model screening matrix ranked strictly by `Q²_cv`.
- Explain the selection of `Hierarchical_MC_MLR` as the champion model and contrast its robust linear behavior with the overfitting tendency of high-capacity RF/SVR models.
- Report MCCV stability, Y-scrambling chance-correlation rejection, Williams Plot applicability-domain interpretation, and the direct Gramatica 2014 benchmark comparison.

Mandatory Phase 5 manuscript protocol:
- Read the full current `docs/Paper_Draft.md` before any manuscript edit.
- Preserve the validated Phase 1-4 baseline and build cleanly from it.
- Apply the strict academic formatting rules recorded in `docs/PROJECT.md`.

Locked final Champion Model for manuscript compilation:
- Model: `Hierarchical_MC_MLR`
- Algorithm: MLR
- Descriptors: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, `SlogP_VSA2`

## Completed Tasks
1. Parsed the QDB.177 logKoc endpoint (`M2.logKoc`) and constructed a unified dataset containing `SMILES` and `logKoc`.
2. Validated and canonicalized SMILES using RDKit.
3. Removed invalid/disconnected structures and duplicate canonical SMILES. The raw set contained 643 entries; 1 duplicate was removed; 642 compounds were secured for modeling.
4. Applied Y-Ranking by sorting compounds according to `logKoc` and assigning every fifth ranked compound to the external Test set.
5. Saved the processed datasets:
   - `data/processed/train.csv` with 514 compounds
   - `data/processed/test.csv` with 128 compounds
6. Confirmed zero data leakage between Training and Test sets based on canonical SMILES overlap.
7. Calculated raw 2D Mordred descriptors using `Calculator(descriptors, ignore_3D=True)`.
8. Extracted 1,613 Mordred descriptors for both splits without preliminary filtering.
9. Saved the raw descriptor matrices:
   - `data/features/mordred_train.csv` with shape `(514, 1615)`
   - `data/features/mordred_test.csv` with shape `(128, 1615)`
10. Completed native PyQSAR3 pre-filtering using `pyqsar.data_tools`:
    - Applied `pyqsar.data_tools.NonNumricFilter`.
    - Applied `FilteringTools.rm_nan`.
    - Applied `FilteringTools.rm_inf`.
    - Applied `FilteringTools.rm_novar`.
    - Applied a strict train-set-only Pearson high-correlation filter (`|r| > 0.90`) because PyQSAR3 does not expose a native high-correlation pre-filter.
11. Saved the PyQSAR3-filtered descriptor matrices:
    - `data/features/filtered_train_pyqsar3.csv` with shape `(514, 379)`
    - `data/features/filtered_test_pyqsar3.csv` with shape `(128, 379)`
12. Corrected Step 3.2 from sample clustering to feature clustering, aligning the workflow with PyQSAR3 GA/MC requirements.
13. Clustered all 377 filtered descriptors using native PyQSAR3 feature-clustering modules:
    - Hierarchical: `pyqsar.model_tools.FeatureCluster`, tuned distance cut, 20 feature clusters
    - K-Means: `pyqsar.model_tools.FeatureCluster_KMeans`, configured to 20 feature clusters
    - SOM: `pyqsar.model_tools.FeatureCluster_Minisom`, 5x5 map, 25 feature clusters
14. Saved descriptor-to-cluster mappings and native sidecar files:
    - `data/features/feature_clusters_pyqsar3.json`
    - `data/features/feature_clusters_pyqsar3_hierarchical.cluster`
    - `data/features/feature_clusters_pyqsar3_kmeans.cluster`
    - `data/features/feature_clusters_pyqsar3_som.cluster`
15. Deleted the obsolete row-wise clustering outputs because they were incompatible with PyQSAR3 feature-selection engines.
16. Completed the strict 8-descriptor 12-model linear baseline matrix:
    - Hierarchical FeatureCluster: GA-MLR, GA-PLS, MC-MLR, MC-PLS
    - K-Means FeatureCluster: GA-MLR, GA-PLS, MC-MLR, MC-PLS
    - SOM FeatureCluster: GA-MLR, GA-PLS, MC-MLR, MC-PLS
17. Calculated comprehensive internal and external metrics for all 12 linear baseline models.
18. Saved Step 3.3 outputs:
    - `data/features/12_model_extended_metrics.csv`
    - `data/features/final_selected_features_8vars.json`
    - `data/features/test_predictions_8vars.csv`
19. Generated the server-side non-linear modeling script:
    - `run_nonlinear_models.py`
20. Completed Step 3.4 server-side non-linear modeling on the Hierarchical feature-cluster track:
    - Hierarchical-GA-SVR: `Q²_cv = 0.811`, `Q²_ext F2 = 0.809`
    - Hierarchical-MC-SVR: `Q²_cv = 0.783`, `Q²_ext F2 = 0.764`
    - Hierarchical-GA-RF: `Q²_cv = 0.794`, `Q²_ext F2 = 0.826`
    - Hierarchical-MC-RF: `Q²_cv = 0.787`, `Q²_ext F2 = 0.803`
21. Saved Step 3.4 outputs:
    - `data/features/nonlinear_search_log.txt`
    - `data/features/best_nonlinear_config.json`
    - `data/features/nonlinear_model_metrics.csv`
22. Completed strict best-model selection across all 16 evaluated models using `Q²_cv` only:
    - Champion Model: `Hierarchical_MC_MLR`
    - Algorithm: MLR
    - `Q²_cv = 0.8119744761107354`
    - Selected descriptors: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, `SlogP_VSA2`
23. Generated the standalone MCCV script:
    - `run_mccv.py`
24. Completed 100-iteration MCCV for the locked Champion Model:
    - Champion Model: `Hierarchical_MC_MLR`
    - Algorithm: MLR
    - Champion `Q²_cv = 0.812`
    - Champion descriptors: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, `SlogP_VSA2`
    - MCCV `R²_ext`: mean = 0.806, SD = 0.032
    - MCCV `RMSE_ext`: mean = 0.518, SD = 0.043
    - MCCV `MAE_ext`: mean = 0.407, SD = 0.033
25. Generated the standalone Y-scrambling script for the locked Champion Model:
    - `run_yscrambling.py`
26. Completed 100-iteration Y-scrambling chance-correlation test:
    - Champion Model: `Hierarchical_MC_MLR`
    - Original `R²_train = 0.820`
    - Average scrambled `R²_y-sc = 0.016`
    - Scrambled `R²_y-sc` SD = 0.007
    - Maximum scrambled `R² = 0.039`
    - Conclusion: chance correlation is mathematically rejected; the model relies on authentic structure-property relationships.
27. Completed Phase 4 predicted-versus-experimental visualization for the Champion Model:
    - Output: `data/features/figure_predicted_vs_experimental.png`
    - Training set: `n = 514`
    - External test set: `n = 128`
    - Diagnostic interpretation: predictions were symmetrically distributed around the `y = x` line, with no systematic bias or evident heteroscedasticity.
28. Completed Phase 4 Williams Plot applicability-domain analysis:
    - Output: `data/features/figure_williams_plot.png`
    - Warning leverage threshold: `h* = 0.0525`
    - Standardized residual bounds: `[-3, +3]`
    - Diagnostic interpretation: the vast majority of training and external validation compounds fell within the accepted response and leverage domain.
29. Completed comprehensive Phase 4 metrics integration:
    - Output: `data/features/final_comprehensive_metrics.csv`
    - Shape: `(16, 33)`
    - Content: 12 linear models, 4 non-linear models, MCCV summaries, and Y-scrambling summaries in one unified master table.
30. Completed Phase 5 - Step 5.1 quantitative benchmark and OECD validation coding:
    - Output: `data/features/benchmark_oecd_comparison.csv`
    - Historical baseline: Gramatica 2014 Model 4, `R²_train = 0.790`, `Q²_cv = 0.780`, external `R²_ext = 0.794`
    - Champion model: `Hierarchical_MC_MLR`, `R²_train = 0.820182`, `Q²_cv = 0.811974`, external `Q²_ext F2 = 0.814252`
31. Completed Phase 5 - Step 5.2 manuscript drafting:
    - Section 1, `Introduction`, was authored.
    - Section 2, `Materials and Methods`, was refined through Sections 2.1-2.5.
    - Sections 3-6 were preserved as placeholders for subsequent manuscript compilation steps.

## Phase 1 Outputs
- `data/processed/train.csv`
- `data/processed/test.csv`

## Phase 2 Outputs
- `data/features/mordred_train.csv`
- `data/features/mordred_test.csv`

## Phase 3 Step 3.1 Outputs
- `data/features/filtered_train_pyqsar3.csv`
- `data/features/filtered_test_pyqsar3.csv`

## Phase 3 Step 3.2 Outputs
- `data/features/feature_clusters_pyqsar3.json`
- `data/features/feature_clusters_pyqsar3_hierarchical.cluster`
- `data/features/feature_clusters_pyqsar3_kmeans.cluster`
- `data/features/feature_clusters_pyqsar3_som.cluster`

## Phase 3 Step 3.3 Outputs
- `data/features/12_model_extended_metrics.csv`
- `data/features/final_selected_features_8vars.json`
- `data/features/test_predictions_8vars.csv`

## Phase 3 Step 3.4 Outputs
- `run_nonlinear_models.py`
- `data/features/nonlinear_search_log.txt`
- `data/features/best_nonlinear_config.json`
- `data/features/nonlinear_model_metrics.csv`

## Phase 3 Step 3.5 Outputs
- `run_mccv.py`
- `data/features/mccv_100_iterations.csv`
- `data/features/mccv_summary.json`

## Phase 3 Step 3.6 Outputs
- `run_yscrambling.py`
- `data/features/yscrambling_100_iterations.csv`
- `data/features/yscrambling_summary.json`

## Phase 4 Outputs
- `data/features/figure_predicted_vs_experimental.png`
- `data/features/figure_williams_plot.png`
- `data/features/final_comprehensive_metrics.csv`

## Phase 5 Expected Outputs
- `data/features/benchmark_oecd_comparison.csv`
- Completed `docs/Paper_Draft.md` with Introduction, Discussion, Conclusion, and References finalized.
- Final benchmark comparison against Gramatica et al. 2014.
- Repository-facing workflow documentation prepared for project closure.

## Environment For Next Step
Recommended environment:
```bash
conda activate qsar_ml
```

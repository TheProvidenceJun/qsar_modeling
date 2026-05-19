# Project Handoff

## Current Status
Phase 3 - Step 3.6 script generation is complete. The Y-scrambling chance-correlation script has been created for the locked Champion Model.
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
- A standalone Y-scrambling script has been generated to test the locked Champion Model against chance correlation.

## Active Project
- **Project:** QSAR logKoc Modeling Project
- **Dataset:** QDB.177, Gramatica et al. 2014
- **Target property:** logKoc
- **Main objective:** Build and validate an ensemble/consensus QSAR model.
- **Original benchmark:** Test R² = 0.794
- **Benchmark descriptor count:** 8 descriptors

## Current Phase
**Phase 3: Advanced Feature Selection & Modeling (Non-linear & Validation)**
This is the next action.

## Next Action
The Project Manager will execute the standalone Y-scrambling script:
`run_yscrambling.py`

Execution command:
```bash
conda activate qsar_ml
python run_yscrambling.py
```

The script will:
- Use the locked Champion Model: `Hierarchical_MC_MLR`.
- Use the locked 8 descriptors: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, `SlogP_VSA2`.
- Fit the original MLR model on the unscrambled training target and calculate `R²_train`.
- Randomly permute the training `logKoc` values across 100 trials.
- Refit the MLR model for each scrambled response vector and calculate `R²_y-sc`.
- Save `data/features/yscrambling_100_iterations.csv`.
- Save `data/features/yscrambling_summary.json`.

Do not report Y-scrambling results until the script has been executed and the actual outputs have been reviewed.

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

## Phase 3 Step 3.6 Script Output
- `run_yscrambling.py`

## Expected Phase 3 Step 3.6 Runtime Outputs
- `data/features/yscrambling_100_iterations.csv`
- `data/features/yscrambling_summary.json`

## Environment For Next Step
Recommended environment:
```bash
conda activate qsar_ml
```

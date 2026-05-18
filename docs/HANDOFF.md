# Project Handoff

## Current Status
Phase 3 - Step 3.2b: Multi-Track Data Partitioning has been completed successfully.
- The QDB.177 `M2.logKoc` endpoint was parsed from `data/raw/properties/M2.logKoc/values`.
- Raw compound identifiers were mapped to SMILES files under `data/raw/compounds/<Compound Id>/daylight-smiles`.
- RDKit validation and canonicalization secured 642 unique valid compounds from 643 raw entries.
- A Y-Ranking 80/20 split was generated with zero SMILES overlap between Training and Test sets.
- Raw 2D Mordred descriptors were extracted independently for the Training and Test sets.
- The raw descriptor feature space was preserved without missing-value filtering, zero-variance filtering, or correlation filtering.
- Strict pre-filtering was fitted only on the Training set and then applied identically to the Test set.
- The filtered feature space now contains 379 robust descriptors for downstream feature selection.
- The clustered Training data has been physically partitioned into three independent track files for stable PyQSAR3 execution.

## Active Project
- **Project:** QSAR logKoc Modeling Project
- **Dataset:** QDB.177, Gramatica et al. 2014
- **Target property:** logKoc
- **Main objective:** Build and validate an ensemble/consensus QSAR model.
- **Original benchmark:** Test R² = 0.794

## Current Phase
**Phase 3: Feature Filtering & PyQSAR3 Modeling**
This is the next action.

## Next Action
Step 3.3: Execute the 12-Model PyQSAR3 Matrix in the Phase 3 notebook:
`notebooks/03_pyqsar3_modeling.ipynb`

Execution plan:
1. Feed the three physically separated cluster-track datasets into the coupled PyQSAR3 feature-selection/modeling engines.
2. Execute the full 12-model matrix:
   - Hierarchical Track: Cluster-GA (MLR/PLS) and Cluster-MC (MLR/PLS)
   - K-Means Track: Cluster-GA (MLR/PLS) and Cluster-MC (MLR/PLS)
   - SOM Track: Cluster-GA (MLR/PLS) and Cluster-MC (MLR/PLS)
3. Pause consensus modeling for now; focus only on exploring and recording the 12 base-model pathways.

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
10. Completed strict Training-set-fitted pre-filtering to prevent data leakage:
    - Removed 491 descriptors containing missing values in the Training set.
    - Removed 165 zero-variance/constant descriptors using `VarianceThreshold(threshold=0.0)` fitted on the Training set.
    - Removed 578 highly correlated descriptors using absolute Pearson correlation > 0.90 computed on the Training set.
    - Retained 379 descriptors for both Training and Test sets.
11. Saved the filtered descriptor matrices:
    - `data/features/filtered_train.csv` with shape `(514, 381)`
    - `data/features/filtered_test.csv` with shape `(128, 381)`
12. Completed Training-set clustering and saved per-compound cluster labels:
    - `data/features/cluster_labels_train.csv`
13. Physically partitioned the clustered Training data into three PyQSAR3-ready track files. Each file contains `SMILES`, `logKoc`, 379 filtered descriptors, and one standardized `Cluster_Label` column:
    - `data/features/train_track_hierarchical.csv` with shape `(514, 382)`
    - `data/features/train_track_kmeans.csv` with shape `(514, 382)`
    - `data/features/train_track_som.csv` with shape `(514, 382)`

## Phase 1 Outputs
- `data/processed/train.csv`
- `data/processed/test.csv`

## Phase 2 Outputs
- `data/features/mordred_train.csv`
- `data/features/mordred_test.csv`

## Phase 3 Step 3.1 Outputs
- `data/features/filtered_train.csv`
- `data/features/filtered_test.csv`

## Phase 3 Step 3.2/3.2b Outputs
- `data/features/cluster_labels_train.csv`
- `data/features/train_track_hierarchical.csv`
- `data/features/train_track_kmeans.csv`
- `data/features/train_track_som.csv`

## Environment For Next Step
Recommended environment:
```bash
conda activate pq3
```

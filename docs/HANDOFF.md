# Project Handoff

## Current Status
Phase 3 - Step 3.1: Native PyQSAR3 Pre-Filtering has been completed successfully.
- The QDB.177 `M2.logKoc` endpoint was parsed from `data/raw/properties/M2.logKoc/values`.
- Raw compound identifiers were mapped to SMILES files under `data/raw/compounds/<Compound Id>/daylight-smiles`.
- RDKit validation and canonicalization secured 642 unique valid compounds from 643 raw entries.
- A Y-Ranking 80/20 split was generated with zero SMILES overlap between Training and Test sets.
- Raw 2D Mordred descriptors were extracted independently for the Training and Test sets.
- The raw descriptor feature space was preserved without missing-value filtering, zero-variance filtering, or correlation filtering.
- Native PyQSAR3 preprocessing reduced the modeling tables to 379 columns while preserving strict train-set-only filtering.

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
Step 3.2: Native PyQSAR3 Clustering & Silhouette Scoring in the Phase 3 notebook:
`notebooks/03_pyqsar3_modeling.ipynb`

Use the PyQSAR3-filtered outputs from Step 3.1 as the starting point and continue with native PyQSAR3-compatible clustering.

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

## Phase 1 Outputs
- `data/processed/train.csv`
- `data/processed/test.csv`

## Phase 2 Outputs
- `data/features/mordred_train.csv`
- `data/features/mordred_test.csv`

## Phase 3 Step 3.1 Outputs
- `data/features/filtered_train_pyqsar3.csv`
- `data/features/filtered_test_pyqsar3.csv`

## Environment For Next Step
Recommended environment:
```bash
conda activate pq3
```

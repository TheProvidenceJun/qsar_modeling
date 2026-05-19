# Paper Draft

## 1. Introduction

## 2. Materials and Methods

### 2.1. Data Collection and Curation

The modeling dataset was derived from QDB.177, the QSAR DataBank archive associated with the logKoc study reported by Gramatica et al. (2014). The target endpoint was the soil sorption coefficient expressed as logKoc, registered in the archive as `M2.logKoc`. Experimental response values were extracted from the corresponding QDB property table, and compound identifiers were used to retrieve the associated Daylight SMILES representations from the compound registry files.

Molecular structures were curated using RDKit prior to descriptor calculation and model development. Each SMILES string was parsed into an RDKit molecular object to verify chemical validity, and valid structures were converted to canonical SMILES to provide a unique representation for duplicate detection. Disconnected structures, salts, and mixtures were excluded from the modeling set to retain a single interpretable molecular graph for each observation. Duplicate records were then removed on the basis of canonical SMILES. The initial QDB endpoint contained 643 raw entries. After RDKit-based standardization and duplicate removal, 642 unique valid substances were retained for subsequent QSAR analysis.

### 2.2. Dataset Splitting

The curated dataset was divided into training and external test subsets using a Y-Ranking strategy based on the response variable. Compounds were first sorted in ascending order of their experimental logKoc values. The ranked list was then sampled systematically so that every fifth compound was assigned to the external test set, while the remaining compounds were retained for model training. This response-sorted allocation was selected to promote uniform coverage of the logKoc domain in both subsets and to reduce the risk that the external test set would be concentrated in a narrow response interval.

The final split produced a training set containing 514 compounds and an external test set containing 128 compounds, corresponding to an approximately 80/20 partition of the curated dataset. Canonical SMILES were compared across the two subsets after splitting, and no overlap was detected. The processed datasets were saved as `data/processed/train.csv` and `data/processed/test.csv` for all subsequent descriptor extraction, feature selection, model development, and external validation steps.

### 2.3. Molecular Descriptor Calculation

Molecular descriptors were calculated using the Mordred descriptor library. Descriptor extraction was performed separately for the training and external test sets using the Mordred `Calculator` initialized with `ignore_3D=True`, thereby restricting the calculation to two-dimensional descriptors and avoiding conformer-dependent three-dimensional descriptor generation. The resulting descriptor tables retained the original `SMILES` and `logKoc` columns together with the calculated molecular descriptors.

A total of 1,613 raw two-dimensional Mordred descriptors were generated for each compound in both subsets. The training descriptor table contained 514 compounds and 1,615 columns, comprising `SMILES`, `logKoc`, and 1,613 descriptor variables. The external test descriptor table contained 128 compounds with the same 1,615-column schema. These raw descriptor matrices were saved as `data/features/mordred_train.csv` and `data/features/mordred_test.csv`.

No preliminary descriptor filtering was applied during this stage. Missing descriptor values were not imputed or removed, zero-variance descriptors were not excluded, and correlation-based filtering was not performed. Retaining the complete initial descriptor space at the extraction stage was a deliberate methodological choice to prevent information from the external test set from influencing feature selection. All filtering and feature selection steps are therefore reserved for the subsequent modeling phase, where they can be fitted strictly on the training set and then applied to the external test set.

### 2.4. Data Pre-Filtering

Prior to feature selection and model development, a strict descriptor pre-filtering protocol was applied to remove invalid, non-informative, and redundant variables. The native preprocessing modules of the `pyqsar3` library (`pyqsar.data_tools`) were used to eliminate non-numeric, infinite, missing, and zero-variance descriptors. Specifically, `pyqsar.data_tools.NonNumricFilter` was used for native numeric descriptor screening, followed by `FilteringTools.rm_nan`, `FilteringTools.rm_inf`, and `FilteringTools.rm_novar` for missing-value, infinite-value, and zero-variance filtering.

To prevent information leakage, all statistical filtering decisions were fitted exclusively on the training set and then applied unchanged to the external test set. Because the installed `pyqsar3` preprocessing utilities do not provide a native high-correlation pre-filter, a Pearson correlation filter was subsequently applied using only the training set correlation matrix. Descriptor pairs with absolute correlation greater than 0.90 were treated as highly redundant, and one descriptor from each such pair was removed before applying the same retained descriptor list to the external test set.

The final filtered datasets consisted of 377 robust molecular descriptors, corresponding to 379 total columns after retaining `SMILES` and `logKoc`. These PyQSAR3-filtered matrices were saved as `data/features/filtered_train_pyqsar3.csv` and `data/features/filtered_test_pyqsar3.csv`, providing a noise-controlled feature space for downstream clustering and model development.

## 3. Results

## 4. Discussion

## 5. Conclusion

## 6. References

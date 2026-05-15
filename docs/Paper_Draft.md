# Paper Draft

## 1. Introduction

## 2. Materials and Methods

### 2.1. Data Collection and Curation

The modeling dataset was derived from QDB.177, the QSAR DataBank archive associated with the logKoc study reported by Gramatica et al. (2014). The target endpoint was the soil sorption coefficient expressed as logKoc, registered in the archive as `M2.logKoc`. Experimental response values were extracted from the corresponding QDB property table, and compound identifiers were used to retrieve the associated Daylight SMILES representations from the compound registry files.

Molecular structures were curated using RDKit prior to descriptor calculation and model development. Each SMILES string was parsed into an RDKit molecular object to verify chemical validity, and valid structures were converted to canonical SMILES to provide a unique representation for duplicate detection. Disconnected structures, salts, and mixtures were excluded from the modeling set to retain a single interpretable molecular graph for each observation. Duplicate records were then removed on the basis of canonical SMILES. The initial QDB endpoint contained 643 raw entries. After RDKit-based standardization and duplicate removal, 642 unique valid substances were retained for subsequent QSAR analysis.

### 2.2. Dataset Splitting

The curated dataset was divided into training and external test subsets using a Y-Ranking strategy based on the response variable. Compounds were first sorted in ascending order of their experimental logKoc values. The ranked list was then sampled systematically so that every fifth compound was assigned to the external test set, while the remaining compounds were retained for model training. This response-sorted allocation was selected to promote uniform coverage of the logKoc domain in both subsets and to reduce the risk that the external test set would be concentrated in a narrow response interval.

The final split produced a training set containing 514 compounds and an external test set containing 128 compounds, corresponding to an approximately 80/20 partition of the curated dataset. Canonical SMILES were compared across the two subsets after splitting, and no overlap was detected. The processed datasets were saved as `data/processed/train.csv` and `data/processed/test.csv` for all subsequent descriptor extraction, feature selection, model development, and external validation steps.

## 3. Results

## 4. Discussion

## 5. Conclusion

## 6. References

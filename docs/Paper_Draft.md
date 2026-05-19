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

### 2.5. Feature Space Clustering

Following descriptor pre-filtering, unsupervised clustering was applied to the descriptor variables rather than to the chemical samples. This feature-space clustering strategy was used to group highly correlated and structurally similar molecular descriptors, thereby mitigating multicollinearity and organizing the retained descriptor space for subsequent PyQSAR3 feature selection. The procedure was conducted on the 377 filtered descriptors from the training set and was implemented using the native feature-clustering modules available in `pyqsar.model_tools`.

Three complementary descriptor-clustering approaches were evaluated to provide alternative topological representations of the filtered feature space. Hierarchical descriptor clustering was performed using `FeatureCluster`, with a tuned distance cut yielding 20 feature clusters. K-Means descriptor clustering was performed using `FeatureCluster_KMeans`, configured to produce 20 feature clusters. Self-Organizing Map descriptor clustering was performed using `FeatureCluster_Minisom`; a 5x5 map was used, resulting in 25 feature clusters after assignment of the filtered descriptors to SOM units.

All 377 retained descriptors were mapped to their corresponding feature-cluster labels for each of the three clustering strategies. These descriptor groupings were saved as `data/features/feature_clusters_pyqsar3.json`, together with native PyQSAR3 `.cluster` sidecar files. The resulting feature-cluster topologies provide the structural basis for the subsequent parallel multi-track genetic algorithm (GA) and Monte Carlo (MC) feature-selection engines, without introducing any sample-wise clustering or external test-set information into the feature-selection workflow.

## 3. Results

### 3.1. Feature Selection and Internal Performance Evaluation

The baseline linear modeling stage was executed as a systematic 12-model matrix in which three descriptor-clustering spaces were coupled with two PyQSAR3 feature-selection engines and two linear regression algorithms. Hierarchical, K-Means, and Self-Organizing Map feature spaces were each supplied to the native PyQSAR3 Genetic Algorithm (GA) and Monte Carlo (MC) feature-selection engines. Each selected descriptor subset was then evaluated using multiple linear regression (MLR) and partial least squares (PLS), yielding twelve baseline linear models. In accordance with the descriptor count reported for the Gramatica et al. (2014) benchmark model, all feature-selection paths were strictly constrained to select exactly eight molecular descriptors.

The internally fitted models exhibited strong goodness-of-fit across the full matrix. The highest training coefficients of determination were observed for the Hierarchical-MC models (`R²_train = 0.820` for MLR and `R²_train = 0.818` for PLS), followed closely by the Hierarchical-GA and SOM-MC pathways. These values indicate that the eight-descriptor subsets retained substantial explanatory capacity despite the strong parsimony constraint imposed on the feature-selection process.

Internal robustness was evaluated using 5-fold cross-validation, with cross-validated predictions summarized through `Q²_cv`, `CCC_cv`, `RMSE_cv`, and `MAE_cv`. The Hierarchical-MC track demonstrated particularly high internal stability, reaching `Q²_cv = 0.812` and `CCC_cv = 0.897` for both MLR and PLS variants. The Hierarchical-GA track also showed strong robustness, with both MLR and PLS models achieving `Q²_cv = 0.807` and `CCC_cv = 0.894`. The consistency between training performance and cross-validated predictivity suggests that the descriptor-cluster-guided search produced stable eight-variable baseline models rather than overparameterized fitting artifacts.

### 3.2. External Predictivity of Baseline Linear Models

The twelve baseline linear models were next evaluated on the external test set using a comprehensive validation framework. External predictivity was quantified using three complementary external predictive coefficients (`Q²_ext F1`, `Q²_ext F2`, and `Q²_ext F3`), together with the external concordance correlation coefficient (`CCC_ext`), `RMSE_ext`, and `MAE_ext`. The full extended performance matrix is reported below, ranked by `Q²_ext F2`, which is equivalent to the standard external test-set coefficient of determination.

| Rank | Model | Q²_ext F2 | Q²_ext F1 | Q²_ext F3 | CCC_ext | RMSE_ext | MAE_ext | Q²_cv | CCC_cv | R²_train |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Hierarchical_GA_PLS | 0.820 | 0.820 | 0.825 | 0.910 | 0.502 | 0.405 | 0.807 | 0.894 | 0.816 |
| 2 | Hierarchical_GA_MLR | 0.820 | 0.820 | 0.824 | 0.910 | 0.503 | 0.404 | 0.807 | 0.894 | 0.817 |
| 3 | SOM_MC_PLS | 0.817 | 0.817 | 0.822 | 0.908 | 0.506 | 0.383 | 0.804 | 0.893 | 0.816 |
| 4 | Hierarchical_MC_PLS | 0.817 | 0.817 | 0.821 | 0.907 | 0.507 | 0.397 | 0.812 | 0.897 | 0.818 |
| 5 | Hierarchical_MC_MLR | 0.814 | 0.814 | 0.819 | 0.906 | 0.510 | 0.397 | 0.812 | 0.897 | 0.820 |
| 6 | KMeans_MC_MLR | 0.812 | 0.812 | 0.817 | 0.902 | 0.513 | 0.396 | 0.787 | 0.882 | 0.800 |
| 7 | KMeans_GA_MLR | 0.810 | 0.810 | 0.815 | 0.898 | 0.515 | 0.413 | 0.761 | 0.867 | 0.784 |
| 8 | SOM_MC_MLR | 0.801 | 0.801 | 0.806 | 0.899 | 0.528 | 0.406 | 0.807 | 0.894 | 0.816 |
| 9 | KMeans_GA_PLS | 0.796 | 0.796 | 0.801 | 0.890 | 0.535 | 0.420 | 0.751 | 0.861 | 0.769 |
| 10 | KMeans_MC_PLS | 0.793 | 0.794 | 0.799 | 0.894 | 0.538 | 0.416 | 0.775 | 0.875 | 0.790 |
| 11 | SOM_GA_PLS | 0.779 | 0.779 | 0.785 | 0.886 | 0.556 | 0.428 | 0.796 | 0.887 | 0.804 |
| 12 | SOM_GA_MLR | 0.779 | 0.779 | 0.785 | 0.887 | 0.556 | 0.427 | 0.796 | 0.887 | 0.805 |

External predictivity varied across descriptor-clustering and feature-selection strategies, indicating that the topology imposed on the descriptor space had a measurable influence on the generalization capacity of the resulting linear models. The Hierarchical-GA pathway produced the strongest external predictivity among the baseline linear models, with both PLS and MLR variants achieving `Q²_ext F2 = 0.820` and `CCC_ext = 0.910`. Several additional models also exceeded or closely approached this performance range, including SOM-MC-PLS (`Q²_ext F2 = 0.817`) and the Hierarchical-MC models (`Q²_ext F2 = 0.817` for PLS and `0.814` for MLR).

Importantly, multiple eight-descriptor linear baselines surpassed the historical Gramatica et al. (2014) external benchmark of `R²_ext = 0.794`. This included both Hierarchical-GA models and several MC-driven pathways, demonstrating that the feature-cluster-guided PyQSAR3 workflow can produce externally predictive, highly parsimonious linear models under the same descriptor-count constraint as the reference model. At the same time, the results are interpreted here as baseline linear-model performance rather than final model selection, because additional non-linear and ensemble modeling expansions remain planned.

Overall, the high external predictive coefficients, strong concordance values, and limited gaps between internal and external metrics indicate that the eight-descriptor linear baselines provide a robust benchmarking foundation. These results justify subsequent investigation of non-linear or ensemble modeling frameworks, such as Support Vector Regression and Random Forests, to determine whether additional non-linear structure-property relationships can be captured beyond the current linear descriptor representations.

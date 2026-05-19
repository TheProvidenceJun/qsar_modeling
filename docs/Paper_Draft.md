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

The twelve baseline linear models were next evaluated using the complete internal and external validation framework. Although external predictivity was quantified using three complementary external predictive coefficients (`Q²_ext F1`, `Q²_ext F2`, and `Q²_ext F3`), together with `CCC_ext`, `RMSE_ext`, and `MAE_ext`, model ranking at this stage was organized strictly by internal cross-validation robustness. This ordering prevents model selection from being influenced by the external test set and preserves the external set for final verification rather than optimization.

The full extended performance matrix is reported below, ranked in descending order by `Q²_cv`. In the event of near-equivalent cross-validation performance, training fit and external predictivity are reported only as descriptive supporting metrics, not as the primary selection criterion.

| Rank | Model | R²_train | CCC_tr | RMSE_train | Q²_cv | CCC_cv | RMSE_cv | MAE_cv | Q²_ext F1 | Q²_ext F2 | Q²_ext F3 | CCC_ext | RMSE_ext | MAE_ext |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Hierarchical_MC_MLR | 0.820182 | 0.901209 | 0.508648 | 0.811974 | 0.896656 | 0.520127 | 0.419998 | 0.814267 | 0.814252 | 0.819122 | 0.906043 | 0.510145 | 0.397140 |
| 2 | Hierarchical_MC_PLS | 0.818294 | 0.900068 | 0.511312 | 0.811959 | 0.896633 | 0.520149 | 0.420112 | 0.816532 | 0.816517 | 0.821328 | 0.906940 | 0.507025 | 0.396826 |
| 3 | Hierarchical_GA_PLS | 0.816350 | 0.898890 | 0.514040 | 0.807286 | 0.894252 | 0.526572 | 0.417462 | 0.819976 | 0.819961 | 0.824682 | 0.909574 | 0.502244 | 0.405073 |
| 4 | Hierarchical_GA_MLR | 0.816775 | 0.899148 | 0.513444 | 0.807124 | 0.894196 | 0.526792 | 0.418151 | 0.819748 | 0.819733 | 0.824460 | 0.909900 | 0.502562 | 0.404386 |
| 5 | SOM_MC_MLR | 0.815659 | 0.898472 | 0.515005 | 0.806598 | 0.893811 | 0.527511 | 0.410281 | 0.800784 | 0.800767 | 0.805991 | 0.899315 | 0.528338 | 0.406454 |
| 6 | SOM_MC_PLS | 0.815863 | 0.898595 | 0.514721 | 0.803669 | 0.892734 | 0.531490 | 0.414283 | 0.817486 | 0.817471 | 0.822257 | 0.908353 | 0.505706 | 0.382805 |
| 7 | SOM_GA_MLR | 0.804537 | 0.891682 | 0.530314 | 0.796103 | 0.887085 | 0.541635 | 0.427792 | 0.779310 | 0.779292 | 0.785079 | 0.886501 | 0.556084 | 0.426606 |
| 8 | SOM_GA_PLS | 0.803691 | 0.891162 | 0.531461 | 0.795522 | 0.886833 | 0.542406 | 0.428920 | 0.779388 | 0.779370 | 0.785155 | 0.886260 | 0.555986 | 0.427970 |
| 9 | KMeans_MC_MLR | 0.799899 | 0.888827 | 0.536569 | 0.787166 | 0.881709 | 0.553378 | 0.441287 | 0.811955 | 0.811940 | 0.816871 | 0.901864 | 0.513310 | 0.395886 |
| 10 | KMeans_MC_PLS | 0.789885 | 0.882610 | 0.549832 | 0.774646 | 0.874989 | 0.569421 | 0.458416 | 0.793505 | 0.793488 | 0.798903 | 0.893517 | 0.537903 | 0.415934 |
| 11 | KMeans_GA_MLR | 0.783595 | 0.878669 | 0.558001 | 0.760978 | 0.867046 | 0.586435 | 0.462205 | 0.810462 | 0.810446 | 0.815417 | 0.898439 | 0.515344 | 0.413021 |
| 12 | KMeans_GA_PLS | 0.768694 | 0.869222 | 0.576892 | 0.751181 | 0.861042 | 0.598333 | 0.468576 | 0.795683 | 0.795666 | 0.801024 | 0.890228 | 0.535059 | 0.420050 |

The reordered linear baseline matrix shows that the Hierarchical-MC track provides the strongest internal robustness among the linear models. Hierarchical-MC-MLR achieved the highest cross-validation score (`Q²_cv = 0.811974`), followed immediately by Hierarchical-MC-PLS (`Q²_cv = 0.811959`). Both models also showed high cross-validated concordance (`CCC_cv = 0.896656` and `0.896633`, respectively), indicating that the MC-guided descriptor selection within the Hierarchical feature-cluster topology produced especially stable linear baselines.

External predictivity remained strong across several models, and multiple eight-descriptor linear baselines exceeded the historical Gramatica et al. (2014) external benchmark of `R²_ext = 0.794`. However, these external results are interpreted as verification evidence rather than model-selection criteria. The distinction is essential: choosing a model because of superior external test performance would constitute data snooping and would compromise the independence of the external validation set.

Overall, the Hierarchical feature-clustering strategies dominated the upper ranks of the cross-validation-ordered linear matrix, with both Hierarchical-MC and Hierarchical-GA pathways producing high internal robustness and strong external verification metrics. These results establish a rigorous linear benchmarking foundation for comparison with the subsequent non-linear SVR and RF models, while preserving the central principle that the final champion model must be selected by `Q²_cv` rather than external test performance.

### 3.3. Non-linear Modeling Performance (SVR and RF)

Following the linear baseline matrix, non-linear modeling was performed using Support Vector Regression (SVR) and Random Forest (RF) algorithms on the Hierarchical feature-clustering track. The same descriptor-count constraint was retained: all non-linear models were restricted to exactly eight molecular descriptors. Genetic Algorithm (GA) and Monte Carlo (MC) feature-selection strategies were used to identify candidate descriptor subsets, and model hyperparameters were optimized by internal 5-fold GridSearchCV for each evaluated feature subset. For SVR, the optimized parameters included the regularization constant `C` and kernel coefficient `gamma`; for RF, the optimized parameters included the number of trees and maximum tree depth.

The four resulting non-linear models are summarized in Table 2. The Hierarchical-GA-SVR model achieved `R²_train = 0.880`, `Q²_cv = 0.811`, `CCC_cv = 0.893`, `Q²_ext F2 = 0.809`, and `CCC_ext = 0.903`, using the optimized parameter combination `C = 100` and `gamma = scale`. The Hierarchical-MC-SVR model yielded lower internal and external predictivity (`R²_train = 0.817`, `Q²_cv = 0.783`, `Q²_ext F2 = 0.764`) with `C = 100` and `gamma = 0.01`. The RF models achieved stronger apparent training fits and high external scores: Hierarchical-GA-RF reached `R²_train = 0.966`, `Q²_cv = 0.794`, `Q²_ext F2 = 0.826`, and `CCC_ext = 0.912` with `n_estimators = 100` and `max_depth = 10`, whereas Hierarchical-MC-RF reached `R²_train = 0.972`, `Q²_cv = 0.787`, `Q²_ext F2 = 0.803`, and `CCC_ext = 0.899` with `n_estimators = 100` and unrestricted tree depth.

| Model | R²_train | Q²_cv | CCC_cv | Q²_ext F2 | CCC_ext | Optimized Hyperparameters |
|---|---:|---:|---:|---:|---:|---|
| Hierarchical_GA_SVR | 0.880 | 0.811 | 0.893 | 0.809 | 0.903 | `C = 100`, `gamma = scale` |
| Hierarchical_MC_SVR | 0.817 | 0.783 | 0.880 | 0.764 | 0.871 | `C = 100`, `gamma = 0.01` |
| Hierarchical_GA_RF | 0.966 | 0.794 | 0.888 | 0.826 | 0.912 | `n_estimators = 100`, `max_depth = 10` |
| Hierarchical_MC_RF | 0.972 | 0.787 | 0.885 | 0.803 | 0.899 | `n_estimators = 100`, `max_depth = None` |

The non-linear results revealed a clear contrast between apparent external performance and internal robustness. The RF models produced slightly higher external predictive coefficients, with `Q²_ext F2` reaching 0.826 for Hierarchical-GA-RF. However, both RF models also showed large gaps between training fit and cross-validated performance, with `R²_train` values greater than 0.96 but `Q²_cv` values below 0.80. This pattern indicates a greater tendency toward overfitting relative to the SVR models. In contrast, the GA-SVR model maintained a more balanced validation profile, combining strong training fit with the highest non-linear cross-validation score (`Q²_cv = 0.811`) and high external concordance.

These findings demonstrate that the non-linear expansion provides useful complementary evidence, but they do not by themselves define the final champion model. A rigorous model-selection criterion must therefore be applied across the complete 16-model candidate set, comprising the 12 linear baseline models and the 4 non-linear models, before proceeding to Monte Carlo cross-validation and chance-correlation analysis.

### 3.4. Champion Model Selection and Monte Carlo Cross-Validation (MCCV)

To prevent data snooping and minimize the risk of overfitting-driven model selection, the final champion model was selected exclusively according to internal cross-validation robustness. External test-set performance was not used as a selection criterion. Across the complete set of 16 evaluated models, comprising 12 linear baselines and 4 non-linear models, the multiple linear regression model derived from Monte Carlo feature selection within the Hierarchical descriptor-clustering space (`Hierarchical_MC_MLR`) achieved the highest cross-validation score (`Q²_cv = 0.812`) and was therefore selected as the champion model.

The champion model retained eight molecular descriptors: `ABC`, `BCUTs-1h`, `C1SP2`, `ETA_shape_y`, `FilterItLogS`, `NdS`, `SlogP_VSA1`, and `SlogP_VSA2`. These descriptors define the locked feature set used for subsequent robustness evaluation. The selection of a highly parsimonious multiple linear regression model is notable because it indicates that strong predictive behavior was achieved without reliance on high-capacity non-linear learners or overparameterized descriptor combinations.

The stability of the locked champion model was then evaluated by 100-iteration Monte Carlo Cross-Validation (MCCV). In each iteration, the complete curated dataset was randomly partitioned into an 80% training subset and a 20% test subset. The locked `Hierarchical_MC_MLR` model was refitted on each random training subset using the same eight descriptors, and predictive performance was evaluated on the corresponding held-out subset. This repeated resampling protocol provides a stringent assessment of whether model performance remains stable across alternative train/test partitions rather than depending on a single favorable split.

Across the 100 MCCV iterations, the champion model achieved an average external coefficient of determination of `R²_ext = 0.806 ± 0.032`, an average external root mean squared error of `RMSE_ext = 0.518 ± 0.043`, and an average external mean absolute error of `MAE_ext = 0.407 ± 0.033`. The low standard deviations across the repeated random splits demonstrate that the model's predictive performance is highly stable and only weakly dependent on the specific composition of the training and test partitions.

These MCCV results confirm the structural generalizability of the champion model. The mean external performance consistently exceeded the historical Gramatica et al. (2014) external benchmark (`R²_ext = 0.794`), while the low dispersion of the MCCV metrics indicates that the observed predictive capability is robust across repeated resampling. The locked `Hierarchical_MC_MLR` model is therefore a statistically stable champion candidate for final chance-correlation analysis by Y-scrambling.

### 3.5. Chance Correlation Test (Y-Scrambling)

To exclude the possibility that the high predictive performance of the champion model resulted from chance correlation, a Y-scrambling analysis was performed. The experimental `logKoc` response values in the training set were randomly permuted 100 times, while the molecular descriptor matrix was kept unchanged. For each permuted dataset, the same multiple linear regression modeling procedure was reapplied using the locked eight-descriptor feature set of the `Hierarchical_MC_MLR` model, and the resulting scrambled coefficient of determination was recorded.

The original champion model achieved a training coefficient of determination of `R²_train = 0.820`. In contrast, the Y-scrambled models showed a near-complete collapse in explanatory performance, with an average scrambled coefficient of determination of `R²_y-sc = 0.016 ± 0.007`. The maximum scrambled `R²` observed across all 100 permutations was only 0.039. This large separation between the original model and the scrambled-response models demonstrates that the predictive behavior of the champion model cannot be reproduced when the response-descriptor relationship is destroyed.

The Y-scrambling results therefore mathematically reject the hypothesis of chance correlation. The pronounced degradation from `R²_train = 0.820` to an average `R²_y-sc = 0.016` confirms that the `Hierarchical_MC_MLR` model is based on authentic structure-property relationships rather than random statistical association. Together with the external validation and MCCV results, this establishes the final champion model as both statistically robust and chemically meaningful.

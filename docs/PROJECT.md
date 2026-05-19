# QSAR logKoc Modeling Project

## Project Overview
This project builds and validates an advanced QSAR model for predicting **logKoc** using the QDB.177 dataset from Gramatica et al. 2014.

The main objective is to develop an ensemble/consensus QSAR workflow, test a wide range of model combinations, identify the best-performing result, and compare the final external validation metrics against the original paper.

- **Target property:** logKoc
- **Reference dataset:** QDB.177, Gramatica et al. 2014
- **Reference benchmark:** Original external Test R² = 0.794
- **Primary goal:** Build the strongest possible consensus QSAR model and report the results in formal review-paper format.

## Roles
- **Parent Engineer & Documentation Manager:** Plans the workflow, manages documentation, reviews reported outputs, and gives step-by-step instructions.
- **Executor:** Runs all Python code locally in Jupyter Lab and reports outputs back for interpretation and next-step planning.

## Conda Environments
This project uses separate environments to avoid package conflicts.

### `mordred`
Used for:
- RDKit-based SMILES cleaning and standardization
- Mordred descriptor extraction (Raw calculation only)

### `pq3`
Used for:
- Data pre-processing based strictly on the Training set (Handling NaNs, Variance Thresholds, Correlation Filtering)
- PyQSAR3 feature selection and modeling
- K-Fold cross-validation & Consensus model construction
- External validation metrics & Applicability Domain analysis

## Master Curriculum

### Phase 1: Data Preparation & Split
Clean the raw SMILES using RDKit and prepare the modeling dataset.
- Load the raw QDB.177 data (Parsing XML/directory structure in `data/raw/`).
- Validate and clean SMILES with RDKit (Remove invalid structures, salts, duplicates).
- Secure approximately 643 valid compounds.
- Split the final dataset into Train/Test using Y-Ranking (80/20 split).

### Phase 2: Mordred Descriptor Extraction
Extract raw molecular descriptors in the `mordred` environment.
- Calculate Mordred descriptors for the Train and Test sets separately.
- **CRITICAL RULE 1:** Do NOT filter columns (e.g., NaNs, zero-variance) in this phase.
- **CRITICAL RULE 2 (Reference):** Strictly refer to `examples/example_mordred/example_mordred.md` for syntax and best practices before writing the notebook.

### Phase 3: Advanced Feature Selection & Modeling (Non-linear & Validation)

Perform strict preprocessing, descriptor clustering, linear baseline modeling, server-side non-linear expansion, and robustness validation under the Gramatica et al. 2014 8-descriptor constraint.

**Feature Constraint:** All GA and MC feature-selection runs in this phase MUST be constrained to select exactly **8 descriptors**, matching the descriptor count used in the Gramatica et al. 2014 benchmark model. Any exploratory run using a different descriptor count is non-production and must not be reported as the primary benchmark comparison.

**OECD Evaluation Metrics Framework:** Every generated production model must report the following metrics:
- **Internal Goodness-of-fit:** $R^2_{train}$, $CCC_{tr}$, $RMSE_{train}$
- **Internal Robustness (5-Fold CV):** $Q^2_{cv}$, $CCC_{cv}$, $RMSE_{cv}$, $MAE_{cv}$
- **External Predictivity (Test Set):** $Q^2_{ext\ F1}$, $Q^2_{ext\ F2}$, $Q^2_{ext\ F3}$, $CCC_{ext}$, $RMSE_{ext}$, $MAE_{ext}$
- **Chance Correlation:** $R^2_{y-sc}$ by Y-scrambling, performed only on the final Best Model

**Champion Model Selection Criteria:** The absolute Best Model must be selected strictly by the highest cross-validation score, $Q^2_{cv}$.
- **Prevention of Overfitting:** Simple training fit ($R^2_{train}$) can be artificially inflated by memorizing the training data. Cross-validation ($Q^2_{cv}$) tests the model on unseen subsets within the training space and therefore reveals its true mathematical capability.
- **Zero Tolerance for Data Leakage (Anti-Data Snooping):** Selecting a model based on external validation metrics ($R^2_{ext}$ or $Q^2_{ext}$) is a statistical fallacy because it allows information from the external Test set to leak into the model selection process. The Test set must remain completely unseen until final verification.
- **Objective Proof of Generalization:** A high $Q^2_{cv}$ is the most honest and robust statistical evidence that the model can handle unknown chemical structures without structural or mathematical volatility.

#### Step 3.1 - Step 3.3: Completed Data Preparation, Feature Clustering, and Linear Baseline Matrix
- Step 3.1 completed native PyQSAR3-compatible pre-filtering using training-set-only filtering decisions.
- Step 3.2 completed native PyQSAR3 descriptor/feature clustering using Hierarchical, K-Means, and SOM tracks.
- Step 3.3 completed the 12-model linear baseline matrix under the strict 8-descriptor constraint:
  - 3 descriptor-clustering tracks: Hierarchical, K-Means, SOM
  - 2 feature-selection engines: GA and MC
  - 2 linear regressors: MLR and PLS

#### Step 3.4: Server-Side Non-linear Modeling Script
Execute SVR and Random Forest models through a standalone Python script, not a Jupyter notebook.
- **Script target:** `run_nonlinear_models.py`
- **Execution environment:** `qsar_ml`
- **Parallelization:** Use multiprocessing with `max(1, total_cores - 2)` workers to preserve server responsiveness.
- **Scope:** Focus strictly on the Hierarchical feature-cluster track.
- **Feature selection:** Couple GA and MC feature-selection strategies with SVR and RF modeling, constrained to exactly 8 descriptors.
- **Hyperparameter tuning:** Use rigorous `GridSearchCV` for SVR and RF hyperparameters.
- **Required outputs:**
  - `data/features/best_nonlinear_config.json`
  - `data/features/nonlinear_model_metrics.csv`
  - `data/features/nonlinear_search_log.txt`
- **Required metric columns for `nonlinear_model_metrics.csv`:**
  - $R^2_{train}$, $CCC_{tr}$, $RMSE_{train}$
  - $Q^2_{cv}$, $CCC_{cv}$, $RMSE_{cv}$, $MAE_{cv}$
  - $Q^2_{ext\ F1}$, $Q^2_{ext\ F2}$, $Q^2_{ext\ F3}$
  - $CCC_{ext}$, $RMSE_{ext}$, $MAE_{ext}$

#### Step 3.5: MCCV (Monte Carlo Cross-Validation)
Select the absolute Best Model among the 16 evaluated production candidates strictly according to the highest $Q^2_{cv}$:
- 12 linear baseline models from Step 3.3
- 4 non-linear models from Step 3.4: GA-SVR, GA-RF, MC-SVR, MC-RF

External validation metrics must not be used to choose the Best Model. After $Q^2_{cv}$-based selection, lock the Best Model's 8 descriptors and optimal hyperparameters. Perform 100 random Train/Test splits to quantify robustness across repeated resampling. Report the full internal and external metric suite for the MCCV distribution, including mean, standard deviation, and relevant confidence intervals where appropriate.

#### Step 3.6: Y-Scrambling
Perform chance-correlation analysis on the final Best Model only.
- Randomly permute the training response vector across repeated scrambling trials.
- Refit the locked modeling pipeline under each scrambled response condition.
- Report $R^2_{y-sc}$ and compare scrambled performance against the non-scrambled Best Model to demonstrate that the final model is not attributable to chance correlation.

#### CRITICAL RULE (Reference)
Strictly refer to `examples/example_pyqsar3/example_pyqsar3.md` to map the exact native syntax for PyQSAR3 preprocessing, clustering integration, and GA/MC engine invocation before writing production notebook cells or scripts that consume PyQSAR3 artifacts.


### Phase 4: Applicability Domain & Visualization
Generate publication-quality visual diagnostics for the final Best Model.

#### Step 4.1: Predicted vs. Experimental Plot
- Generate publication-quality scatter plots comparing experimental and predicted logKoc values for the final Best Model.
- Plot Training and Test predictions with clear visual separation.
- Include the identity line, fitted trend if appropriate, and the final model's key validation metrics.
- Export figures in publication-ready formats suitable for manuscript and repository documentation.

#### Step 4.2: Williams Plot (Applicability Domain)
- Calculate standardized residuals for the final Best Model.
- Calculate hat leverage values and the warning leverage threshold $h^*$.
- Generate a Williams Plot to identify structural outliers and response outliers.
- Interpret the Applicability Domain in relation to model reliability and external predictivity.

### Step 4.3: Comprehensive QMRF Metric Integration & Benchmark (NEW)

Synthesize and align all multi-track validation statistics into a unified, publication-ready master matrix to satisfy rigorous OECD QMRF specifications.

1. **Dual-Architecture Synthesis (Linear vs. Non-Linear):**
   - Import the 12 linear baseline results (MLR/PLS from Step 3.3) and the 4 non-linear advanced estimator results (SVR/RF from Step 3.4).
   - Construct a symmetrical master benchmark table comparing all 16 models across 3 core dimensions: Goodness-of-Fit, Robustness (LMO), and Predictivity.

2. **Champion Lock-in Statistics Extraction:**
   - For the absolute winning combination (Hierarchical_MC_MLR), extract and append the exact 1,000-run Monte Carlo Cross-Validation (MCCV) statistical distributions denoted as Mean ± Standard Deviation ($SD$).
   - Compute and log the final Y-Scrambling coefficient ($R^2_{y-sc}$) to confirm the model's structural resistance against chance correlation.

3. **OECD QMRF Unified Metric Compliance Matrix:**
   - Ensure the integrated DataFrame strictly computes and displays the following columns for cross-architecture verification:
     - **Goodness-of-Fit:** $R^2_{train}$, $RMSE_{train}$, $CCC_{tr}$
     - **Robustness (Internal CV):** $Q^2_{cv}$, $RMSE_{cv}$, $MAE_{cv}$, $CCC_{cv}$
     - **Validation/Predictivity:** $Q^2_{ext\ F1}$, $Q^2_{ext\ F2}$, $Q^2_{ext\ F3}$, $RMSE_{ext}$, $MAE_{ext}$, $CCC_{ext}$
     - **Reliability Safeguards:** $R^2_{y-sc}$ (Champion exclusive), MCCV Stability ($R^2_{ext} \pm SD$)

### Phase 5: Paper Draft & Repository Documentation
Prepare formal documentation and manuscript-style outputs.
- Draft `docs/Paper_Draft.md` (Intro to Supp Info).
- Compare final model performance against the original paper.
- Document the workflow in `README.md`.

## 📂 Directory Architecture
```text
qsar_modeling/
├── data/                       # 📊 데이터 저장소 (Git push 제외 권장)
│   ├── raw/                    # QDB.177 원본 데이터 (SMILES, logKoc)
│   ├── processed/              # Y-Ranking으로 분할된 Train/Test 데이터
│   └── features/               # Mordred 추출 Descriptor CSV
├── envs/                       # 🌐 Conda 환경 설정 (yml 첨부 위치)
├── notebooks/                  # 📓 Jupyter Lab 작업 공간
│   ├── 01_data_split.ipynb           
│   ├── 02_mordred_extract.ipynb      
│   ├── 03_pyqsar3_modeling.ipynb     
│   └── 04_result_metrics.ipynb       
├── docs/                       # 📝 논문 작성 및 LLM 지침 문서
│   ├── literature/             # 참고 문헌
│   ├── Paper_Draft.md          # 1~6장 최종 논문 초안
│   ├── README.md               # GitHub 메인 설명서
│   ├── PROJECT.md              # [현재 파일] 프로젝트 방향성
│   ├── HANDOFF.md              # 세션 전환용 상태 요약본
│   └── AGENTS.md               # LLM 라우팅 지침
├── examples/                   # 💡 참조용 예시 코드
└── .gitignore                  # Git 제외 목록

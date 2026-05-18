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

### Phase 3: Feature Filtering & Multi-Track PyQSAR3 Modeling

Perform strict descriptor clustering, multi-track feature selection, and train linear statistical models natively within the py3 environment.

#### 1. Feature Filtering & Common Preprocessing
- **Basic Filters:** Apply strict preprocessing filters (remove NaNs, zero-variance, high correlation descriptors) strictly using the Training set only. Transform the Test set accordingly to prevent data leakage.
- **3-Track Clustering Evaluation:** Apply three distinct clustering methods to the filtered training features to handle multicollinearity:
  1) Hierarchical Clustering
  2) K-Means Clustering
  3) SOM (Self-Organizing Maps) *[Note: Use native pyqsar3 or sklearn approximation if SOM is constrained]*
- **Silhouette Scoring:** Compute the Silhouette Score for all three clustering tracks to evaluate group structural quality.

#### 2. Multi-Track PyQSAR3 Execution Matrix (12 Models)
Feed each of the 3 cluster datasets into the coupled pyqsar3 selection-modeling engines. Execute both Feature Selection and Model Building simultaneously for all 12 combinations:
- **Hierarchical Track (4 Models):**
  - Cluster-GA ➔ MLR, PLS
  - Cluster-MC ➔ MLR, PLS
- **K-Means Track (4 Models):**
  - Cluster-GA ➔ MLR, PLS
  - Cluster-MC ➔ MLR, PLS
- **SOM Track (4 Models):**
  - Cluster-GA ➔ MLR, PLS
  - Cluster-MC ➔ MLR, PLS

#### 3. Validation & Model Selection
- **Internal Cross-Validation:** Evaluate all 12 coupled pipelines using 5-Fold Cross-Validation ($Q^2$ or $R^2_{cv}$) exclusively on the training folds.
- **External Evaluation:** Validate predictions on the split Test set ($R^2_{ext}$, RMSE).
- **Winning Track Selection:** Compare metrics to find the Absolute Winner among the 12 models.
- **Optional Advanced ML Expansion:** (Optional) Extract the optimal feature indices from the winning clustering track, and feed them into external SVR and Random Forest models with hyperparameter tuning to benchmark linear vs. non-linear performance.

#### CRITICAL RULE (Reference)
Strictly refer to `examples/example_pyqsar3/example_pyqsar3.md` to map the exact native syntax for clustering integration and GA/MC engine invocation before writing the production notebook cells.


### Phase 4: External Validation & Applicability Domain
Evaluate the selected model on the external test set.
- Calculate external validation metrics (R²_ext, RMSE, MAE, CCC).
- Compare performance against the original Gramatica 2014 benchmark.
- Draw a Williams Plot and assess Applicability Domain.

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

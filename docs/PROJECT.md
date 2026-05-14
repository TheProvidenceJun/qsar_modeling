# 🚀 QSAR logKoc Review Project

## 📌 Project Overview
본 프로젝트는 Gramatica et al. (2014)의 환경 오염 물질 QSAR 모델링 연구를 기반으로, 최신 앙상블/컨센서스(Consensus) 머신러닝 기법을 적용하여 모델의 예측력을 극대화하고 이를 엄격하게 검증하는 리뷰 논문 작성 프로젝트입니다.

- **Target Property:** logKoc (Soil sorption coefficient)
- **Dataset:** 원본 데이터 정제를 통해 유효 데이터셋을 확보 
- **Primary Goal:** 다양한 알고리즘 조합을 테스트하여 "Best Result"를 도출하고, 원본 논문의 성능($R^2_{ext} = 0.794$, $\sigma = 0.543$) 및 세계적 QSAR 신뢰도 검증 평가지표(Metric)와 비교 분석합니다.

## 📝 Review Paper Structure
최종 결과물인 `Paper_Draft.md`는 다음의 학술 논문 형식을 엄격히 따릅니다:
1. **Introduction:** logKoc의 중요성 및 기존 연구 한계, 앙상블 기법 도입의 필요성
2. **Methods:** Y-Ranking 데이터 분할, Mordred 디스크립터 추출, PyQSAR3 앙상블 모델링 및 교차 검증 방법론
3. **Results and discussion:** - 자체 추가된 엄격한 Metric (예: CCC, Q^2_F1/F2/F3 등) 결과
    - Original 결과(Gramatica 2014)와의 직접 비교
    - Williams Plot을 통한 Applicability Domain(AD) 검증
4. **Conclusion:** Best 조합의 우수성 및 환경 독성 예측 모델로서의 가치 요약
5. **References:** 최대 10 pages 이내로 마무리
6. **Supplementary Information:** 평가지표(Metric) 산출에 사용된 핵심 Python 코드 첨부

## 🛠 Tech Stack & Conda Environments
환경 간의 패키지 충돌(Data Leakage 및 Dependency 이슈)을 방지하기 위해 두 개의 독립된 Conda 가상환경을 사용합니다.
- **`conda activate mordred`**: 분자 디스크립터 추출 및 데이터 전처리 (`02_mordred_extract.ipynb`)
- **`conda activate py3`**: PyQSAR3 기반 앙상블 모델 학습 및 Metric 평가 (`03_pyqsar3_modeling.ipynb`, `04_result_metrics.ipynb`)

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

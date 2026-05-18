# 🌍 QSAR Modeling for Soil Sorption Coefficient (logKoc)

![Status](https://img.shields.io/badge/Status-In%20Progress-orange)
![Python](https://img.shields.io/badge/Python-3.x-blue)
![Machine Learning](https://img.shields.io/badge/Task-QSAR%20Modeling-success)

## 📌 Project Overview
본 프로젝트는 환경 오염 물질의 토양 흡착 계수(Soil sorption coefficient, logKoc)를 예측하기 위한 **QSAR(Quantitative Structure-Activity Relationship)** 모델링 파이프라인입니다. 
Gramatica et al. (2014)의 연구(QDB.177 데이터셋)를 기반으로 하여, 최신 **앙상블/컨센서스(Ensemble/Consensus) 머신러닝 기법**을 적용해 예측력을 개선하고 엄격한 외부 검증을 수행하는 것을 목표로 합니다.

- **Target Endpoint:** logKoc
- **Original Benchmark:** $R^2_{ext} = 0.794$, $\sigma = 0.543$ (Gramatica 2014)
- **Methodology:** Y-Ranking Data Split, Mordred Descriptor Extraction, PyQSAR3 Ensemble Modeling, Williams Plot (Applicability Domain)

## 🚀 Workflow & Pipeline
1. **Data Preprocessing:** `RDKit`을 이용한 SMILES 표준화 및 불순물/중복 데이터 제거.
2. **Rational Data Splitting:** 화학적 공간(Chemical Space)을 보존하기 위한 **Y-Ranking (Response Sorting)** 기반 Train(80%)/Test(20%) 분할.
3. **Feature Extraction:** `Mordred` 라이브러리를 활용한 2D/3D 분자 표현자(Descriptor) 추출 및 정제.
4. **Ensemble Modeling:** `PyQSAR3`를 이용한 K-Fold 교차 검증 및 최적의 앙상블 조합(Best Combination) 발굴.
5. **Validation & Applicability Domain:** CCC(Concordance Correlation Coefficient) 산출 및 Williams Plot 시각화.

## ⚙️ Environments & Installation
의존성 충돌을 방지하기 위해 데이터 전처리와 모델링 환경을 엄격히 분리하여 운영합니다.

```bash
# 1. 분자 특성 추출 환경 (Mordred)
conda env create -f envs/env_mordred.yml
conda activate mordred

# 2. QSAR 앙상블 모델링 환경 (PyQSAR3)
conda env create -f envs/env_pyqsar3.yml
conda activate pq3

```

## 📂 Directory Structure

```text
qsar_modeling/
├── data/           # 정제된 데이터 및 Descriptor Feature 보관
├── envs/           # Conda 가상환경 설정 파일 (yml)
├── notebooks/      # 단계별 Jupyter Notebook 실행 파일 (01~04)
├── docs/           # 논문 초안 (Paper_Draft.md) 및 상태 관리 문서
└── examples/       # 코드 작성 참조용 예시 스크립트

```

## 📈 Final Results & Conclusion

> 🚧 **진행 중 (Work In Progress)**
> *모델링이 완료되면 이 섹션에 우리 모델의 최종 평가 지표(Metric)와 원본 논문과의 비교 결과를 업데이트할 예정입니다.*

## 📚 References

* Gramatica, P.; Cassani, S.; Chirico, N. QSARINS-chem: Insubria datasets and new QSAR/QSPR models for environmental pollutants in QSARINS. *J. Comput. Chem.* **2014**, 35, 1036–1044. (QDB archive DOI: 10.15152/QDB.177)


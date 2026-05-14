# 🔄 Project Handoff Status

**Current Phase:** Phase 1 - Data Preprocessing & Cleaning
**Active Environment:** `mordred` (RDKit 사용을 위해 추천)
**Status:** 프로젝트 설계 완료 및 원본 데이터 로드 준비 중.

## ✅ Completed Tasks
- [x] 프로젝트 폴더 구조 및 가상환경 전략 수립.
- [x] `PROJECT.md`에 데이터 정제 단계 반영 및 업데이트 완료.

## 🚀 Next Mission (For Next Session)
- **Task: 01_data_prep_and_split.ipynb 실행**
- **목표:** `data/raw/`에 있는 가공되지 않은 데이터를 정제하여 **약 643개의 유효 화합물**을 추출하고, 이를 Y-Ranking 방식으로 분할함.
- **상세 작업:** 1. SMILES 표준화 및 중복/염 제거.
   		 2. 최종 정제된 리스트를 `train.csv`와 `test.csv`로 저장.
- **주의사항:** RDKit 패키지를 사용하여 SMILES의 화학적 무결성을 확인해야 함.

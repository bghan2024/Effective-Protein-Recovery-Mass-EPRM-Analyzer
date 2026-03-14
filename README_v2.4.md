# EPRM Analyzer v2.4
## Effective Protein Recovery Mass Analyzer

단백질 서열의 물리화학적 특성을 분석하여 정제 후 예상 회수 농도를 예측하는 도구입니다.

---

## 📋 목차 (Table of Contents)

1. [소개](#소개)
2. [주요 기능](#주요-기능)
3. [설치 방법](#설치-방법)
4. [빠른 시작](#빠른-시작)
5. [사용 가이드](#사용-가이드)
6. [파라미터 설명](#파라미터-설명)
7. [출력 결과](#출력-결과)
8. [인용 방법](#인용-방법)
9. [라이선스](#라이선스)
10. [문제 해결](#문제-해결)
11. [참고문헌](#참고문헌)

---

## 🎯 소개

EPRM Analyzer는 단백질 정제 실험에서 예상되는 회수 농도를 예측하는 도구입니다. 
단백질 서열만 입력하면 다음과 같은 정보를 제공합니다:

- **예상 회수 농도**: 정제 후 예상되는 유효 단백질 농도
- **물리화학적 특성**: 분자량, 불안정성 지수, GRAVY, 등전점
- **불확실성 정량화**: Monte Carlo 시뮬레이션을 통한 신뢰구간
- **실험 가이드**: 목표 농도(20nM) 달성을 위한 희석 배수

### 왜 EPRM Analyzer를 사용하나요?

- ✅ **시간 절약**: 실험 전에 예상 농도를 미리 알 수 있어 실험 계획이 쉬워집니다
- ✅ **비용 절감**: 불필요한 실험을 줄일 수 있습니다
- ✅ **과학적 근거**: 10개 이상의 주요 논문에 기반한 알고리즘
- ✅ **사용하기 쉬움**: FASTA 파일만 있으면 자동으로 분석합니다

---

## ✨ 주요 기능

### 1. 물리화학적 특성 분석
- 분자량 (Molecular Weight)
- 불안정성 지수 (Instability Index)
- GRAVY (소수성 지수)
- 등전점 (Isoelectric Point, pI)

### 2. 회수율 예측
- 키트 효율 반영
- 시스템 손실 반영 (피펫팅, Dead Volume, 표면 흡착)
- 단백질 특성 기반 손실 반영

### 3. 불확실성 정량화
- Monte Carlo 시뮬레이션 (1000회 반복)
- 평균, 표준편차, 95% 신뢰구간 제공

### 4. 일괄 처리
- 여러 FASTA 파일 자동 처리
- 프로젝트 파일 자동 제외
- 상세한 로그 및 JSON 결과 저장

---

## 📦 설치 방법

### 필수 요구사항

- Python 3.7 이상
- pip (Python 패키지 관리자)

### 1단계: 필수 라이브러리 설치

터미널(Windows: PowerShell 또는 CMD, Mac/Linux: Terminal)에서 다음 명령어를 실행하세요:

```bash
# 필수 라이브러리
pip install biopython numpy

# 선택적 라이브러리 (YAML 설정 파일 사용 시)
pip install pyyaml
```

### 2단계: 파일 다운로드

GitHub에서 `eprm_analyzer_v2.4.py` 파일을 다운로드하거나 클론하세요:

```bash
git clone https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer.git
cd EPRM-Analyzer
```

### 3단계: 설치 확인

Python에서 import가 되는지 확인하세요:

```python
from eprm_analyzer_v2_4 import EPRMAnalyzer
print("설치 완료!")
```

---

## 🚀 빠른 시작

### 예시 1: 단일 서열 분석

```python
from eprm_analyzer_v2_4 import EPRMAnalyzer

# 분석기 생성 (기본 설정)
analyzer = EPRMAnalyzer()

# 단백질 서열 분석
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDPSIHAGHSVEVLELKP"
result = analyzer.calculate_eprm(sequence)

# 결과 확인
print(f"분자량: {result['MW_kDa']:.2f} kDa")
print(f"등전점: {result['pI']:.2f}")
print(f"예상 농도: {result['C_Effective_uM'][0]:.4f} ± {result['C_Effective_uM'][1]:.4f} uM")
```

### 예시 2: FASTA 파일 일괄 처리

```python
from eprm_analyzer_v2_4 import EPRMAnalyzer

# 분석기 생성
analyzer = EPRMAnalyzer()

# 현재 디렉토리의 모든 .fasta 파일 처리
analyzer.process_files()

# 결과는 EPRM_Results_YYYYMMDD_HHMMSS 폴더에 저장됩니다
```

### 예시 3: 사용자 정의 파라미터

```python
from eprm_analyzer_v2_4 import EPRMAnalyzer

# 실험 조건에 맞게 파라미터 설정
analyzer = EPRMAnalyzer(
    initial_conc_um=20.0,      # 초기 농도: 20 uM
    initial_vol_ul=100.0,      # 초기 부피: 100 uL
    final_vol_ul=500.0,        # 최종 부피: 500 uL
    buffer_ph=7.0,             # 버퍼 pH: 7.0
    random_seed=42             # 재현성 보장
)

# 분석 실행
result = analyzer.calculate_eprm("MKTAYIAKQR")
```

---

## 📖 사용 가이드

### FASTA 파일 형식

EPRM Analyzer는 표준 FASTA 형식을 지원합니다:

```
>Protein_1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDPSIHAGHSVEVLELKP

>Protein_2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDPSIHAGHSVEVLELKP
```

**주의사항**:
- 헤더는 `>`로 시작해야 합니다
- 서열은 표준 아미노산 20가지만 포함해야 합니다 (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
- 서열은 대소문자 구분 없이 입력 가능합니다

### 설정 파일 사용

JSON 또는 YAML 형식의 설정 파일을 사용할 수 있습니다:

**config.json 예시**:
```json
{
  "initial_conc_um": 10.0,
  "initial_vol_ul": 90.0,
  "final_vol_ul": 450.0,
  "eta_kit": 0.50,
  "systemic_efficiency": 0.75,
  "buffer_ph": 7.4,
  "instability_threshold": 40.0,
  "instability_penalty_factor": 80.0,
  "gravy_penalty_factor": 0.15
}
```

**사용 방법**:
```python
analyzer = EPRMAnalyzer(config_file="config.json")
```

---

## ⚙️ 파라미터 설명

### 실험 조건 파라미터

| 파라미터 | 설명 | 기본값 | 단위 | 변경 가능 |
|---------|------|--------|------|----------|
| `initial_conc_um` | 초기 단백질 농도 | 10.0 | uM | ✅ |
| `initial_vol_ul` | 초기 부피 | 90.0 | uL | ✅ |
| `final_vol_ul` | 최종 부피 | 450.0 | uL | ✅ |
| `buffer_ph` | 버퍼 pH | 7.4 | - | ✅ |

### 효율 파라미터

| 파라미터 | 설명 | 기본값 | 범위 | 변경 가능 |
|---------|------|--------|------|----------|
| `eta_kit` | 키트 효율 | 0.50 | 0.0 ~ 1.0 | ✅ |
| `systemic_efficiency` | 시스템 효율 | 0.75 | 0.0 ~ 1.0 | ✅ |

**시스템 효율 설명**:
- 0.70 ~ 0.75: 일반적인 실험 환경
- 0.75 ~ 0.80: 깨끗하고 정밀한 실험 환경
- 0.80 이상: 매우 이상적인 조건 (드뭄)

### 알고리즘 파라미터

| 파라미터 | 설명 | 기본값 | 변경 가능 |
|---------|------|--------|----------|
| `instability_threshold` | 불안정성 임계값 | 40.0 | ⚠️ 권장 안 함 |
| `instability_penalty_factor` | 불안정성 페널티 계수 | 80.0 | ⚠️ 권장 안 함 |
| `gravy_penalty_factor` | GRAVY 페널티 계수 | 0.15 | ⚠️ 권장 안 함 |

**주의**: 알고리즘 파라미터는 과학적 근거에 기반하여 설정되었으므로, 
특별한 이유가 없으면 기본값을 사용하는 것을 권장합니다.

---

## 📊 출력 결과

### 결과 파일 구조

분석을 실행하면 `EPRM_Results_YYYYMMDD_HHMMSS` 폴더가 생성되고, 
다음 파일들이 저장됩니다:

```
EPRM_Results_20261201_143022/
├── results.json              # 모든 분석 결과 (JSON 형식)
├── eprm_analysis_detail.log  # 상세 로그
└── config.json               # 사용된 설정 파라미터
```

### results.json 구조

```json
[
  {
    "file": "protein.fasta",
    "header": ">Protein_1",
    "sequence": "MKTAYIAKQR...",
    "results": {
      "MW_kDa": 25.3,
      "Instability": 35.2,
      "GRAVY": -0.15,
      "pI": 7.8,
      "Stab_Factor": 1.0,
      "Ads_Factor": 0.98,
      "pI_Factor": 0.99,
      "Eta_Prot": 0.97,
      "Total_Coeff": 0.36,
      "C_Theo_Max_uM": 2.0,
      "C_Effective_uM": [0.72, 0.05],
      "C_Effective_CI_95": [0.63, 0.82]
    }
  }
]
```

### 결과 값 설명

- **MW_kDa**: 분자량 (단위: kDa)
- **Instability**: 불안정성 지수 (40 이상이면 불안정)
- **GRAVY**: 소수성 지수 (음수=친수성, 양수=소수성)
- **pI**: 등전점
- **C_Effective_uM**: 예상 유효 농도
  - `[평균, 표준편차]` 형식 (불확실성 포함 시)
  - 단일 값 (불확실성 제외 시)
- **C_Effective_CI_95**: 95% 신뢰구간 `[하한값, 상한값]`

---

## 📝 인용 방법

### ⚠️ 중요: 이 도구를 사용할 경우 반드시 인용해야 합니다

### 출판 전 (Pre-publication)

논문이나 보고서에 다음을 포함하세요:

**Methods 섹션**:
```
Protein recovery prediction was performed using EPRM Analyzer v2.4 
(Han, 2026; https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer).
```

**참고문헌**:
```
Han, B. (2026). EPRM Analyzer v2.4. GitHub repository. 
https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer
```

### 출판 후 (Post-publication)

논문이 출판되면 논문을 인용하세요 (논문 정보는 추후 업데이트 예정):

**Methods 섹션**:
```
Protein recovery prediction was performed using EPRM Analyzer 
(Han, 2026; https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer).
```

**참고문헌**:
```
Han, B. (2026). EPRM Analyzer: A tool for predicting effective 
protein recovery mass. [Journal Name], [Volume], [Pages]. 
https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer
```

### BibTeX 형식

```bibtex
@software{eprm_analyzer_2026,
  author = {Han, Byeong-gu},
  title = {EPRM Analyzer v2.4},
  year = {2026},
  url = {https://github.com/bghan2026/Effective-Protein-Recovery-Mass-EPRM-Analyzer},
  version = {2.4.0}
}
```

---

## 📜 라이선스

이 프로그램은 **GNU Affero General Public License v3.0 (AGPL-3.0)** 
조건에 따라 배포됩니다.

### 라이선스 요약

- ✅ **자유로운 사용**: 개인, 학술, 상업적 목적으로 자유롭게 사용 가능
- ✅ **수정 가능**: 코드를 수정하여 사용 가능
- ✅ **재배포 가능**: 수정본도 재배포 가능
- ⚠️ **조건**: 
  - 원작자 표시 필수
  - 수정본도 동일한 AGPL-3.0 라이선스 적용
  - 웹 서비스로 제공 시 소스 코드 공개 필수

### 전체 라이선스 텍스트

전체 라이선스 텍스트는 다음에서 확인할 수 있습니다:
https://www.gnu.org/licenses/agpl-3.0.html

### 개발자 권익 보호

이 라이선스는 개발자의 권익을 최대한 보장합니다:
- 연구실을 떠나더라도 저작권 유지
- 코드 사용 시 개발자 인용 필수
- 상업적 사용 시에도 라이선스 조건 준수

---

## 🔧 문제 해결

### 자주 묻는 질문 (FAQ)

#### Q1: "Biopython 라이브러리가 설치되지 않았습니다" 오류

**해결 방법**:
```bash
pip install biopython
```

#### Q2: "유효하지 않은 서열" 오류

**원인**: 서열에 비표준 아미노산이 포함되어 있음

**해결 방법**:
- 서열에 표준 아미노산 20가지만 포함되어 있는지 확인
- 공백, 숫자, 특수문자 제거
- 서열이 최소 2개 아미노산 이상인지 확인

#### Q3: 결과가 예상과 다릅니다

**확인 사항**:
1. 파라미터 값이 실험 조건과 일치하는지 확인
2. 버퍼 pH 값이 올바른지 확인
3. 키트 효율 값이 키트 제조사 스펙과 일치하는지 확인

#### Q4: 분석이 너무 느립니다

**해결 방법**:
```python
# 불확실성 계산 제외 (빠른 분석)
result = analyzer.calculate_eprm(sequence, include_uncertainty=False)

# 또는 시뮬레이션 반복 횟수 줄이기
result = analyzer.calculate_eprm(sequence, n_iterations=500)
```

#### Q5: FASTA 파일이 인식되지 않습니다

**확인 사항**:
1. 파일 확장자가 `.fasta` 또는 `.txt`인지 확인
2. FASTA 형식이 올바른지 확인 (헤더는 `>`로 시작)
3. 프로젝트 파일(README.md, requirements.txt 등)은 자동으로 제외됩니다

### 에러 메시지 해석

| 에러 메시지 | 의미 | 해결 방법 |
|------------|------|----------|
| "초기 농도는 양수여야 합니다" | 농도 값이 0 이하 | 양수 값 입력 |
| "키트 회수율은 0과 1 사이여야 합니다" | 효율 값이 범위 초과 | 0.0 ~ 1.0 사이 값 입력 |
| "설정 파일을 찾을 수 없습니다" | 설정 파일 경로 오류 | 파일 경로 확인 |

---

## 📚 참고문헌

이 도구는 다음 과학적 논문들에 기반하여 개발되었습니다:

1. **Instability Index**: 
   Guruprasad, K., et al. (1990). Protein Engineering, 3(2), 155-161.

2. **GRAVY**: 
   Kyte, J., & Doolittle, R. F. (1982). Journal of Molecular Biology, 157(1), 105-132.

3. **Protein Adsorption**: 
   Norde, W. (1986). Advances in Colloid and Interface Science, 25(4), 267-340.

4. **Protein Purification**: 
   Janson, J. C. (Ed.). (2011). Protein Purification (3rd ed.). John Wiley & Sons.

5. **pI and Solubility**: 
   Gromiha, M. M., & Selvaraj, S. (2004). Progress in Biophysics and Molecular Biology, 86(2), 235-277.

6. **pI Calculation**: 
   Bjellqvist, B., et al. (1993). Electrophoresis, 14(1), 1023-1031.

7. **Solubility at pI**: 
   Shaw, K. L., et al. (2001). Protein Science, 10(6), 1206-1215.

8. **Experimental Loss Factors**: 
   Scopes, R. K. (2010). Protein Purification (3rd ed.). Springer.

9. **Monte Carlo Methods**: 
   Rubinstein, R. Y., & Kroese, D. P. (2016). Simulation and the Monte Carlo Method (3rd ed.). John Wiley & Sons.

전체 참고문헌 목록은 코드 파일의 헤더 주석을 참고하세요.

---

## 👤 작성자 및 연락처

**Han Byeong-gu**
- Email: hanbyeonggu@gmail.com
- GitHub: [@bghan2024](https://github.com/bghan2024)

### 문의 사항

- **버그 리포트**: GitHub Issues 사용
- **기능 제안**: GitHub Issues 사용
- **학술적 협업**: 이메일로 연락
- **인용 관련 문의**: 이메일로 연락

---

## 📈 버전 히스토리

### v2.4.0 (2026)
- 초보자 친화적인 문서화
- 상세한 인용 가이드라인
- AGPL-3.0 라이선스 적용
- 사용자 변경 가능 파라미터 명확히 표시

### v2.3.0
- Systemic Loss Factor 추가
- pI-pH Interaction 추가
- Instability Threshold 업데이트 (40.0)
- 참고문헌 보강

### v2.1.0
- FASTA 파일 형식 검증 강화
- 프로젝트 파일 자동 제외
- 불확실성 정량화 추가

---

## 🙏 감사의 말

이 도구를 사용해주셔서 감사합니다. 
버그 리포트, 기능 제안, 사용 후기는 언제든 환영합니다!

---

**마지막 업데이트**: 2026년  
**버전**: 2.4.0  
**라이선스**: GNU Affero General Public License v3.0

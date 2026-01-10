"""
Effective Protein Recovery Mass (EPRM) Analyzer v2.4

단백질 서열의 물리화학적 특성을 분석하여 정제 후 예상 회수 농도를 예측하는 도구입니다.
실험적 데이터를 바탕으로 시스템 손실과 등전점(pI) 효과를 보정하여 실제 실험과 유사한 
회수 계수를 산출하도록 개선된 버전입니다.

================================================================================
⚠️ 중요: 인용 및 라이선스 정보 ⚠️
================================================================================

이 코드와 알고리즘을 사용할 경우, 반드시 다음을 준수해야 합니다:

1. 출판 전 (Pre-publication):
   - GitHub 저장소를 인용하세요:
     "EPRM Analyzer was obtained from https://github.com/bghan2024/Effective-Protein-Recovery-Mass-EPRM-Analyzer"
   - 또는 코드 내 주석을 인용:
     "Analysis was performed using EPRM Analyzer v2.4 
      (https://github.com/bghan2024/Effective-Protein-Recovery-Mass-EPRM-Analyzer)"

2. 출판 후 (Post-publication):
   - 논문이 출판되면 논문을 인용하세요 (논문 정보는 추후 업데이트 예정)
   - 선택적으로 GitHub 저장소도 함께 인용 가능:
     "Analysis was performed using EPRM Analyzer (Han, 2024; 
      https://github.com/bghan2024/Effective-Protein-Recovery-Mass-EPRM-Analyzer)"

3. 실험 방법 섹션:
   - Methods 섹션에 다음을 포함하세요:
     "Protein recovery prediction was performed using EPRM Analyzer v2.4 
      (Han, 2024), which predicts effective protein recovery concentration 
      based on physicochemical properties including instability index, 
      GRAVY score, and isoelectric point."

4. 참고문헌:
   - 이 도구를 사용한 모든 논문/보고서에 반드시 인용을 포함해야 합니다.
   - 인용 없이 사용하는 것은 학술 윤리에 위배됩니다.

================================================================================
📜 라이선스: GNU Affero General Public License v3.0 (AGPL-3.0)
================================================================================

이 프로그램은 자유 소프트웨어입니다. GNU Affero General Public License v3.0 
조건에 따라 재배포하거나 수정할 수 있습니다.

이 라이선스는 개발자(작성자)의 권익을 최대한 보장합니다:

1. 저작권 보호:
   - 원작자의 저작권이 명확히 보호됩니다.
   - 코드 사용 시 원작자 표시가 필수입니다.

2. Copyleft 조항:
   - 이 코드를 수정하거나 파생작을 만들 경우, 동일한 AGPL-3.0 라이선스를 
     적용해야 합니다.
   - 웹 서비스로 제공하는 경우에도 소스 코드를 공개해야 합니다.

3. 연구실 이탈 시 보호:
   - 연구실을 떠나더라도 개발자의 저작권은 유지됩니다.
   - 코드 사용 시 개발자 인용이 필수입니다.
   - 상업적 사용 시에도 라이선스 조건을 준수해야 합니다.

4. 상업적 사용:
   - 상업적 사용이 가능하지만, AGPL-3.0 조건을 준수해야 합니다.
   - 수정된 코드도 공개해야 합니다.

전체 라이선스 텍스트는 다음에서 확인할 수 있습니다:
https://www.gnu.org/licenses/agpl-3.0.html

================================================================================
👤 작성자
================================================================================

**Han Byeong-gu**
- Email: hanbyeonggu@gmail.com
- GitHub: [@bghan2024](https://github.com/bghan2024)

개발자 연락처:
- 코드 관련 문의, 버그 리포트, 기능 제안은 GitHub Issues를 통해 
  연락해주세요.
- 학술적 협업이나 인용 관련 문의는 이메일로 연락해주세요.

================================================================================
📚 참고문헌 (References)
================================================================================

이 도구는 다음 과학적 논문들에 기반하여 개발되었습니다. 
이 도구를 사용할 때는 다음 참고문헌도 함께 인용하는 것을 권장합니다:

1. Instability Index: 
   Guruprasad, K., Reddy, B. V., & Pandit, M. W. (1990). 
   Correlation between stability of a protein and its dipeptide composition: 
   a novel approach for predicting in vivo stability of a protein from its primary sequence.
   Protein Engineering, 3(2), 155-161.
   Note: Instability Index > 40 indicates unstable proteins (Guruprasad et al., 1990).

2. GRAVY (Grand Average of Hydropathy): 
   Kyte, J., & Doolittle, R. F. (1982). 
   A simple method for displaying the hydropathic character of a protein.
   Journal of Molecular Biology, 157(1), 105-132.

3. Protein Adsorption: 
   Norde, W. (1986). 
   Adsorption of proteins from solution at the solid-liquid interface.
   Advances in Colloid and Interface Science, 25(4), 267-340.
   Note: Hydrophobic proteins show increased adsorption to solid surfaces.

4. Protein Purification Recovery: 
   Janson, J. C. (Ed.). (2011). 
   Protein Purification: Principles, High Resolution Methods, and Applications (3rd ed.).
   John Wiley & Sons.
   Note: "Losses occur at every step due to unspecific adsorption and aggregate removal."

5. Isoelectric Point (pI) and Solubility: 
   Gromiha, M. M., & Selvaraj, S. (2004). 
   Inter-residue interactions in protein folding and stability.
   Progress in Biophysics and Molecular Biology, 86(2), 235-277.
   Note: Solubility is minimal at pI, increasing aggregation risk when buffer pH ≈ pI.

6. pI Calculation: 
   Bjellqvist, B., et al. (1993). 
   The focusing positions of polypeptides in immobilized pH gradients can be predicted 
   from their amino acid sequences. Electrophoresis, 14(1), 1023-1031.

7. Protein Solubility at pI: 
   Shaw, K. L., et al. (2001). 
   The effect of net charge on the solubility, activity, and stability of ribonuclease Sa.
   Protein Science, 10(6), 1206-1215.
   Note: Proteins exhibit minimal solubility near their isoelectric point.

8. Experimental Loss Factors: 
   Scopes, R. K. (2010). 
   Protein Purification: Principles and Practice (3rd ed.).
   Springer Science & Business Media.
   Note: Typical recovery losses include pipetting errors (2-5%), dead volume (5-10%), 
   and surface adsorption (5-15%), resulting in overall systemic efficiency of 70-80%.

9. Monte Carlo Uncertainty Quantification: 
   Rubinstein, R. Y., & Kroese, D. P. (2016). 
   Simulation and the Monte Carlo Method (3rd ed.). 
   John Wiley & Sons.

10. Empirical Calibration: 
    The systemic efficiency parameter (default 0.75) is calibrated 
    based on empirical observations from protein purification experiments, accounting for 
    cumulative losses from pipetting, dead volumes, and non-specific surface adsorption.

================================================================================
📝 버전 정보
================================================================================

Version: 2.4.0
Release Date: 2024

Changes in v2.4:
- Enhanced documentation for beginners
- Added comprehensive citation requirements
- Updated license to AGPL-3.0 for maximum developer rights protection
- Added user-configurable parameter markers
- Improved code comments and explanations

Changes in v2.3:
- Added Systemic Loss Factor: Accounts for experimental handling losses 
  (pipetting, dead volume, surface adsorption)
- Added pI-pH Interaction: Considers isoelectric point effects on solubility and aggregation
- Updated Instability Threshold: Changed from 20.0 to 40.0 based on Guruprasad et al. (1990)
- Enhanced References: Comprehensive literature support for all model components
- Improved Calibration: Model parameters tuned to match empirical recovery data (0.3-0.4 range)
- Enhanced Logging: More detailed reporting including pI and systemic efficiency factors

Changes in v2.1:
- Enhanced FASTA file format validation
- Automatic exclusion of project files (requirements.txt, README.md, etc.)
- Content-based file format verification
- Additional scientific references
- Improved error handling and logging

================================================================================
🚀 빠른 시작 가이드 (Quick Start Guide)
================================================================================

초보자를 위한 간단한 사용법:

1. 필수 라이브러리 설치:
   pip install biopython numpy
   pip install pyyaml  # 선택사항 (YAML 설정 파일 사용 시)

2. 기본 사용법:
   # Python 스크립트에서
   from eprm_analyzer_v2_4 import EPRMAnalyzer
   
   # 분석기 생성
   analyzer = EPRMAnalyzer()
   
   # 단백질 서열 분석
   result = analyzer.calculate_eprm("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDPSIHAGHSVEVLELKP")
   
   # 결과 확인
   print(f"예상 농도: {result['C_Effective_uM']} uM")

3. FASTA 파일 일괄 처리:
   analyzer = EPRMAnalyzer()
   analyzer.process_files()  # 현재 디렉토리의 .fasta 파일 자동 처리

더 자세한 사용법은 README.md 파일을 참고하세요.

================================================================================
⚙️ 사용자 변경 가능 파라미터 (User-Configurable Parameters)
================================================================================

다음 파라미터들은 사용자가 실험 조건에 맞게 쉽게 변경할 수 있습니다.
코드 내에서 "🔧 USER CONFIGURABLE" 주석으로 표시되어 있습니다.

주요 변경 가능 파라미터:
1. initial_conc_um: 초기 단백질 농도 (기본값: 10.0 uM)
2. initial_vol_ul: 초기 부피 (기본값: 90.0 uL)
3. final_vol_ul: 최종 부피 (기본값: 450.0 uL)
4. eta_kit: 키트 효율 (기본값: 0.50)
5. systemic_efficiency: 시스템 효율 (기본값: 0.75)
6. buffer_ph: 버퍼 pH (기본값: 7.4)
7. instability_threshold: 불안정성 임계값 (기본값: 40.0)
8. instability_penalty_factor: 불안정성 페널티 계수 (기본값: 80.0)
9. gravy_penalty_factor: GRAVY 페널티 계수 (기본값: 0.15)

자세한 설명은 각 파라미터의 독스트링을 참고하세요.

================================================================================
"""

import os
import json
import logging
import numpy as np
from datetime import datetime
from glob import glob
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import warnings

# ============================================================================
# 외부 라이브러리 의존성 확인
# ============================================================================
# Biopython: 단백질 서열 분석을 위한 필수 라이브러리
# 설치 방법: pip install biopython
try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    raise ImportError(
        "Biopython 라이브러리가 설치되지 않았습니다. "
        "'pip install biopython'을 실행해주세요.\n"
        "설치 가이드: https://biopython.org/wiki/Download"
    )

# PyYAML: YAML 설정 파일을 읽기 위한 선택적 라이브러리
# 설치 방법: pip install pyyaml
# 없어도 JSON 설정 파일은 사용 가능
try:
    import yaml
except ImportError:
    yaml = None
    warnings.warn(
        "PyYAML이 설치되지 않았습니다. YAML 설정 파일을 사용할 수 없습니다. "
        "JSON 설정 파일은 계속 사용 가능합니다."
    )


# ============================================================================
# 상수 정의 (Constants)
# ============================================================================

# 표준 아미노산 (20가지)
# A=Alanine, C=Cysteine, D=Aspartic acid, E=Glutamic acid, F=Phenylalanine,
# G=Glycine, H=Histidine, I=Isoleucine, K=Lysine, L=Leucine, M=Methionine,
# N=Asparagine, P=Proline, Q=Glutamine, R=Arginine, S=Serine, T=Threonine,
# V=Valine, W=Tryptophan, Y=Tyrosine
STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

# 제외할 파일명 패턴 (프로젝트 파일 등)
# FASTA 파일이 아닌 프로젝트 관리 파일들을 자동으로 제외합니다.
EXCLUDED_FILES = {
    'requirements.txt', 'readme.md', 'readme.txt', 'license', 'license.txt',
    'config.yaml', 'config.yml', 'config.json', '.gitignore', 'setup.py',
    'pyproject.toml', 'pom.xml', 'package.json', 'cargo.toml'
}


class EPRMAnalyzer:
    """
    Effective Protein Recovery Mass (EPRM) Analyzer v2.4.
    
    단백질 서열의 물리화학적 특성(불안정성 지수, GRAVY, 등전점)과 정제 키트의 효율을 기반으로
    실험 후 회수될 유효 단백질 농도를 예측하고 실험 가이드를 제공하는 클래스입니다.
    
    이 클래스는 다음과 같은 기능을 제공합니다:
    1. 단백질 서열 분석 (분자량, 불안정성 지수, GRAVY, 등전점 계산)
    2. 회수율 예측 (키트 효율, 시스템 효율, 단백질 특성 기반)
    3. 불확실성 정량화 (Monte Carlo 시뮬레이션)
    4. FASTA 파일 일괄 처리
    5. 상세한 로깅 및 결과 저장
    
    주요 개선사항 (v2.4):
    - 초보자 친화적인 문서화
    - 상세한 인용 가이드라인
    - AGPL-3.0 라이선스 적용 (개발자 권익 최대 보장)
    - 사용자 변경 가능 파라미터 명확히 표시
    
    주요 개선사항 (v2.3):
    1. Systemic Loss Factor: 피펫팅, Dead Volume, 튜브 벽면 흡착 등 실험적 손실 반영
       (기본값: 0.75, Scopes 2010 기준)
    2. pI-pH Interaction: 버퍼 pH와 단백질 pI의 차이에 따른 용해도 감소 반영
       (Gromiha et al. 2004, Shaw et al. 2001)
    3. Instability Threshold 업데이트: 20.0 → 40.0 (Guruprasad et al. 1990 기준)
    4. Empirical Calibration: 실험 데이터 기반 파라미터 보정 (회수 계수 0.3-0.4 범위)
    
    사용 예시:
        >>> analyzer = EPRMAnalyzer()
        >>> result = analyzer.calculate_eprm("MKTAYIAKQR")
        >>> print(result['C_Effective_uM'])
    """

    def __init__(
        self,
        initial_conc_um: float = 10.0,          # 🔧 USER CONFIGURABLE
        initial_vol_ul: float = 90.0,           # 🔧 USER CONFIGURABLE
        final_vol_ul: float = 450.0,            # 🔧 USER CONFIGURABLE
        eta_kit: float = 0.50,                 # 🔧 USER CONFIGURABLE
        systemic_efficiency: float = 0.75,     # 🔧 USER CONFIGURABLE
        buffer_ph: float = 7.4,                 # 🔧 USER CONFIGURABLE
        instability_threshold: float = 40.0,   # 🔧 USER CONFIGURABLE
        instability_penalty_factor: float = 80.0, # 🔧 USER CONFIGURABLE
        gravy_penalty_factor: float = 0.15,     # 🔧 USER CONFIGURABLE
        config_file: Optional[str] = None,
        random_seed: Optional[int] = None
    ):
        """
        초기화 메서드: 실험 조건 설정 및 결과 디렉토리 생성.
        
        이 메서드는 분석에 필요한 모든 파라미터를 설정하고, 결과를 저장할 디렉토리를 생성합니다.
        설정 파일을 사용하거나 직접 파라미터를 지정할 수 있습니다.
        
        Args:
            initial_conc_um (float): 
                초기 단백질 농도 (단위: 마이크로몰, uM). 
                기본값: 10.0 uM
                예시: 정제 전 단백질 농도가 10 uM이면 10.0을 입력
                
            initial_vol_ul (float): 
                초기 부피 (단위: 마이크로리터, uL).
                기본값: 90.0 uL
                예시: 정제 전 샘플 부피가 90 uL이면 90.0을 입력
                
            final_vol_ul (float): 
                최종 희석/용출 부피 (단위: 마이크로리터, uL).
                기본값: 450.0 uL
                예시: 정제 후 최종 부피가 450 uL이면 450.0을 입력
                
            eta_kit (float): 
                키트 자체의 이론적 최대 효율 (0.0 ~ 1.0 사이의 값).
                기본값: 0.50 (50%)
                설명: 사용하는 정제 키트의 이론적 최대 회수율입니다.
                예시: 키트 제조사가 50% 회수율을 명시하면 0.50 입력
                
            systemic_efficiency (float): 
                실험적 조작 효율 (0.0 ~ 1.0 사이의 값).
                기본값: 0.75 (75%, Scopes 2010 기준)
                설명: 피펫팅 오차, Dead Volume, 표면 흡착 등으로 인한 손실을 반영합니다.
                일반적으로 0.70 ~ 0.80 범위입니다.
                예시: 실험 환경이 깨끗하고 정밀하면 0.80, 일반적이면 0.75
                
            buffer_ph (float): 
                실험 버퍼의 pH 값.
                기본값: 7.4 (PBS 기준)
                설명: 단백질의 등전점(pI)과 버퍼 pH의 차이에 따라 용해도가 달라집니다.
                예시: PBS 사용 시 7.4, Tris-HCl 사용 시 해당 pH 값 입력
                
            instability_threshold (float): 
                불안정성 지수 임계값.
                기본값: 40.0 (Guruprasad et al. 1990 기준)
                설명: 이 값보다 높으면 불안정한 단백질로 간주됩니다.
                일반적으로 변경할 필요 없습니다.
                
            instability_penalty_factor (float): 
                불안정성 페널티 계수.
                기본값: 80.0
                설명: 불안정성 지수가 임계값을 초과할 때 적용되는 페널티의 강도를 조절합니다.
                값이 작을수록 페널티가 강해집니다.
                
            gravy_penalty_factor (float): 
                GRAVY 페널티 계수.
                기본값: 0.15
                설명: 소수성(GRAVY 값)이 높을수록 표면 흡착이 증가하므로 
                이 계수로 손실을 조절합니다.
                
            config_file (Optional[str]): 
                설정 파일 경로 (YAML 또는 JSON 형식).
                None이면 기본값 사용.
                예시: "config.json" 또는 "config.yaml"
                
            random_seed (Optional[int]): 
                재현성을 위한 랜덤 시드.
                None이면 매번 다른 결과 (재현 불가).
                정수 값을 주면 항상 같은 결과 (재현 가능).
                예시: 42 (재현성을 위해 권장)
        
        Raises:
            ValueError: 파라미터가 유효하지 않은 경우 (예: 음수 값, 범위 초과).
            FileNotFoundError: 설정 파일을 찾을 수 없는 경우.
        
        사용 예시:
            >>> # 기본 설정으로 생성
            >>> analyzer = EPRMAnalyzer()
            
            >>> # 사용자 정의 파라미터로 생성
            >>> analyzer = EPRMAnalyzer(
            ...     initial_conc_um=20.0,  # 20 uM
            ...     initial_vol_ul=100.0,  # 100 uL
            ...     final_vol_ul=500.0,   # 500 uL
            ...     buffer_ph=7.0,         # pH 7.0
            ...     random_seed=42         # 재현성 보장
            ... )
            
            >>> # 설정 파일 사용
            >>> analyzer = EPRMAnalyzer(config_file="my_config.json")
        """
        # 설정 파일 로드 (있는 경우)
        # 설정 파일이 있으면 파일의 값을 우선 사용하고, 없으면 기본값 사용
        if config_file:
            self._load_config(config_file)
        else:
            # 🔧 USER CONFIGURABLE: 실험 파라미터
            # 아래 값들을 실험 조건에 맞게 변경하세요
            self.c_start = initial_conc_um
            self.v_start = initial_vol_ul
            self.v_final = final_vol_ul
            self.eta_kit = eta_kit
            self.systemic_efficiency = systemic_efficiency
            self.buffer_ph = buffer_ph
            self.instability_threshold = instability_threshold
            self.instability_penalty_factor = instability_penalty_factor
            self.gravy_penalty_factor = gravy_penalty_factor

        # 파라미터 검증 (잘못된 값이 입력되면 에러 발생)
        self._validate_parameters()

        # 랜덤 시드 설정 (재현성 보장)
        # 같은 시드를 사용하면 항상 같은 결과를 얻을 수 있습니다.
        if random_seed is not None:
            np.random.seed(random_seed)
            self.random_seed = random_seed
        else:
            self.random_seed = None

        # 결과 저장소 설정 (Timestamp 기반 폴더링)
        # 실행할 때마다 새로운 폴더가 생성되어 결과가 덮어씌워지지 않습니다.
        # 폴더명 형식: EPRM_Results_20241201_143022 (날짜_시간)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = f"EPRM_Results_{self.timestamp}"
        os.makedirs(self.output_dir, exist_ok=True)

        # 로깅 시스템 초기화
        # 분석 과정을 파일과 콘솔에 동시에 기록합니다.
        self.log_path = os.path.join(self.output_dir, "eprm_analysis_detail.log")
        self._setup_logging()

        # 설정 저장 (재현성)
        # 나중에 같은 설정으로 재현할 수 있도록 모든 파라미터를 저장합니다.
        self._save_config()

    def _load_config(self, config_file: str) -> None:
        """
        설정 파일 로드.
        
        YAML 또는 JSON 형식의 설정 파일을 읽어서 파라미터를 설정합니다.
        
        설정 파일 예시 (config.json):
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
        """
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"설정 파일을 찾을 수 없습니다: {config_file}")

        with open(config_path, 'r', encoding='utf-8') as f:
            if config_path.suffix.lower() in ['.yaml', '.yml']:
                if yaml is None:
                    raise ImportError(
                        "YAML 파일을 읽으려면 PyYAML이 필요합니다. "
                        "'pip install pyyaml'을 실행해주세요."
                    )
                config = yaml.safe_load(f)
            else:
                config = json.load(f)

        # 설정 적용 (설정 파일에 없는 값은 기본값 사용)
        self.c_start = config.get('initial_conc_um', 10.0)
        self.v_start = config.get('initial_vol_ul', 90.0)
        self.v_final = config.get('final_vol_ul', 450.0)
        self.eta_kit = config.get('eta_kit', 0.50)
        self.systemic_efficiency = config.get('systemic_efficiency', 0.75)
        self.buffer_ph = config.get('buffer_ph', 7.4)
        self.instability_threshold = config.get('instability_threshold', 40.0)
        self.instability_penalty_factor = config.get('instability_penalty_factor', 80.0)
        self.gravy_penalty_factor = config.get('gravy_penalty_factor', 0.15)

    def _save_config(self) -> None:
        """
        현재 설정을 파일로 저장 (재현성).
        
        분석에 사용된 모든 파라미터를 JSON 파일로 저장하여 
        나중에 같은 조건으로 재현할 수 있도록 합니다.
        """
        config = {
            'initial_conc_um': self.c_start,
            'initial_vol_ul': self.v_start,
            'final_vol_ul': self.v_final,
            'eta_kit': self.eta_kit,
            'systemic_efficiency': self.systemic_efficiency,
            'buffer_ph': self.buffer_ph,
            'instability_threshold': self.instability_threshold,
            'instability_penalty_factor': self.instability_penalty_factor,
            'gravy_penalty_factor': self.gravy_penalty_factor,
            'random_seed': self.random_seed,
            'timestamp': self.timestamp,
            'version': '2.4.0'  # 버전 정보 추가
        }
        config_path = os.path.join(self.output_dir, "config.json")
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)

    def _validate_parameters(self) -> None:
        """
        파라미터 유효성 검증.
        
        입력된 파라미터가 논리적으로 올바른지 확인합니다.
        잘못된 값이 있으면 ValueError를 발생시킵니다.
        """
        if self.c_start <= 0:
            raise ValueError(
                f"초기 농도는 양수여야 합니다: {self.c_start} uM. "
                "현재 값이 0 이하입니다. 올바른 농도 값을 입력해주세요."
            )
        if self.v_start <= 0:
            raise ValueError(
                f"초기 부피는 양수여야 합니다: {self.v_start} uL. "
                "현재 값이 0 이하입니다. 올바른 부피 값을 입력해주세요."
            )
        if self.v_final <= 0:
            raise ValueError(
                f"최종 부피는 양수여야 합니다: {self.v_final} uL. "
                "현재 값이 0 이하입니다. 올바른 부피 값을 입력해주세요."
            )
        if not (0 < self.eta_kit <= 1):
            raise ValueError(
                f"키트 회수율은 0과 1 사이여야 합니다: {self.eta_kit}. "
                "현재 값이 범위를 벗어났습니다. 0.0 ~ 1.0 사이의 값을 입력해주세요."
            )
        if not (0 < self.systemic_efficiency <= 1):
            raise ValueError(
                f"시스템 효율은 0과 1 사이여야 합니다: {self.systemic_efficiency}. "
                "현재 값이 범위를 벗어났습니다. 0.0 ~ 1.0 사이의 값을 입력해주세요."
            )
        # 최종 부피가 초기 부피보다 작으면 농축 과정으로 간주 (경고만 출력)
        if self.v_final < self.v_start:
            logging.warning(
                f"최종 부피({self.v_final} uL)가 초기 부피({self.v_start} uL)보다 작습니다. "
                "농축 과정을 가정합니다. 이는 정상적인 상황일 수 있습니다."
            )

    def _validate_sequence(self, sequence: str) -> Tuple[bool, Optional[str]]:
        """
        단백질 서열 유효성 검증.
        
        입력된 서열이 유효한 단백질 서열인지 확인합니다.
        
        검증 항목:
        1. 서열이 비어있지 않은지
        2. 표준 아미노산(20가지)만 포함하는지
        3. 최소 길이(2개 아미노산) 이상인지
        
        Args:
            sequence (str): 검증할 서열 (예: "MKTAYIAKQR")
        
        Returns:
            Tuple[bool, Optional[str]]: 
                - (True, None): 유효한 서열
                - (False, "에러 메시지"): 유효하지 않은 서열
        """
        if not sequence:
            return False, "서열이 비어있습니다. 단백질 서열을 입력해주세요."

        # 표준 아미노산만 포함하는지 확인
        # 서열을 대문자로 변환하고, 표준 아미노산 집합에 없는 문자가 있는지 확인
        invalid_chars = set(sequence.upper()) - STANDARD_AMINO_ACIDS
        if invalid_chars:
            return False, (
                f"비표준 아미노산이 포함되어 있습니다: {invalid_chars}. "
                f"표준 아미노산은 {''.join(sorted(STANDARD_AMINO_ACIDS))} 입니다."
            )

        # 최소 길이 확인 (너무 짧은 서열은 분석 불가)
        if len(sequence) < 2:
            return False, (
                f"서열이 너무 짧습니다 (현재 길이: {len(sequence)}). "
                "최소 2개 아미노산이 필요합니다."
            )

        return True, None

    def _is_fasta_file(self, file_path: str) -> bool:
        """
        파일이 FASTA 형식인지 확인합니다.
        
        FASTA 형식의 특징:
        1. '>'로 시작하는 헤더 라인이 있음
        2. 헤더 다음에 아미노산 서열이 있음
        3. 서열은 표준 아미노산 문자로만 구성
        
        Args:
            file_path (str): 확인할 파일 경로
        
        Returns:
            bool: FASTA 형식이면 True, 아니면 False
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read(1000)  # 처음 1000자만 읽어서 확인 (성능 향상)
                
            lines = content.split('\n')
            has_header = False
            has_sequence = False
            sequence_chars = 0
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                    
                # FASTA 헤더 확인 ('>'로 시작)
                if line.startswith('>'):
                    has_header = True
                    continue
                
                # 서열 라인 확인
                line_upper = line.upper().replace(' ', '').replace('\t', '')
                if line_upper:
                    # 표준 아미노산 비율 확인
                    valid_chars = sum(1 for c in line_upper if c in STANDARD_AMINO_ACIDS)
                    if len(line_upper) > 0:
                        valid_ratio = valid_chars / len(line_upper)
                        # 80% 이상이 표준 아미노산이면 서열로 간주
                        if valid_ratio > 0.8:
                            has_sequence = True
                            sequence_chars += valid_chars
            
            # FASTA 파일 조건:
            # 1. 헤더가 있거나
            # 2. 서열 문자가 충분히 많고 (최소 10자)
            # 3. 전체 내용의 50% 이상이 표준 아미노산 문자
            total_chars = len(content.replace(' ', '').replace('\n', '').replace('\t', ''))
            if total_chars > 0:
                sequence_ratio = sequence_chars / total_chars
                return (has_header or has_sequence) and sequence_chars >= 10 and sequence_ratio > 0.5
            
            return False
            
        except Exception:
            # 파일 읽기 오류 시 False 반환
            return False

    def _should_exclude_file(self, file_path: str) -> bool:
        """
        파일을 제외해야 하는지 확인합니다.
        
        프로젝트 관리 파일(README, requirements.txt 등)은 자동으로 제외합니다.
        
        Args:
            file_path (str): 확인할 파일 경로
        
        Returns:
            bool: 제외해야 하면 True, 아니면 False
        """
        file_name = os.path.basename(file_path).lower()
        
        # 제외 목록에 있는 파일
        if file_name in EXCLUDED_FILES:
            return True
        
        # 특정 확장자 제외
        excluded_extensions = {'.md', '.py', '.json', '.yaml', '.yml', '.toml', 
                               '.xml', '.log', '.gitignore', '.txt'}
        # 단, .fasta와 .txt는 FASTA 형식일 수 있으므로 내용 확인 필요
        if file_path.endswith('.fasta'):
            return False  # .fasta는 항상 처리
        
        # .txt 파일은 내용 확인
        if file_path.endswith('.txt'):
            # 제외 목록에 있으면 제외
            if file_name in EXCLUDED_FILES:
                return True
            # FASTA 형식이 아니면 제외
            if not self._is_fasta_file(file_path):
                return True
            return False
        
        # 다른 제외 확장자
        for ext in excluded_extensions:
            if file_path.lower().endswith(ext):
                return True
        
        return False

    def _setup_logging(self) -> None:
        """
        로깅 핸들러 설정 (File + Stream).
        
        분석 과정을 파일과 콘솔에 동시에 기록합니다.
        기존 핸들러를 제거하여 로그 중복 출력을 방지합니다.
        """
        logger = logging.getLogger()

        # 기존 핸들러 초기화 (중복 방지)
        if logger.hasHandlers():
            logger.handlers.clear()

        # 로깅 설정
        logging.basicConfig(
            level=logging.INFO,  # INFO 레벨 이상의 로그만 기록
            format='%(asctime)s - [%(levelname)s] - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_path, encoding='utf-8'),  # 파일 출력
                logging.StreamHandler()  # 콘솔 출력
            ]
        )
        logging.info(f"EPRM Analysis v2.4 Started. Output Directory: {self.output_dir}")
        logging.info(
            f"Configuration: c_start={self.c_start} uM, "
            f"v_start={self.v_start} uL, v_final={self.v_final} uL, "
            f"eta_kit={self.eta_kit}, systemic_efficiency={self.systemic_efficiency}, "
            f"buffer_pH={self.buffer_ph}"
        )

    def calculate_eprm(
        self,
        sequence: str,
        include_uncertainty: bool = True,
        n_iterations: int = 1000  # 🔧 USER CONFIGURABLE: 시뮬레이션 반복 횟수
    ) -> Dict[str, Union[float, Tuple[float, float]]]:
        """
        단백질 서열을 분석하여 예측 회수율과 유효 농도를 계산합니다.
        
        이 메서드는 EPRM 분석의 핵심 기능입니다. 단백질 서열을 입력받아
        물리화학적 특성을 분석하고, 회수 농도를 예측합니다.
        
        계산 과정:
        1. Biopython을 이용해 MW, Instability Index, GRAVY, pI 계산
        2. 안정성 계수(stab_factor), 흡착 계수(ads_factor), pI-pH 상호작용(pi_factor) 계산
        3. 키트 효율(eta_kit), 시스템 효율(systemic_efficiency), 단백질 효율(eta_prot) 결합
        4. (선택적) Monte Carlo 시뮬레이션을 통한 불확실성 정량화
        
        Args:
            sequence (str): 
                아미노산 서열 문자열 (예: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDPSIHAGHSVEVLELKP")
                표준 아미노산 20가지만 포함해야 합니다.
                
            include_uncertainty (bool): 
                불확실성 계산 포함 여부. 
                True: Monte Carlo 시뮬레이션 수행 (시간이 더 걸리지만 정확함)
                False: 단순 계산만 수행 (빠름)
                기본값: True
                
            n_iterations (int): 
                Monte Carlo 시뮬레이션 반복 횟수.
                값이 클수록 정확하지만 시간이 더 걸립니다.
                기본값: 1000
                권장 범위: 500 ~ 10000
        
        Returns:
            Dict[str, Union[float, Tuple[float, float]]]: 분석 결과 딕셔너리
                - MW_kDa: 분자량 (단위: kDa)
                - Instability: 불안정성 지수 (40 이상이면 불안정)
                - GRAVY: 소수성 지수 (음수=친수성, 양수=소수성)
                - pI: 등전점 (Isoelectric Point)
                - Stab_Factor: 안정성 계수 (1.0에 가까울수록 안정)
                - Ads_Factor: 흡착 계수 (1.0에 가까울수록 흡착 적음)
                - pI_Factor: pI-pH 상호작용 계수 (1.0에 가까울수록 좋음)
                - Eta_Prot: 단백질 효율 계수 (0.0 ~ 1.0)
                - Total_Coeff: 총 회수 계수 (eta_kit × systemic_efficiency × eta_prot)
                - C_Theo_Max_uM: 이론적 최대 농도 (단순 희석만 고려, 단위: uM)
                - C_Effective_uM: 유효 농도 (보정 적용, 단위: uM)
                    - include_uncertainty=False: float 값
                    - include_uncertainty=True: (평균, 표준편차) 튜플
                - C_Effective_CI_95: 95% 신뢰구간 (include_uncertainty=True인 경우만)
                    - (하한값, 상한값) 튜플
        
        Raises:
            ValueError: 서열이 유효하지 않은 경우 (비표준 아미노산 포함 등)
            RuntimeError: 단백질 분석 중 오류 발생 (Biopython 오류)
        
        사용 예시:
            >>> analyzer = EPRMAnalyzer()
            >>> 
            >>> # 불확실성 포함 분석
            >>> result = analyzer.calculate_eprm("MKTAYIAKQR", include_uncertainty=True)
            >>> print(f"예상 농도: {result['C_Effective_uM'][0]:.4f} ± {result['C_Effective_uM'][1]:.4f} uM")
            >>> 
            >>> # 빠른 분석 (불확실성 제외)
            >>> result = analyzer.calculate_eprm("MKTAYIAKQR", include_uncertainty=False)
            >>> print(f"예상 농도: {result['C_Effective_uM']:.4f} uM")
        """
        # 서열 검증
        is_valid, error_msg = self._validate_sequence(sequence)
        if not is_valid:
            raise ValueError(f"유효하지 않은 서열: {error_msg}")

        # 1. 기초 물성 분석 (Biopython 사용)
        # Biopython의 ProteinAnalysis 클래스를 사용하여 단백질의 물리화학적 특성을 계산합니다.
        try:
            analysis = ProteinAnalysis(sequence)
            mw = analysis.molecular_weight() / 1000.0  # Da -> kDa 변환
            instability_index = analysis.instability_index()  # 불안정성 지수
            gravy = analysis.gravy()  # GRAVY (소수성 지수)
            pI = analysis.isoelectric_point()  # 등전점
        except Exception as e:
            raise RuntimeError(
                f"단백질 분석 중 오류 발생: {str(e)}\n"
                "서열 형식이 올바른지 확인해주세요."
            )

        # 2. 물리화학적 감쇄 로직 (Refined)
        # 단백질의 특성에 따라 회수율이 달라지므로, 이를 계수로 반영합니다.
        
        # A. Stability Factor (Guruprasad et al. 1990)
        # 불안정성 지수가 임계값(40)을 넘으면 페널티를 부여합니다.
        # Instability Index > 40 is unstable. Linear penalty applied for values above threshold.
        # 수식: 1.0 - (초과분 / penalty_factor)
        # 예: instability_index=50, threshold=40, penalty_factor=80
        #     stab_factor = 1.0 - (50-40)/80 = 1.0 - 0.125 = 0.875
        stab_factor = 1.0 - (
            max(0, instability_index - self.instability_threshold) /
            self.instability_penalty_factor
        )

        # B. Adsorption Factor (Norde, 1986)
        # 소수성(GRAVY 값)이 높을수록 튜브 벽면 등에 흡착되어 손실이 증가합니다.
        # Hydrophobic proteins show increased adsorption. 
        # Basal loss is handled by systemic_efficiency; here we calculate additional hydrophobicity-related loss.
        # 수식: 1.0 - (|GRAVY| * penalty_factor)
        # 예: gravy=0.5, penalty_factor=0.15
        #     ads_factor = 1.0 - (0.5 * 0.15) = 1.0 - 0.075 = 0.925
        ads_factor = 1.0 - (abs(gravy) * self.gravy_penalty_factor)

        # C. pI-pH Solubility Factor (Gromiha et al. 2004, Shaw et al. 2001) [NEW]
        # 버퍼 pH가 단백질 pI에 가까우면 순전하가 0에 가까워져 응집 위험이 증가합니다.
        # When pH ≈ pI, net charge approaches zero, increasing aggregation risk.
        # Gaussian-like penalty: maximum 15% loss when pH = pI, decreasing as |pH - pI| increases.
        delta_ph = abs(self.buffer_ph - pI)
        # When delta_ph = 0, maximum 15% loss (0.85); approaches 1.0 as delta_ph increases
        # 수식: 1.0 - (0.15 * exp(-(delta_ph^2) / 2.0))
        # 예: buffer_ph=7.4, pI=7.4 (delta_ph=0)
        #     pi_factor = 1.0 - (0.15 * exp(0)) = 1.0 - 0.15 = 0.85
        # 예: buffer_ph=7.4, pI=9.0 (delta_ph=1.6)
        #     pi_factor = 1.0 - (0.15 * exp(-1.6^2/2)) ≈ 1.0 - 0.15*0.28 ≈ 0.96
        pi_factor = 1.0 - (0.15 * np.exp(-(delta_ph**2) / 2.0))

        # D. Combined Protein Efficiency
        # 세 가지 계수를 곱하여 단백질 고유 효율을 계산합니다.
        # 0 이하가 되지 않도록 max(0.0, ...)로 보호합니다.
        eta_prot = max(0.0, stab_factor * ads_factor * pi_factor)

        # 3. 최종 농도 계산
        # Total Coeff = Kit_Max × Systemic_Handling × Protein_Specifics
        # Systemic efficiency (default 0.75) accounts for tube binding, dead volume, pipetting errors.
        total_recovery_coeff = self.eta_kit * self.systemic_efficiency * eta_prot

        # 이론적 최대 희석 농도 (단순 희석만 고려)
        # 수식: (초기 농도 × 초기 부피) / 최종 부피
        # 예: c_start=10 uM, v_start=90 uL, v_final=450 uL
        #     c_theo_max = (10 × 90) / 450 = 2.0 uM
        c_theo_max = (self.c_start * self.v_start) / self.v_final

        # 보정된 유효 농도 (Effective Concentration)
        # 이론적 최대 농도에 회수 계수를 곱하여 실제 예상 농도를 계산합니다.
        c_eff = c_theo_max * total_recovery_coeff

        # 결과 딕셔너리 구성
        result = {
            "MW_kDa": mw,
            "Instability": instability_index,
            "GRAVY": gravy,
            "pI": pI,
            "Stab_Factor": stab_factor,
            "Ads_Factor": ads_factor,
            "pI_Factor": pi_factor,
            "Eta_Prot": eta_prot,
            "Total_Coeff": total_recovery_coeff,
            "C_Theo_Max_uM": c_theo_max,
            "C_Effective_uM": c_eff
        }

        # 4. 불확실성 정량화 (Monte Carlo 시뮬레이션)
        # 파라미터에 불확실성이 있다고 가정하고, 여러 번 시뮬레이션하여
        # 평균, 표준편차, 신뢰구간을 계산합니다.
        if include_uncertainty:
            uncertainty_result = self._calculate_uncertainty(
                sequence, instability_index, gravy, pI, n_iterations
            )
            result["C_Effective_uM"] = (
                uncertainty_result["mean"],
                uncertainty_result["std"]
            )
            result["C_Effective_CI_95"] = uncertainty_result["ci_95"]

        return result

    def _calculate_uncertainty(
        self,
        sequence: str,
        instability_index: float,
        gravy: float,
        pI: float,
        n_iterations: int = 1000
    ) -> Dict[str, Union[float, Tuple[float, float]]]:
        """
        Monte Carlo 시뮬레이션을 통한 불확실성 정량화.
        
        파라미터에 노이즈를 추가하여 여러 번 시뮬레이션하고 통계를 계산합니다.
        이 방법을 통해 예측값의 불확실성을 정량화할 수 있습니다.
        
        Args:
            sequence (str): 단백질 서열 (현재는 사용하지 않지만 향후 확장 가능)
            instability_index (float): 계산된 불안정성 지수
            gravy (float): 계산된 GRAVY 값
            pI (float): 계산된 등전점 (Isoelectric Point)
            n_iterations (int): 시뮬레이션 반복 횟수 (기본값: 1000)
        
        Returns:
            Dict[str, Union[float, Tuple[float, float]]]: 불확실성 통계
                - mean: 평균값
                - std: 표준편차
                - ci_95: 95% 신뢰구간 (하한값, 상한값) 튜플
        """
        results = []

        # 파라미터 불확실성 가정 (표준편차)
        # 실제 실험에서는 파라미터 값에 불확실성이 있습니다.
        # 이를 정규분포로 모델링합니다.
        instability_std = instability_index * 0.05  # 5% 불확실성
        gravy_std = 0.1  # 10% 불확실성 (고정값)
        eta_kit_std = self.eta_kit * 0.05  # 5% 불확실성
        sys_std = self.systemic_efficiency * 0.05  # 5% 불확실성 (실험적 조작 변동)

        # Monte Carlo 시뮬레이션 반복
        for _ in range(n_iterations):
            # 파라미터에 노이즈 추가 (정규분포에서 샘플링)
            inst_perturbed = max(0, np.random.normal(instability_index, instability_std))
            gravy_perturbed = np.random.normal(gravy, gravy_std)
            # 0~1 범위로 제한
            eta_kit_perturbed = np.clip(
                np.random.normal(self.eta_kit, eta_kit_std), 0, 1
            )
            sys_perturbed = np.clip(
                np.random.normal(self.systemic_efficiency, sys_std), 0, 1
            )

            # 보정 계수 재계산 (노이즈가 추가된 파라미터로)
            stab_factor = 1.0 - (
                max(0, inst_perturbed - self.instability_threshold) /
                self.instability_penalty_factor
            )
            ads_factor = 1.0 - (abs(gravy_perturbed) * self.gravy_penalty_factor)
            
            # pI는 상대적으로 고정값으로 가정 (단백질 고유 특성)
            delta_ph = abs(self.buffer_ph - pI)
            pi_factor = 1.0 - (0.15 * np.exp(-(delta_ph**2) / 2.0))
            
            eta_prot = max(0, stab_factor * ads_factor * pi_factor)
            total_recovery_coeff = eta_kit_perturbed * sys_perturbed * eta_prot
            
            # 농도 계산
            c_theo_max = (self.c_start * self.v_start) / self.v_final
            c_eff = c_theo_max * total_recovery_coeff

            results.append(c_eff)

        # 통계 계산
        results = np.array(results)
        mean = float(np.mean(results))
        std = float(np.std(results))
        ci_lower = float(np.percentile(results, 2.5))   # 2.5 백분위수
        ci_upper = float(np.percentile(results, 97.5))  # 97.5 백분위수

        return {
            "mean": mean,
            "std": std,
            "ci_95": (ci_lower, ci_upper)
        }

    def process_files(
        self,
        input_dir: Optional[str] = None,  # 🔧 USER CONFIGURABLE: 입력 디렉토리
        include_uncertainty: bool = True
    ) -> None:
        """
        지정된 디렉토리의 .fasta 및 .txt 파일을 찾아 분석을 수행하고 로그를 기록합니다.
        
        이 메서드는 여러 FASTA 파일을 한 번에 처리할 때 사용합니다.
        현재 디렉토리 또는 지정된 디렉토리에서 .fasta 및 FASTA 형식의 .txt 파일을
        자동으로 찾아서 분석합니다.
        
        처리 과정:
        1. 디렉토리에서 .fasta 및 .txt 파일 검색
        2. FASTA 형식 검증 (프로젝트 파일 자동 제외)
        3. 각 파일의 모든 서열 분석
        4. 결과를 JSON 파일로 저장
        5. 상세 로그를 파일로 저장
        
        Args:
            input_dir (Optional[str]): 
                입력 파일 디렉토리. 
                None이면 현재 디렉토리(".") 사용.
                예시: "./data", "C:/Users/Desktop/proteins"
                
            include_uncertainty (bool): 
                불확실성 계산 포함 여부.
                True: Monte Carlo 시뮬레이션 수행 (정확하지만 느림)
                False: 단순 계산만 수행 (빠름)
                기본값: True
        
        사용 예시:
            >>> analyzer = EPRMAnalyzer()
            >>> 
            >>> # 현재 디렉토리의 파일 처리
            >>> analyzer.process_files()
            >>> 
            >>> # 특정 디렉토리의 파일 처리
            >>> analyzer.process_files(input_dir="./my_proteins")
            >>> 
            >>> # 빠른 분석 (불확실성 제외)
            >>> analyzer.process_files(include_uncertainty=False)
        """
        if input_dir is None:
            input_dir = "."

        # 대상 파일 검색
        # glob 모듈을 사용하여 .fasta 및 .txt 파일을 찾습니다.
        all_candidate_files = (
            glob(os.path.join(input_dir, "*.fasta")) +
            glob(os.path.join(input_dir, "*.txt"))
        )

        # 파일 필터링: 제외 목록 및 형식 검증
        target_files = []
        excluded_files = []
        
        for file_path in all_candidate_files:
            # 현재 실행 중인 파이썬 스크립트 파일은 제외
            if file_path == os.path.basename(__file__):
                excluded_files.append((file_path, "Python script file"))
                continue
            
            # 제외 목록 확인
            if self._should_exclude_file(file_path):
                excluded_files.append((file_path, "Excluded file or invalid format"))
                continue
            
            # FASTA 형식 검증 (특히 .txt 파일)
            if file_path.endswith('.txt') and not self._is_fasta_file(file_path):
                excluded_files.append((file_path, "Not a valid FASTA format"))
                continue
            
            target_files.append(file_path)

        # 제외된 파일 로그 출력
        if excluded_files:
            logging.info(f"Excluded {len(excluded_files)} files:")
            for file_path, reason in excluded_files[:5]:  # 최대 5개만 로그
                logging.info(f"  - {os.path.basename(file_path)}: {reason}")
            if len(excluded_files) > 5:
                logging.info(f"  ... and {len(excluded_files) - 5} more files")

        # 유효한 파일이 없으면 에러 메시지 출력
        if not target_files:
            logging.error(
                f"No valid FASTA files found in {input_dir}. "
                "Please place sequence files (.fasta or FASTA-formatted .txt) in the directory."
            )
            return

        logging.info(f"Found {len(target_files)} valid FASTA file(s) to process.")

        # 결과 저장용 리스트
        all_results = []

        # 각 파일 처리
        for file_path in target_files:
            logging.info(f"{'='*10} Processing File: {file_path} {'='*10}")

            try:
                with open(file_path, "r", encoding='utf-8') as f:
                    # FASTA 포맷 파싱 ('>' 기준으로 분리)
                    # FASTA 형식: >헤더\n서열\n서열...
                    raw_content = f.read()
                    entries = raw_content.split('>')

                    valid_entries = 0
                    for entry in entries:
                        if not entry.strip():
                            continue  # 빈 항목 건너뛰기

                        lines = entry.strip().split('\n')
                        header = lines[0].strip() if lines else "Unknown"
                        # 줄바꿈 제거 및 공백 제거, 대문자 변환
                        seq = "".join(lines[1:]).replace(" ", "").strip().upper()

                        # 유효하지 않은 서열 건너뛰기
                        is_valid, error_msg = self._validate_sequence(seq)
                        if not is_valid:
                            logging.warning(f"Skipping invalid sequence '{header}': {error_msg}")
                            continue

                        try:
                            # 분석 실행
                            res = self.calculate_eprm(seq, include_uncertainty=include_uncertainty)
                            valid_entries += 1

                            # 결과 저장
                            result_entry = {
                                "file": file_path,
                                "header": header,
                                "sequence": seq,
                                "results": res
                            }
                            all_results.append(result_entry)

                            # --- 결과 리포팅 ---
                            logging.info(f"[Analysis Target: {header}]")
                            logging.info(
                                f"  • Properties: MW={res['MW_kDa']:.1f} kDa, "
                                f"pI={res['pI']:.2f}, GRAVY={res['GRAVY']:.2f}"
                            )
                            logging.info(
                                f"  • Instability Index: {res['Instability']:.2f} "
                                f"(Threshold: {self.instability_threshold})"
                            )
                            logging.info(
                                f"  • Coefficients: Kit({self.eta_kit:.2f}) × "
                                f"Sys({self.systemic_efficiency:.2f}) × "
                                f"Prot({res['Eta_Prot']:.3f})"
                            )
                            logging.info(
                                f"  • Final Recovery Coeff: {res['Total_Coeff']:.4f}"
                            )

                            # 불확실성 정보 출력
                            if include_uncertainty and "C_Effective_CI_95" in res:
                                mean, std = res["C_Effective_uM"]
                                ci_lower, ci_upper = res["C_Effective_CI_95"]
                                logging.info(
                                    f"  • >> Estimated Effective Conc: {mean:.4f} ± {std:.4f} uM"
                                )
                                logging.info(
                                    f"  • >> 95% CI: [{ci_lower:.4f}, {ci_upper:.4f}] uM"
                                )
                            else:
                                logging.info(
                                    f"  • >> Estimated Effective Conc: {res['C_Effective_uM']:.4f} uM"
                                )

                            # 실험 가이드: 20nM 타겟 희석비 계산
                            # 일반적으로 실험에서 20nM 농도를 목표로 하므로,
                            # 예상 농도에서 20nM으로 희석하는 배수를 계산합니다.
                            target_conc = 0.02  # 20 nM = 0.02 uM
                            if include_uncertainty:
                                c_eff_value = res["C_Effective_uM"][0]
                            else:
                                c_eff_value = res["C_Effective_uM"]

                            if c_eff_value > target_conc:
                                dilution_factor = int(c_eff_value / target_conc)
                                logging.info(
                                    f"  • [EXPERIMENTAL GUIDE] For 20nM final: Dilute 1:{dilution_factor}"
                                )
                                logging.info(
                                    f"    (Calculation: {c_eff_value:.4f} uM / 0.02 uM ≈ {dilution_factor})"
                                )
                            else:
                                logging.warning(
                                    "  • [GUIDE] Concentration too low (< 20nM) for standard dilution"
                                )
                            logging.info("-" * 50)

                        except Exception as e:
                            logging.error(f"Error analyzing sequence '{header}': {str(e)}")
                            continue

                    if valid_entries == 0:
                        logging.warning(f"No valid sequences found in {file_path}")

            except Exception as e:
                logging.error(f"Error processing file '{file_path}': {str(e)}")

        # 결과를 JSON 파일로 저장
        # 모든 분석 결과를 구조화된 JSON 형식으로 저장하여
        # 나중에 다른 프로그램에서 읽을 수 있도록 합니다.
        results_path = os.path.join(self.output_dir, "results.json")
        with open(results_path, 'w', encoding='utf-8') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False)

        logging.info(f"{'='*40}")
        logging.info(f"All analysis completed. Check details in: {self.log_path}")
        logging.info(f"Results saved to: {results_path}")


# ============================================================================
# 메인 실행 부분
# ============================================================================
# 이 파일을 직접 실행할 때 (python eprm_analyzer_v2.4.py) 실행되는 코드입니다.
# 다른 파일에서 import해서 사용할 때는 실행되지 않습니다.

if __name__ == "__main__":
    """
    메인 실행 예시.
    
    이 코드는 예시입니다. 사용자의 실험 조건에 맞게 파라미터를 수정하세요.
    """
    # 🔧 USER CONFIGURABLE: 아래 파라미터들을 실험 조건에 맞게 수정하세요
    
    # 인스턴스 생성 및 실행
    analyzer = EPRMAnalyzer(
        # 실험 조건 파라미터
        initial_conc_um=10.0,          # 초기 단백질 농도 (uM)
        initial_vol_ul=90.0,           # 초기 부피 (uL)
        final_vol_ul=450.0,            # 최종 부피 (uL)
        
        # 효율 파라미터
        eta_kit=0.50,                  # 키트 효율 (50%)
        systemic_efficiency=0.75,      # 시스템 효율 (75%, 경험적 값)
        
        # 버퍼 조건
        buffer_ph=7.4,                 # 버퍼 pH (PBS 기준)
        
        # 재현성
        random_seed=42                 # 재현성을 위한 시드 (같은 결과를 얻으려면 같은 값 사용)
    )
    
    # FASTA 파일 일괄 처리
    # 현재 디렉토리의 .fasta 파일을 자동으로 찾아서 분석합니다.
    analyzer.process_files(include_uncertainty=True)
    
    # 결과는 EPRM_Results_YYYYMMDD_HHMMSS 폴더에 저장됩니다.
    # - results.json: 모든 분석 결과 (JSON 형식)
    # - eprm_analysis_detail.log: 상세 로그
    # - config.json: 사용된 설정 파라미터

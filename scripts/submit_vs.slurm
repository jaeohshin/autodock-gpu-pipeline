#!/bin/bash
#SBATCH --job-name=vs_egfr             # Job 이름
#SBATCH --output=vs_%j.out            # 표준 출력 로그 (%j: job ID)
#SBATCH --error=vs_%j.err             # 표준 에러 로그
#SBATCH --partition=normal            # 사용할 파티션 이름
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16             # 사용할 CPU 개수
#SBATCH --time=120:00:00               # 최대 실행 시간

# 환경 설정
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

cd /store/jaeohshin/work/dock/scripts

# 특정 kinase에 대해서 실행 (예: abl1)
python run_vs.py --kinase egfr

# 또는 전체를 대상으로 실행하고 싶다면 아래 라인을 대신 사용:
#python run_vs.py --all


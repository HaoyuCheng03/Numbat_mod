#!/bin/bash
#SBATCH --mail-type=end
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=haoyuch3@andrew.cmu.edu
#SBATCH -e /ix1/alee/LO_LAB/Personal/Rick/Haoyu/numbat/simulation_tm/preprocess_mod/logs/%x.%j.e
#SBATCH -o /ix1/alee/LO_LAB/Personal/Rick/Haoyu/numbat/simulation_tm/preprocess_mod/logs/%x.%j.o
#SBATCH --account ctseng
#SBATCH --cpus-per-task=32
#SBATCH --nodelist=htc-n[32-49]
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH --mem=100g

echo LOG START $(date '+%Y-%m-%d %H:%M:%S')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 0: setting         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID="${sampleID:?}"
script=/ix1/alee/LO_LAB/Personal/Rick/Haoyu/numbat/simulation_tm/preprocess_mod/pileup_and_phase_mod.R

module load anaconda3/2023.09-0-python_3.11.5
source activate /ix1/alee/LO_LAB/Personal/Rick/env_conda  # for cellsnp-lite package

module load gcc/12.2.0 r/4.4.0
module load eagle2/2.4.1

echo PAP START $(date '+%Y-%m-%d %H:%M:%S')

/usr/bin/time -v Rscript "$script" --sampleID "$sampleID"


echo LOG END $(date '+%Y-%m-%d %H:%M:%S')

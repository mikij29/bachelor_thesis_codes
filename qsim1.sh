#!/bin/bash
#
#SBATCH --job-name=qsim1
#SBATCH --output=qsim1.out
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH -p express3
#SBATCH --array=0-199

readarray -t RANKS < ranks200.txt
RANK=${RANKS[$SLURM_ARRAY_TASK_ID]}
export RANK
	
module load R

Rscript --vanilla simulacni_studie.R

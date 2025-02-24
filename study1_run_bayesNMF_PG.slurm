#!/bin/bash

#SBATCH --job-name=study1_PG
#SBATCH --array=1-16
#SBATCH --time=20:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=2G       				
#SBATCH --error=logs/study1_PG.error
#SBATCH --output=logs/study1_PG.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript study1_run_bayesNMF_PG.R $SLURM_ARRAY_TASK_ID
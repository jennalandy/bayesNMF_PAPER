#!/bin/bash

#SBATCH --job-name=study2_SBFI
#SBATCH --array=1-20
#SBATCH --time=10:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=logs/study2_SBFI.error
#SBATCH --output=logs/study2_SBFI.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript study2_run_bayesNMF_PT_SBFI.R $SLURM_ARRAY_TASK_ID
#!/bin/bash

#SBATCH --job-name=study1_NT
#SBATCH --array=1-16
#SBATCH --time=2:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=1G       				
#SBATCH --error=logs/study1_NT.error
#SBATCH --output=logs/study1_NT.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript study1_run_bayesNMF_NT.R $SLURM_ARRAY_TASK_ID
#!/bin/bash

#SBATCH --job-name=pcawg_bayesNMF
#SBATCH --array=1-20
#SBATCH --time=20:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=2G       				
#SBATCH --error=logs/pcawg_bayesNMF.error
#SBATCH --output=logs/pcawg_bayesNMF.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript pcawg_run_bayesNMF.R $SLURM_ARRAY_TASK_ID
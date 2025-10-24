#!/bin/bash

#SBATCH --job-name=study2_PE_MH_SBFI
#SBATCH --array=1-20
#SBATCH --time=5:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=../../logs/study2/PE_MH_SBFI.error
#SBATCH --output=../../logs/study2/PE_MH_SBFI.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study2_PE_MH_SBFI.R $SLURM_ARRAY_TASK_ID
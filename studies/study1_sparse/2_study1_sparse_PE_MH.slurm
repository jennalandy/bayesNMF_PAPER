#!/bin/bash

#SBATCH --job-name=study1_sparse_PE_MH
#SBATCH --array=5-7

#SBATCH --time=2:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=../../logs/study1_sparse/PE_MH.error
#SBATCH --output=../../logs/study1_sparse/PE_MH.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study1_sparse_PE_MH.R $SLURM_ARRAY_TASK_ID
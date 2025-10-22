#!/bin/bash

#SBATCH --job-name=study1_sparse_PE
#SBATCH --array=5-7
#SBATCH --time=20:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=5G       				
#SBATCH --error=../../logs/study1_sparse/PE.error
#SBATCH --output=../../logs/study1_sparse/PE.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study1_sparse_PE.R $SLURM_ARRAY_TASK_ID
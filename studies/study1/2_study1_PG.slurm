#!/bin/bash

#SBATCH --job-name=study1_PG
#SBATCH --array=1-10
#SBATCH --time=20:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=5G       				
#SBATCH --error=../../logs/study1/PG.error
#SBATCH --output=../../logs/study1/PG.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study1_PG.R $SLURM_ARRAY_TASK_ID
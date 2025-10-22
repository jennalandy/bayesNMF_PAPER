#!/bin/bash

#SBATCH --job-name=study1_NE
#SBATCH --array=1-10
#SBATCH --time=2:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=../../logs/study1/NE.error
#SBATCH --output=../../logs/study1/NE.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study1_NE.R $SLURM_ARRAY_TASK_ID
#!/bin/bash

#SBATCH --job-name=study2_PE_MH_BFI
#SBATCH --array=1-20
#SBATCH --time=2:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=../../logs/study2/PE_MH_BFI.error
#SBATCH --output=../../logs/study2/PE_MH_BFI.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study2_PE_MH_BFI.R $SLURM_ARRAY_TASK_ID
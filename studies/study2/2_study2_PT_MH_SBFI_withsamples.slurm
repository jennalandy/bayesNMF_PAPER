#!/bin/bash

#SBATCH --job-name=study2_PT_MH_SBFI_withsamples
#SBATCH --array=1,5,10,15,20
#SBATCH --time=2:00:00       			
#SBATCH --ntasks=1
#SBATCH --mem=3G       				
#SBATCH --error=../../logs/study2/PT_MH_SBFI_withsamples.error
#SBATCH --output=../../logs/study2/PT_MH_SBFI_withsamples.out

module load gcc/9.2.0
module load R/4.3.2
module load openblas
Rscript 2_study2_PT_MH_SBFI_withsamples.R $SLURM_ARRAY_TASK_ID
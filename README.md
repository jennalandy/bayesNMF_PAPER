# bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures

This repository contains all code necessary to recreate all simulation studies and data application results presented in the paper ["bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures"](https://arxiv.org/abs/2502.18674). For the R package, see the [bayesNMF](https://github.com/jennalandy/bayesNMF/blob/master/README.md) repository.

Each section has three steps: data creation/processing, fitting models (running slurm jobs), and analyzing results. Below, we define all files needed to reproduce our results. This GitHub repostory demonstrates the expected directory structure.

## Simulation Study 1

- `study1_data.qmd`: simulates all datasets for study 1, stores data in `data/study1`
- Study 1 jobs, assumes `study1_data.qmd` has been run:
  - `study1_run_bayesNMF_<model>.R`: R script to run `bayesNMF` for a given model (indicated in file name) and N,G combination (indicated by command line argument corresponding to `study1_assignments.csv`), stores results in `output/study1`
  - `study1_run_bayesNMF_<model>.slurm`: slurm job array script that calls R script with all N,G combinations
- Study 1 analysis, assumes all jobs have completed:
  - `study1_analysis.qmd`: analyzes results of all study 1 experiments, creates Figure 1 and Supplementary Figure 1 for paper, loads in `study1_analysis.R`, creates figures in `figures/study1`

## Simulation Study 2

- `study2_data.qmd`: simulates all datasets for study 2, stores data in `data/study2`
- Study 2 jobs, assumes `study2_data.qmd` has been run:
  - `study2_run_bayesNMF_PT_<method>.R`:  R script to run `bayesNMF` for Poisson-Truncated Normal model with a given method to learn rank (indicated in file name) and N,G combination (indicated by command line argument corresponding to `study2_assignments.csv`), stores output in `output/study2`
  - `study2_run_bayesNMF_PT_<method>.slurm`: slurm job array script that calls R script with all N,G combinations
  - `study2_run_SignatureAnalyzer_all.sh`: bash script to run SignatureAnalyzer on all study 2 datasets, calls `study2_run_SignatureAnalyzer.sh`, stores results in `output/study2_SignatureAnalyzer`
- Study 2 analysis, assumes all jobs have completed:
  - `study2_analysis.qmd`: analyzes results of all study 2 experiments, creates Figure 2 for paper, loads in `study2_analysis.R`, stores figures in `figures/study2`

## PCAWG Data Application

- PCAWG data, must be run in the following order
  - `pcawg_data.qmd`: downloads data from UCSC Xena browser and processes it into SBS mutational counts matrices for each histology group, stores data in `data/PCAWG-UCSC`
  - `pcawg_defining_hypermutated.qmd`: defines hypermutation status for each sample in each histology group, saves `output/pcawg_defining_hypermutated/nonhyper_sample_names_updated.rds`
  - `pcawg_defining_hypermutated_EDA.qmd`: creates Supplementary Figure 2
- PCAWG jobs, assumes all `pcawg_data.qmd` and `pcawg_defining_hypermutated.qmd` have been run:
  - `pcawg_run_bayesNMF.R`: R script to run `bayesNMF` for Poisson-Truncated Normal model with SBFI and for a given histology group (indicated by command line argument corresponding to `pcawg_assignments.csv`), results saved in `output/pcawg_bayesNMF`
  - `pcawg_run_bayesNMF.slurm`: slurm job array script that calls R script with all histology groups
  - `pcawg_run_Signatureanalyzer_all.sh`: bash script to run SignatureAnalyzer on all histology groups, calls `pcawg_run_Signatureanalyzer.sh`, results saved in `output/pcawg_SignatureAnalyzer`
- PCAWG analysis, assumes all jobs have completed:
  - `pcawg_analysis.qmd`: creates Figures 3-4, loads in `pcawg_analysis_bayesNMF.R` and `pcawg_analysis_SignatureAnalyzer.R`, stores final figures in `figures/pcawg`

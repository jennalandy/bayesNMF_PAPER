
# bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures

This repository contains all code necessary to recreate all simulation studies and data application results presented in the paper ["bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures"](https://arxiv.org/abs/2502.18674). 

For the R package, see the [bayesNMF](https://github.com/jennalandy/bayesNMF) repository.

## Directory Structure
- `studies`: holds code to run all analyses (quarto notebooks, R scripts, and slurm job files, described in detail below)
- `data`: holds all raw simulated data
- `logs`: slurm job log files
- `output`: holds raw output from bayesNMF or SignatureAnalyzer
- `processed`: holds all processed data and processed output, incuding PCAWG mutational counts matrices matrices and simulation metrics
- `figures`: saved figures from all analyses

## Reproducible Analysis
All analyses must be run in a particular order: generating or processing data, running slurm jobs or bash scripts, processing output, and finally visualizing results. For the user's convenience, files are named in numerical order. Descriptions of each step are below. Any R scripts not numbered / not listed here contain helper functions imported by other files (don't worry about them!).

Core analyses (in paper):
- [`study1`](studies/study1): first set of simulation studies providing correct rank
    - `1_study1_data.qmd` simulates data for study 1
    - `2_study1_<model>.R` and `.slurm` runs bayesNMF with specified model on all simulated datsets in study 1
    - `3_study1_processing.R` loads and processes all study 1 results
    - `4_study1_results.qmd` creates visualizations of study 1 results
- [`study2`](studies/study2): second set of simulation studies, learning rank
    - `1_study2_data.qmd` simulates data for study 2
    - `2_study2_<model>.R` and `.slurm` runs bayesNMF with specified model on all simulated datsets in study 2
        - `PT_MH_SBFI_withsamples` reruns a subset of `PT_MH_SBFI` with saving all samples set to `TRUE` to showcase the label switching diagnostic plot. These are not used in analysis because logging all samples increases compute time.
    - `3_study2_run_SignatureAnalyzer_all_<prior>.sh` runs SignatureAnalyzer with a Poisson objective and L2 penalties on all simulated datsets in study 2
    - `4_study2_processing_SignatureAnalyzer.py` converts SignatureAnalyzer `.h5` output to csvs for processing in R
    - `5_study2_processing<variant>.R` processes the output of steps 3 and 4
    - `6_study2_results.qmd` creates visualizations of study 2 results
- [`pcawg`](studies/PCAWG/): pan cancer analysis of whole genomes (PCAWG) database analysis
    - `1_pcawg_data.qmd` to download and process PCAWG data
    - `2_pcawg_defining_hypermutated.qmd` to define which samples are hypermutated using a Negative Binomial mixture model
    - `3_pcawg_EDA.qmd` creates plots showing sample sizes of each histology group, number of mutations per sample, and identification of hypermutated samples
    - `4_pcawg_run_bayesNMF.R` and `.slurm` runs bayesNMF using the Poisson-Truncated Normal + MH model with SBFI on each histology group
    - `5_pcawg_run_SignatureAnalyzer_all.sh` runs SignatureAnalyzer with a Poisson objective and L2 penalties on each histology group
    - `6_pcwag_results.qmd` processes output from 4 and 5 and creates visualizations of results.
    - `7_capabilities.qmd` creates visualizations of R package capabilities (paper figure 1) using breast adenocarcinoma as an example dataset

Sensitivity analyses:
- [`study1_sparse`](studies/study1_sparse/): sensitivity analysis on study 1 for N = 4 only with 20-40% sparsity (see [this plot of sparsity](figures/study1_sparse/sparsity.png)). These simulations also have low mutation counts below 4 (see [this plot of median counts](figures/study1_sparse/median_count.png)). These analyses have a constant expected number of mutations per sample *per signature* of 100 mutations, so 400 mutations per sample (spread across 96 mutation types). This analysis only compares Poisson-Exponential and Poisson-Exponential + MH (focusing on the case of exact model alignment). With these results, we created a new version of main text Figure 2 comparing MAP point estimate [here](figures/study1_sparse/comparing_MAPs.png) as well as a new version of Appendix Figure C.2 comparing credible interval widths [here](figures/study1_sparse/compare_widths_combo.png) (See Appendix C for details). **This sensitivty analysis shows that even with sparse data, the Poisson-Exponential + MH model is still able to match the standard Poisson-Exponential model, both in terms of MAP estimates and posterior uncertainty.**
- [`study2_sparse`](studies/study2_sparse/): sensitivity analysis on study 2 with 10-75% sparsity (see [this plot of sparsity](figures/study2_sparse/sparsity.png)). These simulations also have low mutation counts below 15 (see [this plot of median counts](figures/study2_sparse/median_count.png)). These analyses have a constant expected number of mutations per sample *per signature* of 100 mutations, meaning N = 1 is far more sparse than N = 20. This was done to maintain a relatively constant power to discover each signature (as was done in the primary analysis study2 with expected 1000 per sample per signature). This analysis only compares bayesNMF using Poisson-Truncated Normal + MH with SBFI and SignatureAnalyzer with L2 priors. With these results, we created a new version of main text Figure 3 comparing rank accuracy, precision, and sensitivity between bayesNMF and SignatureAnalyzer [here](figures/study2_sparse/rank_metrics_aligned.png). **This sensitivity analysis shows that even with sparse data, Our Poisson-Truncated Normal + MH with SBFI performs comparably to SignatureAnalyzer with L2 priors. However, for large ranks, SBFI might induce too much sparsity causing precision and sensitivity to drop below SignatureAnalyzer.**
- [`study2_N10`](studies/study2_N10): sensitivity analysis to understand SBFI behavior when learning rank within 1:10 instead of 1:20. Comparing [study2_N10 results](figures/study2_N10/rank_metrics_aligned.png) to [regular study2 results](figures/study2/rank_metrics_aligned_P_T_MH_SBFI.png), we see the same general results (with a bit of noise) for ranks 1-10 in terms of rank bias, precision, and sensitivity. **This sensitivity analysis shows that decreasing the maximum rank does not bias the learned rank downwards, as long as the true rank is within the specified range.** For ranks over 10, bayesNMF estimates the specified maximum rank of 10. This shows that when there is enough signal in the data, bayesNMF is able to learn that maximum rank if appropriate, reassuring us that our model's underestimation of ranks 16-20 is due to a lack of signal in the data, not a downward bias from the specified maximum rank. Further, this leads us to the conclusion that **if the maximum rank is estimated, it is recommended to rerun the sampler with a higher maximum**.
- [`study2_N40`](studies/study2_N40): sensitivity analysis to understand SBFI behavior when learning rank within 1:40 instead of 1:20. Comparing [study2_N40 results](figures/study2_N40/rank_metrics_aligned.png) to [regular study2 results](figures/study2/rank_metrics_aligned_P_T_MH_SBFI.png), we see the same general results (with a bit of noise) in terms of rank bias, precision, and sensitivity. **This sensitivity analysis shows that increasing the maximum rank does not bias the learned rank upwards.** It also shows that our model's underestimation of ranks 16-20 is due to a lack of signal in the data, not a downward bias from the specified maximum rank.
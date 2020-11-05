# Wuxi
This project is used in collaboration with Wuxi NextCode, though it can be applied to more problems in general. The approach is to map a classification (binary or multinomial) to an Ising model such that it can be run on D-Wave. The mapping adapts a simple logistic regression. Collaborators used R, so code is a bit of a complicated mess of MATLAB, Python, and R :grimacing:.

## Description of files
The code is divided into folders based on the three languages. Roughly speaking:
- Matlab code is used to analyze solutions. Used because ported Code to Python was much slower. 
- R code is used to preprocess and output the data files, run classical comparisons, evaluate performance measures and plot
- Python code is mainly used to interface with the DW device. 

In other words, workflow was as follows: 
1. Get data from collaborators (usually feather files)
2. Preprocess data with R and output as .mat files
3. Run with DW using Python. Output solutions as .mat files
4. Combine solutions with Matlab. Print predictions of labels with matlab, saved as .mat file
5. Loaded predictions of labels with R. Calculate relevant metrics, generate plots as needed.

### R files
Mainly used for preprocessing and output data files, running classical classifiers, evaluating and plotting performance measures.

#### R/dataScripts
- output_serial_dilutions_pca_data.R: Loads data with all features, performs PCA, outputs serial dilutions. Outputs a separate .mat file for each cut of the data. Note: need to run another Matlab script to combine all .mat files for each cut into a single .mat file 
- output_bootstrap_resamples_data_for_DW.R: Analagous script for bootstrap resamples.
- output_lumAB_genes_for_DW.R: Output .mat files for lumAB using top 44 genes from PC1 (load preexisting gene-name file).
- output_6cancer_bootstrap_resamples_for_DW.R: Script for output 6 cancer data. Used because R is not great at memory management. 
- output_6cancer_serial_dilutions_pca.R: output .mat files for 6cancer with serial dilutions (incremental decrease) of data. Takes commandline arguments in order to have better memory management. 
- Notes: 
    - basically all files can be run with commandline arguments. A simple wrapper script can be written to run on the cluster, embarrassingly parallel.
    - to save time for both the conventional and annealing-based approaches, .mat files are used as the inputs with all splits contained within. This seemed preferable to loading the entire data matrix each time, performing PCA, and then moving forward.

#### R/runScripts
- \{6cancer,all\}\_bootstrap_resamples.R: Runs conventional approaches on \{6cancer, all binomial\} datasets. Outputs .RDS files. 
- \{6cancer,\}\_serial_dilutions_cl.R: Runs conventional approaches on \{6cancer, binomial\} fractional decrease of training datasets. Outputs .RDS files
- lumAB_diffexp_bootstrap_resamples.R: Runs conventional approaches on genes from differential expression. Outputs .RDS files.
- \{6cancer\_,\}serial_dilutions_RBM.R: Runs RBM on \{6cancer, all binomial\} fractional decrease of training datasets. Outputs .RDS files 
- \{6cancer\_,\}RBM\_Script.R: Runs RBM on \{6cancer, all binomial\} datasets. Outputs .RDS files.
- Notes: all files can be run with commandline arguments to facilitate running on a cluster. Rather than loading the splits, data is precomputed and saved in .mat files. 

#### R/analysisScripts
- \{6cancer,\}output_bootstrap_resamples_rds.R: From DW (or SA, Rand, Field) .mat files with actual and predicted labels generates .rds file for \{6cancer, binomial\}bootstrap resamples.
- output\_\{6cancer,\}_serial\_dilutions_rds.R: From DW (or SA, Rand Field). mat files with actual and predicted labels generates .rds file for \{6cancer, binomial\} serial dilutions.
- \{bootstrap_resamples\_,lumAB_serial_dilutions\}\_plots.R:
- plot_bootstrap_resamples_stacked.R: Generates stacked plot for binomial bootstrap resamples. Data files for both conventional and annealing-based methods are saved in output/  
- plot_bootstrap_resamples\_\{6cancer,lumAB_gene\}.R: Generates for \{6cancer, lumAB gene\} bootstrap resamples. Data files for both conventional and annealing-based methods are saved in output/  
- plot_serial_dilutions.R: Generates plots for incremental decrease of training data. Can be used for binomial or 6cancer data.
- wilcox_rank_sum_v2.R: Takes the metrics from different splits of data and calculates Wilcoxon rank sum test to generate p-value of hypothesis that two populations are from the same distribution.
- Notes:
    - These scripts can be used after running conventional ML (R/runScripts) and annealing-based (python/runScripts then matlab/analysisScripts) methods
    - annealing based approaches generate .mat files with predictions. \*output\* compiles the predictions into an .RDS file that has the same format as the conventional ML methods
    - plot\*.R generates plots used in the paper

### Python files
Mainly used for running DW

- DW_utils.py: List of methods for running on DW. Includes some functions to automatically do cross-validation. Also includes functions to generate heuristic and balanced embeddings. 
- analyze_utils.py: Analysis methods for solutions. Not used that much because slower than Matlab code.

#### python/runScripts
Fairly self-explanatory: run bootstrap and lumAB serial dilutions. 

### Matlab files
Mainly used for combining DW (and SA) output solutions.

- analyzeLogisticResults.m: Main function to analyze results when using logistic regression as model. 
- analyzeMultinomialResults.m: Main function to analyze results when using multinomial regression as model. Less updated than analyzeLogisticResults.m
- printSACV.m: Function to print instance files for SA.
- getSACVSols.m: Function to read in solutions and save .mat file after running SA.
- the other files are helper functions that are called in other scripts. 

The rest of files are helper functions that might occasionally want to call separately. Probably could be further condensed. 

#### matlab/dataScripts
- wrap_print_SA_\{bootstrap_resamples,serial_dilutions\}.m: Script to provide parameters to print SA files for \{bootstrap resamples, fractional decrease of training data\}. 
- combine_mat_files.m: Script to combine .mat files output from R into a single .mat file 

#### matlab/analysisScripts
- analyze_bootstrap_resamples\_\{simple,SA,rand,field\}.m: Simplified version of analysis script \{DW,SA,Rand,Field\} for bootstrap resamples with select parameters for post_processing.  
- analyze_serial_dilutions\_\{simple,SA,rand,field\}.m: Analysis script for incremental decrease of training data. Uses a predefined set of parameters. 
- analyze_6cancer\_\{bootstrap_resamples,serial_dilutions\}\_\{simple,SA,rand,field\}.m Analysis script for annealing-based methods of 6cancer \{bootstrap resamples, incremental decrease of training data\}. 
- compare_bootstrap_resamples_energies.m: Compare DW and SA energy distributions for bootstrap resamples.
- compare_lumAB_serial_dilutions_energies.m: Compare DW and SA energy distributions for seiral dilutions.
- output_serial_dilutions_preds_for_R.m: Simple wrapper function to output .mat files for predictions for incremental decrease of training data.

#### matlab/SAScripts/
- get_SA_\{bootstrap_resamples,serial_dilutions\}\_sols.m: Wrapper function to read SA output files. 
- the other files are helper functions for reading and printing SA files 

### SA
- https://arxiv.org/abs/1401.1084 for paper and original source code
- Code to run simulated annealing (as a zip). The executable that should be used for complete graph problems is an_ss_ge_fi_vdeg. 
    - To generate the binary with multiple threads, type `make an_ss_ge_fi_vdeg_omp` (omp refers to multi-threading) in the directory with the makefile
    - See README.txt for list of parameters to use
    - See arXiv reference for other versions of SA (e.g., for Chimera)
- runNoReg.sh: script to run SA once the binary is compiled. Includes some default parameters. 

## Dependencies (python, use pip unless otherwise indicated):
- minorminer
- sklearn 
- numpy
- scipy
- chimera_embedding (github.com/USCqserver)
- dwave_sapi (download api from dw website)

## Examples 
Possible workflow:
- Receive file from collaborators containing data
- Run lumAB_output_serial_dilutions_pca.R (set save_directory)
- Run combine_mat_files.m
- Run lumAB_serial_dilutions.py
- (Concurrently, run lumAB_serial_dilutions_cl.R)
- Run run_lumAB_serial_dilutions.py
- Run lumAB_serial_dilutions_output_dw_rds.R
- Run bootstrap_resamples_plots.R

## Authors
* Richard Li 

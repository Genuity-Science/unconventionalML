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
- lumAB_output_serial_dilutions.pca.R: Loads data with all features, performs PCA, outputs serial dilutions. Outputs a separate .mat file for each cut of the data. Note: need to run another Matlab script to combine all .mat files for each cut into a single .mat file 
- output_bootstrap_resamples_data_for_DW.R: Analagous script for bootstrap resamples.
- output_6cancer_bootstrap_resamples_for_DW.R: Script for output 6 cancer data. Used because R is pretty terrible at memory management. 

#### R/runScripts
- lumAB_serial_dilutions.R: Wuxi version of code. Explicitly writes out classical approaches, save mean and std in .txt files
- lumAB_serial_dilutions_cl.R and lumAB_serial_dilutions_cl2.R: Version of code. Writes out classical approaches in loops. Saves .RDS file with more metrics. Consider adding a few lines (taken from lumAB_serial_dilutions.R) to output mean and std in .txt files. Split into two scripts for running on hpc
- all_bootstrap_resamples.R: Runs binomial bootstrap resamples. 

#### R/analysisScripts
- bootstrap_resamples_output_\{dw,sa\}_rds.R: From DW (or SA) .mat files with actual and predicted labels generates .rds file for bootstrap resamples.
- lumAB_serial_dilutions_output_\{dw,sa\}_rds.R: From DW (or SA). mat files with actual and predicted labels generates .rds file for lumAB serial dilutions.
- \{bootstrap_resamples_,lumAB_serial_dilutions\}_plots.R: Generates plots for bootstrap resamples or (lumAB serial dilutions). Classical data is in mean and std text files, SA and DW files are saved as .rds files. Should standardize output format.

#### R/archive
Contains old scripts used to generate data outputs, plots. 

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

The rest of files are helper functions that might occasionally want to call separately. Probably could be further condensed. 

#### matlab/dataScripts
- wrap_print_SA_bootstrap_resamples.m: Script to provide parameters to print SA files
- combine_mat_files.m: Script to combine .mat files output from R into a single .mat file 

#### matlab/analysisScripts
- analyze_bootstrap_resamples.m: Simplified version of analysis script for bootstrap resamples. There are additional parameters that could choose to tune. See analyze_lumAB_serial_dilutions.m
- analyze_lumAB_serial_dilutions.m: Analysis script for bootstrap resamples. Does a grid search over possible parameters to improve performance.
- analyze_lumAB_serial_dilutions_simple.m: Analysis script for bootstrap resamples. Uses a predefined set of parameters. 
- analyze_SA_bootstrap_resamples.m: Analysis script for bootstrap resamples with SA.
- analyze_SA_lumAB_serial_dilutions.m: Analysis script for lumAB serial dilutions with SA.
- compare_bootstrap_resamples_energies.m: Compare DW and SA energy distributions for bootstrap resamples.
- compare_lumAB_serial_dilutions_energies.m: Compare DW and SA energy distributions for seiral dilutions.
- get_SA_\{bootstrap_resamples,lumAB_serial_dilutions\}_sols.m: Read in SA solutions for bootstrap resamples and lumAB serial dilutions
- output_\{bootstrap_resamples,lumAB_serial_dilutions\}_preds_for_R.m: output .mat file with predicted and actual labels for analysis with R. Should output the files in analyze\*.m, but in case want to output predictions from another results file
- wrap_analyze_multinomial_solutions.m: Wrapper script for analyzing multinomial solutions. Might be slightly out of date.

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

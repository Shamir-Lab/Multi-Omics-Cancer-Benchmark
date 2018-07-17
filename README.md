# Multi-Omics-Cancer-Benchmark

This repository contains the code that was used for the benchmark on TCGA data.
All preprocessed datasets, survival data and clinical labels that were used in the analysis are available here: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.
To run the benchmark, some configuration is needed, mainly paths to the datasets, survival and clinical data. Additionally, paths to binaries are required for some methods (MultiNMF, rMKL-LPP).
The expected directory structure is: for the clinical data, a single directory with all clinical data files.
The datasets and survival data are organized as follows: a single root directory, with subdirectories for every cancer type. The directory for each cancer type contains a file for every omic, and a file called "survival" with survival data.
Omic files, survival data files and clinical label files should be formatted as the files in http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.

# Population variation in the immune response to influenza
code necessary to reproduce the main analyses published in:

Genetic ancestry effects on the response to viral infection are pervasive but cell type specific (Science, 2021)

HE Randolph, JK Fiege, BK Thielen, CK Mickelson, M Shiratori, J Barroso-Batista, RA Langlois, and LB Barreiro

# General dependencies
R (tested in version 3.6.3)

gcc or LLMV (tested with Apple LLVM version 8.0.0 (clang-800.0.42.1))

# Usage
1. place the folder `inputs`, available for download on Zenodo (10.5281/zenodo.4273999), in the desired working directory
2. place the folders `main_analyses` and `common_functions` in the same working directory
3. check individual dependencies at the header of each script stored in `main_analyses`, and install the necessary CRAN and Bioconductor packages
4. change the working directory in the script to your desired working directory
```
current = [desired working directory]
```
5. run the scripts stored in `main_analyses` in R -- if `main_analyses` and `inputs` are located in the same directory, results will populate in a folder named `outputs`

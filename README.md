positive-unlabeled learning for dms datasets (pudms)
================

## Description

This package offers a streamlined analysis via PUlasso algorithm for
learning sequence-function relationships using deep mutational scanning
data sets.

## Installation

This step is only needed once on a computer where the `pudms` package has not
been installed before. Start the `R` interpreter and run the lines below to
install from Github using the the **devtools** package:

```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("RomeroLab/pudms")
quit(save="no")
```

## Input files

Two protein sequence text files, one for the labeled (positive) and the other for the unlabeled (unselected).
    
- Text files with one sequence per line
- (Optional) Can be gzip compressed.
- Example
    
    
    ```console
    MYYKEIAHALFSDLFALSELYIAVRY*
    MYYKEIAHALFSALFALSPLYIAVRM*
    MYYKELAHALFSAHFALSELYIAVRY*
    MYYKEIAHALFSALFAEHELYIAVRY*
    MRYKEIAHALFSALFALPELYIAVRY*
    ...
    ```


-----

## Example

To setup this example copy all the following files in the [quickexample
directory](./inst/quickexample) to your working directory. 

This [example script](./inst/quickexample/fit_PU_model.R) demonstrates basic
usage of the package using sample labeled and unlabeled sequences.

``` r
# Run this script from the command line 
# R --vanilla --no-save < fit_PU_model.R
remove(list=ls())
# LOAD THE PUDMS LIBRARY
library(pudms)


# SET THE POSITIVE AND UNLABELED SEQUENCE FILES (can be gzipped or not)
pos_file = 'Rocker_sel_sequences_filtered.txt.gz'
unlab_file = 'Rocker_ref_sequences_filtered.txt.gz'


# VARIOUS RUN OPTIONS 
py = NULL         # Proportion of positive sequences in unlabeled set (i.e. fraction functional).
                  # NULL scans a range of possible py values between 1e-3 and 0.5
order = 1         # Model order: 1 for main effects or 2 for pairwise
refstate = NULL   # Reference state for regression.  
                  # NULL chooses the most common residue at each position (preferable for DMS data).  
                  # In contrast, chimera data should use a fixed reference (e.g. 'A')
nobs_thresh = 10  # Filters out columns in X that sum to less than nobs_thresh
n_eff_prop = 1    # Scales the p-values to account for redundant sequence sampling at the NGS step. 
                  # See more in note below.


# OUTPUT FILES
outroc = 'Rocker_CV_ROC.png'
outcsv = 'Rocker_parameters.csv'


# CREATE A PROTEIN DATA SET
pudata = create_protein_dat(path_l = pos_file, path_u = unlab_file) 


# PERFORM CROSS-VALIDATED FITTING OF PU MODEL
cvfit = v.pudms(protein_dat = pudata,
                py1 = py,
		order = order,
		refstate = refstate,
		nobs_thresh = nobs_thresh,
		n_eff_prop = n_eff_prop,
		nhyperparam = 10, # The number of py values to scan. Log spaced between 1e-3 and 0.5
		nfolds = 5,       # The number of cross-validation folds
		nCores = 10)      # The number of threads to use for CV.  


# PLOT THE PU-CORRECTED ROC CURVE FOR THE CV FIT
rocplot = with(cvfit, rocplot(roc_curve = roc_curves[[which(py1 == py1.opt)]], py1 = py1.opt))
ggsave(filename = outroc, plot = rocplot)


# REFIT ALL THE DATA WITH THE OPTIMAL PY VALUE AND WRITE MODEL PARAMETERS/PVALUES TO CSV
optpy = cvfit$py1.opt
cat("The optimal py value is", optpy, "\nRefitting model on all the data with this py value\n")
fit = pudms(protein_dat = pudata, 
            py1 = optpy,
	    order = order,
	    refstate = refstate,
	    nobs_thresh = nobs_thresh,
	    n_eff_prop = n_eff_prop,
	    outfile = outcsv) 
```


### Run example
Run the edited script on your data. This will save a csv file with the output in your working directory. 
```shell
R --vanilla --no-save < fit_PU_model.R
```

## Output 

This example's curve looks like 
<p align="center">
<img src="inst/quickexample/Rocker_CV_ROC.png" width="400" title="CV ROC " />
</p>

The first few rows of the results look like

|   |coef               |se                |zvalue           |p                   |p.adj               |nobs|eff_nobs|sep|group|chi2value       |p.grp|p.grp.adj|dat.g.nobs.grp|
|------|-------------------|------------------|-----------------|--------------------|--------------------|----|--------|---|-----|----------------|-----|---------|--------------|
|Y0.*  |2.15471956767264   |0.0315308746734723|68.3368155811251 |0                   |0                   |7550|7550    |*  |0    |1307680.35553662|0    |0        |192550        |
|Y0.A  |-0.0813983989127809|0.0262583447608329|-3.09990594053725|0.00193582106653184 |0.00216512105618937 |8129|8129    |*  |0    |1307680.35553662|0    |0        |192550        |
|Y0.C  |-1.40667536321285  |0.0902292785710968|-15.5900100886261|8.51201104105626e-55|1.74870184985716e-54|1507|1507    |*  |0    |1307680.35553662|0    |0        |192550        |



The entire example results are available [here](./inst/quickexample/Rocker_parameters.csv)
 

positive-unlabeled learning for dms datasets (pudms)
================

# Description

This package offers a streamlined analysis via PUlasso algorithm for
learning sequence-function relationships using deep mutational scanning
data sets.

# Installation

Install using **devtools** package:

``` r
# install.packages("devtools")
devtools::install_github("RomeroLab/pudms")
```

# Input files

Two protein sequence text files, one for the labeled (positive) and the other for the unlabeled (unselected).
    
      - Text files with one sequence per line
      - (Optional) Can be gzip compressed.
      - Example
    
    <!-- end list -->
    
    ``` console
    MYYKEIAHALFSDLFALSELYIAVRY*
    MYYKEIAHALFSALFALSPLYIAVRM*
    MYYKELAHALFSAHFALSELYIAVRY*
    MYYKEIAHALFSALFAEHELYIAVRY*
    MRYKEIAHALFSALFALPELYIAVRY*
    ...
    ```


-----

# Example

This example script (Link to fit_PU_model.R) demonstrates basic usage of the package using sample labeled and unlabeled sequences.

``` r
# Run this script on the group server using the command
# R --vanilla --no-save < fit_PU_model.R

# LOAD THE PUDMS LIBRARY
library(pudms)


# SET THE POSITIVE AND UNLABELED SEQUENCE FILES (can be gzipped)
pos_file = 'Rocker_sel_sequences_filtered.txt'
unlab_file = 'Rocker_ref_sequences_filtered.txt'


# VARIOUS RUN OPTIONS 
py = NULL         # Proportion of positive sequences in unlabeled set (i.e. fraction functional).  NULL scans a range of possible py values between 1e-3 and 0.5
order = 1         # Model order: 1 for main effects or 2 for pairwise
refstate = NULL   # Reference state for regression.  NULL chooses the most common residue at each position (preferable for DMS data).  In contrast, chimera data should use a fixed reference (e.g. 'A')
nobs_thresh = 10  # Filters out columns in X that sum to less than nobs_thresh
n_eff_prop = 1    # Scales the p-values to account for redundant sequence sampling at the NGS step. See more in note below.


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
		nhyperparam = 3, # The number of py values to scan. Log spaced between 1e-3 and 0.5
		nfolds = 4,      # The number of cross-validation folds
		nCores = 1,      # The number of cores to use for CV.  Q: nCores > 1 causes this to crash with the error: 2 nodes produced errors; first error: object 'refstate' not found; Calls: v.pudms ... clusterApply -> staticClusterApply -> checkForRemoteErrors
		full.fit = FALSE) # Q: Hyebin going to make this default 


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
	    outfile = outcsv) # Q: this does not return group pvalues

```




# Run
Run the edited script on your data. This will save a csv file with the output in your working directory. 
```shell
R --vanilla --no-save < fit_model_roc.R
```

## Output 

This example's curve looks like 
![ROC curve](roc_curve.png)



The first few rows of the results look like

|    |coef              |se               |zvalue            |p                 |p.adj            |nobs|
|------|------------------|-----------------|------------------|------------------|-----------------|----|
|S1.P  |-1.37592271951636 |0.737876703960768|-1.86470546113016 |0.0622227008752244|0.292995889618509|12  |
|I3.M  |0.0150377890146403|0.72996792761058 |0.0206006160624945|0.98356424902434  |0.991380653652348|13  |
|I3.V  |0.857480027936806 |0.966663090177415|0.887051586690281 |0.375051127465389 |0.691074837775208|16  |
|S4.P  |-0.566703453924022|0.759689634851526|-0.745967073823218|0.455687305312763 |0.74104242419004 |10  |
|L5.P  |0.275912251075667 |0.748047109029646|0.368843416069846 |0.71224444142229  |0.884792247202554|15  |

The entire example results are available [here](./results.csv). 
 

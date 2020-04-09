# Run this script from the command line 
# R --vanilla --no-save < fit_PU_model.R
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



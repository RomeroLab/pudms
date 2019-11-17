## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(pudms)
# recover paths
path_l = system.file("extdata", "sample_positive_mutations.txt", 
                     package = "pudms",mustWork = T)
path_u = system.file("extdata", "sample_unlabeled_mutations.txt", 
                     package = "pudms",mustWork = T)

# WT (reference) sequence
WT = 'MSGISLDNSYKMDYPEMGLCIIINNKNFHKSTGMTSRSGTDVDAANLRETFRNLKYEVRNKNDLTREEIVELMRDVSKEDHSKRSSFVCVLLSHGEEGIIFGTNGPVDLKKITNFFRGDRCRSLTGKPKLFIIQACRGTELDCGIETDSGVDDDMACHKIPVEADFLYAYSTAPGYYSWRNSKDGSWFIQSLCAMLKQYADKLEFMHILTRVNRKVATEFESFSFDATFHAKKQIPCIVSMLTKELYFYHLEHHHHHH*'

## ------------------------------------------------------------------------
smpl_dat = create_protein_dat(path_l = path_l, path_u = path_u, WT = WT)

## ------------------------------------------------------------------------
# we filter sequences containing mutations with total number < nobs_thresh
fit1 = pudms(protein_dat = smpl_dat,py1 = 0.2,nobs_thresh = 10) 

## ----eval=F--------------------------------------------------------------
#  View(fit1$result_table)

## ----eval=F--------------------------------------------------------------
#  require(xlsx)
#  write.xlsx(x = fit1$result_table, sheetName = "result", file = "your path/result_table.xlsx")

## ------------------------------------------------------------------------
cv_smpl = cv_grouped_dat(grouped_dat = smpl_dat,
                         test_idx = 1,
                         nfolds = 10,
                         seed = 1)

tr_smpl = cv_smpl$train_grouped_dat # training set
tt_smpl = cv_smpl$test_grouped_dat # test set

## ------------------------------------------------------------------------
cv.r = pudms(protein_dat = tr_smpl,py1 = 0.2,nobs_thresh = 10)

## ----fig.height=6,fig.width=6--------------------------------------------
roc = corrected_roc_curve(coef = coef(cv.r), test_grouped_dat = tt_smpl,py1 = 0.2)
roc$rocplot


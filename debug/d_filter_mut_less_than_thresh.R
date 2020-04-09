## debug file for filter_mut_less_than_thresh
remove(list=ls())
library(pudms)
data('rocker')

verbose = T
protein_dat = rocker
refstate = NULL
order= 2

if(verbose) cat(" 1. create model frames from an aggregated dataset:\n")
Xprotein = create_model_frame(grouped_dat = protein_dat[1:10000,],
                              order = order, 
                              aggregate = T,
                              refstate = refstate,
                              verbose = verbose)
library(magrittr)

Xprotein$X %>% dim
Xprotein$blockidx %>% length


p0 = 600
nobs_thresh = 10
Xprotein$X = Xprotein$X[,1:p0]
Xprotein$blockidx  = Xprotein$blockidx[1:p0]


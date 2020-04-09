library(magrittr)
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

Xprotein$X = Xprotein$X[,1:p0]
Xprotein$blockidx  = Xprotein$blockidx[1:p0]


p0 = 600
nobs_thresh = 10

filtered_Xprotein = filter_mut_less_than_thresh(Xprotein = Xprotein,
                                                order = order,
                                                nobs_thresh = nobs_thresh,
                                                checkResFullRank = T,
                                                verbose = verbose)


with(filtered_Xprotein, Diagonal(x=wei)%*%X) %>% colSums

lambda =0
initial_coef = NULL
py1 = 0.001
nlambda = 10
maxit = 2
eps = 1e-2
inner_eps = 1e-1
if(verbose) cat("\n\n 2. fit a model\n")
fit = grpPUlasso(X = filtered_Xprotein$X,
                 z = filtered_Xprotein$z,
                 weights = filtered_Xprotein$wei,
                 py1 = py1,
                 verbose = verbose,
                 lambda = lambda, 
                 nlambda = nlambda,
                 maxit = maxit,
                 eps = eps,
                 inner_eps = inner_eps,
                 initial_coef = initial_coef
)


pvalue = T
n_eff_prop  = 1
if(pvalue){
  if(verbose) cat("\n\n 3. compute p-values\n")
  if(lambda>0){stop("currently p-value computation is supported only for lambda == 0")}
  if(n_eff_prop>1 || n_eff_prop <= 0){stop ("n_eff_prop should be between 0 and 1")}
  
  pvalues <-
    pval_pu(
      X = filtered_Xprotein$X,
      z = filtered_Xprotein$z,
      theta = as.numeric(fit$coef),
      py1 = py1,
      weights = filtered_Xprotein$wei,
      effective_n_prop = n_eff_prop
    )
  
}else{
  pvalues = NULL
}

nobs = c(0,colSums(Diagonal(x=filtered_Xprotein$wei)%*%filtered_Xprotein$X))

Xprotein = filtered_Xprotein
p.adjust.method = "BH"




return_table = return_tables(Xprotein = filtered_Xprotein, fit=fit,pvalues = pvalues,nobs = nobs,n_eff_prop = n_eff_prop,p.adjust.method = p.adjust.method)

library(Matrix)
# export results to the excel file



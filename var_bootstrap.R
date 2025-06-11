# VAR_BOOTSTRAP
# function for calculating variance matrices by bootstrap
# bootstrap median and bootstrap mean estimators (as defined in attachment A.1 of the thesis) included

# INPUT: 
# - regressors (without the intercept), response; 
# - R - the number of bootstrap realizations
# - which estimators should be used (T/F for e_lws, e_ls, e_mm, e_MM, e_MMa, e_s, e_lts, e_lms, e_rdl1)
# - if LWS, which parameters should be used (lws_1, lws_2, lws_3, lws_4, lws_5, dd_lws, dd_wls, eps_w, J, iniest, ninit)
# - boot_med, boot_mean - T/F, whether should we estimate beta by the bootstrap mean and/or the median estimator based on R bootstrap realizations

# OUTPUT:
# - beta - an array consisting of (boot_med + boot_mean) matrices (possibly for bootstrap median and for bootstrap mean estimator) of the estimated coefficients (in each column of these matrices there is one estimator of beta)
# - resids - an array consisting of (boot_med + boot_mean) matrices (possibly for bootstrap median and for bootstrap mean estimator) of the model residuals (in each column of these matrices there is one vector of residuals)
# - ecmx - empirical covariance matrices for all the individual estimators, i.e. noest matrices, where noest is the total number of estimators; each covariance matrix is computed from R bootstrap realizations
# - esttimes - the matrix of times we spent computing the quantities for a given estimator (in the rows; the columns correspond to the output of the function system.time)


source("load_names.R")
source("just_estimate.R")

var_bootstrap = function(regressors = NA, response, R, e_lws = TRUE, lws_1 = TRUE, lws_2 = TRUE, lws_3 = TRUE, lws_4 = TRUE, lws_5 = TRUE, dd_lws = TRUE, dd_wls = TRUE, eps_w = 1e-33, J = 10000, iniest = "S", ninit = 555, e_ls = FALSE, e_mm = FALSE, e_MM = F, e_MMa = F, e_s = FALSE, e_lts = FALSE, e_lms = FALSE, e_rdl1 = FALSE, boot_med = TRUE, boot_mean = TRUE) {
  
  loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
  noest = loaded_estnames$number
  estnames = loaded_estnames$names
  
  loaded_bootnames = load_bootnames(boot = R, boot_med = boot_med, boot_mean = boot_mean)
  noboot = loaded_bootnames$number
  bootnames = loaded_bootnames$names
  
  icpt_only = anyNA(regressors)
  
  if(icpt_only)
    p = 1
  else {
    p = NCOL(regressors)+1
    if(!is.data.frame(regressors))
      regressors = data.frame(regressors)
  }
  n = NROW(response)
  
  gamma = array(dim = c(R, p, noest))
  esttimes = array(0, dim = c(noest, 5))
  dimnames(esttimes)[[1]] = estnames
  # dimnames(gamma)[[3]] = estnames
  if(noboot > 0) {
    beta = array(dim = c(p, noest, noboot))
    dimnames(beta)[[2]] = estnames
    dimnames(beta)[[3]] = bootnames
  }
  else
    beta = NA
  for (i in 1:R) {
    rows = sample(1:n, size = n, replace = TRUE)
    
    resp_new = response[rows]
    if(icpt_only)
      regr_new = NA
    else
      regr_new = regressors[rows,]
    
    just_estimated = just_estimate(regressors = regr_new, response = resp_new, e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, eps_w = eps_w, J = J_boot, iniest = iniest, ninit = ninit, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
    gamma[i,,] = just_estimated$coef
    esttimes = esttimes + just_estimated$time
  }
  
  ecmx = apply(gamma, 3, cov)
  
  i = 1
  if(boot_med) {
    beta[,,i] = apply(gamma, c(2, 3), median)
    i = i + 1
  }
  if(boot_mean) {
    beta[,,i] = apply(gamma, c(2, 3), mean)
    i = i + 1
  }
  # next idea: in.which lowest mse

  ret_res = function(beta_est){
    if(icpt_only)
      return(response-beta_est)
    else
      return(response-as.matrix(cbind(1, regressors))%*%beta_est)
  }
  
  if(noboot > 0)
    resids = apply(beta, c(2, 3), ret_res)
  else
    resids = NA
  
  return(list(beta, resids, ecmx, esttimes))
}


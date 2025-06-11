# JUST_ESTIMATE
# calculating estimates via various estimators, without any estimation of their variances

# INPUT: 
# - regressors (without the intercept), response; 
# - which estimators should be used (T/F for e_lws, e_ls, e_mm, e_MM, e_MMa, e_s, e_lts, e_lms, e_rdl1)
# - if LWS, which parameters should be used (lws_1, lws_2, lws_3, lws_4, lws_5, dd_lws, dd_wls, eps_w, J, iniest, ninit)

# OUTPUT:
# - coef - the matrix of estimated coefficients (in each column there is one vector of estimated coefficients; the number of columns is equal to the number of estimators used)
# - resid - the matrix of model residuals (in each column there is one vector of residuals)
# - time - the matrix of times it took to calculate the given estimator (in the rows; the columns correspond to the output of the function system.time)


source("lws.R")
source("load_names.R")
library("MASS")
library("robust")
library("robustbase")

just_estimate = function(response, regressors = NA, e_lws = TRUE, lws_1 = TRUE, lws_2 = TRUE, lws_3 = TRUE, lws_4 = TRUE, lws_5 = TRUE, dd_lws = TRUE, dd_wls = TRUE, eps_w = 1e-33, J = 10000, iniest = "S", ninit = 555, e_ls = FALSE, e_mm = FALSE, e_MM = FALSE, e_MMa = FALSE, e_s = FALSE, e_lts = FALSE, e_lms = FALSE, e_rdl1 = FALSE) {
  loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
  noest = loaded_estnames$number
  estnames = loaded_estnames$names
  lws_names = loaded_estnames$lws_names
  
  icpt_only = anyNA(regressors)
  
  if(icpt_only) 
    p = 1
  else {
    p = NCOL(regressors) + 1
    if(!is.data.frame(regressors))
      regressors = data.frame(regressors)
  }
  n = length(response)
  beta = array(0, c(p, noest))
  resid = array(0, c(n, noest))
  esttimes = array(0, dim = c(noest, 5))
  dimnames(esttimes)[[1]] = estnames
  # dimnames(esttimes)[[2]] = c("user", "system", "elapsed")
  
  k = 1
  if(e_lws) {
    for(lws_name in lws_names) {
      esttimes[k, ] = esttimes[k, ] + system.time({
        model = lws(regressors = regressors, response = response, psi_type = lws_name, eps_w = eps_w, J = J, iniest = iniest, ninit = ninit)
        beta[, k] = model$coefficients
        resid[, k] = model$residuals })
      k = k + 1
    }
  }
  if(e_ls) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lm(response~., data = regressors)
      else  
        model = lm(response~1)
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_mm) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lmrob(response~., data = regressors)
      else
        model = lmrob(response~1)
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_MM) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lmRob(response~., data = regressors)
      else
        model = lmRob(response~1)
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_MMa) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lmRob(response~., data = regressors, control = lmRob.control(final.alg = "Adaptive"))
      else
        model = lmRob(response~1, control = lmRob.control(final.alg = "Adaptive"))
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_s) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lqs(x = regressors, y = response, method = "S")
      else
        model = lqs(response ~ 1, method = "S")
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_lts) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lqs(x = regressors, y = response, method = "lts")
      else
        model = lqs(response ~ 1, method = "lts")
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_lms) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lqs(response~., data = regressors, method = "lms")
      else
        model = lqs(response~1, method = "lms")
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
  }
  if(e_rdl1) {
    esttimes[k, ] = esttimes[k, ] + system.time({
      if(!icpt_only)
        model = lm(response~., data = regressors)
      else
        model = lm(response~1)
      beta[, k] = model$coefficients       
      resid[, k] = model$residuals })
    k = k + 1
    # not yet implemented...
    print("This one should have been turned off.")
  }
  
  return(list("coef" = beta, "resid" = resid, "time" = esttimes))
}


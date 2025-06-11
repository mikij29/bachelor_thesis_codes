# PSI_DEP
# function for creating data-dependent weights

# INPUT: 
# - regressors, response;
# - iniest - the initial estimator - one of the set {"ltsMy", "ls", "S", "lts", "lms"}
# - ninit - the number of iterations for the initial estimator computation - by default 555
# - eps_w - the individual weights lower than this constant are replaced by this constant (in order not to divide by zero); by default 1e-33

# OUTPUT:
# - data-dependent weights


library("MASS")

psi_dep = function(response, regressors = NA, iniest = "S", ninit = 555, eps_w = 1e-33){
  # both forms of regressors are acceptable: data.frame and matrix
  
  icpt_only = anyNA(regressors)
  
  if(icpt_only)
    p = 1
  else
    p = NCOL(regressors)+1
  n = NROW(response)
  
  h = floor(n/2) + floor((p+1)/2)
  # one might choose his own h - possible improvement...
  # if (missing(h) || (h <= ceiling(n/2)) || (n < h))...

  # the initial estimate - we are interested only in the residuals
  if(icpt_only) {
    if (iniest == "ltsMy")
      init.res = lws(response = response, psi_type = 5, J = ninit)$residuals
    else if (iniest == "ls")
      init.res = lm(response ~ 1)$residuals
    else if (iniest == "S")
      init.res = lqs(response ~ 1, method = "S", nsamp = ninit)$residuals
    else if (iniest == "lts")
      init.res = lqs(response ~ 1, method = "lts", nsamp = ninit)$residuals
    else if (iniest == "lms")
      init.res = lqs(response ~ 1, method = "lms", nsamp = ninit)$residuals
  }
  else {
    if (iniest == "ltsMy")
      init.res = lws(regressors = regressors, response = response, psi_type = 5, J = ninit)$residuals
    else if (iniest == "ls")
      init.res = lm(response ~ ., data = regressors)$residuals
    else if (iniest == "S")
      init.res = lqs(response ~ ., data = regressors, method = "S", nsamp = ninit)$residuals
    else if (iniest == "lts")
      init.res = lqs(response ~ ., data = regressors, method = "lts", nsamp = ninit)$residuals
    else if (iniest == "lms")
      init.res = lqs(response ~ ., data = regressors, method = "lms", nsamp = ninit)$residuals
  }
  
  # standardized residuals:
  med.abs.dev = mad(init.res)
  # We standardize in a robust way.
  sort.res = sort(init.res^2) / med.abs.dev^2
  sort.res[sort.res < eps_w] = eps_w
  
  # weights to make squared residuals chi-squared - distributed
  weights = qchisq( ((1:n)-0.5)/n, 1 ) / sort.res
  
  w_out <- weights[rank(init.res^2,ties.method="first")]
  
  return(w_out)
}


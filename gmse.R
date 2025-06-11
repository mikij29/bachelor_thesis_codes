# BETA_MSE
# function for computation of mean square error of coefficient (as defined in the thesis) of the given estimated coefficients

# INPUT:
# - beta_est - a vector / a matrix of estimates of beta (if matrix, the individual estimates must be in the rows)
# - - the vector / matrix selection must be specified:
# - - - type = 1 ... for a vector
# - - - type = 2 ... for a matrix
# - beta_true - true values of parameters beta (by default we take beta_true = c(1, 2, -3) as we don't use anything else in our simulation study)

# OUTPUT:
# - MSEoC of the given estimated coefficients


# MSE, TMSE, WMSE
# INPUT:
# - resids - the model residuals
# - alpha - by default 0.9; the percentage of residuals that is NOT trimmed out (for TMSE and WMSE only)
# OUTPUT:
# - MSE / TMSE / WMSE for the given residuals


source("psi_indep.R")

beta_mse = function(beta_est, type = 1, beta_true = c(1, 2, -3)) {
  if(type == 1) 
    return(sum((beta_est - beta_true)^2))
  else if (type == 2) {
    S = NROW(beta_est)
    p = NCOL(beta_est)
    mse_sum = 0
    for (i in 1:S) {
      for (j in 1:p) {
        mse_sum = mse_sum + (beta_est[i, j] - beta_true[j])^2
        # it should be possible to make this faster via sum(*) or %*%
      }
    }
    return(mse_sum)
  }
  else {
    print("Wait a second!")
    return(1)  # OK, this one row is needed for the function to work...
  }
}


mse = function(resids) {
  return(sum(resids^2)/length(resids))
}

tmse = function(resids, alpha = 0.9) {
  n = length(resids)
  trimming = floor(n*alpha)
  return(sum(sort(resids^2)[1:trimming])/trimming)
}

wmse = function(resids, alpha = 0.9) {
  n = length(resids)
  trimming = floor(n*alpha)
  return(sum( (sort(resids^2)*psi_indep(type = 4, n = n))[1:floor(n*alpha)] )/trimming)
}

# resids = c(1, 4, 3, 2)
# wmse(resids)
# tmse(resids)


# LWS
# implementation of the least weighted squares algorithm

# INPUT: 
# - response, regressors (without the intercept); 
# - psi_type - options "dd_wls" for the data-dependent weights without iterations (dd-WLS), "dd_lws" for the data-dependent weights with iterations (dd-LWS), or types 1-5 for the independent weights LWS-1 to LWS-5
# - iniest - the initial estimator - one of the set {"ltsMy", "ls", "S", "lts", "lms"}
# - ninit - the number of iterations for the initial estimator computation - by default 555
# - J - the number of iterations in the main function cycle - by default 10000
# - eps - a constant; we stop the iterations of the estimator-residuals-weights pattern (for a given j in {1, ..., J}) if the weighted residual sum of squares did not get better from the step k to the step k + 1 at least by this particular eps - default 10^(-6)
# - k_max - the maximum number of iterations of the estimator-residuals-weights pattern - default 30
# - print_k ... should we always print actual iteration "k" in the inner cycle? (i.e. values from 1 to k_max)
# - tau - the trimming constant for the function psi_indep - default 0.75
# - eps_w - the individual weights lower than this constant are replaced by this constant (in order not to divide by zero); by default 1e-33

# OUTPUT:
# - coefficients - vector of estimated coefficients
# - weights - obtained vector of weights
# - residuals - vector of model residuals


source("psi_dep.R")
source("psi_indep.R")

is_full_rank <- function(mx, tol = .Machine$double.eps^0.5) {
  # helping function for assessing correctness of the input
  if (!(is.matrix(mx) || is.data.frame(mx))) {
    stop("Input must be a matrix or data frame.")
  }
  
  n <- NROW(mx)
  p <- NCOL(mx)
  
  if (n < p) {
    stop("Number of data rows must be greater than or equal to number of columns (n >= p).")
  }
  
  # Use SVD to compute the rank
  s <- svd(mx)$d  # singular values
  rank <- sum(s > tol * max(s))  # count significant singular values
  
  return(rank == p)
}


lws = function(response, regressors = NA, psi_type = "dd_wls", iniest = "S", ninit = 555, J = 10000, eps = 10^(-6), k_max = 30, print_k = FALSE, tau = 3/4, eps_w = 1e-33) {

  icpt_only = anyNA(regressors)
  
  if(!icpt_only) {
    if(!is_full_rank(regressors))
       stop("Model matrix must be of a full rank.")
    if(!is.data.frame(regressors))
      regressors = data.frame(regressors)
  }
  
  if(psi_type == "dd_lws" || psi_type == "dd_wls") {
    dep_weights = psi_dep(regressors = regressors, response = response, iniest = iniest, ninit = ninit, eps_w = eps_w)
    if(psi_type == "dd_wls") {
      w_lws = dep_weights
      if(icpt_only)
        model = lm(response~1, weights = w_lws)
      else  
        model = lm(response~., data = regressors, weights = w_lws)
      b_lws = model$coefficients
      u_lws = model$residuals
      m_lws = w_lws%*%u_lws^2
      return(list("coefficients" = b_lws, "weights" = w_lws, "mse" = m_lws, "residuals" = u_lws))
    }
    else {
      output = iterate(regressors = regressors, response = response, icpt_only = icpt_only, weights = sort(dep_weights, decreasing = TRUE), J = J, eps = eps, k_max = k_max, print_k = print_k)
    }
  }
  else {
    if(icpt_only)
      p = 1
    else
      p = NCOL(regressors)+1
    n = NROW(response)
    indep_weights = psi_indep(type = psi_type, n = n, p = p, tau = tau)
    # types LWS-1 - LWS-5
    output = iterate(regressors = regressors, response = response, icpt_only = icpt_only, weights = indep_weights, J = J, eps = eps, k_max = k_max, print_k = print_k)
    # indep_weights are sorted (in non increasing order) by default
  }
  return(output)
}
  
iterate = function(regressors, response, icpt_only, weights, J = 10000, eps = 10^(-6), k_max = 30, print_k = FALSE) {
  # weights must be sorted
  if(icpt_only)
    p = 1
  else
    p = NCOL(regressors)+1
    # This works better than ncol. In this way it assigns the value even when 'regressors' is only a vector.
  n = NROW(response)
  
  # LWS algorithm, by articles from J. Kalina:
  m_lws = Inf
  for (j in 1:J){
    m0 = Inf
    srpp1o = sample(1:n, p+1, replace = FALSE)
    # select randomly p plus 1 observations
    if(icpt_only)
      model = lm(response[srpp1o]~1)
    else
      model = lm(response[srpp1o]~., data = data.frame(regressors[srpp1o,]))
      # necessary to use data.frame here?
    
    b0 = model$coefficients
    if(icpt_only)
      u0 = response - b0
    else
      u0 = response - cbind(1,as.matrix(regressors))%*%b0
      # u0 = model$residuals would be wrong here
    
    k = 0
    while (k <= k_max) {
      w1 = weights[rank(u0^2, ties.method = "random")]
      # every variation of ties.method is possible, except from the default one - "average" (we don't want to allow non-integer here)
      # to achieve a bijection indep_weights -> w1 every time, it's also impossible to use "max" and "min" (to author's hypothesis)
      if(icpt_only)
        model = lm(response~1, weights = w1)
      else
        model = lm(response~., data = regressors, weights = w1)
      b1 = model$coefficients
      u1 = model$residuals
      m1 = w1%*%u1^2
      # if (m1 > m0 + eps){
        # very tiny increment in MSE with weights was allowed in some previous versions
          # it might have slowed down the algorithm in some instances
      if (m1 + eps > m0){
        break
      }
      w0 = w1
      m0 = m1
      b0 = b1
      u0 = u1
      if(print_k)
        print(k)
      k = k+1
    }
    if (m0 < m_lws) {
      b_lws = b0
      w_lws = w0
      m_lws = m0
      u_lws = u0
    }
  }
  return(list("coefficients" = b_lws, "weights" = w_lws, "residuals" = u_lws))
  # "wrss" = m_lws (weighted residual sum of squares) was removed from the output, as it can be always easily calculated from weights and residuals
}


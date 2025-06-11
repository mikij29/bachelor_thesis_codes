# VAR_SANDWICH
# function for calculating variance matrices by the sandwich method, as well as some estimates

# INPUT: 
# - regressors (without the intercept), response; 
# - which estimates should be used (T/F for e_lws, e_ls)
# - if LWS, which parameters should be used (lws_1, lws_2, lws_3, lws_4, lws_5, dd_lws, dd_wls, eps_w, J, iniest, ninit)
# - var_wls - T/F - whether to estimate WLS variance matrix as derived in section 1.8.2 of the thesis
# - var_ls_homo - T/F - should we calculate the basic homoscedastic variance matrix for LS as well?

# OUTPUT
# - coefficients - the matrix of estimated coefficients (in each column there is one vector of estimated coefficients; the number of columns is equal to the number of estimators used)
# - residuals - the matrix of model residuals (in each column there is one vector of residuals)
# - varmx - an array of the estimated variance matrices for all the estimators 
# - time - the matrix of times it took to calculate the given estimator (in the rows; the columns correspond to the output of the function system.time)
# - varmx_W - weighted version of sandwich variance estimator (in thesis called "guess" variance estimator)
# - varmx_wls - WLS variance matrix estimate as derived in 1.8.2


source("lws.R")
source("load_names.R")

var_sandwich = function(response, regressors = NA, e_lws = TRUE, lws_1 = TRUE, lws_2 = TRUE, lws_3 = TRUE, lws_4 = TRUE, lws_5 = TRUE, dd_lws = TRUE, dd_wls = TRUE, var_wls = TRUE, eps_w = 1e-33, J = 10000, iniest = "S", ninit = 555, e_ls = TRUE, var_ls_homo = TRUE) {
  
  {
    ## Calculating some dimensions
    loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = F, e_MM = F, e_MMa = F, e_s = F, e_lts = F, e_lms = F, e_rdl1 = F)
    estnames = loaded_estnames$names
    lws_names = loaded_estnames$lws_names
    noest = loaded_estnames$number
    n_lws = loaded_estnames$n_lws
    
    icpt_only = anyNA(regressors)
    
    if(icpt_only)
      p = 1
    else
      p = NCOL(regressors) + 1
    n = length(response)
  }
  
  {
    ## Estimators and corresponding matrices
    beta = array(dim = c(p, noest))
    resid = array(dim = c(n, noest))
    if(e_lws)
      weights = array(dim = c(n, n_lws))
    esttimes = array(0, dim = c(noest, 5))
    dimnames(esttimes)[[1]] = estnames
    dimnames(beta)[[2]] = estnames
    dimnames(resid)[[2]] = estnames
    
    S_X = array(0, dim = c(p, p, noest))
    V_X = array(0, dim = c(p, p))
    # be aware, that the V_X matrix definition in the text differs: in the text, it's divided by 'n'
      # the same goes for S_X
    varmx = array(dim = c(p, p, noest))
    if(e_lws) {
      S_X_W = array(0, dim = c(p, p, n_lws))
      varmx_W = array(dim = c(p, p, n_lws))
    }
    else {
      varmx_W = NA
    }
    
    sx_time = array(0, dim = c(noest, 5))
    matrix_time = rep(0, 5)
    
    k = 1
    if(e_lws){
      for(lws_name in lws_names) {
        esttimes[k,] = esttimes[k,] + system.time({ 
          lws_out = lws(regressors = regressors, response = response, psi_type = lws_name, eps_w = eps_w, J = J, iniest = iniest, ninit = ninit) })
        beta[,k] = lws_out$coefficients
        resid[,k] = lws_out$residuals
        weights[,k] = lws_out$weights
        if (var_wls && lws_name == "dd_wls") {
          # no timing, this should be negligible (as well)
          sigma_sq_hat_W = sum((resid[,k])^2*weights[,k])/(n-p)
          XTWX = array(0, dim = c(p, p))
        }
        k = k + 1
      }
    }
    if(e_ls){
      esttimes[k,] = esttimes[k,] + system.time({ 
        if(icpt_only)
          ls_out = lm(response~1)
        else
          ls_out = lm(response~., data = data.frame(regressors))  })
          # I believe it's not such a biggie for the software to eventually convert data.frame to data.frame
      beta[,k] = ls_out$coefficients
      resid[,k] = ls_out$residuals
      sigma_sq_hat = sum((resid[,k])^2)/(n-p)
    }
  }
  {
    ## Matrix multiplication etc.
    for (i in 1:n){
      matrix_time = matrix_time + system.time({
        if(icpt_only)
          X_i = 1
        else
          X_i = c(1, unlist(regressors[i,]))
          # data.frame je totiz typu list
        outer_prod = X_i%*%t(X_i)
        V_X = V_X + outer_prod
      })
      k = 1
      for(estname in estnames) {
        sx_time[k,] = sx_time[k,] + system.time({
          u_ik = as.numeric(response[i] - X_i%*%beta[,k])
          S_X[,,k] = S_X[,,k] + outer_prod*(u_ik^2)
          if(k <= n_lws) 
            S_X_W[,,k] = S_X_W[,,k] + outer_prod*(u_ik^2)*weights[i,k]
        })
        if (estname == "d.no.iter") {
          XTWX = XTWX + outer_prod * weights[i,k]
        }
        k = k + 1
      }
      
    }
    matrix_time = matrix_time + system.time({ V_X_inv = solve(V_X) })
  }
  {
    ## Creating output
    for(k in 1:noest) {
      esttimes[k, ] = esttimes[k, ] + matrix_time + sx_time[k,] + system.time({
        varmx[,,k] = V_X_inv%*%S_X[,,k]%*%V_X_inv
        if(k <= n_lws)
          varmx_W[,,k] = V_X_inv%*%S_X_W[,,k]%*%V_X_inv 
      })
    }
    if(e_ls && var_ls_homo) {
      varmx_ls_homo = sigma_sq_hat * V_X_inv 
    }
    else {
      varmx_ls_homo = NA
    }
      
    if(e_lws && dd_wls && var_wls) {
      varmx_wls = sigma_sq_hat_W * solve(XTWX)
      # we didn't implement time measuring here yet (and we perhaps won't)
    }
    else {
      varmx_wls = NA
    }
    
  } 
  return(list("coefficients" = beta, "residuals" = resid, "varmx" = varmx, "time" = esttimes, "varmx_ls_homo" = varmx_ls_homo, "varmx_wls" = varmx_wls, "varmx_W" = varmx_W))
}


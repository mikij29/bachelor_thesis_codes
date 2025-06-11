# FILL_MATRICES
# an auxiliary function to be used inside the main R script "simulacni_studie.R"
# its purpose is to fill matrices with output of functions "var_sandwich" and "var_bootstrap"

fill_matrices = function(meth) {
  if(meth == "sand") {
    # k <<- 1
    while(k <= nosand) {
      # vzdy zaciname sandwichem, jinak zde k0:(k0+nosand-1)
      beta[i,,k,m] <<- swch[[1]][, k]
      xmse[i,,k,m] <<- swch[[2]][, k]
      diva[i,,k,m] <<- diag(as.matrix(swch[[3]][,,k]))
      # as.matrix here serves for dealing with intercept only (p = 1)
      if(e_lws && k <= n_lws)
        diva_W[i,,k,m] <<- diag(as.matrix(swch$varmx_W[,,k]))
      time[t,k,,m ] <<- time[t,k,,m] + swch[[4]][k, ]
      k <<- k + 1
    }
    t <<- t + 1
  }
  if(meth == "boot") {
    kk = 1
    while(kk <= noest) {
      if(boot_est > 0) {
        beta[i,,k:(k+boot_est-1),m] <<- btst[[1]][,kk,]
        xmse[i,,k:(k+boot_est-1),m] <<- btst[[2]][,kk,]
      }
      if(p == 1)
        diva[i,,k0 + (kk - 1),m] <<- btst[[3]][kk]
      else
        diva[i,,k0 + (kk - 1),m] <<- diag(array(btst[[3]][,kk], dim = c(p, p)))
      time[t,kk,,m] <<- time[t,kk,,m] + btst[[4]][kk, ]
      kk = kk + 1
      k <<- k + boot_est
    }
    t <<- t + 1
  }
  return(0)
}


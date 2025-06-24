# various functions for keeping dimensions names in order

# LOAD_ESTNAMES
# INPUT
# - which estimators should be used
# - if LWS, which parameters should be used
# OUTPUT
# - names - vector of names of all the estimates (firstly LWS with all its variants, then other estimators)
# - number - length of 'names'
# - n_lws - how many variants of LWS are to be calculated
# - lws_names - names of all the LWS variants used

# LOAD_INIESTNAMES
# INPUT
# - which initial estimators should be used
# OUTPUT
# - names - vector of names of the initial estimators used
# - number - length of 'names'

# LOAD_BOOTNAMES
# INPUT
# - whether 'bootstrap' should be used and with which of its variant of "bootstrap parameters estimator"
# OUTPUT
# - names - vector of names of the variants used
# - number - length of 'names' (0, 1 or 2)

# LOAD_BMV ("beta, method, variance")
# INPUT
# - which estimators should be used
# - if LWS, which parameters should be used
# - which methods should be used (including the variants of 'bootstrap')
# - - sand ... sandwich estimator (only for LS and LWS)
# - - boot ... bootstrap
# - - jest ... just estimate, no estimation of variance matrix
# OUTPUT
# - nobeta - number of estimates of coefficients
# - nometh - number of methods used (variation estimate by sandwich or bootstrap and/or just estimating the coefficients)
# - novar - number of variance estimators (each coefficient estimator is by means of methods 'sandwich' or 'bootstrap' accompanied by a variance estimation)
# - nosand - number of sandwich variance estimators
# - betanames - vector of names of all the estimates (in the form "method.estimator")
# - methnames - vector of names of the methods used
# - varnames - vector of variance estimators names


load_estnames = function(e_lws = T, lws_1 = TRUE, lws_2 = TRUE, lws_3 = TRUE, lws_4 = TRUE, lws_5 = TRUE, dd_lws = TRUE, dd_wls = TRUE, e_ls = T, e_mm = T, e_MM = T, e_MMa = T, e_s = T, e_lts = T, e_lms = T, e_rdl1 = F) {
  n_lws = e_lws * (lws_1 + lws_2 + lws_3 + lws_4 + lws_5 + dd_lws + dd_wls)
  noest = e_lws * n_lws + e_ls + e_mm + e_MM + e_MMa + e_s + e_lts + e_lms + e_rdl1
  if(noest == 0) {
    print("We always need at least one estimator!")
    return(1)
  }
  estnames = rep(NA, noest)
  if(e_lws) {
    if(n_lws == 0) {
      print("Turn at least one weight function on, if you want to estimate via LWS!")
      return(1)
    }
    else
      lws_names = rep(NA, n_lws)
  }
  else
    lws_names = NA
  i = 1
  if(e_lws) {
    if(lws_1) {
      estnames[i] = "lws-1"
      lws_names[i] = 1
      # we always start with lws (otherwise introduce another index)
      i = i + 1
    }
    if(lws_2) {
      estnames[i] = "lws-2"
      lws_names[i] = 2
      i = i + 1
    }
    if(lws_3) {
      estnames[i] = "lws-3"
      lws_names[i] = 3
      i = i + 1
    }
    if(lws_4) {
      estnames[i] = "lws-4"
      lws_names[i] = 4
      i = i + 1
    }
    if(lws_5) {
      estnames[i] = "lws-5"
      lws_names[i] = 5
      i = i + 1
    }
    if(dd_lws) {
      estnames[i] = "dd_lws"
      lws_names[i] = "dd_lws"
      i = i + 1
    }
    if(dd_wls) {
      estnames[i] = "dd_wls"
      lws_names[i] = "dd_wls"
      i = i + 1
    }
  }
  if(e_ls) {
    estnames[i] = "ls"
    i = i + 1
  }
  if(e_mm) {
    estnames[i] = "mm"
    i = i + 1
  }
  if(e_MM) {
    estnames[i] = "MM"
    i = i + 1
  }
  if(e_MMa) {
    estnames[i] = "MMa"
    i = i + 1
  }
  if(e_s) {
    estnames[i] = "s"
    i = i + 1
  }
  if(e_lts) {
    estnames[i] = "lts"
    i = i + 1
  }
  if(e_lms) {
    estnames[i] = "lms"
    i = i + 1
  }
  if(e_rdl1) {
    estnames[i] = "rdl1"
    i = i + 1
  }
  return(list("names" = estnames, "number" = noest, "n_lws" = n_lws, "lws_names" = lws_names))
  # n_lws not used yet
}


load_iniestnames = function(i_ls = F, i_ltsMy = F, i_lts = F, i_lms = F, i_s = TRUE) {
  noiniest = i_ls + i_ltsMy + i_lts + i_lms + i_s
  if(noiniest == 0) {
    # print("We have no initial estimators! Hope you switched off e_lws as well.")
    iniestnames = NA
  }
  else {
    iniestnames = rep(NA, noiniest)
    i = 1
    if(i_ls) {
      iniestnames[i] = "ls"
      i = i + 1
    }
    if(i_ltsMy) {
      iniestnames[i] = "ltsMy"
      i = i + 1
    }
    if(i_lts) {
      iniestnames[i] = "lts"
      i = i + 1
    }
    if(i_lms) {
      iniestnames[i] = "lms"
      i = i + 1
    }
    if(i_s) {
      iniestnames[i] = "S"
      i = i + 1
    }
  }
  return(list("names" = iniestnames, "number" = noiniest))
}


load_bootnames = function(boot = 0, boot_med = TRUE, boot_mean = TRUE) {
  if(boot <= 0 || boot_med + boot_mean == 0)
    return(list("number" = 0, "names" = NA))
  else {
    boot_est = boot_med + boot_mean
    bootnames = rep(NA, boot_est)
    i = 1
    if(boot_med) {
      bootnames[i] = "med"
      i = i + 1
    }
    if(boot_mean) {
      bootnames[i] = "mean"
      i = i + 1
    }
    return(list("number" = boot_est, "names" = bootnames))
  }
}


load_bmv = function(e_lws = T, lws_1 = TRUE, lws_2 = TRUE, lws_3 = TRUE, lws_4 = TRUE, lws_5 = TRUE, dd_lws = TRUE, dd_wls = TRUE, e_ls = T, e_mm = T, e_MM = T, e_MMa = T, e_s = T, e_lts = T, e_lms = T, e_rdl1 = F, sand = TRUE, boot = 0, boot_med = TRUE, boot_mean = TRUE, jest = TRUE) {
  
loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
  estnames = loaded_estnames$names
  n_lws = loaded_estnames$n_lws
  lws_names = loaded_estnames$lws_names
  
  loaded_boot = load_bootnames(boot = boot, boot_med = boot_med, boot_mean = boot_mean)
  boot_est = loaded_boot$number
  if(boot_est > 0) {
    bootnames = loaded_boot$names
  }
  
  nosand = e_lws * n_lws + e_ls
  if(sand && nosand == 0) {
    print("Sandwich won't be performed, we have no suitable estimators!")
    sand = FALSE
  }
  
  nometh = sand + (boot > 0) + jest
  if(nometh == 0) {
    print("What do you want to execute when I don't have any methods to execute it with?!")
    return(1)
  }
  methnames = rep(NA, nometh)
  
  nobeta = nosand * (sand + boot_est) + (e_mm + e_MM + e_MMa + e_s + e_lts + e_lms + e_rdl1) * boot_est + noest * jest
  if(nobeta == 0) {
    print("We need some estimators!")
    return(1)
  }
  betanames = rep(NA, nobeta)
  
  novar = nosand * (nometh - jest) + (e_mm + e_MM + e_MMa + e_s + e_lts + e_lms + e_rdl1) * (boot > 0)
  if(novar > 0)
    varnames = rep(NA, novar)
  else
    varnames = NA
  
  i = 1
  t = 1
  
  if(sand) {
    methnames[t] = "sandwich"
    t = t + 1
    if(e_lws) {
      for(lws_name in lws_names) {
        betanames[i] = paste0("sand.", lws_name)
        varnames[i] = paste0("sand.", lws_name)
        i = i + 1
      }
    }
    if(e_ls) {
      betanames[i] = "sand.ls"
      varnames[i] = "sand.ls"
      i = i + 1
    }
  }
  
  if(boot > 0) {
    methnames[t] = "bootstrap"
    t = t + 1
    j = i
    for(estname in estnames) {
      varnames[j] = paste0("boot.", estname)
      j = j + 1
      if(boot_est > 0) {
        for(k in 1:boot_est)
          betanames[i + k - 1] = paste0("boot.", estname, ".", bootnames[k])
        i = i + boot_est
      }
    }
  }
  
  if(jest) {
    methnames[t] = "jest"
    t = t + 1
    for(estname in estnames) {
      betanames[i] = paste0("jest.", estname)
      i = i + 1
    }
  }
  
  return(list("nobeta" = nobeta, "nometh" = nometh, "novar" = novar, "nosand" = nosand, "betanames" = betanames, "methnames" = methnames, "varnames" = varnames))
}


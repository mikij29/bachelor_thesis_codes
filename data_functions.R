# functions for processing simulation output from 'simulacni_studie.R' and creating tables based on these data

# LOAD_DATA
# loads available data saved during runs of computational part of the script simulacni_studie.R

# INPUT:
# - lenora - LENgth Of RAnks - how many elements contained the 'rank' list
# - simnum - SIMulation NUMber

# OUTPUT:
# - betab, xmseb, divab, timeb, diva_W_b, swnsw - matrices with all of the computed data available; they are further to be supplied in function 'process_data'
# - nos - Number Of Samples
# - lenora, simnum - for convenience purposes


# METHODOLOGICAL NOTE
# When simulating via cluster with multiple nodes, we run the same script with a different value of one parameter. That may serve for two disjoint reasons:
# - 1) Identification only - when storing computed values, for each run of the script we create a file (R Workspace type) with this unique value of the parameter in its name. Except of that, no object in these resulting R Workspaces depends on the value of that parameter.
# - 2) Else - each value of the parameter provides objects which we want to study in dependence on this parameter. In this case, after merging all the R Workspaces in one 'big' matrix, we need to divide it back to its original segments when performing further analysis (typically by means of for cycle).
# - - For producing tables as in thesis.pdf, only 1) is relevant.


# PROCESS_DATA
# via this function, the array of coefficient estimates is converted into an array of error metrics of interest

# INPUT:
# - betab, xmseb - from the function 'load_data'
# - MSEoC, MSE, TMSE, WMSE - T/F, whether this particular error metric should be evaluated
# - alphas - for which values of alpha should TMSE and/or WMSE be calculated
# - nos - Number Of Samples
# - lenora - LENght Of RAnks 
# - simnum - SIMulation NUMber
# - reason - value 1 or 2 - see the methodological note above

# OUTPUT:
# - valies - array of error metrics of interest


# PRODUCE_TABLE_VAR

# INPUT: 
# - data - as output from the function 'load_data'
# - simnum - SIMulation NUMber
# - permute - T/F - do we want to permute the estimators on the output?
# - - if yes, further parameters have to be provided - resulting permutations of all estimates and of sandwich variance estimates
# - kick_out - while simulating, we evaluated quite too many MM type estimators compared to others. Hence, one of them was decided not to be included in tables in thesis.pdf, i.e. it was "kicked out".
# - latex - T/F - do you want to produce a latex table? If not, the table is printed only in console

# OUTPUT:
# - latex file OR console output, which contain estimate of variance of each element of the vector of estimated coefficients


# PRODUCE_TABLE_EST

# INPUT: 
# - data - as output from the function 'load_data'
# - simnum - SIMulation NUMber
# - permute - T/F - do we want to permute the estimators on the output?
# - - if yes, further parameter has to be provided - the resulting permutation of all estimates
# - MSEoC, MSE, TMSE, WMSE - T/F - whether this particular error metric should be evaluated
# - alphas - for which values of alpha should TMSE and/or WMSE be calculated
# - mean_btst - should we include the mean bootstrap estimator results in the table as well? Note: we have no implementation for the median bootstrap estimator
# - latex - T/F - do you want to produce a latex table? If not, the table is printed only in console
# - min_bold - T/F - should the minimal value in each column be printed in bold?
# - kick_out - while simulating, we evaluated quite too many MM type estimators compared to others. Hence, one of them was decided not to be included in tables in thesis.pdf, i.e. it was "kicked out".
  
# OUTPUT:
# - latex file OR console output, which contain values of selected metrics of errors for given coefficient estimators


source("gmse.R")
library(xtable)    # for creating latex tables

load_data = function(lenora = 200, simnum = 71) {
  setwd(paste0("./qsim",simnum))
  load(paste0("qsim", simnum,  "_", 1, ".RData"))
  
  # gaining dimensions needed:
  nos = dim(beta)[1]
  p = dim(beta)[2]
  ss = dim(xmse)[2]
  noiniest = dim(beta)[4]
  nobeta = dim(beta)[3]
  nometh = dim(time)[1]
  noest = dim(time)[2]
  novar = dim(diva)[3]
  
  # b stands for 'big' - as we are storing all the values (possibly returned from more computing cores) in one matrix:
  betab = array(dim = c(lenora*nos, p, nobeta, noiniest))
  xmseb = array(dim = c(lenora*nos, ss, nobeta, noiniest))
  timeb = array(0, dim = c(nometh, noest, 5, noiniest))
  I_diva_W = !anyNA(diva_W)
  I_swch = !anyNA(swch_nonswch)
  
  # assigning names to corresponding dimensions:
  betanames = dimnames(beta)[[3]]
  estnames = dimnames(time)[[2]]
  iniestnames = dimnames(beta)[[4]]

  dimnames(betab)[[3]] = betanames
  dimnames(betab)[[4]] = iniestnames
  dimnames(xmseb)[[3]] = betanames
  dimnames(xmseb)[[4]] = iniestnames
  dimnames(timeb)[[1]] = dimnames(time)[[1]]
  dimnames(timeb)[[2]] = estnames
  dimnames(timeb)[[4]] = iniestnames
  
  if(novar > 0) {
    divab = array(dim = c(lenora*nos, p, novar, noiniest))
    dimnames(divab)[[3]] = dimnames(diva)[[3]]
    dimnames(divab)[[4]] = iniestnames
  }
  else
    divab = NA
  
  if(I_diva_W) {
    n_lws = dim(diva_W)[3]
    diva_W_b = array(dim = c(lenora*nos, p, n_lws, noiniest))
    dimnames(diva_W_b)[[3]] = dimnames(diva_W)[[3]]
    dimnames(diva_W_b)[[4]] = iniestnames
  }
  else
    diva_W_b = NA
  
  if(I_swch) {
    noswnsw = dim(swch_nonswch)[3]
    swnsw = array(dim = c(lenora*nos, p, noswnsw, noiniest))
    dimnames(swnsw)[[4]] = iniestnames
  }
  else
    swnsw = NA
  
  # main data loading cycle
  for(i in 1:lenora) {
    load(paste0("qsim", simnum, "_", i, ".RData"))
    
    if(novar > 0) {
      divab[(nos*(i-1)+1):(nos*i),,,] = diva
    }
    if(I_diva_W) {
      diva_W_b[(nos*(i-1)+1):(nos*i),,,] = diva_W
    }
    if(I_swch) {
      swnsw[(nos*(i-1)+1):(nos*i),,,] = swch_nonswch
    }
    
    betab[(nos*(i-1)+1):(nos*i),,,] = beta
    xmseb[(nos*(i-1)+1):(nos*i),,,] = xmse
    timeb = timeb + time
  }
  
  setwd("..")
  
  # in older simulated data, these values were not saved (they were not modified) - hence, here we handle the case whey they're missing:
  if(!exists("nosand"))
    nosand = 8
  if(!exists("beta_true"))
    beta_true = c(1, 2, -3)
  
  return(list("betab" = betab, "xmseb" = xmseb, "divab" = divab, "timeb" = timeb, "diva_W_b" = diva_W_b, "swnsw" = swnsw, "nos" = nos, "lenora" = lenora, "simnum" = simnum, "nosand" = nosand, "beta_true" = beta_true))
}


process_data = function(betab, xmseb, MSEoC = T, beta_true = c(1, 2, -3), MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), nos = 1, lenora = 1, simnum = 1, reason = 1) {
  # function concerning various characteristics of estimated coefficients themselves (not the variance of these estimates)
  
  noalpha = length(alphas)
  nometrics = MSEoC + MSE + noalpha * (TMSE + WMSE)
  metricsnames = array(dim = nometrics)
  
  betanames = dimnames(betab)[[3]]
  iniestnames = dimnames(betab)[[4]]
  nobeta = length(betanames)
  noiniest = length(iniestnames)
  
  if(reason == 1) {
    # when all of the data from different files are equivalent (i.e. there is no principal difference for different values of 'rank' in objects in resulting R workspaces)
    valies = array(dim = c(nometrics, lenora*nos, nobeta, noiniest))
    # if(noiniest == 1) {
      j = 1
      if(MSEoC) {
        valies[j,,,] = apply(betab, c(1, 3, 4), beta_mse, beta_true = beta_true)
        metricsnames[j] = "MSEoC"
        j = j + 1
      }
      if(MSE) {
        valies[j,,,] = apply(xmseb, c(1, 3, 4), mse)
        metricsnames[j] = "MSE"
        j = j + 1
      }
      if(TMSE) {
        for (prop in alphas) {
          valies[j,,,] = apply(xmseb, c(1, 3, 4), tmse, alpha = prop)
          metricsnames[j] = paste0("TMSE ", prop)
          j = j + 1
        }
      }
      if(WMSE) {
        for (prop in alphas) {
          valies[j,,,] = apply(xmseb, c(1, 3, 4), wmse, alpha = prop)
          metricsnames[j] = paste0("WMSE ", prop)
          j = j + 1
        }
      }
    # }
    # else {
    #   valies[1,,,] = apply(betab, c(1, 3, 4), beta_mse)
    #   valies[2,,,] = apply(xmseb, c(1, 3, 4), tmse, alpha = 0.75)
    #   valies[3,,,] = apply(xmseb, c(1, 3, 4), wmse, alpha = 0.75)
    #   valies[4,,,] = apply(xmseb, c(1, 3, 4), tmse, alpha = 0.9)
    #   valies[5,,,] = apply(xmseb, c(1, 3, 4), wmse, alpha = 0.9)
    #   valies[6,,,] = apply(xmseb, c(1, 3, 4), mse)
    # }
    dimnames(valies)[[1]] = metricsnames
    dimnames(valies)[[3]] = betanames
  }
  
  else if (reason == 2) {
    valies = array(dim = c(6, lenora, nos, nobeta, noiniest))
    # as we expect to obtain different quantities from different values of 'rank', we need to 'cut big matrices' back to corresponding blocks again; on those blocks we run 6 metrics of errors of our interest 
    for(i in 1:lenora) {
      betab_cut = betab[(nos*(i-1)+1):(nos*i), , ,]
      xmseb_cut = xmseb[(nos*(i-1)+1):(nos*i), , ,]
    
      # the following rows of code are yet to be updated
      if(noiniest == 1) {
        valies[1,i,,,] = apply(betab_cut, c(1, 3), beta_mse)
        valies[2,i,,,] = apply(xmseb_cut, c(1, 3), tmse, alpha = 0.75)
        valies[3,i,,,] = apply(xmseb_cut, c(1, 3), wmse, alpha = 0.75)
        valies[4,i,,,] = apply(xmseb_cut, c(1, 3), tmse, alpha = 0.9)
        valies[5,i,,,] = apply(xmseb_cut, c(1, 3), wmse, alpha = 0.9)
        valies[6,i,,,] = apply(xmseb_cut, c(1, 3), mse)
      }
      else {
        valies[1,i,,,] = apply(betab_cut, c(1, 3, 4), beta_mse)
        valies[2,i,,,] = apply(xmseb_cut, c(1, 3, 4), tmse, alpha = 0.75)
        valies[3,i,,,] = apply(xmseb_cut, c(1, 3, 4), wmse, alpha = 0.75)
        valies[4,i,,,] = apply(xmseb_cut, c(1, 3, 4), tmse, alpha = 0.9)
        valies[5,i,,,] = apply(xmseb_cut, c(1, 3, 4), wmse, alpha = 0.9)
        valies[6,i,,,] = apply(xmseb_cut, c(1, 3, 4), mse)
      }
    }
    dimnames(valies)[[1]] = rnms
    dimnames(valies)[[4]] = betanames
  }
  return(valies)
}

# Unluckily, this piece of data still has to be processed manualy:
  # swnswM = data$swnsw
  # swnswM
  # apply(swnswM, c(2,3), mean)


# These functions are used inside functions 'produce_table_est' and 'produce_table_var' from the script 'data_functions.R'.
permute_ests = function(ests, res_permutation_est) {
  if(length(ests) != length(res_permutation_est)) {
    print("No permuting is performed, the provided permutation does not fit the number of estimators used.")
    return(NA)
  }
  return(ests[res_permutation_est])
}


permute_vars = function(vars, res_permutation_sand) {
  if(length(vars) != length(res_permutation_sand)) {
    print("No permuting is performed, the provided permutation does not fit the number of variance estimators used.")
    return(NA)
  } 
  return(res_permutation_sand)
}

  
produce_table_var = function(data, simnum = 71, permute = T, res_permutation_est = NA, res_permutation_sand = NA, kick_out = T, latex = F) {
  # no defaults for res_permutation_est and res_permutation_sand - when permute = F, it does not matter they are not provided
  
  if(anyNA(data$divab))
    stop("We need some variance estimators to produce a table of these!")
  divabM = data$divab
  I_diva_W_b = !anyNA(data$diva_W_b)
  
  if(I_diva_W_b)
    divabMW = data$diva_W_b
    # notation: diva - diagonal variance, b - big, M - matrix, W - weighted :)
  
  iniestnames = dimnames(data$betab)[[4]]
  est_print_copy = est_print
  sand_print_copy = sand_print
  noest = dim(data$timeb)[2]
  
  if(permute) {
    if(anyNA(res_permutation_est)) {
      print("permute = TRUE and no resulting permutation of estimates is provided -> I'm switching permute to FALSE")
      permute = FALSE
    }
    else if(anyNA(res_permutation_sand)) {
      print("permute = TRUE and no resulting permutation of sandwich variance estimates is provided -> I'm switching permute to FALSE")
      permute = FALSE
    }
    else {
      perm_ests = permute_ests(est_print_copy, res_permutation_est)
      perm_vars = permute_vars(sand_print_copy, res_permutation_sand)
      if(anyNA(perm_ests) || anyNA(perm_vars))
        permute = FALSE
      else {
        est_print_copy = perm_ests
        sand_print_copy = perm_vars
        nosand = data$nosand
      }
    }
  }
  
  m = 1
  for(iniest in iniestnames) {
    if(data$lenora * data$nos == 1) {
      tab_var = divabM[,,,m]
      if(I_diva_W_b)
        tab_impr = divabMW[,,,m]
    }
    else {
      tab_var = apply(divabM[,,,m], c(2,3), mean)
      if(I_diva_W_b)
        tab_impr = apply(divabMW[,,,m], c(2,3), mean)
    }
    if(I_diva_W_b)
      tab_var[,1:7] = tab_impr
    # improved variance estimate (by using the weighted version of sandwich variance estimator)
    
    if(permute) {
      tab_var = cbind(tab_var[,permute_vars(1:nosand, res_permutation_sand)], tab_var[,nosand+permute_ests(1:noest, res_permutation_est)])
      if(dim(tab_var)[2] == length(sand_print_copy) + length(est_print_copy))
        dimnames(tab_var)[[2]] = c(sand_print_copy, paste0("B-", est_print_copy)) # "B" stands for Bootstrap
    }
    
    if(noest > 1 && kick_out) {
      save_colnames = colnames(tab_var)[!(colnames(tab_var) == "B-MMko")]
      tab_var = tab_var[,!(colnames(tab_var) == "B-MMko")]
      colnames(tab_var) = save_colnames
    }
    if(latex) {
      tab_var_final = t(formatC(tab_var, digits = 3, format = "f"))
      print(xtable(tab_var_final, type = "latex"), file = paste0("./tab_var_", simnum, "_ini_", iniest, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
    }
    else {
      tab_var_final = t(round(tab_var, digits = 3))
      print(paste0("Initial estimator - ", iniest))
      print(tab_var_final)
    }
    m = m + 1
  }
  return(0)
}

  
produce_table_est = function(data, simnum = 71, permute = T, res_permutation_est = NA, MSEoC = T, MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), mean_btst = T, latex = F, min_bold = T, kick_out = T) {
  # data... in format of an output from the function 'load_data'
  {
    valies = process_data(betab = data$betab, xmseb = data$xmseb, MSEoC = MSEoC, beta_true = data$beta_true, MSE = MSE, TMSE = TMSE, WMSE = WMSE, alphas = alphas, nos = data$nos, lenora = data$lenora, simnum = data$simnum)
    metricsnames = dimnames(valies)[[1]]
    nometrics = length(metricsnames)
    iniestnames = dimnames(data$betab)[[4]]
    noest = dim(data$timeb)[2]
    
    betanames = dimnames(valies)[[3]]
    nobeta = length(betanames)
    betanames_cut = betanames[!(grepl(".med", betanames, fixed=TRUE) | grepl("sand", betanames, fixed=TRUE))]
    if(!mean_btst)
      betanames_cut = betanames_cut[!(grepl("boot", betanames_cut, fixed=TRUE))]
    bn_cut_ind = match(betanames_cut, betanames)
    finals = length(bn_cut_ind)
    
    est_print_copy = est_print
    if(permute) {
      if(anyNA(res_permutation_est)) {
        print("permute = TRUE and no resulting permutation of estimates is provided -> I'm switching permute to FALSE")
        permute = FALSE
      }
      else {
        perm_ests = permute_ests(est_print_copy, res_permutation_est)
        if(anyNA(perm_ests))
          permute = FALSE
        else {
          est_print_copy = perm_ests
        }
      }
    }
  }
  
  m = 1
  for(iniest in iniestnames) {
    tab_prep = array(dim = c(nometrics, finals))
    dimnames(tab_prep)[[1]] = metricsnames
    dimnames(tab_prep)[[2]] = betanames_cut

    for(i in 1:nometrics) {
      if(data$lenora * data$nos == 1)
        tab_prep[i,] = valies[i,,betanames_cut,m]
      else
        tab_prep[i,] = apply(valies[i,,betanames_cut,m], 2, mean)
    }
    if(mean_btst) {
      # both regular estimate & bootstrap mean estimate in the table... not much useful, even though implemented
      if(permute) {
        tab_prep = cbind(tab_prep[,permute_ests(1:noest, res_permutation_est)], tab_prep[,permute_ests(noest+1:noest, res_permutation_est)])
        if(dim(tab_prep)[2] == 2*length(est_print_copy))
          dimnames(tab_prep)[[2]] = c(est_print_copy, paste0("B-",est_print_copy))
          # "B-" worked for identification purposes... there might have been some problem without it
      }
    }
    else
      if(permute) {
        tab_prep = tab_prep[,permute_ests(1:noest, res_permutation_est)]
        if(dim(tab_prep)[2] == length(est_print_copy))
          dimnames(tab_prep)[[2]] = est_print_copy
      }
    if(noest > 1 && kick_out) {
      save_colnames = colnames(tab_prep)[!(colnames(tab_prep) %in% c("MMko", "B-MMko"))]
      tab_prep = tab_prep[,!(colnames(tab_prep) %in% c("MMko", "B-MMko"))]
      colnames(tab_prep) = save_colnames
      # dimnames(c(2,3))[[2]] = "W"
    }
    if(latex) {
      argmins = apply(tab_prep, 1, which.min)
      tab_est_final = t(formatC(tab_prep, digits = 3, format = "f"))
      if(min_bold) {
        i = 1
        for (argmin in argmins) {
          tab_est_final[argmin, i] = paste0("\\", "multicolumn{1}{B{.}{.}{2.3}}{", tab_est_final[argmin, i], "}")
          i = i + 1
        }
      }
      print(xtable(tab_est_final, type = "latex"), file = paste0("./tab_est_", simnum, "_ini_", iniest, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL, sanitize.text.function = identity)
    }
    else {
      tab_est_final = t(round(tab_prep, digits = 3))
      print(paste0("Initial estimator - ", iniest))
      print(tab_est_final)
    }
    m = m + 1
  }
  return(0)
}


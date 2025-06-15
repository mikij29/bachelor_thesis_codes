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
# - - if yes, further parameters have to be provided (nosand - Number Of SANDwich variance estimators, noest - Number Of Estimators in total)
# - - - as well as the permutation itself - see the corresponding code chunk in the script 'simulacni_studie.R'
# - latex - T/F - do you want to produce a latex table? If not, the table is printed only in console

# OUTPUT:
# - latex file OR console output


# PRODUCE_TABLE_EST

# INPUT: 
# - data - as output from the function 'load_data'
# - simnum - SIMulation NUMber
# - MSEoC, MSE, TMSE, WMSE - T/F - whether this particular error metric should be evaluated
# - alphas - for which values of alpha should TMSE and/or WMSE be calculated
# - permute - T/F - do we want to permute the estimators on the output?
# - - if yes, further parameter has to be provided (noest - Number Of Estimators in total)
# - - - as well as the permutation itself - see the corresponding code chunk in the script 'simulacni_studie.R'
# - mean_btst - should we include the mean bootstrap estimator results in the table as well? Note: we have no implementation for the median bootstrap estimator
# - kick_out - while simulating, we evaluated quite too many MM type estimators compared to others. Hence, one of them was decided not to be included in the tables in thesis.pdf, i.e. it was "kicked out"
# - latex - T/F - do you want to produce a latex table? If not, the table is printed only in console
# - min_bold - T/F - should the minimal value in each column be printed in bold?
  
# OUTPUT:
# - latex file OR console output


source("gmse.R")

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
  if(novar > 0)
    divab = array(dim = c(lenora*nos, p, novar, noiniest))
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
    dimnames(divab)[[3]] = dimnames(diva)[[3]]
    dimnames(divab)[[4]] = iniestnames
  }
  
  if(I_diva_W) {
    n_lws = dim(diva_W)[3]
    diva_W_b = array(dim = c(lenora*nos, p, n_lws, noiniest))
    dimnames(diva_W_b)[[3]] = dimnames(diva_W)[[3]]
    dimnames(diva_W_b)[[4]] = iniestnames
  }
  
  if(I_swch) {
    noswnsw = dim(swch_nonswch)[3]
    swnsw = array(dim = c(lenora*nos, p, noswnsw, noiniest))
    dimnames(swnsw)[[4]] = iniestnames
  }
  
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
  
  return(list("betab" = betab, "xmseb" = xmseb, "divab" = divab, "timeb" = timeb, "diva_W_b" = diva_W_b, "swnsw" = swnsw, "nos" = nos, "lenora" = lenora, "simnum" = simnum))
}


process_data = function(betab, xmseb, MSEoC = T, MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), nos = 10, lenora = 200, simnum = 61, reason = 1) {
  # function concerning various characteristics of estimated coefficients themselves (not the variance of these estimates)
  
  # alphas = c(0.7, 0.75, 0.8, 0.9)
  # {MSEoC = T; MSE = T; TMSE = T; WMSE = T}
  # valies
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
        valies[j,,,] = apply(betab, c(1, 3, 4), beta_mse)
        # valies[j,,,] = apply(betab, c(1, 3, 4), beta_mse, simplify = F)
        # replacing c(1, 3) by c(1, 3, 4) suddenly works... for the case it would not again, I think "simplify = F" could have been a solution
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
  
  
produce_table_var = function(data, simnum = 71, permute = T, nosand = 8, noest = 14, latex = F) {
  
  divabM = data$divab
  divabMW = data$diva_W_b
  # notation: diva - diagonal variance, b - big, M - matrix, W - weighted :)
  
  tab_var = apply(divabM, c(2,3), mean)
  tab_impr = apply(divabMW, c(2,3), mean)
  # improved variance estimate (by using the weighted version of sandwich variance estimator)
  tab_var[,1:7] = tab_impr
  
  if(permute) {
    tab_var = cbind(tab_var[,permute_vars(1:nosand)], tab_var[,nosand+permute_ests(1:noest)])
    est_print = permute_ests(est_print)
    sand_print = permute_vars(sand_print)
    # tab_var_semifinal = cbind(tab_var[,permute_vars(1:nosand)], tab_var[,nosand+permute_ests(1:noest)])
    # est_print_p = permute_ests(est_print)
    # sand_print_p = permute_vars(sand_print)
  }
  # dimnames(tab_var_semifinal)[[2]] = c(sand_print_p, paste0("B", est_print_p))
  dimnames(tab_var)[[2]] = c(sand_print, paste0("B-", est_print)) # "B" stands for Bootstrap
  
  if(latex) {
    tab_var_final = t(formatC(tab_var, digits = 3, format = "f"))
    # print(xtable(tab_var_final, type = "latex"), file = paste0("D:/skola/bc_work/thesis-en/tab/new_var_", simnum, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
    print(xtable(tab_var_final, type = "latex"), file = paste0("./tab_var_", simnum, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
  }
  else {
    tab_var_final = t(round(tab_var, digits = 3))
    print(tab_var_final)
  }
  return(0)
}
  

produce_table_est = function(data, simnum = 71, permute = T, noest = 14, MSEoC = T, MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), mean_btst = T, latex = F, min_bold = T, kick_out = T) {
  # data... in format of an output from the function 'load_data'
  {
    valies = process_data(betab = data$betab, xmseb = data$xmseb, MSEoC = MSEoC, MSE = MSE, TMSE = TMSE, WMSE = WMSE, alphas = alphas, nos = data$nos, lenora = data$lenora, simnum = data$simnum)
    metricsnames = dimnames(valies)[[1]]
    nometrics = length(metricsnames)
    
    betanames = dimnames(valies)[[3]]
    nobeta = length(betanames)
    betanames_cut = betanames[!(grepl(".med", betanames, fixed=TRUE) | grepl("sand", betanames, fixed=TRUE))]
    if(!mean_btst)
      betanames_cut = betanames_cut[!(grepl("boot", betanames_cut, fixed=TRUE))]
    bn_cut_ind = match(betanames_cut, betanames)
    finals = length(bn_cut_ind)
    
    # est_print_p = permute_ests(est_print)
    if(permute)
      est_print = permute_ests(est_print)
  }
  tab_prep = array(dim = c(nometrics, finals))
  dimnames(tab_prep)[[1]] = metricsnames
  dimnames(tab_prep)[[2]] = betanames_cut
  for(i in 1:nometrics) {
    tab_prep[i,] = apply(valies[i,,betanames_cut,], 2, mean)
  }
  if(mean_btst) {
    # regular estimate, btst mean estimate
    tab_prep = cbind(tab_prep[,permute_ests(1:noest)], tab_prep[,permute_ests(noest+1:noest)])
    # tab_est = cbind(tab_prep[,permute_ests(1:14)], tab_prep[,permute_ests(15:28)])
    # dimnames(tab_est)[[2]] = c(est_print_p, paste0("B",est_print_p))
    dimnames(tab_prep)[[2]] = c(est_print, paste0("B-",est_print))
    # "B-" worked for identification purposes... there might have been some problem without it
  }
  else {
    # tab_est = tab_prep[,permute_ests(1:14)]
    # dimnames(tab_est)[[2]] = est_print_p
    tab_prep = tab_prep[,permute_ests(1:14)]
    dimnames(tab_prep)[[2]] = est_print
  }
  if(kick_out)
    # tab_est = tab_est[,!(colnames(tab_est) %in% c("MMko", "BMMko"))]
    tab_prep = tab_prep[,!(colnames(tab_prep) %in% c("MMko", "BMMko"))]
  if(latex) {
    argmins = apply(tab_prep, 1, which.min)
    # argmins = apply(tab_est, 1, which.min)
    # tab_est_final = t(formatC(tab_est, digits = 3, format = "f"))
    tab_est_final = t(formatC(tab_prep, digits = 3, format = "f"))
    if(min_bold) {
      i = 1
      for (argmin in argmins) {
        # if(i == 4 || i == 7)
        tab_est_final[argmin, i] = paste0("\\", "multicolumn{1}{B{.}{.}{2.3}}{", tab_est_final[argmin, i], "}")
        # else
        # tab_est_final[argmin, i] = paste0("\\", "multicolumn{1}{B{.}{.}{1.3}}{", tab_est_final[argmin, i], "}")
        i = i + 1
      }
    }
    print(xtable(tab_est_final, type = "latex"), file = paste0("./tab_est_", simnum, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL, sanitize.text.function = identity)
  }
  else {
    # tab_est_final = t(round(tab_est, digits = 3))
    tab_est_final = t(round(tab_prep, digits = 3))
    print(tab_est_final)
  }
  return(0)
}


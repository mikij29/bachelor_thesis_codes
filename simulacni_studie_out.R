# script for processing output from 'simulacni_studie.R' and creating tables based on these data

setwd("D:/skola/bc_work")
source("gmse.R")
source("load_names.R")
library(xtable)

# to do maybe: include some special cases - when diva not computed etc.

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
  divab = array(dim = c(lenora*nos, p, novar, noiniest))
  timeb = array(0, dim = c(nometh, noest, 5, noiniest))
  
  n_lws = dim(diva_W)[3]
  if(!exists("swch_nonswch"))
    swch_nonswch = swch_special
  noswnsw = dim(swch_nonswch)[3]
  
  diva_W_b = array(dim = c(lenora*nos, p, n_lws, noiniest))
  dimnames(diva_W_b)[[3]] = dimnames(diva_W)[[3]]
  
  swnsw = array(dim = c(lenora*nos, p, noswnsw, noiniest))
  
  
  # assigning names for corresponding dimensions: (check for effectivity)
  betanames = dimnames(beta)[[3]]
  estnames = dimnames(time)[[2]]
  iniestnames = dimnames(beta)[[4]]

  dimnames(timeb)[[1]] = dimnames(time)[[1]]
  dimnames(divab)[[3]] = dimnames(diva)[[3]]
  
  dimnames(betab)[[3]] = betanames
  dimnames(xmseb)[[3]] = betanames
  dimnames(timeb)[[2]] = estnames
  dimnames(betab)[[4]] = iniestnames
  dimnames(xmseb)[[4]] = iniestnames
  dimnames(divab)[[4]] = iniestnames
  dimnames(timeb)[[4]] = iniestnames
  
  if(1 == 1) {
    dimnames(diva_W_b)[[4]] = iniestnames
    dimnames(swnsw)[[4]] = iniestnames
  }
  
  for(i in 1:lenora) {
    load(paste0("qsim", simnum, "_", i, ".RData"))
    if(1 == 1) {
      diva_W_b[(nos*(i-1)+1):(nos*i),,,] = diva_W
      swnsw[(nos*(i-1)+1):(nos*i),,,] = swch_nonswch
    }
    divab[(nos*(i-1)+1):(nos*i),,,] = diva
    timeb = timeb + time
    betab[(nos*(i-1)+1):(nos*i),,,] = beta
    xmseb[(nos*(i-1)+1):(nos*i),,,] = xmse
  }
  setwd("..")
  return(list("betab" = betab, "xmseb" = xmseb, "divab" = divab, "timeb" = timeb, "diva_W_b" = diva_W_b, "swnsw" = swnsw, "nos" = nos, "lenora" = lenora, "simnum" = simnum))
}

# When simulating via cluster with multiple nodes, we run the same script with a different value of one parameter. That may serve for two different reasons:
# # 1) Identification only - when storing computed values, for each run of the script we create a file (R Workspace type) with this unique value of the parameter in its name. Except of that, no object in these resulting R Workspaces depends on the value of that parameter.
# # 2) Else - each value of the parameter provides objects which we want to study in dependence on this parameter. Hence, after merging all the R Workspaces in one 'big' matrix, we need to divide it back to its original segments when performing further analysis (typically by means of a for cycle).

process_data = function(betab, xmseb, MSEoC = T, MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), nos = 10, lenora = 200, simnum = 61, reason = 1) {
  
  noalpha = length(alphas)
  nometrics = MSEoC + MSE + noalpha * (TMSE + WMSE)
  metricsnames = array(dim = nometrics)
  
  # function concerning various characteristics of estimated coefficients themselves (not the variance of these estimates)
  betanames = dimnames(betab)[[3]]
  iniestnames = dimnames(betab)[[4]]
  nobeta = length(betanames)
  noiniest = length(iniestnames)
  
  if(reason == 1) {
    # when all of the data from different files are equivalent (i.e. there is no principal difference for different values of 'rank' in terms of objects in resulting R workspaces)
    valies = array(dim = c(nometrics, lenora*nos, nobeta, noiniest))
    if(noiniest == 1) {
      j = 1
      if(MSEoC) {
        valies[j,,,] = apply(betab, c(1, 3), beta_mse)
        metricsnames[j] = "MSEoC"
        j = j + 1
      }
      if(MSE) {
        valies[j,,,] = apply(xmseb, c(1, 3), mse)
        metricsnames[j] = "MSE"
        j = j + 1
      }
      if(TMSE) {
        for (prop in alphas) {
          valies[j,,,] = apply(xmseb, c(1, 3), tmse, alpha = prop)
          metricsnames[j] = paste0("TMSE ", prop)
          j = j + 1
        }
      }
      if(WMSE) {
        for (prop in alphas) {
          valies[j,,,] = apply(xmseb, c(1, 3), wmse, alpha = prop)
          metricsnames[j] = paste0("WMSE ", prop)
          j = j + 1
        }
      }
    }
    else {
      valies[1,,,] = apply(betab, c(1, 3, 4), beta_mse)
      valies[2,,,] = apply(xmseb, c(1, 3, 4), tmse, alpha = 0.75)
      valies[3,,,] = apply(xmseb, c(1, 3, 4), wmse, alpha = 0.75)
      valies[4,,,] = apply(xmseb, c(1, 3, 4), tmse, alpha = 0.9)
      valies[5,,,] = apply(xmseb, c(1, 3, 4), wmse, alpha = 0.9)
      valies[6,,,] = apply(xmseb, c(1, 3, 4), mse)
    }
    dimnames(valies)[[1]] = metricsnames
    dimnames(valies)[[3]] = betanames
  }
  else if (reason == 2) {
    valies = array(dim = c(6, lenora, nos, nobeta, noiniest))
    # as we expect to obtain different quantities from different values of 'rank', we need to 'cut big matrices' back to corresponding blocks again; on those blocks we run 6 metrics of errors of our interest 
    for(i in 1:lenora) {
      betab_cut = betab[(nos*(i-1)+1):(nos*i), , ,]
      xmseb_cut = xmseb[(nos*(i-1)+1):(nos*i), , ,]
    
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

  
produce_table_var = function(data, simnum = 71, nosand = 8, latex = F) {
  # data... in format of an output from the function 'load_data'
  divabM = data$divab
  divabMW = data$diva_W_b
  
  tab_var = apply(divabM, c(2,3), mean)
  tab_impr = apply(divabMW, c(2,3), mean)
  # improved variation estimate
  tab_var[,1:7] = tab_impr
  
  tab_var_semifinal = cbind(tab_var[,permute_vars(1:nosand)], tab_var[,nosand+permute_ests(1:14)])
  # 14 nahrad noest ci necim takym
  est_print_p = permute_ests(est_print)
  sand_print_p = permute_vars(sand_print)
  
  dimnames(tab_var_semifinal)[[2]] = c(sand_print_p, paste0("B", est_print_p))
  
  if(latex) {
    tab_var_final = t(formatC(tab_var_semifinal, digits = 3, format = "f"))
    print(xtable(tab_var_final, type = "latex"), file = paste0("D:/skola/bc_work/thesis-en/tab/new_var_", simnum, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
  }
  else {
    tab_var_final = t(round(tab_var_semifinal, digits = 3))
    print(tab_var_final)
  }
  return(0)
}
  
produce_table_est = function(data, simnum = 71, MSEoC = T, MSE = T, TMSE = T, WMSE = T, alphas = c(0.7, 0.75, 0.8, 0.9), mean_btst = T, latex = F, kick_out = T) {
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
    
    est_print_p = permute_ests(est_print)
  }
  tab_prep = array(dim = c(nometrics, finals))
  dimnames(tab_prep)[[1]] = metricsnames
  dimnames(tab_prep)[[2]] = betanames_cut
  for(i in 1:nometrics) {
    tab_prep[i,] = apply(valies[i,,betanames_cut,], 2, mean)
  }
  if(mean_btst) {
    # regular estimate, btst mean estimate
    tab_est = cbind(tab_prep[,permute_ests(1:14)], tab_prep[,permute_ests(15:28)])
    dimnames(tab_est)[[2]] = c(est_print_p, paste0("B",est_print_p))
    # "B" worked for identification purposes... and maybe there was some problem without it
  }
  else {
    tab_est = tab_prep[,permute_ests(1:14)]
    dimnames(tab_est)[[2]] = est_print_p
  }
  if(kick_out)
    tab_est = tab_est[,!(colnames(tab_est) %in% c("MMko", "BMMko"))]
  if(latex) {
    argmins = apply(tab_est, 1, which.min)
    tab_est_final = t(formatC(tab_est, digits = 3, format = "f"))
    i = 1
    for (argmin in argmins) {
      # if(i == 4 || i == 7)
        tab_est_final[argmin, i] = paste0("\\", "multicolumn{1}{B{.}{.}{2.3}}{", tab_est_final[argmin, i], "}")
      # else
        # tab_est_final[argmin, i] = paste0("\\", "multicolumn{1}{B{.}{.}{1.3}}{", tab_est_final[argmin, i], "}")
      i = i + 1
    }
    print(xtable(tab_est_final, type = "latex"), file = paste0("D:/skola/bc_work/thesis-en/tab/new_est_", simnum, ".tex"), only.contents = TRUE, include.colnames = FALSE, hline.after = NULL, sanitize.text.function = identity)
  }
  else {
    tab_est_final = t(round(tab_est, digits = 3))
    print(tab_est_final)
  }
  return(0)
}
swnswM = data$swnsw
swnswM
apply(swnswM, c(2,3), mean)

# rnms = c("beta_mse","tmse75","wmse75","tmse","wmse","mse")
est_print = c("LWS-1", "LWS-2", "LWS-3", "LWS-4", "LWS-5", "dd-LWS", "dd-WLS", "LS", "MMko", "MM", "REWLS", "S", "LTS", "LMS")
sand_print = c("LWS-1", "LWS-2", "LWS-3", "LWS-4", "LWS-5", "dd-LWS", "dd-WLS", "LS")
# default orderings out of the simulation study; however, we want another order on the output:
permute_ests = function(ests) return(ests[c(8, 1:5, 13, 11, 6, 7, 10, 14, 12, 9)])
permute_vars = function(ests) return(ests[c(8,1:7)])

# getwd()
# setwd("..")

simnum = 93
data = load_data(lenora = 200, simnum = simnum)
# produce_table_est(data, simnum = simnum, MSE = F, alphas = c(0.7, 0.8, 0.9), mean_btst = F, latex = F, kick_out = T)
produce_table_est(data, simnum = simnum, MSE = F, alphas = c(0.7, 0.8, 0.9), mean_btst = F, latex = T, kick_out = T)
# produce_table_var(data, simnum = simnum, nosand = 8, latex = T)
produce_table_var(data, simnum = simnum, nosand = 8, latex = T)

# setwd("..")
# getwd()
# betanames[order(betanames)]
estnames
nosand

dim(data$betab)


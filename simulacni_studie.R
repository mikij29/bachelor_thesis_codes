#### SOURCE
{
  source("generate_data.R")
  source("var_sandwich.R")
  source("var_bootstrap.R")
  source("just_estimate.R")
  source("load_names.R")
  source("fill_matrices.R")
  source("data_functions.R")
}

#### INPUT
{
  data_compute = T    # TRUE - data are simulated and saved / FALSE - no data will be simulated (e.g. when we only want to process already simulated data)
  simnum = 2406    # simulation number... helps in saving the data
  
  data_process = T    # if TRUE, saved data are loaded and processed to tables
  simnum_out = 2406    # SIMulation NUMber
  
  ### GENERAL SIMULATING PARAMETERS
  {
    rank = 1    # option for no parallel computation
    # rank = Sys.getenv("RANK")    # option for parallel computation - instead of running this same script multiple times (or setting number of samples 'nos' higher than 1), via this command and via batch .sh script submitted to a computational cluster, this script might be performed simultaneously on multiple nodes - each run is then distinguished from others via this variable 'rank'
    ## similarly with other variables to modify on cluster
    set.seed(rank)    # comment/uncomment for setting seed
    ss = 200    # sample size
    nos = 1    # number of samples
  }
  
  ### MODEL SETTINGS 
  {
    moding = 25   # model setting shortcut... see a code chunk below
    
    ## PARAMETERS FOR THE FUNCTION 'generate_data' (see the script generate_data.R for descriptions and details)
    {c0 = 1; px1 = T; dx1_type = "n"; dx1_mean = 0; dx1_par = 1; c1 = 2; px2 = T; dx2_type = "unif"; dx2_par1 = 0; dx2_par2 = 2; c2 = -3} 
    {
      ## Distribution of error terms to be considered - only one parameter at a time can be larger than zero.
      er_n = 1
      er_laplace = 0  # this one, however, requires package "nimble" to work...
      er_t = 0
    }
    {
      ## Data contamination
      c_perc = 0.2  # Which proportion of responses should be contaminated. More than one choice of this percentage of contamination might be achieved via running this script multiple times with seed set to a constant value.
      {
        ## Type, distribution and parameters of contamination
        c_vertical = F
        cdy_type = "unif_space"    # options are "n" "unif" "unif_space"
        {cdy_par1 = -50; cdy_par2 = -10}    # for "n" (normal distribution) ... these encode its mean & variance; "unif" ... defining points of the interval for uniform distribution; "unif_space" ... defining points of the first interval of uniform distribution on two disjoint intervals
        {cdy_par3 = 10; cdy_par4 = 50}    # only relevant for unif_space - defining points of the second interval
        
        c_leverage = F
        {c_x1 = TRUE; cdx1_type = "n"; cdx1_par1 = 6; cdx1_par2 = 1}
        {c_x2 = FALSE; cdx2_type = "unif"; cdx2_par1 = -5; cdx2_par2 = 1}
        
        c_hetero = FALSE
        
        {c_locmod = T; epsilon = 0.2}
      }
    }
  }
  
  ### METHODS TO BE USED
  {
    ## PARAMETERS FOR THE FUNCTION 'var_sandwich' - coefficient estimates & variance estimates via sandwich methods (see the script var_sandwich.R for descriptions and details).
    {
      sand = T    # turns on/off the whole sandwich method
      J_sand = 2    # J with lower index always means the number of LWS iterations for the corresponding method.
      # J_sand = 10000
      var_wls = T    # whether to estimate WLS variance matrix as derived in section 1.8.2 of the thesis 
      var_ls_homo = T    # should we calculate the basic homoscedastic variance matrix for LS as well?
    }

    ## PARAMETERS FOR THE FUNCTION 'var_bootstrap' (see the script var_bootstrap.R for descriptions and details)
    {
      boot = 2
      # boot = 100    # number of bootstrap replicates; if boot is set to zero, it turns off the whole bootstrap method
      J_boot = 2    
      # J_boot = 1000    # for computational purposes, by default J is of one order lower
      {boot_med = F; boot_mean = T}    # special bootstrap coefficient estimators... however, they usually do not produce interesting results; hence to be found only in attachment of the thesis
    }
    
    ## PARAMETERS FOR THE FUNCTION 'just_estimate', which just estimates coefficients with no corresponding estimate of variance (see the script just_estimate.R for descriptions and details).
    {
      jest = T    # whether to use 'just_estimate' method
      J_jest = 2
      # J_jest = 10000    # J with lower index always means the number of LWS iterations for the corresponding method.
      # Note: if (sand && jest && (e_lws || e_ls)), then we get 2 estimates of the same quantity (but each possibly with different value of J)
    }
  }
  
  ### ESTIMATES TO BE USED
  {
    ## LEAST WEIGHTED SQUARES
    {
      {e_lws = T; lws_1 = T; lws_2 = T; lws_3 = T; lws_4 = T; lws_5 = T; dd_lws = T; dd_wls = T}    # least weighted squares switch & which weight functions to use
      eps_w = 10^(-9)    # One parameter of iterative LWS procedure to be modified here... others (yet) have to be modified in the source script lws.R.
      {i_ls = F; i_ltsMy = F; i_lts = F; i_lms = F; i_s = TRUE}    # Which initial estimators to use. (ltsMy stands for author's own implemented least trimmed squares estimator (via the function lws)) Applicable only when dd_lws or dd_wls is TRUE.
      ninit = 555    # number of iterations of the initial estimator
    }
      
    ## OTHER ESTIMATORS
    {
      e_ls = T    # ordinary least squares estimator
      e_mm = T    # MM estimator via library lmrob
      e_MM = T    # MM estimator via library lmRob
      e_MMa = T    # REWLS estimator (via library lmRob - a version of MM estimation)
      
      # following 3 estimators have been implemented using the function 'lqs' from library MASS:
      e_s = T    # S estimator
      e_lts = T    # least trimmed squares estimator
      e_lms = T    # least median squares estimator
      
      e_rdl1 = FALSE  # whatever this one should have been... not implemented yet
    }
  }
  
  ### OUTPUT CHOICES
  {
    lenora = 1    # LENgth Of RAnks (it needs to be provided)
    
    {out_est = T; out_var = T}    # Whether to produce coefficient estimation and/or variance estimation tables.
    est_print = c("LWS-1", "LWS-2", "LWS-3", "LWS-4", "LWS-5", "dd-LWS", "dd-WLS", "LS", "MMko", "MM", "REWLS", "S", "LTS", "LMS")    # This vector follows the default order of estimators in tables, if they are all 'switched on'. The purpose of this vector is to give estimators in the tables other (more suitable) names.
    sand_print = c("LWS-1", "LWS-2", "LWS-3", "LWS-4", "LWS-5", "dd-LWS", "dd-WLS", "LS")    # (relevant only for variance table) The same as est_print, however here only the estimators for the sandwich method are included.
    # So these are default orders of estimators, as they come out of simulations. However, we might want to permute them...
    {permute_est = F; permute_sand = F}
    # ...hence the following parameters of resulting permutations:
    res_permutation_est = c(8, 1:5, 13, 11, 6, 7, 10, 14, 12, 9)
    res_permutation_sand = c(8, 1:7)
    # If the corresponding number of estimators used does not match the length of this particular permutation vector, no permutation is performed.
    
    kick_out = T    # if it was computed - to get rid of "MMko" estimator (as we might have too much of MM type estimators for our table...)
    latex = F    # Do we want to create latex file output? If not, tables will be printed into console.
    
    ## OPTIONS REGARDING COEFFICIENT ESTIMATES TABLE ONLY ('out_est'):
    {
      {MSEoC = T; MSE = T; TMSE = T; WMSE = T}    # which metrics of errors do we want to consider
      alphas = c(0.7, 0.75, 0.8, 0.9)    # for which values of alpha to evaluate TMSE and WMSE
      mean_btst = F    # if they were computed - include mean bootstrap estimators in the resulting table?
      min_bold = T    # should the first minimal value in each column (a column corresponds to one evaluated error metric) be printed in bold? (only if latex = T)
    }
  }    
}

{
  #### SIMULATING DATA
  if(data_compute) {
    ## Loading - of model design, numbers and names
    {
      {
        ## Loading model by shortcuts
        if(moding == 25) {
          c_perc = 0.2
          cdy_type = "unif_space"
          cdy_par1 = -50; cdy_par2 = -10; cdy_par3 = 10; cdy_par4 = 50
        }
        else if(moding == 21) {
          c_perc = 0.2
          cdy_type = "unif_space"
          cdy_par1 = -10; cdy_par2 = -5; cdy_par3 = 5; cdy_par4 = 10
        }
        else if(moding == 15) {
          c_perc = 0.1
          cdy_type = "unif_space"
          cdy_par1 = -50; cdy_par2 = -10; cdy_par3 = 10; cdy_par4 = 50
        }
        else if(moding == 11) {
          c_perc = 0.1
          cdy_type = "unif_space"
          cdy_par1 = -10; cdy_par2 = -5; cdy_par3 = 5; cdy_par4 = 10
        }
        else
          print("Moding is off.")
      }
      
      ## Loading numbers and names
      {
        if(!e_lws || !(dd_lws || dd_wls)) {i_ls = F; i_ltsMy = F; i_lts = F; i_lms = F; i_s = F}    # switching off initial estimators, when they're not needed (for their names not to be printed in output)
        loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
        noest = loaded_estnames$number
        estnames = loaded_estnames$names
        n_lws = loaded_estnames$n_lws
        lws_names = loaded_estnames$lws_names
        
        loaded_iniestnames = load_iniestnames(i_ls = i_ls, i_ltsMy = i_ltsMy, i_lts = i_lts, i_lms = i_lms, i_s = i_s)
        noiniest = loaded_iniestnames$number
        if(noiniest == 0) {
          if(e_lws && (dd_lws || dd_wls))
            stop("Wait, in this setting of estimates we need some initial estimator!")
          else {
            # for the simulation cycle to work (and for output to make sense) even without dd_lws or dd_wls
            noiniest = 1
            iniestnames = "not needed in this case."
          }
        }
        else
          iniestnames = loaded_iniestnames$names
        
        loaded_bootnames = load_bootnames(boot = boot, boot_med = boot_med, boot_mean = boot_mean)
        boot_est = loaded_bootnames$number
        # we don't need bootnames (if they exist) here 
        
        loaded_bmv = load_bmv(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1, sand = sand, boot = boot, boot_med = boot_med, boot_mean = boot_mean, jest = jest)
        nobeta = loaded_bmv$nobeta
        betanames = loaded_bmv$betanames
        nometh = loaded_bmv$nometh
        methnames = loaded_bmv$methnames
        novar = loaded_bmv$novar
        if(novar > 0)
          varnames = loaded_bmv$varnames
        
        nosand = loaded_bmv$nosand
        
        c_number = floor(c_perc*ss)
      }
    }
    
    {  
      total_time = system.time( {
        
        ## Creating arrays for storing values
        p = 1 + px1 + px2
        beta = array(0, dim = c(nos, p, nobeta, noiniest))
        # matrix for storing estimated coefficients
        xmse = array(0, dim = c(nos, ss, nobeta, noiniest))
        # matrix for storing corresponding residuals
        time = array(0, dim = c(nometh, noest, 5, noiniest))
        # not perfect, when only sandwich is performed and there are some non-sandwich estimators "= T"
        # while simulating for thesis, output of this function was not considered much...
        
        if (novar > 0) {
          diva = array(0, dim = c(nos, p, novar, noiniest))
          # diva stands for "DIagonal VAriance" (because variance of an estimator is found on the diagonal of the variance matrix)
          if(sand && nosand > 0) {
            if(var_ls_homo || var_wls)
              swch_nonswch = array(dim = c(nos, p, var_ls_homo + var_wls, noiniest))
            else
              swch_nonswch = NA
            if(e_lws) 
              diva_W = array(0, dim = c(nos, p, n_lws, noiniest))
            # W stands for "Weighted version"
            else
              diva_W = NA
          }
        }
        
        ## Main cycle
        for(i in 1:nos) {
          
          ## Data generation
          data = generate_data(ss = ss, c0 = c0, px1 = px1, dx1_type = dx1_type, dx1_mean = dx1_mean, dx1_par = dx1_par, c1 = c1, px2 = px2, dx2_type = dx2_type, dx2_par1 = dx2_par1, dx2_par2 = dx2_par2, c2 = c2, er_n = er_n, er_laplace = er_laplace, er_t = er_t, c_perc = c_perc, c_vertical = c_vertical, cdy_type = cdy_type, cdy_par1 = cdy_par1, cdy_par2 = cdy_par2, cdy_par3 = cdy_par3, cdy_par4 = cdy_par4, c_leverage = c_leverage, c_x1 = c_x1, cdx1_type = cdx1_type, cdx1_par1 = cdx1_par1, cdx1_par2 = cdx1_par2, c_x2 = c_x2, cdx2_type = cdx2_type, cdx2_par1 = cdx2_par1, cdx2_par2 = cdx2_par2, c_hetero = c_hetero, c_locmod = c_locmod, epsilon = epsilon)
          
          regressors = data$regressors
          response = data$response
          if(i == 1)    # enough to be saved only once
            beta_true = data$beta_true
          
          m = 1
          for(iniest in iniestnames) {
            
            ## Filling matrices
            k = 1
            t = 1
            
            if(sand && nosand > 0) {
              swch = var_sandwich(regressors = regressors, response = response, e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, var_wls = var_wls, eps_w = eps_w, J = J_sand, iniest = iniest, ninit = ninit, e_ls = e_ls, var_ls_homo = var_ls_homo)
              fill_matrices(meth = "sand")
              if(var_ls_homo || var_wls) {
                s = 1
                if(var_ls_homo == TRUE) {
                  swch_nonswch[i,,s,m] = diag(as.matrix(swch$varmx_ls_homo))
                  s = s + 1
                }
                if(var_wls == TRUE) {
                  swch_nonswch[i,,s,m] = diag(as.matrix(swch$varmx_wls))
                  # as.matrix ... for cases with intercept only - some weird behaviour was observed otherwise
                  s = s + 1
                }
              }
            }
            
            if(boot > 0) {
              btst = var_bootstrap(regressors = regressors, response = response, R = boot, e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, eps_w = eps_w, J = J_boot, iniest = iniest, ninit = ninit, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1, boot_med = boot_med, boot_mean = boot_mean)
              k0 = k
              fill_matrices(meth = "boot")
            }
            
            if(jest) {
              just_estimated = just_estimate(regressors = regressors, response = response, e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, eps_w = eps_w, J = J_jest, iniest = iniest, ninit = ninit, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
              beta[i,,k:nobeta,m] = just_estimated$coef
              xmse[i,,k:nobeta,m] = just_estimated$resid
              time[t,,,m] = just_estimated$time
            }
            m = m + 1
          }
        }
      })
      
      ## Assigning names to the corresponding dimensions
      {
        dimnames(beta)[[1]] = paste0("sample", seq(1, nos))
        dimnames(xmse)[[1]] = paste0("sample", seq(1, nos))
        if(sand && nosand > 0) {
          if(e_lws) {
            dimnames(diva_W)[[1]] = paste0("sample", seq(1, nos))
            dimnames(diva_W)[[2]] = paste0("var.beta", seq(1, p))
            dimnames(diva_W)[[3]] = lws_names
            dimnames(diva_W)[[4]] = iniestnames
          }
          if(var_ls_homo && var_wls)  # :)
            dimnames(swch_nonswch)[[3]] = c("var_ls_homo", "varmx_wls")
        }
        
        dimnames(beta)[[2]] = paste0("beta", seq(1, p))
        dimnames(xmse)[[2]] = paste0("res.beta", seq(1, ss))
        dimnames(beta)[[3]] = betanames
        dimnames(xmse)[[3]] = betanames
        dimnames(beta)[[4]] = iniestnames
        dimnames(xmse)[[4]] = iniestnames
        
        dimnames(time)[[1]] = methnames
        dimnames(time)[[2]] = estnames
        dimnames(time)[[3]] = c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
        dimnames(time)[[4]] = iniestnames
        
        if(novar > 0) {
          dimnames(diva)[[1]] = paste0("sample", seq(1, nos))
          dimnames(diva)[[2]] = paste0("var.beta", seq(1, p))
          dimnames(diva)[[3]] = varnames
          dimnames(diva)[[4]] = iniestnames
        }
      }
      
      dir.create(paste0("./qsim", simnum))
      save(beta, xmse, diva, time, diva_W, swch_nonswch, beta_true, nosand, file = paste0("./qsim", simnum, "/qsim", simnum, "_", as.integer(rank), ".RData"))
      # my way of saving the data - compatible with 'simulacni_studie_out.R'
    }
  }
  
  #### PROCESSING DATA
  if(data_process) {
    simulated_data = load_data(lenora = lenora, simnum = simnum_out)
    
    if(out_est) {
      if(!exists("res_permutation_est"))
        res_permutation_est = NA
      produce_table_est(simulated_data, simnum = simnum_out, permute = permute_est, res_permutation_est = res_permutation_est, MSEoC = MSEoC, MSE = MSE, TMSE = TMSE, WMSE = WMSE, alphas = alphas, mean_btst = mean_btst, latex = latex, min_bold = min_bold, kick_out = kick_out)
    }
    if(out_var) {
      if(!exists("res_permutation_est"))
        res_permutation_est = NA
      if(!exists("res_permutation_sand"))
        res_permutation_sand = NA
      produce_table_var(simulated_data, simnum = simnum_out, permute = permute_sand, res_permutation_est = res_permutation_est, res_permutation_sand = res_permutation_sand, kick_out = kick_out, latex = latex)
    }
  }
}


# warnings()
# beta
# xmse
# diva
# time
# diva_W
# swch_nonswch


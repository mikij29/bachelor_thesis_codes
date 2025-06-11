# setwd("D:/skola/bc_work")
{
  #### SOURCE
  source("generate_data.R")
  source("var_sandwich.R")
  source("var_bootstrap.R")
  source("just_estimate.R")
  source("load_names.R")
  source("lws.R")
  source("fill_matrices.R")
}

{
  #### INPUT
  {
    ### GENERAL SIMULATING PARAMETERS
    
    simnum = 1   # simulation number... helps in saving the data
    rank = 1  # option for no parallel computation
    # rank = Sys.getenv("RANK")   # option for parallel computation - instead of running this same script multiple times (or setting number of samples 'nos' higher than 1), via this command and via batch .sh script submitted to computational cluster this script might be performed simultaneously on multiple nodes - each run is distinguished from others via this variable 'rank'
    set.seed(rank)  # comment/uncomment for setting seed
    
    ss = 200  # sample size
    nos = 1   # number of samples
  }
  
  {
    ### MODEL SETTINGS 
    
    moding = 25   # model shortcuts... see a code chunk below
    
    ## PARAMETERS FOR THE FUNCTION 'generate_data' (see the script generate_data.R for descriptions and details)
    {c0 = 1; px1 = T; dx1_type = "n"; dx1_mean = 0; dx1_par = 1; c1 = 2; px2 = T; dx2_type = "unif"; dx2_par1 = 0; dx2_par2 = 2; c2 = -3} 
    {
      ## Distribution of error terms to be considered - only one at a time
      er_n = 1
      er_laplace = 0  # this one, however, requires package "nimble" to work...
      er_t = 0
    }
    {
      ## Contamination & distribution of contaminants
      # used to be: a cycle over the "main cycle"
      # now: only one contam. at a time - the cycle over multiple contamination is inactive
      {
        ## Percentage / number of contaminated observations
        # c_perc_set = c(0.1, 0.2, 0.3, 0.4)
        c_perc = 0.2  # which proportion of responses should be contaminated
        
        # c_number = FALSE
        # if(!c_number)
        #   c_number = floor(c_perc*ss)
      }
      
      {
        ## Type and distribution of contamination
        c_vertical = F
        {
          cdy_type = "unif_space"
          # "n" "unif" "unif_space"
          cdy_par1 = -50
          cdy_par2 = -10
          # "n" ... mean & variance
          # "unif" ... defining points of the interval
          # "unif_space" ... defining points of the first interval
          # for unif_space:
          cdy_par3 = 10
          cdy_par4 = 50
        }
        
        c_leverage = F
        {
          c_x1 = TRUE
          cdx1_type = "n"
          cdx1_par1 = 6
          cdx1_par2 = 1
          
          c_x2 = FALSE
          cdx2_type = "unif"
          cdx2_par1 = -5
          cdx2_par2 = 1
        }
        
        c_hetero = FALSE
        
        c_locmod = T
        {
          epsilon = 0.2
        }
      }
    }
  }
  
  {
    ### METHODS TO BE USED
    
    # J with lower index... number of corresponding LWS iterations
    
    
    ## PARAMETERS FOR THE FUNCTION 'var_sandwich' (see the script var_sandwich.R for descriptions and details)
    
    # coefficient estimates + variance estimates via sandwich methods
    
    sand = T    # turns on/off whole sandwich method
    J_sand = 2
    # J_sand = 10000
    
    var_ls_homo = T
    var_wls = T
    # whether to produce special variance matrices
    
    ## PARAMETERS FOR THE FUNCTION 'var_bootstrap' (see the script var_bootstrap.R for descriptions and details)
    
    boot = 2
    # boot = 100    # number of bootstrap replicates
    J_boot = 2
    # J_boot = 1000   # for computational purposes, J is of one order lower
    
    boot_med = F
    boot_mean = T
    # special bootstrap coefficient estimators... however, they usually do not produce interesting results; only in attachment of the thesis
    
    
    ## PARAMETERS FOR THE FUNCTION 'just_estimate' (see the script just_estimate.R for descriptions and details)
    
    # just estimate coefficients - with no corresponding estimate of variance
    
    jest = T
    J_jest = 2
    # J_jest = 10000
    # if (sand && jest && (e_lws || e_ls)), then we get 2 estimates of the same quantity (but possibly with different values of J)
  }
  
  {
    ### ESTIMATES TO BE USED
    
    # all functions obey this ordering
    
    # least weighted squares
    {
      e_lws = T
      lws_1 = T
      lws_2 = T
      lws_3 = T
      lws_4 = T
      lws_5 = T
      dd_lws = T
      dd_wls = T
      
      eps_w = 10^(-9) # one parameter of LWS procedure to be modified here... others (yet) have to be modified in the source script lws.R itself
    
      # initial estimator
      sim_oie = F   # Should we simulate the initial estimator only - without the following LWS/WLS procedure?
      
      i_ls = F
      i_ltsMy = F   # my own implemented least trimmed squares estimator (via the function lws)
      i_lts = F
      i_lms = F
      i_s = TRUE
      ninit = 555
    }
    
    # Other estimators
    e_ls = T
    e_mm = T   # MM estimator via library lmrob
    e_MM = T   # MM estimator via library lmRob
    e_MMa = T   # REWLS estimator (via library lmRob - a version of MM estimation)
    e_s = T
    e_lts = T
    e_lms = T
    e_rdl1 = FALSE  # whatever this one should have been... not implemented yet
  }
}

{
  #### LOADING
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
    loaded_estnames = load_estnames(e_lws = e_lws, lws_1 = lws_1, lws_2 = lws_2, lws_3 = lws_3, lws_4 = lws_4, lws_5 = lws_5, dd_lws = dd_lws, dd_wls = dd_wls, e_ls = e_ls, e_mm = e_mm, e_MM = e_MM, e_MMa = e_MMa, e_s = e_s, e_lts = e_lts, e_lms = e_lms, e_rdl1 = e_rdl1)
    noest = loaded_estnames$number
    estnames = loaded_estnames$names
    n_lws = loaded_estnames$n_lws
    lws_names = loaded_estnames$lws_names
    
    loaded_iniestnames = load_iniestnames(i_ls = i_ls, i_ltsMy = i_ltsMy, i_lts = i_lts, i_lms = i_lms, i_s = i_s)
    noiniest = loaded_iniestnames$number
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
  }
}

# for(c_perc in c_perc_set) {
# this used to be a cycle, but the simulation studies for thesis happened to run without it

{  
  #### MAIN CODE
  
  c_number = floor(c_perc*ss)
  
  total_time = system.time( {
    
    ## Creating arrays for storing values
    p = 1 + px1 + px2
    time = array(0, dim = c(nometh, noest, 5, noiniest))
      # not perfect, when only sandwich is performed and there are some non-sandwich estimators "= T"
      # while simulating for thesis, outputs of this function were not considered much...
    beta = array(0, dim = c(nos, p, nobeta, noiniest))
    xmse = array(0, dim = c(nos, ss, nobeta, noiniest))
    
    if (novar > 0) {
      diva = array(0, dim = c(nos, p, novar, noiniest))
      if(sand && nosand > 0) {
        if(var_ls_homo || var_wls)
          swch_nonswch = array(dim = c(nos, p, var_ls_homo + var_wls, noiniest))
        else
          swch_nonswch = NA
        if(e_lws) 
          diva_W = array(0, dim = c(nos, p, n_lws, noiniest))
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

      m = 1
      for(iniest in iniestnames) {
        
        if(sim_oie) {
          # elsewhere than here, we don't need the function lws(...) in this script
          if (iniest == "ltsMy")
            imod = lws(regressors = regressors, response = response, psi_type = 5, J = ninit, eps_w = eps_w)
          else if (iniest == "lts")
            imod = lqs(response ~ ., data = regressors, method = "lts", nsamp = ninit)
          else if (iniest == "lms")
            imod = lqs(response ~ ., data = regressors, method = "lms", nsamp = ninit)
          else if (iniest == "s")
            imod = lqs(response ~ ., data = regressors, method = "S", nsamp = ninit)
          else if (iniest == "ls")
            imod = lm(response ~ ., data = regressors)
            
          beta[i,,1,m] = imod$coefficients
          xmse[i,,1,m] = imod$residuals
        }
        
        else {
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
        }
        m = m + 1
      }
    }
  })
  
  {
    ## Assigning names to the corresponding dimensions
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
  save(beta, xmse, diva, time, diva_W, swch_nonswch, file = paste0("./qsim", simnum, "/qsim", simnum, "_", as.integer(rank), ".RData"))
  # my way of saving the data - compatible with 'simulacni_studie_out.R'
}

beta
xmse
diva
time
diva_W
swch_nonswch


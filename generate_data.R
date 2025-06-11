# GENERATE_DATA
# generates data in the form Y = c0 + c1 * x1 + c2 * x2 + eps

# OUTPUT: 
# - regressors in the form c(x1, x2)
# - response Y

# INPUT parameters:
# - ss - sample size
# - c0, c1, c2 - optional constant values

# - px1 - T/F - should the variable x1 be present?
# - dx1_type - the distribution of x1 - "n" for normal, "exp" for exponential, "t" for t-distribution
# - dx1_mean - the location parameter for the normal distribution
# - dx1_par
# - - the scale parameter for the normal distribution
# - - the rate for the exponential distribution
# - - degrees of freedom for the t-distribution

# - px2 - T/F - should the variable x2 be present?
# - dx2_type - the distribution of x2 - "unif" for uniform
# - dx2_par1 - the left boundary point of the interval on which the uniform distribution will be performed
# - dx2_par2 - the right boundary point of the interval on which the uniform distribution will be performed

# - er_n, er_laplace, er_t - exactly one should be set non-zero
# - - er_n - variance for normally distributed errors
# - - er_laplace - scale for errors with the Laplace distribution
# - - er_t - degrees of freedom for errors with t-distribution

# - c_perc - the percentage of contaminated observations
# - c_vertical - T/F - should the response be contaminated vertically?
# - - cdy_type = "n" for normal "N(cdy_par1, cdy_par2)", "unif" for the uniform distribution on (cdy_par1, cdy_par2), "unif_space" for the uniform distribution on the union of intervals (cdy_par1, cdy_par2) and (cdy_par3, cdy_par4)
# - - c_leverage - T/F - should the regressors be contaminated as well?
# - - - c_x1 - T/F - should the regressor x1 be contaminated? 
# - - - - cdx1_type - "n" for normal distribution N(cdx1_par1, cdx1_par2), "unif" for uniform distribution on (cdx1_par1, cdx1_par2)
# - - - c_x2 - T/F - should the regressor x2 be contaminated? 
# - - - - cdx2_type - "n" for normal distribution N(cdx2_par1, cdx2_par2), "unif" for uniform distribution on (cdx2_par1, cdx2_par2)
# - c_hetero - T/F - should the heteroscedasticity be present in our data? If yes, the random errors have the exponential distribution with rate (x1+x2)
# - c_locmod - T/F - should the local modification of regressors be performed?
# - - epsilon - a positive number determining uniform distribution on the interval (-epsilon, epsilon), by which regressors x1 and x2 will be locally modified


generate_data = function(ss = 200, c0 = 1, px1 = TRUE, dx1_type = "n", dx1_mean = 0, dx1_par = 1, c1 = 2, px2 = TRUE, dx2_type = "unif", dx2_par1 = 0, dx2_par2 = 2, c2 = -3, er_n = 1, er_laplace = 0, er_t = 0, c_perc = 0.1, c_vertical = TRUE, cdy_type = "unif_space", cdy_par1 = -10, cdy_par2 = -5, cdy_par3 = 5, cdy_par4 = 10, c_leverage = FALSE, c_x1 = FALSE, cdx1_type = "n", cdx1_par1 = 6, cdx1_par2 = 1, c_x2 = FALSE, cdx2_type = "n", cdx2_par1 = -5, cdx2_par2 = 1, c_hetero = FALSE, c_locmod = FALSE, epsilon = 0.2) {
  
  if(px1) {
    if(dx1_type == "n")
      x1 = rnorm(ss, dx1_mean, dx1_par)
    else if(dx1_type == "exp") {
      x1 = rexp(ss, dx1_par)
      # dx1_par... rate (mean = 1/rate)
    }
    else if(dx1_type == "t") {
      x1 = rt(ss, dx1_par)
      # centered t dist.
      # dx1_par... df
    }
  }
  else
    x1 = 0
  
  if(px2) {
    if(dx2_type == "unif")
      x2 = runif(ss, dx2_par1, dx2_par2)
  }
  else
    x2 = 0
      
  if (er_n > 0)
    eps = rnorm(ss, 0, er_n)
  else if (er_laplace > 0)
    eps = rdexp(ss, scale = er_laplace)
  else if (er_t > 0)
    eps = rt(ss, df = er_t)
  else
    print("Wait, we need some distribution!")
  
  c_number = floor(c_perc*ss)
  if(c_vertical + c_leverage + c_hetero + c_locmod > 0) {
    rows = sample(1:ss, size=c_number, replace = FALSE)
    
    if((c_vertical || c_leverage) + c_hetero + c_locmod > 1)
      print("Are U sure what you're doing?")
  
    if(c_vertical) {
      if(cdy_type == "n")
        eps[rows] = rnorm(c_number, cdy_par1, cdy_par2)
      else if(cdy_type == "unif")
        eps[rows] = runif(c_number, cdy_par1, cdy_par2)
      else if(cdy_type == "unif_space") {
        for (row in rows) {
          coin = rbinom(1, 1, 1/2)
          if (coin == 0)
            eps[row] = runif(1, cdy_par1, cdy_par2)
          else if (coin == 1)
            eps[row] = runif(1, cdy_par3, cdy_par4)
        }
      }
      
      if(c_leverage) {
        if(px1 && c_x1) {
          if(cdx1_type == "n")
            x1[rows] = rnorm(c_number, cdx1_par1, cdx1_par2)
          else if(cdy_type == "unif")
            x1[rows] = runif(c_number, cdx1_par1, cdx1_par2)
        }
        if(px2 && c_x2) {
          if(cdx2_type == "unif")
            x2[rows] = runif(c_number, cdx2_par1, cdx2_par2)
          else if(cdx2_type == "n")
            x2[rows] = rnorm(c_number, cdx2_par1, cdx2_par2)
        }
      }
    }
    
    else if(c_hetero) {
      # eps[rows] = rnorm(c_number, 0, exp(x1+x2))
      eps = rnorm(ss, 0, exp(x1+x2))
    }
    
    else if(c_locmod) {
      if(px1) {
        z1 = runif(ss, -epsilon, epsilon)
        x1 = x1 + z1
      }
      if(px2) {
        z2 = runif(ss, -epsilon, epsilon)
        x2 = x2 + z2
      }
    }
  }
  
  y = c0 + c1*x1 + c2*x2 + eps
  
  if(!px2 && !px1)
    regressors = NA
  else if(!px2)
    regressors = x1
  else if(!px1)
    regressors = x2
  else
    regressors = data.frame(cbind(x1, x2))
  
  return(list("regressors" = regressors, "response" = y))
}


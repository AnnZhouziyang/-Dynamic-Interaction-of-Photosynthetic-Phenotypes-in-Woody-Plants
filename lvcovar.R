
SAD1.get_mat <- function (par0, times, traits = 1, options = list()) {
  
  par <- par0
  if (class(par0) == "list") 
    par <- unlist(par0)
  t_len <- length(times)
  SAD.1 <- array(0, dim = c(t_len * traits, t_len * traits))
  for (i0 in 1:traits) for (i1 in 1:traits) {
    if (i0 == i1) 
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        SAD.1[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                k1] <- abs(par[i0 * 2]) * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
        SAD.1[(i0 - 1) * t_len + k1, (i1 - 1) * t_len + 
                k0] <- abs(par[i0 * 2]) * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
      }
  }
  return(SAD.1)
}

SAD1.get_inv_mat <- function (par, times, traits = 1, options = list()) {
  
  SAD.1 <- SAD1.get_mat(par, times, traits, options)
  return(solve(SAD.1))
}

SAD1.get_mat_det <- function (par, times, traits = 1, options = list()) {
  
  SAD.1 <- SAD1.get_mat(par, times, traits, options)
  return(det(SAD.1))
}


SAD3.get_mat <- function (par0, times, traits = 1, options = list()) {
  
  par <- par0
  if (class(par0) == "list") 
    par <- unlist(par0)
  t_len <- length(times)
  SAD.3 <- array(0, dim = c(t_len * traits, t_len * traits))
  for (i0 in 1:traits) for (i1 in 1:traits) {
    if (i0 == i1) 
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        if(k0==k1){
          SAD.3[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                  k1] <- abs(par[i0 * 2])^2 * ((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
        }
        if(k0<k1){
          SAD.3[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                  k1] <- abs(par[i0 * 2])^2 *par[i0 * 2 - 1]^(k1 - k0)*(1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2)
          SAD.3[(i0 - 1) * t_len + k1, (i1 - 1) * t_len + 
                  k0] <- abs(par[i0 * 2])^2 *par[i0 * 2 - 1]^(k1 - k0)* (1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2)
        }
      }
    if (i0<i1)
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        if(k0==k1){
          SAD.3[(i1-1) * t_len + k0, (i0-1 ) * t_len + 
                  k1] <- par[2*traits+1]*par[i0*2]*par[i1*2]
          SAD.3[(i0-1) * t_len + k0, (i1-1) * t_len + 
                  k1] <- par[2*traits+1]*par[i0*2]*par[i1*2]
        }
        if(k0<k1){
          SAD.3[(i1-1) * t_len + k0, (i0-1 ) * t_len + 
                  k1] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))*par[2*traits+1]
          SAD.3[(i1-1) * t_len + k1, (i0-1 ) * t_len + 
                  k0] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))*par[2*traits+1]
          SAD.3[(i0-1) * t_len + k0, (i1-1 ) * t_len + 
                  k1] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))*par[2*traits+1]
          SAD.3[(i0-1) * t_len + k1, (i1-1 ) * t_len + 
                  k0] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))*par[2*traits+1]
        }
      }
  }
  return(SAD.3)
}





SAD3.get_inv_mat <- function (par, times, traits = 1, options = list()) {
  
  SAD.3 <- SAD3.get_mat(par, times, traits, options)
  return(solve(SAD.3))
}

SAD3.get_mat_det <- function (par, times, traits = 1, options = list()) {
  
  SAD.3 <- SAD3.get_mat(par, times, traits, options)
  return(det(SAD.3))
}
setrqn.mle <- function(s.par,s.y,s.l){
  A <- sum((s.y - lvetrqn.get_mu(s.par,s.l))^2 )
  A
}
sany.mle <- function(s.par,s.y,s.l,x1,x2){
  A <- sum((s.y - lvany.get_mu(s.par,s.l,x1,x2))^2 )
  A
}
lv.f <-function( parameters, state, light ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*X*(1-(X/k1))+r1*X*Y^theta21*(log(k2)-log(Y))^((b2+1)/b2)
            dY <- a2*b2*Y*(log(k2)-log(Y))^((b2+1)/b2)+r2*X^theta12*Y*(1-(X/k1))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = light, func = Lorenz, parms = parameters,method="rk4") );
  out;
}
lv.f1 <-function( parameters, state, light ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*X*(1-(X/k1))
            dY <- a2*b2*Y*(log(k2)-log(Y))^((b2+1)/b2)
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = light, func = Lorenz, parms = parameters,method="rk4") );
  out;
}
lvetrqn.get_mu <- function(par, light, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      r1 = par[3],
      theta12 = par[4],
      a2 = par[5],
      b2 = par[6],
      k2 = par[7],
      r2 = par[8],
      theta21 = par[9]);
  }
  
  state0 <- c(X=0.0217102, Y=0.011);
  y <- lv.f( par0, state0, light );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.light <- 1:length(light)
  return ( c(y[index.light, 2:3] ) );
}

lvetrqn.get_mu1 <- function(par, light, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      a2 = par[3],
      b2 = par[4],
      k2 = par[5]);
  }
  
  state0 <- c(X=0.0217102, Y=0.011);
  y <- lv.f1( par0, state0, light );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.light <- 1:length(light)
  return ( c(y[index.light, 2:3] ) );
}
lvany.get_mu<-function(par,light,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      r1 = par[3],
      theta12 = par[4],
      a2 = par[5],
      b2 = par[6],
      k2 = par[7],
      r2 = par[8],
      theta21 = par[9]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- lv.f( par0, state0, light );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.light <- 1:length(light)
  return ( c(y[index.light, 2:3] ) ); 
}
lvany.get_mu1<-function(par,light,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      k1 = par[2],
      a2 = par[3],
      b2 = par[4],
      k2 = par[5]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- lv.f1( par0, state0, light );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.light <- 1:length(light)
  return ( c(y[index.light, 2:3] ) ); 
}
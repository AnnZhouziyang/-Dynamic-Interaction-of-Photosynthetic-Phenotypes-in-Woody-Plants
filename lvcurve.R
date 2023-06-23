LV.EQ.H0 <- function(y11,y22){
  
  allpheno <- cbind(y11,y22)
  
  
  parin2 <- c( 0.6860704 ,0.3882337,0.7157615,0.2485133 ,0.2758537,
               0.88979811,0.99885548,2.05712064,-0.27184401,0.04035019,9.34159464,0.91603154,0.20793071,2.22845128)
  mpheno <- as.numeric(colMeans(allpheno))
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin2);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin2[6:14])
      AA <- curve.mle(nnpar,y=allpheno,light2 =light2,x1=mpheno[1],x2=mpheno[12])
      AA
    }
    r1.covar <- optim(parin2[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=allpheno,light2 =light2,x1=mpheno[1],x2=mpheno[12])
      AA
    }
    r1 <- optim(c(parin2[6:14]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin2,y=allpheno,light2 =light2,x1=mpheno[1],x2=mpheno[12])
  return(c(LL,parin2))
}



curve.mle <-function( par,y,light2,x1,x2)
{
  len.cov <- 5
  par.covar <- par[1:len.cov]
  n  <- length(y[,1])
  light2<-light2
  sigma <- SAD3.get_mat(par.covar,light2, 2)
  
  curve.par <- par[(len.cov+1):(len.cov+ 9)]
  mu <- lvany.get_mu(curve.par,light2,x1,x2)
  
  yy <- y
  fy <- dmvnorm(yy,mu,sigma)
  A <- -sum(log(fy))
  return(A)
}
LV.EQ.H1 <- function(y11,y22,SNP1,par1,times){
  
  SNP1<-as.numeric(SNP1)
  index <- table(SNP1)
  snp.type <- as.numeric(names(index))
  
  g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(SNP1==snp.type[i])
    SNP.p <- as.numeric(c(colMeans(cbind(y11,y22)[SNP.n,],na.rm=T)))
    r1 <- optim(par1[6:14],sany.mle,s.y=SNP.p,s.l=light2,x1=SNP.p[1],x2=SNP.p[12],
                method="BFGS",control=list(maxit=32000))
    par <- r1$par
    g.par <- c(g.par,par)
    SNP.index[[i]] <- SNP.n
  }
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  parin <- c(par1[1:5],g.par)
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[-(1:5)])
      AA <- mle.fun(nnpar,y11=y11,y22=y22,light2=light2,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- mle.fun(nnpar,y11=y11,y22=y22,light2=light2,SNP.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[-(1:5)]),mle1.g,method = "BFGS",control=list(maxit=32000))
    
    newpar <- c(new.covar1,r1.g$par)
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}

mle.fun <- function(par,y11=y11,y22=y22,light2=light2,
                    SNP.index=SNP.index,snp.type=snp.type){
  
  Y1 <- cbind(y11,y22)
  n  <- length(y11[,1])
  
  len.cov <- 5
  par.covar <- par[1:len.cov]
  
  sigma <- SAD3.get_mat(par.covar,light2, 2)
  len.gen <- 9
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    mu.g <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    yy1 <- Y1[SNP.index[[i]],]
    
    nyy1 <- yy1
    Myy1 <- as.numeric(colMeans(nyy1,na.rm=T))
    mu <- lvany.get_mu(mu.g, light2,x1=Myy1[1],x2=Myy1[12])
    
    fy1 <- dmvnorm( nyy1, mu, sigma)
    A1 <- c(A1,-sum(log(fy1)))
    len <- len + len.gen
  }
  A <- sum(A1)
  return (A);
}
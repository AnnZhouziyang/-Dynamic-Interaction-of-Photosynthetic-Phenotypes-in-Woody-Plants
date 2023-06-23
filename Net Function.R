###FUNMAP
Logistic <- function(par,t){
  
  y <- (par[1]/(1+par[2]*exp(-par[3]*t))-(par[4]*exp(-par[5]*t)))
  return(y)
  
} 
SAD1 <- function(par, times = t) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi<- par[1]   
  v2 <- par[2]
  
  tmp <- (1-phi^2)   
  sigma <- array(1, dim=c(n,n))   
  
  for(i in 1:n)   
  {     
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp     
    sigma[i:n,i] <- sigma[i,i:n]   
  }   
  sigma <- sigma * abs(v2)  
  
  return(sigma)
}
Likelihood = function(pheno,t,par){
  miu = Logistic(par[1:5],t)
  sigma = SAD1(par = c(par[6],par[7]),times = t)
  L0 = c()
  L0 = sum(dmvnorm(pheno,miu,sigma,log = T))
  return(-L0)
}
optim_deal <- function(pheno_deal,t){
  library(mvtnorm)
  itime <- 100
  itimes <- 1
  #par0 <- c(get_initial_par(pheno_ck,t),get_initial_par(pheno_salt,t),0.1,0.1,0.1,0.1)
  
  #par0 <- c( 19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,1.0660809,1.8476307,1.0659458,3.7009838 )
  par0 <- c(24.52088643,11.64758911,0.06179968,8.74285591,0.09998942,0.1,0.1)
  #par0 <- c(19.3898481,2.9533945,0.1775538,1.86229345  ,0.01507654,0.1,0.1)
  repeat{
    a <- optim(par=par0,Likelihood,pheno=pheno_deal,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood,pheno=pheno_deal,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H0 Finished",'\n')
  b
} 
Likelihood_loss = function(pheno,t,par){
  miu = Logistic(par,t)
  L0 = c()
  L0 = sum((apply(pheno, 2, mean)-miu)^2)
  return(L0)
}
optim_loss <- function(pheno,t){
  library(mvtnorm)
  itime <- 100
  itimes <- 1
  par0 <- c(16.19012753, -28.48257768,  41.94282944,  22.09364587 ,  0.02529833 ) ##Stress
  repeat{
    a <- optim(par=par0,Likelihood_loss,pheno=pheno,t=t,control=list(maxit=20000),method = "BFGS")
    
    b <- optim(a$par,Likelihood_loss,pheno=pheno,t=t,control=list(maxit=20000),method = "BFGS")
    
    cat("Logistic",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) == 0 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  b
  
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,-1],2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),0.1,0.1)
}  

###FUNMAP_diff
#t需要与实际数据时间点一一对应
Logistic_diff <- function(par,t){
  
  y <- (par[1]/(1+par[2]*exp(-par[3]*t))-(par[4]*exp(-par[5]*t))) - (par[6]/(1+par[7]*exp(-par[8]*t))-(par[9]*exp(-par[10]*t)))
  return(y)
  
} 
SAD1_diff <- function(par, times = t) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi_ck<- par[1]   
  v2_ck <- par[2]
  phi_salt<- par[3]   
  v2_salt <- par[4] 
  
  tmp_ck <- (1-phi_ck^2)   
  sigma_ck <- array(1, dim=c(n,n))   
  
  tmp_salt <- (1-phi_salt^2)   
  sigma_salt <- array(1, dim=c(n,n)) 
  
  
  for(i in 1:n)   
  {     
    sigma_ck[i,i:n] <- phi_ck^( c(i:n) - i ) * (1-phi_ck^(2*i))/tmp_ck     
    sigma_ck[i:n,i] <- sigma_ck[i,i:n]   
    
    sigma_salt[i,i:n] <- phi_salt^( c(i:n) - i ) * (1-phi_salt^(2*i))/tmp_salt     
    sigma_salt[i:n,i] <- sigma_salt[i,i:n] 
    
  }   
  sigma_ck <- sigma_ck * abs(v2_ck)  
  sigma_salt <- sigma_salt * abs(v2_salt)
  sigma <- sigma_ck + sigma_salt
  
  return(sigma); 
}

Likelihood_diff = function(pheno,t,par){
  miu = Logistic_diff(par,t)
  sigma = SAD1_diff(par = c(par[11],par[12],par[13],par[14]),times = t)
  L0 = c()
  L0 = sum(dmvnorm(pheno,miu,sigma,log = T))
  return(-L0)
}##pheno每一行为一个系号的表型数据
Likelihood_diff_H1_2 = function(pheno_AA,pheno_aa,t,par){
  
  miu1 = Logistic_diff(par[1:10],t)
  miu1.1 = Logistic_diff(par[11:20],t)
  
  sigma = SAD1_diff(par = c(par[21],par[22],par[23],par[24]),times = t)
  L1 = c()
  L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
  L1.1 = c()
  L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
  
  L0 <- L1 + L1.1 
  
  return(-L0)
}
Likelihood_diff_H1_3 = function(pheno_AA,pheno_aa,pheno_Aa,t,par){
  
  miu1 = Logistic_diff(par[1:10],t)
  miu1.1 = Logistic_diff(par[11:20],t)
  miu1.2 = Logistic_diff(par[21:30],t)
  
  sigma = SAD1_diff(par = c(par[31],par[32],par[33],par[34]),times = t)
  L1 = c()
  L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
  L1.1 = c()
  L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
  L1.2 = c()
  L1.2 = sum(dmvnorm(pheno_Aa,miu1.2,sigma,log = T))
  
  L0 <- L1 + L1.1 + L1.2
  return(-L0)
}

test_par0 <- function(pheno_diff,par0,t){
  
  par1 <- optim_diff(pheno_diff[,-1],par0,t)
  par_test <<- par1$par
  aa <- c()
  bb <- c()
  for (n in 1:dim(pheno_diff)[1]) {
    aa <- cbind(t,as.numeric(pheno_diff[n,-1]),rep(n,dim(pheno_diff[,-1])[2]))
    bb <- rbind(bb,aa)
  }
  bb <- as.data.frame(bb)
  colnames(bb) <- c("t","value","group")
  
  
  library(ggplot2)
  cc <- ggplot() + geom_line(data = bb,aes(t,value,group=group),cex=0.7,col="skyblue")
  cc <- cc + geom_line(aes(t,Logistic_diff(par1$par[1:10],t)),cex=1,col="blue")
  cc <- cc + geom_line(aes(t,apply(pheno_diff[,-1], 2, mean)),cex=1,col="red")
  cc <- cc + scale_x_continuous(limits = c(min(t),max(t)),breaks = t,labels = t,expand = c(0,0) ) + xlab("") 
  cc <- cc + theme_bw()
  cc <- cc + 
    theme(
      axis.title.x =element_text(size=14), 
      axis.title.y=element_text(size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  cc
}
test_par1 <- function(pheno_diff,par0,t){
  
  optim_diff_BFGS <- function(pheno_diff,par0,t){
    library(mvtnorm)
    itime <- 20
    itimes <- 1
    
    repeat{
      a <- optim(par=par0,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000),method = "BFGS")
      
      b <- optim(a$par,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000),method = "BFGS")
      
      #cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
        break
      }else{
        par0 <- b$par
      }
    }
    cat("H0 Finished",'\n')
    b
  } 
  
  par1 <- optim_diff_BFGS(pheno_diff[,-1],par0,t)
  par_test <<- par1$par
  aa <- c()
  bb <- c()
  for (n in 1:dim(pheno_diff)[1]) {
    aa <- cbind(t,as.numeric(pheno_diff[n,-1]),rep(n,dim(pheno_diff[,-1])[2]))
    bb <- rbind(bb,aa)
  }
  bb <- as.data.frame(bb)
  colnames(bb) <- c("t","value","group")
  
  
  library(ggplot2)
  cc <- ggplot() + geom_line(data = bb,aes(t,value,group=group),cex=0.7,col="skyblue")
  cc <- cc + geom_line(aes(t,Logistic_diff(par1$par[1:10],t)),cex=1,col="blue")
  cc <- cc + geom_line(aes(t,apply(pheno_diff[,-1], 2, mean)),cex=1,col="red")
  cc <- cc + scale_x_continuous(limits = c(min(t),max(t)),breaks = t,labels = t,expand = c(0,0) ) + xlab("") 
  cc <- cc + theme_bw()
  cc <- cc + 
    theme(
      axis.title.x =element_text(size=14), 
      axis.title.y=element_text(size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  cc
}

optim_diff <- function(pheno_diff,par0,t){
  library(mvtnorm)
  itime <- 20
  itimes <- 1
  
  repeat{
    a <- optim(par=par0,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H0 Finished",'\n')
  b
}  ##H0情况
optim_diff_H1_2 <- function(pheno_AA_diff,pheno_aa_diff,par0,t){
  
  itime <- 20
  itimes <- 1
  
  par0 <- c(par0[1:10],par0[1:10],par0[11:14])
  repeat{
    a <- optim(par=par0,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H1_1 Finished",'\n')
  b
  
}  ##两基因型
optim_diff_H1_3 <- function(pheno_AA_diff,pheno_aa_diff,pheno_Aa_diff,par0,t){
  
  itime <- 10
  itimes <- 1
  par0 <- c(par0 <- c(par0[1:10],par0[1:10],par0[1:10],par0[11:14]))
  repeat{
    a <- optim(par=par0,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H1_3 Finished",'\n')
  b
  
}  ##三基因型

FunMap_Diff <- function(pheno_diff,marker,par0,t){
  
  library(mvtnorm)
  diff_LR <- c() 
  diff_par <- c()
  T_marker <- marker[,colnames(marker)%in%c(pheno_diff[,1])]
  for (a in 1:dim(marker)[1]) {
    
    AA <- as.numeric(names(which(T_marker[a,]==1)))
    aa <- as.numeric(names(which(T_marker[a,]==0)))
    Aa <- as.numeric(names(which(T_marker[a,]==2)))
    all <- as.numeric(names(which(T_marker[a,]!=9)))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    all_diff <- pheno_diff[pheno_diff[,1]%in%all,]
    
    AA_diff <- pheno_diff[pheno_diff[,1]%in%AA,]
    
    aa_diff <- pheno_diff[pheno_diff[,1]%in%aa,]
    
    optim_all_diff <- optim_diff(pheno_diff=all_diff[,-1],par0,t)
    
    
    if(NAa==0){ 
      
      optim_AA_aa_diff <- optim_diff_H1_2(AA_diff[,-1],aa_diff[,-1],par0,t)
      
      LR_diff <- -2*(abs(optim_all_diff$value) - abs(optim_AA_aa_diff$value))
      par_diff <- c(optim_all_diff$par,optim_AA_aa_diff$par,rep(NA,10))
      
      
    } else{
      
      Aa_diff <- pheno_diff[pheno_diff[,1]%in%Aa,]
      
      optim_AA_aa_Aa_diff <- optim_diff_H1_3(AA_diff[,-1],aa_diff[,-1],Aa_diff[,-1],par0,t)
      
      LR_diff <- -2*(abs(optim_all_diff$value) - abs(optim_AA_aa_Aa_diff$value))
      par_diff <- c(optim_all_diff$par,optim_AA_aa_Aa_diff$par)
    }
    diff_LR <- c(diff_LR,-LR_diff)
    diff_par <- rbind(diff_par,par_diff)
    cat("SNP",a,"finished","\n")
    
  }
  ##参数结果，两基因型1：14为H0情况下10个差值逻辑斯蒂参数+4个矩阵参数，15：24为AA差值逻辑斯蒂参数，
  ##                25：34为aa差值逻辑斯蒂参数，最后35：38为H1情况下矩阵参数。
  
  ###参数结果，三基因型1：14为H0情况下10个差值逻辑斯蒂参数+4个矩阵参数，15：24为AA差值逻辑斯蒂参数，
  ###                25：34为aa差值逻辑斯蒂参数，35：44为Aa情况下差值逻辑斯蒂参数，最后45：48为H1-1情况下矩阵参数。
  
  return(list(diff_LR,diff_par) ) 
} ##FM基本运算
FunMap_Diff_parallel <- function(pheno_diff,marker,markerColName,par0,t,method){
  
  Logistic_diff <- function(par,t){
    
    y <- (par[1]/(1+par[2]*exp(-par[3]*t))-(par[4]*exp(-par[5]*t))) - (par[6]/(1+par[7]*exp(-par[8]*t))-(par[9]*exp(-par[10]*t)))
    return(y)
    
  } 
  SAD1_diff <- function(par, times = t) {   
    n <- ifelse (is.vector(times), length(times), NCOL(times) )   
    phi_ck<- par[1]   
    v2_ck <- par[2]
    phi_salt<- par[3]   
    v2_salt <- par[4] 
    
    tmp_ck <- (1-phi_ck^2)   
    sigma_ck <- array(1, dim=c(n,n))   
    
    tmp_salt <- (1-phi_salt^2)   
    sigma_salt <- array(1, dim=c(n,n)) 
    
    
    for(i in 1:n)   
    {     
      sigma_ck[i,i:n] <- phi_ck^( c(i:n) - i ) * (1-phi_ck^(2*i))/tmp_ck     
      sigma_ck[i:n,i] <- sigma_ck[i,i:n]   
      
      sigma_salt[i,i:n] <- phi_salt^( c(i:n) - i ) * (1-phi_salt^(2*i))/tmp_salt     
      sigma_salt[i:n,i] <- sigma_salt[i,i:n] 
      
    }   
    sigma_ck <- sigma_ck * abs(v2_ck)  
    sigma_salt <- sigma_salt * abs(v2_salt)
    sigma <- sigma_ck + sigma_salt
    
    return(sigma); 
  }
  
  Likelihood_diff = function(pheno,t,par){
    miu = Logistic_diff(par,t)
    sigma = SAD1_diff(par = c(par[11],par[12],par[13],par[14]),times = t)
    L0 = c()
    L0 = sum(dmvnorm(pheno,miu,sigma,log = T))
    return(-L0)
  }##pheno每一行为一个系号的表型数据
  Likelihood_diff_H1_2 = function(pheno_AA,pheno_aa,t,par){
    
    miu1 = Logistic_diff(par[1:10],t)
    miu1.1 = Logistic_diff(par[11:20],t)
    
    sigma = SAD1_diff(par = c(par[21],par[22],par[23],par[24]),times = t)
    L1 = c()
    L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
    L1.1 = c()
    L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
    
    L0 <- L1 + L1.1 
    
    return(-L0)
  }
  Likelihood_diff_H1_3 = function(pheno_AA,pheno_aa,pheno_Aa,t,par){
    
    miu1 = Logistic_diff(par[1:10],t)
    miu1.1 = Logistic_diff(par[11:20],t)
    miu1.2 = Logistic_diff(par[21:30],t)
    
    sigma = SAD1_diff(par = c(par[31],par[32],par[33],par[34]),times = t)
    L1 = c()
    L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
    L1.1 = c()
    L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
    L1.2 = c()
    L1.2 = sum(dmvnorm(pheno_Aa,miu1.2,sigma,log = T))
    
    L0 <- L1 + L1.1 + L1.2
    return(-L0)
  }
  
  test_par0 <- function(pheno_diff,par0,t){
    
    par1 <- optim_diff(pheno_diff[,-1],par0,t)
    par_test <<- par1$par
    aa <- c()
    bb <- c()
    for (n in 1:dim(pheno_diff)[1]) {
      aa <- cbind(t,as.numeric(pheno_diff[n,-1]),rep(n,dim(pheno_diff[,-1])[2]))
      bb <- rbind(bb,aa)
    }
    bb <- as.data.frame(bb)
    colnames(bb) <- c("t","value","group")
    
    
    library(ggplot2)
    cc <- ggplot() + geom_line(data = bb,aes(t,value,group=group),cex=0.7,col="skyblue")
    cc <- cc + geom_line(aes(t,Logistic_diff(par1$par[1:10],t)),cex=1,col="blue")
    cc <- cc + geom_line(aes(t,apply(pheno_diff[,-1], 2, mean)),cex=1,col="red")
    cc <- cc + scale_x_continuous(limits = c(min(t),max(t)),breaks = t,labels = t,expand = c(0,0) ) + xlab("") 
    cc <- cc + theme_bw()
    cc <- cc + 
      theme(
        axis.title.x =element_text(size=14), 
        axis.title.y=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    cc
  }
  optim_diff <- function(pheno_diff,par0,t,method){
    library(mvtnorm)
    itime <- 20
    itimes <- 1
    
    repeat{
      a <- optim(par=par0,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000),method = method)
      
      b <- optim(a$par,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000),method = method)
      
      #cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
        break
      }else{
        par0 <- b$par
      }
    }
    cat("H0 Finished",'\n')
    b
  }  ##H0情况
  optim_diff_H1_2 <- function(pheno_AA_diff,pheno_aa_diff,par0,t,method){
    
    itime <- 20
    itimes <- 1
    
    par0 <- c(par0[1:10],par0[1:10],par0[11:14])
    repeat{
      a <- optim(par=par0,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000),method = method)
      
      b <- optim(a$par,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000),method = method)
      
      #cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
        break
      }else{
        par0 <- b$par
      }
    }
    cat("H1_1 Finished",'\n')
    b
    
  }  ##两基因型
  optim_diff_H1_3 <- function(pheno_AA_diff,pheno_aa_diff,pheno_Aa_diff,par0,t,method){
    
    itime <- 10
    itimes <- 1
    par0 <- c(par0 <- c(par0[1:10],par0[1:10],par0[1:10],par0[11:14]))
    repeat{
      a <- optim(par=par0,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000),method = method)
      
      b <- optim(a$par,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000),method = method)
      
      #cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
        break
      }else{
        par0 <- b$par
      }
    }
    cat("H1_3 Finished",'\n')
    b
    
  }  ##三基因型
  library(parallel)
  library(mvtnorm)
  diff_LR <- c() 
  diff_par <- c()
  T_marker <- marker[markerColName%in%c(pheno_diff[,1])]
  
  AA <- as.numeric(names(which(T_marker==1)))
  aa <- as.numeric(names(which(T_marker==0)))
  Aa <- as.numeric(names(which(T_marker==2)))
  all <- as.numeric(names(which(T_marker!=9)))
  
  NAA <- length(AA)
  Naa <- length(aa)
  NAa <- length(Aa)
  
  all_diff <- pheno_diff[pheno_diff[,1]%in%all,]
  
  AA_diff <- pheno_diff[pheno_diff[,1]%in%AA,]
  
  aa_diff <- pheno_diff[pheno_diff[,1]%in%aa,]
  
  optim_all_diff <- optim_diff(pheno_diff=all_diff[,-1],par0,t,method)
  
  
  if(NAa==0){ 
    
    optim_AA_aa_diff <- optim_diff_H1_2(AA_diff[,-1],aa_diff[,-1],par0,t,method)
    
    LR_diff <- -2*(abs(optim_all_diff$value) - abs(optim_AA_aa_diff$value))
    par_diff <- c(optim_all_diff$par,optim_AA_aa_diff$par,rep(NA,10))
    
    
  } else{
    
    Aa_diff <- pheno_diff[pheno_diff[,1]%in%Aa,]
    
    optim_AA_aa_Aa_diff <- optim_diff_H1_3(AA_diff[,-1],aa_diff[,-1],Aa_diff[,-1],par0,t,method)
    
    LR_diff <- -2*(abs(optim_all_diff$value) - abs(optim_AA_aa_Aa_diff$value))
    par_diff <- c(optim_all_diff$par,optim_AA_aa_Aa_diff$par)
  }
  diff_LR <- c(diff_LR,-LR_diff)
  diff_par <- rbind(diff_par,par_diff)
  
  
  
  ##参数结果，两基因型1：14为H0情况下10个差值逻辑斯蒂参数+4个矩阵参数，15：24为AA差值逻辑斯蒂参数，
  ##                25：34为aa差值逻辑斯蒂参数，最后35：38为H1情况下矩阵参数。
  
  ###参数结果，三基因型1：14为H0情况下10个差值逻辑斯蒂参数+4个矩阵参数，15：24为AA差值逻辑斯蒂参数，
  ###                25：34为aa差值逻辑斯蒂参数，35：44为Aa情况下差值逻辑斯蒂参数，最后45：48为H1-1情况下矩阵参数。
  
  return(list(diff_LR,diff_par)) 
} ##并行运算FM
out_parallel_deal <- function(out_parallel,marker_name){
  
  aa <- c()
  bb <- c()
  cc <- c()
  dd <- c()
  nn <- length(out_parallel)
  for (n in 1:nn) {
    
    aa <- out_parallel[[n]][[1]]
    cc <- out_parallel[[n]][[2]]
    bb <- c(bb,aa)
    dd <- rbind(dd,cc)
    
  }
  
  ee <- manhat_data(abs(bb),marker_name)
  
  return(list(bb,dd,ee))
  
}  ##将并行运算结果整理成标准格式
manhat_data <- function(FM_LR,marker_name){
  
  chr1 <- c(1:1000);chr2 <- c(1001:1634);chr3 <- c(1635:2205);chr4 <- c(2206:2773);chr5 <- c(2774:3289);chr6 <- c(3290:3779)
  chr7 <- c(3780:4229);chr8 <- c(4230:4633);chr9 <- c(4634:5019);chr10 <- c(5020:5391);chr11 <- c(5392:5749);chr12 <- c(5750:6120)
  chr13 <- c(6121:6492);chr14 <- c(6493:6849);chr15 <- c(6850:7185);chr16 <- c(7186:7488);chr17 <- c(7489:7783);chr18 <- c(7784:8072)
  chr19 <- c(8073:8305)
  
  ee <- marker_name[,3]
  ff <- as.integer(gsub("[^0-9]","",marker_name[,1]))
  gg <- c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19)
  hh <- abs(FM_LR)
  ii <- data.frame(ee,ff,gg,hh)
  colnames(ii) <- c("SNP","CHR","BP","P")
  ii
}
#####?????????????????????????????????
get_VG <- function(FunMap_par,pheno_diff,marker_data,t){
  
  T_marker <- marker_data[,colnames(marker_data)%in%c(pheno_diff[,1])]
  
  diff_vg <- c() 
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(names(which(T_marker[a,]==1)))
    aa <- as.numeric(names(which(T_marker[a,]==0)))
    Aa <- as.numeric(names(which(T_marker[a,]==2)))
    all <- as.numeric(names(which(T_marker[a,]!=9)))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A基因频率
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a基因频率
    
    mean_AA <- Logistic_diff(FunMap_par[a,][15:24],t)
    mean_aa <- Logistic_diff(FunMap_par[a,][25:34],t)
    AE <- (mean_AA - mean_aa)
    
    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{
      
      mean_Aa <- Logistic_diff(FunMap_par[a,][35:44],t)
      AE <- (mean_AA - mean_aa)/2
      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE
      
    }
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  rownames(diff_vg) <- c(1:dim(marker_data)[1])
  return(sqrt(diff_vg)) 
}###根据求得的参数计算遗传标准差

get_VG_Mean <- function(pheno_diff,marker_data,t){
  
  T_marker <- marker_data[,colnames(marker_data)%in%c(pheno_diff[,1])]
  
  diff_vg <- c() 
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(names(which(T_marker[a,]==1)))
    aa <- as.numeric(names(which(T_marker[a,]==0)))
    Aa <- as.numeric(names(which(T_marker[a,]==2)))
    all <- as.numeric(names(which(T_marker[a,]!=9)))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A基因频率
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a基因频率
    
    mean_AA <- colMeans(pheno_diff[pheno_diff[,1]%in%AA,-1])
    mean_aa <- colMeans(pheno_diff[pheno_diff[,1]%in%aa,-1])
    AE <- (mean_AA - mean_aa)
    
    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{
      
      mean_Aa <- colMeans(pheno_diff[pheno_diff[,1]%in%Aa,-1])
      AE <- (mean_AA - mean_aa)/2
      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE
      
    }
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  rownames(diff_vg) <- c(1:dim(marker_data)[1])
  return(sqrt(diff_vg)) 
}###根据求得的参数计算遗传标准差

get_cluster <- function(cutree,cluster_num,plast_Value){
  
  allcluster <- list()    
  for (a in 1:cluster_num) {
    
    allcluster[[a]] <- as.data.frame(plast_Value)[which(cutree==c(1:cluster_num)[a]),]
    
  }
  allcluster
} #将每一类的样本表型数据整合在一起

##变量选择与微分方程用
Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}
dLegendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1]*0 + 1*mu[2];
  if (np.order>=2)
    L <- L + 0.5 * (6 * ti)* mu[3] ;
  if (np.order>=3)
    L <- L +0.5 * (15 * ti ^ 2 - 3)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)* mu[5];
  if (np.order>=5)
    L <- L + 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)*mu[6];
  if (np.order>=6)
    L <- L + (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *ti)* mu[7];
  if (np.order>=7)
    L <- L + (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)* mu[8];
  return(L);
}
Likehood_Legendre_ind <- function(times,para,E){
  
  sum( (E -Legendre.model(t=times,mu=para))^2)
  
}
Likehood_Legendre_cluster <- function(times,para,E){
  
  sum((apply(E, 2, mean)-Legendre.model(t=times,mu=para))^2)
  
} #类,先平均再最小二乘
smooth.optim_ind <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- optim(para,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS")
  L <- optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS")
  L <- optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method = "BFGS")
  mean_par <- L$par
  mean_par
}
smooth.optim_cluster <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- optim(para,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  mean_par <- L$par
  mean_par
}
get_cluster <- function(cutree,cluster_num,plast_Value){
  
  allcluster <- list()    
  for (a in 1:cluster_num) {
    
    allcluster[[a]] <- as.data.frame(plast_Value)[which(cutree==c(1:cluster_num)[a]),]
    
  }
  allcluster
} #将每一类的样本表型数据整合在一起
varsel1 <- function(X,Y,tt,order=6){
  
  nobs = nrow(X)
  ndim = ncol(X)
  dfo = rep(order-1,ndim)
  index = rep(1:ndim,times=dfo)
  aa2 <- c()
  for(k in 1:ndim){
    aa1 <- c()
    for(j in 1:(order-1)){
      aa <- c()
      for(i in 1:length(tt)){
        aa <- c(aa,Legendre.modelve(tt[i],np.order = j,tmin = min(tt), tmax = max(tt)))
      }
      aa1 <- cbind(aa1,aa*X[,k])
    }
    aa2 <- cbind(aa2,aa1)
  }
  
  
  Xc = scale(aa2,center=T,scale=T)
  n = nrow(Xc)
  
  connect = matrix(0,nrow=ndim,ncol=ndim)
  coefest = matrix(0,nrow=sum(dfo),ncol=ndim)
  regfun = vector("list",length=ndim)
  for(i in 1:ndim)
  {
    inde <- index[-which(index==i)]
    XXc <- Xc[,-which(index==i)]
    
    yc <- Y[,i]-mean(Y[,i])
    
    out1 <- GrpLasso(X=XXc,y=yc,index=inde,lambda=30,crit="BIC")
    var.grp <- out1$var.select  # genes selected
    coef.grp <- out1$coef
    
    ### Adaptive Group Lasso
    index.adp <- index[is.element(index,var.grp)]
    W.adp = sapply(1:length(var.grp),function(j) sqrt(sum(coef.grp[index.adp==var.grp[j]]^2)))
    Xc.adp = Xc[,is.element(index,var.grp)]
    Xcs.adp = scale(Xc.adp,center=F,scale=rep(1/W.adp,times=dfo[var.grp]))
    init.adp = coef.grp/rep(W.adp,times=dfo[var.grp])
    lambda = lambdamax(Xcs.adp,yc,index=index.adp,coef.ini=init.adp,
                       penscale=sqrt,center=F,standardize=F,model=LinReg())*0.7^(1:18)
    out2 = GrpLasso(X=Xc.adp,y=yc,W=W.adp,index=index.adp,ini=coef.grp,
                    lambda=lambda,crit="BIC")
    var.adp = out2$var.select
    coef.adp = out2$coef
    connect[i,var.adp] <-  1
    coefest[is.element(index,var.adp),i] <-  coef.adp
    regfun[[i]] <-  sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    connect[i,var.adp] = 1
    coefest[is.element(index,var.adp),i] = coef.adp
    regfun[[i]] = sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    cat("var=",i,var.adp,"\n")
  }
  return(list(connect=connect,regfun=regfun,coefest=coefest))
}
Legendre.modelve <- function(t, np.order, tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- NA
  L <- 1;
  if (np.order >= 2)
    L <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L <- 0.5 * (15 * ti ^ 2 - 3)
  if (np.order >= 4)
    L <- 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)
  if (np.order >= 5)
    L <- 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L <- (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                       ti)
  if (np.order >= 7)
    L <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                       ti ^ 2 - 35) 
  return(L);
}
GrpLasso <- function(X,y,W=NULL,index,ini=rep(0,ncol(X)),lambda=NULL,crit=c("BIC","EBIC"),center=F){
  if(center==T){
    y = y-mean(y)
    X = scale(X,center=T,scale=F)
  }
  n = nrow(X)
  ind = unique(index)
  p = length(ind)
  dfo = sapply(1:p,function(j) sum(index==ind[j]))
  
  # fit model for a sequence of penalty parameters
  if(!is.null(W)){
    X = scale(X,center=F,scale=rep(1/W,times=dfo))
    ini = ini/rep(W,times=dfo)
  }
  
  # set up the candidates for penalty parameter
  if(is.null(lambda)){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:20)
  }else if(length(lambda)==1){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:lambda)
  }
  
  fit = grplasso(X,y,index=index,lambda=lambda,model=LinReg(),center=F,
                 standardize=F,coef.ini=ini,penscale=sqrt,
                 control=grpl.control(update.hess="lambda",tol=10^-8,trace=0))
  
  # calculate BIC/EBIC
  nlambda = length(lambda)
  rss = sapply(1:nlambda,function(j) sum((y-fit$fitted[,j])^2))
  var.select = sapply(1:nlambda,function(j) unique(index[abs(fit$coef[,j])>0]))
  dfo.lambda = sapply(1:nlambda,function(j) sum(dfo[is.element(ind,var.select[[j]])]))
  if(crit!="BIC" & crit!="EBIC"){
    cat("Error: Criterion not implemented. Reset to BIC!\n")
    crit = "BIC"
  }
  if(crit=="BIC"){
    bic = log(rss)+dfo.lambda*log(n)/n
  }else if(crit == "EBIC"){
    bic = log(rss)+dfo.lambda*log(n)/n+0.5*dfo.lambda*log(p)/n
  }
  
  # select model with smallest value of selection criterion
  id.ss = which.min(bic)
  var.ss = var.select[[id.ss]]
  fit.ss = fit$fitted[,id.ss]
  coef.ss = fit$coef[,id.ss]
  if(!is.null(W)){
    coef.ss = coef.ss*rep(W,times=dfo)
  }
  coef.ss = coef.ss[is.element(index,var.ss)]
  
  return(list(var.select=var.ss,coefficients=coef.ss,fitted=fit.ss,BIC=bic,
              lambda=lambda,fit.org=fit))
}
optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){
  
  diag(connect) <- 1
  
  LL <- L_derivative(nt=times,nstep=nstep,order=order)
  
  nx <- dim(effect)[2]
  
  grp <- floor(nx/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nx))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nx-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }
  
  grp.ii <- unique(grp.i)
  
  res.list <- mclapply(grp.ii, function(i)
  {
    y.c <- 	which(grp.i==i)
    A <- sapply(y.c, proc, connect=connect,effect=effect,LL=LL,nstep=nstep,order=order,times=times);
    return (unlist(A));
  }, mc.cores=n.cores )
  
  res1 <- do.call("c", res.list)
  res2 <- parallel.data.optim(res1,connect,times)
  return(res2)
}
L_derivative <- function(nt,nstep,order){
  
  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){
    
    tg1 <- Legendre.model1((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model1(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}
Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- 1;
  if (np.order >= 2)
    L[2] <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L[3] <- 0.5 * (15 * ti ^ 2 - 3) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * 4 * ti ^ 3 - 60 * ti) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                         ti) 
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                          ti ^ 2 - 35)
  return(L);
}
Legendre.model1 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- ti;
  if (np.order >= 2)
    L[2] <- 0.5 * (3 * ti^2 - 1) 
  if (np.order >= 3)
    L[3] <- 0.5 * (5 * ti ^ 3 - 3 * ti) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * ti ^ 4 - 30 * ti^2 + 3) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * ti ^ 5 - 70 * ti ^ 3 + 15*ti)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * ti ^ 6 - 315 * ti ^ 4 + 105 * ti^2 -5) 
  
  return(L);
}
ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  self <- y.c 
  indexx <- which(connect[y.c,]==1)
  #para <- rep(0.00001,length(indexx)*(order-1))
  para <- rep(0.001,length(indexx)*(order-1))
  #res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
  #            LL=LL,method="BFGS",control=list(maxit=2000,trace=T))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,control=list(maxit=4000,trace=F))
  cat("Gene=",y.c," ",res$value,"\n")
  A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL,self=self)
  return(A)
}
fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL,self = self)
  sum((NG[,self]-rowSums(odes))^2) ##最小二乘
  
}
ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL,self){
  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)
  
  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep <- matrix(rep(0,length(index)),nrow=1)
 # allrep[which(index==self)] <- NG[1,self]
  nn <- 1
  for(j in 1:nstep){
    tg1 <- (rowSums(t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index])
    tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] +stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
  
  
  self_name <- colnames(NG)[self]
  no_name <- which(colnames(allrep)==self_name)
  allrep[,no_name] <- allrep[,no_name]
  allrep
}

parallel.data.optim <- function(rd,nm,ntt){
  
  nrd <- matrix(rd,nrow=length(ntt))
  nn <- dim(nm)[1]
  ki <- 0
  allist <- list()
  for(i in 1:nn){
    iii <- (which(nm[i,]==1))
    iiil <- length(iii)
    tmp.d <- nrd[,(ki+1):(ki+iiil)]
    if(is.matrix(tmp.d)){
      colnames(tmp.d) <- iii
    }else{
      names(tmp.d) <- iii
    }
    
    allist[[i]] <- tmp.d
    ki <- ki + iiil
  }
  
  return(allist)
}
get_eageInform <- function(eage,cluster){
  
  eage_sati <- c()
  for (n in rownames(cluster)) {
    
    
    if(any(which(eage[,1]==n))){
      a1 <- eage[which(eage[,1]==n),]
      a1.0 <- 0
      a1.1 <- 0
      if(!is.matrix(a1)){
        if(a1[3]==0){
          a1.0 <- 1
          a1.1 <- 0}
        else{
          a1.0 <- 0
          a1.1 <- 1}}
      else{
        for (q in 1:dim(a1)[1]) {
          a2 <- as.matrix(a1)[q,]
          if(a2[3]==0){
            a1.0 <- a1.0+1
            a1.1 <- a1.1+0
          }
          else{
            a1.0 <- a1.0+0
            a1.1 <- a1.1+1
          }
        }
      }
      a3 <- c(a1.0,a1.1)
      
    }
    else{a3 <- c(0,0)}
    
    
    b1.0 <- 0
    b1.1 <- 0
    b1 <- eage[which(eage[,2]==n),]
    if(!is.matrix(b1)){
      if(b1[3]==0){
        b1.0 <- 1
        b1.1 <- 0
      }
      else{
        b1.0 <- 0
        b1.1 <- 1
      }
    }
    else{
      for (w in 1:dim(b1)[1]) {
        b2 <- b1[w,]
        if(b2[3]==0){
          b1.0 <- b1.0+1
          b1.1 <- b1.1+0
        }
        else{
          b1.0 <- b1.0+0
          b1.1 <- b1.1+1}
      }
    }
    b3 <- c(b1.0,b1.1)
    cc <- c(a3,b3)
    eage_sati <- rbind(eage_sati,cc)
    
  }
  eage_sati <- cbind(rownames(cluster),eage_sati)
  colnames(eage_sati) <- c("no","sour_con","sour_pro","tar_con","tar_pro")
  rownames(eage_sati) <- c(1:dim(cluster)[1])
  eage_sati
  
  
}
get_integrate <- function(out,t){
  s <- t[2] - t[1]
  c <- c()
  for (n in 1:length(out[1,])) {
    a <- out[,n]
    b <- c()
    for (i in 1:(length(out[,1])-1)) {
      b <- c(b,a[i] * s)
    }
    c <- cbind(c,b)
  }
  colnames(c) <- colnames(out)
  c/100
}
regasso_moudule <- function(connect1,gene,interaction){
  aaa <- colnames(connect1)
  diag(connect1) <- 0
  ng <- dim(gene)[1]
  allcor <- list()
  for(i in 1:ng){
    a1 <- gene[i,]
    nng1 <- (interaction[[i]])
    if(!is.matrix(nng1)){
      next
    }else{
      nng <- as.matrix(nng1[,-which(colnames(nng1)==aaa[i])])
      corr <- c()
      for(j in 1:dim(nng)[2]){
        corr <- c(corr,cor(as.numeric(a1),as.numeric(nng[,j])))
      }
      allcor[[i]] <- corr
    }
  }
  
  for (q in 1:length(allcor)) {
    connect1[q,which(connect1[q,]==1)] <- allcor[[q]]
  }
  
  return(connect1)
}


##解微分方程，求得SNP互作数据
optim_inter <- function(all_cluster_value,ind_Lpar11,norder){
  
  library(splines)
  library(orthogonalsplinebasis)
  library(MASS)
  library(grplasso)
  library(parallel)
  
  ind_Lpar1 <- ind_Lpar11
  clusterAA <- all_cluster_value
  
  
  ind_LparAA <- ind_Lpar1[as.numeric(rownames(clusterAA)),]
  rownames(ind_LparAA) <-  rownames(clusterAA) #调出先前求出的该类个体的勒让德拟合参数
  
  ttt <- seq(13,78,length.out = (dim(ind_LparAA)[1]+5))
  cAA_mean <- apply(ind_LparAA,1, Legendre.model,t=ttt)
  cAA_d <- apply(ind_LparAA,1, dLegendre.model,t=ttt)
  
  clusterAA_50 <- varsel1(X=cAA_mean,Y=cAA_d,tt=ttt)
  colnames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  rownames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  
  
  clusterAA_list <- list()
  for (n in 1:dim(ind_LparAA)[1]) {
    
    ttemp <- list()
    ttemp[[1]] <-  rownames(clusterAA)[n]
    ttemp[[2]] <- names(which(clusterAA_50$connect[n,]==1))
    clusterAA_list[[n]] <- ttemp
  }
  
  
  tryyAA <- optim.parallel(connect=clusterAA_50$connect,effect=t(all_cluster_value),
                           n.cores=1,proc=ode.optim,order=norder,times=seq(13,78,length=30),nstep=29)
  
  for (e in 1:length(tryyAA)) {
    
    colnames(tryyAA[[e]]) <- rownames(clusterAA)[as.numeric(colnames(tryyAA[[e]]))]
    
  } 
  
  
  
  ##计算相关系数
  inter_effect <- regasso_moudule(connect1=clusterAA_50$connect,gene=all_cluster_value,interaction=tryyAA)
  
  
  
  for (n in 1:length(clusterAA_list)) {
    
    
    nono <- which(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0])[order(as.numeric(names(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0]))))]==0)
    
    clusterAA_list[[n]][[3]] <- colSums(get_integrate(out=tryyAA[[n]],t=seq(13,78,length=30)))
    
    clusterAA_list[[n]][[3]][nono] <- 0
  } 
  
  
  
  ##assign(paste("optim_cluster_",which_cluster,"list",sep=""),clusterAA_list, envir = .GlobalEnv)
  optim_cluster <- list()
  optim_cluster[[1]] <- clusterAA_list
  optim_cluster[[2]] <- tryyAA
  optim_cluster
}


##根据上述求得的数据输出画图用的数据
out_Netdata_cytoscape <- function(optim_cluster,which_cluster,clusterAA){
  
  clusterAA_list <- optim_cluster[[1]]
  afterAA <- c()
  effect_selfAA <- c()
  for (i in 1:length(clusterAA_list)){
    dep <- clusterAA_list[[i]][[1]]
    ind <- clusterAA_list[[i]][[2]]
    effectAA <- clusterAA_list[[i]][[3]]
    effect_selfAA <- rbind(effect_selfAA,effectAA[which(names(effectAA)==dep)])
    effectAA <-effectAA[-which(names(effectAA)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effectAA[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effectAA[j]),type))
    }
    afterAA <- rbind(afterAA,one)
    
  }
  eageAA <- cbind(afterAA[,-3],c(0:(length(afterAA[,1])-1)),afterAA[,3])
  colnames(eageAA) <- c("source","target","colour","Id","Weight")
  
  new <- eageAA
  new[,1] <- paste("S",new[,1],sep = "")
  new[,2] <- paste("S",new[,2],sep = "")
  write.csv(new,paste("M",which_cluster,".csv",sep = ""), row.names = F,fileEncoding = "UTF-8")
  
  get_eageInform(eageAA,clusterAA) 
  
}





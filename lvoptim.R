lv.EQ.est <- function(expheno,renewgen,interval=c(1,100)){
  
  y1 <- as.matrix(expheno[,12:22]/100)
  y2 <- as.matrix(expheno[,45:55])
  
  light<- light2
  geno_table <- renewgen
  nm <- dim(geno_table)[1]
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nm)
    n2 <- nm
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=50)
  for(i in n1:n2){
    SNP <- geno_table[i,]
    missing <- which(is.na(SNP))
    if ( length(missing) > 0)
    {
      SNP1 <- SNP[-(missing)]
      y11 <- y1[ -(missing), ]
      y22 <- y2[ -(missing), ]
    }else{
      SNP1 <- SNP
      y11 <- y1
      y22 <- y2
    }
    
    
    h01 <- try(LV.EQ.H0(y11,y22),TRUE)
    if (class(h01) == "try-error") 
      h01 <- NA
    par1 <- h01[-1]
    h02 <- try(LV.EQ.H1(y11,y22,SNP1,par1=par1,times),TRUE)
    if (class(h02) == "try-error") 
      h02 <- NA
    LR <- 2*(h01[1]-h02[1])
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h02)
    }
    res[(i-(n1-1)),(1:length(allpar))] <- allpar
  }
  return(res)
}
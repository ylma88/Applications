rm(list=ls())
library(MASS)
library(splines)
library(magrittr) #or library(dplyr)
library(ismev)

#Application 1
data(fremantle)
SeaLevel = fremantle$SeaLevel
Year = fremantle$Year
SOI = fremantle$SOI
plot(Year,SeaLevel) #Annual maximum sea levels at Fremantle
plot(SOI,SeaLevel)

#Approximate ML estimate 
loglike=function(Z,n0,nknot){
  o = order(Z$SeaLevel,Z$SOI,Z$Year)
  sort_Z = Z[o,]
  y = sort_Z[,3]
  u_n = y[n-n0]
  Z.samp = sort_Z[(n-n0+1):n,]
  W1 = Z.samp[,1]
  X2 = Z.samp[,2]
  Y = Z.samp[,3]
  
  num.knot = seq(0,0.85,length.out=nknot)
  knots.1 = quantile(W1, num.knot) #节点个数,取为w1的样本分位数
  B.basis.1 = bs(W1,knots=knots.1,degree=D,intercept=FALSE) #B-spline basis function
  v = cbind(X2,B.basis.1)
  loglik = function(d){
    -sum(as.vector(v%*%d)-exp(as.vector(v%*%d))*log(Y/u_n))
  }
  res = optim(c(0.1,rep(0.01,dim(v)[2]-1)),loglik)
  U.hat = exp(-exp(v%*%res$par)*log(Y/u_n))
  Fn = ecdf(U.hat)
  D.hat = (1/n0)*sum((U.hat-Fn(U.hat))^2)
  #choose the number of knots
  #AIC = 2*res$value+2*(nknot+D+1)
  Gamma = (1/n0)*lapply(1:n0,function(j) as.matrix(v[j,])%*%t(as.matrix(v[j,])))%>%Reduce('+',.)
  return(list(par.est=res$par,log.value=res$value,D_hat=D.hat,Gamma=Gamma))
}

  #case1:
  n=dim(fremantle)[1];
  st=0.05*n;en=0.5*n;n0=seq(st,en,length.out=100);
  knot=3:7;m=35;D=3;

  X=scale(cbind(SOI,Year))
  Z=as.data.frame(cbind(X,SeaLevel))
  out.Dk=matrix(0,nrow=length(knot),ncol=length(n0))
  for(b in 1:length(knot))
  {
    nknot=knot[b]
    for(j in 1:length(n0))
    {
      e.size=n0[j] 
      out.Dk[b,j]=loglike(Z,e.size,nknot)$D_hat
    }
  }
  #Best n0  
  sub=as.vector(which(out.Dk==min(out.Dk),arr.ind=TRUE)) 
  nknot.opt=knot[sub[1]]
  n0.opt=n0[sub[2]]
  print(c(nknot.opt,n0.opt))
  output=loglike(Z,n0.opt,nknot.opt)
  Par = output$par.est
  Gamma.hat = output$Gamma
  par.est = Par[1]
  #estimate error
  np=1  #number of parameters
  Cov.est=solve(Gamma.hat)/n0.opt 
  var.beta.1=Cov.est[1,1]
  theta.Cov.1=Cov.est[(np+1):(nknot.opt+D+np),(np+1):(nknot.opt+D+np)]
  SE_0=qnorm(0.975)*sqrt(var.beta.1);SE_0
  
  #estimated function
  plotrange=c(round(min(SOI)),floor(max(SOI)))          
  t.1=seq(plotrange[1],plotrange[2],length.out=m)
  n.knot = seq(0,1,length.out=nknot.opt)
  eval.knots.1 = quantile(t.1, n.knot)
  bs.1 = bs(t.1, knots=eval.knots.1, degree=D, intercept=FALSE)
  nonpar.est.1=Par[-1]
  phi.1=as.vector(bs.1%*%nonpar.est.1)
  #estimate se
  phi.1.se=unlist(lapply(1:m,function(j) ((bs.1%*%theta.Cov.1)%*%t(bs.1))[j,j]))%>%as.vector()

  #Interval estimation 
  lower.1 = phi.1-qnorm(0.975)*sqrt(phi.1.se)
  upper.1 = phi.1+qnorm(0.975)*sqrt(phi.1.se)

  #fitted model
  plot(t.1,phi.1,type="l",ylim=c(-6, 6),lwd=2,lty=1,
       xlab="SOI",ylab="g^hat")
  lines(t.1,lower.1,type="l",col='red',lwd=2,lty=3)
  lines(t.1,upper.1,type="l",col='red',lwd=2,lty=3)
  
  
  
  

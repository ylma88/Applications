##model:log(alpha(x,w))=beta_1*SOI+g(t)
rm(list=ls())
library(ismev)
library(splines)
data(fremantle)
SeaLevel <- fremantle$SeaLevel
Year <- fremantle$Year
SOI <- fremantle$SOI
plot(Year,SeaLevel) #Annual maximum sea levels at Fremantle
plot(SOI,SeaLevel)  #Annual maximum sea levels versus SOI

m <- length(Year)
plotrange <- c(1,m)    
grid.num <- m  
T = seq(plotrange[1],plotrange[2],length.out=grid.num)
equa.knots <- c(T[1],T[18],T[36],T[54],T[72],T[86])
bs <- bs(T,knots=equa.knots,degree=D)
D <- 3 #cubic spline

W1 <- 1:m
X1 <- fremantle$SOI
Y1 <- fremantle$SeaLevel
n <- length(Y1)
z <- cbind(X1,W1,Y1) #covariate
o <- order(Y1,X1,W1)
sort_z <- z[o, ]
V <- sort_z[,3]

#Define k function 
k_fun <- function(k){
  u <- V[n-k]
  z.samp <- sort_z[(n-k+1):n,]
  Xi <- z.samp[,1]
  Wi <- z.samp[,2]
  Yi <- z.samp[,3]
  knots <- c(quantile(Wi,0.1),quantile(Wi,0.25),quantile(Wi,0.5), 
             quantile(Wi,0.65),quantile(Wi,0.8),quantile(Wi,0.95)) 
  Bsf <- bs(Wi,knots = knots,degree = D) #B-spline basis
  v <- cbind(Xi,Bsf) #covariate
  #Define the objective function
  loglik <- function(d){
    -sum(as.vector(v %*% d)-exp(as.vector(v%*%d))*log(Yi/u) )
    }
  result <- optim(c(0.1,rep(0.01,dim(v)[2]-1)),loglik)$par
  #sample fraction
  alpha_hat <- exp(v%*%result)
  U_hat <- exp(-alpha_hat*log(Yi/u))
  U_hat_k <- sort(U_hat)
  su <- rep(0,k)
  for(j in 1:k){
    su[j] <- (U_hat_k[j]-(j/(k+1)))^2
    }
  return(c(mean(su),result))
  }
k0 <- seq(5,45,by=1)
k_sf <- rep(0,length(k0))
  for(i in 1:length(k0)){
    kk <- i+4
    k_sf[i] <- k_fun(kk)[1]
  }
k_opt <- which.min(k_sf)+4

#use k_opt to estimate the beta and g
z.samp2 <- sort_z[(n-k_opt+1):n,]
X2 <- z.samp2[,1]
W2 <- z.samp2[,2]
knots_2 <- c(quantile(W2,0.1),quantile(W2,0.25),quantile(W2,0.5), 
             quantile(W2,0.65),quantile(W2,0.8),quantile(W2,0.95)) 
Bsf2 <- bs(W2,knots=knots_2,degree=D) 
v2 <- cbind(X2,Bsf2) #covariate
res <- k_fun(k_opt)[2:length(k_fun(k_opt))]
#estimate Omega
Omega <- matrix(0,length(res),length(res))
for(j in 1:k_opt){
    S <- matrix(,length(res),length(res))
    S <- as.matrix(v2[j,])%*%t(as.matrix(v2[j,]))
    Omega <- Omega+S
    }
Omega_hat <- Omega/k_opt   
Omega_inv <- solve(Omega_hat)/k_opt  
theta_inv <- Omega_inv[2:dim(Omega_inv)[2],2:dim(Omega_inv)[2]]  
g_SE <- rep(0,length(T))
for(j in 1:length(T)){
    g_SE[j] <- ((bs%*%theta_inv)%*%t(bs))[j,j] 
}
#estimate
g_hat <- as.vector(bs%*%res[2:length(res)]) 
lower <- g_hat-qnorm(0.975)*sqrt(g_SE)
upper <- g_hat+qnorm(0.975)*sqrt(g_SE)
#plot
plot(Year,g_hat,type="l",col='blue',lwd=2,ylim=c(-4,5),lty=1)
lines(Year,lower,type="l",col='red',lwd=2,lty=2)
lines(Year,upper,type="l",col='red',lwd=2,lty=2)
#standard error of beta
beta_hat <- res[1];beta_hat
SE_0 <- qnorm(0.975)*Omega_inv[1,1];SE_0




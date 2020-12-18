## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(stats)
data(Orange) 
summary(Orange)
aov.circumference <- aov(sqrt(age) ~ circumference, data = Orange)
summary(aov.circumference) 

## ----echo=T-------------------------------------------------------------------
  knitr::kable(head(Orange),caption="Orange trees Data")

## -----------------------------------------------------------------------------
opar <- par() 
plot(aov.circumference) 

termplot(aov.circumference, se=TRUE, partial.resid=TRUE, rug=TRUE)


## -----------------------------------------------------------------------------

set.seed(2342)
pareto<-function(a,b,n){
   u <- runif(n) 
   x <- b/(1-u)^{1/a}
   return(x)
}

a<-pareto(2,2,1000)
hist(a[a<15], breaks=100,prob = TRUE, xlab='x',main = 'Pareto(2,2)')
y <- seq(2, 15, .01)
lines(y, 8*y^{-3})


## -----------------------------------------------------------------------------
set.seed(1234)
kernel<-function(n){
   u1 <- runif(n,-1,1) 
   u2 <- runif(n,-1,1)
   u3 <- runif(n,-1,1)
   u=vector(mode="numeric",length=n)
   for(i in 1:n){
      if(abs(u3[i])<abs(u2[i])||abs(u3[i])<abs(u1[i])) u[i]=u3[i]
      else u[i]=u2[i]
   }
   return(u)
}

a<-kernel(2000)


hist(a, breaks=100,prob = TRUE, xlab='x',main = 'kernel')
y <- seq(-1, 1, .01)
lines(y, 3/4*(1-y^{2}))


## -----------------------------------------------------------------------------
set.seed(1645)
n <- 1e3
r <- 4
beta <- 2 
lambda <- rgamma(n, r, beta) 
y <- rexp(n, lambda)
hist(y[y<10], breaks=100,prob = TRUE, xlab='y',main = 'mixture')

t <- seq(0, 10, .01)
lines(t, 64*(t+2)^{-5})

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1e4; 
x <- runif(m, min=0, max=pi/3) 
theta.hat <- mean(sin(x)) * pi/3 
print(theta.hat)

## -----------------------------------------------------------------------------

set.seed(7912)

myanti<-function(n=1000, anti=TRUE){
   t<-runif(n/2,0,1)
   if (!anti) s<-runif(n/2,0,1) 
   else s<-1-t
   x<-c(t, s)
   g<-exp(x)
   mean(g)
}

m<-1000
mcint1=mcint2=numeric(m)
for(i in 1:m){
   mcint1[i]<-myanti(n=1000,anti=FALSE)
   mcint2[i]<-myanti(n=1000)
}

c(var(mcint1),var(mcint2),(var(mcint1)-var(mcint2))/var(mcint1)) 


## ----fig.width=10-------------------------------------------------------------
set.seed(22320)
n <- 1000
u <- runif(n)
x1 <- 1/(1-u) 

x2=rnorm(n,1,1)
x2[x2<1]=2-x2[x2<1]

a=exp(-x1*x1/2)*x1*x1/sqrt(2*pi)*x1*x1
b=exp(-x2*x2/2)*x2*x2/sqrt(2*pi)/(exp(-(x2-1)*(x2-1)/2)*sqrt(2/pi))
   

## ----fig.width=10-------------------------------------------------------------
c(mean(a),mean(b))

## ----fig.width=10-------------------------------------------------------------
c(sd(a),sd(b))

## ----fig.width=10-------------------------------------------------------------
x <- seq(1, 10, 0.1)
w <- 2
g <- exp(-x*x/2)*x*x/sqrt(2*pi)
f1 <- 1/(x*x)
f2 <- exp(-(x-1)*(x-1)/2)*sqrt(2/pi)
gs <- c(expression(g(x)==e^{-x^{2}/2}*x^{2}/sqrt(2*pi)),
        expression(f[1](x)==1/x^{2}),
        expression(f[2](x)==e^{-(x-1)^{2}/2}*sqrt(2/pi)))
#for color change lty to col
par(mfrow=c(1,2))
#figure (a)
plot(x, g, type = "l", ylab = "",
      ylim=c(0,1),lwd = w,col=1,main='pdf')
lines(x, f1, lty = 2, lwd = w,col=2)
lines(x, f2, lty = 3, lwd = w,col=3)

legend("topright", legend = gs,
       lty = 1:3, lwd = w, inset = 0.02,col=1:3)

#figure (b)
plot(x, g/f1, type = "l", ylab = "",
     ylim = c(0,1.0), lwd = w, lty = 2,col=2,main='g(x)/f(x)')
lines(x, g/f2, lty = 3, lwd = w,col=3)

legend("topright", legend = gs[-1],
       lty = 2:3, lwd = w, inset = 0.02,col=2:3)


## ----fig.width=10-------------------------------------------------------------
M=1000
k=5
r=M/k
N=100
theta<-numeric(N)
T=numeric(k)
g=function(x)exp(-x)/(1+x^2)
f=function(x)exp(-x)/(1-exp(-1))
a=numeric(k)
for(i in 1:N){
  for(j in 1:k) a[j]=exp(-(j-1)/5)-exp(-j/5)
  #a=(1-exp(-1))/5
  for(j in 1:k){
    u=runif(r,0,1)
    x=-log(exp(-(j-1)/5)-a[j]*u)
    T[j]=mean(a[j]*g(x)/(f(x)*(1-exp(-1))))
  }
  theta[i]=sum(T)
}

## ----fig.width=10-------------------------------------------------------------
c(mean(theta),sd(theta))

## ----fig.width=10-------------------------------------------------------------
mu <- 0
sigma <- 1 # Null hypothesis!
m<-1e4
n<-100
set.seed(123)
mu_hat<-sigma_hat<-pval<-numeric(m)
for(i in 1:m){
  x<-rnorm(n,mu,sqrt(sigma))
  mu_hat[i]<-mean(x)
  sigma_hat[i]<-sd(x)
  pval[i]<-2*(1-pt(abs(sqrt(n)*(mu_hat[i]-mu)/sigma_hat[i]),n-1))
}

## ----fig.width=10-------------------------------------------------------------
print(1-mean(pval<=0.05))

## ----fig.width=10-------------------------------------------------------------
mu <- 2
sigma <- 2 
m<-1e5
n<-20
set.seed(1234)
mu_hat<-sigma_hat<-pval<-numeric(m)
for(i in 1:m){
  x<-rchisq(n, df=2)
  mu_hat[i]<-mean(x)
  sigma_hat[i]<-sd(x)
  pval[i]<-2*(1-pt(abs(sqrt(n)*(mu_hat[i]-mu)/sigma_hat[i]),n-1))
}

## ----fig.width=10-------------------------------------------------------------
print(1-mean(pval<=0.05))

## ----fig.width=10-------------------------------------------------------------
alpha=0.05
m<-1e4
n<-20
set.seed(1234)
UCL<-numeric(m)
for(i in 1:m){
  x<-rchisq(n, df=2)
  UCL[i]<- (n-1)*var(x)/qchisq(alpha, df=n-1)
}

## ----fig.width=10-------------------------------------------------------------
print(mean(UCL<4))

## ----include=FALSE------------------------------------------------------------
n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv <- qnorm(.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
p.reject <- numeric(length(n)) #to store sim. results
m <- 1000 #num. repl. each sim.
for (i in 1:length(n)) {
sktests <- numeric(m) #test decisions
for (j in 1:m) {
x <- rbeta(n[i],2,2)
#test decision is 1 (reject) or 0
sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
}
p.reject[i] <- mean(sktests) #proportion rejected
}
p.reject

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
n <- 30
m <- 200
nu <- seq(1, 20, 1)
N <- length(nu)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
   e <- nu[j]
   sktests <- numeric(m)
   for (i in 1:m) { #for each replicate
      x <- rbeta(n,e,e)
      sktests[i] <- as.integer(abs(sk(x)) >= cv)
   }
   pwr[j] <- mean(sktests)
}
#plot power vs nu
plot(nu, pwr, type = "b",xlab = 'a', ylim = c(0,0.2),main='Power of the skewness test against the Beta(a,a) distribution')
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(nu, pwr+se, lty = 3)
lines(nu, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
n <- 30
m <- 200
nu <- seq(1, 20, 1)
N <- length(nu)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
   e <- nu[j]
   sktests <- numeric(m)
   for (i in 1:m) { #for each replicate
      x <- rt(n,e)
      sktests[i] <- as.integer(abs(sk(x)) >= cv)
   }
   pwr[j] <- mean(sktests)
}
#plot power vs nu
plot(nu, pwr, type = "b",xlab = 'v', ylim = c(0,1),main='Power of the skewness test against the t(v) distribution')
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(nu, pwr+se, lty = 3)
lines(nu, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
#m: number of replicates; n: sample size
varcount<-function(sigma1,sigma2,m,n){
  mean(replicate(m, expr={
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    count5test(x, y)
  }))
}

varftest<-function(sigma1,sigma2,m,n){
  mean(replicate(m, expr={
    x <- rnorm(n, 0, sigma1)
    y <- rnorm(n, 0, sigma2)
    a <- var.test(x, y)$p.value
    as.integer(a < 0.055)
  }))
}

## -----------------------------------------------------------------------------
print(c(varcount(1,1.5,200,20),varftest(1,1.5,200,20)))

## -----------------------------------------------------------------------------
print(c(varcount(1,1.5,200,50),varftest(1,1.5,200,50)))

## -----------------------------------------------------------------------------
print(c(varcount(1,1.5,200,100),varftest(1,1.5,200,100)))

## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=200
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=200
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap) #for the law data
theta_hat<-cor(law$LSAT, law$GPA)
n<-nrow(law)
theta_jack <- numeric(n)
for(i in 1:n){
  x<-law[-i,]
  theta_jack[i]<-cor(x$LSAT, x$GPA)
}

bias_jack <- (n-1)*(mean(theta_jack)-theta_hat)
se_jack <- sqrt((n-1)*mean((theta_jack-theta_hat)^2))
round(c(original=theta_hat,bias.jack=bias_jack,se.jack=se_jack),3)

## -----------------------------------------------------------------------------
library(boot) #for the law data
B=1000
set.seed(1234)
theta_boot <- numeric(B)
thetaboot <- function(data, index) {
  #function to compute the statistic
  mean(data[index]) 
}
air=aircondit[,1]
boot.result <- boot(air, statistic = thetaboot, R = B)

## -----------------------------------------------------------------------------
print(boot.result)

## -----------------------------------------------------------------------------
boot.ci(boot.out = boot.result, conf = 0.95, type = c("norm","basic","perc","bca"))

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
n<-nrow(scor)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
print(theta_hat)

## -----------------------------------------------------------------------------
#Jackknife
theta_jack <- numeric(n)
for (i in 1:n){
  x<-scor[-i,]
  lambda_hat <- eigen(cov(x))$values
  theta_jack[i] <- lambda_hat[1] / sum(lambda_hat)
}
biasj<-(n-1)*(mean(theta_jack)-theta_hat)
sej<-sqrt((n-1)*(n-1)/n)*sd(theta_jack)
print(c(biasj, sej))

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)
# for leave-two-out cross validation
for (j in 2:n){
  for (i in 1:(j-1)){
    y <- magnetic[c(-i,-j)]
    x <- chemical[c(-i,-j)]
    #Linear model
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i]
    yhat2 <- J1$coef[1] + J1$coef[2] * chemical[j]
    a <- magnetic[i] - yhat1
    b <- magnetic[j] - yhat2
    e1[(j-1)*(j-2)/2 + i] <- (a*a+b*b)/2
    #Quadratic model
    J2 <- lm(y~ x + I(x^2))
    yhat1 <- J2$coef[1] + J2$coef[2] * chemical[i] + J2$coef[3] * chemical[i]^2
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[j] + J2$coef[3] * chemical[j]^2
    a <- magnetic[i] - yhat1
    b <- magnetic[j] - yhat2
    e2[(j-1)*(j-2)/2 + i] <- (a*a+b*b)/2
    #Exponential model
    J3 <- lm(log(y)~x)
    logyhat1 <- J3$coef[1] + J3$coef[2] * chemical[i]
    logyhat2 <- J3$coef[1] + J3$coef[2] * chemical[j]
    yhat1 <- exp(logyhat1)
    yhat2 <- exp(logyhat2)
    a <- magnetic[i] - yhat1
    b <- magnetic[j] - yhat2
    e3[(j-1)*(j-2)/2 + i] <- (a*a+b*b)/2
    #Log-Log model
    J4 <- lm(log(y)~log(x))
    logyhat1 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
    logyhat2 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
    yhat1 <- exp(logyhat1)
    yhat2 <- exp(logyhat2)
    a <- magnetic[i] - yhat1
    b <- magnetic[j] - yhat2
    e4[(j-1)*(j-2)/2 + i] <- (a*a+b*b)/2
  }
}
c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
set.seed(123)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1000
tests <- replicate(m, expr = {
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x - mean(x) #centered by sample mean
  y <- y - mean(y)
  count5test(x, y)
} )
alphahat <- mean(tests)
print(alphahat)

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)

maxoutp<-function(x,ix,sizes){
  data<-x[ix]
  n1 <- sizes[1]
  n2 <- sizes[2]
  data1<-data[1:n1]
  data2<-data[(n1+1):(n1+n2)]
  X<-data1-mean(data1)
  Y<-data2-mean(data2)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}
set.seed(123)
m=1000
sizes=c(n1,n2)
testp<-numeric(m)
for(i in 1:m){
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  data<-c(x,y)
  boot.obj <- boot(data = data, statistic = maxoutp, sim="permutation", R=999, sizes=sizes)
  e <- boot.obj$t0
  E <- c(boot.obj$t, e)
  a=mean(E > e)
  testp[i]<-as.integer(a>0.95)
}
mean(testp)

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0,0,0)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
  mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
rbimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1), mu1, sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2), mu2, sd2)
  return(x)
}
for(i in 1:m){
  mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
  mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=10
n2=30
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
f<-function(x) 0.5*exp(-abs(x))
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
set.seed(2001)
N <- 3000
sigma <- c(.05, .3, 1.5, 10)
x0 <- 10
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
index=1:N
plot(index, rw1$x, type="l", main="sigma=0.05", ylab="X")
plot(index, rw2$x, type="l", main="sigma=0.3", ylab="X")
plot(index, rw3$x, type="l", main="sigma=1.5", ylab="X")
plot(index, rw4$x, type="l", main="sigma=10", ylab="X")
print(c(rw1$k, rw2$k, rw3$k, rw4$k))

## -----------------------------------------------------------------------------
print(c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N))

## -----------------------------------------------------------------------------
y<-rw3$x[501:N]
hist(y, breaks=50, main="Laplace distribution", xlab="X", freq=FALSE)
a<-seq(-6,6,0.01)
lines(a, f(a))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}
set.seed(2001)
sigma <- 0.3 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 10000 #length of chains
b <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- rw.Metropolis(sigma, x0[i], n)$x
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
#plot psi for the four chains
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l", main=paste('sigma=0.3, X[0]=', x0[i]), xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="length-1000", ylab="R", main="R-hat statistics")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
set.seed(2001)
sigma <- 0.5 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 8000 #length of chains
b <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- rw.Metropolis(sigma, x0[i], n)$x
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
#plot psi for the four chains
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l", main=paste('sigma=0.5, X[0]=', x0[i]), xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="length-1000", ylab="R", main="R-hat statistics")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
f1<-function(x,k){
  a<-sqrt((k-1)*x^{2}/(k-x^{2}))
  return(1-pt(a,k-1))
}
f2<-function(x,k){
  a<-sqrt(k*x^{2}/(k+1-x^{2}))
  return(1-pt(a,k))
}

f<-function(x) f1(x,4)-f2(x,4)
curve(f, 0.01, 1.99, bty="l", xlab="x", ylab="y",main='k=4')
abline(h=0)
f<-function(x) f1(x,25)-f2(x,25)
curve(f, 0.01, sqrt(25)-0.01, bty="l", xlab="x", ylab="y",main='k=25')
abline(h=0)
f<-function(x) f1(x,100)-f2(x,100)
curve(f, 0.01, sqrt(100)-0.01, bty="l", xlab="x", ylab="y",main='k=100')
abline(h=0)
f<-function(x) f1(x,500)-f2(x,500)
curve(f, 0.01, sqrt(500)-0.01, bty="l", xlab="x", ylab="y",main='k=500')
abline(h=0)

## -----------------------------------------------------------------------------
root<-function(k){
  f1<-function(x){
    a<-sqrt((k-1)*x^{2}/(k-x^{2}))
    return(1-pt(a,k-1))
  }
  f2<-function(x){
    a<-sqrt(k*x^{2}/(k+1-x^{2}))
    return(1-pt(a,k))
  }
  f<-function(x) f1(x)-f2(x)
  return(uniroot(f,lower = 0.001, upper = 1.999))
}
k<-c(4:25,100,500,1000)
a<-numeric(length(k))
for (i in 1:length(k)){
  a[i]=root(k[i])$root
}
result<-matrix(a, nrow=1,dimnames=list(c('intersections'),c(4:25,100,500,1000)))
print(result)

## -----------------------------------------------------------------------------
iteration=15
n_A=444
n_B=132
n_OO=361
n_AB=63
p<-numeric(iteration+1)
q<-numeric(iteration+1)
n_AA<-numeric(iteration)
n_BB<-numeric(iteration)
loglh<-numeric(iteration)
p[1]<-q[1]<-0.5

for (i in 1:iteration) {
  # E-step: update the value of missing data (n_AA,n_BB)
  n_AA[i]=n_A*p[i]/(2-p[i]-2*q[i])
  n_BB[i]=n_B*q[i]/(2-q[i]-2*p[i])
  # M-step: update p and q  
  a<-n_AA[i]+n_A+n_AB
  b<-2*n_OO+2*n_A+n_B+n_AB-n_BB[i]
  c<-n_BB[i]+n_B+n_AB
  d<-2*n_OO+2*n_B+n_A+n_AB-n_AA[i]
  pt<-p[i+1]<-(a*d-a*c)/(b*d-a*c)
  qt<-q[i+1]<-(b*c-a*c)/(b*d-a*c)
  r<-1-pt-qt
  loglh[i]<-n_AA[i]*log(pt/r)+n_BB[i]*log(qt/r)+2*n_OO*log(r)+n_A*log(pt*r)+
    n_B*log(qt*r)+n_AB*log(pt*qt)+(n_A+n_B+n_AB-n_AA[i]-n_BB[i])*log(2)
}

## -----------------------------------------------------------------------------
print(matrix(c(p[-1],q[-1],loglh),nrow=3, byrow=T, dimnames = list(c('p','q','llh'),1:15)))

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
for(x in formulas){
  print(lm(x, mtcars))
}

## -----------------------------------------------------------------------------
lapply(formulas, function(x) lm(x,mtcars))

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(trials, function(x) x$p.value)

## -----------------------------------------------------------------------------
datalist <- list(mtcars, faithful)
lapply(datalist, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE) return(simplify2array(out))
  unlist(out)
}
mylapply(datalist, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp) 

cppFunction('NumericVector rwMetropolis(double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0;
  NumericVector u=runif(N);
  for (int i=1;i<N;i++) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1])))){
      x[i] = y[0];
    } 
    else {
        x[i] = x[i-1];
    }
  }
  return(x);
}')

set.seed(2001)
N <- 5000
sigma <- c(.05, .3, 1.5, 10)
x0 <- 10
rw1 <- rwMetropolis(sigma[1], x0, N)
rw2 <- rwMetropolis(sigma[2], x0, N)
rw3 <- rwMetropolis(sigma[3], x0, N)
rw4 <- rwMetropolis(sigma[4], x0, N)
index=1:N
plot(index, rw1, type="l", main="sigma=0.05 (C)", ylab="X")
plot(index, rw2, type="l", main="sigma=0.3 (C)", ylab="X")
plot(index, rw3, type="l", main="sigma=1.5 (C)", ylab="X")
plot(index, rw4, type="l", main="sigma=10 (C)", ylab="X")

## -----------------------------------------------------------------------------
f<-function(x) 0.5*exp(-abs(x))
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
      }
  }
  return(x)
}
set.seed(2001)
N <- 5000
sigma <- c(.05, .3, 1.5, 10)
x0 <- 10
rw.1 <- rw.Metropolis(sigma[1], x0, N)
rw.2 <- rw.Metropolis(sigma[2], x0, N)
rw.3 <- rw.Metropolis(sigma[3], x0, N)
rw.4 <- rw.Metropolis(sigma[4], x0, N) #result of the R function
index=1:N
plot(index, rw.1, type="l", main="sigma=0.05 (R)", ylab="X")
plot(index, rw.2, type="l", main="sigma=0.3 (R)", ylab="X")
plot(index, rw.3, type="l", main="sigma=1.5 (R)", ylab="X")
plot(index, rw.4, type="l", main="sigma=10 (R)", ylab="X")

## -----------------------------------------------------------------------------
qqplot(rw1,rw.1,main='sigma=0.05',xlab='Rcpp',ylab='R')
qqplot(rw2,rw.2,main='sigma=0.3',xlab='Rcpp',ylab='R')
qqplot(rw3,rw.3,main='sigma=1.5',xlab='Rcpp',ylab='R')
qqplot(rw4,rw.4,main='sigma=10',xlab='Rcpp',ylab='R')

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(rwR=rw.Metropolis(sigma[1], x0, N),rwC=rwMetropolis(sigma[1], x0, N))
summary(ts)[,c(1,3,5,6)]


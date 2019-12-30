## ------------------------------------------------------------------------
n <- 10000
#par(mfrow = c(2, 2))
set.seed(75)
for (i in 1:4) {
sigma <- 10^(i-2)
u  <- runif(n)
x <- (-2*sigma^2*log(1-u))^(1/2)
h <- (-2*sigma^2*log(1/2))^(1/2)
h <- h/sigma^2*exp(-h^2/(2*sigma^2)) # to find the maximum value of f(x)
#par(mar=c(5.1, 4.1, 4, 2.1)) # to let the formula show complete
hist(x, prob = TRUE, main = expression(f(x) == frac(x,sigma^2)*e^-frac(x^2,2*sigma^2)), sub = paste("sigma","=",sigma), ylim = c(0, 1.1*h))
y <- seq(0, (-2*sigma^2*log(0.001))^(1/2), .01) # to let the x-range of hist and line be similar
lines(y, y/sigma^2*exp(-y^2/(2*sigma^2)))
}

## ------------------------------------------------------------------------
n <- 1000
set.seed(100)
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 3, 1)
u <- runif(n)
p1 <- 0.75
k <- as.integer(u < p1)
z1 <- k*x1 + (1 - k)*x2
hist(z1, prob = T, main = paste("mixture","with","p1","=",p1))
#par(mfrow = c(3,3))
for (i in 1:9) {
u <- runif(n)
p1 <- 0.1*i
k <- as.integer(u < p1)
z1 <- k*x1 + (1 - k)*x2
hist(z1, prob = T, main = paste("mixture","with","p1","=",p1))
}

## ------------------------------------------------------------------------
n = 4
Sigma <- matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2)
set.seed(1000)
rwish.bart <- # a function to generate samples from a Wishart distribution Wd(n,Sigma) based on Bartlett's decomposition
  function(n, Sigma){ #n>d+1>=1
    d = sqrt(length(Sigma))
    T <- matrix(rep(0,d^2), nrow = d, ncol = d)
    for (i in 1:d) {
      for (j in 1:d) {
        if (i>j){T[i,j] = rnorm(1)} 
        if (i==j){T[i,j] = rchisq(1, n-i+1)}
      }
    }
    A <- T %*% t(T)
    L <- t(chol(Sigma))
    X <- L %*% A %*% t(L)
    X
  }
rwish.bart(n, Sigma)

## ------------------------------------------------------------------------
m <- 1e4
set.seed(1000)
x <- runif(m, min = 0, max = pi/3)
theta.hat <- mean(sin(x))*pi/3
print(c(theta.hat,-cos(pi/3)+cos(0)))#Compare the estimte and the true value

## ------------------------------------------------------------------------
MC <- function( R , antithetic = TRUE) { # Write a function to generate the estimate with antithetic variables or not
  u <- runif(R/2)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)
  theta <- mean(exp(-u)/(1+u^2))
  theta
}
m <- 1000
MC1 <- MC2 <- numeric(m)
set.seed(100)
for (i in 1:m) {
  MC1[i] <- MC(R = 1000, anti = FALSE)# the estimate without variance reduction
  MC2[i] <- MC(R = 1000)# the estimate with variance reduction
}
c(sd(MC1),sd(MC2),sd(MC2)/sd(MC1))

## ------------------------------------------------------------------------
M <- 10000
k <- 5
r <- M/k # replicates per stratum
N <- 100 # number of times to repeat the estimation
T2 <- numeric(k)
estimate <- matrix(0, N, 2)
g <- function(x)exp(-x)/(1+x^2)
f <- function(x)exp(-x)/(1-exp(-1))
F <- function(x)-log(1-x*(1-exp(-1))) # the inverse function of the cdf of f(x)
set.seed(10)
for (i in 1:N) 
{
  u <- runif(M)
  x <- F(u)
  estimate[i,1] <- mean(g(x)/f(x))
  for (j in 1:k ) 
  {
    u <- runif( M/k, F((j-1)/k),F(j/k))
    x<- F(u)
    T2[j] <- mean(g(x)/(f(x)))
  }
  estimate[i,2] <- mean(T2)
}
knitr::kable(rbind(apply(estimate,2,mean),apply(estimate,2,sd)), formate = "html",col.names = c("5.10","5.13"))

## ------------------------------------------------------------------------
n <- 20
alpha <- .05
N <- 1000
TL1 <- TL2 <- UCL <-numeric(N)
set.seed(12345)
for (i in 1:N) {
  x <- rchisq(n, df = 2)
  y <- qt(1-alpha/2, df = n-1)*sqrt(var(x)/n)
  TL1[i] <- mean(x)-y# the lower confidence limit of t-interval
  TL2[i] <- mean(x)+y# the upper confidence limit of t-interval
  UCL[i] <- (n-1)*var(x)/qchisq(alpha, df = n-1)# the upper confidence limit of the variance
}
mean(TL1 < 2 & TL2 > 2)# the coverage probability of the t-interval
mean(UCL > 4)# the coverage probability of the one side confidence interval of the variance
rm(list = ls())

## ------------------------------------------------------------------------
n <- 100#the size of normal samples
m <- 100#the size of skewness samples
N <- 100#the number of simulations
qua <- c(0.025,0.05,0.95,0.975)#the vecter of quantiles
set.seed(1000)
sk <- function(x){# compute the sample skewness coeff
  xmean <- mean(x)
  m3 <- mean((x-xmean)^3)
  m2 <- mean((x-xmean)^2)
  m3/m2^1.5
}
q1 <- q2 <- q3 <- q4 <-numeric(N)# to storage the 0.025,0.05,0.95 and 0.975 quantiles of the skewness
skw <- numeric(m)# to storage the skewness of the samples
for (i in 1:N) {
  for (j in 1:m) {
    x <- rnorm(n)
    skw[j] <- sk(x)
  }
  q <- quantile(skw,qua,names = F)#find the quantiles
  q1[i] <- q[1]
  q2[i] <- q[2]
  q3[i] <- q[3]
  q4[i] <- q[4]
}
q.hat <- numeric(4)
q.hat[1] <- mean(q1)#0.025 quantile estimate
q.hat[2] <- mean(q2)#0.05 quantile estimate
q.hat[3] <- mean(q3)#0.95 quantile estimate
q.hat[4] <- mean(q4)#0.975 quantile estimate
q.hat#the estimates of such quantiles

qua.var <- function(q, n){#the 2.14 formula
  xq <- qnorm(q, sd = sqrt(6*(n-2)/((n+1)*(n+3))))
  q*(1-q)/(n*dnorm(xq, sd = sqrt(6*(n-2)/((n+1)*(n+3))))^2)
}
q.sd <- numeric(4)# to storage the standard error
for (i in 1:4) {
  q.sd[i] <- sqrt(qua.var(qua[i], n))
}
q.sd# the standard error of the estimates from 2.14

q.large <- numeric(4)
for (i in 1:4) {
  q.large[i] <- qnorm(qua[i], sd = sqrt(6/n))
}
q.large# the estimates of quantiles with the large sample approximation

## ---- echo=FALSE---------------------------------------------------------
n <- 1000#the size of normal samples
m <- 100#the size of skewness samples
N <- 100#the number of simulations
qua <- c(0.025,0.05,0.95,0.975)#the vecter of quantiles
set.seed(1000)
sk <- function(x){# compute the sample skewness coeff
  xmean <- mean(x)
  m3 <- mean((x-xmean)^3)
  m2 <- mean((x-xmean)^2)
  m3/m2^1.5
}
q1 <- q2 <- q3 <- q4 <-numeric(N)# to storage the 0.025,0.05,0.95 and 0.975 quantiles of the skewness
skw <- numeric(m)# to storage the skewness of the samples
for (i in 1:N) {
  for (j in 1:m) {
    x <- rnorm(n)
    skw[j] <- sk(x)
  }
  q <- quantile(skw,qua,names = F)#find the quantiles
  q1[i] <- q[1]
  q2[i] <- q[2]
  q3[i] <- q[3]
  q4[i] <- q[4]
}
q.hat <- numeric(4)
q.hat[1] <- mean(q1)#0.025 quantile estimate
q.hat[2] <- mean(q2)#0.05 quantile estimate
q.hat[3] <- mean(q3)#0.95 quantile estimate
q.hat[4] <- mean(q4)#0.975 quantile estimate
q.hat#the estimates of such quantiles

q.large <- numeric(4)
for (i in 1:4) {
  q.large[i] <- qnorm(qua[i], sd = sqrt(6/n))
}
q.large# the estimates of quantiles with the large sample approximation

## ------------------------------------------------------------------------
alpha <- 0.1 #the significance level
n <- 100 #the size of samples
m <- 2000 #the number of replicates
a <- seq(0, 10, 0.5) #the alphas of beta distributions
N <- length(a)
pwr.beta <- numeric(N)
set.seed(99)

sk <- function(x){# compute the sample skewness coeff
  xmean <- mean(x)
  m3 <- mean((x-xmean)^3)
  m2 <- mean((x-xmean)^2)
  m3/m2^1.5
}

cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2)/((n+1)*(n+3)))) #critical value for the skewness test

for (i in 1:N) {
  ai <- a[i]
  sktests <- numeric(m)
  for (j in 1:m) {
    x <- rbeta(n, ai, ai)
    sktests[j] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr.beta[i] <- mean(sktests)
}

plot(a, pwr.beta, type = "b", xlab = bquote(alpha), ylim = c(0, 1))
abline(h =0.1, lty = 3)
se <- sqrt(pwr.beta*(1-pwr.beta)/m)
lines(a, pwr.beta+se, lty = 3)
lines(a, pwr.beta-se, lty = 3)

v <- 1:20 #the vs of t distributions
N <- length(v)
pwr.t <- numeric(N)
for (i in 1:N) {
  vi <- v[i]
  sktests <- numeric(m)
  for (j in 1:m) {
    x <- rt(n, vi)
    sktests[j] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr.t[i] <- mean(sktests)
}

plot(v, pwr.t, type = "b", xlab = bquote(v), ylim = c(0, 1))
abline(h =0.1, lty = 3)
se <- sqrt(pwr.t*(1-pwr.t)/m)
lines(v, pwr.t+se, lty = 3)
lines(v, pwr.t-se, lty = 3)

## ------------------------------------------------------------------------
n <- 20 #the size of samples
alpha <- 0.05 #the significance level
mu0 <- 1 #the mean of the chisq(1)
m <- 10000 #the number of replicates
p <- numeric(m) #the storage for p-values
set.seed(100)
for (i in 1:m) {
  x <- rchisq(n, mu0)
  ttest <- t.test(x,mu = mu0)
  p[i] <- ttest$p.value
}

p.hat <- mean(p<alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
n <- 20 #the size of samples
alpha <- 0.05 #the significance level
mu0 <- 1 #the mean of the U(0, 2)
m <- 10000 #the number of replicates
p <- numeric(m) #the storage for p-values
set.seed(19)
for (i in 1:m) {
  x <- runif(n, 0, 2)
  ttest <- t.test(x, mu = mu0)
  p[i] <- ttest$p.value
}

p.hat <- mean(p<=alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
n <- 20 #the size of samples
alpha <- 0.05 #the significance level
mu0 <- 1 #the mean of the E(1)
m <- 10000 #the number of replicates
p <- numeric(m) #the storage for p-values
set.seed(100)
for (i in 1:m) {
  x <- rexp(n)
  ttest <- t.test(x, mu = mu0)
  p[i] <- ttest$p.value
}

p.hat <- mean(p<=alpha)
se.hat <- sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat, se.hat))

## ------------------------------------------------------------------------
library(bootstrap)
score <- scor #store the scor data into score
pairs(score[,1:5]) #the scatter plots for each pair of test scores
cormatrix <- cor(score) #the sample correlation matrix
round(cormatrix, 2)
r <- function(x, i){ #the correlation of columns 1 and 2
  cor(x[i, 1], x[i, 2])
}
library(boot)
set.seed(100)
boot.sd <- function(i, j){ #bootstrap estimates of the standard error
  obj <- boot(data = score[,c(i,j)], statistic = r, R = 2000)
  sd(obj$t)
}
rho12 <- boot.sd(1, 2)
rho34 <- boot.sd(3, 4)
rho35 <- boot.sd(3, 5)
rho45 <- boot.sd(4, 5)
round(c(rho12, rho34, rho35, rho45), 3)

## ------------------------------------------------------------------------
library(boot)
n <- 20 #sample size
N <- 1000
sk <- 0
set.seed(100)

ski <- function(x, i){# compute the sample skewness coeff
  xmean <- mean(x[i])
  m3 <- mean((x[i]-xmean)^3)
  m2 <- mean((x[i]-xmean)^2)
  m3/m2^1.5
}

ci.norm <- ci.basic <- ci.perc <- matrix(NA, N, 2)

for (i in 1:N) {
  sam <- rnorm(n)
  de <- boot(sam, statistic = ski, R = 1000)
  ci <- boot.ci(de, type = c("norm", "basic", "perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}
cat('norm: cover rate = ', mean(ci.norm[,1]<=sk & ci.norm[,2]>=sk), ', missing rate left = ', mean(ci.norm[,1]>sk), ', missing rate right = ', mean(ci.norm[,2]<sk)) #norm

cat('basic: cover rate = ', mean(ci.basic[,1]<=sk & ci.basic[,2]>=sk), ', missing rate left = ', mean(ci.basic[,1]>sk), ', missing rate right = ', mean(ci.basic[,2]<sk)) #basic

cat('perc: cover rate = ', mean(ci.perc[,1]<=sk & ci.perc[,2]>=sk), ', missing rate left = ', mean(ci.perc[,1]>sk), ', missing rate right = ', mean(ci.perc[,2]<sk)) #perc

## ------------------------------------------------------------------------
n <- 20 #sample size
N <- 1000
sk <- sqrt(8/5)
set.seed(100)

ski <- function(x, i){# compute the sample skewness coeff
  xmean <- mean(x[i])
  m3 <- mean((x[i]-xmean)^3)
  m2 <- mean((x[i]-xmean)^2)
  m3/m2^1.5
}

ci.norm <- ci.basic <- ci.perc <- matrix(NA, N, 2)

for (i in 1:N) {
  sam <- rchisq(n, 5)
  de <- boot(sam, statistic = ski, R = 1000)
  ci <- boot.ci(de, type = c("norm", "basic", "perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}
cat('norm: cover rate = ', mean(ci.norm[,1]<=sk & ci.norm[,2]>=sk), ', missing rate left = ', mean(ci.norm[,1]>sk), ', missing rate right = ', mean(ci.norm[,2]<sk)) #norm

cat('basic: cover rate = ', mean(ci.basic[,1]<=sk & ci.basic[,2]>=sk), ', missing rate left = ', mean(ci.basic[,1]>sk), ', missing rate right = ', mean(ci.basic[,2]<sk)) #basic

cat('perc: cover rate = ', mean(ci.perc[,1]<=sk & ci.perc[,2]>=sk), ', missing rate left = ', mean(ci.perc[,1]>sk), ', missing rate right = ', mean(ci.perc[,2]<sk)) #perc

## ------------------------------------------------------------------------
library(bootstrap)
score <- scor #store the scor data into score
pca1 <- function(x){
  m <- nrow(x)
  covmatrix <- (m-1)/m*cov(x) #MLE of Sigma
  eigen(covmatrix)$values[1]/sum(eigen(covmatrix)$values) #the proportion of first principle component
}
theta.hat <- pca1(score)
n <- nrow(score)
theta.jack <- numeric(n)
for (i in 1:n) {
  theta.jack[i] <- pca1(score[-i,])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat) #bias
se.jack <- sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2)) #standard error
round(c(original = theta.hat, bias.jack = bias.jack, se.jack = se.jack),3)

## ------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
  J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3

  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +
  J4$coef[3] * chemical[k]^2 +J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
summary(J1)
summary(J2)
summary(J3)
summary(J4)
detach(ironslag)

## ------------------------------------------------------------------------
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 100 #times of replicates
R <- 1000 #times of permutations
set.seed(10)
stat <- numeric(m)
K <- 1:50
reps <- numeric(R)

for (i in 1:m) {
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  z <- c(x, y)
  for (j in 1:R) {
    k <- sample(K, size = n1, replace = FALSE)
    x <- z[k];y <-z[-k]
    reps[j] <- maxout(x,y)
  }
  stat[i] <- quantile(reps,0.95)
}
mean(stat)

count7test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 7))
}
tests <- replicate(m, expr = {
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
count7test(x, y)
} )
alphahat <- mean(tests) #type I error
print(alphahat)

## ------------------------------------------------------------------------
library(Ball)
library(MASS)
library(boot)
dCov <- function(x, y) {
x <- as.matrix(x); y <- as.matrix(y)
n <- nrow(x); m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d); M <- mean(d)
a <- sweep(d, 1, m); b <- sweep(a, 2, m)
b + M
}
A<- Akl(x); B <- Akl(y)
sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q <- dims[2]
d <- p + q
x <- z[ , 1:p] #leave x as is
y <- z[ix, -(1:p)] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
set.seed(100)
n <- c(20, 50, 100, 200)
m <- 50
power.cor1 <- power.ball1 <- power.cor2 <- power.ball2 <- numeric(4)
p.cor1 <- p.ball1 <- p.cor2 <- p.ball2 <-numeric(m)
alpha <- 0.1
for (j in 1:4) {
for (i in 1:m) {
I <- matrix(c(1,0,0,1),2,2)
X <- mvrnorm(n[j],rep(0,2),I)
E <- mvrnorm(n[j],rep(0,2),I)
X <- X/4
Y1 <- X+E #model1
Y2 <- X * E #model2
Z1 <- cbind(X, Y1)
Z2 <- cbind(X, Y2)
boot.obj1 <- boot(data = Z1, statistic = ndCov2, R=100, sim = "permutation",dims = c(2,2))
boot.obj2 <- boot(data = Z2, statistic = ndCov2, R=100, sim = "permutation",dims = c(2,2))
tb1 <- c(boot.obj1$t0, boot.obj1$t)
p.cor1[i] <- mean(tb1>=tb1[1])
tb2 <- c(boot.obj2$t0, boot.obj2$t)
p.cor2[i] <- mean(tb2>=tb2[1])
p.ball1[i] <- bcov.test(X, Y1, R = 100, seed = i*123)$p.value
p.ball2[i] <- bcov.test(X, Y2, R = 100, seed = i*123)$p.value
}
power.cor1[j] <- mean(p.cor1<alpha)
power.ball1[j] <- mean(p.ball1<alpha)
power.cor2[j] <- mean(p.cor2<alpha)
power.ball2[j] <- mean(p.ball2<alpha)
}

plot(n,power.cor1,type = "b",ylab = "power1", xlim=c(0,200), ylim=c(0,1) , col = "1")
lines(n,power.ball1,type = "b", col="2")
legend("bottomright", pch = c(15, 15), legend = c("power.cor1","power.ball1"), col = c(1, 2), bty = "n")

plot(n,power.cor2,type = "b",ylab = "power2", xlim=c(0,200),ylim=c(0.5,1) , col = "1")
lines(n,power.ball2,type = "b", col="2")
legend("bottomright", pch = c(15, 15), legend = c("power.cor2","power.ball2"), col = c(1, 2), bty = "n")
rm(list = ls())

## ------------------------------------------------------------------------
library(GeneralizedHyperbolic)
rw.Metropolis <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if(u[i] <= (dskewlap(y)/dskewlap(x[i-1])))
      x[i] <- y
    else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x, k=k))
}

N <- 2000
sigma <- c(0.05, 0.5, 2, 16)

x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

no.reject <- data.frame(sigma=sigma,no.reject=c(rw1$k, rw2$k, rw3$k, rw4$k),accept.rate=c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N))
    knitr::kable(no.reject,format='html')
    
#par(mfrow=c(2,2))  #display 4 graphs together
refline <- qskewlap(c(.025, .975))
rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
    plot(rw[,j], type="l",
          xlab=bquote(sigma == .(round(sigma[j],3))),
          ylab="X", ylim=range(rw[,j]))
    abline(h=refline)
}
#par(mfrow=c(1,1)) #reset to default

## ------------------------------------------------------------------------
x <- 10
y <- log(exp(x))
z <- exp(log(x))
x==y
x==z
y==z
all.equal(x,y)
all.equal(x,z)
all.equal(y,z)

## ------------------------------------------------------------------------
k <- c(4:25,100,500,1000) #values of ks
s <- numeric(25) # to store the values of roots
for (i in 1:25) {
g <- function(a)
{
  q1 <- sqrt(a^2*(k[i]-1)/(k[i]-(a^2)))
  p1 <- pt(q1, k[i]-1,lower.tail = F, log.p = TRUE)
  q2 <- sqrt(a^2*k[i]/(k[i]+1-(a^2)))
  p2 <- pt(q2, k[i],lower.tail = F, log.p = TRUE)
  return(p1-p2)
}
s[i] <- uniroot(g,c(-1,sqrt(k[i])-0.01))$root
}
r <- cbind(k,s)
r

## ------------------------------------------------------------------------
k <- c(4:25,100,500,1000)
s <- numeric(25)
for (i in 1:25) {
  f <- function(u, k) (1+u^2/(k-1))^(-k/2)
  ck <- function(a, k) sqrt(a^2*(k-1)/(k-a^2))
  g <- function(a){
    in1 <- integrate(f, lower = ck(a, k[i]), upper = Inf, k = k[i])$value
    in2 <- integrate(f, lower = ck(a, k[i]+1), upper = Inf, k = k[i]+1)$value
    p1 <- log(2) + lgamma(k[i]/2) - 0.5*log(pi*(k[i]-1)) - lgamma((k[i]-1)/2) + log(in1)
    p2 <- log(2) + lgamma((k[i]+1)/2) - 0.5*log(pi*k[i]) - lgamma(k[i]/2) + log(in2)
    return(p1-p2)
  }
  s[i] <- uniroot(g,c(-1, 1.9))$root
}
r <- cbind(k,s)
r

## ------------------------------------------------------------------------
library(nloptr)
# Mle 
eval_f0 = function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 = function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  return(sum(x)-0.999999)
}

opts = list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle = NULL
r = matrix(0,1,2)
r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j = 2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res = nloptr( x0=c(0.3,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=28,n.B=24,nOO=41,nAB=70 )
j = j+1
r = rbind(r,res$solution)
mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
#the result of EM algorithm
r 
#the max likelihood values
plot(mle,type = 'l')


## ------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

la1 <- lapply(formulas, lm, data = mtcars) #lapply

lo1 <- vector("list", length(formulas)) #for loop
for (i in seq_along(formulas)) {
  lo1[[i]] <- lm(formulas[[i]], data = mtcars)
}

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

la2 <- lapply(bootstraps, lm, formula = mpg ~ disp) #apply

lo2 <- vector("list", length(bootstraps)) #for loop
for (i in seq_along(bootstraps)){
  lo2[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}

## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

sapply(la1, rsq) #3 la1
sapply(lo1, rsq) #3 lo1

sapply(la2, rsq) #4 la2
sapply(lo2, rsq) #4 lo2

## ------------------------------------------------------------------------
trials <- replicate(
  100, 
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

sapply(trials, function(x) x[["p.value"]])
rm(list = ls())

## ------------------------------------------------------------------------
library(parallel)
parsapply <- function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, 
  chunk.size = NULL) 
{
  FUN <- match.fun(FUN)
  answer <- parLapply(cl = cl, X = as.list(X), fun = FUN, 
    ..., chunk.size = chunk.size)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}
rm(list = ls())

## ------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

rwM<-function(x0,sigma,N){ # R function
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,mean=x[i-1],sd=sigma)
    if(u[i]<=exp(-abs(y))/exp(-abs(x[i-1]))) {
      x[i]<-y
      k<-k+1}
    else x[i]<-x[i-1]
  }
  return(list(x=x,k=k))
}

#Rcpp function
cppFunction('List rwMc(double x0, double sigma, int N){
NumericVector x(N);
x[0]=x0;
int k=0;
for(int i=1;i<N;i++){

double y=as<double>(rnorm(1,x[i-1], sigma));
double u=as<double>(runif(1));

if (u<=exp(-abs(y))/exp(-abs(x[i-1]))) {
x[i]=y;
k=k+1;  
}
else x[i]=x[i-1];
}

return(List::create(Named("x")=x,Named("k")=k));
}')

sigma<-c(0.05,0.5,2,16)
x0<-25
N<-2000
set.seed(100)

rw1<-rwM(x0,sigma[1],N) #R samples
rw2<-rwM(x0,sigma[2],N)
rw3<-rwM(x0,sigma[3],N)
rw4<-rwM(x0,sigma[4],N)

rwc1<-rwMc(x0,sigma[1],N) #Rcpp samples
rwc2<-rwMc(x0,sigma[2],N)
rwc3<-rwMc(x0,sigma[3],N)
rwc4<-rwMc(x0,sigma[4],N)
#par(mfrow=c(2,2))

qqplot(rw1$x,rwc1$x) #compare qqplot
qqplot(rw2$x,rwc2$x)
qqplot(rw3$x,rwc4$x)
qqplot(rw4$x,rwc4$x)

#compare computing time
ts1<-microbenchmark(rwM(x0,sigma[1],N),rwMc(x0,sigma[1],N)) 
ts2<-microbenchmark(rwM(x0,sigma[2],N),rwMc(x0,sigma[2],N))
ts3<-microbenchmark(rwM(x0,sigma[3],N),rwMc(x0,sigma[3],N))
ts4<-microbenchmark(rwM(x0,sigma[4],N),rwMc(x0,sigma[4],N))
summary(ts1)[,c(1,3,5,6)]
summary(ts2)[,c(1,3,5,6)]
summary(ts3)[,c(1,3,5,6)]
summary(ts4)[,c(1,3,5,6)]


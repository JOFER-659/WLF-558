###############################################################################
# POLR in Bayesian Framework
###########
# BY: Jacob Hofer
###############################################################################

# See POLR.R for details on the POLR model and exploration of the 
# likelihood function

# First lets library our necessary packages
library(R2jags)
library(FSAdata)
library(coda)
library(MASS)

#Loading in dataset
data("Jonubi1")
jonubi<-Jonubi1[c(-410,-409,-408),c(2,1)]
rm(Jonubi1)
names(jonubi)<-c('Age', 'Length')

##########################################################################
# Population Wide Model
##########################################################################

# Prior values for alphas
# Find a good guess for beta using a frequentist POLR
# Use average length as x
jonubi$Age <- factor(jonubi$Age)
summary(freqpolr<-polr(Age ~ Length, data = jonubi))
jonubi$Age <- as.integer(jonubi$Age)
p1<-(sum(jonubi[,1] == 1)/nrow(jonubi))
p2<-(sum(jonubi[,1] == 2)/nrow(jonubi))
p3<-(sum(jonubi[,1] == 3)/nrow(jonubi))
p4<-(sum(jonubi[,1] == 4)/nrow(jonubi))
p5<-(sum(jonubi[,1] == 5)/nrow(jonubi))
Bx <- 1.986 * mean(jonubi$Length)
#Alpha 1
a1mean<- Bx-log((1/p1)-1)
#Alpha 2
a2mean<-Bx-log((1/(p1+p2))-1)
#Alpha 3
a3mean<- Bx-log((1/(p1+p2+p3))-1)
#Alpha 4
a4mean<- Bx-log((1/(p1+p2+p3+p4))-1)

# Now lets define our likelihood function and priors to be passed to JAGS
POLR.model <- function(){
  #likelihood
  for (i in 1:N){
    eta[i] <- beta*x[i]
    y[i] ~ dcat(p[i, 1:k])
    
    p[i,1] <- ilogit(alpha[1] - eta[i])
    for(k in 2:(k-1)){
      p[i, k] <- ilogit(alpha[k] - eta[i]) - ilogit(alpha[k-1] - eta[i])
    }
    p[i,k] <- 1 - ilogit(alpha[k-1] - eta[i])
  }
  
  #priors
  alpha[1] ~ dnorm(a1mean, 1) # informative priors
  alpha[2] ~ dnorm(a2mean, 1)
  alpha[3] ~ dnorm(a3mean, 1)
  alpha[4] ~ dnorm(a4mean, 1)
  beta ~ dunif(0,10) # uninformative prior, beta > 0, though
}

# House Keeping
# passing data
POLR.data <- list(y = jonubi$Age, x = jonubi$Length, N = length(jonubi$Length),
                  a1mean = a1mean, a2mean = a2mean, a3mean = a3mean,
                  a4mean = a4mean, k = length(unique(jonubi$Age)))
# parameters to keep track off
params <- c('alpha[1]', 'alpha[2]', 'alpha[3]', 'alpha[4]','beta')
# starting values of parameters
Inits<-function(){list('alpha[1]' = runif(1,5,25),'alpha[2]' = runif(1,7,27),
                       'alpha[3]' = runif(1,8,28), 'alpha[4]' = runif(1,10,30),
                       'beta' = runif(1,0,10))}
nc <- 3 #Number of chains per parameters
ni <- 50000 #Number of iterations of each MCMC chain
nb <- 5000 # burn-in
nt <- 5 # Thinning rate
post.samples <- (ni-nb)/(nt*nc) #Number of MCMC samples used to estimate the posterior of each parameter.

# Calling JAGS
jags.out.pop <- jags(model.file = POLR.model, data = POLR.data, inits = Inits,
                     parameters.to.save = params, n.chains = nc,
                     n.iter = ni, n.burnin = nb, n.thin = nt)

# Check R hat values for convergence
jags.out.pop

# Check trace plots for convergence
plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,1],type='l',xlab='MCMC Iteration',
     ylab='alpha1',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,1],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,1],col=4)

plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,2],type='l',xlab='MCMC Iteration',
     ylab='alpha2',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,2],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,2],col=4)

plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,3],type='l',xlab='MCMC Iteration',
     ylab='alpha3',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,3],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,3],col=4)

plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,4],type='l',xlab='MCMC Iteration',
     ylab='alpha4',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,4],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,4],col=4)

plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,5],type='l',xlab='MCMC Iteration',
     ylab='beta',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,5],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,5],col=4)

plot(1:dim(jags.out.pop$BUGSoutput$sims.array)[1],
     jags.out.pop$BUGSoutput$sims.array[,1,6],type='l',xlab='MCMC Iteration',
     ylab='deviance',col=2)
lines(jags.out.pop$BUGSoutput$sims.array[,2,6],col=3)
lines(jags.out.pop$BUGSoutput$sims.array[,3,6],col=4)

# Credible intervals for parameters of interest

# All the alpha's
HPDinterval(mcmc(jags.out.pop$BUGSoutput$sims.list$alpha),0.95)

# Beta 
HPDinterval(mcmc(jags.out.pop$BUGSoutput$sims.list$beta),0.95)


# Plot the posteriors of the Parameters
plot(density(jags.out.pop$BUGSoutput$sims.list$alpha[,1]),xlab='Alpha 1',main='')
plot(density(jags.out.pop$BUGSoutput$sims.list$alpha[,2]),xlab='Alpha 2',main='')
plot(density(jags.out.pop$BUGSoutput$sims.list$alpha[,3]),xlab='Alpha 3',main='')
plot(density(jags.out.pop$BUGSoutput$sims.list$alpha[,4]),xlab='Alpha 4',main='')
plot(density(jags.out.pop$BUGSoutput$sims.list$beta),xlab='Beta',main='')

##########################################################################
# Individual Pred Model ( NOT WORKING D: )
##########################################################################

#Alpha 1
a1mean<- 20
#Alpha 2
a2mean<-21
#Alpha 3
a3mean<- 22
#Alpha 4
a4mean<- 23

# Now lets define our likelihood function and priors to be passed to JAGS
POLR.model <- function(){
  #likelihood
  for (i in 1:N){
    eta[i] <- beta*x[i]
    y[i] ~ dcat(p[i, 1:k])
    
    p[i,1] <- ilogit(alpha[1] - eta[i])
    for(k in 2:(k-1)){
      p[i, k] <- ilogit(alpha[k] - eta[i]) - ilogit(alpha[k-1] - eta[i])
    }
    p[i,k] <- 1 - ilogit(alpha[k-1] - eta[i])
  }
  
  #priors
  alpha[1] ~ dnorm(a1mean, 1) # informative priors
  alpha[2] ~ dnorm(a2mean, 1)
  alpha[3] ~ dnorm(a3mean, 1)
  alpha[4] ~ dnorm(a4mean, 1)
  beta ~ dunif(0,10) # uninformative prior, beta > 0, though
}

# House Keeping
# passing data
POLR.data <- list(y = jonubi$Age, x = jonubi$Length, N = length(jonubi$Length),
                  a1mean = a1mean, a2mean = a2mean, a3mean = a3mean,
                  a4mean = a4mean, k = length(unique(jonubi$Age)))
# parameters to keep track off
params <- c('alpha[1]', 'alpha[2]', 'alpha[3]', 'alpha[4]','beta')
# starting values of parameters
Inits<-function(){list('alpha[1]' = runif(1,5,25),'alpha[2]' = runif(1,7,27),
                       'alpha[3]' = runif(1,8,28), 'alpha[4]' = runif(1,10,30),
                       'beta' = runif(1,0,10))}
nc <- 3 #Number of chains per parameters
ni <- 50000 #Number of iterations of each MCMC chain
nb <- 5000 # burn-in
nt <- 5 # Thinning rate
post.samples <- (ni-nb)/(nt*nc) #Number of MCMC samples used to estimate the posterior of each parameter.

# Calling JAGS
jags.out.ind <- jags(model.file = POLR.model, data = POLR.data, inits = Inits,
                     parameters.to.save = params, n.chains = nc,
                     n.iter = ni, n.burnin = nb, n.thin = nt)

# Check R hat values for convergence
jags.out.ind

# Check trace plots for convergence
plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,1],type='l',xlab='MCMC Iteration',
     ylab='alpha1',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,1],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,1],col=4)

plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,2],type='l',xlab='MCMC Iteration',
     ylab='alpha2',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,2],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,2],col=4)

plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,3],type='l',xlab='MCMC Iteration',
     ylab='alpha3',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,3],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,3],col=4)

plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,4],type='l',xlab='MCMC Iteration',
     ylab='alpha4',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,4],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,4],col=4)

plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,5],type='l',xlab='MCMC Iteration',
     ylab='beta',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,5],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,5],col=4)

plot(1:dim(jags.out.ind$BUGSoutput$sims.array)[1],
     jags.out.ind$BUGSoutput$sims.array[,1,6],type='l',xlab='MCMC Iteration',
     ylab='deviance',col=2)
lines(jags.out.ind$BUGSoutput$sims.array[,2,6],col=3)
lines(jags.out.ind$BUGSoutput$sims.array[,3,6],col=4)

# Credible intervals for parameters of interest

# All the alpha's
HPDinterval(mcmc(jags.out.ind$BUGSoutput$sims.list$alpha),0.95)

# Beta 
HPDinterval(mcmc(jags.out.ind$BUGSoutput$sims.list$beta),0.95)


# Plot the posteriors of the Parameters
plot(density(jags.out.ind$BUGSoutput$sims.list$alpha[,1]),xlab='Alpha 1',main='')
plot(density(jags.out.ind$BUGSoutput$sims.list$alpha[,2]),xlab='Alpha 2',main='')
plot(density(jags.out.ind$BUGSoutput$sims.list$alpha[,3]),xlab='Alpha 3',main='')
plot(density(jags.out.ind$BUGSoutput$sims.list$alpha[,4]),xlab='Alpha 4',main='')
plot(density(jags.out.ind$BUGSoutput$sims.list$beta),xlab='Beta',main='')

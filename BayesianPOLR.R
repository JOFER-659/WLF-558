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

#Loading in dataset
data("Jonubi1")
jonubi<-Jonubi1[c(-410,-409,-408),c(2,1)]
rm(Jonubi1)
names(jonubi)<-c('Age', 'Length')

# Finding values to use as priors
# Center of two age groups average lengths
#Alpha 1
a1mean<-((mean(jonubi[jonubi$Age == 2,2])-mean(jonubi[jonubi$Age == 1,2]))/2) + mean(jonubi[jonubi$Age == 1,2])
#Alpha 2
a2mean<-((mean(jonubi[jonubi$Age == 3,2])-mean(jonubi[jonubi$Age == 2,2]))/2) + mean(jonubi[jonubi$Age == 2,2])
#Alpha 3
a3mean<-((mean(jonubi[jonubi$Age == 4,2])-mean(jonubi[jonubi$Age == 3,2]))/2) + mean(jonubi[jonubi$Age == 3,2])
#Alpha 4
a4mean<-((mean(jonubi[jonubi$Age == 5,2])-mean(jonubi[jonubi$Age == 4,2]))/2) + mean(jonubi[jonubi$Age == 4,2])

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
  alpha[1] ~ dnorm(a1mean, 0.01)
  alpha[2] ~ dnorm(a2mean, 0.01)
  alpha[3] ~ dnorm(a3mean, 0.01)
  alpha[4] ~ dnorm(a4mean, 0.01)
  beta ~ dnorm(0,0.001)
}

# House Keeping
# passing data
POLR.data <- list(y = jonubi$Age, x = jonubi$Length, N = length(jonubi$Length),
                  a1mean = a1mean, a2mean = a2mean, a3mean = a3mean,
                  a4mean = a4mean, k = length(unique(jonubi$Age)))
# parameters to keep track off
params <- c('alpha[1]', 'alpha[2]', 'alpha[3]', 'alpha[4]','beta')
# starting values of parameters
Inits<-function(){list('alpha[1]' = runif(1,5,15),'alpha[2]' = runif(1,7,17),
                       'alpha[3]' = runif(1,8,18), 'alpha[4]' = runif(1,10,20),
                       'beta' = runif(1,0,5))}
nc <- 3 #Number of chains per parameters
ni <- 50000 #Number of iterations of each MCMC chain
nb <- 5000 # burn-in
nt <- 5 # Thinning rate
post.samples <- (ni-nb)/(nt*nc) #Number of MCMC samples used to estimate the posterior of each parameter.

# Calling JAGS
jags.out <- jags(model.file = POLR.model, data = POLR.data, inits = Inits,
                 parameters.to.save = params, n.chains = nc,
                 n.iter = ni, n.burnin = nb, n.thin = nt)

# Check R hat values for convergence
jags.out

# Check trace plots for convergence
plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,1],type='l',xlab='MCMC Iteration',
     ylab='alpha1',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,1],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,1],col=4)

plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,2],type='l',xlab='MCMC Iteration',
     ylab='alpha2',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,2],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,2],col=4)

plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,3],type='l',xlab='MCMC Iteration',
     ylab='alpha3',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,3],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,3],col=4)

plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,4],type='l',xlab='MCMC Iteration',
     ylab='alpha4',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,4],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,4],col=4)

plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,5],type='l',xlab='MCMC Iteration',
     ylab='beta',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,5],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,5],col=4)

plot(1:dim(jags.out$BUGSoutput$sims.array)[1],
     jags.out$BUGSoutput$sims.array[,1,6],type='l',xlab='MCMC Iteration',
     ylab='deviance',col=2)
lines(jags.out$BUGSoutput$sims.array[,2,6],col=3)
lines(jags.out$BUGSoutput$sims.array[,3,6],col=4)

# Credible intervals for parameters of interest

# All the alpha's
HPDinterval(mcmc(jags.out$BUGSoutput$sims.list$alpha),0.95)

# Beta 
HPDinterval(mcmc(jags.out$BUGSoutput$sims.list$beta),0.95)


# Plot the posteriors of the Parameters
plot(density(jags.out$BUGSoutput$sims.list$alpha[,1]),xlab='Alpha 1',main='')
plot(density(jags.out$BUGSoutput$sims.list$alpha[,2]),xlab='Alpha 2',main='')
plot(density(jags.out$BUGSoutput$sims.list$alpha[,3]),xlab='Alpha 3',main='')
plot(density(jags.out$BUGSoutput$sims.list$alpha[,4]),xlab='Alpha 4',main='')
plot(density(jags.out$BUGSoutput$sims.list$beta),xlab='Beta',main='')


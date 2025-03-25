#############################################################
## Continuation Ratio Logit Models
#########
## BY: Jacob Hofer
#############################################################
# Continuation Ratio Logit Models (CRLM's) are another method of 
# modeling ordinal data, similar to POLR. The main difference is
# that instead of using cumulative logits, like
# logit[P(Y <= j|x)], we use contuation ratio logits where
# logit[P(Y = j|Y >= j, x)]

# Again we have c categories. We only generate c-1 equations
# as the last equation is automatically defined by 1-P(Y>=c-1).

#############################################################
#load necessary packages
library(FSAdata)

#############################################################

#Writing the likelihood function

#function to calculate survival probabilities up to level j
myprod<-function(y_i, x_i, alpha, beta){
    prod(1 - (plogis(alpha[1:(y_i - 1)] + beta[1:(y_i - 1)] * x_i)))
}


crlm.log.lik <- function(data, alpha, beta){
  if (length(alpha) != length(beta)) {
    stop("alpha and beta must be the same length.")
  }
  value <- numeric(length = nrow(data))
  for (i in 1:nrow(data)){
    y_i <- data[i,1]
    x_i <- data[i,2]
    
    if (y_i == 1){ # probability if first category
      p_i <- plogis(alpha[1] + beta[1] * x_i)
    } else if (y_i == length(alpha) + 1){ #probability of last cat
      p_i <- myprod(y_i = y_i, x_i = x_i, alpha = alpha, beta = beta)
    } else {
      p_i = plogis(alpha[y_i] + beta[y_i] * x_i)*myprod(y_i = y_i,
                                                        x_i = x_i,
                                                        alpha = alpha,
                                                        beta = beta) #Prob of other cats
    }
    value[i] <- log(pmax(p_i, 1e-40))
  }
  return(sum(value))
}

data <- data.frame(Age = c(1,2,3,4,5), Length = c(2,4,5,7,9))
#Specify alphas (note: must be increasing)
alpha <- c(1,1,1,1)
#Specify Beta
beta<-c(1,2,2.3,3)

#Testing the function
crlm.log.lik(data,alpha,beta)


#### Now lets load the jonubi dataset and estimate the parameters of interest
# using optim

# Alphas go in first half of parms vector, betas in second
# Recall: length(alpha) = length(beta)
crlm.neg.log.lik <- function(parms, data){
  alpha <- parms[1:((length(parms))/2)] # Extract alphas
  beta <- parms[(((length(parms))/2)+1):length(parms)]
  value <- numeric(length = nrow(data))
  for (i in 1:nrow(data)){
    y_i <- data[i,1]
    x_i <- data[i,2]
    
    if (y_i == 1){ # probability if first category
      p_i <- plogis(alpha[1] + beta[1] * x_i)
    } else if (y_i == length(alpha) + 1){ #probability of last cat
      p_i <- myprod(y_i = y_i, x_i = x_i, alpha = alpha, beta = beta)
    } else {
      p_i = plogis(alpha[y_i] + beta[y_i] * x_i)*myprod(y_i = y_i,
                                                        x_i = x_i,
                                                        alpha = alpha,
                                                        beta = beta) #Prob of other cats
    }
    value[i] <- log(pmax(p_i, 1e-40))
  }
  return(-sum(value))
}

data("Jonubi1")
jonubi<-Jonubi1[c(-410,-409,-408),c(2,1)]
rm(Jonubi1)
names(jonubi)<-c('Age', 'Length')
jonubi$Age <- as.integer(jonubi$Age)

parms<-c(1,1,1,1,1,1,1,1)

optimresult <- optim(par = parms,
                     fn = crlm.neg.log.lik,
                     data = jonubi,
                     method = 'BFGS')

optimresult

### Work checking
library(numDeriv)
grad(func = crlm.neg.log.lik, x = parms, data = jonubi)


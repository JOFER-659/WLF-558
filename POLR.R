#########
# POLR (Proportional Odds Logistic Regression)
#########
# BY: Jacob Hofer
#########

###############################################################################
# POLR is a method of modeling ordinal data within a regression frame work. It
# is similar to multinomial regression but takes advantage of extra information
# stored within the data, given by its ordinal form. 

# Data typically comes in the form of 1,2,3,...,c where n is the number of
# caterogies. From the ordinal structure, it is assumed 1 < 2 < 3 < ... < c. 
# For example the likert scale data, we know disagree is less than neutral, and
# neutral is less than agree. 

# The model is generalized linear model, making use of the logit kernel.
# Cumulative probabilities of each category are found for an observation. 

# The model is typically written as 
# logit[P(Y <= j|x)] = a_j - B'x, where j = 1,2,3,...,c-1, a_j is an intercept,
# B is a vector of parameters, and x is a vector of predictor values.

# Notice that j only goes to c-1, this is because the last categories equation
# is already determined by the prior c-1 equations due to the constraint that
# p(Y <= j) = 1

# Lets think about how we will get probabilities from this model. The model in 
# its current form will give us log odds of being in catergory j or lower. So to
# turn this into probabilities, we must use the inverse logit function. 
# That is: p(Y <= j) = 1/(1 + exp(-(a_j - B'x))

# If we want to find the probability of a specific catergory, as opposed to 
# the cumulative distribution, we use basic probability definitions
# p(Y = 1) = p(Y <= 1), p(Y = 2) = p(Y<=2)-P(Y<=1), and
# p(Y = j) = p(Y <= j) - p(Y <= j-1) in general,
# p(Y = j) = 1 - p(Y <= j-1) for the last category

# The likelihood function for this model is not simple, but can be written out.
# The way to define the likelihood is to multiply the probability of the data
# being in its respective category. Lets write it out now (note we are only
# using one predictor variable)

# data should be in the form of a dataframe where
# the first column is response, second is predictor.
# alpha should be c-1 in length, and beta a scalar

###############################################################################

# Lets write out a function to calculate the log-likelihood if we know the data
# and the parameters
log.lik <- function(data, alpha, beta){
  value <- numeric(length = nrow(data))
  for (i in 1:nrow(data)){
  y_i <- data[i,1]
  x_i <- data[i,2]
  
  if (y_i == 1){ # probability if first category
   p_i <- plogis(alpha[1] - beta * x_i)
  } else if (y_i == length(alpha) + 1){ # Probability if second catergory
   p_i <- 1 - plogis(alpha[length(alpha)] - beta * x_i)
  } else { #Probability if any other category
   p_i <- plogis(alpha[y_i] - beta * x_i) - plogis(alpha[y_i-1] - beta * x_i)
  }
  value[i] <- log(p_i)
  }
  return(sum(value))
}

#Create example dataframe
data <- data.frame(Age = c(4,5,3,4,2,1,2), Length = c(10,12,7,9,5,4,6))
#Specify alphas (note: must be increasing)
alpha <- c(1,3,4,9)
#Specify Beta
beta<-1

#Testing the function
log.lik(data,alpha,beta)

###############################################################################

# Now lets write the log likelihood in a form where we can pass it to the 
# optim function, essentially finding the MLE's of the alphas and beta
# parameters.

#Note, parms should be of the form (alpha_1, alpha_2,..., alpha_n, beta)
neg.log.lik <- function(parms, data){
  alpha <- parms[1:(length(parms)-1)] # Extract alphas
  beta <- parms[length(parms)] #Extract Beta
  value <- numeric(length = nrow(data))
  for (i in 1:nrow(data)){
    y_i <- data[i,1]
    x_i <- data[i,2]
    
    if (y_i == 1){ # probability if first category
      p_i <- plogis(alpha[1] - beta * x_i)
    } else if (y_i == length(alpha) + 1){ # Probability if second catergory
      p_i <- 1 - plogis(alpha[length(alpha)] - beta * x_i)
    } else { #Probability if any other category
      p_i <- plogis(alpha[y_i] - beta * x_i) - plogis(alpha[y_i-1] - beta * x_i)
    }
    value[i] <- log(p_i)
  }
  return(-sum(value)) # Made negative because optim minimizes functions
}

#lets load a dataset to try this out
dat<-read.csv("SL_csv.csv")
dat <-dat[,c(3,2)]
names(dat) <- c('learn', 'age')
dat$learn <- factor(dat$learn, ordered = T)
dat$age <- as.numeric(dat$age)
str(dat)
head(dat)

# fit polr using MASS package to compare results
polrresult <- polr(learn ~ age, data = dat)

#turn learn into interger so likelihood function can handle it
dat$learn <- as.integer(dat$learn)

# Set starting values (4 alphas and 1 beta)
parms <- c(-4,-2,-1,1,0)

optimresult <- optim(par = parms,
      fn = neg.log.lik,
      data = dat,
      method = 'BFGS')

polrresult
optimresult

# Notice the results are near identical :)
# WOOHOO
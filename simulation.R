---
title: "ZINB simulation"
author: "Jian Huang"
date: "2023-08-25"
---
        
#load packages

 library(car) 
 library(MASS)  
 library(pscl)  
library(runjags)
library(rjags)

#Set your work path
setwd("your work path") 

set.seed(123)  # Random seed
n <- 1000  # Sample size


age <- runif(n,min =4, max = 12)
ethnic <- factor(sample(1:2,n, replace = TRUE), labels = c("W", "K"))
BMIAZ = rnorm(n, mean = 0, sd = 0.5)

# Negative binomial distribution parameter
theta <- 2
age_mean = mean(age)
age_sd = sd(age)
z_age = age-age_mean

B_beta = z_age*(-0.5) + as.numeric(ethnic)*1 + 1.5*BMIAZ

# Negative binomial mean
mu <- exp( -1 + B_beta) 
# zero-inflated prob
logit_p = exp( -1  - 0.5 * B_beta)
zero_prob <- logit_p/(1+logit_p)

zero_prob <- 0.25
# y count
count <- rnegbin(n, mu = mu, theta = theta)

p_x = sum(count == 0)/n
prob = (p_x*zero_prob)/((1-p_x)*(1-zero_prob))
prob[prob>1] =1
count_pro_zero = rbinom(sum(count != 0), 1, 1-prob)
count[count != 0] = count[count != 0] * count_pro_zero

age = z_age
data <- data.frame(count,age, ethnic, BMIAZ)



# hurdle function was used to fit the zero expansion negative binomial regression model
fit <- hurdle(count ~ age + ethnic + BMIAZ, dist = "negbin", link = "logit")


predicted_counts <- as.data.frame(predict(fit, newdata = data, type = "response"))
predicted_counts$counts = data$count
accc = sum(abs(predicted_counts[,1] - predicted_counts$counts))/n
 summary(fit)

# hurdle_zero_prob <- 1/(1 + exp(-fit$coefficients$zero[1] - fit$coefficients$zero[2] * age - fit$coefficients$zero[3] * as.numeric(ethnic) -fit$coefficients$zero[4] *BMIAZ))
# mean(hurdle_zero_prob)
# sd(hurdle_zero_prob)
 
 

##########################################jags model#########


#Defining initial values for two chains
N1=length(data$count)
data1=list("N" = N1,"count"=count,"age"=age,
           "ethnic"=ethnic,"BMIAZ"=BMIAZ)
cl <- makeCluster(2)
ethnic_effect <-  c(NA, 1)
zero_ethnic_effect <-  c(NA, 1)
non_zero_group <- rep(1, N1)

inits=list("k"=1,"intercept"=0,"ethnic_effect"=ethnic_effect,
           "age_coefficient"= 1,"BMIAZ_coefficient"=1,
           "inflation_intercept"=0,"zero_ethnic_effect"=zero_ethnic_effect,
           "zero_age_coefficient"= -1,"zero_BMIAZ_coefficient"=1,
           "zero_inflation_intercept"=-1,"tau" =1,
           "non_zero_group"=non_zero_group)


inits1=rep(list(inits),2)

# The main function
jags1=run.jags("ZINB_tau_simlu_model.txt",burnin =2000,sample =6000,adapt =2000,n.chains=2,
               inits=inits1,monitor=c("k","intercept",
                                      "ethnic_effect", "age_coefficient",
                                      "BMIAZ_coefficient",
                                      "zero_ethnic_effect",
                                      "zero_inflation_intercept",
                                      "zero_age_coefficient",
                                      "zero_BMIAZ_coefficient", "tau",
                                      "non_zero_propotion","resid.sum.sq"),
               data=data1,summary=runjags.getOption("force.summary"))



# Export the summary of the estimations from the main function
S1 = add.summary(jags1,summary.iters=10000)
S1

summary(S1$mcmc)

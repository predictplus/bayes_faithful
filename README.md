# bayes_faithful
The capstone project for the UCSC on the Bayesian Statistics

# Executive Summary
There has been an original problem of the eruption time estimation given the waiting time over 100 minutes for the Old Faithful geyser in Yellowstone National Park, Wyoming, USA. Two MCMC models have been created to estimate the linear regression model parameters: eruption time = b0 + b1 * waiting time. The best performing model has been selected out of two. The selected model considers the fact that the geyser operates in 2 modes: long wait and short wait, and thus has bi-modal distribution of eruption time and waiting time. With the selected model it has been estimated that for the waiting time of 100 minutes the eruption time mean will be around 5.3 minutes.

# Introduction
The problem I’d like to research is the dependence of a geyser eruption time on the time between eruptions. Specifically, I’m interested in understanding eruption time, given we succeed in keeping the eruption delayed over 100 minutes (beyond usual waiting time)?

# Data
An existing dataset has been used for these purposes, called “faithful”. The dataset connects a waiting time between eruptions and the duration of the eruption for the Old Faithful geyser in Yellowstone National Park, Wyoming, USA. The data format is the following:

272 observations on 2 variables.
[,1]
eruptions
numeric
Eruption time in mins
[,2]
waiting
numeric
Waiting time to next eruption (in mins)

The dataset is perfect to answer the question since it has all necessary variables and may help to create required models and answer the questions. The dataset was apparently collected through visual observations when both eruption time and waiting time were measured by a person.
Based on the data analysis it appears that there are no missing values
> sum(is.na(dat))
[1] 0

Let’s take a look at the barplot of the histograms for the eruption time and waiting time to see how their distributions may potentially look like:
> hist(dat$eruptions)
> hist(dat$waiting)
> plot(density(dat$eruptions))
> plot(density(dat$waiting))
And then let’s check if there is any correlation between these values:
> plot(dat$eruptions~dat$waiting)
Both distributions seem to be bimodal and there is likely a linear correlation between these two variables from the faithful dataset.

We can also check the statistical summary for this dataset:
> summary(dat)
   	eruptions        waiting    
 Min.   :1.600   Min.   :43.0  
 1st Qu.:2.163   1st Qu.:58.0  
 Median :4.000   Median :76.0  
 Mean   :3.488   Mean   :70.9
 3rd Qu.:4.454   3rd Qu.:82.0  
 Max.   :5.100   Max.   :96.0 
 Model

Linear regression may be a great fit to answer the question we have in Clause 1, since there is a good linear correlation between both variables. Inference for regression Bayesian statistical model will allow us to essentially calculate any conditional probability between both variables given any input data.
We can first check the parameters of a simple linear regression model
> rm = lm(eruptions ~ waiting, data = dat)
> summary(rm)
Coefficients:
             Estimate Std. Error t value
(Intercept) -1.874016   0.160143  -11.70
waiting      0.075628   0.002219   34.09
            Pr(>|t|)    
(Intercept)   <2e-16 ***
waiting       <2e-16 ***
This gives us a good idea about potential means and standard deviations for each coefficient that we might further use for sampling.
The original model that we will build will assume the eruption time variable (y) comes from normal distribution with mean mu[i], and the mean comes from linear model of waiting time. Betas (coefficients of the linear model) will come from a normal distribution, with mean 0 and low variance. There will also be a precision parameter for the normal distribution which is inverse of sigma.
mod_string = "model {
     for (i in 1:n) {
         y[i] ~ dnorm(mu[i],prec)
         mu[i] = b[1] + b[2]*waiting[i]
     }
     for (j in 1:2) {
         b[j] ~ dnorm(0.0, 1.0/1.0e-4)
     }
     prec ~ dgamma(5.0/2.0, 5.0*10.0/2)
     sig2 = 1.0/prec
     sig = sqrt(sig2)
 }"

The original model looks like this

Running simulations for 5000 iterations with the 30000 iterations burn-out gives us the following plots. One can notice that overall there is a great convergence.

Let’s check the diagnostics for this model.
> gelman.diag(mod_sim)
Potential scale reduction factors:
     Point est. Upper C.I.
b[1]          1          1
b[2]          1          1
sig           1          1

Multivariate psrf
1
> autocorr.diag(mod_sim)
                b[1]        b[2]          sig
Lag 0   1.0000000000 1.000000000  1.000000000
Lag 1   0.0539659455 0.040366671  0.004508573
Lag 5  -0.0204763129 0.001507012 -0.003019463
Lag 10  0.0007562073 0.005117313 -0.004754853
Lag 50 -0.0055873482 0.012087695  0.009611340

> effectiveSize(mod_sim)
    b[1]     b[2]      sig 
13779.15 13833.54 15000.00 

> summary(mod_sim)

Iterations = 30001:35000
Thinning interval = 1 
Number of chains = 3 
Sample size per chain = 5000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

         Mean        SD  Naive SE Time-series SE
b[1] -0.00264 0.0099846 8.152e-05      8.508e-05
b[2]  0.04998 0.0006329 5.168e-06      5.381e-06
sig   0.73981 0.0314688 2.569e-04      2.569e-04

2. Quantiles for each variable:
         2.5%       25%       50%      75%   97.5%
b[1] -0.02223 -0.009357 -0.002633 0.004155 0.01684
b[2]  0.04873  0.049557  0.049987 0.050399 0.05120
sig   0.68150  0.717999  0.738788 0.760508 0.80507
Let’s check the residuals for the model.
> X = cbind(rep(1.0,data_jags$n), data_jags$waiting)
> pm_params = colMeans(mod_csim)
> yhat = drop(X %*% pm_params[1:2])
> resid = data_jags$y - yhat
> plot(resid)


Looks not too bad but there is still a high spread of residuals and also a clear dependence between yhat and resid. Let’s see if we can do better. One thing that can be improved is the model itself. Instead of sampling y[i] from a single distribution, we can sample from 2 different distributions with different mu[i], where bettas in mu[i] depend on a category of the eruption time. From the dat$waiting distribution we can see that it is a mixture of two normal distributions. The border between two is roughly at 66 minutes, and we can consider all waiting times below that as ones that belong to cat1 and the ones above that as ones that belong to cat2. 
> dat$cat <- ifelse(dat$waiting > 66, 1, 2)
The modified model will look like:
mod_string2 = "model {
     for (i in 1:n) {
         y[i] ~ dnorm(mu[i],prec)
         mu[i] = b0[cat[i]] + b1[cat[i]]*waiting[i]
     }
    
     for (j in 1:2) {
     	   b0[j] ~ dnorm(0.0, 1.0/1.0e-4)
     	   b1[j] ~ dnorm(0.0, 1.0/1.0e-4)
     }
    
     prec ~ dgamma(5.0/2.0, 5.0*10.0/2)
     sig2 = 1.0/prec
     sig = sqrt(sig2)
}"

Running simulations for 5000 iterations with the 30000 iterations burn-out gives us the following plots. One can notice that there is a great convergence.


> gelman.diag(mod_sim2)
Potential scale reduction factors:
      Point est. Upper C.I.
b0[1]          1          1
b0[2]          1          1
b1[1]          1          1
b1[2]          1          1
sig            1          1

Multivariate psrf
1

> autocorr.diag(mod_sim2)
              b0[1]        b0[2]        b1[1]        b1[2]
Lag 0   1.000000000  1.000000000  1.000000000  1.000000000
Lag 1   0.056251694  0.032628255  0.048036996  0.036583247
Lag 5  -0.001103798  0.005962964  0.003064703  0.002319884
Lag 10  0.001661488 -0.005541951  0.005795987 -0.006343716
Lag 50 -0.011053228 -0.002139406 -0.010210652 -0.005099828
                sig
Lag 0   1.000000000
Lag 1   0.000516131
Lag 5  -0.005622860
Lag 10 -0.011329846
Lag 50  0.003164589

> effectiveSize(mod_sim2)
   b0[1]    b0[2]    b1[1]    b1[2]      sig 
13695.11 14358.52 13950.23 13939.28 15000.00

> summary(mod_sim2)
Iterations = 30001:35000
Thinning interval = 1 
Number of chains = 3 
Sample size per chain = 5000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:
           Mean        SD  Naive SE Time-series SE
b0[1] 0.0013021 0.0100410 8.198e-05      8.591e-05
b0[2] 0.0008544 0.0100566 8.211e-05      8.393e-05
b1[1] 0.0531115 0.0005726 4.675e-06      4.851e-06
b1[2] 0.0377775 0.0010929 8.924e-06      9.258e-06
sig   0.5947281 0.0253338 2.068e-04      2.069e-04

2. Quantiles for each variable:
          2.5%       25%       50%      75%   97.5%
b0[1] -0.01863 -0.005462 0.0012445 0.008002 0.02060
b0[2] -0.01854 -0.005971 0.0008735 0.007639 0.02046
b1[1]  0.05197  0.052729 0.0531170 0.053491 0.05425
b1[2]  0.03565  0.037036 0.0377873 0.038508 0.03991
sig    0.54713  0.577199 0.5939694 0.611297 0.64689
Let’s check the residuals for the model now.
> X2 = cbind(rep(1.0,data_jags2$n), data_jags2$waiting)
> pm_params2 = colMeans(mod_csim2)
> 
> yhat2 <- rep(0, nrow(dat))
> 
> for (k in 1:nrow(dat)) {
+   if (dat$cat[k] == 1) {
+     yhat2[k] = as.numeric(X2[k,] %*% pm_params2[c(1,3)])
+     resid2[k] = data_jags2$y[k] - yhat2[k]
+   } else {
+     yhat2[k] = as.numeric(X2[k,] %*% pm_params2[c(2,4)])
+     resid2[k] = data_jags2$y[k] - yhat2[k]
+   }
+ }
> 
> plot(resid2)
> points(resid, col = "red")
> plot(yhat2, resid2)

The results are way better now. Firstly, the residuals have become way smaller. Secondly, the linear pattern between yhat2 and resid2 has almost disappear

 Results and Conclusions
Now let’s try to estimate the eruption time, given we succeed in keeping the eruption delayed over 100 minutes.

> wtg = 100
> err = mod_csim2[,"b0[1]"] + mod_csim2[,"b1[1]"] * wtg
> plot(density(err))
> mean(err)
[1] 5.31245
The eruption time is expected to be around 5.3 minutes. Below is the PDF for the eruption time given the waiting time is 100 minutes.


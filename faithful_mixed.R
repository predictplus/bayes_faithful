library("rjags")

dat = faithful

# mod_string1 = "model {
#     for (i in 1:n) {
#         y[i] ~ dnorm(mu[i], prec)
#         z[i] ~ dcat(omega)
#         mu[i] = b[z[i]*2-1] + b[z[i]*2]*waiting[i]
#     }
#     for (j in 1:4) {
#         b[j] ~ dnorm(0.0, 1.0/1.0e6)
#     }
#     prec ~ dgamma(5.0/2.0, 5.0*10.0/2)
#     sig2 = 1.0 / prec
#     sig = sqrt(sig2)
#     omega ~ ddirch(c(1.0, 1.0))
# }"

set.seed(1)

mod_string1 = "model {
    for (i in 1:n) {
        y[i] ~ dnorm(mu[i], prec)
        z[i] ~ dcat(omega)
        mu[i] = b[z[i]*2-1] + b[z[i]*2]*theta[z[i]]
    }
    
    theta[1] ~ dnorm(mu_theta[1], prec_theta[1])
    theta[2] ~ dnorm(mu_theta[2], prec_theta[2])
    
    sigma1 = 5
    mu1 = 53.5
    sigma2 = 7
    mu2 = 80
    mu_theta[1] ~ dnorm(mu1, 1.0/sigma1^2) 
    mu_theta[2] ~ dnorm(mu2, 1.0/sigma2^2)  
    prec_theta[1] ~ dgamma(1.0/2.0, 1.0*sigma1^2/2) 
    prec_theta[2] ~ dgamma(1.0/2.0, 1.0*sigma2^2/2)
    
    for (j in 1:4) {
        b[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    prec ~ dgamma(1.0/2.0, 1.0*1.0/2)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
    omega ~ ddirch(c(1.0, 1.0))
}"

data_jags1 = list(y = dat$eruptions, n = nrow(dat), waiting = dat$waiting)
params1 = c("b", "sig", "omega")
inits1 = list("b" = rnorm(4, 0.0, 100.0), "prec" = rgamma(1, 1.0, 1.0))

mod1 = jags.model(textConnection(mod_string1), data = data_jags1, inits = inits1, n.chains = 3)
update(mod1, 1000)
mod_sim1 = coda.samples(model=mod1, variable.names = params1, n.iter = 3000)
mod_csim1 = do.call(rbind, mod_sim1)
plot(mod_sim1)

autocorr.diag(mod_sim1)
effectiveSize(mod_sim1)
summary(mod_sim1)
gelman.diag(mod_sim1)

X = cbind(rep(1.0,data_jags1$n), data_jags1$waiting)
pm_params1 = colMeans(mod_csim1)
yhat1 = drop(X %% pm_params1[1:2])
resid1 = data_jags1$y - yhat1
plot(resid1)

wtg = 100
err = mod_csim1[,"b[1]"] + mod_csim1[,"b[2]"] * wtg
print(err)
```

```R

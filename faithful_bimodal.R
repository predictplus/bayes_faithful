library("rjags")

dat = faithful

dat$cat <- ifelse(dat$waiting > 66, 1, 2)

datcat1 <- subset(dat, cat == 1)
datcat2 <- subset(dat, cat == 2)

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

set.seed(11)

data_jags2 = list(y = dat$eruptions, waiting = dat$waiting, n = nrow(dat), cat = dat$cat)
params2 = c("b0", "b1","sig")
inits2 = list("prec" = rgamma(1, 1.0, 1.0))

mod2 = jags.model(textConnection(mod_string2), data = data_jags2, inits = inits2, n.chains = 3)
update(mod2, 30000)
mod_sim2 = coda.samples(model=mod2, variable.names = params2, n.iter = 5000)
mod_csim2 = do.call(rbind, mod_sim2)
plot(mod_sim2)

gelman.diag(mod_sim2)
autocorr.diag(mod_sim2)
effectiveSize(mod_sim2)
raftery.diag(mod_sim2)
summary(mod_sim2)

X2 = cbind(rep(1.0,data_jags2$n), data_jags2$waiting)
pm_params2 = colMeans(mod_csim2)

yhat2 <- rep(0, nrow(dat))

for (k in 1:nrow(dat)) {
  if (dat$cat[k] == 1) {
    yhat2[k] = as.numeric(X2[k,] %*% pm_params2[c(1,3)])
    resid2[k] = data_jags2$y[k] - yhat2[k]
  } else {
    yhat2[k] = as.numeric(X2[k,] %*% pm_params2[c(2,4)])
    resid2[k] = data_jags2$y[k] - yhat2[k]
  }
}

plot(resid2)
points(resid, col = "red")
plot(yhat2, resid2)

wtg = 100
err = mod_csim2[,"b0[1]"] + mod_csim2[,"b1[1]"] * wtg

quantile(err, c(0.025, 0.975))

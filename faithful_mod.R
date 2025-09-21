library("rjags")

dat = faithful

dat$cat <- ifelse(dat$waiting > 66, 1, 2)

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

data_jags = list(y = dat$eruptions, n = nrow(dat), waiting = dat$waiting)
params = c("b", "sig")
inits = list("b" = rnorm(2, 0.0, 100.0), "prec" = rgamma(1, 1.0, 1.0))

set.seed(11)

mod = jags.model(textConnection(mod_string), data = data_jags, inits = inits, n.chains = 3)
update(mod, 30000)
mod_sim = coda.samples(model=mod, variable.names = params, n.iter = 5000)
mod_csim = do.call(rbind, mod_sim)
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)
summary(mod_sim)

X = cbind(rep(1.0,data_jags$n), data_jags$waiting)
pm_params = colMeans(mod_csim)
yhat = drop(X %*% pm_params[1:2])
resid = data_jags$y - yhat
plot(resid)
library(tidyverse)
library(readxl)
library(cowplot)
library(rstan)
library(loo)
library(survival)
library(xtable)
rstan_options(auto_write = TRUE)

logmeanexp <- function(x){
  tmp <- max(x)
  tmp + log(sum(exp(x-tmp))) - log(length(x))
}

mod_wei <- stan_model("weibull.stan")
mod_gam <- stan_model("gamma.stan")
mod_lnorm <- stan_model("lognormal.stan")

#Data is available in Backer et al (2020)  
data_raw <- read_tsv(file = "20-00062_BAKER_SupplementMaterial/BAKER_Supplementary material S1_data.tsv")

data2 <- data_raw %>%
  mutate(exposure = as.integer(as.Date(exposure_end,format="%m/%d/%Y")-as.Date(exposure_start,format="%m/%d/%Y")),
         incubation = as.integer(as.Date(symptom_onset,format="%m/%d/%Y") - as.Date(exposure_end,format="%m/%d/%Y"))) %>% 
  dplyr::filter(incubation > 0)

sf1 <- survfit(Surv(time=incubation,time2=exposure+incubation, type="interval2")~1, data=data2)

dat4stan <- list(
  N=nrow(data2),
  d=ifelse(is.na(data2$exposure),0,ifelse(data2$exposure==0,1,2)),
  tE=ifelse(is.na(data2$exposure),0,data2$exposure),
  tI=data2$incubation
)

smp_wei <- sampling(mod_wei,dat4stan,seed=1)
traceplot(smp_wei)
summary(smp_wei)$summary
ex_wei <- rstan::extract(smp_wei)

smp_gam <- sampling(mod_gam,dat4stan,seed=1)
traceplot(smp_gam)
summary(smp_gam)$summary
ex_gam <- rstan::extract(smp_gam)

smp_lnorm <- sampling(mod_lnorm,dat4stan,seed=1)
traceplot(smp_lnorm)
summary(smp_lnorm)$summary
ex_lnorm <- rstan::extract(smp_lnorm)


r_eff_wei <- relative_eff(exp(ex_wei$log_lik),chain_id = rep(1:4, each = 1000))
r_eff_gam <- relative_eff(exp(ex_gam$log_lik),chain_id = rep(1:4, each = 1000))
r_eff_lnorm <- relative_eff(exp(ex_lnorm$log_lik),chain_id = rep(1:4, each = 1000))

loo_wei <- loo(extract_log_lik(smp_wei), r_eff = r_eff_wei)
loo_gam <- loo(extract_log_lik(smp_gam), r_eff = r_eff_gam)
loo_lnorm <- loo(extract_log_lik(smp_lnorm), r_eff = r_eff_lnorm)

dists <- c("Weibull","gamma","log-normal")
loo3 <- c(loo_wei$estimates[3,1],
loo_gam$estimates[3,1],
loo_lnorm$estimates[3,1])

#show Table 1
print(xtable(data.frame(distribution = dists, looic=loo3)),
      include.rownames = FALSE)

p_wei_pred <- function(t){
  exp(logmeanexp(pweibull(t,ex_wei$m,ex_wei$eta,lower.tail = FALSE, log.p = TRUE)))
}

p_gam_pred <- function(t){
  exp(logmeanexp(pgamma(t,ex_gam$a,ex_gam$b,lower.tail = FALSE, log.p = TRUE)))
}

p_lnorm_pred <- function(t){
  exp(logmeanexp(plnorm(t,ex_lnorm$mu,ex_lnorm$sigma,lower.tail = FALSE, log.p = TRUE)))
}

p_wei_pred_v <- Vectorize(p_wei_pred)
p_gam_pred_v <- Vectorize(p_gam_pred)
p_lnorm_pred_v <- Vectorize(p_lnorm_pred)

#produce Figure S3 
pdf("survfit.pdf")
plot(sf1,conf.int=FALSE,lty=2,xlim=c(0,30))
curve(p_wei_pred_v(x),add=TRUE, col="orange2", lwd=1.5, lty=3)
curve(p_gam_pred_v(x),add=TRUE, col="forestgreen", lwd=1.5, lty=4)
curve(p_lnorm_pred_v(x),add=TRUE, col="royalblue", lwd=1.5, lty=5)
legend("topright",c("non-parametric","Weibull","gamma","log-normal"),lty = c(2:5),lwd="1.5",
       col=c("black","orange2","forestgreen","royalblue"))
dev.off()

dfpred <- data.frame(dist=rep(factor(dists,levels = dists),each=4000),
                     pred=c(pred_wei,pred_gam,pred_lnorm))

#show Table 2
print(xtable(
t(simplify2array(tapply(dfpred$pred, dfpred$dist, quantile, prob=c(0.025,0.5,0.975))))
))

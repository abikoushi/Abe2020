library(rstan)
library(loo)
library(cowplot)
rstan_options(auto_write = TRUE)
mod_weibull <- stan_model("weibull_intcens.stan") #marginalized
mod_weibull2 <- stan_model("weibull_intcens2.stan") #Backer's method

logmeanexp <- function(x){
  tmp <- max(x)
  tmp + log(sum(exp(x-tmp))) - log(length(x))
}
logdiffexp <- function(x,y){
  x + log1p(-exp(y-x))
}

simfunc <- function(N,iter=100){
  loores <- numeric(iter)
  loores2 <- numeric(iter)
  for (i in 1:iter) {
    set.seed(i)
    x <- rweibull(N,2,2)
    y <- floor(x)
    dat <- list(N=N,y=y)
    smp_weibull <- sampling(mod_weibull,dat,cores=4,seed=1)
    smp_weibull2 <- sampling(mod_weibull2,dat,cores=4,seed=1)
    ex <- rstan::extract(smp_weibull)
    r_eff <- relative_eff(exp(ex$log_lik),chain_id = rep(1:4, each = 1000))
    est_loo <- loo::loo(ex$log_lik,r_eff=r_eff)$estimates[3,1]/(2*N)
    ex2 <- rstan::extract(smp_weibull2)
    r_eff2 <- relative_eff(exp(ex2$log_lik),chain_id = rep(1:4, each = 1000))
    est_loo2 <- loo::loo(ex2$log_lik,r_eff=r_eff2)$estimates[3,1]/(2*N)
    logdist_pred <- function(t){
      logmeanexp(logdiffexp(pweibull(t+1,ex$m,ex$eta,log.p = TRUE),
                            pweibull(t,ex$m,ex$eta,log.p = TRUE)))
    }
    logdist_pred_v <- Vectorize(logdist_pred)
    logdist_pred2 <- function(t){
      logmeanexp(dweibull(t,ex$m,ex$eta,log = TRUE))
    }
    logdist_pred_v2 <- Vectorize(logdist_pred2)
    dist_true <- function(t){
      pweibull(t+1,2,2)-pweibull(t,2,2)
    }
    dist_true2 <- function(t){
      dweibull(t,2,2)
    }
    f_ge <- function(t){
      -logdist_pred_v(t)*dist_true(t)
    }
    f_ge2 <- function(t){
      -logdist_pred_v2(t)*dist_true2(t)
    }
    ge_v <- f_ge(0:100)
    GE2 <- integrate(f_ge2,0,Inf)
    loores[i] <- est_loo-(sum(ge_v, na.rm = TRUE))
    loores2[i] <- est_loo2-(GE2$value)
  }
  return(list(loo=loores,loo2=loores2))
}
out25 <- simfunc(25)
out50 <- simfunc(50)
out100 <- simfunc(100)
# save(out25,out50,out100,file = "sim_Weibull_intcens.Rdata")

df <- data.frame(method = c(rep("marginalized",100*3),rep("Backer",100*3)),
                 n = factor(rep(c(rep(25,100),rep(50,100),rep(100,100)),2)),
                 value =c(out25$loo,out50$loo,out100$loo,
                          out25$loo2,out50$loo2,out100$loo2))

mean_sd <- function (x, mult = 1){
  x <- stats::na.omit(x)
  sd <- mult * stats::sd(x)
  mean <- mean(x)
  ggplot2:::new_data_frame(list(y = mean, ymin = mean - sd, ymax = mean + 
                        sd), n = 1)
}

#show Figure S3
ggplot(df,aes(x=n,y=value,colour=method,linetype=method,
              group=interaction(method,n)))+
  geom_violin()+
  stat_summary(geom = "pointrange", fun.data  = mean_sd, position = position_dodge(width=0.9))+
  geom_hline(yintercept = 0,linetype=2)
# ggsave("simGEplot.pdf")

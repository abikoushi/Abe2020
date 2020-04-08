data{
  int<lower = 1> N;
  int<lower = 0, upper=2> d[N];
  real tI[N];
  real tE[N];
}

parameters{
  real<lower = 0> m;
  real<lower = 0> eta;
}
model{
  for(i in 1:N){
    if(d[i]==1){
     target += weibull_lpdf(tI[i]|m,eta); 
    }else if(d[i]==2){
      target += log_diff_exp(weibull_lcdf(tE[i]+tI[i]|m,eta),weibull_lcdf(tI[i]|m,eta));
    }else{
     target += weibull_lccdf(tI[i]|m,eta); 
    }
  }
}
generated quantities{
  real meanhat;
  real pred;
  real log_lik[N];
  for(i in 1:N){
    if(d[i]==1){
      log_lik[i] = weibull_lpdf(tI[i]|m,eta); 
    }else if(d[i]==2){
      log_lik[i] = log_diff_exp(weibull_lcdf(tE[i]+tI[i]|m,eta),weibull_lcdf(tI[i]|m,eta));
    }else{
      log_lik[i] = weibull_lccdf(tI[i]|m,eta); 
    }
  }
  meanhat = eta*tgamma(1+1/m);
  pred = weibull_rng(m,eta);
}

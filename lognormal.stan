data{
  int<lower = 1> N;
  int<lower = 0, upper=2> d[N];
  real tI[N];
  real tE[N];
}

parameters{
  real<lower = 0> sigma;
  real mu;
}
model{
  for(i in 1:N){
    if(d[i]==1){
     target += lognormal_lpdf(tI[i]|mu,sigma); 
    }else if(d[i]==2){
      target += log_diff_exp(lognormal_lcdf(tE[i]+tI[i]|mu,sigma),lognormal_lcdf(tI[i]|mu,sigma));
    }else{
     target += lognormal_lccdf(tI[i]|mu,sigma); 
    }
  }
}
generated quantities{
  real meanhat;
  real pred;
  real log_lik[N];
  for(i in 1:N){
    if(d[i]==1){
      log_lik[i] =  lognormal_lpdf(tI[i]|mu,sigma); 
    }else if(d[i]==2){
      log_lik[i] = log_diff_exp(lognormal_lcdf(tE[i]+tI[i]|mu,sigma),lognormal_lcdf(tI[i]|mu,sigma));
    }else{
      log_lik[i] = lognormal_lccdf(tI[i]|mu,sigma); 
    }
  }
  meanhat = exp(mu+sigma^2/2);
  pred = lognormal_rng(mu,sigma);
}

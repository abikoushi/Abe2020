data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real<lower=0> m;
  real<lower=0> eta;
}
model {
  for(i in 1:N){
    if(y[i]==0.0){
     target += weibull_lcdf(1.0 | m, eta);
    }else{
     target += log_diff_exp(weibull_lcdf(y[i]+1.0| m, eta),weibull_lcdf(y[i]| m, eta)); 
    }
  }
}
generated quantities{
  vector[N] log_lik;
    for(i in 1:N){
    if(y[i]==0.0){
      log_lik[i] = weibull_lcdf(1.0 | m, eta);
    }else{
      log_lik[i] = log_diff_exp(weibull_lcdf(y[i]+1.0| m, eta),weibull_lcdf(y[i]| m, eta));
    }
  }
}

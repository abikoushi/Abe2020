data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real<lower=0> m;
  real<lower=0> eta;
  vector<lower=0,upper=1>[N] ep;
}
model {
  for(i in 1:N){
     target += weibull_lpdf(y[i]+ep[i] | m, eta);
  }
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = weibull_lpdf(y[i]+ep[i] | m, eta);
  }
}

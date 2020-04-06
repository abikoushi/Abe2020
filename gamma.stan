data{
  int<lower = 1> N;
  int<lower = 0, upper=2> d[N];
  real tI[N];
  real tE[N];
}

parameters{
  real<lower = 0> a;
  real<lower = 0> b;
}
model{
  for(i in 1:N){
    if(d[i]==1){
     target += gamma_lpdf(tI[i]|a,b); 
    }else if(d[i]==2){
      target += log_diff_exp(gamma_lcdf(tE[i]+tI[i]|a,b),gamma_lcdf(tE[i]|a,b));
    }else{
     target += gamma_lccdf(tI[i]|a,b); 
    }
  }
}
generated quantities{
  real meanhat;
  real pred;
  real log_lik[N];
  for(i in 1:N){
    if(d[i]==1){
     log_lik[i] = gamma_lpdf(tI[i]|a,b); 
    }else if(d[i]==2){
      log_lik[i]  = log_diff_exp(gamma_lcdf(tE[i]+tI[i]|a,b),gamma_lcdf(tE[i]|a,b));
    }else{
     log_lik[i] = gamma_lccdf(tI[i]|a,b); 
    }
  }
  meanhat = a/b;
  pred = gamma_rng(a,b);
}

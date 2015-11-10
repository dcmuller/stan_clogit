## clogit.stan
## Conditional logistic regression 
## (the same likelihood calculation as used by survival::clogit in R and 
## clogit in Stata)
## David C Muller

functions {
  ## function to return the number of observatios in a group
  int group_size(vector ref, int value) {
    int toselect[rows(ref)] ;
    for (ii in 1:rows(ref))
      toselect[ii] <- ref[ii]==value;
    return sum(toselect);
  }
  
  ## function to return the number of cases or events in a set
  int n_cases(vector y) {
    int toselect[rows(y)] ;
    for (ii in 1:rows(y))
      toselect[ii] <- y[ii]==1;
    return sum(toselect);
  }
  
  ## function to subset a vector (return just those observations in a given group)
  vector subset_vector(vector y, vector ref, int value) {
    int jj;
    int toselect[rows(y)];
    vector[group_size(ref, value)] res;
    for(ii in 1:rows(y))
      toselect[ii] <- ref[ii] == value;
    jj <- 1;
    for(ii in 1:rows(y)) {
      if (toselect[ii] == 1) {
        res[jj] <- y[ii];
        jj <- jj+1;
      }
    }
    return res;
  }
  
  ## recursive function to evaluate the denominator of the conditional likelihood 
  real cl_denom(int N_g, int D_g, vector xb); 
  real cl_denom(int N_g, int D_g, vector xb) {
    real res;
    if (N_g < D_g) {
      return 0;
    }
    if (D_g == 0) {
      return 1;
    }
    res <- cl_denom(N_g-1, D_g, xb) + cl_denom(N_g-1, D_g-1, xb)*exp(xb[N_g]);
    return res;
  }
}

data {
  int<lower=0> N;
  int<lower=1> n_grp;
  int<lower=1> n_coef;
  int<lower=1, upper=n_grp> grp[N];
  vector[N] y;
  matrix[N, n_coef] x;
}
transformed data {
 int n_group[n_grp]; # number of observations in the group
 int n_case[n_grp]; # number of cases/events in the group
 for (ii in 1:n_grp) {
   n_group[ii] <- group_size(to_vector(grp), ii);
   n_case[ii] <- n_cases(subset_vector(y, to_vector(grp), ii));
 }
}
parameters {
  vector[n_coef] b;
}

transformed parameters {
  vector[n_coef] oddsratio;
  oddsratio <- exp(b);
}
model {
  vector[N] xb; # observation level linear predictor
  real ll; # log likelihood
  
  ## diffuse normal prior for log odds ratios
  b ~ normal(0, 3);
  
  ## log likelihood is a sum over each group
  xb <- x * b;
  for (ii in 1:n_grp) {
    vector[n_group[ii]] y_g;
    vector[n_group[ii]] xb_g;
    y_g <- subset_vector(y, to_vector(grp), ii);
    xb_g <- subset_vector(xb, to_vector(grp), ii);
    ll <- dot_product(y_g, xb_g) - log(cl_denom(n_group[ii], n_case[ii], xb_g));
    increment_log_prob(ll);
  }
}

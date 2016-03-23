## clogit.stan
## Conditional logistic regression
## (the same likelihood calculation as used by survival::clogit in R and
## clogit in Stata)
## David C Muller

functions {
  ## function to return the number of observations in a group
  int group_size(int[] ref, int value) {
    int count;
    count <- 0;
    for (ii in 1:size(ref))
      if (ref[ii]==value)
        count <- count + 1;
    return count;
  }

  ## function to subset a vector (return just those observations in a given group)
  vector subset_vector(vector y, int[] ref, int value) {
    int jj;
    vector[group_size(ref, value)] res;
    if (size(ref) != rows(y))
      reject("illegal input: non-matching dimensions")
    jj <- 1;
    for(ii in 1:size(ref)) {
      if (ref[ii] == value) {
        res[jj] <- y[ii];
        jj <- jj+1;
      }
    }
    return res;
  }

  ## function to subset an integer array (return just those observations in a given group)
  int[] subset_intarray(int[] y, int[] ref, int value) {
    int jj;
    int res[group_size(ref, value)];
    if (size(ref) != size(y))
      reject("illegal input: non-matching dimensions")
    jj <- 1;
    for(ii in 1:size(ref)) {
      if (ref[ii] == value) {
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
    res <- cl_denom(N_g-1, D_g, xb) + exp(log(cl_denom(N_g-1, D_g-1, xb)) + xb[N_g]);
    return res;
  }
}

data {
  int<lower=0> N; # Number of observations
  int<lower=1> n_grp; # Number of groups
  int<lower=1> n_coef; # Number of coefficients (log odds ratios) to estimate
  int<lower=1, upper=n_grp> grp[N]; # stratum/group identifier
  int<lower=0, upper=1> y[N]; # array of 0/1 outcomes
  matrix[N, n_coef] x; # Matrix of regressors
}

transformed data {
 int n_group[n_grp]; # number of observations in the group
 int n_case[n_grp]; # number of cases/events in the group
 for (ii in 1:n_grp) {
   n_group[ii] <- group_size(grp, ii);
   {
     int subset_y[n_group[ii]];
     subset_y <- subset_intarray(y, grp, ii);
     n_case[ii] <- group_size(subset_y, 1);
   }
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
    int y_g[n_group[ii]];
    vector[n_group[ii]] xb_g;
    y_g <- subset_intarray(y, grp, ii);
    xb_g <- subset_vector(xb, grp, ii);
    ll <- dot_product(to_vector(y_g), xb_g) - log(cl_denom(n_group[ii], n_case[ii], xb_g));
    increment_log_prob(ll);
  }
}

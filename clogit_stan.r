## example conditional logistic regression using Stan
## David C Muller

library(survival)
library(rstan)

set.seed(77834)

## use the infertility data from the survival package
datlist <- list(N=nrow(infert), 
                n_grp=max(infert[, "stratum"]), 
                n_coef=2,
                x=infert[,c("spontaneous", "induced")], 
                y=infert[, "case"],
                grp=infert[, "stratum"])

## fit using default parameters
clogit_stan <- stan("clogit.stan", data=datlist)
clogit_stan

## fit using survival::clogit for comparison
clogit_survival <- clogit(case ~ spontaneous + induced + strata(stratum), 
                          data=infert)
summary(clogit_survival)

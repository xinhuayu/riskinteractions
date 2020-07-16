#####################
# exploring risk interactions of covid-19 in South Korea
# study period: March 10 to April 30, 2020
# revised date: July 14, 2020
# author: Xinhua Yu
#####################
rm(list = ls())

library(foreign)
library(dplyr)
library(data.table)
library(incidence)
library(mgcv)
library(ggplot2)


# data source, updated monthly;
# https://www.kaggle.com/kimjihoo/coronavirusdataset

rootpath = "C:/Users/xinhuayu/Desktop/COVID19/southkorea/"

covid_kr=data.table(read.csv(paste0(rootpath,"TimeAge.csv")))

covid_kr[,agegroup:=recode(age, "0s"=0,"10s"=0,"20s"=1,"30s"=1,"40s"=2,"50s"=2,"60s"=3,"70s"=3,"80s"=3)]
# counts are already cumulative sum, pool over age groups by date
covid_kr[,c("cases","deaths"):=list(sum(confirmed),sum(deceased)), by=c("agegroup","date")]

# obtain the last obs for each age and date group
covid_kr2<-covid_kr[, .SD[.N], by=c("agegroup","date"),.SDcols=c("cases","deaths")]
covid_kr2[,time:=seq_len(.N),by="agegroup"]
# sort the data
setkeyv(covid_kr2,c("agegroup","date"))
# new cases and deaths;
covid_kr2[,c("newcases","newdeaths"):=list(cases-shift(cases),deaths-shift(deaths)), by=c("agegroup")]

# remove new case in the first date in each age group (missing)
# create a date variable: case_date
covid_kr2<-covid_kr2[!is.na(newcases),][,case_date:=as.Date(date)]
setkeyv(covid_kr2,c("agegroup","case_date"))

####################################
# plot the epidemic curve by age groups
# South Korea epidemic starts on Jan 20, but real boost about Feb 20
# the age date starts at March 2, 2020
agedaily = covid_kr2[case_date >=as.Date("2020-03-10") & case_date <as.Date("2020-05-01"),]

age10daily = agedaily[agegroup==0,]
age30daily = agedaily[agegroup==1,]
age50daily = agedaily[agegroup==2,]
age60daily = agedaily[agegroup==3,]

# padding NA to make the same length to all age group data
maxlen<-max(table(agedaily$agegroup))
age10daily = age10daily[1:maxlen,]
age30daily = age30daily[1:maxlen,]
age50daily = age50daily[1:maxlen,]
age60daily = age60daily[1:maxlen,]

########################################
# function for obtaining the predicted cases
myfitgam <- function(grpdata){
  grpdata[,time:=seq_len(.N)]
  setnames(grpdata,"time","datevalue")
  
  Mixmodel <- gam(newcases~s(datevalue,bs="tp",k=16,m=2),family=poisson,
                  data=grpdata,method="REML")
  
  print(summary(Mixmodel))
  print(anova(Mixmodel))
  
  newd <-data.frame(datevalue=grpdata$datevalue)
  
  # prediction + SE from NB model 
  fitTN <- predict( Mixmodel ,newd, type="response",se = TRUE )
  fitTN2 <- data.table(cbind(fitTN,newd))
  return(fitTN2)
}

# separate model
gfitTN10 = myfitgam(age10daily)
gfitTN30 = myfitgam(age30daily)
gfitTN50 = myfitgam(age50daily)
gfitTN60 = myfitgam(age60daily)

# dates from one age group 
gfitdata <- data.frame(newdate=as.Date(age30daily$case_date),
                       fitTN10 = gfitTN10$fit,
                       seTN10 = gfitTN10$se.fit,
                       fitTN30 = gfitTN30$fit,
                       seTN30 = gfitTN30$se.fit,
                       fitTN50 = gfitTN50$fit,
                       seTN50 = gfitTN50$se.fit,
                       fitTN60 = gfitTN60$fit,
                       seTN60 = gfitTN60$se.fit,
                       Dailycase10=age10daily$newcases,
                       Dailycase30=age30daily$newcases,
                       Dailycase50=age50daily$newcases,
                       Dailycase60=age60daily$newcases)

# export figure at 800 x 436
ggplot(gfitdata) + 
  geom_point(aes(newdate, Dailycase10, color="0-19")) + geom_line(aes(newdate, Dailycase10,color="0-19"),linetype="dashed") +
  geom_point(aes(newdate, Dailycase30, color="20-39")) + geom_line(aes(newdate, Dailycase30,color="20-39"),linetype="dashed") +
  geom_point(aes(newdate, Dailycase50, color="40-59")) + geom_line(aes(newdate, Dailycase50,color="40-59"),linetype="dashed") +
  geom_point(aes(newdate, Dailycase60, color="60+")) + geom_line(aes(newdate, Dailycase60,color="60+"),linetype="dashed") +
  geom_line(aes(newdate, fitTN10,color="0-19"),size = 1,linetype="solid") +
  geom_line(aes(newdate, fitTN30,color="20-39"),size = 1,linetype="solid") +
  geom_line(aes(newdate, fitTN50,color="40-59"),size = 1,linetype="solid") +
  geom_line(aes(newdate, fitTN60,color="60+"),size = 1,linetype="solid") +
  scale_color_manual("Age groups",values=c("0-19" = "green","20-39"="blue","40-59"="orange","60+"="red")) + 
  labs(x = "Date",y="Daily New Cases") + scale_x_date(date_breaks="1 week",date_labels = "%b %d") + 
  theme(legend.position=c(.72,.65)) +
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 


############### time series analysis ##################
# starting from the end of first wave epidemic, March 10, 2020
covid_kr4<-covid_kr2[case_date>=as.Date("2020-03-10") & case_date< as.Date("2020-05-01"),]
setkeyv(covid_kr4,c("agegroup","case_date"))

# age distribution of total cases, the last max obs
covid_summ = covid_kr4[,max(cases),by="agegroup"]
covid_summ[,sum:=sum(V1)][,pct:=V1/sum]
covid_summ[]

# deaths
covid_summ = covid_kr4[,max(deaths),by="agegroup"]
covid_summ[,sum:=sum(V1)][,pct:=V1/sum]
covid_summ[]


# reshape data to create counts for each age group;
covid_kr5<-dcast(covid_kr4,case_date~agegroup,value.var = "newcases")
setnames(covid_kr5,c("0","1","2","3"),c("age10","age30","age50","age60p"))

covid_kr5[]

# library(foreign)
# plot the epidemic curve in stata;
# write.dta(covid_kr5,"C:/Users/xinhuayu/Desktop/COVID19/southkorea/covid_kr5.dta")

# library(ggplot2)
# plot the epidemic curves
# covid_kr4 %>% 
#   ggplot( aes(x=date, y=newcases, group=agegroup, color=as.factor(agegroup), linetype=as.factor(agegroup))) +
#   geom_line()


# need to use baysian methods with rstan
library(rstan)
library(bayesplot)

rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

# old package, forced multicore chain
library(rstanmulticore)

options(max.print = 4000)
options(scipen = 9999)

######full model with multiplicative scale ####################
# assume the correlation between series is latent (random effect in mixed model)
# also normal priors for regression coefficients;
##############################################################

##############################################################
# Standard Poisson model with covariance between series
##############################################################

ptcode <- '
data {
  int<lower=1> N;   // length of time series
  int<lower=1> K;   // K  of lags
  int<lower=1> J;   // J time series, J equations
  int y[N,J];       // N obs, J time series
}

transformed data {
  vector[J] zeros;  // J zeros for multinormal prior for random effects
  matrix[N,J] lny;
  int ones[N,J];
  
  zeros = to_vector(rep_array(0,J));
  
  ones = rep_array(1,N,J);
  lny = log(to_matrix(y)+to_matrix(ones)); //add ones to avoid log(0)
}

parameters {
  vector[J] alpha;   // J equations
  vector[J] bt;    // J random effects
  real beta[J,J,K]; // J equations, J cross predictors, and K lags 
  corr_matrix[J] Omega; 
  vector<lower=0>[J] sigma; //SD
}

transformed parameters {
  cov_matrix[J] Sigma;  //covariance of random effects
  Sigma = quad_form_diag(Omega, sigma);   
}

// add covariance between these processes
model {
  real lambda[N,J];
  
  // Weakly informative prior, around means or zero etc
    for (j in 1:J) {
      alpha[j] ~ student_t(5,0,2.5);
      for (m in 1:J) {
        for (k in 1:K) 
          beta[j,m,k] ~ student_t(5,0,2.5);
      }
    }

  // covariance of latent variables (similar to random effect in a mixed model)
  sigma ~ cauchy(0,5); // prior on the standard deviations
  Omega ~ lkj_corr(2); // LKJ prior on the correlation matrix, diagnal modality density
  
   
  bt ~ multi_normal(zeros,Sigma);
  
  // could be vectorized  
  for (j in 1:J) {
    for (n in K+1:N){
      lambda[n,j] = alpha[j] + bt[j];
      for (m in 1:J) {
        for (k in 1:K) 
          lambda[n,j] += beta[j,m,k]*lny[n-k,m];
      }
    }
    y[K+1:N,j] ~ poisson_log(lambda[K+1:N,j]);
  }
}

generated quantities {
  real exp_beta[J,J,K];
  exp_beta = exp(beta);
}
'
# default 5 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 5,        # no of lags, 3,5,7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
                )

fit5 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,         # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50             # sampling slice per chain
)

fit5

saveRDS(fit5,paste0(rootpath,"PO_lag5_fit5.RDS"))
#total elapse fit1
print(get_elapsed_time(fit5))

# model diagnostic check 
# library(shinystan)
# launch_shinystan(fit5)

# 3 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 3,        # no of lags, 3,5,7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)


fit5 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,         # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50             # sampling slice per chain
)

fit5
saveRDS(fit5,paste0(rootpath,"PO_lag3_fit5.RDS"))

#total elapse fit1
print(get_elapsed_time(fit5))

# 7 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 7,        # no of lags, 3,5,7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit5 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,         # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50             # sampling slice per chain
)

fit5
saveRDS(fit5,paste0(rootpath,"PO_lag7_fit5.RDS"))

#total elapse fit1
print(get_elapsed_time(fit5))



################negative binomial model##################;
# slower mixing, need more iters, wider confidence interval
##########################################################
ptcode <- '
data {
  int<lower=1> N;   // length of time series
  int<lower=1> K;   // K  of lags
  int<lower=1> J;   // J time series, J equations
  int y[N,J];       // N obs, J time series
}

transformed data {
  vector[J] zeros;  // J zeros for multinormal prior for random effects
  matrix[N,J] lny;
  int ones[N,J];
  
  zeros = to_vector(rep_array(0,J));
  ones = rep_array(1,N,J);
  lny = log(to_matrix(y)+to_matrix(ones));  //add ones to log to avoid log(0)
}

parameters {
  vector[J] alpha;   // J equations
  vector[J] bt;    // J random effects
  real beta[J,J,K]; // J equations, J cross predictors, and K lags 
  corr_matrix[J] Omega; 
  vector<lower=0>[J] sigma; //SD
  real<lower=0> invsqrtphi[J];  //negative binomial, NB(mu, phi), scale parameter, var(x) = mu^2/phi 
}

transformed parameters {
  cov_matrix[J] Sigma;  //Variance of random effects
  Sigma = quad_form_diag(Omega, sigma);   
}

// add covariance between these processes
model {
  real lambda[N,J];
  real phi[J];
  
  // weakly informative prior, could be vectorized
    for (j in 1:J) {
      alpha[j] ~ student_t(5,0,2.5);
      invsqrtphi[j] ~ normal(0,0.3); // half normal prior for NB scale,  make it closer to Poisson 
      phi[j] = 1/(invsqrtphi[j]^2);  //convert inverse of square root of phi for proper prior;
      for (m in 1:J){
        for (k in 1:K) 
          beta[j,m,k] ~ student_t(5,0,2.5);
      }
    }
    
  // covariance of latent variables (similar to random effect in a mixed model)
  sigma ~ cauchy(0,5); // prior on the standard deviations
  Omega ~ lkj_corr(2); // LKJ prior on the correlation matrix, diagnal modal distributed
  
  // combine alpha, bt  
  bt ~ multi_normal(zeros,Sigma);
  
  // could be vectorized  
  for (j in 1:J) {
    for (n in K+1:N){
      lambda[n,j] = alpha[j] + bt[j];
      for (m in 1:J) {
        for (k in 1:K) 
          lambda[n,j] += beta[j,m,k]*lny[n-k,m];
      }
    }
    // alternative parameterization, nb(mu, phi) where phi is a scale parameter for overdispersion
    y[K+1:N,j] ~ neg_binomial_2_log(lambda[K+1:N,j],phi[j]);
  }
}

generated quantities {
  real exp_beta[J,J,K];
  real disp[J];
  exp_beta = exp(beta);
  for (j in 1:J) disp[j] = invsqrtphi[j]^2;
}
'
# default 5 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 5,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit7 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 82000,        # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # forced initial values
)

fit7

saveRDS(fit7,paste0(rootpath,"NB_lag5_fit7.RDS"))

#library(shinystan)
#launch_shinystan(fit7)

# 3 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 3,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit7 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 82000,        # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # forced initial values
)

fit7

saveRDS(fit7,paste0(rootpath,"NB_lag3_fit7.RDS"))


# 7 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 7,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit7 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 82000,        # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # forced initial values
)

fit7

saveRDS(fit7,paste0(rootpath,"NB_lag7_fit7.RDS"))


###############Generialized poisson distribution #############;
# the pdf of standard poisson is:  (note: gamma(y+1) = y!) 
#          {lambda^y*exp(-lambda)}/gamma(y+1)
# or       {lambda*lambda^(y-1)*exp(-lambda)}/gamma(y+1)     
#
# the pdf of generalized poisson is (comparing the above poisson pdf, introducing a scale factor xi):
#          {lambda*(lambda+xi*y)^(y-1) *exp(-(lambda+xi*y))}/gamma(y+1)
#
# them mean is mu = lambda/(1-xi); and variance is var = lambda/(1-xi)^3 = mu/(1-xi)^2
# the dispersion factor is 1/(1-xi)^2  
#
# reparametralize the pdf of generalized poisson in terms of mu and xi
#          {mu*(1-xi) * (mu - xi*(mu-y))^(y-1)*exp(-(mu-xi*(mu-y)))}/gamma(y+1)
# the log likelihood of new pdf is:
#          log(mu*(1-xi)) + (y-1)*log(mu-xi*(mu-y)) - (mu - xi*(mu-y)) - log(gamma(y+1))
#
# we can still model log(mu) = beta*X, and let xi be the scale parameter, and 1/(1-xi)^2 the dispersion factor
# when xi = 0, then the generalized poisson is standard poisson
# when 0 < xi < 1, then overdispersion
# when xi < 0, then underdispersion  (uncommon, but for autocorrelated counts, maybe)
# we will write our own pdf and loglikelihood
# they are slow mixing, not optimized and not in CPP
##############################################################;

ptcode <- '
functions {
  real genpoisson_lg_log(int[] y,real[] lgmu,real xi) {
      // log likelihood, log mu input
      real gplp;
      vector[num_elements(y)] mu;
      vector[num_elements(y)] muy;
    
      gplp = 0;
      mu = exp(to_vector(lgmu));
      muy = mu-xi*(mu-to_vector(y));

      // no need the lgamma(y+1) in the log likelihood, they are constant;
      for (i in 1:num_elements(y)){
         gplp += log(mu[i]*(1-xi)) + (y[i]-1)*log(muy[i]) - muy[i] ;
      }
      return gplp;     
    }
}

data {
  int<lower=1> N;   // length of time series
  int<lower=1> K;   // K  of lags
  int<lower=1> J;   // J time series, J equations
  int y[N,J];       // N obs, J time series
}

transformed data {
  vector[J] zeros;  // J zeros for multinormal prior for random effects
  matrix[N,J] lny;
  int ones[N,J];
  
  zeros = to_vector(rep_array(0,J));
  
  ones = rep_array(1,N,J);
  lny = log(to_matrix(y)+to_matrix(ones)); //add ones to avoid log(0)
}

parameters {
  vector[J] alpha;   // J equations
  vector[J] bt;    // J random effects
  real beta[J,J,K]; // J equations, J cross predictors, and K lags 
  corr_matrix[J] Omega; 
  vector<lower=0>[J] sigma; //SD
  real<lower=-0.9, upper=0.9> xi[J];  // -0.9 is about 25% underdisperson, and 0.9 is about 100 times overdispersion
}

transformed parameters {
  cov_matrix[J] Sigma;  //covariance of random effects
  Sigma = quad_form_diag(Omega, sigma);   
}

// add covariance between these processes
model {
  real lambda[N,J];

  // Weakly informative prior, around means or zero etc
    for (j in 1:J) {
      alpha[j] ~ student_t(5,0,2.5);
      xi[j] ~ normal(0,0.3);  //xi can be negative
      for (m in 1:J) {
        for (k in 1:K) 
          beta[j,m,k] ~ student_t(5,0,2.5);
      }
    }

  // covariance of latent variables (similar to random effect in a mixed model)
  sigma ~ cauchy(0,5); // prior on the standard deviations
  Omega ~ lkj_corr(2); // LKJ prior on the correlation matrix, diagnal modality density
  
  bt ~ multi_normal(zeros,Sigma);
  
  for (j in 1:J) {
    for (n in K+1:N){
      lambda[n,j] = alpha[j] + bt[j];
      for (m in 1:J) {
        for (k in 1:K) 
          lambda[n,j] += beta[j,m,k]*lny[n-k,m];
      }
    }
    y[K+1:N,j] ~ genpoisson_lg(lambda[K+1:N,j],xi[j]);  // functions take vectors
  }
}

generated quantities {
  real exp_beta[J,J,K];
  real disp[J];
  exp_beta = exp(beta);
  for (j in 1:J) 
    disp[j] = 1/((1-xi[j])^2);  //dispersions
}
'
# default 5 lag models, reported in the paper

kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 5,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
                )

fit9 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # set initial values
)

fit9

saveRDS(fit9,paste0(rootpath,"GP_lag5_fit9.RDS"))

#total elapse fit1
print(get_elapsed_time(fit9))

# check models
# library(shinystan)
# launch_shinystan(fit9)

# 3 lag model
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 3,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit9 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # set initial values
)

fit9

saveRDS(fit9,paste0(rootpath,"GP_lag3_fit9.RDS"))

# 7 lag file
kr_age_data<-list(N = length(covid_kr5$age60p),
                  K = 7,        # no of lags, 3, 5, 7
                  J = 4,        # no of time series
                  y = covid_kr5[,c("age60p","age50","age30","age10")]
)

fit9 <- pstan(
  model_code = ptcode,  # Stan program code
  data = kr_age_data,   # named list of data
  chains = 5,           # number of Markov chains
  warmup = 2000,        # number of warmup iterations per chain
  iter = 52000,         # total number of iterations per chain
  thin = 50,            # sampling slice per chain
  init = 0              # set initial values
)

fit9

saveRDS(fit9,paste0(rootpath,"GP_lag7_fit9.RDS"))


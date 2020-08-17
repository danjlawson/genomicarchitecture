## Usage: "Rscript sim_bayesinf.R compile" to make the require stan models
## Compile once and make sure "model.RDS" and "model3.RDS" are in your directory
## Then "Rscript sim_bayesinf.R <S> <Fst> <seed>" for integer seeds to run the stan models
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Date: August 2020
## Licence: GPLv3

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
compile=FALSE

args=commandArgs(TRUE)
if(length(args)>0) {
    if(args[[1]]=="compile"){
        compile=TRUE
    }else if(length(args)==3){
        S=as.numeric(args[[1]])
        Fst=as.numeric(args[[2]])
        seed=as.numeric(args[[3]])
        set.seed(seed)
        print(seed)
    }else{
      stop("Usage: Rscript sim_bayesinf.R <S> <Fst> <seed>")
    }
}
#stop()

## Define our model
## This model accounts for genetic architecture but not uncertainty
simple_model_code <-
'
data {
  int<lower=0> J; // number of snps
  real beta[J]; // effect sizes
  real<lower=0> f[J]; // minor allele frequences
}
parameters {
  real S;
  real<lower=0> sigma;
}
transformed parameters {
  real snpsd[J];
  for (j in 1:J)
    snpsd[j] = pow(f[j]*(1-f[j]),S/2)*sigma;
}
model {
    S ~ uniform(-2,2);
    sigma ~ uniform(0,2);
    for(j in 1:J)
        beta[j] ~ normal(0,snpsd[j]);
}
'

## Define our model
## This model accounts for genetic architecture and uncertainty using Fst
drift_model_code <-
'
data {
  int<lower=0> J; // number of snps
  real beta[J]; // effect sizes
  real<lower=0> f[J]; // minor allele frequences
  real<lower=0> Fst; // Known value of Fst
}
parameters {
  real S; 
  real<lower=0> sigma_b;
  real<lower=0,upper=0.5> p[J]; // true minor allele frequences
}
transformed parameters {
}
model {
    S ~ uniform(-2,2);
    sigma_b ~ uniform(0,2);
    for(j in 1:J){
        f[j] ~ beta( ((1-Fst)/Fst)*p[j], ((1-Fst)/Fst)*(1-p[j]) );
        beta[j] ~ normal(0,pow(p[j]*(1-p[j]),S/2)*sigma_b);
    }
}
'


## Create the stan model compiled object
if(compile){
    sm <- stan_model(model_code = simple_model_code)
    saveRDS(sm,"simple_model.RDS")
    ##sm2 <- stan_model(model_code = model2_code)
    smd <- stan_model(model_code = drift_model_code)
    saveRDS(smd,"drift_model.RDS")
    stop("Finished compiling code")
}else{
    sm=readRDS("simple_model.RDS")
    smd=readRDS("drift_model.RDS")
}
#####################

## create data in the form stan likes, from a 2 column data frame of beta and f
data_frame_to_stan_list=function(data,use="obs",Fst=NULL){
    if(use=="obs"){
        ret=list(J=dim(data)[1],beta=data[,"beta"],f=data[,"f"])
    }else if(use=="true"){
        ret=list(J=dim(data)[1],beta=data[,"b"],f=data[,"p"])
    }else stop("invalid \"use\"")
    if(!is.null(Fst)) ret$Fst=Fst
    ret
}

## Make example test data
make_test_data=function(N,sigma_b=0.01,sigma_beta=0,sigma_f=0.1,S=-1){
    p=runif(2*N,0.01,0.5)
    b=rnorm(2*N,0,sigma_b*(p*(1-p))^(S/2))
    if(sigma_f==0){ f=p
    }else f=rbeta(length(p),(1-sigma_f)*p/sigma_f,(1-sigma_f)*(1-p)/sigma_f) # Balding and Nichols model

    ok=which((f>0.01)&(p>0.01))
    if(sum(ok)>=N) {
        ok=sample(ok,N)
    }else ok=sample(ok,N,replace=TRUE)
    f[f>0.5]=1-f[f>0.5]
    f[f<0.01]=0.01
    if(sigma_beta==0){
      beta=b
    }else{
      beta=rnorm(N,b,sigma_beta)
    }
    r=data.frame(f=f,beta=beta,p=p,b=b)
    r[ok,]
}
## Extract MCMC samples from a stan results object: get the one called S
getS=function(stanres){
    svals=as.numeric(extract(stanres,"S")[[1]])
    c(mean=mean(svals),quantile(svals,c(0.05,0.25,0.5,0.75,0.95)))
}


#########
#########
#########
## START OF DATA GENERATION

## Some parameters we don't need to change
thin=500 # Report mcmc samples after this many steps
N=10000 # Number of SNPs
iter=10000 # Number of MCMC iterations (increase if you have convergence problems)

## Make appropriate data
test=make_test_data(N,sigma_b=0.01,sigma_f=Fst,S=S)
## Where the f and beta are taken from the "true" simulated f and beta's, "in Africa" These are for when we have "direct access" to the data.
data_direct=data_frame_to_stan_list(test,"true")
## Where the f and beta are taken from the drifted f.
data_obs=data_frame_to_stan_list(test,"obs",Fst=Fst)


#########
## START OF INFERENCE

## Infer in the drift model, observed and direct datasets
smd_obs<-sampling(smd,data=data_obs,chains=2,iter=iter,thin=thin)
test_obs<-sampling(sm,data=data_obs,chains=2,iter=iter,thin=thin)
test_direct<-sampling(sm,data=data_direct,chains=2,iter=iter,thin=thin)

## Report output
allres=list(
    res=cbind(direct=getS(test_direct),
              full=getS(smd_obs),
              obs=getS(test_obs)),
        data=test,
        obs=test_obs,
        direct=test_direct,
        smd=smd_obs
)


############
## Write results to disk
save(allres, file=paste0('s_',S,'_Fst',Fst,'_seed',seed,'.RData'))

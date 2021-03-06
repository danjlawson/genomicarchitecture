## Usage: "Rscript sim_bayesinf.R compile" to make the require stan models
## Compile once and make sure "model.RDS" and "model3.RDS" are in your directory
## Make sure that "all_maf.RData" is in your directory, for generating data from 1kg
## Then "Rscript sim_bayesinf_1kg.R <S> <pop> <Fst> <seed>" for integer seeds to run the stan models
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
    }else if(length(args)==4){
        S=as.numeric(args[[1]])
        target=args[[2]]
        Fst=as.numeric(args[[3]])
        seed=as.numeric(args[[4]])
        set.seed(seed)
        print(seed)
    }else{
      stop("Usage: Rscript sim_bayesinf.R <S> <seed>")
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
    saveRDS(sm4,"drift_model.RDS")
    stop("Finished compiling code")
}else{
    sm=readRDS("simple_model.RDS")
    smd=readRDS("drift_model.RDS")
}

##########################################
## Load the SNP frequencies created from hapgen
snpfreqs=load("all_maf.RData")
snpfreqs=get(snpfreqs)
ok=apply(snpfreqs[,2:6]>0,1,all) ## Restrict to those SNPs that are not at zero frequency

## create data in the form stan likes, from a 2 column data frame of beta and f
data_frame_to_stan_list=function(data,use="obs",Fst=NULL){
    if(use=="obs"){
        data=data[data[,"f"]>0,]
        ret=list(J=dim(data)[1],beta=data[,"beta"],f=data[,"f"])
    }else if(use=="true"){
        data=data[data[,"p"]>0,]
        ret=list(J=dim(data)[1],beta=data[,"b"],f=data[,"p"])
    }else stop("invalid \"use\"")
    if(!is.null(Fst)) ret$Fst=Fst
    ret
}

## Make example test data
make_test_data=function(N,snpfreqs0,snpfreqs1,sigma_b=0.01,sigma_beta=0,S=-1){
    tdf=data.frame(p=snpfreqs0,f=snpfreqs1)
    tdf=tdf[tdf$p>0.01,]
    mysnps=sample(1:dim(tdf)[1],N)
    tdf=tdf[mysnps,]
    tdf$b=rnorm(N,0,sigma_b*(tdf$p*(1-tdf$p))^(S/2))
    if(sigma_beta==0){
      tdf$beta=tdf$b
    }else{
      tdf$beta=rnorm(N,tdf$b,sigma_beta)
    }
    tdf
}

## Extract MCMC samples from a stan results object: get the one called S
getS=function(stanres){
    svals=as.numeric(extract(stanres,"S")[[1]])
    c(mean=mean(svals),quantile(svals,c(0.05,0.25,0.5,0.75,0.95)))
}


#########
thin=500
N=10000
iter=10000
test=make_test_data(N,snpfreqs[,"afr"],snpfreqs[,target],sigma_b=0.01,S=S)
data_obs=data_frame_to_stan_list(test,"obs",Fst=Fst)
data_direct=data_frame_to_stan_list(test,"true")


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

save(allres, file=paste0('s_',S,'_target',target,'_Fst',Fst,'_seed',seed,'_1kg.RData'))

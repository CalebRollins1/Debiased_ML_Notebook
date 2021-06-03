########
# set up
########

rm(list=ls())

library("foreign")
library("dplyr") 
library("ggplot2")
library("quantreg") #used for rq.fit.sfn
library("nnet")	#used for mulitnom
library("randomForest")
library("keras")
library("mvtnorm")
library("Rlab")

setwd("~/Downloads/rrr_401k_PCR")

#######################
# clean and format data
#######################

source('specifications.R')




df  <- read.dta("sipp1991.dta")
spec=3 #spec in (1-3)
#quintile=0 #quintile in (1-5). 0 means all quintiles
data<-get_data(df,spec,0) #trimming like Farrell; different than Chernozhukov et al. 

Y.vH=data[[1]]
T.vH=data[[2]]
X.vH=data[[3]] #no intercept

##################
# helper functions
##################

source('primitives.R')
source('stage1.R')





###########
# algorithm
###########

set.seed(1) # for sample splitting

alpha_estimator=5
gamma_estimator=7
bias=0
#alpha_estimator: 0 dantzig, 1 lasso, 2 PCR, 3 PCR2, 4 PCR3 (no dict), 5 PCR.missing
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn, 4 PCR, 5 PCR2, 6 PCR3 (no dict), 7 PCR.missing

source('stage2.R')

coverage_experiment = function(data,n,noise.func,noise.params,iter = 500,rs,target,m = NULL,biased = FALSE,alt_params = NULL,CIs,alpha_est = 2, gamma_est = 4){
  ###Inputs
  #data: either a data generating function that takes as input n or list(Y,T,X)
  #n: dimension of X  if data is a function(X is n x n) (irrelevant is data is a list)
  #noise.func: function that takes as inputs (X and an element from noise.params to add
    #noise to X
  #noise.params: list of parameters to itererate over. 
  #iter: number of iterations each coverage experiment is run over
  #rs: list of pcs to try when running PCR
  #target: true value (H0) you're testing to see if it falls in the confidence interval
  #m: not ijmplemented
  #biased: indicates whether or not to run a trial of PCR with biasing
  #alt_params: list of alpha and gamma values to try
  #CIs: list of confidence intervals to test for (as percents, so 95% CI = 95)
  #alpha_est: alpha for PCR (number 2-5)
  #gamma_est: gamma for PCR (number 4-7)
  
  ###Outputs
  #results: data frame where each row gives data for 1 coverage experiment
  
  
  ###Designed to run coverage experiments that vary across a number of dimensions. These include:
  ###1. Amount of noise of added to X
  ###2. number of principal components kept for PCR
  ###3. biased and unbiased models
  ###4. other models that aren't PCR (such as lasso or nn)
  
  ###The code does the following:
  ###run a coverage experiment for each combination of noise, principal components, and biases:
  ### For iter rounds, randomly draw data (if data is a function) and add random noise
  ###   Run PCR with alpha_est and gamma_est as your alpha and gamma estimators
  ###   Record from the trial 1.ATE 2.SE 3.which confidence intervals had target in them
  ###   Run all of your alternate parameters and record same as above
  ###Each row in the output corresponds to 1 combination of noise, principal components, and biases
  ###As 1 row in the data frame, record:
  ###   1. the ATE averaged over all trials
  ###   2. the SE average over all trials
  ###   3. the coverage of each confidence interval in CIs
  
  #checks if data is a function or list
  if (length(data)!=1){
    Y = data[[1]]
    T = data[[2]]
    X = data[[3]]
  } else {
    D = data(n)
    Y = D[[1]]
    T = D[[2]]
    X = D[[3]]
  }
  
  #Set up initial paremeters for debiased ML
  # dictionary
  dict=b # b for partially linear model, b2 for interacted model. note that b2 appears in stage1.R for NN
  p=length(b(T[1],X[1,]))
  
  #p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  p0=ceiling(p/4) 
  if (p>60){
    p0=ceiling(p/40)
    
  }
  
  
  D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
  D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
  max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights
  
  #num_of_specs is the number of different specifications for each noise value
  num_of_specs = length(rs)*(biased+1)+length(alt_params)
  #print(num_of_specs)
  params = c()
  biases = c()
  CI_mat = matrix(0,length(CIs),0)
  ATEs = c()
  SEs = c()
  algs = c()
  if (is.null(m)){
    m = n
  }
  
  #first loop; iterates over noise values
  for (param in noise.params){
    
    #Uses to keep track of ATE, SE, and how ofter CI contains target
    CI = matrix(0,length(CIs),num_of_specs)
    ATE_mat = matrix(0,iter,num_of_specs)
    SE_mat = matrix(0,iter,num_of_specs)
    #Coverage Experiment
    for (j in 1:iter){
      print(paste(param,j))
      
      #generates data if data is a function
      if (length(data)==1){
        D = data(n)
        Y = D[[1]]
        T = D[[2]]
        X = D[[3]]
      }
      
      X.noisy = noise.func(X,param) #add noise
      #loop over different ranks
      for (i in 1:length(rs)){
        r = rs[i]
        
        
        #results with debiasing
        unbiased_pcr<-rrr(Y,T,X.noisy,p0,D_LB,D_add,max_iter,dict,alpha_est,gamma_est,0,pcs = r)
        for (k in 1:length(CIs)){
          if (abs(unbiased_pcr[3]-target)<(qnorm(0.5*CIs[k]/100+0.5)*unbiased_pcr[4])){
            #CI[1] = CI[1]+1
            CI[k,i] = CI[k,i]+1
          }
        
        }
        ATE_mat[j,i] = unbiased_pcr[3]
        SE_mat[j,i] = unbiased_pcr[4]
        printer(unbiased_pcr)
        #biased results
        if (biased){
          biased_pcr<-rrr(Y,T,X.noisy,p0,D_LB,D_add,max_iter,dict,alpha_est,gamma_est,1,pcs = r)
          for (k in 1:length(CIs)){
            if (abs(biased_pcr[3]-target)<qnorm(0.5*CIs[k]/100+0.5)*biased_pcr[4]){
              #CI[1] = CI[1]+1
              CI[k,i+length(rs)] = CI[k,i+length(rs)]+1
            }
          }
          ATE_mat[j,i+length(rs)] = biased_pcr[3]
          SE_mat[j,i+length(rs)] = biased_pcr[4]
          printer(biased_pcr)
        }
        
        
      }
      #loop over all alternative parameters
      if (!is.null(alt_params)){
        i = length(rs)*(biased+1)
        for (spec in alt_params){
          i = i+1
          alpha = spec[[1]]
          gamma = spec[[2]]
          res = rrr(Y,T,X.noisy,p0,D_LB,D_add,max_iter,dict,alpha,gamma,0,pcs = 0)
          for (k in 1:length(CIs)){
            if (abs(res[3]-target)<qnorm(0.5*CI[k]/100+0.5)*res[4]){
              #CI[1] = CI[1]+1
              CI[k,i] = CI[k,i]+1
            }
            
          }
          ATE_mat[j,i] = res[3]
          SE_mat[j,i] = res[4]
          #if (j==1){algs = c(algs,paste(alpha,gamma,collapse = ' '))}
          printer(res)
        }
      }
      
    } #end of loop for coverage experiment
    
    #updating lists that will make up the result data frame
    params =c(params,rep(paste(param,collapse = ' '),num_of_specs))
    
    
    for (col in 1:num_of_specs){
      ATEs = c(ATEs,sum(ATE_mat[,col])/iter)
    }
    
    for (col in 1:num_of_specs){
      SEs = c(SEs,sum(SE_mat[,col]^2)^0.5/iter)
    }
    CI_mat = cbind(CI_mat,CI)
    if (biased){
      algs = c(algs,rep(paste("2 4",rs),2))
      biases = c(biases,c(rep(0,length(rs)),rep(1,length(rs)),rep(0,length(alt_params))))
    } else {
      algs = c(algs,rep(paste("2 4",rs),1))
      biases = c(biases,rep(1,num_of_specs))
    }
    
    for (elem in alt_params){
      algs = c(algs,paste(elem,collapse = ' '))
    }
  
  } #end of loop over noise parameters
  
  results = data.frame('noise' = params,
                          'ATE' = ATEs,
                          'SE' = SEs,
                          'bias'= biases,
                          'alg' = algs)
  #print(CI_mat)
  #add columns for each confidence level tested
  for (i in 1:length(CIs)){
    results[[paste("CI",CIs[i])]] = CI_mat[i,]/iter
  }
  
  return(results)
}

rd = function(mat,p){
  return(random_dropout(mat,p)[[1]])
}


#copy of simulate_data. Used to make is easier to change rank
sim = function(n){
  return(simulate_data(n,rank = 8))
}


#data: either a data generating function that takes as input n or list(Y,T,X)
#n: dimension of X  if data is a function(X is n x n) (irrelevant is data is a list)
#noise.func: function that takes as inputs (X and an element from noise.params to add
#noise to X
#noise.params: list of parameters to itererate over. 
#iter: number of iterations each coverage experiment is run over
#rs: list of pcs to try when running PCR
#target: true value (H0) you're testing to see if it falls in the confidence interval
#m: not ijmplemented
#biased: indicates whether or not to run a trial of PCR with biasing
#alt_params: list of alpha and gamma values to try
#CIs: list of confidence intervals to test for (as percents, so 95% CI = 95)
#alpha_est: alpha for PCR (number 2-5)
#gamma_est: gamma for PCR (number 4-7)

results_noise = coverage_experiment(data = sim, 
                    n = 500, 
                    noise.func = rd, 
                    noise.params = c(0.1,0.2,0.3), 
                    iter = 500,
                    rs = c(8,10,13),
                    target = 2.2,
                    m = NULL,
                    biased = TRUE,
                    alt_params = c(),
                    CIs = c(80,95),
                    alpha_est = 5,
                    gamma_est = 7)

results_dropout = coverage_experiment(data = list(Y.vH,T.vH,X.vH), #sim
                                    n = 500, 
                                    noise.func = rd, #add_noise
                                    noise.params = c(0.1,0.2,0.3), 
                                    iter = 100,
                                    rs = c(15),
                                    target = 8319.94,
                                    m = NULL,
                                    biased = FALSE,
                                    alt_params = c(), #list(list(1,1),list(1,3))
                                    CIs = c(80,95),
                                    alpha_est = 5,
                                    gamma_est = 7)



#write.csv(dropout_coverage,"file")










library(ggplot2)
library(dplyr)
library(foreign)
library(Rlab)
library(nnet)
library(tidyr)


setwd("~/Downloads/rrr_401k_blackbox")


df  <- read.dta("sipp1991.dta")


get_data<-function(df,spec,quintile){
  ###  Inputs
  ###  df: dataset that is being used. In this case a .dta but could be another format
  ###  spec: value of 1,2,or 3.
  ###  quintile: integer in [0,5]. 0 means all data, 1-5 indicates which quintile of income distribution to filter for
  ###  Output
  ###  List of (Y,D,X)
  Y <- df[,"net_tfa"]
  D <- df[,"e401"]
  #  All three specifications incorporate the same variables, but differ in the number of higher order terms. X.L has no higher order terms. X.H has higher order terms up to the 8th degree for some variables. X.vH incorporates higher order and interaction terms.
  ## low dim specification
  X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
  
  ## high dim specification.
  X.H <- cbind(poly(df[,"age"], 6, raw=TRUE),
               poly(df[,"inc"], 8, raw=TRUE),
               poly(df[,"educ"], 4, raw=TRUE),
               poly(df[,"fsize"], 2, raw=TRUE),
               df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
  
  ## very high dim specification
  X.vH=model.matrix(~(poly(df[,"age"], 6, raw=TRUE) +
                        poly(df[,"inc"], 8, raw=TRUE) +
                        poly(df[,"educ"], 4, raw=TRUE) +
                        poly(df[,"fsize"], 2, raw=TRUE) +
                        df[,"marr"] +
                        df[,"twoearn"] +
                        df[,"db"] +
                        df[,"pira"] +
                        df[,"hown"])^2)
  X.vH=X.vH[,-1]
  
  #Selects for desires specification
  if (spec==1){
    X=X.L
  } else if (spec==2) {
    X=X.H
  } else {
    X=X.vH
  }
  
  X <- scale(X,center=TRUE,scale=TRUE) #Centers and scales X
  n=nrow(X)
  #Impose common support. More discussion in Appendix B.
  p.1 <- multinom(D~X-1, trace=FALSE)$fitted.values
  indexes.to.drop <- which(p.1 < min(p.1[D==1]) | max(p.1[D==1]) < p.1)
  if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
  n.per.treatment <- as.vector(table(D[-indexes.to.drop]))
  n.trim <- n.per.treatment[1]+n.per.treatment[2]
  
  Y.trimmed=Y[-indexes.to.drop]
  D.trimmed=D[-indexes.to.drop]
  X.trimmed=X[-indexes.to.drop,]
  
  if (spec==1){
    inc=X.trimmed[,2]
  } else if (spec==2) {
    inc=X.trimmed[,7]
  } else {
    inc=X.trimmed[,7]
  }
  #Filters for desired quintile
  if (quintile>0){
    q <- ntile(inc, 5)
    Y.q=Y.trimmed[q==quintile]
    D.q=D.trimmed[q==quintile]
    X.q=X.trimmed[q==quintile,]
  } else {
    Y.q=Y.trimmed
    D.q=D.trimmed
    X.q=X.trimmed
  }
  
  return(list(Y.q,D.q,X.q))
  
}

X.L = get_data(df,1,0)[[3]]
X.H = get_data(df,2,0)[[3]]
X.vH = get_data(df,3,0)[[3]]

add_noise <- function(mat,var){
  #Adds gaussian noise to every entry in a matrix
  #
  #Inputs
  #mat: matrix that you want to add noise to
  #var: variance of the noise
  #
  #Ouput: matrix of same dimensions as mat with noise added
  N = dim(mat)[1]
  M = dim(mat)[2]
  noise = matrix(rnorm(N*M,mean=0,sd=var^0.5), N, M)
  return(mat+noise)
}

random_dropout <- function(mat,p){
  #Randomly makes each entry in a matrix NA with probability p
  #
  #Inputs
  #mat: matrix that you're changing
  #p: probability of an entry becoming NA
  #
  #Output[[1]]: mat with some entries as NA
  #Output[[2]]: matrix with same NAs as output[[1]] but all other entries 1
  N = dim(mat)[1]
  M = dim(mat)[2]
  dropout_matrix = matrix(rbern(N*M,1-p),N,M)
  dropout_matrix[dropout_matrix == 0] = NA
  
  return(list(dropout_matrix*mat,dropout_matrix))
}

fill.na <- function(X,return_scale = FALSE){
  #Replaces NAs with zeroes and scales columns
  #First step in matrix completion for the dropour matrix
  #
  #Inputs
  #X: matrix that is being filled in
  #return scale: returns vector that is used to scale each row
  #
  #Output: filled in matrix
  #Output[[2]] (if return_scale = TRUE): 1XN matrix with scaling ratios
  N = dim(X)[1]
  M = dim(X)[2]
  scale.ratios = 1-(t(is.na(X))%*%matrix(1,N,1))/N
  X.solution = X/t(matrix(scale.ratios,M,N))
  X.solution[is.na(X.solution)] = 0
  if (return_scale){
    return(list(X.solution,scale.ratios))
  }
  return(X.solution)
}

fill.na.custom = function(X,scale.ratios){
  #Replaces NAs with zeroes and scales by scale.ratios
  #Used on test data when you want to use training scale ratios
  #
  #Inputs
  #X: matrix being filled in
  #scale.ratios: ratios to scale each column by
  #
  #Output: filled in matrix
  N = dim(X)[1]
  M = dim(X)[2]
  X.solution = X/t(matrix(scale.ratios,M,N))
  X.solution[is.na(X.solution)] = 0
  return(X.solution)
}

denoise <- function(X,pcs,return_U = FALSE){
  #By truncating D matrix, attempts to clean data with gaussian noise
  #
  #Inputs
  #X: matrix that is being corrected
  #pcs: number of principal components to keep in D
  #return_U: returns 1:pcs columns of U from svd
  #
  #Output: corrected matrix
  #Output[[2]] (if return_U = TRUE): truncated U matrix
  decomp = svd(X)
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  X.fixed = U%*%D_trunc%*%t(V)
  
  if (return_U){
    return(list(X.fixed,U[,1:pcs]))
  }
  
  return(X.fixed)
  
}


MSE <- function(X,var,pcs){
  #Returns the square error for a noisy matrix
  #Replaced by MSE.dropout.fast
  #
  #Inputs
  #X: matrix that will have noise added to it
  #var: variance of noise added
  #pcs: number of pcs used to denoise X
  #
  #Output: Sqaured error between X and the denoised version of X
  Z = add_noise(X,var)
  decomp = svd(Z)
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  X.fixed = U%*%D_trunc%*%t(V)
  
  return(norm(X-X.fixed,type = "F"))
  
}

MSE.dropout <- function(X,p,pcs,na_only = FALSE){
  #Returns the square error for a matrix that has had some values dropped out
  #Replaced by MSE.dropout.fast
  #
  #Inputs
  #X: matrix that will have values dropped out
  #p: probability an entry becomes NA
  #pcs: number of pcs used to denoise X
  #
  #Output: Sqaured error between X and the denoised version of X
  new = random_dropout(X,p)
  na_matrix = new[[2]]
  Z = fill.na(new[[1]])
  decomp = svd(Z)
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  X.fixed = U%*%D_trunc%*%t(V)
  if (na_only){
    na_matrix[na_matrix == 1] = 0
    na_matrix[is.na(na_matrix)] = 1
    return(norm((X-X.fixed)*na_matrix,type = "F")^2/(norm((X)*na_matrix,type = "F")^2))
  }
  return(norm(X-X.fixed,type = "F")^2/norm(X,type = "F")^2)
  
}

NMSE_1 = function(X,X.fixed,na_matrix,na_only){
  #Error function that divides by square error by n*(max(X)-min(X))
  #
  #Inputs
  #X: Original clean matrix
  #X.fixed: fixed version of X after noise was added
  #na_matrix: matrix with NAs where values were dropped from X
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: error
  if(na_only){
    na_matrix.copy = na_matrix
    na_matrix.copy[na_matrix == 1] = 0
    na_matrix.copy[is.na(na_matrix)] = 1
    n = sum(na_matrix.copy)
    return(norm((X-X.fixed)*na_matrix.copy,type = "F")^2/(n*(max(abs(X*na_matrix),na.rm = TRUE)-min(abs(X*na_matrix),na.rm = TRUE))))
  } else {
    n = dim(X)[1]*dim(X)[2]
    return(norm(X-X.fixed,type = "F")^2/(n*(max(abs(X))-min(abs(X)))))
  }
}

NMSE_2 = function(X,X.fixed,na_matrix,na_only){
  #Error function that divides square error by norm(X)
  #
  #Inputs
  #X: Original clean matrix
  #X.fixed: fixed version of X after noise was added
  #na_matrix: matrix with NAs where values were dropped from X
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: error
  if(na_only){
    na_matrix[na_matrix == 1] = 0
    na_matrix[is.na(na_matrix)] = 1
    return(norm((X-X.fixed)*na_matrix,type = "F")^2/(norm((X)*na_matrix,type = "F")^2))
  } else {
    return(norm(X-X.fixed,type = "F")^2/norm(X,type = "F")^2)
  }
}

R2 = function(X,X.fixed,na_matrix,na_only){
  #Error function that divides square error by norm(X-bar_y) where bar_y is mean of X.fixed
  #
  #Inputs
  #X: Original clean matrix
  #X.fixed: fixed version of X after noise was added
  #na_matrix: matrix with NAs where values were dropped from X
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: error
  if(na_only){
    na_matrix[na_matrix == 1] = 0
    na_matrix[is.na(na_matrix)] = 1
    #print(X.fixed)
    #print(na_matrix)
    bar_y = sum(X.fixed*na_matrix)/sum(na_matrix)
    return(norm((X-X.fixed)*na_matrix,type = "F")^2/(norm((X-bar_y)*na_matrix,type = "F")^2))
  } else {
    bar_y = mean(X.fixed)
    return(norm(X-X.fixed,type = "F")^2/norm(X-bar_y,type = "F")^2)
  }
}

MSE.dropout.fast <- function(X,decomp,na_matrix,p,pcs,error_metric = "NMSE_1",na_only = FALSE){
  #Faster version of MSE.dropout
  #Sped up by passing in SVD instead of calculating it each time
  #
  #Inputs
  #X: matrix that had noise added to it
  #decomp: SVD of X after values were dropped out and columns scaled
  #na_matrix: matrix with NAs where values were dropped from X
  #p: probability an entry becomes NA
  #pcs: number of pcs used to denoise X
  #error_metric: One of 3 metrics NMSE_1, NMSE_2, R^2
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: error
  
  
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0 #could be sped up doing rank one updates
  
  
  if (error_metric == "NMSE_1"){
    error_func = NMSE_1
  } else if (error_metric == "NMSE_2"){
    error_func = NMSE_2
  } else if (error_metric == "R^2") {
    error_func = R2
  } else {
    error_func = NMSE_1
  }
  
  X.fixed = U%*%D_trunc%*%t(V)
  return(error_func(X,X.fixed,na_matrix,na_only))
  
}

xval.fast <- function(X,Y,Y.complete,decomp,na_matrix.X,na_matrix.Y,p,pcs,error_metric = "NMSE_1",na_only = FALSE){
  #Uses SVD from training set to complete test matrix and return errors for both
  #Sped up by passing in SVD instead of calculating it each time
  #
  #Inputs
  #X: matrix that had noise added to it
  #Y: test matrix with no noise
  #Y.complete: test matrix after noise has been added and columns rescaled
  #decomp: SVD of X after values were dropped out and columns scaled
  #na_matrix.X: matrix with NAs where values were dropped from X
  #na_matrix.Y: matrix with NAs where values were dropped from Y
  #p: probability an entry becomes NA
  #pcs: number of pcs used to denoise X
  #error_metric: One of 3 metrics NMSE_1, NMSE_2, R^2
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output[[1]]: error for training data
  #Output[[2]]: error for test data
  
  U = decomp$u
  D = decomp$d
  V = decomp$v
  U_trunc = U[,1:pcs] #first pcs columns of U
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0 #could be sped up doing rank one updates
  
  #selects error function and assigns to error_func
  if (error_metric == "NMSE_1"){
    error_func = NMSE_1
  } else if (error_metric == "NMSE_2"){
    error_func = NMSE_2
  } else if (error_metric == "R^2") {
    error_func = R2
  } else {
    error_func = NMSE_1
  }
  #svd.Y = svd(Y.complete)
  #U.Y = svd.Y$u
  #D.Y = svd.Y$d
  #V.Y = svd.Y$v
  X.fixed = U%*%D_trunc%*%t(V) #uses first pcs principa components of noisy X
  #X.fixed = U.Y[,1:pcs]%*%t(U.Y[,1:pcs])%*%Y.complete
  Y.fixed = U_trunc%*%t(U_trunc)%*%Y.complete 
  #print(Y.complete)
  #print("-------")
  #print(Y)
  #print(Y.fixed)
  return(list(error_func(X,X.fixed,na_matrix.X,na_only),error_func(Y,Y.fixed,na_matrix.Y,na_only)))
  
}


split_data = function(data,ratio = 0.5){
  #Splits data matrix randomly into to seperate matrices
  train.idx = sample.int(dim(data)[1],round(ratio*dim(data)[1]))
  train = data[train.idx,]
  test = data[-train.idx,]
  return(list(train,test))
}



xval_to_df.dropout = function(p_list,max_pcs,data,error_metric = "NMSE_1",na_only = FALSE){
  #Uses xval.dropout.fast to make a data frame with errors for different assignments of p
  #
  #Inputs
  #p_list: list of p values to use when dropping out values of X
  #max_pcs: maximum number of principal components to try when completing X
  #data: matrix that will have values dropped out
  #error_metric: one of NMSE_1, NMSE_2, R^2
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: data frame with 4 columns
  #p: probability of dropout
  #pcs: principal components used to approximate X
  #error.train: error on training data
  #error.test: error on test data
  
  #splits data into X and Y with equal number of rows
  split = split_data(data[1:(2*round(dim(data)[1]/2)),])
  X = split[[1]]
  Y = split[[2]]
  
  #print("X")
  #print(X)
  
  #print("Y")
  #print(Y)
  
  
  #X = matrix(c(rep(-1,100),rep(1,100),rep(1,100)),100)
  #variance_list = c(0.2,0.4,0.6,0.8,1)
  #max_zeros = 3
  
  #Initializes lists to store arguments passed to xval.dropout.fast
  #These each only need computed once for each p, so to speed up the process
  #we dropout entries and calculate svd once for each value of p
  vars = c() #holds values of p
  svds = list() #holds svd of X after dropour for each p in p_list
  Y.completes = list() #holds completed Y matrices using new scale ratios for each p
  na_matrices.X = list() #holds NA matrics for X after dropout for each p
  na_matrices.Y = list() #holds NA matrices for Y after dropour for each p
  for (i in 1:length(p_list)){
    vars = append(vars,rep(p_list[i],max_pcs)) #values from p_list repeated max_pcs times
    new.X = random_dropout(X,p_list[i])
    new.Y = random_dropout(Y,p_list[i])
    
    #print("new.X")
    #print(new.X)
    
    #print("new.Y")
    #print(new.Y)
    na_matrix.X = new.X[[2]]
    na_matrix.Y = new.Y[[2]]
    Z = fill.na(new.X[[1]],return_scale = TRUE)
    scaling = Z[[2]]
    Z = Z[[1]]
    decomp = svd(Z)
    Y.complete = fill.na.custom(new.Y[[1]],scale.ratios = scaling)
    
    #print("Y.comlpete")
    #print(scaling)
    #print(Y.complete)
    svds[[as.character(p_list[i])]] = decomp
    na_matrices.X[[as.character(p_list[i])]] = na_matrix.X
    na_matrices.Y[[as.character(p_list[i])]] = na_matrix.Y
    Y.completes[[as.character(p_list[i])]] = Y.complete
  }
  #print(svds[["0.1"]]$d)
  pcs = rep(1:max_pcs,length(p_list)) #repeats 1:pcs length of p_list times
  error.train = c()
  error.test = c()
  #length pcs = max_pcs*length(p_list)
  #This loop calculates the errors for each p and pcs combo
  for (i in 1:length(pcs)){
    output = xval.fast(X,
                       Y,
                       Y.completes[[as.character(vars[i])]],
                       svds[[as.character(vars[i])]],
                       na_matrices.X[[as.character(vars[i])]],
                       na_matrices.Y[[as.character(vars[i])]],
                       vars[i],
                       pcs[i],
                       error_metric = error_metric)
    error.train = c(error.train,output[[1]])
    error.test = c(error.test,output[[2]])
  }
  
  results = data.frame("p" = vars,
                       "pcs" = pcs,
                       "error.train" = error.train,
                       "error.test" = error.test)
  
  return(results)
}


xval_to_df.variance = function(var_list,max_pcs,data,error_metric = "NMSE_1"){
  #Uses xval.dropout.fast to make a data frame with errors for different assignments of var
  #very similar to xval_to_df.dropout accept adds noise instead of dropping values
  #
  #Inputs
  #var_list: list of var values to use when adding noise
  #max_pcs: maximum number of principal components to try when completing X
  #data: matrix that will have noise added
  #error_metric: one of NMSE_1, NMSE_2, R^2
  #
  #Output: data frame with 4 columns
  #p: variance of noise added
  #pcs: principal components used to approximate X
  #error.train: error on training data
  #error.test: error on test data
  data = data[1:(2*round(dim(data)[1]/2)),]
  train.idx = sample.int(dim(data)[1],round(0.5*dim(data)[1]))
  X = data[train.idx,]
  Y = data[-train.idx,]
  print(dim(X))
  print(dim(Y))
  #print("X")
  #print(X)
  
  #print("Y")
  #print(Y)
  
  
  #X = matrix(c(rep(-1,100),rep(1,100),rep(1,100)),100)
  #variance_list = c(0.2,0.4,0.6,0.8,1)
  #max_zeros = 3
  vars = c()
  svds = list()
  Y.completes = list()
  na_matrices.X = list() #na_matrices are filled with 0s. These are irrelevant arguments but
  na_matrices.Y = list() #xval.dropout.fast needs them to be passed through
  for (i in 1:length(var_list)){
    vars = append(vars,rep(var_list[i],max_pcs))
    data.noisy = add_noise(data,var_list[i])
    X.noisy = data.noisy[train.idx,]
    Y.noisy = data.noisy[-train.idx,]
    
    #print("new.X")
    #print(new.X)
    
    #print("new.Y")
    #print(new.Y)
    decomp = svd(X.noisy)
    
    #print("Y.comlpete")
    #print(scaling)
    #print(Y.complete)
    svds[[as.character(var_list[i])]] = decomp
    na_matrices.X[[as.character(var_list[i])]] = 0
    na_matrices.Y[[as.character(var_list[i])]] = 0
    Y.completes[[as.character(var_list[i])]] = Y.noisy
  }
  #print(svds[["0.1"]]$d)
  pcs = rep(1:max_pcs,length(var_list))
  error.train = c()
  error.test = c()
  for (i in 1:length(pcs)){
    output = xval.fast(X,
                       Y,
                       Y.completes[[as.character(vars[i])]],
                       svds[[as.character(vars[i])]],
                       na_matrices.X[[as.character(vars[i])]],
                       na_matrices.Y[[as.character(vars[i])]],
                       vars[i],
                       pcs[i],
                       error_metric = error_metric)
    error.train = c(error.train,output[[1]])
    error.test = c(error.test,output[[2]])
  }
  
  results = data.frame("p" = vars,
                       "pcs" = pcs,
                       "error.train" = error.train,
                       "error.test" = error.test)
  
  return(results)
}
xval.noisy = xval_to_df.variance(c(0.5),275,X.vH)

graph_xval(xval.noisy,"Gausian Noise","NMSE_1")
test = add_noise(matrix(1:300,30)+matrix(300:1,30),1)

xval.test2 = xval_to_df.dropout(c(0.3),9,X.,error_metric = "R^2")

graph_xval(xval.test2)

MSE_to_df = function(variance_list,max_pcs,spec){
  #make a data frame with errors for different assignments of var
  #
  #Inputs
  #variance_list: list of var values to use when adding noise
  #max_pcs: maximum number of principal components to try when completing X
  #data: matrix that will have noise added
  #error_metric: one of NMSE_1, NMSE_2, R^2
  #
  #Output: data frame with 3 columns
  #var: variance of noise added
  #pcs: principal components used to approximate X
  #MSE: error
X = get_data(df,spec,0)[[3]]
#X = matrix(c(rep(-1,100),rep(1,100),rep(1,100)),100)
#variance_list = c(0.2,0.4,0.6,0.8,1)
#max_zeros = 3
vars = c()
for (i in 1:length(variance_list)){
  vars = append(vars,rep(variance_list[i],max_pcs))
}

results = data.frame("var" = vars,
                     "pcs" = rep(1:max_pcs,length(variance_list)))


results = results%>%
  rowwise()%>%
  mutate("MSE" = MSE(X,var,pcs))

return(results)
}




MSE_to_df.dropout = function(p_list,max_pcs,data,error_metric = "NMSE_1",na_only = FALSE){
  #make a data frame with errors for different assignments of p
  #
  #Inputs
  #p_list: list of var values to use when dropping values out
  #max_pcs: maximum number of principal components to try when completing X
  #data: matrix that will have noise added
  #error_metric: one of NMSE_1, NMSE_2, R^2
  #na_only: return results for only NA values that dropped out of X when noise was added
  #
  #Output: data frame with 3 columns
  #var: variance of noise added
  #pcs: principal components used to approximate X
  #MSE: error
  X = data
  #X = matrix(c(rep(-1,100),rep(1,100),rep(1,100)),100)
  #variance_list = c(0.2,0.4,0.6,0.8,1)
  #max_zeros = 3
  vars = c()
  svds = list()
  na_matrices = list()
  for (i in 1:length(p_list)){
    vars = append(vars,rep(p_list[i],max_pcs))
    new = random_dropout(X,p_list[i])
    na_matrix = new[[2]]
    Z = fill.na(new[[1]])
    decomp = svd(Z)
    svds[[as.character(p_list[i])]] = decomp
    na_matrices[[as.character(p_list[i])]] = na_matrix
  }
  #print(svds[["0.1"]]$d)
  
  results = data.frame("p" = vars,
                       "pcs" = rep(1:max_pcs,length(p_list)))
  
  
  results = results%>%
    rowwise()%>%
    mutate("error" = MSE.dropout.fast(X,
                                    svds[[as.character(p)]],
                                    na_matrices[[as.character(p)]],
                                    p,
                                    pcs,
                                    error_metric = error_metric,
                                    na_only = na_only))
  
  return(results)
}

graph_MSE = function(results,title){
  #graphs a data frame from MSE_to_df
results%>%
  ggplot(aes(x = pcs,y = MSE,color = var))+
  geom_point(aes(color = var))+
  geom_line(aes(group = var))+
  ggtitle(title)
  
}

graph_MSE.dropout = function(results,title,error_metric){
  #graphs a data from from MSE_to_df.dropout
  results%>%
    ggplot(aes(x = pcs,y = error,color = p))+
    geom_point(aes(color = p))+
    geom_line(aes(group = p))+
    ggtitle(title)+
    ylab(error_metric)+
    xlab("Principal Components")
  
}


graph_xval = function(results,title,error_metric){
  #graphs a data frame from xval_to_df.dropout or xval_to_df.variance
  results.long = results%>%
    pivot_longer(c(error.train,error.test))
  
  results.long %>%
    ggplot(aes(x = pcs, y = value, color = name))+
    geom_point()+
    geom_line(aes(group = name))+
    ggtitle(title)+
    ylab(error_metric)+
    xlab("Principal Components")
}

xval.lowp = xval_to_df.dropout(c(0.01),275,X.vH,error_metric = "R^2")
graph_xval(xval.lowp,"dropout p=0.01","R^2")

xval.zero = xval_to_df.dropout(c(0),275,X.vH,error_metric = "R^2")
graph_xval(xval.zero,"dropout p=0","R^2")


#MSE with noise
results1 = MSE_to_df(c(0.1,0.2,0.3,0.5,0.8,1),9,1)
results2 = MSE_to_df(c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),9,1)
results3 = MSE_to_df(c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),9,2)
graph_MSE(results1,"X.L")
graph_MSE(results2,"X.L")
graph_MSE(results3,"X.H")


#MSE with dropout
results4 = MSE_to_df.dropout(c(0.2,0.4,0.6,0.8),9,1)
graph_MSE(results4,"dropout")
results4 = MSE_to_df.dropout(c(0.1,0.2,0.3,0.4,0.5),3,1)
graph_MSE(results4,"dropout")
results5 = MSE_to_df.dropout(c(0.2,0.4,0.6,0.8),9,3,na_only = TRUE)
graph_MSE(results5,"dropout na only")

results5 = MSE_to_df.dropout(c(0.5),20,3,na_only = TRUE)
graph_MSE(results5,"dropout na only")

results5 = MSE_to_df.dropout(p_list = c(0.1),
                             max_pcs = 275,
                             data  = X.vH,
                             error_metric = "R^2",
                             na_only = TRUE)
graph_MSE.dropout(results5,"dropout na only","R^2")

full_pc.R2 = MSE_to_df.dropout(p_list = c(0.1),
                             max_pcs = 275,
                             data  = X.vH,
                             error_metric = "R^2",
                             na_only = TRUE)



full_pc.full_matrix = MSE_to_df.dropout(p_list = c(0.1),
                             max_pcs = 275,
                             data  = X.vH,
                             error_metric = "NMSE_2",
                             na_only = FALSE)

full_pc.multi_prob = MSE_to_df.dropout(p_list = c(0.1,0.2,0.3),
                                   max_pcs = 275,
                                   data  = X.vH,
                                   error_metric = "NMSE_2",
                                   na_only = TRUE)

graph_MSE.dropout(full_pc.R2,"dropout na only R^2","R^2")
graph_MSE.dropout(full_pc.multi_prob,"dropout na only","NMSE_2")
graph_MSE.dropout(full_pc.full_matrix,"dropout full matrix","NMSE_2")

ptm = proc.time()
MSE_to_df.dropout(c(0.1,0.2,0.3),5,X.vH)
print(proc.time()-ptm)

new = random_dropout(X.vH,0.1)
na_matrix = new[[2]]
Z = fill.na(new[[1]])
decomp = svd(Z)

ptm = proc.time()
MSE.dropout.fast(X.vH,
                 decomp,
                 na_matrix,
                 0.1,
                 3,
                 na_only = FALSE)
print(proc.time()-ptm)


ptm = proc.time()
MSE.dropout(X.vH,0.1,3,na_only = FALSE)
print(proc.time()-ptm)

XY = data.frame(X = c(1,2,3),Y = c(1,2,3))

ggplot(XY,aes(x = X, y = Y))+geom_point()

#write.csv(full_pc.full_matrix,"full_pc_full_matrix.csv")
#write.csv(full_pc.multi_prob,"full_pc_multi_prob.csv")
#write.csv(full_pc.R2,"full_pc_R2.csv")
write.csv(xval.test,"xval_test.csv")


test = matrix(1:9,3)

svd = svd(test)
U = svd$u
D = svd$d
V = svd$v

D_trunc = diag(D)
D_trunc[,-(1:2)] = 0

U_trunc = U[,1:2]

U
U_trunc
D_trunc

U%*%D_trunc%*%t(V)

U_trunc%*%t(U_trunc)%*%test

xval.test.long = xval.test%>%
  pivot_longer(c(error.train,error.test))

xval.test.long %>%
  ggplot(aes(x = pcs, y = value, color = name))+
  geom_point()+
  geom_line(aes(group = name))







 


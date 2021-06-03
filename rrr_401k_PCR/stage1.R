#stage 1 function
l=0.1


RMD_dantzig <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  
  p <- ncol(G)
  zp <- rep(0, p)
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,T,X,TX)
  
  A <- solve(diag(D),G)
  R <- rbind(A, -A)
  
  a <- solve(diag(D),M)
  r <- c(a - lambda*L, -a - lambda*L)
  
  if(sparse) {
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  } else {
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  
  return(f)
}

RMD_lasso <- function(M, G, D, lambda=0, control = list(maxIter = 1000, optTol = 10^(-5), 
                                                        zeroThreshold = 10^(-6)), beta.start = NULL) {
  

  
  p <- ncol(G)
  
  Gt<-G
  Mt<-M
  
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,...)
  lambda_vec=lambda*L*D #v3: insert D here
  
  if (is.null(beta.start)) {
    beta <- rep(0,p) #vs low-dimensional initialization
  }
  else {
    beta <- beta.start
  }
  wp <- beta
  mm <- 1
  while (mm < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      rho=Mt[j]-Gt[j,]%*%beta+Gt[j,j]*beta[j]
      z=Gt[j,j]
      
      if (sum(is.na(rho)) >= 1) {
        beta[j] <- 0
        next
      }
      if (rho < -1 * lambda_vec[j]) 
        beta[j] <- (rho+lambda_vec[j])/z
      if (abs(rho) <= lambda_vec[j]) 
        beta[j] <- 0
      if (rho > lambda_vec[j]) 
        beta[j] <- (rho-lambda_vec[j])/z
    }
    wp <- cbind(wp, beta)
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    mm <- mm + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = mm))
}


inv.D <- function(D,pcs){
  #Used for the psuedoinverse
  #Returns 1/D for first pcs on diagonal and 0 for all past pcs
  #if pcs is zero, sets all pcs that are too small to zero (can change threshold in line 87)
  if (pcs == 0){
    D.inv = 1/D
    for (i in 1:length(D.inv)){
      if (D.inv[i]>100){
        D.inv[i] = 0
      } 
    }
    return(diag(D.inv))
  } else {
  D.inv = 1/D
  D.inv[-(1:pcs)] = 0
  return(diag(D.inv))
  }
}

PCR<-function(Y,T,X,b,pcs,is_alpha){
  #Do the SVD of X and keep first pcs principal components
  decomp = svd(X)
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  #Combine new X with T into a matrix W
  W_trunc = cbind(T,U%*%D_trunc%*%t(V))
  
  n=nrow(X)
  p=length(b(T[1],X[1,]))
  
  x.cols = length(X[1,])+1
  
  B=matrix(0,n,p)
  
  #Apply the dictionary b to W
  for (i in 1:n){
    B[i,]=b(W_trunc[i,1],W_trunc[i,2:x.cols])
  }
  
  #Do the SVD of B
  decomp.B = svd(B)
  U.B = decomp.B$u
  D.B = decomp.B$d
  V.B = decomp.B$v
  
  if (is_alpha ==1){
    #calculate G^-1 and M using the psuedoinverse
    G = V.B%*%inv.D(D.B,0)^2%*%t(V.B)
    M = rep(0,p)
    for (i in 1:n){
      M = M + (b(1,W_trunc[i,2:x.cols])-b(0,W_trunc[i,2:x.cols]))
    }
    
    
    return(G%*%M)
  }
  
  if (is_alpha==0){
    #calculate beta using the psuedoinverse
    print(paste("pcs:",pcs))
    beta = V.B%*%inv.D(D.B,0)%*%t(U.B)%*%Y
    
    return(beta)
  }
}


PCR.missing<-function(Y,T,X,b,pcs,is_alpha){
  
  #fill in NA values and find delta
  N = dim(X)[1]
  M = dim(X)[2]
  scale.ratios = 1-(t(is.na(X))%*%matrix(1,N,1))/N
  delta = diag(c(1/scale.ratios))
  X = X/t(matrix(scale.ratios,M,N))
  X[is.na(X)] = 0
  #Do the SVD of X and keep first pcs principal components
  decomp = svd(X)
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  #Combine new X with T into a matrix W
  W_trunc = cbind(T,U%*%D_trunc%*%t(V))
  
  n=nrow(X)
  p=length(b(T[1],X[1,]))
  
  x.cols = length(X[1,])+1
  
  #Apply the dictionary b to W
  B=matrix(0,n,p)
  for (i in 1:n){
    B[i,]=b(W_trunc[i,1],W_trunc[i,2:x.cols])
  }
  
  decomp.B = svd(B)
  U.B = decomp.B$u
  D.B = decomp.B$d
  V.B = decomp.B$v
  
  if (is_alpha ==1){
    G = V.B%*%inv.D(D.B,0)^2%*%t(V.B)
    M = rep(0,p)
    for (i in 1:n){
      M = M + (b(1,W_trunc[i,2:x.cols])-b(0,W_trunc[i,2:x.cols]))
    }
    
    
    return(list(G%*%M,delta))
  }
  
  if (is_alpha==0){
    print(paste("pcs:",pcs))
    beta = V.B%*%inv.D(D.B,0)%*%t(U.B)%*%Y
    
    return(list(beta,delta))
  }
}


PCR2<-function(Y,T,X,b,pcs,is_alpha){
  #Same as PCR but does SVD on T and X together. Performed poorly in trials
  decomp = svd(cbind(T,X))
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  W_trunc = U%*%D_trunc%*%t(V)
  
  n=nrow(X)
  p=length(b(T[1],X[1,]))
  
  x.cols = length(X[1,])+1
  
  B=matrix(0,n,p)
  for (i in 1:n){
    B[i,]=b(W_trunc[i,1],W_trunc[i,2:x.cols])
  }
  
  decomp.B = svd(B)
  U.B = decomp.B$u
  D.B = decomp.B$d
  V.B = decomp.B$v
  
  if (is_alpha ==1){
    G = V.B%*%inv.D(D.B,0)^2%*%t(V.B)
    M = rep(0,p)
    for (i in 1:n){
      M = M + (b(1,W_trunc[i,2:x.cols])-b(0,W_trunc[i,2:x.cols]))
    }
    
    return(G%*%M)
  }
  
  if (is_alpha==0){
    print(paste("pcs:",pcs))
    beta = V.B%*%inv.D(D.B,0)%*%t(U.B)%*%Y
    
    return(beta)
  }
}


PCR3<-function(Y,T,X,b,pcs,is_alpha){
  #Same as PCR2 but with no dictionary mapping. Also performed poorly
  decomp = svd(cbind(T,X))
  U = decomp$u
  D = decomp$d
  V = decomp$v
  D_trunc = diag(D)
  D_trunc[,-(1:pcs)] = 0
  
  W_trunc = U%*%D_trunc%*%t(V)
  
  n=nrow(X)
  p=length(c(T[1],X[1,]))
  
  x.cols = length(X[1,])+1
  
  if (is_alpha ==1){
    G = V%*%inv.D(D,pcs)^2%*%t(V)
    M = rep(0,p)
    for (i in 1:n){
      M = M + (c(1,W_trunc[i,2:x.cols])-c(0,W_trunc[i,2:x.cols]))
    }
    
    return(G%*%M)
  }
  
  if (is_alpha==0){
    print(paste("pcs:",pcs))
    beta = V%*%inv.D(D,pcs)%*%t(U)%*%Y
    
    return(beta)
  }
}

get_D <- function(Y,T,X,m,rho_hat,b){
  n=nrow(X)
  p=length(b(T[1],X[1,]))
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=b(T[i],X[i,])*as.vector(rho_hat %*% b(T[i],X[i,]))-m(Y[i],T[i],X[i,],b)
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

c=0.5
alpha=0.1
tol=1e-6

RMD_stable<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,is_alpha,is_lasso){
  
  k=1
  
  p=length(b(T[1],X[1,]))
  n=length(T)
  
  # low-dimensional moments
  X0=X[,1:p0]
  MNG0<-get_MNG(Y,T,X0,b)
  M_hat0=MNG0[[1]]
  N_hat0=MNG0[[2]]
  G_hat0=MNG0[[3]]
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0=solve(G_hat0,N_hat0)
  beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
  
  # moments
  MNG<-get_MNG(Y,T,X,b)
  M_hat=MNG[[1]]
  N_hat=MNG[[2]]
  G_hat=MNG[[3]]
  
  # penalty
  lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) # snippet
  
  if(is_alpha){ 
    ###########
    # alpha_hat
    ###########
    diff_rho=1
    while(diff_rho>tol & k<=max_iter){
      
      # previous values
      rho_hat_old=rho_hat+0
      
      # normalization
      D_hat_rho=get_D(Y,T,X,m,rho_hat_old,b)
      D_hat_rho=pmax(D_LB,D_hat_rho)
      D_hat_rho=D_hat_rho+D_add
      
      # RMD estimate
      if(is_lasso){
        rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }else{
        rho_hat=RMD_dantzig(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }
      
      # difference
      diff_rho=two.norm(rho_hat-rho_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(rho_hat)
    
  } else { 
    ###########
    # gamma_hat
    ###########
    diff_beta=1
    while(diff_beta>tol & k<=max_iter){
      
      # previous values
      beta_hat_old=beta_hat+0
      
      # normalization
      D_hat_beta=get_D(Y,T,X,m2,beta_hat_old,b)
      D_hat_beta=pmax(D_LB,D_hat_beta)
      D_hat_beta=D_hat_beta+D_add
      
      # RMD estimate
      if(is_lasso){
        beta_hat=RMD_lasso(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }else{
        beta_hat=RMD_dantzig(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }
      
      # difference
      diff_beta=two.norm(beta_hat-beta_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(beta_hat)
    
  }
}

arg_Forest<- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
arg_Nnet<- list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)

get_stage1<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,alpha_estimator,gamma_estimator,pcs = 10){
  
  p=length(b(T[1],X[1,]))
  n=length(T)
  MNG<-get_MNG(Y,T,X,b)
  B=MNG[[4]]
  
  ###########
  # alpha hat
  ###########
  if(alpha_estimator==0){ # dantzig
    
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,0)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  } else if(alpha_estimator==1){ # lasso
    
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,1)
    
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  } else if(alpha_estimator==2){ # PCR
    
    rho_hat=PCR(Y,T,X,b,pcs,1)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  } else if(alpha_estimator==3){ # PCR2
    
    rho_hat=PCR2(Y,T,X,b,pcs,1)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  }else if(alpha_estimator==4){ # PCR3 no dict
    
    rho_hat=PCR3(Y,T,X,b,pcs,1)
    alpha_hat<-function(d,z){
      return(c(d,z)%*%rho_hat)
    }
    
  }else if(alpha_estimator==5){ # PCR.missing 
      
    res=PCR.missing(Y,T,X,b,pcs,1)
    
    rho_hat = res[[1]]
    delta_hat = res[[2]]
    alpha_hat<-function(d,z){
      z[is.na(z)] = 0
      z = z%*%delta_hat
      return(b(d,z)%*%rho_hat)
    }
    
  }
  
  ###########
  # gamma hat
  ###########
  if(gamma_estimator==0){ # dantzig
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,0)
    gamma_hat<-function(d,z){
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==1){ # lasso
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,1)
    gamma_hat<-function(d,z){ 
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==2){ # random forest
    
    forest<- do.call(randomForest, append(list(x=B,y=Y), arg_Forest))
    gamma_hat<-function(d,z){
      return(predict(forest,newdata=b(d,z), type="response"))
    }
    
  } else if(gamma_estimator==3){ # neural net
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    
    maxs_Y<-max(Y)
    mins_Y<-min(Y)
    
    # hack to ensure that constant covariates do not become NA in the scaling
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])
    
    NN_Y<-scale(Y, center = mins_Y, scale = maxs_Y - mins_Y)
    
    nn<- do.call(nnet, append(list(x=NN_B,y=NN_Y), arg_Nnet))
    gamma_hat<-function(d,z){
      
      test<-t(as.vector(b2(d,z)))
      NN_b<-test
      NN_b[,keep]<-scale(t(NN_b[,keep]), 
                         center = mins_B[keep], 
                         scale = maxs_B[keep] - mins_B[keep])
      
      NN_Y_hat<-predict(nn,newdata=NN_b)
      Y_hat=NN_Y_hat*(maxs_Y-mins_Y)+mins_Y
      
      return(Y_hat)
    }} else if (gamma_estimator == 4){ #PCR
      beta = PCR(Y,T,X,b,pc = pcs,0)
      gamma_hat <- function(d,z){
        return(b(d,z)%*%beta)
      }
    }  else if (gamma_estimator == 5){ #PCR2
      beta = PCR2(Y,T,X,b,pc = pcs,0)
      gamma_hat <- function(d,z){
        return(b(d,z)%*%beta)
      }
    } else if (gamma_estimator == 6){ #PCR3 no dict
      beta = PCR3(Y,T,X,b,pc = pcs,0)
      gamma_hat <- function(d,z){
        return(c(d,z)%*%beta)
      }
    } else if (gamma_estimator == 7){ #PCR.missing
      res = PCR.missing(Y,T,X,b,pc = pcs,0)
      
      beta = res[[1]]
      delta_hat = res[[2]]
      gamma_hat <- function(d,z){
        z[is.na(z)] = 0
        z = z%*%delta_hat
        
        return(b(d,z)%*%beta)
      }
    }
  
  return(list(alpha_hat,gamma_hat))
  
}

get_data<-function(df,spec,quintile){
  
  Y <- df[,"net_tfa"]
  T <- df[,"e401"]
  
  ## low dim specification
  X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
  
  ## high dim specification. NOTE: original paper is this spec squared (pairwise interactions thereof?)
  X.H <- cbind(poly(df[,"age"], 6, raw=TRUE),
               poly(df[,"inc"], 8, raw=TRUE),
               poly(df[,"educ"], 4, raw=TRUE),
               poly(df[,"fsize"], 2, raw=TRUE),
               df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"]) 
  
  # copied from EJ: "(poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown)^2"
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
  
  if (spec==1){
    X=X.L
  } else if (spec==2) {
    X=X.H
  } else {
    X=X.vH
  }
  
  X <- scale(X,center=TRUE,scale=TRUE)
  
  #impose common support
  p.1 <- multinom(T~X-1, trace=FALSE)$fitted.values
  indexes.to.drop <- which(p.1 < min(p.1[T==1]) | max(p.1[T==1]) < p.1)
  if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
  n.per.treatment <- as.vector(table(T[-indexes.to.drop]))
  n.trim <- n.per.treatment[1]+n.per.treatment[2]
  
  Y.trimmed=Y[-indexes.to.drop]
  T.trimmed=T[-indexes.to.drop]
  X.trimmed=X[-indexes.to.drop,]
  
  if (spec==1){
    inc=X.trimmed[,2]
  } else if (spec==2) {
    inc=X.trimmed[,7]
  } else {
    inc=X.trimmed[,7]
  }
  
  if (quintile>0){
    q <- ntile(inc, 5)
    Y.q=Y.trimmed[q==quintile]
    T.q=T.trimmed[q==quintile]
    X.q=X.trimmed[q==quintile,]
  } else {
    Y.q=Y.trimmed
    T.q=T.trimmed
    X.q=X.trimmed
  }
  
  return(list(Y.q,T.q,X.q))
  
}


simulate_data = function(n,method = 3,rank = 5){
  ###designed to return list(Y,T,X) with ATE 2.2
  ###Inputs
  #n: dimensions of data Y:length n vec T:length n vec, X: n x n
  #method: method for generating X. 
  #  1 Method specified by Rahul
  #  2 Adjustments to Rahul's method to have faster decay in singular values
  #  3 Low rank meht
  #rank: indicates rank of X when method = 3
  ###Output: list(Y,T,X)
  
  X = matrix(0,n,100)
  Y = c()
  T = c()
  beta = c()
  for (i in 1:n){
    beta = c(beta,1/(i^2))
  }
  if (method ==1){
    sigma = diag(100)
    for (i in 1:99){
      sigma[i,i+1] = 0.5
      sigma[i+1,i] = 0.5
    }
  } else if (method == 2){
    #can adjust power of abs(i-j) to change decay of singular values (smaller = faster decay)
    sigma = matrix(0,100,100)
    for (i in 1:100){
      for (j in 1:100){
        sigma[i,j] = 1/(abs(i-j)^0.1+1)
        sigma[j,i] = 1/(abs(i-j)^0.1+1)
      }
    }
  } else if (method ==3){
    #create two matrices of random gaussians
    U = matrix(rnorm(n*rank),n,rank)
    V = matrix(rnorm(n*rank),n,rank)
    X = U%*%t(V)
    v = rnorm(n,0,1)
    eps = rnorm(n,0,1)
    T = as.vector(ifelse(3*X%*%beta+0.75*v>=0,1,0))
    Y = 2.2*T+X%*%beta+X[,1]*T+eps
    return(list(Y,T,X))
  }
  for (i in 1:n){
    X_i = rmvnorm(1,rep(0,100),sigma)
    v = rnorm(1,0,1)
    eps = rnorm(1,0,1)
    if (3*X_i%*%beta+0.75*v>=0){
      T_i = 1
    }else{
      T_i = 0
    }
    Y = c(Y,1.2*T_i+1.2*X_i%*%beta+T_i^2+T_i*X_i[1]+eps)
    T = c(T,T_i)
    X[i,] = X_i
  }
  
  return(list(Y,T,X))
  
  
  
}






  
  
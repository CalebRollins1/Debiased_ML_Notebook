#alpha_avg = c(0,0)

two.norm <- function(x){
  return(sqrt(x %*% x))
} 

one.norm<-function(x){
  return(sum(x%*%sign(x)))
}

one.norm.grad<-function(x){
  return(sign(x))
}

# intercept
b<-function(d,z){
  return(c(1,d,z))
}

b.PCR <- function(d,z){
  return(c(d,z))
}

# intercept and interaction
b2<-function(d,z){
  return(c(1,d,z,d*z))
}

obj<-function(rho,G,M,r){
  return(t(rho)%*%G%*%rho-2*t(M)%*%rho+2*r*one.norm(rho))
}

m<-function(y,d,z,gamma){ #all data arguments to make interchangeable with m2
  return(gamma(1,z)-gamma(0,z))
}

m2<-function(y,d,z,gamma){
  return(y*gamma(d,z))
}

psi_tilde<-function(y,d,z,m,alpha,gamma){
  
  #alpha_avg <<- c((alpha_avg[2]*alpha_avg[1]+alpha(d,z))/(alpha_avg[2]+1),alpha_avg[2]+1)
  return(m(y,d,z,gamma)+alpha(d,z)*(y-gamma(d,z)))
}

psi_tilde_bias<-function(y,d,z,m,alpha,gamma){
  return(m(y,d,z,gamma))
}

get_MNG<-function(Y,T,X,b){
  
  p=length(b(T[1],X[1,]))
  n.nl=length(T)
  
  B=matrix(0,n.nl,p)
  M=matrix(0,p,n.nl)
  N=matrix(0,p,n.nl)
  
  for (i in 1:n.nl){
    B[i,]=b(T[i],X[i,])
    M[,i]=m(Y[i],T[i],X[i,],b)
    N[,i]=m2(Y[i],T[i],X[i,],b)  # this is a more general formulation for N
  }
  
  M_hat=rowMeans(M)
  N_hat=rowMeans(N)
  G_hat=t(B)%*%B/n.nl
  
  return(list(M_hat,N_hat,G_hat,B))
}

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
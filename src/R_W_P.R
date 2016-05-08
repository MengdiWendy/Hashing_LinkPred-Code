# 6-18-2015
# code for the publish version of "Non-transitive Hashing with Latent Similarity Components".
# update June 18,2015
# author: Mengdi Wang
# this is the regularization part for objective function

R_W_P <- function(W_matrix,X,H_matrix,r1,r2){
  N = ncol(X)
  L = nrow(X)
  Km = ncol(W_matrix)
  WW = diag(tcrossprod(W_matrix, W_matrix)) # product of W_i and t(W_i)
  
  R_W = mat.or.vec(L, Km)
  
  for(k in 1:Km){
    R_W[,k]= r1*4*(sum(WW)-WW[k])*W_matrix[,k]+r2*2*W_matrix[,k]     
  }
  
  return (R_W)
}






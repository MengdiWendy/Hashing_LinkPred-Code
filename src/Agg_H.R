# July 4th, 2015
# update July 4th, 2015
# @author: Mengdi Wang

# H_matrix - N*Km matrix, N: number of nodes, Km= K*M
# K - number of bits in each hashing table
# M - number of hashing tables
library(e1071)
Agg_H <- function(H_matrix, K, M){
  Km = ncol(H_matrix)
  N = nrow(H_matrix)
 
  H = data.frame(H_matrix)
  start <- seq(1, by = K, length = ncol(H) / K)
  K = K-1
  H = lapply(start, function(i, H) H[i:(i+K)], H = H)
 
  aH = mat.or.vec(N,N)
 
  for(i in 1:M){
    H_m = as.matrix(H[[i]])
    temp = hamming.distance(H_m)
    temp = (K+1)^(i-1)*(temp+1)
    temp = 1/temp
    aH = aH + temp
  }
 
 aH = 1/aH
 
 return (aH)
}
 
 

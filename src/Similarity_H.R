# 
# Date: Dec.06, 2015
# @author: Mengdi Wang

# H_matrix - N*Km matrix, N: number of nodes, Km= K*M
# K - number of bits in each hashing table
# M - number of hashing tables

Similarity_H <- function(H_matrix, K, M){
  Km = ncol(H_matrix)
  N = nrow(H_matrix)
  
  H = data.frame(t(H_matrix))
  start <- seq(1, by = K, length = ncol(H) / K)
  K = K-1
  H = lapply(start, function(i, H) H[i:(i+K)], H = H)
  
  upper = lapply(H, function(x) cal_upper(x))
  upper = Reduce("+", upper)
  down = lapply(H, function(x) cal_down(x))
  down = Reduce("+", down)
  
  return (upper/down)
}

cal_upper <- function(H){
  x = as.matrix(H)
  S = cosine(t(x))
  return (S*exp(S))
}
cal_down <- function(H){
  x = as.matrix(H)
  S = cosine(t(x))
  return (exp(S))
}

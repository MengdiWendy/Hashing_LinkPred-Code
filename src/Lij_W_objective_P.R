# June 27, 2015 
# objective function in published paper (without regularization)
# Mij : objective function for similairty pairs
# H is a list [H_1', H_2',..., H_m'], where H_1' is [N,K], H_m is [K,N]
# H_matrix is Km*N (Km is the whole hashing code length)
# W is a list [W_1, W_2, ..., W_m], where W_1 [L,K]
# X is a matrix [L,N]
# R is -1 when disimilarity; 1 for similairity 
# cor is the location of the non-zero value in the relationship matrix. (-1 means disimilarity, 1 means similarity)
# return value is equation (16), which is [W_1, W_2,..., W_m], with size [L, KM] 

Lij_W_objective_P <- function(cor, X, W, H, H_matrix, r, mFirstTerm,WX_list){
  
  L = nrow(W[[1]])
  K = ncol(W[[1]])
  i = cor[1]
  j = cor[2]
  num_Tables = length(W)
  # eSij_sum is the sum of e(S^n_ij) (n = 1,...,M) in equation (17) in that paper
  

  # get W_m and concatenate them together to get W_M = Lij_W_ob 
  # W_M is [W_1, W_2,..., W_m], with size [L, KM]  
  Lij_W_ob = NULL
  for (m in 1:num_Tables){
    Lij_W_ob = cbind(Lij_W_ob, Sij_Wk_ob(i,j,X, W[[m]], t(H[[m]]), H_matrix, mFirstTerm[[m]],WX_list[[m]]))
  } 
  # t is the first term of the right part of equation (16)  
  # Lij_W_ob is the second term of the right part of equation (16)  
  t = -r/(1+exp(r*crossprod(H_matrix[,i], H_matrix[,j])))
  t = as.numeric(t)
  Lij_W_ob = t*Lij_W_ob
  return (Lij_W_ob)  
}

# this function is for equation (17)
Sij_Wk_ob <- function(i, j,X, Wm,  Hm, H_matrix, FirstTerm, WX){
  L = nrow(X)
  K = ncol(Wm)  
  Sij_Wk_o = mat.or.vec(L, K)  
  #  concatenate K [L 1] vector to get the equation below (17)  
  for(k in 1:K){       
    Sij_Wk_o[,k] = FirstTerm[i,j]*(Hm[k,j]*WX[k,i]*X[,i]+Hm[k,i]*WX[k,j]*X[,j])
  }  
  return (Sij_Wk_o)
}


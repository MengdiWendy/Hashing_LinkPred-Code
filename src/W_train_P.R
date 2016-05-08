# 6-28-2015
# update June 28, 2015
# update May 2nd, 2016
# @author: Mengdi Wang
## this is the main function to train W.

source("R_W_P.R")
source("Lij_W_objective_P.R")
source("Lij_W_objective_P_nonbinary.R")


# This function return the list(W_matrix=W_matrix, converge = 1, steps = i )

W_train_P <- function(W_matrix, K, X, similarity, disimilarity, r1, r2, beta, a, alpha, M, N){
  i=1
  thro =10


  while(thro>0.001){
    print(i)
    #print(W_matrix)
    
   
    W_new_matrix = W_matrix-alpha*L_W_P(W_matrix, K, X, similarity, disimilarity, r1, r2, beta, a, M, N)
   
    i = i+1
    
    Wm_new= cal_2_norm(W_new_matrix)    
    cat("2 norm of W_new_matrix: ",Wm_new)
    Wm = cal_2_norm(W_matrix)
    thro = abs((Wm_new-Wm))
#    inter = c(inter, thro)
    cat("The threshold is:",thro) 
    cat("    ")
    W_matrix = W_new_matrix 
####
## if step size is larger than 50, stop
## if step size is larger than 50 and thro > 0.005, set converge = 0
## if step size is larger than 50 and thro < 0.005, set converge = 2, meaning that it close to converge (which is 0.001)
## if thro < 0.001 within 50 steps, set converge = 1.
    if(i>70){
      if(thro>0.005){return (list(W_matrix=W_matrix, converge = 0, steps=i ))
      }else{return (list(W_matrix=W_matrix, converge = 2, steps = i ))
      }
    }
  }
  cat('\n')
  return (list(W_matrix=W_matrix, converge = 1, steps = i ))
}

# calucate 2 norm, because norm(X, type = c("2")) does not work in this R version
cal_2_norm <- function(W){
  s = svd(W)
  d = s$d
  norm = d[1]
  return (norm)
}

# L_W_P calculate one step for W 
L_W_P <- function(W_matrix, K, XF, similarity, disimilarity, r1, r2, beta, a, M, N ){
  
  L = nrow(XF)
  Km = ncol(W_matrix)
  
  L_W = matrix(0,L, Km)
  M_W = matrix(0,L, Km)
  
  H_matrix = 2/(1+exp(-crossprod(W_matrix,XF)))-1
  
  W =data.frame( W_matrix)
  H = data.frame(t(H_matrix))
  
  start <- seq(1, by = K, length = ncol(W) / K)
  K = K-1
  W <- lapply(start, function(i, W) W[i:(i+K)], W = W)
  H <- lapply(start, function(i, H) H[i:(i+K)], H = H)
  
  ############################
  # pre-calculate 
  # calculate mTristTerm for latter use 
  num_M = length(W)
  allSimilarity = crossprod(H_matrix, H_matrix)  
  
  mSimilarity = lapply(H, function(x) tcrossprod(as.matrix(x),as.matrix(x)))
  
  expSimilarity = lapply(mSimilarity, function(x) exp(x))  
  eSij_sum = Reduce('+',expSimilarity )
    
  mFirstTerm = lapply(mSimilarity, function(x, allSimilarity, eSij_sum) (1+x-allSimilarity)*exp(x)/eSij_sum,allSimilarity=allSimilarity,eSij_sum=eSij_sum )
  
  
  ##############
  # pre-calculate WX_list
  
  WX_list = lapply(W, function(x, XFeature) 2*exp(-crossprod(as.matrix(x),XFeature))/((1+exp(-crossprod(as.matrix(x),XFeature)))^2), XFeature = XF)

  library(parallel)
  
  source("Lij_W_objective_P.R")
  test <- Lij_W_objective_P(c(similarity[1,1], similarity[1,2]), XF, W, H,H_matrix, 1, mFirstTerm, WX_list)
  
  if(N > 0){
#    M_W <- lapply(similarity, function(x, XF, W, H,H_matrix, r,mFirstTerm,WX_list) 
#               Lij_W_objective_P(x, XF, W, H,H_matrix, r,mFirstTerm, WX_list), XF = XF, W = W, H = H,H_matrix= H_matrix, r = 1,mFirstTerm = mFirstTerm, WX_list = WX_list) 
    M_W <- mclapply(similarity, function(x, XF, W, H,H_matrix, r,mFirstTerm,WX_list) 
                     Lij_W_objective_P(x, XF, W, H,H_matrix, r,mFirstTerm, WX_list), XF = XF, W = W, H = H,H_matrix= H_matrix, r = 1,
                    mFirstTerm = mFirstTerm, WX_list = WX_list,mc.cores=6)    
    M_W <- Reduce("+", M_W)
  }
  #  registerDoMC(2)
  #  if(N>0){
  #    M_W <- foreach(i=1:N,.combine='+')%dopar%Mij_W(similarity[i,1], similarity[i,2], K, X, W, H, a)
  #  }
  #  if(N > 0){
  #    for(i in 1:N){
  #      M_W = M_W + Mij_W(similarity[i,1], similarity[i,2], K, X, W, H, a)
  #    }
  #  }
  
  if(M > 0){
    L_W <- mclapply(disimilarity, function(x, XF, W, H,H_matrix, r,mFirstTerm,WX_list) 
      Lij_W_objective_P(x, XF, W, H,H_matrix, r,mFirstTerm, WX_list), XF = XF, W = W, H = H,H_matrix= H_matrix, r = -1,mFirstTerm = mFirstTerm, WX_list = WX_list,mc.cores=6)    
    L_W <- Reduce("+", L_W)
  }
  #  if(N>0){
  #   L_W <- foreach(i=1:N,.combine='+')%dopar%Lij_W(disimilarity[i,1],disimilarity[i,2], X, W, H)
  # }
  
  #  if(N > 0){
  #    for(i in 1:N){
  #      L_W = L_W + Lij_W(disimilarity[i,1],disimilarity[i,2], X, W, H)
  #    }
  #  }
  source("R_W_P.R")
  R_W_matrix = R_W_P (W_matrix,XF,H_matrix,r1,r2)
  
  L_W = M_W + beta*L_W+R_W_matrix
  
  return (L_W)
}


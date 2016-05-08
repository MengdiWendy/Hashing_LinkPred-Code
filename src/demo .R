
source("get.AUC.R")
source("Agg_H.R")
source("W_train_P.R")
source("get.select.node.ID.R") # optional 

# change the direction to the data folder
setwd("~/data")
# load testing and training data
load("AA.similarity.RData")
load("testing.1178.RData")
load("training.1178.RData")

testing = R.2011.testing.1178
training = R.2011.training.1178

# change to the src folder
setwd("~/src")
# set the parameters 
# here we use 2 hashing tables, in each hashing table we use 16bits
K = 16
r1 = 10^(-6)
r2 = 10^2
beta = 1
a = 1
alpha = 0.001

# select some nodes to focus, this is optional
# Here we use nodes with more than 2 degrees in both testing and training 
negative.positive.matrix <- get.negative.positive.matrix(testing, training)
positiveMatrix <- negative.positive.matrix$positiveMatrix
negativeMatrix <- negative.positive.matrix$negativeMatrix
select.node.ID <- get.select.node.ID(testing, training, positiveMatrix, degree=2)

# prepare feature matrix 
bits = 100 # Here we use 100 bits for AA similarity matrix
AA_feature <- svd(AA.similarity, nu = bits)$u
X = AA_feature # X is our feature matrix 
# Here we only use one feature to build our feature matrix in this demo, but we can use more feature matrix
# and the feature matrix could be X = cbind(as.matrix(feature[1], feature[2],..., feature[N])


data.pca = prcomp(X, scale=TRUE) 
W_matrix_original = data.pca$rotation[,1:32]
X = as.matrix(t(X))
  
positiveMatrix = testing-training
similarity <-  which(positiveMatrix[select.node.ID,]>0, arr.ind=TRUE)
disimilarity <- which(negativeMatrix[select.node.ID,]>0, arr.ind=TRUE)
n.disimilarity <- nrow(disimilarity)
disimi.s <- sample(n.disimilarity,nrow(similarity)*1)
disimilarity <- disimilarity[disimi.s,]

M = nrow(disimilarity)
N = nrow(similarity)
        
similarity <- t(similarity)
similarity <- data.frame(similarity)
#disimilarity <- NULL
disimilarity <- t(disimilarity)
disimilarity <- data.frame(disimilarity)

W_matrix = W_matrix_original
        
W<- try(W_train_P(W_matrix, K, X, similarity, disimilarity, r1, r2, beta, a, alpha, M, N))
if (class(W) =="try-error"){
  print("Error in training W ")
}else{
   W_matrix=W$W_matrix
   
   H = sign(crossprod(W_matrix,as.matrix(X))) # create hashing code -1 or +1
   Y = t((H+1)/2) # covert H to Y, -1 in H is 0 in Y and +1 in H is +1 in Y 
aH <- Agg_H (Y, K=16, M=2)
AUC <- get.AUC(testing, training, 1/aH,select.node.ID )
cat("AUC is ", AUC)
}



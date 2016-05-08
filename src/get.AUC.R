# 7-2-2015
# @author: Mengdi Wang
# update May 2nd, 2016
# This return the AUC of the proposed method
# testing is the testing matrix
# training is training matrix
# metricMatrix is the similairity matrix caculated from hashing tables
# select.node.ID are node id which we forcus, it can be the list of id of all nodes or some particular nodes.


get.AUC <- function(testing, training, metricMatrix,select.node.ID ){
  CorPositive = which(testing[select.node.ID,]-training[select.node.ID,]==1,arr.ind=TRUE) 
  CorPositive = data.frame(t(CorPositive))
  CorNegative= which(testing[select.node.ID,]+training[select.node.ID,]==0,arr.ind=TRUE) 
  CorNegative = data.frame(t(CorNegative))
  
  library(parallel)
  
  metricPositive <- mclapply(CorPositive, function(x, metricMatrix) get.metric.entry(x, metricMatrix), metricMatrix=metricMatrix[select.node.ID,], mc.cores=6)
  metricNegative<- mclapply(CorNegative, function(x, metricMatrix) get.metric.entry(x, metricMatrix), metricMatrix=metricMatrix[select.node.ID,], mc.cores=6)  
  
  n1.n2.list <- mclapply(metricPositive, function(x,metricNegative) get.n1.n2(x,metricNegative),metricNegative=metricNegative, mc.cores=6)
  
  test = matrix(unlist(n1.n2.list),ncol = 2, byrow = T)
  #n1.n2 <- Reduce("+", n1.n2.list)
  #AUC = (n1.n2[1]+0.5*n1.n2[2])/(length(metricPositive)*length(metricNegative))
  
  AUC = ((mean(test[,1])+0.5*(mean(test[,2]))))/length(metricNegative)
  return (AUC)
}

get.metric.entry <- function(cor,matrixMatrix){
  return (matrixMatrix[cor[1], cor[2]])
}

get.n1.n2 <- function(positve, metricNegative){
  n1 = length(which(metricNegative<positve))
  n2 = length(which(metricNegative==positve))
  
  return (c(n1,n2))
}



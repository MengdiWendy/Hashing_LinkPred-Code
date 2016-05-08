# This is a function to select nodes
# This function is optional 
get.select.node.ID <- function(testing, trainng, positiveMatrix, degree=2){
  degree_list <- rowSums(training)
  names(degree_list) <- rownames(training)
  select_node_ID_1 <- which(degree_list>degree)
  
  degree_list <- rowSums(positiveMatrix)
  names(degree_list) <- rownames(positiveMatrix)
  select_node_ID_2 <- which(degree_list>degree)
  
  select_node_ID = intersect(select_node_ID_1,select_node_ID_2) 
  
  select_node_ID
}

get.negative.positive.matrix <- function(testing, training){
  cor = which((testing-training)==1,arr.ind=TRUE) 
  testing.new = sparseMatrix(i = cor[,1], j = cor[,2], x=1, 
                             dims = dim(training), symmetric = F )
  
  rownames(testing.new) = rownames(training)
  colnames(testing.new) = rownames(training)
  
  cor = which((testing+training)==0,arr.ind=TRUE) 
  testing.negative = sparseMatrix(i = cor[,1], j = cor[,2], x=1, 
                                  dims = dim(training), symmetric = F )
  
  rownames(testing.new) = rownames(training)
  colnames(testing.new) = rownames(training)
  
  return (list(positiveMatrix =testing.new, negativeMatrix = testing.negative ))
}
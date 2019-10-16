
###AMOPLS
#R function implementing AMOPLS
#Based on Boccard & Rudaz (2016) as well as MatLab code kindly provided by J. Boccard

#KOPLS library
#install.package(kopls) #main package used in the analysis - Kernel OPLS
library(kopls)
library(plyr)
library(matrixcalc)


###Main function
AMOPLS_function<-function(
  Yall, #decomposed matrices
  npred, #number of predictive components to compute
  ortho, #number of orthogonal components to compute
  X = NULL,# Model matrix 
  UseX=FALSE ){# Use or not model matrix as response matrix

#Data preparation
Xres <- Yall[[length(Yall)]]
Y <- Yall[1:(length(Yall)-1)]

#Identify the number of matrices
f <- length(Y)

#Identify the number of observations
n <- nrow(Y[[1]])


#Create the response matrix
#############################
#Apply SVD to each submatrix and extract all non-zero eigenvectors = extract level barycentres
#Multiply U by the eigenvalues
#Then concatenate in responsematrix
#dec=number of decimals to round and decide if 0
responsematrix<-c()
Ysvd <- lapply(Y,svd)
for (i in 1:f){
  for(j in 1:n){
    if(round(Ysvd[[i]]$d[j],4)!=0){
      d <- Ysvd[[i]]$d[j]
      # d <- Ysvd[[i]]$d[j]/sum(Ysvd[[i]]$d)
    responsematrix <- cbind(responsematrix,Ysvd[[i]]$u[,j]*d)
      # responsematrix <- cbind(responsematrix,Y[[i]])
            }
        }
    }

# Use the model matrix as response
if(UseX==TRUE) 
{svdX=svd(X)
  responsematrix <- svdX$u%*%diag(svdX$d)}


###Create the predictor kernel matrix W_mat
#############################################
#Add the pure effects and residual submatrices = predictors Xi
predictors <- mapply(function(a) a + Xres, a = Y, SIMPLIFY = FALSE)

#Append the residual matrix to the list
predictors <- c(predictors,list(Xres))

#Create objects

AMat <- vector("list", length=(f+1))
W_mat <- matrix(0,n,n)

#Loop for f matrices
for (i in 1:(f+1)) {
#Create the kernel matrix with the augmented matrix as first matrix and no second matrix
  #p=polynomial kernel, 1=parameter of the kernel function

temp <- (predictors[[i]]%*%t(predictors[[i]]))

#Calculate the spectral norm of the matrix and normalise
AMat[[i]] <- temp/norm(x = temp, type = "F")

#Add all kernel matrices = predictor matrix W_mat
W_mat <- W_mat + AMat[[i]]
}


#Run KOPLS to explain the response by the predictor
#######################################################
#Nr of OPLS predictive components = Nr of responses in the Y matrix
#Number of orthogonal components identified by permutation

## Model
#Number of predictive components = length of Y
#Number of Y-orthogonal components = ortho
finalmodel <- koplsModel(W_mat, responsematrix, npred,ortho)
summary(finalmodel)


#Compute lambda
########################
lambda <- matrix(NA,nrow=length(Yall),ncol=npred+ortho)
for (j in 1:length(Yall)){

  for (k in 1:npred){
    lambda[j,k] <- finalmodel$T[,k] %*% AMat[[j]] %*% finalmodel$T[,k]
  }
  for (l in 1:ortho){
    lambda[j,k+l] <- finalmodel$To[,l] %*% AMat[[j]] %*% finalmodel$To[,l]
  }
}

#Normalise lambda
for (i in 1:ncol(lambda)){
  lambda[,i] <- lambda[,i]/sum(lambda[,i])
}

return(list(finalmodel,lambda, responsematrix, W_mat))
#End of AMOPLS function
}



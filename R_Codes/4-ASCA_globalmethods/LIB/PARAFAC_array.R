
##
#R code creating an array for use in PARAFASCA (Jansen et al. (2008))


PARAFAC_array<-function(
  data #decomposed matrices to be transformed into an array)
  ){
  #Create an array with the decomposed matrices: rows and columns = same as the matrices
  #third dimension = number of matrices
  PARAFAC_array<-array(data=NA,dim=c(nrow(data[[1]]),ncol(data[[1]]),length(data)))
    for (i in 1:length(data)) {
        PARAFAC_array[,,i]<-data[[i]]
    }
  return(PARAFAC_array)
}




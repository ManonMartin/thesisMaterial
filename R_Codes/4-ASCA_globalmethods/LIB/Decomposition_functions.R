

##########Wrapping function bringing together the GLM and ANOVA decompositions
#in one function and based on the matrixDecomposition code written by mthiel & rvanoirbeek


decomposition<-function(
  data, #matrix to decompose
  design, # experimental design
  decomp.meth, #"ANOVA" or "GLM"
  formula=FALSE,#formula for the GLM decomposition
  perm.calc=FALSE, #permutation calculation for the GLM decomposition
  graphs=FALSE, #"ASCA" or "APCA" for the ANOVA decomposition
  decomp.name="dataset" #name of decomposition for ANOVA decomposition
){
  
  if(decomp.meth=="ANOVA"){
    Ydecomp<-ASCA_APCA_crossed(data,design,decomposition=decomp.meth,method=graphs,name=decomp.name)
    #Get the row names for all matrices
    rownames(Ydecomp$alpha)=rownames(Ydecomp$beta)=rownames(Ydecomp$alphabeta)<-rownames(Ydecomp$epsilon)
    #Y and Xres objects with the original data order
    Ypure<-list(Yalpha=Ydecomp$alpha,Ybeta=Ydecomp$beta,Yalphabeta=Ydecomp$alphabeta)
    Xres<-Ydecomp$epsilon
  }
  
  if(decomp.meth=="GLM"){
    Ydecomp<-matrixDecomposition(formula,data,design)
    Ypure<-list(Yalpha=Ydecomp$effectMatrices[[2]],Ybeta=Ydecomp$effectMatrices[[3]],Yalphabeta=Ydecomp$effectMatrices[[4]])
    Xres<-Ydecomp$residuals
  }
  
  #Augmented matrices
  Yaugmented<-mapply(function(a) a+Xres, a = Ypure, SIMPLIFY = FALSE)
  #Append the residual matrix to the list
  Ypure<-c(Ypure,list(Xres))
  Yaugmented<-c(Yaugmented,list(Xres))
  #3 objects: pure effect matrices, augmented effect matrices, residual matrix
  return(list(Ypure=Ypure,Yaugmented=Yaugmented,Yres=Xres))
}



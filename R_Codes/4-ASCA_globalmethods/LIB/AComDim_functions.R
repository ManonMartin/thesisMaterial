
#####ComDim
#R function implementing ComDim
#Based on Jouan-Rimbaud Bouveresse et al. (2011) 
#and the MatLab ComDim function from the saisir package (https://www.chimiometrie.fr/saisir_webpage.html)

#ComDim Function
#Function to compute r common dimensions
ComDim<-function(
  Y, #original data as a list of matrices
  R, #number of common dimensions to compute
  thresh=1E-20, #convergence threshold
  stand=TRUE #Standardisation or not
  ) 
  {
  #Identify the number of matrices
  K<-length(Y)
  #Identify the number of observations
  n<-nrow(Y[[1]])
  # Center and normalize data by Frobenius Norm
  YS=Y
  for(k in 1:K)
  {
  YS[[k]]=  scale(Y[[k]], center = TRUE, scale = FALSE)
  YS[[k]]=YS[[k]]/(norm(YS[[k]], type = "F"))
  }  
  #Create the objects needed to store intermediary and final results
  lambda<-matrix(1,nrow=K,ncol=R)
  q<-matrix(NA,nrow=n,ncol=R)
  muvec=rep(0,R)
  muvec2=rep(0,R)
  expl<-rep(0,R)
  iter<-rep(0,R)
  firstCD<-TRUE
  Sigma0<-list()
  TotInerty=0
  for(k in 1:K){
    #Computation starting SigmaK
    Sigma0[[k]]<-YS[[k]]%*%t(YS[[k]])
    TotInerty<- TotInerty+sum((Sigma0[[k]]*Sigma0[[k]]))
  }
  
  for(r in 1:R) {
    #Create the objects needed
    unfit<-999999
    Aux<-list()
    Auxsum<-0
    conv<-0
    #########Iteration of the computation of SigmaG, q and lambda (see AComDim article)
    repeat {
      #Keep the last value of unfit to check the convergence
      lastunfit<-unfit
      ###Compute the global sample variance covariance matrix SigmaG
      SigmaG<-0
      if(firstCD==TRUE) 
      {Sigma<-Sigma0
        YR<-YS} 
      else 
      {Sigma<-newSigma}
      for(k in 1:K){
        #Computation of SigmaG
        SigmaG<-SigmaG+lambda[k,r]*Sigma[[k]]
      }
      #Singular Value Decomposition of SigmaG
      SigmaSVD<-svd(SigmaG)
      #Allocate the value of the first singular vector to q
      q[,r]<-SigmaSVD$u[,1]
      #Update the value of lambda and compute Auxsum
      for(k in 1:K){
        lambda[k,r]<-as.numeric(t(q[,r])%*%Sigma[[k]]%*%q[,r])
        Aux[[k]]<-Sigma[[k]]-lambda[k,r]*(q[,r]%*%t(q[,r]))
        Auxsum[k]<-sum(Aux[[k]]*Aux[[k]])
      }
      #Compute unfit
      unfit<-sum(Auxsum)
      #Compare new unfit and last unfit to check the convergence
      conv<-lastunfit-unfit
      iter[r]<-iter[r]+1
      if (conv<thresh|iter[r]>200) {break}
    } 
    
    #Compute the explained variance
    muvec2[r]=q[,r]%*%SigmaG %*%q[,r]
    muvec[r]=sum(lambda[,r]^2)
    expl[r]=100*muvec[r]/TotInerty
    #Run ComDim1
    newq<-q[,r]
    
    #Update the deflated matrix
    #(I-qq')X[(I-qq')X]' => (I-qq')Sigma(I-qq')'
    newSigma<-list()
    newSigma2<-rep()
    for (k in 1:K) {
      YR[[k]]=YR[[k]]-newq%*%t(newq)%*%YR[[k]]
      newSigma[[k]]<-YR[[k]]%*%t(YR[[k]])}
    #Indication that not first CD anymore
    firstCD<-FALSE
    #End of the r loop
  }
  #Outputs of the ComDim analysis
  output<-list(n,K,lambda,q,expl,muvec,muvec2,iter,TotInerty,Sigma0,YS)
  names(output)<-c("n","K","lambda","q","expl","muvec","muvec2","iter","TotInerty","Sigma0","YS")
  return(output)
  
  }


#############
#AComDim plots

AComDim_plots<-function(
  res #result of a ComDim analysis
  ){
#Results of the ComDim function for the decomposed matrices
saliences<-res$lambda
common_dimensions<-res$q #also called loadings
explainedvar<-res$expl
gammavec=res$gammavec
K<-res$K
colnames(saliences)=colnames(common_dimensions)<-paste("CD", 1:ncol(saliences), sep="")
rownames(saliences)<-paste("Matrix",1:nrow(saliences),sep="")

#Plots
#Explained variance
par(mfrow=c(1,1),mar=c(3,2,3,2))
a<-barplot((explainedvar),ylim=c(0,110),names.arg=1:R,yaxt='n',main="Explained variance",col="coral2")
text(x=a,y=explainedvar+4,labels=as.character(round(explainedvar)))

#Cumulative Explained variance
par(mfrow=c(1,1),mar=c(3,2,3,2))
explainedvar=cumsum(explainedvar)
a<-barplot(explainedvar,ylim=c(0,110),names.arg=1:R,yaxt='n',main="Cumulative Explained variance",col="coral2")
text(x=a,y=explainedvar+4,labels=as.character(round(explainedvar,1)))
abline(h=100)

#Saliences by matrix for all common dimensions computed
par(mfrow=c(R/2,2),mar=c(2.5,2,2.5,2))
for(i in 1:R) {
  barplot(saliences[,i],main=paste("CD",i),names.arg=c("A","B","AB","E"),ylim=c(0,1),yaxt='n',col="coral2")
  axis(2,at=c(0,1),labels=c("0","1"))
}
}




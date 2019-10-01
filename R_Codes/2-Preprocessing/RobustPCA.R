
##########################################################################
## Robust PCA and ROBPCA for outlier analysis in a multivariate asset
##########################################################################
#install.packages("rrcov")

# k1 and k2 enable to perform ROBPCA on distinct groups of observations
# if group is not NULL

RobustPCA <- function(x, dataname="data", group=NULL, k=0, kgroup=0, impG=FALSE, Out.path="Graphics") { 

# Arguments
  
  # k :   nCP for raw data , if k=0, the algorithm itself will determine the number of components
  # kgroup : vector (or unique value) of nCP to keep for each group


######################
Groups = levels(as.factor(group))
nGroup=length(unique(group))
n=dim(x)[1]

if(length(kgroup==1)){
  cat("\n Chosen default kgroup for all groups: ", kgroup)
  kgroup=rep(kgroup, length(unique(group)))
}


# nbre spectre par groupe
ng=numeric()
for (i in 1:nGroup){
  ng[i]=dim(x[which(group==unique(group)[i]),])[1]
}

# Raw data
X.HubPCA <- PcaHubert(x, k = k,  kmax=min(dim(x)), mcd=FALSE, scale = FALSE) # n CP choisi par l'algorithme
nCP.X=X.HubPCA@k
 #  # number of CP
summary(X.HubPCA)



xgroup=vector("list", nGroup)
Xg.HubPCA=vector("list", nGroup)
subgroup=c()
nsubgroup=numeric()
nCP.Xg=numeric()

if(!is.null(group)) {
  for (i in 1:nGroup) { 
    xgroup[[i]]=x[which(group==unique(group)[i]),]
    subgroup[i]=as.character(unique(group)[i])
#     x1=x[group==unique(group)[1],]
#     group1=as.character(unique(group)[1])
#   x2=x[group==unique(group)[2],]
#   group2=as.character(unique(group)[2])
#    n[i]=dim(xgroup[[i]])[1]
#   n2=dim(x2)[1]

    Xg.HubPCA[[i]] <- PcaHubert(xgroup[[i]], k = kgroup[i],  kmax=min(dim(xgroup[[i]])), scale = FALSE, mcd=FALSE)
    nCP.Xg[i]=Xg.HubPCA[[i]]@k

#     X2.HubPCA <- PcaHubert(x2, k = k2,  kmax=min(dim(x2)), scale = FALSE, mcd=FALSE)
#     nCP.X2=X2.HubPCA@k
  
#   if (impG==TRUE) { 
#   pca.scoreplot(X1.HubPCA, main=paste0("Group1 - nPC: ", nCP.X1),id.n=5)
#   pca.scoreplot(X2.HubPCA, main=paste0("Group2 - nPC: ", nCP.X2),id.n=5)
#   }
}
}



## Good, bad leverage points and Orthogonal outlier
# SD = score distance, corresponds to mahalanobis distance

X.glp=X.blp=X.oo=numeric()

X.glp=which(X.HubPCA@sd>=X.HubPCA@cutoff.sd) 
X.oo=which(X.HubPCA@od>=X.HubPCA@cutoff.od)
X.blp=X.glp[X.glp %in% X.oo]


Xg.glp=vector("list", nGroup)
Xg.oo=vector("list", nGroup)

if(!is.null(group)) { 
  for (i in 1:nGroup) { 
  Xg.glp[[i]]=which(Xg.HubPCA[[i]]@sd>=Xg.HubPCA[[i]]@cutoff.sd)
  Xg.oo[[i]]=which(Xg.HubPCA[[i]]@od>=Xg.HubPCA[[i]]@cutoff.od) 

# X1.glp=which(X1.HubPCA@sd>=X1.HubPCA@cutoff.sd)
# X1.oo=which(X1.HubPCA@od>=X1.HubPCA@cutoff.od) 
# X2.glp=which(X2.HubPCA@sd>=X2.HubPCA@cutoff.sd) 
# X2.oo=which(X2.HubPCA@od>=X2.HubPCA@cutoff.od) 
}
}


X.GLP=names(X.glp)
X.OO=names(X.oo)


Xg.GLP=vector("list", nGroup)
Xg.OO=vector("list", nGroup)

if(!is.null(group)) { 
  for (i in 1:nGroup) { 
Xg.GLP[[i]]=names(Xg.glp[[i]])
Xg.OO[[i]]=names(Xg.oo[[i]])
}
# X2.GLP=names(X2.glp)
# X2.OO=names(X2.oo)
res=list(X.GLP=X.GLP, Xg.GLP=Xg.GLP,  X.OO=X.OO, Xg.OO=Xg.OO)

} else {res=list(X.GLP, X.OO)}






if (impG==TRUE) { 

  
# Outlier map based on ROBPCA model.
  
pdf(paste0(Out.path,"/", dataname, "OD_OutlierMapsRaw.pdf"), width = 6, height = 6)
 
xlim=c(0, max(X.HubPCA@sd, X.HubPCA@cutoff.sd)*1.1)
ylim=c(0, max(X.HubPCA@od, X.HubPCA@cutoff.od)*1.1)
 
name=names(X.HubPCA@od)
namelab=substr(name, 1,5)

col=rep(1,n)
col[X.HubPCA@flag==F]=2

plot(X.HubPCA@sd, X.HubPCA@od, ylim=ylim, xlim=xlim, col=col, xlab="Score distance", main=paste0("Robust PCA (Raw data) - n PC: ", nCP.X), 
      ylab="Orthogonal distance")
abline(h=X.HubPCA@cutoff.od, lty=1, col=2)
abline(v=X.HubPCA@cutoff.sd, lty=1, col=2)
text(X.HubPCA@sd, X.HubPCA@od, labels=namelab, pos=c(1,3),cex = 0.8, col=col)
  
dev.off()
  
## OD_OutlierMapsGroup
###########################

if(!is.null(group)) { 
  pdf(paste0(Out.path,"/", dataname, "OD_OutlierMapsGroup.pdf"), width = ceiling(nGroup)*6, height = 6)
  par(mfrow=c(ceiling(nGroup/2),2))

  for (i in 1:nGroup) {     
  #windows()
  xlimg=c(0, max(Xg.HubPCA[[i]]@sd, Xg.HubPCA[[i]]@cutoff.sd)*1.3)
  ylimg=c(0, max(Xg.HubPCA[[i]]@od, Xg.HubPCA[[i]]@cutoff.od)*1.3)
    
  name=names(Xg.HubPCA[[i]]@od)
  namelab=substr(name, 1,5)
  
  col=rep(1,ng[i])
  col[Xg.HubPCA[[i]]@flag==F]=2
 
  plot(Xg.HubPCA[[i]]@sd, Xg.HubPCA[[i]]@od, ylim=ylimg, xlim=xlimg, col=col, xlab="Score distance", 
       main=paste0("Robust PCA (Group ",unique(group)[i],") - n PC: ", nCP.Xg[i] ), ylab="Orthogonal distance")
  abline(h=Xg.HubPCA[[i]]@cutoff.od, lty=1, col=2)
  abline(v=Xg.HubPCA[[i]]@cutoff.sd, lty=1, col=2)
  text(Xg.HubPCA[[i]]@sd, Xg.HubPCA[[i]]@od, labels=namelab, pos=c(1,3),cex = 0.8, col=col)   
  }
  dev.off()
}



## OD_Outlier_PCAScoreplotRaw
####################################

# pdf(paste0(Out.path,"/", dataname, "OD_Outlier_PCAScoreplotRaw.pdf"), width = 6, height = 6)
# pca.scoreplot(X.HubPCA, main=paste0("Raw data - nPC: ", nCP.X))
# dev.off()
#   


# Plot des outliers "good leverage points" et "orthogonal outlier" avec le spectre moyen par groupe
##############################################################################
# LEs spectres sont selectionnes sur base des PCA robustes globales et par groupe  

#BLP=unique(unlist(res))

# Raw data
#############
MeanSpect=apply(x,2, mean)

# X.glp
# X.oo
# 
# X.GLP=names(X.glp)
# X.OO=names(X.oo)

## GLP

if(!length(X.glp)==0) { 

pdf(paste0(Out.path,"/", dataname, "glpVSMeanSpecRaw.pdf"), width=7, height=7*length(X.glp))
par(mfrow=c(length(X.glp),1))

for (i in 1:length(X.glp)) {
  plot(MeanSpect, type="l", main=paste0("Good leverage points", "\n ", X.GLP[i]))
  
  lines(as.numeric(x[X.glp[i],]), col="red", type="l")
  legend("topright", legend=c("Ref spectrum", X.GLP[i]), col=c(1, "red"), lty=1)
  
} 
dev.off()
}

## OO
if(!length(X.oo)==0) { 
  pdf(paste0(Out.path,"/", dataname, "ooVSMeanSpecRaw.pdf"), width=7, height=7*length(X.oo))
  par(mfrow=c(length(X.oo),1))
  
  for (i in 1:length(X.oo)) {
    plot(MeanSpect, type="l", main=paste0("Orthogonal outlier", "\n ", X.OO[i]))
    
    lines(as.numeric(x[X.oo[i],]), col="red", type="l")
    legend("topright", legend=c("Ref spectrum", X.OO[i]), col=c(1, "red"), lty=1)
    
  } 
  dev.off()
}



# By group
#############
xg=vector("list", nGroup)
groupp=vector("list", nGroup)
MeanSpectG=vector("list", nGroup)

if(!is.null(group)) { 
  for (i in 1:nGroup) {  

    xg[[i]]=x[which(group==unique(group)[i]),]
    groupp[[i]]=as.character(unique(group)[i])
    }  
}
#MeanSpectG[[i]]=apply(xg[[i]],2, mean)


## GLP and OO
sum=length(unlist(Xg.glp))+length(unlist(Xg.oo))
if(!sum==0) { 
  
  pdf(paste0(Out.path,"/", dataname, "glpVSMeanSpecGroup.pdf"), width=8, height=6*sum)
  par(mfrow=c(sum,1))
  
  for (i in 1:nGroup) {  

    if(!length(Xg.glp[[i]])==0) { 
  for (j in 1:length(Xg.glp[[i]])) {
    Moyg=apply(xg[[i]][-Xg.glp[[i]][j],],2, mean)
    plot(Moyg, type="l", main=paste0("Good leverage point Group",Groups[i],"\n ", Xg.GLP[[i]][j]))
    
    lines(as.numeric(xg[[i]][Xg.glp[[i]][j],]), col="red", type="l")
    legend("topright", legend=c("Ref spectrum", Xg.GLP[[i]][j]), col=c(1, "red"), lty=1) 
  } 
}

    if(!length(Xg.oo[[i]])==0) { 
  for (j in 1:length(Xg.oo[[i]])) {
    Moyg=apply(xg[[i]][-Xg.oo[[i]][j],],2, mean)
    plot(Moyg, type="l", main=paste0("OO Group",unique(group)[i], "\n ", Xg.OO[[i]][j]))
    
    lines(as.numeric(xg[[i]][Xg.oo[[i]][j],]), col="red", type="l")
    legend("topright", legend=c("Ref spectrum", Xg.OO[[i]][j]), col=c(1, "red"), lty=1)
    
  } 
}



}
dev.off()
}

}

return(res)
}

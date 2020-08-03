# repret/reprod mixed

######################
# A. Balanced case
######################
out_path <- "/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments"
# save(res.parlmer, outcomes, design,
#    file=file.path(out_path,"res_parlmer_metabo_balanced.RData"))

rm(list=ls())
# load("/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments/sensoryData.RData")
out_path <- "/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments"
load(file.path(out_path,"res_parlmer_metabo_balanced.RData"))
require(lme4)
n <- nrow(outcomes)

MM_full <- res.parlmer
### all fixed estimates and random predictions
# cof_PC <- rbind(fixef_PC, ranef_PC)

# Model matrices ------------------------------
# Fixed effects 
fixef_PC <- sapply(MM_full$merMod_obj, fixef)
FixedModMatlist <- res.parlmer$FixedModMatlist
X   <- do.call(cbind, FixedModMatlist)
X  <- X [,rownames(fixef_PC)] # reorder colnames of X 

# random effects 
length(unique(design$Volunteer))
length(unique(design$Sampling))

RanModMatlist <- res.parlmer$RanModMatlist
dim1RandModMad <- sapply(RanModMatlist, function(x) dim(x)[2])
names_randomEffects <- names(RanModMatlist)
shortRandNames <- gsub("[^A-z]", "", names_randomEffects)
Z <- do.call(cbind, RanModMatlist)
colnames(Z) <- paste0(rep(shortRandNames, dim1RandModMad),colnames(Z))


# Variance components ------------------------------
### Residuals and random judge and CJ variance
varcor_random_full <- sapply(MM_full$merMod_obj, VarCorr) # variances
varcor_random_full <- matrix(unlist(varcor_random_full), nrow=nrow(varcor_random_full), 
                             ncol=ncol(varcor_random_full), dimnames=dimnames(varcor_random_full))

# sqrt(varcor_random_full)
# VarCorr(lmer_res)


# Residuals sd error
Res_std_error_PC <- sapply(MM_full$merMod_obj, sigma) # Std.Dev.
varcor_resid <- Res_std_error_PC^2 # variances 


# ranef --------------
ranef_PC <- sapply(MM_full$merMod_obj, 
                   function(x) unlist(ranef(x, condVar=FALSE)))

#######################################
#          pour 1 reponse #############
#######################################
dev.off()

lmer_res <- res.parlmer$merMod_obj[[1]]
y <- outcomes[,1]

# yhat 
yhat <- predict(lmer_res)

varcor_random_full
varcor_resid 
sigma2r <- varcor_resid[1]
sigma2v <- varcor_random_full["Volunteer" ,1]
sigma2vs <- varcor_random_full["Volunteer:Sampling",1]


q = ncol(Z)
G <- diag(q)*0
G[1:27,1:27] = diag(27)*sigma2vs
G[28:36,28:36] = diag(9)*sigma2v

R <- diag(n)*sigma2r
Rinv <- diag(n)*(1/sigma2r)

plot(yhat,y-yhat,main="Residus")

# Equ (2) Eilers 2018
XRX <- t(X)%*%Rinv%*%X
ZRX <- t(Z)%*%Rinv%*%X
XRZ <- t(X)%*%Rinv%*%Z
ZRZ=t(Z)%*%Rinv%*%Z
ZRZG <- ZRZ + solve(G)
XRY=t(X)%*%Rinv%*%y
ZRY=t(Z)%*%Rinv%*%y
Q <- solve(rbind(cbind(XRX, XRZ), cbind(ZRX, ZRZG)))
vecbc=Q%*%rbind(XRY,ZRY)
cat("\n Estimation des paramètres et effets aléatoires")
print(vecbc)

# cat("Comparaison aux estimateurs R lmer")
cat("beta")
print(cbind(vecbc[1:3],lmer_res@beta))

cat("C - judges")
print(cbind(vecbc[4:30],ranef_PC[1:27,1])) # Volunteer:Sampling
cat("C - candies*judges")
print(cbind(vecbc[31:39],ranef_PC[28:36,1])) # Volunteer
# plot(vecbc[6:16],lmer_res@u)
# print(lm(vecbc[6:16]~lmer_res@u))

# Equ (7) Eilers 2018: compute H
H <- cbind(X,Z)%*%Q%*%rbind(t(X),t(Z))%*%Rinv
cat("Effective dimension : ", sum(diag(H)),"\n")

# Comparaison yhat lmer et Hy
yhat2=H%*%y
# print(lm(yhat2~yhat))
plot(yhat,yhat2,main="Y predit par lmer et par Hy")
plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")


# dimensions effectives sur base de K
K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
# dim(K)
K_diag <- diag(K)
# sum(K_diag)

K11=K[1:3,1:3] # Tube+time+intercept
K22=K[4:30,4:30] # # Volunteer:Sampling
K33=K[31:39,31:39] # # Volunteer
cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")

ED_Tube_time <- sum(diag(K11))
ED_Samp <- sum(diag(K22))
ED_Vol <- sum(diag(K33))

# ED bien calculées ??
sum(ED_Tube_time,ED_Samp, ED_Vol) - sum(K_diag)


#########################################
# Pour toutes les réponses ##############
#########################################

ED_Tube_time <- c()
ED_Samp <- c()
ED_Vol <- c()

for (i in c(1:12, 14:15)){
  
  lmer_res <- res.parlmer$merMod_obj[[i]]
  y <- outcomes[,i]
  
  # yhat 
  yhat <- predict(lmer_res)
  
  sigma2r <- varcor_resid[i]
  sigma2v <- varcor_random_full["Volunteer" ,i]
  sigma2vs <- varcor_random_full["Volunteer:Sampling",i]
  
  
  q = ncol(Z)
  G <- diag(q)*0
  G[1:27,1:27] = diag(27)*sigma2vs
  G[28:36,28:36] = diag(9)*sigma2v
  R <- diag(n)*sigma2r
  Rinv <- diag(n)*(1/sigma2r)
  
  # Equ (2) Eilers 2018
  XRX <- t(X)%*%Rinv%*%X
  ZRX <- t(Z)%*%Rinv%*%X
  XRZ <- t(X)%*%Rinv%*%Z
  ZRZ=t(Z)%*%Rinv%*%Z
  ZRZG <- ZRZ + solve(G)
  XRY=t(X)%*%Rinv%*%y
  ZRY=t(Z)%*%Rinv%*%y
  Q <- solve(rbind(cbind(XRX, XRZ), cbind(ZRX, ZRZG)))
  vecbc=Q%*%rbind(XRY,ZRY)
  # cat("\n Estimation des paramètres et effets aléatoires")
  # print(vecbc)
  # cat("Comparaison aux estimateurs R lmer")
  # cat("beta")
  
  # print(cbind(vecbc[1:5],lmer_res@beta))
  # cat("C")
  # print(cbind(vecbc[6:16],lmer_res@u)) # FIXME
  # plot(vecbc[6:16],lmer_res@u)
  # print(lm(vecbc[6:16]~lmer_res@u))
  
  # Equ (7) Eilers 2018: compute H
  H <- cbind(X,Z)%*%Q%*%rbind(t(X),t(Z))%*%Rinv
  cat("Effective dimension : ", sum(diag(H)),"\n")
  
  # Comparaison yhat lmer et Hy
  yhat2=H%*%y
  # print(lm(yhat2~yhat))
  plot(yhat,yhat2,main="Y predit par lmer et par Hy")
  plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")
  
  
  # dimensions effectives sur base de K
  K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
  # dim(K)
  K_diag <- diag(K)
  # sum(K_diag)
  
  K11=K[1:3,1:3] # Tube+time+intercept
  K22=K[4:30,4:30] # # Volunteer:Sampling
  K33=K[31:39,31:39] # # Volunteer
  cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
  cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
  cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
  cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")

  
  ED_Tube_time[i] <- sum(diag(K11))
  ED_Samp[i] <- sum(diag(K22))
  ED_Vol[i] <- sum(diag(K33))
}

#######################################
#       PC13 : only Sampling ##########
#######################################

lmer_res <- res.parlmer$merMod_obj[[13]]
y <- outcomes[,13]

# yhat 
yhat <- predict(lmer_res)

varcor_random_full
varcor_resid 
sigma2r <- varcor_resid[13]
sigma2v <- varcor_random_full["Volunteer" ,13]
sigma2vs <- varcor_random_full["Volunteer:Sampling",13]

ZZ <- Z[,1:27]

q = ncol(ZZ)
G <- diag(q)*0
G = diag(q)*sigma2vs

R <- diag(n)*sigma2r
Rinv <- diag(n)*(1/sigma2r)

plot(yhat,y-yhat,main="Residus")

# Equ (2) Eilers 2018
XRX <- t(X)%*%Rinv%*%X
ZRX <- t(ZZ)%*%Rinv%*%X
XRZ <- t(X)%*%Rinv%*%ZZ
ZRZ=t(ZZ)%*%Rinv%*%ZZ
ZRZG <- ZRZ + solve(G)
XRY=t(X)%*%Rinv%*%y
ZRY=t(ZZ)%*%Rinv%*%y
Q <- solve(rbind(cbind(XRX, XRZ), cbind(ZRX, ZRZG)))
vecbc = Q %*% rbind(XRY,ZRY)
cat("\n Estimation des paramètres et effets aléatoires")
print(vecbc)

# cat("Comparaison aux estimateurs R lmer")
cat("beta")
print(cbind(vecbc[1:3],lmer_res@beta))

cat("C - Volunteer:Sampling")
print(cbind(vecbc[4:30], ranef_PC[1:27,1])) # Volunteer:Sampling

# plot(vecbc[6:16],lmer_res@u)
# print(lm(vecbc[6:16]~lmer_res@u))

# Equ (7) Eilers 2018: compute H
H <- cbind(X,ZZ)%*%Q%*%rbind(t(X),t(ZZ))%*%Rinv
cat("Effective dimension : ", sum(diag(H)),"\n")

# Comparaison yhat lmer et Hy
yhat2=H%*%y
# print(lm(yhat2~yhat))
plot(yhat,yhat2,main="Y predit par lmer et par Hy")
plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")


# dimensions effectives sur base de K
K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
# dim(K)
K_diag <- diag(K)
# sum(K_diag)

K11=K[1:3,1:3] # Tube+time+intercept
K22=K[4:30,4:30] # # Volunteer:Sampling
K33=0# # Volunteer
cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")

ED_Tube_time[13] <- sum(diag(K11))
ED_Samp[13] <- sum(diag(K22))
ED_Vol[13] <- sum(diag(K33))


####################
ED <- rbind(ED_Tube_time, ED_Samp, ED_Vol)
colnames(ED) <- names(MM_full$merMod_obj)
colSums(ED)


pdf("/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments/ED_metabo.pdf",
    width = 12, height = 7, pointsize = 16)
par(mar=c(5,4,4,4))
plot(1:ncol(ED),ED_Tube_time, ylim=range(ED), xaxt="n",
     type="b", col="red", lty=4, xlab="PC", ylab = "ED")
lines(1:ncol(ED),ED_Samp, type="b", col="blue", lty=2)
lines(1:ncol(ED),ED_Vol, type="b", col="blue", lty=1)
legend("top",legend = c("Tube+time", "Sampling", "Volunteer"),
       col = c("red", "blue", "blue"), lty=c(4,2,1), ncol = 3,
       bty = "n", xpd=TRUE, inset = c(0,-0.2))
axis(side = 1, at = 1:ncol(ED), labels = 1:ncol(ED))
mtext("EDtot:", side=1, at = 0, line=4)
mtext(format(round(colSums(ED),1), nsmall=1), side=1, at = 1:ncol(ED), line=4, cex=0.8)
dev.off()

ED_metabo <- ED
save(ED_metabo, file=file.path(out_path, "ED_metabo.RData"))
write.csv(ED, file.path(out_path,"ED_metabo.csv"))


##############################
#   B. Unbalanced case #######
##############################

out_path <- "/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments"
# save(res.parlmer, outcomes, design,
#      file=file.path(out_path,"res_parlmer_metabo_unbalanced.RData"))

rm(list=ls())
# load("/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments/sensoryData.RData")
out_path <- "/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments"
load(file.path(out_path,"res_parlmer_metabo_unbalanced.RData"))
require(lme4)
n <- nrow(outcomes)

MM_full <- res.parlmer
### all fixed estimates and random predictions
# cof_PC <- rbind(fixef_PC, ranef_PC)

# Model matrices ------------------------------
# Fixed effects 
fixef_PC <- sapply(MM_full$merMod_obj, fixef)
FixedModMatlist <- res.parlmer$FixedModMatlist
X   <- do.call(cbind, FixedModMatlist)
X  <- X [,rownames(fixef_PC)] # reorder colnames of X 

# random effects 
length(unique(design$Volunteer))
length(unique(design$Sampling))

RanModMatlist <- res.parlmer$RanModMatlist
dim1RandModMad <- sapply(RanModMatlist, function(x) dim(x)[2])
names_randomEffects <- names(RanModMatlist)
shortRandNames <- gsub("[^A-z]", "", names_randomEffects)
Z <- do.call(cbind, RanModMatlist)
colnames(Z) <- paste0(rep(shortRandNames, dim1RandModMad),colnames(Z))


# Variance components ------------------------------
### Residuals and random judge and CJ variance
varcor_random_full <- sapply(MM_full$merMod_obj, VarCorr) # variances
varcor_random_full <- matrix(unlist(varcor_random_full), nrow=nrow(varcor_random_full), 
                             ncol=ncol(varcor_random_full), dimnames=dimnames(varcor_random_full))

# sqrt(varcor_random_full)
# VarCorr(lmer_res)


# Residuals sd error
Res_std_error_PC <- sapply(MM_full$merMod_obj, sigma) # Std.Dev.
varcor_resid <- Res_std_error_PC^2 # variances 


# ranef --------------
ranef_PC <- sapply(MM_full$merMod_obj, 
                   function(x) unlist(ranef(x, condVar=FALSE)))

#######################################
#          pour 1 reponse #############
#######################################
dev.off()

lmer_res <- res.parlmer$merMod_obj[[1]]
y <- outcomes[,1]

# yhat 
yhat <- predict(lmer_res)

varcor_random_full
varcor_resid 
sigma2r <- varcor_resid[1]
sigma2v <- varcor_random_full["Volunteer" ,1]
sigma2vs <- varcor_random_full["Volunteer:Sampling",1]


q = ncol(Z)
G <- diag(q)*0
G[1:36,1:36] = diag(36)*sigma2vs # Volunteer:Sampling
G[37:48,37:48] = diag(12)*sigma2v # Volunteer

R <- diag(n)*sigma2r
Rinv <- diag(n)*(1/sigma2r)

plot(yhat,y-yhat,main="Residus")

# Equ (2) Eilers 2018
XRX <- t(X)%*%Rinv%*%X
ZRX <- t(Z)%*%Rinv%*%X
XRZ <- t(X)%*%Rinv%*%Z
ZRZ=t(Z)%*%Rinv%*%Z
ZRZG <- ZRZ + solve(G)
XRY=t(X)%*%Rinv%*%y
ZRY=t(Z)%*%Rinv%*%y
Q <- solve(rbind(cbind(XRX, XRZ), cbind(ZRX, ZRZG)))
vecbc=Q%*%rbind(XRY,ZRY)
cat("\n Estimation des paramètres et effets aléatoires")
print(vecbc)

# cat("Comparaison aux estimateurs R lmer")
cat("beta")
print(cbind(vecbc[1:3],lmer_res@beta))

cat("C - Volunteer:Sampling")
print(cbind(vecbc[4:39],ranef_PC[1:36,1])) # Volunteer:Sampling
cat("C - Volunteer")
print(cbind(vecbc[40:51],ranef_PC[37:48,1])) # Volunteer
# plot(vecbc[6:16],lmer_res@u)
# print(lm(vecbc[6:16]~lmer_res@u))

# Equ (7) Eilers 2018: compute H
H <- cbind(X,Z)%*%Q%*%rbind(t(X),t(Z))%*%Rinv
cat("Effective dimension : ", sum(diag(H)),"\n")

# Comparaison yhat lmer et Hy
yhat2=H%*%y
# print(lm(yhat2~yhat))
plot(yhat,yhat2,main="Y predit par lmer et par Hy")
plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")


# dimensions effectives sur base de K
K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
# dim(K)
K_diag <- diag(K)
# sum(K_diag)

K11=K[1:3,1:3] # Tube+time+intercept
K22=K[4:39,4:39] # # Volunteer:Sampling
K33=K[40:51,40:51] # # Volunteer
cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")

ED_Tube_time <- sum(diag(K11))
ED_Samp <- sum(diag(K22))
ED_Vol <- sum(diag(K33))

# ED bien calculées ??
sum(ED_Tube_time,ED_Samp, ED_Vol) - sum(K_diag)


#########################################
# Pour toutes les réponses ##############
#########################################


ED_Tube_time <- c()
ED_Samp <- c()
ED_Vol <- c()

for (i in c(1:15)){
  
  lmer_res <- res.parlmer$merMod_obj[[i]]
  y <- outcomes[,i]
  
  # yhat 
  yhat <- predict(lmer_res)
  
  sigma2r <- varcor_resid[i]
  sigma2v <- varcor_random_full["Volunteer" ,i]
  sigma2vs <- varcor_random_full["Volunteer:Sampling",i]
  
  
  q = ncol(Z)
  G <- diag(q)*0
  G[1:36,1:36] = diag(36)*sigma2vs # Volunteer:Sampling
  G[37:48,37:48] = diag(12)*sigma2v # Volunteer
  
  R <- diag(n)*sigma2r
  Rinv <- diag(n)*(1/sigma2r)
  
  # Equ (2) Eilers 2018
  XRX <- t(X)%*%Rinv%*%X
  ZRX <- t(Z)%*%Rinv%*%X
  XRZ <- t(X)%*%Rinv%*%Z
  ZRZ=t(Z)%*%Rinv%*%Z
  ZRZG <- ZRZ + solve(G)
  XRY=t(X)%*%Rinv%*%y
  ZRY=t(Z)%*%Rinv%*%y
  Q <- solve(rbind(cbind(XRX, XRZ), cbind(ZRX, ZRZG)))
  vecbc=Q%*%rbind(XRY,ZRY)
  # cat("\n Estimation des paramètres et effets aléatoires")
  # print(vecbc)
  # cat("Comparaison aux estimateurs R lmer")
  # cat("beta")
  
  # print(cbind(vecbc[1:5],lmer_res@beta))
  # cat("C")
  # print(cbind(vecbc[6:16],lmer_res@u)) # FIXME
  # plot(vecbc[6:16],lmer_res@u)
  # print(lm(vecbc[6:16]~lmer_res@u))
  
  # Equ (7) Eilers 2018: compute H
  H <- cbind(X,Z)%*%Q%*%rbind(t(X),t(Z))%*%Rinv
  cat("Effective dimension : ", sum(diag(H)),"\n")
  
  # Comparaison yhat lmer et Hy
  yhat2=H%*%y
  # print(lm(yhat2~yhat))
  plot(yhat,yhat2,main="Y predit par lmer et par Hy")
  plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")
  
  
  # dimensions effectives sur base de K
  K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
  # dim(K)
  K_diag <- diag(K)
  # sum(K_diag)
  
  K11=K[1:3,1:3] # Tube+time+intercept
  K22=K[4:39,4:39] # # Volunteer:Sampling
  K33=K[40:51,40:51] # # Volunteer
  cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
  cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
  cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
  cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")
  
  
  ED_Tube_time[i] <- sum(diag(K11))
  ED_Samp[i] <- sum(diag(K22))
  ED_Vol[i] <- sum(diag(K33))
}

####################
ED <- rbind(ED_Tube_time, ED_Samp, ED_Vol)
colnames(ED) <- names(MM_full$merMod_obj)
colSums(ED)
ED_Tube <- rep(1, ncol(ED))
ED_Time <- rep(1, ncol(ED))
pdf("/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments/ED_metabo_unbal.pdf",
    width = 12, height = 7, pointsize = 16)
par(mar=c(5,4,4,4), cex=1.2)
plot(1:ncol(ED),ED_Tube, ylim=range(ED), xaxt="n",
     type="b", col="red", xlab="PC", ylab = "ED",lty=4, pch=3)
lines(1:ncol(ED),ED_Time, type="b", col="ForestGreen", lty=2, pch=1)
lines(1:ncol(ED),ED_Samp, type="b", col="mediumvioletred", lty=2, pch=4)
lines(1:ncol(ED),ED_Vol, type="b", col="blue", lty=1, pch=2)
legend("top",legend = c("Tube","Time", "Sampling", "Volunteer"),
       col = c("red", "ForestGreen", "mediumvioletred", "blue"), lty=c(4,2,2,1), ncol = 4,
       bty = "n", xpd=TRUE, inset = c(0,-0.2), pch = c(3,1,4,2))
axis(side = 1, at = 1:ncol(ED), labels = 1:ncol(ED))
# mtext("EDtot:", side=1, at = 0, line=4)
# mtext(format(round(colSums(ED),1), nsmall=1), side=1, at = 1:ncol(ED), line=4, cex=0.8)
dev.off()

ED_metabo <- ED
save(ED_metabo, file=file.path(out_path, "ED_metabo_unbal.RData"))
write.csv(ED, file.path(out_path,"ED_metabo_unbal.csv"))

####################
rm(list=ls())
out_path <- "/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments"
load(file.path(out_path,"res_parlmer_metabo_unbalanced.RData"))
load(file=file.path(out_path, "ED_metabo_unbal.RData"))
# require(lme4)
n <- nrow(outcomes)
ED_Tube <- rep(1, ncol(ED_metabo))
ED_Time <- rep(1, ncol(ED_metabo))
ED_Samp <- ED_metabo["ED_Samp", ]
ED_Vol <- ED_metabo["ED_Vol", ]

pdf("/Users/manon/Documents/Doctorat/Manuscript_these/manuscript/Manuscript_review/Paul_comments/ED_metabo_unbal.pdf",
    width = 12, height = 7, pointsize = 13)
par(mar=c(4,4,1,8.5), cex=1.2)
plot(1:ncol(ED_metabo),ED_Tube, ylim=c(0,40), xaxt="n",
     type="p", col="red", xlab="PC", ylab = "ED",pch=3)
lines(1:ncol(ED_metabo),ED_Time, type="l", col="red", lty=2, lwd=2)

lines(1:ncol(ED_metabo),ED_Time, type="p", col="ForestGreen", pch=1, lwd=2)
lines(1:ncol(ED_metabo),ED_Time, type="l", col="ForestGreen", lty=3, lwd=2)


lines(1:ncol(ED_metabo),ED_Samp, type="p", col="mediumvioletred", pch=4, lwd=2)
max_ed <- length(unique(design$Sampling))*length(unique(design$Volunteer))
lines(1:ncol(ED_metabo),rep(max_ed, ncol(ED_metabo)), 
      type="l", col="mediumvioletred", lty=1)


lines(1:ncol(ED_metabo),ED_Vol, type="p", col="blue",pch=2)
max_ed <- length(unique(design$Volunteer))
lines(1:ncol(ED_metabo),rep(max_ed, ncol(ED_metabo)), 
      type="l", col="blue", lty=4)


legend("topright",legend = c("Tube","Time", "Sampling", "Volunteer",
                             "Max Tube","Max Time", "Max Sampling", "Max Volunteer"),
       col = c("red", "ForestGreen", "mediumvioletred", "blue"), 
       lty=c(NA,NA,NA,NA,2,3,1,4), 
       bty = "n", xpd=TRUE, inset = c(-0.25,0), 
       pch = c(3,1,4,2,NA,NA,NA,NA), lwd=2)
axis(side = 1, at = 1:ncol(ED_metabo), labels = 1:ncol(ED_metabo))
dev.off()

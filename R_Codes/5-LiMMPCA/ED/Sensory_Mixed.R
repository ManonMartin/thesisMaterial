    
# Sensory data

load(file.path("ED/res_parlmer_sensory.RData"))
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



#########################################
# For all the responses  ##############
#########################################


ED_candies <- c()
ED_cj <- c()
ED_judges <- c()

for (i in 1:length(res.parlmer$merMod_obj)){
  
  lmer_res <- res.parlmer$merMod_obj[[i]]
  y <- outcomes[,i]
  
  # yhat 
  yhat <- predict(lmer_res)
  
  sigma2r <- varcor_resid[i]
  sigma2g <- varcor_random_full["Judges" ,i]
  sigma2gg <- varcor_random_full["CandiesJudges",i]
  
  
  q = ncol(Z)
  G <- diag(q)*0
  G[1:44,1:44] = diag(44)*sigma2gg
  G[45:55,45:55] = diag(11)*sigma2g
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
  # cat("Effective dimension : ", sum(diag(H)),"\n")
  
  # Comparaison yhat lmer et Hy
  yhat2=H%*%y
  # print(lm(yhat2~yhat))
  # plot(yhat,yhat2,main="Y predit par lmer et par Hy")
  # plot(yhat-yhat2,main="Y predit par lmer Moins par Hy")
  
  
  # dimensions effectives sur base de K
  K=Q%*%cbind(rbind(XRX,ZRX),rbind(XRZ,ZRZ))
  # dim(K)
  K_diag <- diag(K)
  # sum(K_diag)
  
  K11=K[1:5,1:5] # C
  K22=K[6:49,6:49] # CJ
  K33=K[50:60,50:60] # J
  # cat("\n Effective dimension K11 : ", sum(diag(K11)),"\n")
  # cat("\n Effective dimension K22 : ", sum(diag(K22)),"\n")
  # cat("\n Effective dimension K22 : ", sum(diag(K33)),"\n")
  # cat("\n Effective dimension (totale) : ", sum(K_diag),"\n")

  ED_candies[i] <- sum(diag(K11))
  ED_cj[i] <- sum(diag(K22))
  ED_judges[i] <- sum(diag(K33))
}

ED_candies <- ED_candies-1
ED <- rbind(ED_candies, ED_cj, ED_judges)
colnames(ED) <- names(MM_full$merMod_obj)
colSums(ED)
max_ED <- 5 + 11 + 44



ED_sensory <- ED
save(ED_sensory, file=file.path(out_path, "ED_sensory.RData"))
write.csv(ED_sensory, file.path(out_path,"ED_candies.csv"))

n <- nrow(outcomes)
ED_candies = ED_sensory["ED_candies", ]
ED_cj = ED_sensory["ED_cj", ]
ED_judges = ED_sensory["ED_judges", ]

# pdf("ED_candies.pdf",
#     width = 8, height = 5, pointsize = 10)
par(mar=c(4,4,1,8), cex=1.1)
plot(1:8,ED_candies, ylim = c(0,45),
     type="p", col="red",  xlab="PC", ylab = "ED", pch=1, lwd=1.5)
lines(1:8,ED_candies, type="l", col="red", lty=4, lwd=1.5)

lines(1:8,ED_cj, type="p", col="blue", pch=2, lwd=1.5)
max_ed <- (length(unique(designInter$Candies))-1)*length(unique(designInter$Judges))
lines(1:8,rep(max_ed,8), type="l", col="blue", lty=2, lwd=1.5)

lines(1:8,ED_judges, type="p", col="forestgreen", pch=4, lwd=1.5)
max_ed <-length(unique(designInter$Judges))
lines(1:8,rep(max_ed,8), type="l", col="forestgreen", lty=1, lwd=1.5)


legend("topright",legend = c("ED Candies", "ED Judges", "ED CJ", 
                             "Max Candies",
                             "Max Judges", "Max CJ"),
       col = c("red", "forestgreen", "blue","red", "forestgreen", "blue"), 
       lty=c(NA,NA,NA,4,2,1), pch=c(1,2,4, NA, NA, NA),
       bty = "n", xpd=TRUE, inset = c(-0.25,0), lwd=1.5)
# dev.off()



setwd("/Users/manon/Documents/Doctorat/Projects/ASCA_modeles_mixtes/Data/SensoryDataNaes")
rm(list=ls())
X <- read.csv(file = "X.csv", header = FALSE)
Xa <- read.csv(file = "Xa.csv", header = FALSE)
Y <- read.csv(file = "Y.csv", header = FALSE)


colnames(Xa) <- c("Assessors", "Candies")


# PCA
require(MBXUCL)
colnames(Y) <- c("Transp","Acid" ,"Sweet" , "Raspb.", "Sugar", "Bites","Hard" ,"Elastic" ,"Sticky")

res <- SVDforPCA(Y)
res$var

plot(-1*res$loadings[,1], -1*res$loadings[,2], xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.6))
text(-1*res$loadings[,1], -1*res$loadings[,2], rownames(res$loadings), pos=4)


save(X, Xa, Y, file="Candies.RData")

  

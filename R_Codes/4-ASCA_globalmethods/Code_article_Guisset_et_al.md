---
title: "Code article Guisset et Al - PARAFASCA - ACOMDIM - AMOPLS"
authors: "Manon Martin - Severine Guisset - Bernadette Govaerts"
date: 'October 16, 2019,21:37'
output:
  html_document:
    keep_md: yes
    smart: FALSE
    code_folding: hide
    collapsed: yes
    fig_caption: yes
    fig_height: 6
    fig_width: 9
    highlight: tango
    number_sections: yes
    theme: united
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '3'
  word_document: default
editor_options: 
  chunk_output_type: console
references:
  - id: Guisset2019
    title: Comparison of PARAFASCA, AComDim, and AMOPLS approaches in the multivariate GLM modelling of multi-factorial designs
    author:
    - family: Guisset
      given: S.
    - family: Martin
      given: M.
    - family: Govaerts
      given: B.
    container-title: Chemometrics and Intelligent Laboratory Systems
    volume: 184
    URL: 'https://doi.org/10.1016/j.chemolab.2018.11.006'
    DOI: 10.1016/j.chemolab.2018.11.006
    page: 44 - 63
    type: article-journal
    issued:
      year: 2019
      month: 1
---




This code generates the useful results for the article @Guisset2019 (Cf. the Reference section below).

# Load libraries and functions 

Different librairies are needed to run this code.





Specify where are located the functions and data


```r
# Choix du répertoire de travail

# directoryMain <- file.path('/Users/manon/Desktop/CodeArticleGuissetetAl')
# directoryMain <-
# file.path('~/Dropbox/PartageDiversProfessionel/Partage_Package_ASCA/CodeGuissetetAlMinimum')

directoryData <- "../../Datasets/3_citrate_hippurate_urine"

# lecture des routines R
DirectoryFun <- file.path("LIB")

# install manually kopls (downloaded from mac or linux
kopls.path <- "Lib_KOPLS_for_mac/kopls_1.1.2.tar.gz"
# windows kopls.path <- 'Lib_KOPLS_for_mac/kopls_1.1.2.tar.gz'
# install.packages(kopls.path, repos = NULL, type = 'source')


source(file.path(DirectoryFun, "Decomposition_functions.R"))
source(file.path(DirectoryFun, "AComDim_functions.R"))
source(file.path(DirectoryFun, "PARAFAC_array.R"))
source(file.path(DirectoryFun, "AMOPLS_functions.R"))
source(file.path(DirectoryFun, "permutationTest.R"))
source(file.path(DirectoryFun, "matrixDecomposition.R"))


# output path
directoryOut <- "Outputs"

nperm <- 1000
```

# Data

## Uploading data

`D_UCH.RData` contains the complete design of the urine database and `X_UCH.RData` the complete spectral matrix. `matrixDecomposition.R` decomposes the spectral matrix according to the design


```r
load(file.path(directoryData, "UCH.Rdata"))
```


```r
pander("outcomes")
```

outcomes

```r
str(head(outcomes))
```

```
##  num [1:6, 1:600] 0.0312 0.0581 0.027 0.0341 0.0406 ...
##  - attr(*, "dimnames")=List of 2
##   ..$   : chr [1:6] "M2C00D2R1" "M2C00D2R2" "M2C02D2R1" "M2C02D2R2" ...
##   ..$ X1: chr [1:600] "9.9917004" "9.9753204" "9.9590624" "9.9427436" ...
```

```r
pander("design")
```

design

```r
str(design)
```

```
## 'data.frame':	34 obs. of  5 variables:
##  $ Hippurate: Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 2 2 2 2 ...
##  $ Citrate  : Factor w/ 3 levels "0","2","4": 1 1 2 2 3 3 1 1 2 2 ...
##  $ Dilution : Factor w/ 1 level "diluted": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Day      : Factor w/ 2 levels "2","3": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Time     : Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...
```

```r
pander("formula")
```

formula

```r
print(formula)
```

```
## outcomes ~ Hippurate + Citrate + Time + Hippurate * Citrate + 
##     Time * Hippurate + Time * Citrate + Hippurate * Citrate * 
##     Time
```


```r
attach(design)

n <- dim(outcomes)[1]
m <- dim(outcomes)[2]

pander("table(Hippurate, Citrate, Time)")
```

table(Hippurate, Citrate, Time)

```r
table(Hippurate, Citrate, Time)
```

```
## , , Time = 1
## 
##          Citrate
## Hippurate 0 2 4
##         0 2 1 1
##         1 2 2 2
##         2 2 2 2
## 
## , , Time = 2
## 
##          Citrate
## Hippurate 0 2 4
##         0 2 2 2
##         1 2 2 2
##         2 2 2 2
```

```r
Variables = as.numeric(dimnames(outcomes)[[2]])  # Vector with variable names 
Samples = dimnames(outcomes)[[1]]  # Vector with sample names
```

 

## Figure 1 spectrum

```r
### Spectral profile
plot(outcomes["M2C24D2R1", ], type = "l", xaxt = "n", ylab = "Intensity", xlab = "ppm", 
    cex.lab = 1, mgp = c(2.3, 1, 0))
axis(side = 1, at = seq(1, m, 55), labels = round(as.numeric(colnames(outcomes))[seq(1, 
    m, 55)], 1))
title(expression(paste("Example of ", a^{
    1
}, "H NMR spectral profile")), line = 1)
```

<img src="Code_article_Guisset_et_al_files/figure-html/Spectralprofile-1.png" width="800px" />


# Define a nice Scatterplot Matrix function 

```r
#------- define the Scatterplot Matrix function
# if spotPoints == TRUE, will circle specific outlying observations
pairsBG = function(MatScores, titre = "Scores", spotPoints = FALSE, ...) {
    panelup = function(x, y) {
        col <- rep("blue", length(Citrate))
        col[Citrate == "2"] <- "forestgreen"
        col[Citrate == "4"] <- "red"
        pch <- rep(4, length(Hippurate))
        pch[Hippurate == "1"] <- 16
        pch[Hippurate == "2"] <- 2
        
        points(x, y, col = col, pch = pch)
        
        if (spotPoints) {
            text(x, y, textlabel, pos = c(3), col = "darkturquoise")
            points(x, y, pch = pch2, col = "darkturquoise")
        }
        
    }
    paneldown = function(x, y) {
        
        pch = rep(1, length(Time))
        pch[Time == 2] <- 3
        
        col <- rep("orange", length(Time))
        col[Time == 2] <- "black"
        
        points(x, y, col = col, pch = pch)
    }
    pairs(MatScores, main = titre, upper.panel = panelup, lower.panel = paneldown, 
        gap = 0.3, ...)
}


#------- define a function to draw the legends next to the scatterplot

legendsScatterMatrix <- function() {
    legend("topright", title = expression(bold("Hippurate")), legend = c("0", 
        "1", "2"), bty = "n", col = c(1, 1), pch = c(4, 16, 2), inset = c(0.11, 
        0.1), cex = 0.7)
    
    legend("topright", title = expression(bold("Citrate")), legend = c("0", 
        "2", "4"), bty = "n", seg.len = 0.4, pch = 15, col = c("blue", "forestgreen", 
        "red"), inset = c(0.11, 0.4), cex = 0.7)
    
    legend("topright", title = expression(bold("Time")), legend = c("1", "2"), 
        bty = "n", seg.len = 0.4, pch = c(1, 3), col = c("orange", "black"), 
        inset = c(0.11, 0.7), cex = 0.7)
}
```

# PCA on the spectral matrix 

PCA decomposition is first applied on the spectral matrix `outcomes` by applying the function `SVDforPCA` from the package MBXUCL. Scores plots and loadings plots are obtained with `DrawScores` and `DrawLoadings` functions respectively. 

## PCA and % var explained


```r
pcaOutcomes = SVDforPCA(outcomes)
eig.res = rbind(pcaOutcomes$var, pcaOutcomes$var * 100/sum(pcaOutcomes$var), 
    pcaOutcomes$cumvar)[, 1:6]
rownames(eig.res) = c("Variances", "Prop Var", "Cum Eigen Values")
pander(eig.res)
```


-----------------------------------------------------------------------
        &nbsp;           PC1     PC2     PC3     PC4     PC5     PC6   
---------------------- ------- ------- ------- ------- ------- --------
    **Variances**       45.25   24.83   16.75   6.716   3.412   0.9151 

     **Prop Var**       45.25   24.83   16.75   6.716   3.412   0.9151 

 **Cum Eigen Values**   45.25   70.08   86.82   93.54   96.95   97.87  
-----------------------------------------------------------------------

## scores plot


```r
col <- c(`0` = "blue", `2` = "forestgreen", `4` = "red")
pch <- c(`0` = 4, `1` = 1, `2` = 2)

responsematrix_Scores <- DrawScores(obj = pcaOutcomes, type.obj = "PCA", drawNames = F, 
    createWindow = FALSE, main = "Reponse matrix scores plot", color = Citrate, 
    pch = Hippurate, size = 2, cex.lab = 3, axes = c(1, 2), xlab = NULL, ylab = NULL, 
    drawEllipses = FALSE, typeEl = "norm", levelEl = 0.9, drawPolygon = FALSE, 
    noLegend = FALSE, legend_color_manual = col, legend_shape_manual = pch) + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50")) + 
    geom_text(aes(label = Time), hjust = c(0), vjust = -0.9, col = "black", 
        cex = 3)

responsematrix_Scores
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-8-1.png" width="800px" />

```r
Time <- paste0("time ", Time)

DrawScores(obj = pcaOutcomes, type.obj = "PCA", drawNames = F, createWindow = FALSE, 
    main = "Reponse matrix scores plot", color = Citrate, pch = Time, size = 2, 
    cex.lab = 3, axes = c(1, 2), xlab = NULL, ylab = NULL, drawEllipses = FALSE, 
    typeEl = "norm", levelEl = 0.9, drawPolygon = FALSE, noLegend = FALSE, legend_color_manual = col) + 
    scale_shape_manual(name = "Time", breaks = c("time 1", "time 2"), values = c("1", 
        "2"))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-8-2.png" width="800px" />

```r
Time <- design$Time

DrawScores(pcaOutcomes, type.obj = "PCA", drawNames = TRUE, createWindow = F, 
    main = "Reponse matrix score plot", color = Citrate, pch = Hippurate, axes = c(1, 
        2), size = 2.5)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-8-3.png" width="800px" />

```r
DrawScores(pcaOutcomes, type.obj = "PCA", drawNames = TRUE, createWindow = F, 
    main = "Reponse matrix score plot", color = Citrate, pch = Hippurate, axes = c(3, 
        4), size = 2.5)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-8-4.png" width="800px" />

```r
DrawScores(pcaOutcomes, type.obj = "PCA", drawNames = TRUE, createWindow = F, 
    main = "Reponse matrix score plot", color = Citrate, pch = Hippurate, axes = c(5, 
        6), size = 2.5)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-8-5.png" width="800px" />


## scatter plot des scores

```r
par(xpd = TRUE)
pairsBG(pcaOutcomes$scores[, 1:6], "Scores PCA")
legendsScatterMatrix()
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-9-1.png" width="800px" />

# GLM decomposition

## Model definition


```r
ModelTerms_sansE <- attr(terms(formula), "term.labels")  # without residuals
ModelTerms <- c(ModelTerms_sansE, "residuals")  # with residuals

ModelTerms_abbrev <- ModelTerms  # Abbreviated model terms
index <- gregexpr(":", ModelTerms)
for (i in 1:length(ModelTerms)) {
    if (index[[i]][1] != -1) {
        ModelTerms_abbrev[i] <- substr(ModelTerms[i], 1, 1)
        index2 <- index[[i]]
        for (k in 1:length(index2)) {
            ModelTerms_abbrev[i] <- paste(ModelTerms_abbrev[i], substr(ModelTerms[i], 
                index2[k] + 1, index2[k] + 1), sep = "x")
        }
    }
}

ModelTerms_abbrev[length(ModelTerms_abbrev)] = "Residuals"

p <- length(ModelTerms)  # N model terms
```



## GLM decomposition


```r
# Décomposition de Y en matrice des effets par GLM

resGLM <- matrixDecomposition(formula, outcomes, design)
modelMatrix <- sapply(resGLM$modelMatrixByEffect, function(x) x)
nparam <- sum(sapply(modelMatrix, function(x) dim(x)[2])) - 1


# construction de la liste des matrices des effets purs
EffectMatGLM <- resGLM$effectMatrices[-1]  # minus intercept
res <- vector(mode = "list")
res[[1]] <- resGLM$residuals
EffectMatGLM <- c(EffectMatGLM, residuals = res)  # plus residuals
pander("names(EffectMatGLM)")
```

names(EffectMatGLM)

```r
names(EffectMatGLM)
```

```
## [1] "Hippurate"              "Citrate"               
## [3] "Time"                   "Hippurate:Citrate"     
## [5] "Hippurate:Time"         "Citrate:Time"          
## [7] "Hippurate:Citrate:Time" "residuals"
```

```r
# Calcul des matrices augmentées
EffectMatGLMAug <- resGLM$effectMatrices[-1]  # effectMatrices minus intercept
EffectMatGLMAug <- lapply(EffectMatGLMAug, function(x) x + resGLM$residuals)
EffectMatGLMAug <- c(EffectMatGLMAug, residuals = res)  # plus residuals
pander("names(EffectMatGLMAug)")
```

names(EffectMatGLMAug)

```r
names(EffectMatGLMAug)
```

```
## [1] "Hippurate"              "Citrate"               
## [3] "Time"                   "Hippurate:Citrate"     
## [5] "Hippurate:Time"         "Citrate:Time"          
## [7] "Hippurate:Citrate:Time" "residuals"
```

```r
Xmat <- resGLM$modelMatrix[, -1]  # modelMatrix minus intercept
```


## Variation percentages

How to get effects matrices for all the effects.
The sum is equal to 100% for balanced designs


```r
VariationPercentages <- unlist(resGLM$variationPercentages)
pander(round(VariationPercentages, 2))
```


------------------------------------------------------------------
 Hippurate   Citrate   Time    Hippurate:Citrate   Hippurate:Time 
----------- --------- ------- ------------------- ----------------
   39.31      29.91    16.24         1.54               6.23      
------------------------------------------------------------------

Table: Table continues below

 
---------------------------------------------------
 Citrate:Time   Hippurate:Citrate:Time   residuals 
-------------- ------------------------ -----------
     0.54                1.68               4.3    
---------------------------------------------------

```r
barplot(t(unlist(resGLM$variationPercentages)), las = 2)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-12-1.png" width="800px" />

```r
sum(VariationPercentages)
```

```
## [1] 99.74193
```

## Permutation tests
### Result permut

```r
load(file.path(directoryOut, "PermutationResASCA.RData"))

names(pval) <- names(ASCAFSTAT) <- names(ASCAFSTATperm) <- ModelTerms_sansE

pander(pval)
```


--------------------------------------------------------------------------------
 Hippurate   Citrate   Time   Hippurate:Citrate   Hippurate:Time   Citrate:Time 
----------- --------- ------ ------------------- ---------------- --------------
     0          0       0           0.146               0             0.448     
--------------------------------------------------------------------------------

Table: Table continues below

 
------------------------
 Hippurate:Citrate:Time 
------------------------
         0.104          
------------------------

```r
for (i in ModelTerms_sansE) {
    hist(ASCAFSTATperm[[i]], xlim = c(0, max(ASCAFSTAT[i], ASCAFSTATperm[[i]])), 
        breaks = 20, main = i, freq = FALSE)
    points(ASCAFSTAT[i], 0, col = "red", pch = 19)
    rangex = range(ASCAFSTAT[i], ASCAFSTATperm[[i]])
    xf = seq(0, rangex[2], by = (rangex[2] - rangex[1])/1000)
    yf = df(xf, n - 1, n - 1)
    lines(xf, yf)
}
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-1.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-2.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-3.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-4.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-5.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-6.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-13-7.png" width="50%" />


# ASCA-GLM 

## scores/loadings


```r
listgraphs <- list()
varASCA <- list()

for (i in 1:length(ModelTerms)) {
    # PCA on the non augmented effect matrices
    ascaSVD = SVDforPCA(EffectMatGLM[[i]])
    ascaSVD$scores = round(ascaSVD$scores, 5)
    varASCA[[i]] <- ascaSVD$var
    # Scores plots with Citrate/hippurate colors
    listgraphs[[paste0(ModelTerms[i], "-CH")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = TRUE, createWindow = F, main = paste0(ModelTerms[i], "-CH", 
            " 
                                                                     scores plot - ASCA-GLM"), 
        color = as.factor(Citrate), pch = as.factor(Hippurate), axes = c(1, 
            2), size = 2.5)
    # Scores plots with Time/Day colors
    listgraphs[[paste0(ModelTerms[i], "-DM")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = TRUE, createWindow = F, main = paste0(ModelTerms[i], "-DM", 
            " scores plot - ASCA-GLM"), color = as.factor(Day), pch = as.factor(Time), 
        axes = c(1, 2), size = 2.5)
    
    # Loadings plots
    listgraphs[[paste0(ModelTerms[i], "-Loadings")]] <- DrawLoadings(ascaSVD, 
        type.obj = "PCA", createWindow = F, main = paste0(ModelTerms[i], " loadings plot - ASCA-GLM"), 
        axes = c(1:2), loadingstype = "s", num.stacked = 2, xlab = NULL, ylab = NULL, 
        ang = "0", xaxis_type = "numerical", nxaxis = 10)
}

listgraphs
```

```
## $`Hippurate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-1.png" width="800px" />

```
## 
## $`Hippurate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-2.png" width="800px" />

```
## 
## $`Hippurate-Loadings`
## $`Hippurate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-3.png" width="800px" />

```
## 
## 
## $`Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-4.png" width="800px" />

```
## 
## $`Citrate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-5.png" width="800px" />

```
## 
## $`Citrate-Loadings`
## $`Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-6.png" width="800px" />

```
## 
## 
## $`Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-7.png" width="800px" />

```
## 
## $`Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-8.png" width="800px" />

```
## 
## $`Time-Loadings`
## $`Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-9.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-10.png" width="800px" />

```
## 
## $`Hippurate:Citrate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-11.png" width="800px" />

```
## 
## $`Hippurate:Citrate-Loadings`
## $`Hippurate:Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-12.png" width="800px" />

```
## 
## 
## $`Hippurate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-13.png" width="800px" />

```
## 
## $`Hippurate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-14.png" width="800px" />

```
## 
## $`Hippurate:Time-Loadings`
## $`Hippurate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-15.png" width="800px" />

```
## 
## 
## $`Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-16.png" width="800px" />

```
## 
## $`Citrate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-17.png" width="800px" />

```
## 
## $`Citrate:Time-Loadings`
## $`Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-18.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-19.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-20.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-Loadings`
## $`Hippurate:Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-21.png" width="800px" />

```
## 
## 
## $`residuals-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-22.png" width="800px" />

```
## 
## $`residuals-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-23.png" width="800px" />

```
## 
## $`residuals-Loadings`
## $`residuals-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-14-24.png" width="800px" />


## Variance decomposition for PC1 and PC2

```r
names(varASCA) <- ModelTerms

pander("ASCA variance decomposition for PC1 and 2")
```

ASCA variance decomposition for PC1 and 2

```r
# PC1
round(sapply(varASCA, function(x) x[1]), 2)
```

```
##              Hippurate.PC1                Citrate.PC1 
##                      97.71                      98.22 
##                   Time.PC1      Hippurate:Citrate.PC1 
##                     100.00                      44.01 
##         Hippurate:Time.PC1           Citrate:Time.PC1 
##                      93.92                      90.76 
## Hippurate:Citrate:Time.PC1              residuals.PC1 
##                      47.23                      48.54
```

```r
# PC2
round(sapply(varASCA, function(x) x[2]), 2)
```

```
##              Hippurate.PC2                Citrate.PC2 
##                       2.29                       1.78 
##                   Time.PC2      Hippurate:Citrate.PC2 
##                       0.00                      38.51 
##         Hippurate:Time.PC2           Citrate:Time.PC2 
##                       6.08                       9.24 
## Hippurate:Citrate:Time.PC2              residuals.PC2 
##                      27.49                      16.90
```

```r
varASCAPC12 <- rbind(sapply(varASCA, function(x) x[1]), sapply(varASCA, function(x) x[2]))


pander("sum over PC1 and PC2")
```

sum over PC1 and PC2

```r
round(apply(varASCAPC12, 2, sum), 2)
```

```
##              Hippurate.PC1                Citrate.PC1 
##                     100.00                     100.00 
##                   Time.PC1      Hippurate:Citrate.PC1 
##                     100.00                      82.52 
##         Hippurate:Time.PC1           Citrate:Time.PC1 
##                     100.00                     100.00 
## Hippurate:Citrate:Time.PC1              residuals.PC1 
##                      74.72                      65.44
```


## Hippurate scores and loadings

```r
# Hippurate scores and loadings
ascaSVD <- SVDforPCA(EffectMatGLM$Hippurate)
ascaSVD$scores <- round(ascaSVD$scores, 5)

ASCAScoresHippurate <- DrawScores(ascaSVD, type.obj = "PCA", drawNames = F, 
    createWindow = F, main = "ASCA scores plot - Hippurate effect", axes = c(1, 
        2), size = 2.5, pch = Hippurate, color = Hippurate, noLegend = TRUE) + 
    scale_colour_manual(name = "Hippurate:", breaks = c("0", "1", "2"), values = c("blue", 
        "forestgreen", "red")) + scale_shape_manual(name = "Hippurate:", breaks = c("0", 
    "1", "2"), values = c(4, 1, 2))


ASCAScoresHippurate
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-16-1.png" width="800px" />

```r
ASCALoadingHippurate <- DrawLoadings(ascaSVD, type.obj = "PCA", createWindow = F, 
    main = "ASCA loading plot - Hippurate effect", axes = c(1), loadingstype = "s", 
    num.stacked = 2, xlab = "ppm", ylab = "", ang = "0", xaxis_type = "numerical", 
    nxaxis = 10)

ASCALoadingHippurate
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-16-2.png" width="800px" />


# ASCA-E-GLM 

## Overview of scores and loadings

PCA decomposition is applied on each non-augmented spectral matrix. Residuals are then projected on the dimensions 2 first dimension for each effect et scores plots prepared. The loadings are just those of the original PCA. 


```r
listgraphs <- list()
for (i in 1:(length(ModelTerms))) {
    # PCA on the non-augmented effect matrices
    ascaSVD <- SVDforPCA(EffectMatGLM[[i]])
    # scores 1 and 2 correction by adding a projection of the residuals matrix
    ascaSVD$scores[, 1:2] <- (EffectMatGLM[[i]] + EffectMatGLM[[p]]) %*% ascaSVD$loadings[, 
        1:2]
    # Scores plot with Citrate/hippurate colors
    listgraphs[[paste0(ModelTerms[i], "-CH")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = F, createWindow = F, main = paste0(ModelTerms[i], "-CH", 
            " score plot"), color = Citrate, pch = Hippurate, axes = c(1, 2), 
        size = 2.5)
    # Scores plot with Media/Time colors
    listgraphs[[paste0(ModelTerms[i], "-RM")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = F, createWindow = F, main = paste0(ModelTerms[i], "-DM", 
            " score plot"), color = Time, pch = Dilution, axes = c(1, 2), size = 2.5)
    # Loadings plot
    listgraphs[[paste0(ModelTerms[i], "-Loadings")]] <- DrawLoadings(ascaSVD, 
        type.obj = "PCA", createWindow = F, main = paste0(ModelTerms[i], " loading plot"), 
        axes = c(1:2), loadingstype = "s", num.stacked = 2, xlab = NULL, ylab = NULL, 
        ang = "0", xaxis_type = "numerical", nxaxis = 10)
}
listgraphs
```

```
## $`Hippurate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-1.png" width="800px" />

```
## 
## $`Hippurate-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-2.png" width="800px" />

```
## 
## $`Hippurate-Loadings`
## $`Hippurate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-3.png" width="800px" />

```
## 
## 
## $`Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-4.png" width="800px" />

```
## 
## $`Citrate-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-5.png" width="800px" />

```
## 
## $`Citrate-Loadings`
## $`Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-6.png" width="800px" />

```
## 
## 
## $`Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-7.png" width="800px" />

```
## 
## $`Time-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-8.png" width="800px" />

```
## 
## $`Time-Loadings`
## $`Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-9.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-10.png" width="800px" />

```
## 
## $`Hippurate:Citrate-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-11.png" width="800px" />

```
## 
## $`Hippurate:Citrate-Loadings`
## $`Hippurate:Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-12.png" width="800px" />

```
## 
## 
## $`Hippurate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-13.png" width="800px" />

```
## 
## $`Hippurate:Time-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-14.png" width="800px" />

```
## 
## $`Hippurate:Time-Loadings`
## $`Hippurate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-15.png" width="800px" />

```
## 
## 
## $`Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-16.png" width="800px" />

```
## 
## $`Citrate:Time-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-17.png" width="800px" />

```
## 
## $`Citrate:Time-Loadings`
## $`Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-18.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-19.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-20.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-Loadings`
## $`Hippurate:Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-21.png" width="800px" />

```
## 
## 
## $`residuals-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-22.png" width="800px" />

```
## 
## $`residuals-RM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-23.png" width="800px" />

```
## 
## $`residuals-Loadings`
## $`residuals-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-17-24.png" width="800px" />



## PCA on the effect matrices


```r
# PC var for each effect matrix PCA
pcVar <- vector("list", length = length(ModelTerms))
Y_Loadings <- vector("list", length = length(ModelTerms))
ascaSVD_scores <- vector("list", length = (length(ModelTerms) + 1))
names(pcVar) <- names(Y_Loadings) <- ModelTerms

names(ascaSVD_scores) <- c(ModelTerms[-p], "residuals1", "residuals2")

for (i in 1:(p - 1)) {
    ascaSVD = SVDforPCA(EffectMatGLM[[i]])
    pcVar[[ModelTerms[i]]] <- ascaSVD$var/100
    Y_Loadings[[ModelTerms[i]]] <- ascaSVD$loadings
    ascaSVD$scores[, 1:2] = (EffectMatGLM[[i]] + EffectMatGLM$residuals) %*% 
        ascaSVD$loadings[, 1:2]
    ascaSVD_scores[[ModelTerms[i]]] <- ascaSVD$scores[, 1]
}

# Error
ascaSVD = SVDforPCA(EffectMatGLM$residuals)
ascaSVD_error <- ascaSVD
pcVar[["residuals"]] <- ascaSVD$var/100
ascaSVD$scores = round(ascaSVD$scores, 5)
ascaSVD_scores[["residuals1"]] <- ascaSVD$scores[, 1]
ascaSVD_scores[["residuals2"]] <- ascaSVD$scores[, 2]
```

## Saliences plot


```r
saliences <- matrix(0, nrow = p, ncol = p + 1)
colnames(saliences) <- paste(c(ModelTerms[-p], "residuals1", "residuals2"), 
    "PC1")

colnames(saliences)[p + 1] <- "residuals2 PC2"

rownames(saliences) <- ModelTerms
for (i in 1:(p - 1)) {
    saliences[ModelTerms[i], colnames(saliences)[i]] <- pcVar[[i]][1]
}

saliences["residuals", "residuals1 PC1"] <- pcVar$residuals[1]
saliences["residuals", "residuals2 PC2"] <- pcVar$residuals[2]


# plot --------------------------------------------------------------------

ylim <- c(0, 1)
angle1 <- rep(c(45, 45, 135), length.out = p)
angle2 <- rep(c(45, 135, 135), length.out = p)
density1 <- seq(5, 35, length.out = p)
density2 <- seq(5, 35, length.out = p)
col <- rainbow((p))


par(xpd = TRUE, mar = c(5, 4, 4, 4))
x <- barplot(saliences, beside = TRUE, ylim = ylim, xaxt = "n", col = col, angle = angle1, 
    density = density1, main = "Effects explained by ASCA-E components")
barplot(saliences, add = TRUE, beside = TRUE, ylim = ylim, xaxt = "n", col = col, 
    angle = angle2, density = density2)

labs <- as.vector(colnames(saliences))
text(cex = 1, x = colMeans(x) - 0.25, y = -0.1, labs, xpd = TRUE, srt = 45)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-19-1.png" width="800px" />




## ASCA variance decomposition

```r
pcexp_pereffect_Resid <- (ascaSVD_error$var[1:2]/100) * VariationPercentages["residuals"]

pcexp_pereffect_Effects <- VariationPercentages * round(sapply(varASCA, function(x) x[1]), 
    2)/100

pcexp_pereffect_perPC <- c(pcexp_pereffect_Effects[-p], pcexp_pereffect_Resid)

pander("pcexp_pereffect")
```

pcexp_pereffect

```r
pander(pcexp_pereffect_perPC)
```


------------------------------------------------------------------
 Hippurate   Citrate   Time    Hippurate:Citrate   Hippurate:Time 
----------- --------- ------- ------------------- ----------------
   38.41      29.37    16.24        0.6789              5.85      
------------------------------------------------------------------

Table: Table continues below

 
--------------------------------------------------------
 Citrate:Time   Hippurate:Citrate:Time    PC1     PC2   
-------------- ------------------------ ------- --------
    0.4889              0.7954           2.086   0.7263 
--------------------------------------------------------



## Scores plots


```r
ascaSVD <- SVDforPCA(EffectMatGLM$Hippurate)
ascaSVD$scores[, 1:2] <- (EffectMatGLM$Hippurate + EffectMatGLM$residuals) %*% 
    ascaSVD$loadings[, 1:2]

col <- c(`0` = "blue", `1` = "forestgreen", `2` = "red")

ASCAEScoresHippurate <- DrawScores(ascaSVD, type.obj = "PCA", drawNames = F, 
    createWindow = F, main = "ASCA-E scores plot - Hippurate effect", axes = c(1, 
        2), size = 2.5, noLegend = FALSE, pch = Hippurate, color = Hippurate) + 
    scale_colour_manual(name = "Hippurate:", breaks = c("0", "1", "2"), values = c("blue", 
        "forestgreen", "red")) + scale_shape_manual(name = "Hippurate:", breaks = c("0", 
    "1", "2"), values = c(4, 1, 2))

ASCAEScoresHippurate
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-21-1.png" width="800px" />

```r
# Citrate scores
ascaSVD <- SVDforPCA(EffectMatGLM$Citrate)
ascaSVD$scores[, 1:2] <- (EffectMatGLM$Citrate + EffectMatGLM$residuals) %*% 
    ascaSVD$loadings[, 1:2]

col <- c(`0` = "black", `2` = "orange", `4` = "red")
DrawScores(ascaSVD, type.obj = "PCA", drawNames = F, createWindow = F, main = "ASCA-E scores plot - Citrate effect", 
    axes = c(1, 2), size = 2.5, noLegend = TRUE, pch = Citrate, color = Citrate) + 
    scale_colour_manual(name = "Citrate:", breaks = c("0", "2", "4"), values = c("black", 
        "orange", "red")) + scale_shape_manual(name = "Citrate:", breaks = c("0", 
    "2", "4"), values = c(25, 15, 4))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-21-2.png" width="800px" />


### Scores scatter plots


```r
# scores plot
mat_ascaSVD <- do.call(cbind, ascaSVD_scores)
names(ascaSVD_scores)
```

```
## [1] "Hippurate"              "Citrate"               
## [3] "Time"                   "Hippurate:Citrate"     
## [5] "Hippurate:Time"         "Citrate:Time"          
## [7] "Hippurate:Citrate:Time" "residuals1"            
## [9] "residuals2"
```

```r
ascaSVD_scores <- c(ascaSVD_scores[p:(p + 1)], ascaSVD_scores[1:(p - 1)])

ascaSVD_scores$Time <- ascaSVD_scores$Time * -1
ascaSVD_scores$`Hippurate:Time` <- ascaSVD_scores$`Hippurate:Time` * -1


labels <- paste0(c(paste(ModelTerms_abbrev, "\n PC 1"), "Residuals \n PC 2"), 
    "\n ", round(pcexp_pereffect_perPC, 2), "%")

# plot --------------------------------------------------------------------
par(xpd = TRUE)
pairsBG(mat_ascaSVD, titre = "ASCA-E Scores plots", oma = c(3, 3, 5, 15), labels = labels, 
    spotPoints = FALSE)

legendsScatterMatrix()
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-22-1.png" width="800px" />


## Loadings plots


```r
# loadings plot
loadingsAPCAGLM <- do.call(cbind, sapply(Y_Loadings, function(x) x[, 1]))

colnames(loadingsAPCAGLM) <- paste(colnames(loadingsAPCAGLM), "PC1")
colnames(loadingsAPCAGLM)[colnames(loadingsAPCAGLM) == "Hippurate:Time PC1"] <- "HxT PC1"

loadingsAPCAGLM[, "Time PC1"] <- loadingsAPCAGLM[, "Time PC1"] * -1
loadingsAPCAGLM[, "HxT PC1"] <- loadingsAPCAGLM[, "HxT PC1"] * -1

# plot --------------------------------------------------------------------
LinePlot(t(loadingsAPCAGLM[, c(1, 2, 3, 5)]), main = "ASCA-GLM Loadings", type = "s", 
    num.stacked = 4, xaxis_type = "numerical", nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-23-1.png" width="800px" />




# APCA-GLM scores and loadings


```r
listgraphs <- list()
for (i in 1:length(ModelTerms)) {
    # PCA on the augmented effect matrices
    ascaSVD <- SVDforPCA(EffectMatGLMAug[[i]])
    ascaSVD$scores <- round(ascaSVD$scores, 5)
    # Scores plot with Citrate/hippurate color
    listgraphs[[paste0(ModelTerms[i], "-CH")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = F, createWindow = F, main = paste0(ModelTerms[i], "-CH", 
            " score plot"), color = Citrate, pch = Hippurate, axes = c(1, 2), 
        size = 2.5)
    # Scores plot with Media/Day color
    listgraphs[[paste0(ModelTerms[i], "-DM")]] <- DrawScores(ascaSVD, type.obj = "PCA", 
        drawNames = F, createWindow = F, main = paste0(ModelTerms[i], "-DM", 
            " score plot"), color = Day, pch = Dilution, axes = c(1, 2), size = 2.5)
    # Loadings plot
    listgraphs[[paste0(ModelTerms[i], "-Loadings")]] <- DrawLoadings(ascaSVD, 
        type.obj = "PCA", createWindow = F, main = paste0(ModelTerms[i], " loading plot"), 
        axes = c(1:2), loadingstype = "s", num.stacked = 2, xlab = NULL, ylab = NULL, 
        ang = "0", xaxis_type = "numerical", nxaxis = 10)
}
listgraphs
```

```
## $`Hippurate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-1.png" width="800px" />

```
## 
## $`Hippurate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-2.png" width="800px" />

```
## 
## $`Hippurate-Loadings`
## $`Hippurate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-3.png" width="800px" />

```
## 
## 
## $`Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-4.png" width="800px" />

```
## 
## $`Citrate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-5.png" width="800px" />

```
## 
## $`Citrate-Loadings`
## $`Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-6.png" width="800px" />

```
## 
## 
## $`Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-7.png" width="800px" />

```
## 
## $`Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-8.png" width="800px" />

```
## 
## $`Time-Loadings`
## $`Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-9.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-10.png" width="800px" />

```
## 
## $`Hippurate:Citrate-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-11.png" width="800px" />

```
## 
## $`Hippurate:Citrate-Loadings`
## $`Hippurate:Citrate-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-12.png" width="800px" />

```
## 
## 
## $`Hippurate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-13.png" width="800px" />

```
## 
## $`Hippurate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-14.png" width="800px" />

```
## 
## $`Hippurate:Time-Loadings`
## $`Hippurate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-15.png" width="800px" />

```
## 
## 
## $`Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-16.png" width="800px" />

```
## 
## $`Citrate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-17.png" width="800px" />

```
## 
## $`Citrate:Time-Loadings`
## $`Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-18.png" width="800px" />

```
## 
## 
## $`Hippurate:Citrate:Time-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-19.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-20.png" width="800px" />

```
## 
## $`Hippurate:Citrate:Time-Loadings`
## $`Hippurate:Citrate:Time-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-21.png" width="800px" />

```
## 
## 
## $`residuals-CH`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-22.png" width="800px" />

```
## 
## $`residuals-DM`
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-23.png" width="800px" />

```
## 
## $`residuals-Loadings`
## $`residuals-Loadings`[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-24-24.png" width="800px" />



```r
# Hippurate scores and loadings
ascaSVD <- SVDforPCA(EffectMatGLMAug$Hippurate)
ascaSVD$scores <- round(ascaSVD$scores, 5)

APCAScoresHippurate <- DrawScores(ascaSVD, type.obj = "PCA", drawNames = F, 
    createWindow = F, main = "APCA scores plot - Hippurate effect", axes = c(1, 
        2), noLegend = TRUE, size = 2.5, pch = Hippurate, color = Hippurate) + 
    scale_colour_manual(name = "Hippurate:", breaks = c("0", "1", "2"), values = c("blue", 
        "forestgreen", "red")) + scale_shape_manual(name = "Hippurate:", breaks = c("0", 
    "1", "2"), values = c(4, 1, 2)) + guides(colour = guide_legend(override.aes = list(shape = c(25, 
    15, 4), color = c("blue", "red", "red"))))

APCAScoresHippurate
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-25-1.png" width="800px" />

```r
APCALoadingHippurate <- DrawLoadings(ascaSVD, type.obj = "PCA", createWindow = F, 
    main = "APCA loading plot - Hippurate effect", axes = c(1), loadingstype = "s", 
    num.stacked = 2, xlab = "ppm", ylab = "", ang = "0", xaxis_type = "numerical", 
    nxaxis = 10)
```


# Hippurate scores and loadings for ASCA(-E) and APCA


```r
scoresASCAhip <- plot_grid(ASCAScoresHippurate, APCAScoresHippurate, ASCAEScoresHippurate, 
    align = "none", nrow = 1, rel_widths = c(0.3, 0.32, 0.39))

loadingsASCAhip <- plot_grid(ASCALoadingHippurate[[1]] + ylim(-0.2, 0.5), APCALoadingHippurate[[1]] + 
    ylim(-0.2, 0.5), align = "none", nrow = 1, rel_heights = c(1/2, 1/2))


gridExtra::grid.arrange(scoresASCAhip, loadingsASCAhip, heights = c(1.5, 1))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-26-1.png" width="800px" />




# PARAFASCA total

## Construct the array

```r
predictedValues <- unique(resGLM$predictedValues) - matrix(rep(resGLM$parameters[1, 
    ], nrow(unique(resGLM$predictedValues))), ncol = ncol(unique(resGLM$predictedValues)), 
    nrow = nrow(unique(resGLM$predictedValues)), byrow = TRUE)

modMat <- do.call(cbind, modelMatrix)

sel_rownames <- rownames(unique(modMat))
index_rownames <- match(sel_rownames, rownames(modMat))

# rownames(modMat[index_rownames,]) == rownames(predictedValues)

des <- design[index_rownames, c("Time", "Hippurate", "Citrate")]


PARAFAC_array <- array(data = NA, dim = c(3, 3, 600, 2))

funDataArray <- function(x) {
    matrix(data = x, ncol = 3, nrow = 3, dimnames = list(c("Cit1", "Cit2", "Cit3"), 
        c("Hip1", "Hip2", "Hip3")), byrow = TRUE)
}

###### by Time
for (j in 1:nlevels(des$Time)) {
    
    subdes_medium <- des[des$Time == levels(des$Time)[j], ]
    subpredictedValues_medium <- predictedValues[des$Time == levels(des$Time)[j], 
        ]
    
    modmatunique <- modMat[index_rownames, ]
    modmatunique <- modmatunique[des$Time == levels(des$Time)[j], ]
    
    rownames(subpredictedValues_medium) == rownames(modmatunique)
    
    rownames(modmatunique) <- paste0("Hip", subdes_medium$Hippurate, "Cit", 
        subdes_medium$Citrate)
    rownames(subpredictedValues_medium) <- rownames(modmatunique)
    
    subpredictedValues_medium <- subpredictedValues_medium[sort(rownames(subpredictedValues_medium)), 
        ]
    
    
    predictedValues_list <- lapply(seq_len(ncol(subpredictedValues_medium)), 
        function(i) subpredictedValues_medium[, i])
    
    decompMatrices <- lapply(predictedValues_list, funDataArray)
    
    for (i in 1:length(decompMatrices)) {
        PARAFAC_array[, , i, j] <- decompMatrices[[i]]
    }
    
}

# Select the number of components
```




## Run PARAFAC on the array

### Scree plot


```r
Rsq <- c()
for (i in 1:nparam) {
    Res_PARAFAC <- parafac(PARAFAC_array, nfac = i, nstart = 10, const = c("uncons", 
        "uncons", "orthog", "uncons"))
    # orthogonal constraint for the variables mode
    Rsq[i] = Res_PARAFAC$Rsq
}
```



```r
plot(c(0, 1:nparam), c(0, Rsq), type = "b", ylab = "R2", xlab = "n Comp", main = "PARAFASCA Scree plot", 
    mgp = c(2, 1, 0))
points(c(0, 1:nparam), c(0, Rsq), col = "blue", pch = 20)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-29-1.png" width="800px" />

```r
pdf(file.path(directoryOut, "PARAFASCAScreeplot.pdf"), width = 5, height = 4)
plot(c(0, 1:nparam), c(0, Rsq), type = "b", ylab = "R2", xlab = "n Comp", main = "PARAFASCA Scree plot", 
    mgp = c(2, 1, 0))
points(c(0, 1:nparam), c(0, Rsq), col = "blue", pch = 20)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### Estimation for a given R

```r
R <- 4
Res_PARAFAC <- parafac(PARAFAC_array, nfac = R, nstart = 100, const = c("uncons", 
    "uncons", "orthog", "uncons"), verbose = FALSE)
# orthogonal constraint for the variables mode
```


### Save the results

```r
save(Res_PARAFAC, R, file = file.path(directoryOut, "PARAFASCA_results.RData"))
```

## Load the results


```r
load(file = file.path(directoryOut, "PARAFASCA_results.RData"))
```


## Pc var per component


```r
res <- Res_PARAFAC
# Res_PARAFAC$SSE

mydim <- dim(PARAFAC_array)

# fitted array
Xhat1 <- fitted(res)

# save the modes
A <- res$A
B <- res$B
C <- res$C
D <- res$D

# SST
SST <- sumsq(PARAFAC_array)  # SST

## Compute R2

# A. approach
outprod_percomp <- vector(mode = "list", length = R)
SSE_percomp <- c()
for (i in 1:R) {
    outprod_percomp[[i]] <- outer(outer(outer(A[, i], B[, i]), C[, i]), D[, 
        i])
    E_parafac <- PARAFAC_array - outprod_percomp[[i]]
    SSE_percomp[i] <- sumsq(E_parafac)
}


R2_percomp <- 1 - (SSE_percomp/SST)


# B. approach

outprod_cumul <- vector(mode = "list", length = R)
outprod_cumul[[1]] <- outprod_percomp[[1]]
SSE_cumul <- c()
SSE_cumul[1] <- sumsq(PARAFAC_array - outprod_cumul[[1]])

for (i in 2:R) {
    outprod_cumul[[i]] <- outprod_cumul[[i - 1]] + outprod_percomp[[i]]
    E_parafac <- PARAFAC_array - outprod_cumul[[i]]
    SSE_cumul[i] <- sumsq(E_parafac)
}

# outprod_cumul[[6]][1,1,1,1] Xhat1[1,1,1,1] #=> OK!

R2_cumul <- 1 - (SSE_cumul/SST)

# comparison between A and B
plot(R2_cumul, type = "l", main = "R2 per comp and cumulated", xlab = "n components")
points(cumsum(R2_percomp), col = "red")
legend("topright", legend = c("per comp", "cumul"), lty = c(1, NA), pch = c(NA, 
    1), col = c("black", "red"))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-33-1.png" width="800px" />

```r
pander("R2 per component")
```

R2 per component

```r
names(R2_percomp) <- paste("comp", 1:R)
round(R2_percomp * 100, 2)
```

```
## comp 1 comp 2 comp 3 comp 4 
##  38.19  28.26  16.53  11.61
```

```r
pander("Cumulated R2")
```

Cumulated R2

```r
R2_cumul[R]
```

```
## [1] 0.9457823
```

```r
# Res_PARAFAC$Rsq # ==> OK
```

## all loadings/scores graphs

```r
res <- Res_PARAFAC
par(mfrow = c(1, 1))
nstack <- 1
for (i in 1:R) {
    A = as.matrix(res$A[, i])
    B = as.matrix(res$B[, i])
    C = as.matrix(res$C[, i])
    D = as.matrix(res$D[, i])
    # dimnames(A)=list(unique(Hippurate),paste0('Comp',i))
    # dimnames(B)=list(unique(Citrate),paste0('Comp',i))
    # dimnames(C)=list(Variables,paste0('Comp',i))
    a = LinePlot(t(A), main = "Hippurate")
    b = LinePlot(t(B), createWindow = FALSE, main = "Citrate", rows = 1, type = "l")
    c = LinePlot(t(C), createWindow = FALSE, main = "Variables", rows = 1, type = "s")
    d = LinePlot(t(D), createWindow = FALSE, main = "Time", rows = 1)
    cat("PARAFASCA Component ", i)
    do.call(gridExtra::grid.arrange, c(a, b, c, d, list(nrow = 4, ncol = 1)))
}
```

```
## PARAFASCA Component  1
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-34-1.png" width="800px" />

```
## PARAFASCA Component  2
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-34-2.png" width="800px" />

```
## PARAFASCA Component  3
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-34-3.png" width="800px" />

```
## PARAFASCA Component  4
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-34-4.png" width="800px" />


## Loadings modes A B D

```r
mainEffects <- list(Hippurate = Hippurate, Citrate = Citrate, Time = Time)
namesmainEffects <- names(mainEffects)
modes = c("A", "B", "D")
titles <- paste("mode", c("A", "B", "C"))
res <- Res_PARAFAC

# Manipulations sur les modes

res$A[, 1] <- res$A[, 1] * -1
res$D[, 1] <- res$D[, 1] * -1
res$B[, 1] <- res$B[, 1] * -1
res$C[, 1] <- res$C[, 1] * -1

# res$B[,2] <- res$B[,2]*-1 res$C[,2] <- res$C[,2]*-1 res$A[,3] <-
# res$A[,3]*-1 res$B[,3] <- res$B[,3]*-1 res$C[,3] <- res$C[,3]*-1 res$D[,3]
# <- res$D[,3]*-1 res$A[,4] <- res$A[,4]*-1 res$C[,4] <- res$C[,4]*-1

res$D <- res$D * 10
res$C <- res$C/10

plots <- vector("list", length(namesmainEffects))
# color <- rainbow(R)
color <- c("blue", "forestgreen", "red", "yellow3")
for (k in 1:length(modes)) {
    colnames(res[[modes[k]]]) <- paste0("comp", 1:R)
    res[[modes[k]]] <- data.frame(res[[modes[k]]], lev = levels(mainEffects[[k]]))
    res[[modes[k]]] <- as.data.frame(res[[modes[k]]])
    data = res[[modes[k]]]
    
    for (i in 1:R) {
        data[, i] <- as.numeric(data[, i])
    }
    
    data_melted <- data.frame(Components = rep(paste("comp.", 1:R), each = nrow(data)), 
        lev = rep(data$lev, times = R), value = unlist(c(data[, 1:R])))
    data_melted[, "value"] <- as.numeric(data_melted[, "value"])
    rownames(data_melted) <- NULL
    
    plots[[k]] <- ggplot(data = data_melted, aes(x = lev, y = value, colour = Components, 
        group = Components)) + geom_point(aes(shape = Components)) + geom_line(aes(linetype = Components), 
        size = 0.6) + xlab(paste(namesmainEffects[k], "levels")) + ggtitle(label = titles[k]) + 
        scale_color_manual(values = color) + theme(legend.key.size = unit(2.2, 
        "line")) + ylab("Scores")
    
    
}


plots[[1]] <- plots[[1]] + theme(legend.position = "none")

plots[[2]] <- plots[[2]] + theme(legend.position = "none") + theme(axis.title.y = element_blank())

plots[[3]] <- plots[[3]] + theme(axis.title.y = element_blank())

scoresPARAFASCA <- plot_grid(plots[[1]], plots[[2]], plots[[3]], align = "none", 
    nrow = 1, rel_widths = c(0.33, 0.3, 0.45))

scoresPARAFASCA
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-35-1.png" width="800px" />

## Loadings plot for C mode (variables)


```r
C <- res$C
colnames(C) <- paste("comp.", 1:R)
rownames(C) <- colnames(outcomes)

LinePlot(t(C[, c(1, 2, 3, 4)]), createWindow = FALSE, main = "PARAFASCA mode D loadings", 
    type = "s", num.stacked = 4, xaxis_type = c("numerical"), nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-36-1.png" width="800px" />


## Interactions visualisation

### Component 4


```r
col1 <- c("black", "orange", "blue")


A4 <- res$A[, 4]
B4 <- res$B[, 4]
C4 <- res$D[, 4]
names(A4) <- paste0("A", 1:3)
names(B4) <- paste0("B", 1:3)
names(C4) <- paste0("C", 1:2)


### interaction Plots
par(mfrow = c(2, 2))

# AB

vect <- kronecker(A4, B4, make.dimnames = TRUE)
mat <- vect
dim(mat) <- c(3, 3)
dimnames(mat) <- list(names(B4), names(A4))
# lines = mode B (Citrate)
matplot(t(mat), type = "l", ylab = "Scores", xlab = "mode A (Hippurate)", main = "Interaction plot Scores H*C", 
    col = col1)
abline(h = 0, lty = 5, col = "grey")
legend("topright", legend = c("Citrate0", "Citrate2", "Citrate4"), col = col1, 
    lty = c(1:3))

# AC
vect <- kronecker(A4, C4, make.dimnames = TRUE)
mat <- vect
dim(mat) <- c(2, 3)
dimnames(mat) <- list(names(C4), names(A4))
# lines = mode C (Time)
matplot(t(mat), type = "l", ylab = "Scores", xlab = "mode A (Hippurate)", main = "Interaction plot Scores H*T", 
    col = col1)
abline(h = 0, lty = 5, col = "grey")
legend("bottomright", legend = c("Time1", "Time2"), col = col1, lty = c(1:2))


# BC

vect <- kronecker(B4, C4, make.dimnames = TRUE)
mat <- vect
dim(mat) <- c(2, 3)
dimnames(mat) <- list(names(C4), names(B4))
matplot(t(mat), type = "l", ylab = "Scores", xlab = "mode B (Citrate)", main = "Interaction plot Scores C*T", 
    col = col1, ylim = c(-1.3, 0.5))
abline(h = 0, lty = 5, col = "grey")
legend("bottomleft", legend = c("Time1", "Time2"), col = col1, lty = c(1:2))


# ABC
vect <- kronecker(A4, B4, make.dimnames = TRUE)
vectC1 <- vect * C4[1]
vectC2 <- vect * C4[2]
dim(vectC1) <- c(3, 3)
dimnames(vectC1) <- list(names(B4), names(A4))
dim(vectC2) <- c(3, 3)
dimnames(vectC2) <- list(names(B4), names(A4))

matplot(t(vectC1), type = "l", ylab = "Scores", xlab = "mode A (Hippurate)", 
    lty = 1:3, col = col1[1], main = "Interaction plot Scores (H*C*T)")
matplot(t(vectC2), type = "l", lty = 1:3, add = TRUE, col = col1[2])
legend("topright", legend = c("Time1", "Time2", "Citrate1", "Citrate2", "Citrate3"), 
    col = c(col1[1], col1[2], "gray50", "gray50", "gray50"), lty = c(NA, NA, 
        1:3), lwd = c(NA, NA, 2, 2, 2), pch = c(15, 15, NA, NA, NA))
abline(h = 0, lty = 5, col = "grey")
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-37-1.png" width="800px" />



### All Components


```r
# Interactions
plots <- vector(mode = "list", length = 4)
dfs <- vector(mode = "list", length = 4)
names(dfs) <- paste0("comp", 1:4)

par(mfrow = c(2, 2))

for (i in 1:4) {
    
    Ai <- res$A[, i]
    Bi <- res$B[, i]
    Ci <- res$D[, i]
    names(Ai) <- paste0("A", 1:3)
    names(Bi) <- paste0("B", 1:3)
    names(Ci) <- paste0("C", 1:2)
    
    
    # ABC
    
    if (i %in% c(1, 4)) {
        col1 <- c("black", "orange", "blue")
        
        vect <- kronecker(Ai, Bi, make.dimnames = TRUE)
        vect1 <- vect * Ci[1]
        vect2 <- vect * Ci[2]
        dim(vect1) <- c(3, 3)
        dimnames(vect1) <- list(names(Bi), names(Ai))
        dim(vect2) <- c(3, 3)
        dimnames(vect2) <- list(names(Bi), names(Ai))
        
        dfs[[i]] <- list(vect1 = vect1, vect2 = vect2)
        
        if (i == 1) {
            par(mar = c(4, 4, 2, 2))
            matplot(t(dfs[[i]]$vect1), type = "l", ylab = "Mean effect scores", 
                cex.lab = 1.2, xlab = "mode A (Hippurate)", lty = 2:4, xaxt = "n", 
                col = col1[1], main = paste("Component", i))
            matplot(t(dfs[[i]]$vect2), type = "l", lty = 2:4, add = TRUE, col = col1[2])
            axis(side = 1, at = c(1, 2, 3), labels = c(0, 1, 2))
            
            legend("bottomright", legend = c("T=1", "T=2", "C=0", "C=2", "C=4"), 
                col = c(col1[1], col1[2], "gray50", "gray50", "gray50"), lty = c(NA, 
                  NA, 2:4), lwd = c(NA, NA, 2, 2, 2), pch = c(15, 15, NA, NA, 
                  NA))
        } else {
            par(mar = c(4, 2, 2, 2))
            matplot(t(dfs[[i]]$vect1), type = "l", cex.lab = 1.2, ylab = "", 
                xlab = "mode A (Hippurate)", lty = 2:4, xaxt = "n", col = col1[1], 
                main = paste("Component", i))
            matplot(t(dfs[[i]]$vect2), type = "l", lty = 2:4, add = TRUE, col = col1[2])
            axis(side = 1, at = c(1, 2, 3), labels = c(0, 1, 2))
            legend("topright", legend = c("T=1", "T=2", "C=0", "C=2", "C=4"), 
                col = c(col1[1], col1[2], "gray50", "gray50", "gray50"), lty = c(NA, 
                  NA, 2:4), lwd = c(NA, NA, 2, 2, 2), pch = c(15, 15, NA, NA, 
                  NA))
        }
        
        
        
    } else if (i == 2) {
        col1 <- c("black", "orange", "blue")
        par(mar = c(4, 2, 2, 2))
        vect <- kronecker(Bi, Ai, make.dimnames = TRUE)
        vect1 <- vect * Ci[1]
        vect2 <- vect * Ci[2]
        dim(vect1) <- c(3, 3)
        dimnames(vect1) <- list(names(Ai), names(Bi))
        dim(vect2) <- c(3, 3)
        dimnames(vect2) <- list(names(Ai), names(Bi))
        
        
        dfs[[i]] <- list(vect1 = vect1, vect2 = vect2)
        
        matplot(t(dfs[[i]]$vect1), type = "l", cex.lab = 1.2, ylab = "", xlab = "mode B (Citrate)", 
            lty = 2:4, xaxt = "n", col = col1[1], main = paste("Component", 
                i))
        matplot(t(dfs[[i]]$vect2), type = "l", lty = 2:4, add = TRUE, col = col1[2])
        
        axis(side = 1, at = c(1, 2, 3), labels = c(0, 2, 4))
        legend("bottomright", legend = c("T=1", "T=2", "H=0", "H=1", "H=2"), 
            col = c(col1[1:2], "gray50", "gray50", "gray50"), lty = c(NA, NA, 
                2:4), lwd = c(NA, NA, 2, 2, 2), pch = c(15, 15, NA, NA, NA))
        
        
    } else {
        col1 <- c("black", "orange", "blue")
        par(mar = c(4, 2, 2, 2))
        vect <- kronecker(Bi, Ci, make.dimnames = TRUE)
        vect1 <- vect * Ai[1]
        vect2 <- vect * Ai[2]
        vect3 <- vect * Ai[3]
        dim(vect1) <- c(2, 3)
        dimnames(vect1) <- list(names(Ci), names(Bi))
        dim(vect2) <- c(2, 3)
        dimnames(vect2) <- list(names(Ci), names(Bi))
        dim(vect3) <- c(2, 3)
        dimnames(vect3) <- list(names(Ci), names(Bi))
        
        dfs[[i]] <- list(vect1 = t(vect1), vect2 = t(vect2), vect3 = t(vect3))
        
        matplot(t(dfs[[i]]$vect1), type = "l", cex.lab = 1.2, ylab = "", xlab = "mode C (Time)", 
            lty = 2:4, xaxt = "n", col = col1[1], main = paste("Component", 
                i))
        matplot(t(dfs[[i]]$vect2), type = "l", lty = 2:4, add = TRUE, col = col1[2])
        matplot(t(dfs[[i]]$vect3), type = "l", lty = 2:4, add = TRUE, col = col1[3])
        axis(side = 1, at = c(1, 2), labels = c(1, 2))
        legend("topright", legend = c("H=0", "H=1", "H=2", "C=0", "C=2", "C=4"), 
            col = c(col1, "gray50", "gray50", "gray50"), lty = c(NA, NA, NA, 
                2:4), lwd = c(NA, NA, NA, 2, 2, 2), pch = c(15, 15, 15, NA, 
                NA, NA))
        
    }
    
}
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-38-1.png" width="800px" />



# PARAFASCA submodels with H and T


```r
predictedValues <- vector(mode = "list", length = 3)
names(predictedValues) <- c("H+HXT", "T+HXT", "HXT")
predictedValues[["H+HXT"]] <- unique(resGLM$effectMatrices$Hippurate + resGLM$effectMatrices$`Hippurate:Time`)
predictedValues[["T+HXT"]] <- unique(resGLM$effectMatrices$Time + resGLM$effectMatrices$`Hippurate:Time`)
predictedValues[["HXT"]] <- unique(resGLM$effectMatrices$`Hippurate:Time`)

des <- design[index_rownames, c("Time", "Hippurate", "Citrate")]
# table(des$Time, des$Hippurate)


PARAFAC_array_list <- vector(mode = "list", length = 3)
names(PARAFAC_array_list) <- c("H+HXT", "T+HXT", "HXT")

for (j in 1:3) {
    des[rownames(predictedValues[[j]]), c("Time", "Hippurate")]
    PARAFAC_array <- array(data = NA, dim = c(3, 2, 600))
    
    funDataArray <- function(x) {
        matrix(data = x, ncol = 2, nrow = 3, dimnames = list(c("Hip1", "Hip2", 
            "Hip3"), c("Time1", "Time2")), byrow = TRUE)
    }
    
    predictedValues_list <- lapply(seq_len(ncol(predictedValues[[j]])), function(i) predictedValues[[j]][, 
        i])
    decompMatrices <- lapply(predictedValues_list, funDataArray)
    
    for (i in 1:length(decompMatrices)) {
        PARAFAC_array[, , i] <- decompMatrices[[i]]
    }
    
    PARAFAC_array_list[[j]] <- PARAFAC_array
}

par(mfrow = c(4, 1), mar = c(2, 2, 2, 2))
plot(PARAFAC_array_list[[1]][1, 1, ], type = "l", main = "[1,1,]")
plot(PARAFAC_array_list[[1]][3, 1, ], type = "l", main = "[3,1,]")
plot(PARAFAC_array_list[[1]][1, 2, ], type = "l", main = "[1,2,]")
plot(PARAFAC_array_list[[1]][3, 2, ], type = "l", main = "[3,2,]")
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-39-1.png" width="800px" />

## Scree plot




```r
for (j in 1:2) {
    
    Rsq <- c()
    for (i in 1:nparam) {
        
        Res_PARAFAC <- parafac(PARAFAC_array_list[[j]], nfac = i, nstart = 10, 
            const = c("uncons", "uncons", "orthog"), verbose = FALSE)
        # orthogonal constraint for the variables mode
        Rsq[i] = Res_PARAFAC$Rsq
    }
    
    plot(c(0, 1:nparam), c(0, Rsq), type = "b", ylab = "R2", xlab = "n Comp", 
        main = paste0("PARAFASCA Scree plot\n", names(PARAFAC_array_list)[j]), 
        mgp = c(2, 1, 0))
    points(c(0, 1:nparam), c(0, Rsq), col = "blue", pch = 20)
    
}
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-40-1.png" width="800px" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-40-2.png" width="800px" />


## Estimation with a given R

```r
R <- 2
Res_PARAFAC <- vector(mode = "list", length = 3)
for (j in 1:2) {
    Res_PARAFAC[[j]] <- parafac(PARAFAC_array_list[[j]], nfac = R, nstart = 100, 
        const = c("uncons", "uncons", "orthog"), verbose = FALSE)
}

names(Res_PARAFAC) <- names(PARAFAC_array_list)
```


### Save the results

```r
save(Res_PARAFAC, R, file = file.path(directoryOut, "PARAFASCA_Sub_results.RData"))
```

## Load the results

```r
load(file = file.path(directoryOut, "PARAFASCA_Sub_results.RData"))
```


## all loadings/scores graphs


```r
for (j in 1:2) {
    cat(names(Res_PARAFAC)[j], "\n")
    res <- Res_PARAFAC[[j]]
    par(mfrow = c(1, 1))
    nstack <- 1
    for (i in 1:R) {
        A = as.matrix(res$A[, i])
        B = as.matrix(res$B[, i])
        C = as.matrix(res$C[, i])
        
        a = LinePlot(t(A), main = "Hippurate")
        b = LinePlot(t(B), createWindow = FALSE, main = "Time", rows = 1, type = "l")
        c = LinePlot(t(C), createWindow = FALSE, main = "Variables", rows = 1, 
            type = "s")
        
        cat("PARAFASCA Component ", i)
        do.call(gridExtra::grid.arrange, c(a, b, c, list(nrow = 3, ncol = 1)))
    }
}
```

```
## H+HXT 
## PARAFASCA Component  1
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-44-1.png" width="800px" />

```
## PARAFASCA Component  2
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-44-2.png" width="800px" />

```
## T+HXT 
## PARAFASCA Component  1
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-44-3.png" width="800px" />

```
## PARAFASCA Component  2
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-44-4.png" width="800px" />



```r
res <- Res_PARAFAC[[1]]

A1 <- res$A[, 1]
A2 <- res$A[, 2]
B1 <- res$B[, 1]
B2 <- res$B[, 2]
C1 <- res$C[, 1]
C2 <- res$C[, 2]

y11 <- A1[1] * B1[1] * C1 + A2[1] * B2[1] * C2
plot(y11, type = "l")

y31 <- A1[3] * B1[1] * C1 + A2[3] * B2[1] * C2
plot(y31, type = "l")

y12 <- A1[1] * B1[2] * C1 + A2[1] * B2[2] * C2
plot(y12, type = "l")

y31 <- A1[3] * B1[2] * C1 + A2[3] * B2[2] * C2
plot(y31, type = "l")
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-45-1.png" width="49%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-45-2.png" width="49%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-45-3.png" width="49%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-45-4.png" width="49%" />


## Loadings


```r
load <- vector(mode = "list", length = 2)
for (j in 1:2) {
    cat(names(Res_PARAFAC)[j], "\n")
    res <- Res_PARAFAC[[j]]
    nstack <- 1
    c <- vector(mode = "list", length = R)
    for (i in 1:R) {
        C = as.matrix(res$C[, i])
        c[[i]] = LinePlot(t(C), createWindow = FALSE, main = paste(names(Res_PARAFAC)[j], 
            "- Mode C", "- Comp", i), rows = 1, type = "s")
        # d = LinePlot(t(D), createWindow = FALSE, main = 'Time', rows = 1)
        cat("PARAFASCA Component ", i)
        
    }
    load[[j]] <- plot_grid(c[[1]][[1]], c[[2]][[1]], align = "none", nrow = 2, 
        rel_widths = c(0.5, 0.5))
}
```

```
## H+HXT 
## PARAFASCA Component  1PARAFASCA Component  2T+HXT 
## PARAFASCA Component  1PARAFASCA Component  2
```

```r
plot_grid(load[[1]], load[[2]], align = "none", nrow = 2, rel_widths = c(0.5, 
    0.5))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-46-1.png" width="800px" />


## Loadings modes A B C

```r
mainEffects <- list(Hippurate = Hippurate, Time = Time)
namesmainEffects <- names(mainEffects)
modes = c("A", "B")
scoresPARAFASCA <- vector(mode = "list", length = 2)
LoadPARAFASCA <- vector(mode = "list", length = 2)

color_HHT <- c("blue", "yellow3")
color_THT <- c("red", "yellow3")

for (j in 1:2) {
    
    titles <- paste(names(Res_PARAFAC[j]), "- mode", c("A", "B"), "scores")
    res <- Res_PARAFAC[[j]]
    
    colnames(res$C) <- paste("comp.", 1:R)
    rownames(res$C) <- colnames(outcomes)
    
    if (j == 1) {
        res$C[, 1] <- res$C[, 1] * -1
        res$B[, 1] <- res$B[, 1] * -1
        
        res$A[, 2] <- res$A[, 2] * -1
        res$C[, 2] <- res$C[, 2] * -1
        
        color <- color_HHT
    } else {
        res$A[, 1] <- res$A[, 1] * -1
        res$C[, 1] <- res$C[, 1] * -1
        res$C[, 2] <- res$C[, 2] * -1
        res$A[, 2] <- res$A[, 2] * -1
        
        color <- color_THT
    }
    
    
    plots <- vector("list", length(namesmainEffects))
    
    for (k in 1:length(modes)) {
        colnames(res[[modes[k]]]) <- paste0("comp", 1:R)
        res[[modes[k]]] <- data.frame(res[[modes[k]]], lev = levels(mainEffects[[k]]))
        res[[modes[k]]] <- as.data.frame(res[[modes[k]]])
        data = res[[modes[k]]]
        
        for (i in 1:R) {
            data[, i] <- as.numeric(data[, i])
        }
        
        data_melted <- data.frame(Components = rep(paste("comp.", 1:R), each = nrow(data)), 
            lev = rep(data$lev, times = R), value = unlist(c(data[, 1:R])))
        data_melted[, "value"] <- as.numeric(data_melted[, "value"])
        rownames(data_melted) <- NULL
        
        plots[[k]] <- ggplot(data = data_melted, aes(x = lev, y = value, colour = Components, 
            group = Components)) + geom_point(aes(shape = Components)) + geom_line(aes(linetype = Components), 
            size = 0.6) + xlab(paste(namesmainEffects[k], "levels")) + ggtitle(label = titles[k]) + 
            scale_color_manual(values = color) + theme(legend.key.size = unit(2.2, 
            "line")) + ylab("Scores")
        
        
    }
    
    
    plots[[1]] <- plots[[1]] + theme(legend.position = "none")
    
    plots[[2]] <- plots[[2]] + theme(axis.title.y = element_blank())
    
    
    scoresPARAFASCA[[j]] <- plot_grid(plots[[1]], plots[[2]], align = "none", 
        nrow = 1, rel_widths = c(0.5, 0.5))
    
    LoadPARAFASCA[[j]] <- LinePlot(t(res$C[, c(1, 2)]), createWindow = FALSE, 
        main = "mode C loadings", type = "s", num.stacked = 4, xaxis_type = c("numerical"), 
        nxaxis = 10)
    
    ALLscoresPARAFASCA <- plot_grid(scoresPARAFASCA[[1]], scoresPARAFASCA[[2]], 
        align = "none", nrow = 2, rel_widths = c(0.5, 0.5))
    
    ALLscoresPARAFASCA
}


scoresPARAFASCA[[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-47-1.png" width="800px" height="50%" />

```r
LoadPARAFASCA[[1]][[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-47-2.png" width="800px" height="50%" />

```r
scoresPARAFASCA[[2]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-47-3.png" width="800px" height="50%" />

```r
LoadPARAFASCA[[2]][[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-47-4.png" width="800px" height="50%" />

## Interaction plots

### H+HxT and T+HxT

```r
res <- Res_PARAFAC[[1]]

# remettre les scores dans le bon sens
res$C[, 1] <- res$C[, 1] * -1
res$B[, 1] <- res$B[, 1] * -1

res$A[, 2] <- res$A[, 2] * -1
res$C[, 2] <- res$C[, 2] * -1

# Comp 2 = interation
A2 <- res$A[, 2]
B2 <- res$B[, 2]
A1 <- res$A[, 1]
B1 <- res$B[, 1]

names(A1) <- names(A2) <- paste0("A", 1:3)
names(B1) <- names(B2) <- paste0("B", 1:2)


col1 <- c("black", "orange")

# AB comp 1
vect <- kronecker(A1, B1, make.dimnames = TRUE)
mat <- vect
dim(mat) <- c(2, 3)
data <- c(mat[1, ], mat[2, ])
data_melted <- data.frame(Hippurate = as.factor(c(0, 1, 2, 0, 1, 2)), Time = as.factor(c(rep("Time1", 
    3), rep("Time2", 3))), value = data)

col1 <- c("black", "orange")
ABComp1 <- ggplot(data = data_melted, aes(x = Hippurate, y = value, colour = Time, 
    group = Time)) + geom_point(aes(shape = Time)) + geom_line(aes(linetype = Time), 
    size = 0.6) + xlab("Hippurate levels") + ylab("Scores") + ggtitle(label = "Interaction plot scores (H*T) \nComponent 1") + 
    scale_color_manual(values = col1) + theme(legend.key.size = unit(2.2, "line"))


# AB comp 2
vect <- kronecker(A2, B2, make.dimnames = TRUE)
mat <- vect
dim(mat) <- c(2, 3)
data <- c(mat[1, ], mat[2, ])
data_melted <- data.frame(Hippurate = as.factor(c(0, 1, 2, 0, 1, 2)), Time = as.factor(c(rep("Time1", 
    3), rep("Time2", 3))), value = data)

col1 <- c("black", "orange")
ABComp2 <- ggplot(data = data_melted, aes(x = Hippurate, y = value, colour = Time, 
    group = Time)) + geom_point(aes(shape = Time)) + geom_line(aes(linetype = Time), 
    size = 0.6) + xlab("Hippurate levels") + ylab("Scores") + ggtitle(label = "Interaction plot scores (H*T) \nComponent 2") + 
    scale_color_manual(values = col1) + theme(legend.key.size = unit(2.2, "line"))


# load C comp 1 et 2
C <- t(res$C[, c(1, 2)])
rownames(C) <- c("Comp. 1", "Comp. 2")
load <- LinePlot(C, createWindow = FALSE, main = "mode C loadings", type = "s", 
    num.stacked = 4, xaxis_type = c("numerical"), nxaxis = 10)[[1]]



plot_grid(ABComp1, ABComp2, align = "none", nrow = 1)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-48-1.png" width="800px" />

```r
load
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-48-2.png" width="800px" />

### HxT

```r
Res_PARAFAC_HXT <- parafac(PARAFAC_array_list$HXT, nfac = 1, nstart = 100, const = c("uncons", 
    "uncons", "orthog"), verbose = FALSE)

A <- Res_PARAFAC_HXT$A[, 1] * -1
B <- Res_PARAFAC_HXT$B[, 1] * -1

names(A) <- paste0("A", 1:3)
names(B) <- paste0("B", 1:2)

vect <- kronecker(A, B, make.dimnames = TRUE)

mat <- vect
dim(mat) <- c(2, 3)
dimnames(mat) <- list(names(B), names(A))

matplot(t(mat), type = "l", ylab = "Mean effect scores", lty = 2:3, xaxt = "n", 
    lwd = 2.3, xlab = "mode A (Hippurate)", main = "", col = col1, cex.lab = 1.2)
axis(side = "1", at = c(1, 2, 3), labels = c(0, 1, 2))
legend("topright", legend = c("T=1", "T=2"), lwd = 2, col = col1, lty = c(2:3))
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-49-1.png" width="800px" />

# AComDim

## Scree plot


```r
pcvarAComDim <- function(R) {
    resComDimGLM <- ComDim(Y = EffectMatGLMAug, R = R, thresh = 1e-20, stand = TRUE)
    saliences <- resComDimGLM$lambda
    ### New computation by ComDim
    mu_r <- apply(saliences^2, 2, sum)
    ### Centering and Normalisation of the augmented effect matrices
    EffectMatGLMAugN = EffectMatGLMAug
    
    for (i in 1:dim(saliences)[1]) {
        EffectMatGLMAugN[[i]] <- scale(EffectMatGLMAug[[i]], center = TRUE, 
            scale = FALSE)
        EffectMatGLMAugN[[i]] <- EffectMatGLMAugN[[i]]/(norm(EffectMatGLMAugN[[i]], 
            type = "F"))
    }
    
    ## denominator computation
    denom <- 0
    for (i in 1:dim(saliences)[1]) {
        denom <- denom + norm(EffectMatGLMAugN[[i]] %*% t(EffectMatGLMAugN[[i]]), 
            type = "F")^2
    }
    
    PCexp <- mu_r * 100/denom
    CSPCexp <- cumsum(PCexp)
    names(PCexp) <- paste("CD", 1:R)
    
    CSPCexp[R]
}

pcvar <- c()
for (i in 1:nparam) {
    pcvar[i] <- pcvarAComDim(i)
}
pcvar <- pcvar/100


plot(c(0, 1:nparam), c(0, pcvar), type = "b", yaxt = "n", ylab = "R2", xlab = "n CD", 
    main = "AComDim Scree plot", mgp = c(2, 1, 0))
axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
points(c(0, 1:nparam), c(0, pcvar), col = "blue", pch = 20)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-50-1.png" width="800px" />

## Model

```r
R <- 6  # number of Common Components

### AComDim Parameters Y, #original data as a list of matrices R, #number of
### common dimensions to compute thresh=1E-20, #convergence threshold
### stand=TRUE #Standardisation or not

resComDimGLM <- ComDim(Y = EffectMatGLMAug, R = R, thresh = 1e-20, stand = TRUE)


# Eigen value by dimension (2 different computations)
pander("Eigen value by dimension \n")
```

Eigen value by dimension 

```r
resComDimGLM$muvec
```

```
## [1] 0.74479407 0.78476808 0.74291403 0.64212824 0.08870175 0.31122156
```

```r
# resComDimGLM$muvec2


pander("Saliences \n")
```

Saliences 

```r
pander(resComDimGLM$lambda)
```


--------- ----------- ----------- ----------- --------- -----------
 0.04485    0.8858      2.8e-05    0.0002705   0.01487    9.7e-05  

 0.05754    0.00213     0.8619     0.0002914   0.01875   0.0001003 

 0.09757   0.005143    0.001763     0.8013     0.02408   0.0001461 

 0.3848    0.0008275   0.0004839   0.002862    0.1216    0.001822  

 0.1862    0.007833    0.002061    0.007854    0.06546    0.5579   

 0.4336    0.0002973   0.0002514   0.001874    0.1449    0.0008781 

 0.3652    0.0007693   0.0005302   0.002561    0.1472    0.001756  

 0.4751    5.65e-05    0.0002314   0.002182    0.1607    0.0007311 
--------- ----------- ----------- ----------- --------- -----------

```r
pander("Explained variance")
```

Explained variance

```r
resComDimGLM$expl
```

```
## [1] 20.443497 21.540724 20.391893 17.625471  2.434732  8.542572
```

```r
pander("Total explained variance")
```

Total explained variance

```r
sum(resComDimGLM$expl)
```

```
## [1] 90.97889
```



## Saliences plots

```r
saliences <- resComDimGLM$lambda
colnames(saliences) <- paste0("CC", 1:R)
rownames(saliences) <- ModelTerms_abbrev

ylim <- c(0, max(saliences))
angle1 <- rep(c(45, 45, 135), length.out = p)
angle2 <- rep(c(45, 135, 135), length.out = p)
density1 <- seq(5, 35, length.out = p)
density2 <- seq(5, 35, length.out = p)
col <- rainbow(p)

# plot --------------------------------------------------------------------
par(xpd = TRUE, mar = c(5, 4, 4, 4))
x <- barplot(saliences, beside = TRUE, ylim = ylim, yaxt = "n", xaxt = "n", 
    col = col, angle = angle1, density = density1, main = "AComDim Saliences", 
    ann = FALSE, axes = FALSE)
barplot(saliences, add = TRUE, beside = TRUE, ylim = ylim, xaxt = "n", col = col, 
    angle = angle2, density = density2, ann = FALSE, axes = FALSE)
usr <- par("usr")
par(usr = c(usr[1:2], 0, ylim[2]))
axis(2, at = seq(0, 1, 0.1))

labs <- as.vector(colnames(saliences))
text(cex = 1, x = colMeans(x) - 0.25, y = -0.06, labs, xpd = TRUE, srt = 45)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-52-1.png" width="800px" />


## Percentage explained by effect and by dimension

```r
# by effect
tab <- rowSums(resComDimGLM$lambda)
names(tab) <- ModelTerms
pander("Pourcentage explained by effect")
```

Pourcentage explained by effect

```r
pander(tab * 100)
```


--------------------------------------------------------------------------------
 Hippurate   Citrate   Time   Hippurate:Citrate   Hippurate:Time   Citrate:Time 
----------- --------- ------ ------------------- ---------------- --------------
   94.59      94.07     93          51.24             82.73           58.19     
--------------------------------------------------------------------------------

Table: Table continues below

 
------------------------------------
 Hippurate:Citrate:Time   residuals 
------------------------ -----------
          51.8              63.9    
------------------------------------

```r
### New computation by ComDim
mu_r <- apply(saliences^2, 2, sum)

### Centering and Normalisation of the augmented effect matrices
EffectMatGLMAugN = EffectMatGLMAug

for (i in 1:dim(saliences)[1]) {
    EffectMatGLMAugN[[i]] = scale(EffectMatGLMAug[[i]], center = TRUE, scale = FALSE)
    EffectMatGLMAugN[[i]] = EffectMatGLMAugN[[i]]/(norm(EffectMatGLMAugN[[i]], 
        type = "F"))
}

## Denominator computation
denom = 0
for (i in 1:dim(saliences)[1]) {
    denom = denom + norm(EffectMatGLMAugN[[i]] %*% t(EffectMatGLMAugN[[i]]), 
        type = "F")^2
}

PCexp <- mu_r * 100/denom
CSPCexp = cumsum(PCexp)
names(PCexp) <- paste("CD", 1:R)
pander("Percentage explained by dimension")
```

Percentage explained by dimension

```r
pander(rbind(PCexp, CSPCexp))
```


-------------------------------------------------------------
   &nbsp;      CD 1    CD 2    CD 3    CD 4    CD 5    CD 6  
------------- ------- ------- ------- ------- ------- -------
  **PCexp**    20.44   21.54   20.39   17.63   2.435   8.543 

 **CSPCexp**   20.44   41.98   62.38    80     82.44   90.98 
-------------------------------------------------------------


## Scores scatterplot matrix


```r
# Score scatterplot matrix
dimnames(resComDimGLM$q)[[2]] = paste0("ComDim ", 1:R)
dimnames(resComDimGLM$q)[[1]] <- rownames(outcomes)

par(xpd = TRUE)
dimnames(resComDimGLM$q)[[2]] = paste0("ComDim ", 1:R)
# labels <- paste0(colnames(resComDimGLM$q), '\n ',round(PCexp,2), '%')

labels <- paste0(colnames(resComDimGLM$q), " \n", round(PCexp, 2), "%", c("\n \"Residuals\" ", 
    "\n \"Hippurate\"", "\n \"Citrate\" ", "\n \"Time\"", "\n \"Residuals\" ", 
    "\n \"HxT\""))


resComDimGLM$q[, 4] <- resComDimGLM$q[, 4] * -1
resComDimGLM$q[, 6] <- resComDimGLM$q[, 6] * -1


par(xpd = TRUE)
pairsBG(resComDimGLM$q, titre = "AComDim Scores plots", oma = c(3, 3, 5, 15), 
    labels = labels, spotPoints = FALSE)
legendsScatterMatrix()
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-54-1.png" width="800px" />



## Loadings (weighted)

```r
matload <- matrix(nrow = R, ncol = m)
dimnames(matload)[[2]] <- Variables
dimnames(matload)[[1]] <- paste("CC", 1:R)

for (i in 1:R) {
    matload[i, ] <- resComDimGLM$q[, i] %*% Reduce("+", Map("*", EffectMatGLMAugN, 
        sqrt(resComDimGLM$lambda[, i])))
}

# graphique
LinePlot(matload, createWindow = FALSE, type = "s", main = "Loadings", rows = 1:R, 
    num.stacked = R, xaxis_type = c("numerical"), nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-55-1.png" width="800px" />

```r
# hipurate, citrate and dilution main effects loadings
LinePlot(matload[c(2:4, 6), ], createWindow = FALSE, main = "AComDim Loadings", 
    type = "s", num.stacked = 4, xaxis_type = c("numerical"), nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-55-2.png" width="800px" />

## F-test Acomdim


```r
# Prepare the results matrix
AComDimFtest = matrix(nrow = p - 1, ncol = 3)
dimnames(AComDimFtest)[[1]] = ModelTerms[1:(p - 1)]
dimnames(AComDimFtest)[[2]] = c("Salience Dim 1", "F Stat", "P-value")
# Saliences for the first dimension
AComDimFtest[, 1] = resComDimGLM$lambda[1:(p - 1), 1]
# Test stat computation and p-values
AComDimFtest[, 2] = resComDimGLM$lambda[p, 1]/resComDimGLM$lambda[1:(p - 1), 
    1]
AComDimFtest[, 3] = pf(AComDimFtest[, 2], n - 1, n - 1, lower.tail = F)
pander(AComDimFtest)
```


------------------------------------------------------------------
           &nbsp;             Salience Dim 1   F Stat    P-value  
---------------------------- ---------------- -------- -----------
       **Hippurate**             0.04485       10.59    4.395e-10 

        **Citrate**              0.05754       8.258    1.27e-08  

          **Time**               0.09757        4.87    8.279e-06 

   **Hippurate:Citrate**          0.3848       1.235     0.2741   

     **Hippurate:Time**           0.1862       2.552    0.004339  

      **Citrate:Time**            0.4336       1.096     0.3972   

 **Hippurate:Citrate:Time**       0.3652       1.301      0.227   
------------------------------------------------------------------

```r
round(AComDimFtest[, 3], 4)
```

```
##              Hippurate                Citrate                   Time 
##                 0.0000                 0.0000                 0.0000 
##      Hippurate:Citrate         Hippurate:Time           Citrate:Time 
##                 0.2741                 0.0043                 0.3972 
## Hippurate:Citrate:Time 
##                 0.2270
```


## Permutations BG 3

Approche où toutes les matrices des effets pour les données permutées ont déjà été calculées plus haut. 


```r
################################ Recalcul des vraie valeur de la stat F
TRUE_F_stat <- resComDimGLM$lambda[p, 1]/resComDimGLM$lambda[1:(p - 1), 1]
names(TRUE_F_stat) <- ModelTerms_sansE
outcomes_true <- outcomes

################################ Permutations main effects

namesmainEffects <- ModelTerms_sansE[!grepl(":", ModelTerms_sansE)]
mainEffects <- design[, namesmainEffects]
nMainEff <- length(mainEffects)
F_stat <- vector("list", nMainEff)
names(F_stat) <- namesmainEffects

ptm <- proc.time()
kn = 0
for (k in namesmainEffects) {
    cat("\n k = ", k)
    kn = kn + 1
    for (i in 1:nperm) {
        cat("\n k = ", k, " i ", i)
        # construction of the augmented effects matrix list
        EffectMatGLMAug_perm <- DecompEffectMatrixPermut[[k]][[i]][-1]  # minus intercept
        MatRes = DecompEffectMatrixPermut[[k]][[i]]$residuals
        EffectMatGLMAug_perm <- lapply(EffectMatGLMAug_perm, function(x) x + 
            MatRes)
        EffectMatGLMAug_perm$residuals <- EffectMatGLMAug_perm$residuals/2
        
        # AComDim Application
        resComDimGLM_p <- ComDim(EffectMatGLMAug_perm, R)
        
        # computation of the F-stat for permuted data and a given k
        Fstatp = resComDimGLM_p$lambda[p, 1]/resComDimGLM_p$lambda[kn, 1]
        
        # adding it to the vector
        F_stat[[k]] <- c(F_stat[[k]], Fstatp)
    }
    
}

proc.time() - ptm
mainEffects_F_stat <- F_stat

################################ Interactions
nIntEff = length(namesinteractionEffects)
IndIntEff = nMainEff + 1:nIntEff
F_stat <- vector("list", nIntEff)
names(F_stat) <- namesinteractionEffects


ptm <- proc.time()
for (i in 1:nperm) {
    cat("\n i= ", i)
    # Computation of the augmented effect matrix
    EffectMatGLMAug_perm <- DecompEffectMatrixPermut[["Interaction"]][[i]][-1]  # minus intercept
    MatRes <- DecompEffectMatrixPermut[["Interaction"]][[i]]$residuals  # Residuals matrix
    EffectMatGLMAug_perm <- lapply(EffectMatGLMAug_perm, function(x) x + MatRes)  # Augmentation of all the effect matrices (residuals included)
    EffectMatGLMAug_perm$residuals = EffectMatGLMAug_perm$residuals/2  #Dividing the residuals by 2 (because taken 2 times)
    # ComDim application
    resComDimGLM_p <- ComDim(EffectMatGLMAug_perm, R)
    # Computation of the stat F
    Fstatp <- resComDimGLM_p$lambda[p, 1]/resComDimGLM_p$lambda[IndIntEff, 1]
    # Adding the statistics to the vector
    for (ki in 1:nIntEff) {
        F_stat[[namesinteractionEffects[ki]]] <- c(F_stat[[namesinteractionEffects[ki]]], 
            Fstatp[ki])
    }
}

proc.time() - ptm

# Gathering all the computed F stat (main effects and interactions)
interEffects_F_stat <- F_stat
F_statPermuted <- c(mainEffects_F_stat, interEffects_F_stat)

# p-values computation
pval <- c()
for (i in ModelTerms_sansE) {
    pval[i] <- sum(TRUE_F_stat[i] <= F_statPermuted[[i]])/nperm
}
names(pval) <- ModelTerms_sansE

# Saving permutations
save(F_statPermuted, pval, TRUE_F_stat, nperm, file = file.path(directoryOut, 
    "PermutationResAComDim.RData"))

# printing p-values
pander(pval)
rm(list = ls())
```

## Result permut

```r
load(file.path(directoryOut, "PermutationResAComDim.RData"))

names(pval) <- names(TRUE_F_stat) <- names(F_statPermuted) <- ModelTerms_sansE

pander(pval)
```


--------------------------------------------------------------------------------
 Hippurate   Citrate   Time   Hippurate:Citrate   Hippurate:Time   Citrate:Time 
----------- --------- ------ ------------------- ---------------- --------------
     0          0       0           0.366               0             0.499     
--------------------------------------------------------------------------------

Table: Table continues below

 
------------------------
 Hippurate:Citrate:Time 
------------------------
         0.179          
------------------------

```r
for (i in ModelTerms_sansE) {
    hist(F_statPermuted[[i]], xlim = c(0, max(TRUE_F_stat[i], F_statPermuted[[i]])), 
        breaks = 20, main = i, freq = FALSE)
    points(TRUE_F_stat[i], 0, col = "red", pch = 19)
    rangex = range(TRUE_F_stat[i], F_statPermuted[[i]])
    xf = seq(0, rangex[2], by = (rangex[2] - rangex[1])/1000)
    yf = df(xf, n - 1, n - 1)
    lines(xf, yf)
}
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-1.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-2.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-3.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-4.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-5.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-6.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-57-7.png" width="50%" />

# AMOPLS 

## Scree plot


```r
# Choose the number of predictive and orthogonal components (generally 1 or
# 2)
ortho <- 1

### Scree plot R2Yhat
R2Yhat_perPC <- c()
for (i in 1:nparam) {
    res <- AMOPLS_function(EffectMatGLM, npred = i, ortho, X = Xmat)
    R2Yhat_perPC[i] <- res[[1]]$R2Yhat[ortho + 1]
}


plot(c(0, 1:nparam), c(0, R2Yhat_perPC), type = "b", yaxt = "n", ylab = "R2", 
    xlab = "n PrC", main = "AMOPLS Scree plot", mgp = c(2, 1, 0))
axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
points(c(0, 1:nparam), c(0, R2Yhat_perPC), col = "blue", pch = 20)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-58-1.png" width="800px" />


## Model


```r
npred <- 4

### Scree plot R2XO
res <- AMOPLS_function(EffectMatGLM, npred = npred, ortho = 10, X = Xmat, UseX = FALSE)
pander("R2XO")
```

R2XO

```r
res[[1]]$R2XO
```

```
##  [1] 0.0000000 0.2684476 0.3915285 0.4504768 0.4862545 0.5059795 0.5311376
##  [8] 0.5488920 0.5612256 0.5722903 0.5874261
```

```r
plot(c(0, 1:10), c(res[[1]]$R2XO), type = "b", main = "scree plot AMOPLS R2XO", 
    xlab = "n pred comp")
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-59-1.png" width="800px" />

```r
#### run AMOPLS with appropriate n pred and ortho comp
ortho <- 2
# Run AMOPLS on the pure matrices
res <- AMOPLS_function(EffectMatGLM, npred, ortho, X = Xmat, UseX = FALSE)

# Store the results
AMOPLS_model <- res[[1]]
AMOPLS_lambda <- res[[2]]  #saliences, rows = matrices, cols = components

A <- AMOPLS_model$A  #Number of explanatory components

# Show the contribution lambda
pander("the contribution lambda")
```

the contribution lambda

```r
pander(round(AMOPLS_lambda, 2))
```


------ ------ ------ ------ ------ ------
 0.87   0.03   0.01    0     0.01   0.01 

 0.06   0.91   0.01    0     0.02   0.02 

 0.03   0.04   0.88   0.01   0.03   0.03 

  0      0      0     0.01   0.2    0.2  

 0.03   0.02   0.1    0.97   0.08   0.08 

  0      0      0      0     0.23   0.22 

  0      0      0     0.01   0.2    0.21 

  0      0      0      0     0.24   0.23 
------ ------ ------ ------ ------ ------

```r
# Run AMOPLS on the pure matrices
R2X_perPC <- c()
for (i in 1:npred) {
    res <- AMOPLS_function(EffectMatGLM, npred = i, ortho, X = Xmat, UseX = FALSE)
    R2X_perPC[i] <- res[[1]]$R2X[1]
}

R2X <- c(R2X_perPC[1], diff(R2X_perPC))
R2X <- round(R2X * 100, 2)

pander("R2X")
```

R2X

```r
names(R2X) <- c(paste0("PrC", 1:npred))
R2X
```

```
## PrC1 PrC2 PrC3 PrC4 
## 8.13 7.79 7.85 6.62
```

```r
# explained variation for Y-orthogonal model components.
pander("explained variation for Y-orthogonal model components. \n")
```

explained variation for Y-orthogonal model components. 

```r
R2XOcum <- res[[1]]$R2XO[ortho + 1]
R2XO <- diff(res[[1]]$R2XO)
R2XO <- round(R2XO * 100, 2)
names(R2XO) <- c(paste0("PrO", 1:ortho))
R2XO
```

```
##  PrO1  PrO2 
## 26.84 12.31
```

```r
R2X_perPC <- c(R2X, R2XO)
pander("R2X per component")
```

R2X per component

```r
names(R2X_perPC) <- c(paste0("PrC", 1:npred), paste0("PrO", 1:ortho))
R2X_perPC
```

```
##  PrC1  PrC2  PrC3  PrC4  PrO1  PrO2 
##  8.13  7.79  7.85  6.62 26.84 12.31
```

```r
R2Xtot <- sum(R2X_perPC)
pander("R2X tot")
```

R2X tot

```r
R2Xtot
```

```
## [1] 69.54
```

```r
### R2Yhat
R2Yhat_tot <- res[[1]]$R2Yhat[ortho + 1]
roundedR2Yhat <- round(R2Yhat_tot * 100, 2)
pander("R2Yhat")
```

R2Yhat

```r
roundedR2Yhat
```

```
## [1] 94.09
```



## Contribution plot


```r
saliences <- AMOPLS_lambda
colnames(saliences) <- c(paste0("PrC", 1:A), paste0("OC", 1:ortho))
rownames(saliences) <- ModelTerms_abbrev
saliences <- cbind(saliences[, A + (1:ortho)], saliences[, c(1:A)])

ylim <- c(0, 1)
angle1 <- rep(c(45, 45, 135), length.out = p)
angle2 <- rep(c(45, 135, 135), length.out = p)
density1 <- seq(5, 35, length.out = p)
density2 <- seq(5, 35, length.out = p)
col <- rainbow(p)

# plot --------------------------------------------------------------------
par(xpd = TRUE, mar = c(5, 4, 4, 10))
x <- barplot(saliences, beside = TRUE, ylim = ylim, xaxt = "n", col = col, angle = angle1, 
    density = density1, main = "AMOPLS effects contributions")
barplot(saliences, add = TRUE, beside = TRUE, ylim = ylim, xaxt = "n", col = col, 
    angle = angle2, density = density2)

labs <- as.vector(colnames(saliences))
text(cex = 1, x = colMeans(x) - 0.22, y = -0.08, labs, xpd = TRUE, srt = 45)

legend("topright", legend = ModelTerms_abbrev, fill = col, border = col, col = col, 
    angle = angle1, density = density1, inset = c(-0.3, 0), cex = 0.8)
par(bg = "transparent")
legend("topright", legend = ModelTerms_abbrev, col = col, fill = col, border = col, 
    angle = angle2, density = density2, inset = c(-0.3, 0), cex = 0.8)
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-60-1.png" width="800px" />


## Scores plots


```r
# Scores scatterplot matrix
scores <- cbind(AMOPLS_model$To, AMOPLS_model$T)
dimnames(scores)[[2]] <- c(paste0("OC ", 1:ortho), paste0("PrC ", 1:npred))

dimnames(scores)[[2]] <- paste0(dimnames(scores)[[2]], " \n", round(R2X_perPC[c((npred + 
    1):(npred + ortho), 1:npred)], 2), "%", c("", "", "\n \"Hippurate\" ", "\n \"Citrate\" ", 
    "\n \"Time\" ", "\n \"HxT\" "))

dimnames(scores)[[1]] <- rownames(outcomes)



# define textlabel and pch2 for spotPoints
textlabel <- matrix("", ncol = dim(scores)[2], nrow = dim(scores)[1], dimnames = list(rownames(scores), 
    NULL))
textlabel[rownames(scores) == c("M2C04D2R1"), ] <- rep(c("a"), dim(scores)[2])
textlabel[rownames(scores) == c("M1C04D2R2"), ] <- rep(c("b"), dim(scores)[2])
textlabel[rownames(scores) == c("M2C02D2R1"), ] <- rep(c("c"), dim(scores)[2])

pch2 <- rep("", dim(scores)[1])
pch2[rownames(scores) %in% c("M2C04D2R1", "M1C04D2R2", "M2C02D2R1")] <- "O"

# colnames(scores)
scores[, 3] <- scores[, 3] * -1
scores[, 4] <- scores[, 4] * -1

par(xpd = TRUE)
pairsBG(scores, titre = "AMOPLS Scores plots", oma = c(3, 3, 5, 15), spotPoints = FALSE)
legendsScatterMatrix()
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-61-1.png" width="800px" />


## Weighted Loadings


```r
###----------Weighted Loadings
scores <- cbind(AMOPLS_model$T, AMOPLS_model$To)
R <- dim(scores)[2]

matload <- matrix(nrow = R, ncol = m)  # loading matrix
dimnames(matload)[[2]] <- Variables
dimnames(matload)[[1]] <- c(paste("PrC", 1:npred), paste("OC", 1:ortho))

for (i in 1:R) {
    matload[i, ] <- scores[, i] %*% Reduce("+", Map("*", EffectMatGLM, sqrt(AMOPLS_lambda[, 
        i])))
}


# plot
LinePlot(as.matrix(matload), createWindow = FALSE, type = "s", main = "Loadings", 
    rows = 1:R, num.stacked = R, xaxis_type = c("numerical"), nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-62-1.png" width="800px" />

```r
# main effects loadings
matload[c(1), ] <- matload[c(1), ] * -1  # hippurate
matload[c(2), ] <- matload[c(2), ] * -1  # citrate

LinePlot(as.matrix(matload[c(1, 2, 3, 4), ]), createWindow = FALSE, main = "AMOPLS Loadings", 
    type = "s", num.stacked = 4, xaxis_type = c("numerical"), nxaxis = 10)
```

```
## [[1]]
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-62-2.png" width="800px" />

## True RSR


```r
# RSR
RSR_AMOPLS <- AMOPLS_lambda[length(ModelTerms), npred + 1]/AMOPLS_lambda[1:length(ModelTerms_sansE), 
    npred + 1]

pander("RSR")
```

RSR

```r
round(RSR_AMOPLS, 2)
```

```
## [1] 17.97 13.39  7.68  1.20  2.99  1.05  1.19
```


## Permutations BG3


```r
################################
# true F stats
################################

#Run AMOPLS on the pure matrices
res <- AMOPLS_function(EffectMatGLM,npred, ortho) 
AMOPLS_lambda <- res[[2]] #saliences, rows = matrices, cols = components

# RSR computation
A <- res[[1]]$A #Number of explanatory components
RSR_AMOPLS <- AMOPLS_lambda[length(ModelTerms),A+1]/AMOPLS_lambda[1:length(ModelTerms_sansE),A+1]
names(RSR_AMOPLS) <- ModelTerms_sansE

TRUE_RSR_stat  <- RSR_AMOPLS 
outcomes_true <- outcomes 

################################
# Permutations
################################


################################
#  main effects
################################

namesmainEffects <- ModelTerms_sansE[!grepl(":",ModelTerms_sansE)]
mainEffects <- design[,namesmainEffects]

RSR_stat <- vector("list", length(mainEffects))
names(RSR_stat) <- namesmainEffects

ptm <- proc.time()

for (k in namesmainEffects) { 
  for (i in 1:nperm) { 
    cat("\n k = ",k," i = ",i)

# contruction of the entry matrix for AMOPLS  
    EffectMatGLM_perm <-DecompEffectMatrixPermut[[k]][[i]][-1] # minus intercept

# Run AMOPLS on the pure matrices
    res<-AMOPLS_function(EffectMatGLM_perm, npred, ortho) 
    A <- res[[1]]$A #Number of explanatory components
    
# Store the results
    AMOPLS_lambda <- res[[2]] #saliences, rows = matrices, cols = components
    
    RSR_AMOPLS <- AMOPLS_lambda[length(ModelTerms),A+1]/AMOPLS_lambda[1:length(ModelTerms_sansE),A+1]
    names(RSR_AMOPLS) <- ModelTerms_sansE
# Ading the statistics to the vector   
    
    RSR_stat[[k]] <- c(RSR_stat[[k]], RSR_AMOPLS[[k]])
  }

}
proc.time() - ptm

mainEffects_RSR_stat <- RSR_stat

###################
#  Interactions
###################
nIntEff <- length(namesinteractionEffects)
IndIntEff <- nMainEff+1:nIntEff
RSR_stat <- vector("list", nIntEff)
names(RSR_stat) <- namesinteractionEffects

ptm <- proc.time()

for (i in 1:nperm) { 
    cat("\n i = ",i)

# Récupération des matrice des effets ec contruction input AMOPLS 
 EffectMatGLM_perm <- DecompEffectMatrixPermut[["Interaction"]][[i]][-1] # minus intercept
# Application de la fonction AMOPLS
res2 <- AMOPLS_function(EffectMatGLM_perm, npred, ortho) 
A <- res2[[1]]$A #Number of explanatory components

# Recuperation des lambda
AMOPLS_lambda <- res2[[2]] #saliences, rows = matrices, cols = components
# Calcul des RSR    
RSR_AMOPLS <- AMOPLS_lambda[dim(AMOPLS_lambda)[1],A+1]/AMOPLS_lambda[1:((dim(AMOPLS_lambda)[1]-1)),A+1]
names(RSR_AMOPLS) <- ModelTerms_sansE
for (ki in 1:nIntEff)
{k=namesinteractionEffects[ki]
RSR_stat[[k]] <- c(RSR_stat[[k]], RSR_AMOPLS[k])}
}
}

proc.time() - ptm

interEffects_RSR_stat <- RSR_stat
RSR_stat <- c(mainEffects_RSR_stat,interEffects_RSR_stat)
pval <- c()

for (i in ModelTerms_sansE){
pval[i] <- sum(TRUE_RSR_stat[i] <= RSR_stat[[i]])/nperm
}

names(pval) <- ModelTerms_sansE
pval
save(RSR_stat,pval, TRUE_RSR_stat,nperm, file = file.path(directoryOut, "PermutationResAMOPLS.RData"))
```


## Result Permut


```r
load(file = file.path(directoryOut, "PermutationResAMOPLS.RData"))
names(pval) <- names(RSR_stat) <- names(TRUE_RSR_stat) <- ModelTerms_sansE

pander("p-values")
```

p-values

```r
pander(pval)
```


--------------------------------------------------------------------------------
 Hippurate   Citrate   Time   Hippurate:Citrate   Hippurate:Time   Citrate:Time 
----------- --------- ------ ------------------- ---------------- --------------
     0          0       0           0.499               0              0.68     
--------------------------------------------------------------------------------

Table: Table continues below

 
------------------------
 Hippurate:Citrate:Time 
------------------------
         0.548          
------------------------

```r
for (i in ModelTerms_sansE) {
    hist(RSR_stat[[i]], xlim = range(TRUE_RSR_stat[i], RSR_stat[[i]]), breaks = 20, 
        main = i)
    points(TRUE_RSR_stat[i], 0, col = "red", pch = 19)
}
```

<img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-1.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-2.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-3.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-4.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-5.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-6.png" width="50%" /><img src="Code_article_Guisset_et_al_files/figure-html/unnamed-chunk-64-7.png" width="50%" />


# All screeplots


```r
SC_AComDim <- pcvar  # AComDim
SC_AMOPLS <- R2Yhat_perPC  # AMOPLS
SC_PARAFASCA <- Rsq  # PARAFASCA

par(mfrow = c(1, 3), mar = c(3, 4, 2, 2))

# PARAFASCA
plot(c(0, 1:nparam), c(0, Rsq), type = "b", ylab = "R2", yaxt = "n", xlab = "n comp.", 
    main = "(a) PARAFASCA Scree plot", mgp = c(2, 1, 0), ylim = c(0, 1))
axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
points(c(0, 1:nparam), c(0, Rsq), col = "blue", pch = 20)

# ACOMDIM
plot(c(0, 1:nparam), c(0, pcvar), type = "b", yaxt = "n", ylab = "", xlab = "n CC", 
    main = "(b) AComDim Scree plot", mgp = c(2, 1, 0), ylim = c(0, 1))
axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
points(c(0, 1:nparam), c(0, pcvar), col = "blue", pch = 20)

# AMOPLS
plot(c(0, 1:nparam), c(0, SC_AMOPLS), type = "b", yaxt = "n", ylab = "", xlab = "n PrC", 
    main = "(c) AMOPLS Scree plot", mgp = c(2, 1, 0), ylim = c(0, 1))
axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
points(c(0, 1:nparam), c(0, R2Yhat_perPC), col = "blue", pch = 20)
```

# References

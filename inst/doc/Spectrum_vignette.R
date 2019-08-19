## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test1 <- Spectrum(blobs,showpca=TRUE,fontsize=8,dotsize=2)

## ------------------------------------------------------------------------
names(test1)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
RNAseq <- brain[[1]]
test2 <- Spectrum(RNAseq,showdimred=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
RNAseq <- brain[[1]]
test3 <- Spectrum(RNAseq,showres=FALSE,runrange=TRUE,krangemax=10)

## ------------------------------------------------------------------------
head(test3[[2]]$assignments)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test4 <- Spectrum(brain,showdimred=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
brain1 <- brain[[1]]
brain2 <- brain[[2]]
brain3 <- brain[[3]]
brain1 <- brain1[,-5:-10]
brain_m <- list(brain1,brain2,brain3)
test4 <- Spectrum(brain_m,missing=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test5 <- Spectrum(circles,showpca=TRUE,method=2,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test6 <- Spectrum(spirals,showpca=TRUE,method=2,tunekernel=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test7 <- Spectrum(blobs,FASP=TRUE,FASPk=300,fontsize=8,dotsize=2)

## ------------------------------------------------------------------------
names(test7)

## ------------------------------------------------------------------------
head(test7[[1]])

## ------------------------------------------------------------------------
## 1. run my clustering algorithm yielding assignments in vector, e.g. 1,2,2,1,2,2...
## 2. reorder data according to assignments
ind <- sort(as.vector(test2$assignments),index.return=TRUE)
datax <- RNAseq[,ind$ix] ## order the original data 
#annonx <- meta[ind$ix,] ## order the meta data
#annonx$cluster <- ind$x ## add the cluster to the meta data
## 3. do heatmap 
# insert your favourite heatmap function


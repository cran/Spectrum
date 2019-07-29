## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test1 <- Spectrum(blobs,showdimred=TRUE,fontsize=8,dotsize=2)

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


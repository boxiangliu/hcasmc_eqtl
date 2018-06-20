#!/usr/bin/Rscript

args <- commandArgs(TRUE);
DAT <- args[1]
NAME <- args[2]
CORE <- args[3]
SIGNATURE <- args[4]

library(xCell)
## load manually edited xCell function to set number of cores
source("xCell_parallel.R")

dat.gs <- read.table(DAT, head=T, as.is=T, row.names=2) # keep ensemlID but read gene symbol as rownames
dim(dat.gs)

## read custom signature file
load(SIGNATURE) # name is gs.p

## run xCell
res <- xCellAnalysis(dat.gs[,-1], signature=gs.p, return.raw=TRUE, microenv=F, combine.celltype=F, save.raw = T, file.name = NAME, parallel.sz=as.numeric(CORE)) 
res.t <- as.data.frame(t(res))
save(res.t, file=paste0(NAME, ".Rda")) # redundant as save.raw will save the same data as txt

library(coloc)
library(snpStats)


setClass("simdata",
     representation(df1="data.frame",df2="data.frame"))
setValidity("simdata", function(object) {
  n <- nrow(object@df1)
  if(nrow(object@df2)!=n)
    return("nrow of '@df1' should equal nrow of '@df2'")
})
setMethod("show", signature="simdata", function(object) {
  cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n")
})

sim.data <- function(nsnps=100,nsamples=200,causals=c(1,50),nsim=1) {
  cat("Generate",nsim,"small sets of data\n")
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y2 <- rnorm(nsamples,rowSums(X2[,causals]),2)
  colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
  df1 <- cbind(Y=Y1,X1)
  df2 <- cbind(Y=Y2,X2)
  if(nsim==1) {
    return(new("simdata",
           df1=as.data.frame(df1),
           df2=as.data.frame(df2)))
  } else {
    index <- split(1:(nsamples * nsim), rep(1:nsim, nsamples))
    objects <- lapply(index, function(i) new("simdata", df1=as.data.frame(df1[i,]),
                         df2=as.data.frame(df2[i,])))
    return(objects)
  }
}

## simulate some data and load the coloc library
set.seed(12345)
data <- sim.data(nsamples=1000,nsim=1)
library(coloc)


# run coloc:
library(snpStats)

Y1 <- data@df1$Y
Y2 <- data@df2$Y

X1 <- new("SnpMatrix",as.matrix(data@df1[,-1]))
X2 <- new("SnpMatrix",as.matrix(data@df2[,-1]))

p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
p2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
maf <- col.summary(X2)[,"MAF"]

my.res <- coloc.abf(dataset1=list(pvalues=p1,N=nrow(X1),type="quant"),
            dataset2=list(pvalues=p2,N=nrow(X2),type="quant"),
            MAF=maf)
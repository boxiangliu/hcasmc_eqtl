install.packages('GenOrd')
library(GenOrd)

set.seed(1)
# Sets the marginals.
# The values are cumulative so for the first variable the first marginal will be .1, the second is .2, the third is .3, and the fourth is .4
marginal <- list(c(0.25,0.75),c(0.25,0.75))
# Checks the lower and upper bounds of the correlation coefficients.
corrcheck(marginal)
# Sets the correlation coefficients
R <- matrix(c(1,0.2,0.2,1),2,2) # Correlation matrix
n <- 100
##Selects and ordinal sample with given correlation R and given marginals.
m <- ordsample(n, marginal, R)
##compare it with the pre-defined R
cor(m)
table(m[,1],m[,2])
m=data.frame(m)
m
ggplot(m,aes(x=X1-1,y=X2-1))+geom_jitter(width=0.3,height=0.3)
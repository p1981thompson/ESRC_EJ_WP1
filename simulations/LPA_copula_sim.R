

library(copula)
library(KScorrect)


#=========================================================================#

#Define copula
#=========================================================================#
n_ppts <- 1000

set.seed(100)
# constructs an elliptical copula (this has inputs of the correlation matrix - overall, pooled estimates not classones.)
myCop <- normalCopula(param=c(0.5,0.4,0.4,0.6,
                                  0.6,0.5,0.2,
                                      0.7,0.2,
                                          0.2), 
                      dim=5, dispstr="un")

# creates a multivariate distribution via copula (each line below is the marginal distribution for each variable, with means and sds for each profile)
myMvd <- mvdc(copula=myCop, margins=c("mixnorm", "mixnorm", "mixnorm", "mixnorm", "mixnorm"),
              paramMargins=list(list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), 
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), 
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), 
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), 
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2))))

#Generate data based on the copula.
z2 <- rMvdc(n_ppts, myMvd)
colnames(z2) <- c("nara_comp", "nara_acc", "word_acc", "nonword_acc", "wold_comp")

#plot to see if we get something sensible.
pairs.panels(z2)



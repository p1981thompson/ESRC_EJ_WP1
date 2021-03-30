#==============================================================================================#
# Alternate data simulation for LPA model using copula    
#==============================================================================================#

#17-11-2020

#load R packages
library(copula)
library(KScorrect)


#==============================================================================================#

#Define copula: links the marginal distributions together to form the joint distribution
#==============================================================================================#
n_ppts <- 1000

set.seed(100)
# constructs a copula (this has inputs of the correlation matrix - overall, pooled estimates not class ones.) We use a normal or gaussian copula here but there are others. I just seeing whether we might get better results using a different copula. 
myCop <- normalCopula(param=c(0.5,0.4,0.4,0.6, # param is the correlation matrix between the marginal variables
                                  0.6,0.5,0.2,
                                      0.7,0.2,
                                          0.2), 
                      dim=5, dispstr="un") #dim is number of marginal variables; dispstr is the structure of the correlation matrix, so in our case we have unstructured as the variables are correlated but not structured unlike if we had repeated measures.
#==============================================================================================#

# creates a multivariate distribution using the copula (link) - each line below is the marginal distribution for each variable, with means and sds for each profile.

# We define each individual variable's distribution within the multivariate distribution. So each variable has its own distribution, in our case we want a mixture normal distribution; see here: https://en.wikipedia.org/wiki/Multimodal_distribution#Mixture_of_two_normal_distributions

# The mvdc function sets up the multivariate distribution, specifying the copula (link), then the marginals and finally the parameters for setting the location (mean) and spread (sd) and mixing proportions of the marginal distributions.

myMvd <- mvdc(copula=myCop, margins=c("mixnorm", "mixnorm", "mixnorm", "mixnorm", "mixnorm"), #These are the mixture distributions for each variable.
              paramMargins=list(list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), # NARA - comprehension
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), # NARA - accuracy
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), # Word reading
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)), # Nonword reading
                                list(mean=c(105,90,80),  sd=c(5,7,10),pro=c(0.6,0.2,0.2)))) # WOLD listening

#Generate data based on the multivariate distribution. Same as if you did rnorm(), you are saying that you want to take X number of samples from that distribution.
simdat2 <- rMvdc(n_ppts, myMvd)
colnames(simdat2) <- c("nara_comp", "nara_acc", "word_acc", "nonword_acc", "wold_comp")

#==============================================================================================#
#plot to see if we get something sensible.
pairs.panels(simdat2)

#==============================================================================================#
#We can run a basic version of your script below (no missing data). It still does not pick up the structure, so I'm thinking on that. 

sim_dat2 %>% 
  estimate_profiles(1:10) %>% 
  compare_solutions(statistics = c("AIC", "BIC", "BLRT_val", "BLRT_p"))

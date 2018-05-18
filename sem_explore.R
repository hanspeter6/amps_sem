# libraries
library(tidyverse)
library(nFactors)
library(psych)
library(sem)

# reading in the mins (nationals) datasets:
set02_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set02_min.rds")
# set05_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set05_min.rds")
set08_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set08_min.rds")
set10_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set10_min.rds")
set12_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set12_min.rds")
set14_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set14_min.rds")

# try to create single frame:
# function to add year as variable
add_year <- function(set, year) {
        set %>%
                mutate(year = year) %>%
                dplyr::select(qn, pwgt, year, everything())
}

# combine in single dataset, excluding lifestages (since NA in 2008),
#  lifestyle and attitudes (since missing in 2002); 
#  2005 excluded since has only one value for internet
set_min <- rbind.data.frame(add_year(set02_min[,-11], 2002),
                            # add_year(set05_min, 2005),
                            add_year(set08_min[,-c(11,14,15)], 2008),
                            add_year(set10_min[,-c(11,14,15)], 2010),
                            add_year(set12_min[,-c(11,14,15)], 2012),
                            add_year(set14_min[,-c(11,14,15)], 2014))

# pca method (package "nFactors") to estimate optimum number of factors to extract:
## Determine Number of Factors to Extract
ev <- eigen(cor(set_min[,20:ncol(set_min)]))
ap <- parallel(subject=nrow(set_min[,20:ncol(set_min)]),var=ncol(set_min[,20:ncol(set_min)]),
               rep=100,cent=.02)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
# jpeg("nScree_10_min")
plotnScree(nS, main = "All Sets") # optimal = 7

# dev.off()

# extract factors (using pa method and common factor extraction)
factors_set_min_pa <- fa(r = set_min[,20:ncol(set_min)], nfactors = 7, fm = "pa")
factors_set_min_ml <- fa(r = set_min[,20:ncol(set_min)], nfactors = 7, fm = "ml")

# consider loadings:
write.csv(round(factors_set_min_pa$loadings, 3), "loadings_pa.csv")
write.csv(round(factors_set_min_ml$loadings, 3), "loadings_ml.csv")

# consider a first SEM model (maybe confirmatory???)
# Using structure as shown for ml. Interstingly,very little difference between pa and ml


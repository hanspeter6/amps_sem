# libraries
library(tidyverse)
library(nFactors)
library(psych)
library(sem)
# library(ctsem)

# reading in the mins (nationals) datasets:
set02_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set02_min.rds")
# set05_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set05_min.rds")
set08_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set08_min.rds")
set10_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set10_min.rds")
set12_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set12_min.rds")
set14_min <- readRDS("/Users/HansPeter/Dropbox/Statistics/UCTDataScience/Thesis/amps_nationals/set14_min.rds")

# saving them in this directory for sharing:
save(set02_min, set08_min, set10_min, set12_min, set14_min, file = "min_sets.RData")

load("min_sets.RData") # NB for HP: be careful to first read, since details may have changed in data processing

# create single frame:
# function to add year as variable
add_year <- function(set, year) {
        set %>%
                mutate(year = year) %>%
                dplyr::select(qn, pwgt, year, everything())
}

# combine in single dataset:
# -> excluding 'lifestages' (since NA in 2008),
#  -> excluding 'lifestyle' and 'attitudes' (since missing in 2002); 
#  -> already excluded 2005 since has only one value for internet engagement
set_min <- rbind.data.frame(add_year(set02_min[,-11], 2002),
                            add_year(set08_min[,-c(11,14,15)], 2008),
                            add_year(set10_min[,-c(11,14,15)], 2010),
                            add_year(set12_min[,-c(11,14,15)], 2012),
                            add_year(set14_min[,-c(11,14,15)], 2014))

# scaling the media vehicle variables to mean = 0 and sd = 1:
set_min <- cbind.data.frame(set_min[,1:19], scale(set_min[,20:ncol(set_min)]))

## For now consider entire dataset...
## 
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

### consider SEM for individual years:
## first only 2010, and only confirmatory model
# 1st Consider correlations of media vehicles set10_min:
media_vehicles_10 <- set10_min[,22:ncol(set10_min)]

cov_mat_10 <- cov(media_vehicles_10)

# checking covariance matrix is invertible
f <- function(m) class(try(solve(m),silent=T)) == "matrix"
f(cov_mat_10) # TRUE....so that not the problem

# prepare text file 

# try first basic confirmatory model, unit variances of factors ("basic_measurement_model.txt")
basic_model_10 <- specifyModel(file = "basic_measurement_model.txt")
basic_sem_10 <- sem(model = basic_model_10,
                 S = cov_mat_10,
                 # data = media_vehicles_10,
                 N = 25160)

summary(basic_sem_10)
# pathDiagram(basic_sem_10) # cant get this to work...

# adding between factor correlations, rest the same ("factor_correlations_model")
factor_corr_model_10 <- specifyModel(file = "factor_correlations_model.txt")
factor_corr_sem_10 <- sem(model = factor_corr_model_10,
                    S = cov_mat_10,
                    data = set10_min,
                    N = 25160)

summary(factor_corr_sem_10, digits = 3)

pathDiagram(factor_corr_sem_10) # not working
# would like to add demographics to factors.... not sure how to do that..
# also, now what, copare sem for different years??? continuous time wont work...

# Try lavaan package: can deal with ordinal and bivariate (dummies)
library(lavaan)

# basic model cfa model
myModel1 <- '# latent variable definitions:
popPrint =~ Business.Day + Mail.n.Guardian + The.Sunday.Independent + Sunday.Times + You + Car + Cosmopolitan + Getaway + Topcar
afrikaans =~ Rapport + Huisgenoot + Sarie
soccer =~ Soccer.Laduma + Kickoff
african =~ Drum + Bona + Metro.FM
social =~ X5FM + DSTV + int_social + int_radio + int_search
freeTV =~ e.tv + SABC.1 + SABC.2 + SABC.3
news =~ int_print + int_news
# variances and covariances of factors
# variances
popPrint ~~ popPrint
african ~~ african
social ~~ social
afrikaans ~~ afrikaans
soccer ~~ soccer
news ~~ news
freeTV ~~ freeTV
#covariances
freeTV ~~ social
freeTV ~~ african
freeTV ~~ afrikaans
freeTV ~~ news
freeTV ~~ popPrint
freeTV ~~ soccer
social ~~ soccer
social ~~ popPrint
social ~~ news
social ~~ african
social ~~ afrikaans
african ~~ afrikaans
african ~~ soccer
african ~~ news
african ~~ popPrint
afrikaans ~~ popPrint
afrikaans ~~ news
afrikaans ~~ soccer
news ~~ soccer
news ~~ popPrint
popPrint ~~ soccer
'
fit1 <- lavaan::sem(myModel1, data = set10_min, fit.measures = TRUE)
summary(fit1, standardized = TRUE)

# try to add some demographics as regression formulae:
# # duplicate set and turn ordered into numerics
set10_min_cats <- set10_min
set10_min_cats$age <- as.numeric(set10_min_cats$age)
set10_min_cats$edu <- as.numeric(set10_min_cats$edu)
set10_min_cats$hh_inc <- as.numeric(set10_min_cats$hh_inc)
set10_min_cats$lsm <- as.numeric(set10_min_cats$lsm)

# turn sex and race into dummies
library(dummies)
set10_min_dummies <- dummy.data.frame(data = set10_min_cats, names = c("sex", "race"))
names(set10_min_dummies)[4] <- "male"
names(set10_min_dummies)[5]  <- "female"
names(set10_min_dummies)[8]  <- "black"
names(set10_min_dummies)[9]  <- "coloured"
names(set10_min_dummies)[10]  <- "indian"
names(set10_min_dummies)[11]  <- "white"

fullModel <- '# latent variable definitions:
popPrint =~ Business.Day + Mail.n.Guardian + The.Sunday.Independent + Sunday.Times + You + Car + Cosmopolitan + Getaway + Topcar
afrikaans =~ Rapport + Huisgenoot + Sarie
soccer =~ Soccer.Laduma + Kickoff
african =~ Drum + Bona + Metro.FM
social =~ X5FM + DSTV + int_social + int_radio + int_search
freeTV =~ e.tv + SABC.1 + SABC.2 + SABC.3
news =~ int_print + int_news
# variances and covariances of factors
# variances
popPrint ~~ popPrint
african ~~ african
social ~~ social
afrikaans ~~ afrikaans
soccer ~~ soccer
news ~~ news
freeTV ~~ freeTV
#covariances
freeTV ~~ social
freeTV ~~ african
freeTV ~~ afrikaans
freeTV ~~ news
freeTV ~~ popPrint
freeTV ~~ soccer
social ~~ soccer
social ~~ popPrint
social ~~ news
social ~~ african
social ~~ afrikaans
african ~~ afrikaans
african ~~ soccer
african ~~ news
african ~~ popPrint
afrikaans ~~ popPrint
afrikaans ~~ news
afrikaans ~~ soccer
news ~~ soccer
news ~~ popPrint
popPrint ~~ soccer
# regressions
popPrint ~ age  + edu + hh_inc + lsm + male + black
african ~ age  + edu + hh_inc + lsm + male + black
afrikaans ~ age  + edu + hh_inc + lsm + male + black
social ~ age  + edu + hh_inc + lsm + male + black
news ~ age  + edu + hh_inc + lsm + male + black
freeTV ~ age  + edu + hh_inc + lsm + male + black
soccer ~ age  + edu + hh_inc + lsm + male + black
'
fit2 <- lavaan::sem(fullModel, data = set10_min_dummies, fit.measures = TRUE)
summary(fit2, standardized = TRUE)

# worry about distributions (especially binary, so se not accurate:)
fit2_boot <- bootstrapLavaan(object = fit2)

# try another year and compare parameter values


# for longitudinal stuff: options: go back to what I did before; use these new SEM parameters in some way
# try working with single value per category level per year and use longitudinal SEM ctsem
# 
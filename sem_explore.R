# libraries
library(tidyverse)
library(nFactors)
library(psych)
library(lavaan)
library(dummies)
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
# 1st Consider correlations of media vehicles":

cov_mat_10 <- cov(set10_min[,22:ncol(set10_min)])

# checking covariance matrix is invertible
f <- function(m) class(try(solve(m),silent=T)) == "matrix"
f(cov_mat_10) # TRUE....so that not the problem

# function to prep sets into numeric for ordered variables and dummies otherwise
prep_sets_sem <- function(set) {
        
        # ordered to numeric
        set$age <- as.numeric(set$age)
        set$edu <- as.numeric(set$edu)
        set$hh_inc <- as.numeric(set$hh_inc)
        set$lsm <- as.numeric(set$lsm)
        
        # 'sex' and 'race' to dummies
        require(dummies)
        set2 <- dummy.data.frame(data = set, names = c("sex", "race"))
        names(set2)[which(names(set2) == 'sex1')] <- "male"
        names(set2)[which(names(set2) == 'sex2')]  <- "female"
        names(set2)[which(names(set2) == 'race1')]  <- "black"
        names(set2)[which(names(set2) == 'race2')]  <- "coloured"
        names(set2)[which(names(set2) == 'race3')]  <- "indian"
        names(set2)[which(names(set2) == 'race4')]  <- "white"
        return(set2)
        
}

# prep sets
set02_sem <- prep_sets_sem(set02_min)
set08_sem <- prep_sets_sem(set08_min)
set10_sem <- prep_sets_sem(set10_min)
set12_sem <- prep_sets_sem(set12_min)
set14_sem <- prep_sets_sem(set14_min)

# define 'lavaan' SEM model:
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
popPrint ~ age  + edu + hh_inc + lsm + male + black + coloured + white
african ~ age  + edu + hh_inc + lsm + male + black + coloured + white
afrikaans ~ age  + edu + hh_inc + lsm + male + black + coloured + white
social ~ age  + edu + hh_inc + lsm + male + black + coloured + white
news ~ age  + edu + hh_inc + lsm + male + black + coloured + white
freeTV ~ age  + edu + hh_inc + lsm + male + black + coloured + white
soccer ~ age  + edu + hh_inc + lsm + male + black + coloured + white
'

# fit full model to each year seperately
fit_sem_02 <- lavaan::sem(fullModel, data = set02_sem, fit.measures = TRUE)
# summary(fit_sem_02, standardized = TRUE)
fit_sem_08 <- lavaan::sem(fullModel, data = set08_sem, fit.measures = TRUE)
# summary(fit_sem_08, standardized = TRUE)
fit_sem_10 <- lavaan::sem(fullModel, data = set10_sem, fit.measures = TRUE)
# summary(fit_sem_10, standardized = TRUE)
fit_sem_12 <- lavaan::sem(fullModel, data = set12_sem, fit.measures = TRUE)
# summary(fit_sem_12, standardized = TRUE)
fit_sem_14 <- lavaan::sem(fullModel, data = set14_sem, fit.measures = TRUE)
# summary(fit_sem_14, standardized = TRUE)

# fitMeasures(fit_sem_10)

# write parameters to .csv
write.csv(parameterEstimates(fit_sem_02), file = "paraEst_02.csv")
write.csv(parameterEstimates(fit_sem_08), file = "paraEst_08.csv")
write.csv(parameterEstimates(fit_sem_10), file = "paraEst_10.csv")
write.csv(parameterEstimates(fit_sem_12), file = "paraEst_12.csv")
write.csv(parameterEstimates(fit_sem_14), file = "paraEst_14.csv")

# # worry about distributions (especially binary, so se not accurate:)
# fit2_boot <- bootstrapLavaan(object = fit2) # GREAT, BUT TAKES TIME to process

# save objects:
save(fit_sem_02, fit_sem_08, fit_sem_10, fit_sem_12, fit_sem_14, fit2_boot, file = "fit_objects.RData")
# try another year and compare parameter values
load("fit_objects.RData")

# consider correlation matrix between different years' parameter estimates
est_matrix <- cbind(parameterEstimates(fit_sem_02)[,4],
                    parameterEstimates(fit_sem_08)[,4],
                    parameterEstimates(fit_sem_10)[,4],
                    parameterEstimates(fit_sem_12)[,4],
                    parameterEstimates(fit_sem_14)[,4])
colnames(est_matrix) <- c("2002", "2008", "2010", "2012", "2014")
round(cor(est_matrix), 2) # strong correlation of estimate values

# Options for longitudinal stuff:
# 1:  Use these SEM parameters in some way. Eg.,  model the parameter values over time...?
# 2:  Try working with single value per category level per year again (pseudo panels) and use longitudinal SEM ctsem...?
# 3:  Try suggestion from Francesca: use full set and model "year" as a categorical fixed effect in the structural model

# Option 3:use full set and model "year" as a categorical fixed effect in the structural model

# dataset already created: "set_min"
# will set up as before, except also including "year" as ordered - numeric, to limit the number of dummies
set_min_prepped <- prep_sets_sem(set_min)

# changing "year" to start from 0
set_min_prepped$year <- set_min_prepped$year - 2002

# including "year" in structural part:
fullModel_all <- '# latent variable definitions:
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
popPrint ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
african ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
afrikaans ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
social ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
news ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
freeTV ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
soccer ~ year + age + edu + hh_inc + lsm + male + black + coloured + white
'

# fit sem:
fit_sem_all <- lavaan::sem(fullModel_all, data = set_min_prepped, fit.measures = TRUE)
summary(fit_sem_all, standardized = TRUE) 

# not too sure this is best.... 'year' kind of gets lost in the other regression coeff...?? Difficult to see change over time
# --- maybe more dummies?.... and also consider bootstrapping options... for defensible std errs...,
#  will take lots of processing time given large all-set and checking one year seemed to give very similar estimates. Maybe ML OK for large sample...

# Option 2: pseudo panels + longitudinal SEM R package: "ctsem":

# define factor levels
set_min$age <- factor(set_min$age, labels = c("15-24","25-44", "45-54","55+"), ordered = TRUE)
set_min$race <- factor(set_min$race,labels = c("black", "coloured", "indian", "white"), ordered = TRUE)
set_min$edu <- factor(set_min$edu, labels = c("<matric", "matric",">matric" ) ,ordered = TRUE)
set_min$lsm <- factor(set_min$lsm, labels = c("LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10"), ordered = TRUE) #"LSM1-2", 
set_min$sex <- factor(set_min$sex, labels = c("male", "female"), ordered = TRUE)
set_min$hh_inc <- factor(set_min$hh_inc, labels = c("<R2500","R2500-R6999","R7000-R11999",">=R12000"), ordered = TRUE) # NB 2012 levels

# function to calculate se
se <- function(x) sqrt(var(x)/length(x))

# getting dataset into shape
set_fin <- set_min %>%
        gather(key = type, value = category, age, sex, edu, hh_inc, race, lsm) %>%
        dplyr::select(year, category, Business.Day:int_search) %>%
        as_tibble() %>%
        group_by(year, category) %>%
        summarise_all(c("mean"))

# setting year as absolute time, and changing variable to "time"
set_fin$year <- as.numeric(as.character(set_fin$year)) - 2002
names(set_fin)[1] <- "time"
# setting id column as required by the package: named id and numeric
id <- rep(1:22, 5)
# extracting "category" column for later reference:
category <- set_fin$category
set_fin <- set_fin[,-2]

set_fin_fin <- as.data.frame(set_fin) %>%
        mutate(id = id) %>%
        dplyr::select(time, id, everything())

# variable names cant have fullstops, change them to underscores:
library(stringr)
media_vehicles <- names(set_min[,20:ncol(set_min)])
nu_media_vehicles <- str_replace_all(media_vehicles, "\\.", "_")

names(set_fin_fin)[3:ncol(set_fin_fin)] <- nu_media_vehicles

# setting up model for ctsem:
library(ctsem)
mod_mat <- matrix(c(0,0,0,0,1,0,0, 
                    0,0,0,0,1,0,0,
                    0,0,0,1,0,0,0,
                    0,0,0,0,1,0,0,
                    0,0,0,0,1,0,0,
                    0,0,0,0,0,1,0,
                    0,0,0,0,0,0,1,
                    0,0,0,1,0,0,0,
                    0,0,0,0,1,0,0,
                    0,0,0,0,0,1,0,
                    0,0,0,0,0,0,1,
                    0,0,0,0,1,0,0,
                    0,0,0,0,1,0,0,
                    0,0,0,0,1,0,0,
                    0,0,0,1,0,0,0,
                    0,0,0,0,1,0,0,
                    0,1,0,0,0,0,0,
                    0,0,0,0,0,0,1,
                    1,0,0,0,0,0,0,
                    1,0,0,0,0,0,0,
                    1,0,0,0,0,0,0,
                    1,0,0,0,0,0,0,
                    0,1,0,0,0,0,0,
                    0,0,1,0,0,0,0,
                    0,1,0,0,0,0,0,
                    0,0,1,0,0,0,0,
                    0,1,0,0,0,0,0,
                    0,1,0,0,0,0,0), nrow = 28, ncol = 7, byrow = TRUE)

mod_ctsem <- ctModel(n.manifest = 28,
                     n.latent = 7,
                     LAMBDA = mod_mat,
                     Tpoints = 5,
                     manifestNames = nu_media_vehicles)

fit_ctsem <- ctFit(dat = set_fin_fin, mod_ctsem, dataform = "long")

fit_ctsem


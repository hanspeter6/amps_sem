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
        
        # 'sex', 'race', 'year' to dummies
        require(dummies)
        set2 <- dummy.data.frame(data = set, names = c("sex", "race", "year"))
        
        # change names
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

# # adding interaction effects as variables: CAN DO THIS WITH OPERATOR ':'
# set_min_prepped <- set_min_prepped %>%
#         mutate(year_age = year * age,
#                year_male = year * male,
#                year_edu = year * edu,
#                year_hh_inc = year * hh_inc,
#                year_white = year * white,
#                year_coloured = year * coloured,
#                year_indian = year * indian,
#                year_lsm = year * lsm)

# including "year" in structural part:
fullModel_all <- '

# latent variable definitions:
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
popPrint ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
african ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
afrikaans ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
social ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
news ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
freeTV ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
soccer ~ year2008 + year2010 + year2012 + year2014 + age + edu + hh_inc + lsm + male + white + coloured + indian + year2008:age + year2008:edu + year2008:hh_inc + year2008:lsm + year2008:male + year2008:white + year2008:coloured + year2008:indian + year2010:age + year2010:edu + year2010:hh_inc + year2010:lsm + year2010:male + year2010:white + year2010:coloured + year2010:indian + year2012:age + year2012:edu + year2012:hh_inc + year2012:lsm + year2012:male + year2012:white + year2012:coloured + year2012:indian + year2014:age + year2014:edu + year2014:hh_inc + year2014:lsm + year2014:male + year2014:white + year2014:coloured + year2014:indian
'


# fit sem:
fit_sem_all <- lavaan::sem(fullModel_all, data = set_min_prepped, fit.measures = TRUE)
summary(fit_sem_all, standardized = TRUE, fit.measures = TRUE) 

##  try to do plots comparing fitted with actual by category and by year:

# define the model matrix (NB order of columns: freeTV, social, news, afrikaans, popPrint, soccer, african)
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

# construct dataframe of standardised observed by the model:
obs <-as.matrix(set_min[,20:ncol(set_min)]) %*% mod_mat
colnames(obs) <- c("freeTV", "social", "news", "afrikaans", "popPrint","soccer", "african")
comb <- cbind.data.frame(set_min[1:19],scale(obs))

# function to create frames
frames_factors_all <- function(set, category) {
        require(dplyr)
        
        # set$cluster <- factor(set$cluster, labels = c("cluster1", "cluster2", "cluster3", "cluster4"))
        set$age <- factor(set$age, labels = c("15-24","25-44", "45-54","55+"), ordered = TRUE)
        set$race <- factor(set$race,labels = c("black", "coloured", "indian", "white"), ordered = TRUE)
        set$edu <- factor(set$edu, labels = c("<matric", "matric",">matric" ) ,ordered = TRUE)
        set$lsm <- factor(set$lsm, labels = c("LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10"), ordered = TRUE) #"LSM1-2", 
        set$sex <- factor(set$sex, labels = c("male", "female"), ordered = TRUE)
        set$hh_inc <- factor(set$hh_inc, labels = c("<R2500","R2500-R6999","R7000-R11999",">=R12000"), ordered = TRUE) # NB 2012 levels
        
        
        set %>%
                group_by_(year = "year", category = category) %>%
                summarise(popPrint = mean(popPrint),
                          afrikaans = mean(afrikaans),
                          soccer = mean(soccer),
                          african = mean(african),
                          social = mean(social),
                          freeTV = mean(freeTV),
                          news = mean(news)
                          # up_f1 = mean(ML1) + (2 * sd(ML1)/sqrt(length(ML1))),
                          # low_f1 = mean(ML1) - (2 * sd(ML1)/sqrt(length(ML1))),
                          # up_f2 = mean(ML2) + (2 * sd(ML2)/sqrt(length(ML2))),
                          # low_f2 = mean(ML2) - (2 * sd(ML2)/sqrt(length(ML2))),
                          # up_f3 = mean(ML3) + (2 * sd(ML3)/sqrt(length(ML3))),
                          # low_f3 = mean(ML3) - (2 * sd(ML3)/sqrt(length(ML3))),
                          # up_f4 = mean(ML4) + (2 * sd(ML4)/sqrt(length(ML4))),
                          # low_f4 = mean(ML4) - (2 * sd(ML4)/sqrt(length(ML4))),
                          # up_f5 = mean(ML5) + (2 * sd(ML5)/sqrt(length(ML5))),
                          # low_f5 = mean(ML5) - (2 * sd(ML5)/sqrt(length(ML5))),
                          # up_f6 = mean(ML6) + (2 * sd(ML6)/sqrt(length(ML6))),
                          # low_f6 = mean(ML6) - (2 * sd(ML6)/sqrt(length(ML6))),
                          # up_f7 = mean(ML7) + (2 * sd(ML7)/sqrt(length(ML7))),
                          # low_f7 = mean(ML7) - (2 * sd(ML7)/sqrt(length(ML7)))
                          )
        
}

# function to bind the frames by year
frame_bind_factor_all <- function(set) {
        rbind.data.frame(#frames_factors(set,"cluster"),
                frames_factors_all(set,"sex"),
                frames_factors_all(set,"age"),
                frames_factors_all(set,"edu"),
                frames_factors_all(set,"race"),
                frames_factors_all(set, "hh_inc"),
                frames_factors_all(set,"lsm")) %>% 
                dplyr::select(year,category, everything())
        
}
factor_obs_all <- data.frame(frame_bind_factor_all(comb))

all_plots_factors <- function(data, title = "All Scores Observed") {
        ggplot(data = data, title = title) +
                geom_line(aes(year, social, group = category, colour = "Social")) +
                geom_line(aes(year, freeTV, group = category, colour = "Free TV")) +
                geom_line(aes(year, afrikaans, group = category, colour = "Afrikaans")) +
                geom_line(aes(year, soccer, group = category, colour = "Soccer")) +
                geom_line(aes(year, news, group = category, colour = "News")) +
                geom_line(aes(year, african, group = category, colour = "African")) +
                geom_line(aes(year, popPrint, group = category, colour = "PopPrint")) +
                scale_colour_discrete(name="Factors") +
                facet_grid(. ~ category) +
                theme(axis.text.x = element_text(size = 6)) +
                labs(y = "aggregate scores observed", title = title)
        
}

all_plots_factors(factor_obs_all)

vector_row1 <- c("male", "female","15-24","25-44", "45-54","55+","black", "coloured", "indian", "white")
vector_row2 <- c("<matric", "matric",">matric", "<R2500","R2500-R6999","R7000-R11999",">=R12000", "LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10")

# now focussing on the predicted (fitted values):
# predicted values from sem model
pred_sem_all <- predict(fit_sem_all)

# first create similar from to observed values:
first <- cbind.data.frame(set_min[1:19],scale(pred_sem_all))

second <- frame_bind_factor_all(first)

colnames(second) <- c("year",
                      "category",
                      "popPrint_pred",
                      "afrikaans_pred",
                      "soccer_pred",
                      "african_pred",
                      "social_pred",
                      "freeTV_pred",
                      "news_pred")

test_join <- left_join(factor_obs_all, second, by = c("year", "category"))

# function for plotting fitted models
plot_fitted_factors <- function(data, factor) { # factor:one of: popPrint afrikaans soccer  african  social   freeTV    news
        
        if(factor == "social") {
                a <- "social"
                b <- "social_pred"
                # c <- "up_f1"
                # d <- "low_f1"
                e <- "Social"
                f <- "Social with Fitted Values"
        }
        if(factor == "freeTV") {
                a <- "freeTV"
                b <- "freeTV_pred"
                # c <- "up_f2"
                # d <- "low_f2"
                e <- "Free TV"
                f <- "Free TV with Fitted Values"
        }
        if(factor == "afrikaans") {
                a <- "afrikaans"
                b <- "afrikaans_pred"
                # c <- "up_f3"
                # d <- "low_f3"
                e <- "Afrikaans"
                f <- "Afrikaans with Fitted Values"
        }
        if(factor == "soccer") {
                a <- "soccer"
                b <- "soccer_pred"
                # c <- "up_f4"
                # d <- "low_f4"
                e <- "Soccer"
                f <- "Soccer with Fitted Values"
        }
        if(factor == "news") {
                a <- "news"
                b <- "news_pred"
                # c <- "up_f5"
                # d <- "low_f5"
                e <- "News"
                f <- "News with Fitted Values"
        }
        if(factor == "african") {
                a <- "african"
                b <- "african_pred"
                # c <- "up_f6"
                # d <- "low_f6"
                e <- "African"
                f <- "African with Fitted Values"
        }
        if(factor == "popPrint") {
                a <- "popPrint"
                b <- "popPrint_pred"
                # c <- "up_f7"
                # d <- "low_f7"
                e <- "popPrint"
                f <- "PopPrint with Fitted Values"
        }
        
        #plot
        ggplot(data, aes_string("year", a, group = "category")) +
                geom_point(color = "blue", size = 1, fill = "white", alpha = 0.5) +
                geom_line(size = 0.2) +
                geom_line(aes_string("year", b, group = "category"), colour = "red", size = 0.3, linetype = 2 ) +
                facet_grid(.~ category) + theme(axis.text.x = element_text(size = 6)) +
                # geom_errorbar(aes_string(ymax = c, ymin = d), size = 0.3, width = 0.4, alpha = 0.5) +
                labs(y = e, title = f)
        
}
library(gridExtra)
## social
pf_social_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                       factor = "social")
pf_social_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                         factor = "social")
jpeg("social_fitted.jpeg", quality = 100)
grid.arrange(pf_social_up, pf_social_down, nrow = 2)
dev.off()

## freeTV
pf_freeTV_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                             factor = "freeTV")
pf_freeTV_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                               factor = "freeTV")
jpeg("freeTV_fitted.jpeg", quality = 100)
grid.arrange(pf_freeTV_up, pf_freeTV_down, nrow = 2)
dev.off()

## afrikaans
pf_afrikaans_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                                factor = "afrikaans")
pf_afrikaans_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                                  factor = "afrikaans")
jpeg("afrikaans_fitted.jpeg", quality = 100)
grid.arrange(pf_afrikaans_up, pf_afrikaans_down, nrow = 2)
dev.off()

## soccer
pf_soccer_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                                factor = "soccer")
pf_soccer_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                                  factor = "soccer")
jpeg("soccer_fitted.jpeg", quality = 100)
grid.arrange(pf_soccer_up, pf_soccer_down, nrow = 2)
dev.off()

## african
pf_african_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                                factor = "african")
pf_african_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                                  factor = "african")
jpeg("african_fitted.jpeg", quality = 100)
grid.arrange(pf_african_up, pf_african_down, nrow = 2)
dev.off()

## popPrint
pf_popPrint_up <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row1),],
                                                factor = "popPrint")
pf_popPrint_down <- plot_fitted_factors(data = test_join[which(test_join$category %in% vector_row2),],
                                                  factor = "popPrint")
jpeg("popPrint_fitted.jpeg", quality = 100)
grid.arrange(pf_popPrint_up, pf_popPrint_down, nrow = 2)
dev.off()

#
b <- fitted(fit_sem_all)
c <- BIC(fit_sem_all)
show(fit_sem_all)
mi <- modificationIndices(fit_sem_all)
subset(mi, mi > 10)
# not too sure this is best.... 'year' kind of gets lost in the other regression coeff...?? Difficult to see change over time
#  maybe more dummies?.... and also consider bootstrapping options... for defensible std errs...,
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

# # worry about distributions (especially binary, so se not accurate:)
# fit2_boot <- bootstrapLavaan(object = fit2)



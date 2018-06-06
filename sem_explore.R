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

# for sem, want to change variable names since package doesn't take underscores
name_change <- function(set) {
        names(set)[which(names(set)  %in% c('age_actual','hh_inc','lsm_full'))] <- c("age.actual", "hh.inc", "lsm.full")
        return(set)
}

set02_min <- name_change(set02_min)
set08_min <- name_change(set08_min)
set10_min <- name_change(set10_min)
set12_min <- name_change(set12_min)
set14_min <- name_change(set14_min)



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
#  -> exclude 'age' (ie in brackets); and 'lsm' to go with lsm_full
set_min <- rbind.data.frame(add_year(set02_min[,!names(set02_min) %in% c("lifestages", "lifestyle", "attitudes")], 2002),
                            add_year(set08_min[,!names(set08_min) %in% c("lifestages", "lifestyle", "attitudes")], 2008),
                            add_year(set10_min[,!names(set10_min) %in% c("lifestages", "lifestyle", "attitudes")], 2010),
                            add_year(set12_min[,!names(set12_min) %in% c("lifestages", "lifestyle", "attitudes")], 2012),
                            add_year(set14_min[,!names(set14_min) %in% c("lifestages", "lifestyle", "attitudes")], 2014))

# indicate start and end of media vehicles:
strt <- which(names(set_min) == "Business.Day")
lst <- ncol(set_min)

# scaling the media vehicle variables to mean = 0 and sd = 1 for pooled dataset
set_min <- cbind.data.frame(set_min[,1:strt-1], scale(set_min[,strt:lst]))

# save it
saveRDS(set_min, "set_min.rds")
set_min <- readRDS("set_min.rds")

# pca method (package "nFactors") to estimate optimum number of factors to extract:
## Determine Number of Factors to Extract
ev <- eigen(cor(set_min[,strt:lst]))
ap <- parallel(subject=nrow(set_min[,strt:lst]),var=ncol(set_min[,strt:lst]),
               rep=100,cent=.02)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
jpeg("nScree_min_all.jpeg")
plotnScree(nS, main = "Pooled Sets (N = 126 726)") # optimal = 7
dev.off()

library(mvnormtest)
mvnormtest::mshapiro.test(as.matrix(set_min[,strt:lst]))

# to compare ML and PA methods of factor extraction on the full set...
# extract factors (using pa method and common factor extraction)
factors_set_min_pa <- fa(r = set_min[,22:ncol(set_min)], nfactors = 7, fm = "pa")
factors_set_min_ml <- fa(r = set_min[,22:ncol(set_min)], nfactors = 7, fm = "ml")

# consider loadings:
write.csv(round(factors_set_min_pa$loadings, 3), "loadings_pa.csv")
write.csv(round(factors_set_min_ml$loadings, 3), "loadings_ml.csv")

# function to prep sets into numeric for ordered variables and dummies otherwise
prep_sets_sem <- function(set, year = "numeric") { # year = 'numeric' or 'categorical'
        
        # ordered to numeric (only lsm, since year already)
        set$lsm.full <- as.numeric(set$lsm.full)
        
        # To dummies
        require(dummies)
        if(year == 'numeric') {
                set2 <- dummy.data.frame(data = set, names = c("sex", "race", "edu", "hh.inc"))
        }
        if(year == 'categorical') {
                set2 <- dummy.data.frame(data = set, names = c("sex", "race", "edu", "hh.inc", "year"))
        }
        
        # change names for all except year, which is Ok either way
        names(set2)[which(names(set2) == 'sex1')] <- "male"
        names(set2)[which(names(set2) == 'sex2')]  <- "female"
        names(set2)[which(names(set2) == 'race1')]  <- "black"
        names(set2)[which(names(set2) == 'race2')]  <- "coloured"
        names(set2)[which(names(set2) == 'race3')]  <- "indian"
        names(set2)[which(names(set2) == 'race4')]  <- "white"
        names(set2)[which(names(set2) == 'edu1')]  <- "minus.matric"
        names(set2)[which(names(set2) == 'edu2')]  <- "matric"
        names(set2)[which(names(set2) == 'edu3')]  <- "matric.plus"
        names(set2)[which(names(set2) == 'hh.inc1')]  <- "less.5k"
        names(set2)[which(names(set2) == 'hh.inc2')]  <- "x5k.11k"
        names(set2)[which(names(set2) == 'hh.inc3')]  <- "x11k.20k"
        names(set2)[which(names(set2) == 'hh.inc4')]  <- "more.20k"
        return(set2)
        
}

# Options for longitudinal stuff:
# 1:  Use these SEM parameters in some way. Eg.,  model the parameter values over time...?
# 2:  Try working with single value per category level per year again (pseudo panels) and use longitudinal SEM ctsem...?
# 3:  Try suggestion from Francesca: use full set and model "year" as a fixed effect in the structural model (
# year = 'categorical') or (year = 'numeric)??? for now try both

# Option 3a: model "year" as a categorical fixed effect with interaction effects between year categories and all others
# Option 3b: model "year" as numeric fixed effect with interaction effects between year and all others

# set up data
set_prepped_yrCat <- prep_sets_sem(set_min, year = "categorical")
set_prepped_yrNum <- prep_sets_sem(set_min, year = "numeric")

# for numeric, changing "year" to start from 0
set_prepped_yrNum$year <- set_prepped_yrNum$year - 2002

#NBNBNBNBNBNBNBN takes time to run....check to load saved versions ...
# # fit sem: (NBNBNBNBN need to load models into workspace: model_year_categorical.R and model_year_numeric.R)
fit_sem_yrCat <- lavaan::sem(model_year_categorical, data = set_prepped_yrCat, fit.measures = TRUE, meanstructure = TRUE)
fit_sem_yrNum <- lavaan::sem(model_year_numeric, data = set_prepped_yrNum, fit.measures = TRUE, meanstructure = TRUE)




# testing boot (great, but very slow... here just 10 bootstrapped samples)
fit_sem_bootTest <- lavaan::sem(model_test, data = set_prepped_yrNum, fitMeasures = TRUE)#, se = "standard", test = "bootstrap", bootstrap = 5)

summary(fit_sem_bootTest)

boot_cat <- bootstrapLavaan(fit_sem_bootTest, R = 5)
boot_cat # gives free parameter estimates for each of 10 bootstraps
summary(boot_cat) # gives summary stats for each free parameter. So can get for example se of the mean:
str(boot_cat) # named matrix
boot_cat[1:10,1:3] # gets the 10 bootstrapped estimates for first three parameters
str(summary(boot_cat)) # table rows = 6 stats, cols = parameters
summary(boot_cat)[,1]

# function to calculate std errors of the mean:
sems <- function(x) mean(x) + c(-1,1)*2*sqrt(var(x)/length(x))
# calculating confidence intervals (95%) for the bootstrapped estimates of the free parameters
boot_cat_sems <- apply(boot_cat, 2, sems)




# save these to lavaan model objects for 
saveRDS(fit_sem_yrCat, "fit_sem_yrCat.rds")
saveRDS(fit_sem_yrNum, "fit_sem_yrNum.rds")
fit_sem_yrCat <- readRDS("fit_sem_yrCat.rds")
fit_sem_yrNum <- readRDS("fit_sem_yrNum.rds")

# considering the two models:

# # try bootstrapping:
# boot_cat <- bootstrapLavaan(fit_sem_yrCat, R = 10, type = "bollen.stine", 
#                             FUN = "coef")


# ??lavaan::summary # NB: worth checking the help pages on this!!

sum_sem_yrCat <- summary(fit_sem_yrCat, standardized = TRUE, fit.measures = TRUE)
sum_sem_yrNum <- summary(fit_sem_yrNum, standardized = TRUE, fit.measures = TRUE)


write.csv(parameterEstimates(fit_sem_yrCat), file = "paraEst_fit_sem_yrCat.csv")
write.csv(parameterEstimates(fit_sem_yrNum), file = "paraEst_fit_sem_yrNum.csv")

fitMeasures(fit_sem_yrCat)
fitMeasures(fit_sem_yrCat, c("cfi", "tli", "aic", "bic", "rmsea", "srmr")) # would suggest better fit here....,
fitMeasures(fit_sem_yrNum, c("cfi", "tli", "aic", "bic", "rmsea", "srmr")) #  but both actually pretty good.

residuals_yrCat <- residuals(fit_sem_yrCat)
residuals_yrCat$cov # not sure what to do about this?

# Want to consider prediction and plots against aggregated values by demographic category
#  including error bars for confidence intervals of predictions

# will use predictions from year-as-category (they appear better for predictions, although would want to check coeficients of year-as-number later also)

# creating a "yhat" version for use in prediction and also determine confidence intervals for predicted values
yhat <- lavPredict(fit_sem_yrCat, type = "yhat")
X <- yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")]

# extracting the coefficients
coeff_fit <- coef(fit_sem_yrCat)

# extracting the latent factor scores
scores_fit <- data.frame(lavPredict(fit_sem_yrCat, method = "Bartlett")) # not standardised scores!. So, will need to do same for  the predictions also

# create a data frame for preds and scores for each of the latent factors
# use only complete.cases
set_min_complete <- set_min[complete.cases(set_min),]

# sample size
n <- dim(set_min_complete)[1]

# tvalue
t.val <- qt(0.075, n-2)


# 1st inverse of (X'X)^-1
inv_xtx <- solve(t(X)%*%X)

# run for different rows of standard errors
vec_se <- vector()
for(i in 1:n) {
        vec_se[i] <- t(X[i,])%*%inv_xtx%*%X[i,]
        
}

##### freeTV
# observed values
scores_freeTV <- scores_fit[,'freeTV'] 
# predicted values from regression, using model coefficients
pred_freeTV <-  X %*% coeff_fit[which(names(coeff_fit) == "freeTV~year2008"):which(names(coeff_fit) == "freeTV~year2014:indian")]
# correlation for interest sake
corr_freeTV <- cor(pred_freeTV, scores_freeTV)
# SSE (Sum of Squared Errors of prediction)
sse_freeTV <- sum((scores_freeTV - pred_freeTV)^2)
# MSE  (Mean Squared Errors)
mse_freeTV <- sse_freeTV/(n-2)
# standard error of the mean
se_freeTV <- sqrt(mse_freeTV * vec_se)

# upper and lower confidence values
upr_freeTV <- pred_freeTV + t.val * se_freeTV
lwr_freeTV <- pred_freeTV - t.val * se_freeTV

##### popPrint
# observed values
scores_popPrint <- scores_fit[,'popPrint'] 
# predicted values from regression, using model coefficients
pred_popPrint <-  X %*% coeff_fit[which(names(coeff_fit) == "popPrint~year2008"):which(names(coeff_fit) == "popPrint~year2014:indian")]
# correlation for interest sake
corr_popPrint <- cor(pred_popPrint, scores_popPrint)
# SSE (Sum of Squared Errors of prediction)
sse_popPrint <- sum((scores_popPrint - pred_popPrint)^2)
# MSE  (Mean Squared Errors)
mse_popPrint <- sse_popPrint/(n-2)
# standard error of the mean
se_popPrint <- sqrt(mse_popPrint * vec_se)

# upper and lower confidence values
upr_popPrint <- pred_popPrint + t.val * se_popPrint
lwr_popPrint <- pred_popPrint - t.val * se_popPrint

##### news
# observed values
scores_news <- scores_fit[,'news'] 
# predicted values from regression, using model coefficients
pred_news <-  X %*% coeff_fit[which(names(coeff_fit) == "news~year2008"):which(names(coeff_fit) == "news~year2014:indian")]
# correlation for interest sake
corr_news <- cor(pred_news, scores_news)
# SSE (Sum of Squared Errors of prediction)
sse_news <- sum((scores_news - pred_news)^2)
# MSE  (Mean Squared Errors)
mse_news <- sse_news/(n-2)
# standard error of the mean
se_news <- sqrt(mse_news * vec_se)

# upper and lower confidence values
upr_news <- pred_news + t.val * se_news
lwr_news <- pred_news - t.val * se_news

##### afrikaans
# observed values
scores_afrikaans <- scores_fit[,'afrikaans'] 
# predicted values from regression, using model coefficients
pred_afrikaans <-  X %*% coeff_fit[which(names(coeff_fit) == "afrikaans~year2008"):which(names(coeff_fit) == "afrikaans~year2014:indian")]
# correlation for interest sake
corr_afrikaans <- cor(pred_afrikaans, scores_afrikaans)
# SSE (Sum of Squared Errors of prediction)
sse_afrikaans <- sum((scores_afrikaans - pred_afrikaans)^2)
# MSE  (Mean Squared Errors)
mse_afrikaans <- sse_afrikaans/(n-2)
# standard error of the mean
se_afrikaans <- sqrt(mse_afrikaans * vec_se)

# upper and lower confidence values
upr_afrikaans <- pred_afrikaans + t.val * se_afrikaans
lwr_afrikaans <- pred_afrikaans - t.val * se_afrikaans

##### african
# observed values
scores_african <- scores_fit[,'african'] 
# predicted values from regression, using model coefficients
pred_african <-  X %*% coeff_fit[which(names(coeff_fit) == "african~year2008"):which(names(coeff_fit) == "african~year2014:indian")]
# correlation for interest sake
corr_african <- cor(pred_african, scores_african)
# SSE (Sum of Squared Errors of prediction)
sse_african <- sum((scores_african - pred_african)^2)
# MSE  (Mean Squared Errors)
mse_african <- sse_african/(n-2)
# standard error of the mean
se_african <- sqrt(mse_african * vec_se)

# upper and lower confidence values
upr_african <- pred_african + t.val * se_african
lwr_african <- pred_african - t.val * se_african

##### soccer
# observed values
scores_soccer <- scores_fit[,'soccer'] 
# predicted values from regression, using model coefficients
pred_soccer <-  X %*% coeff_fit[which(names(coeff_fit) == "soccer~year2008"):which(names(coeff_fit) == "soccer~year2014:indian")]
# correlation for interest sake
corr_soccer <- cor(pred_soccer, scores_soccer)
# SSE (Sum of Squared Errors of prediction)
sse_soccer <- sum((scores_soccer - pred_soccer)^2)
# MSE  (Mean Squared Errors)
mse_soccer <- sse_soccer/(n-2)
# standard error of the mean
se_soccer <- sqrt(mse_soccer * vec_se)

# upper and lower confidence values
upr_soccer <- pred_soccer + t.val * se_soccer
lwr_soccer <- pred_soccer - t.val * se_soccer

##### social
# observed values
scores_social <- scores_fit[,'social'] 
# predicted values from regression, using model coefficients
pred_social <-  X %*% coeff_fit[which(names(coeff_fit) == "social~year2008"):which(names(coeff_fit) == "social~year2014:indian")]
# correlation for interest sake
corr_social <- cor(pred_social, scores_social)
# SSE (Sum of Squared Errors of prediction)
sse_social <- sum((scores_social - pred_social)^2)
# MSE  (Mean Squared Errors)
mse_social <- sse_social/(n-2)
# standard error of the mean
se_social <- sqrt(mse_social * vec_se)

# upper and lower confidence values
upr_social <- pred_social + t.val * se_social
lwr_social <- pred_social - t.val * se_social


# create single dataset:
set_fin <- set_min_complete %>%
        mutate(pred_popPrint = as.vector(pred_popPrint)) %>%
        mutate(pred_social = as.vector(pred_social)) %>%
        mutate(pred_afrikaans = as.vector(pred_afrikaans)) %>%
        mutate(pred_african = as.vector(pred_african)) %>%
        mutate(pred_news = as.vector(pred_news)) %>%
        mutate(pred_freeTV = as.vector(pred_freeTV)) %>%
        mutate(pred_soccer = as.vector(pred_soccer)) %>%
        mutate(lwr_popPrint = as.vector(lwr_popPrint)) %>%
        mutate(upr_popPrint = as.vector(upr_popPrint)) %>%
        mutate(lwr_social = as.vector(lwr_social)) %>%
        mutate(upr_social = as.vector(upr_social)) %>%
        mutate(lwr_afrikaans = as.vector(lwr_afrikaans)) %>%
        mutate(upr_afrikaans = as.vector(upr_afrikaans)) %>%
        mutate(lwr_african = as.vector(lwr_african)) %>%
        mutate(upr_african = as.vector(upr_african)) %>%
        mutate(lwr_news = as.vector(lwr_news)) %>%
        mutate(upr_news = as.vector(upr_news)) %>%
        mutate(lwr_freeTV = as.vector(lwr_freeTV)) %>%
        mutate(upr_freeTV = as.vector(upr_freeTV)) %>%
        mutate(lwr_soccer = as.vector(lwr_soccer)) %>%
        mutate(upr_soccer = as.vector(upr_soccer))


# function to create frames
frames_factors_all <- function(set, category) {
        require(dplyr)
        
        # set$cluster <- factor(set$cluster, labels = c("cluster1", "cluster2", "cluster3", "cluster4"))
        set$age <- factor(set$age, labels = c("15-24","25-44", "45-54","55+"), ordered = TRUE)
        set$race <- factor(set$race,labels = c("black", "coloured", "indian", "white"), ordered = TRUE)
        set$edu <- factor(set$edu, labels = c("<matric", "matric",">matric" ) ,ordered = TRUE)
        set$lsm <- factor(set$lsm, labels = c("LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10"), ordered = TRUE) #"LSM1-2", 
        set$sex <- factor(set$sex, labels = c("male", "female"), ordered = TRUE)
        set$hh.inc <- factor(set$hh.inc, labels = c("<R2500","R2500-R6999","R7000-R11999",">=R12000"), ordered = TRUE) # NB 2012 levels
        
        
        set %>%
                group_by_(year = "year", category = category) %>%
                summarise(pred_popPrint = mean(pred_popPrint),
                          pred_social = mean(pred_social),
                          pred_soccer = mean(pred_soccer),
                          pred_news = mean(pred_news),
                          pred_freeTV = mean(pred_freeTV),
                          pred_african = mean(pred_african),
                          pred_afrikaans = mean(pred_afrikaans),
                          lwr_popPrint = mean(lwr_popPrint),
                          upr_popPrint = mean(upr_popPrint),
                          lwr_social = mean(lwr_social),
                          upr_social = mean(upr_social),
                          lwr_soccer = mean(lwr_soccer),
                          upr_soccer = mean(upr_soccer),
                          lwr_news = mean(lwr_news),
                          upr_news = mean(upr_news),
                          lwr_freeTV = mean(lwr_freeTV),
                          upr_freeTV = mean(upr_freeTV),
                          lwr_african = mean(lwr_african),
                          upr_african = mean(upr_african),
                          lwr_afrikaans = mean(lwr_afrikaans),
                          upr_afrikaans = mean(upr_afrikaans))
                
}

# function to bind the frames by year
frame_bind_factor_all <- function(set) {
        rbind.data.frame(#frames_factors(set,"cluster"),
                frames_factors_all(set,"sex"),
                frames_factors_all(set,"age"),
                frames_factors_all(set,"edu"),
                frames_factors_all(set,"race"),
                frames_factors_all(set, "hh.inc"),
                frames_factors_all(set,"lsm")) %>% 
                dplyr::select(year,category, everything())
        
}
all_preds_ci <- data.frame(frame_bind_factor_all(set_fin))


vector_row1 <- c("male", "female","15-24","25-44", "45-54","55+","black", "coloured", "indian", "white")
vector_row2 <- c("<matric", "matric",">matric", "<R2500","R2500-R6999","R7000-R11999",">=R12000", "LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10")

# now focussing on the predicted (fitted values):

# function for plotting fitted models
plot_fitted_factors <- function(data, factor) { # factor:one of: popPrint afrikaans soccer  african  social   freeTV    news
        
        if(factor == "social") {
                a <- "pred_social"
                # b <- "pred_social"
                c <- "upr_social"
                d <- "lwr_social"
                e <- "Social"
                f <- "Social: Fitted Values"
        }
        if(factor == "freeTV") {
                a <- "pred_freeTV"
                # b <- "pred_freeTV"
                c <- "upr_freeTV"
                d <- "lwr_freeTV"
                e <- "Free TV"
                f <- "Free TV: Fitted Values"
        }
        if(factor == "afrikaans") {
                a <- "pred_afrikaans"
                # b <- "pred_afrikaans"
                c <- "upr_afrikaans"
                d <- "lwr_afrikaans"
                e <- "Afrikaans"
                f <- "Afrikaans: Fitted Values"
        }
        if(factor == "soccer") {
                a <- "pred_soccer"
                # b <- "pred_soccer"
                c <- "upr_soccer"
                d <- "lwr_soccer"
                e <- "Soccer"
                f <- "Soccer: Fitted Values"
        }
        if(factor == "news") {
                a <- "pred_news"
                # b <- "pred_news"
                c <- "upr_news"
                d <- "lwr_news"
                e <- "News"
                f <- "News: Fitted Values"
        }
        if(factor == "african") {
                a <- "pred_african"
                # b <- "pred_african"
                c <- "upr_african"
                d <- "lwr_african"
                e <- "African"
                f <- "African: Fitted Values"
        }
        if(factor == "popPrint") {
                a <- "pred_popPrint"
                # b <- "pred_popPrint"
                c <- "upr_popPrint"
                d <- "lwr_popPrint"
                e <- "popPrint"
                f <- "PopPrint: Fitted Values"
        }
        
        #plot
        ggplot(data, aes_string("year", a, group = "category")) +
                geom_point(color = "blue", size = 1, fill = "white", alpha = 0.5) +
                geom_line(size = 0.2) +
                # geom_line(aes_string("year", b, group = "category"), colour = "red", size = 0.3, linetype = 2 ) +
                facet_grid(.~ category) +
                theme(axis.text.x = element_text(size = 6)) +
                geom_errorbar(aes_string(ymax = c, ymin = d), size = 0.3, width = 0.4, alpha = 0.5) +
                labs(y = e, title = f)
        
}


library(gridExtra)
## social
pf_social_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                       factor = "social")
pf_social_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                         factor = "social")
jpeg("social_fitted.jpeg", quality = 100)
grid.arrange(pf_social_up, pf_social_down, nrow = 2)
dev.off()

## freeTV
pf_freeTV_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                             factor = "freeTV")
pf_freeTV_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                               factor = "freeTV")
jpeg("freeTV_fitted.jpeg", quality = 100)
grid.arrange(pf_freeTV_up, pf_freeTV_down, nrow = 2)
dev.off()

## afrikaans
pf_afrikaans_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                                factor = "afrikaans")
pf_afrikaans_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                                  factor = "afrikaans")
jpeg("afrikaans_fitted.jpeg", quality = 100)
grid.arrange(pf_afrikaans_up, pf_afrikaans_down, nrow = 2)
dev.off()

## soccer
pf_soccer_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                                factor = "soccer")
pf_soccer_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                                  factor = "soccer")
jpeg("soccer_fitted.jpeg", quality = 100)
grid.arrange(pf_soccer_up, pf_soccer_down, nrow = 2)
dev.off()

## african
pf_african_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                                factor = "african")
pf_african_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                                  factor = "african")
jpeg("african_fitted.jpeg", quality = 100)
grid.arrange(pf_african_up, pf_african_down, nrow = 2)
dev.off()

## popPrint
pf_popPrint_up <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row1),],
                                                factor = "popPrint")
pf_popPrint_down <- plot_fitted_factors(data = all_preds_ci[which(all_preds_ci$category %in% vector_row2),],
                                                  factor = "popPrint")
jpeg("popPrint_fitted.jpeg", quality = 100)
grid.arrange(pf_popPrint_up, pf_popPrint_down, nrow = 2)
dev.off()



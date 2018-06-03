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
# save it
saveRDS(set_min, "set_min.rds")

# indicate start and end of media vehicles:
strt <- which(names(set_min) == "Business.Day")
lst <- ncol(set_min)

# scaling the media vehicle variables to mean = 0 and sd = 1 for pooled dataset
set_min <- cbind.data.frame(set_min[,1:strt-1], scale(set_min[,strt:lst]))

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

# save these to lavaan model objects for 
saveRDS(fit_sem_yrCat, "fit_sem_yrCat.rds")
saveRDS(fit_sem_yrNum, "fit_sem_yrNum.rds")
fit_sem_yrCat <- readRDS("fit_sem_yrCat.rds")
fit_sem_yrNum <- readRDS("fit_sem_yrNum.rds")

# considering the two models:

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
#  in a way that allows for comparisons with pseudo panel approach:

# will use predictions from year-as-category (they appear better for predictions, although would want to check coeficients of year-as-number later also)
# may want to get back to lavPredict(object, type = "yhat")...??
# creating a "yhat" version for use in prediction:

yhat <- lavPredict(fit_sem_yrCat, type = "yhat")
# extracting the coefficients
coeff_fit <- coef(fit_sem_yrCat)
# extracting the latent factor scores
scores_fit <- data.frame(lavPredict(fit_sem_yrCat, method = "Bartlett")) # not standardised scores!. So, will need to do same for  the predictions also

# create a data frame for preds and scores for each of the latent factors
# use only complete.cases
set_min_complete <- set_min[complete.cases(set_min),]

# vectors for popPrint
pred_popPrint <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
        coeff_fit[which(names(coeff_fit) == "popPrint~year2008"):which(names(coeff_fit) == "popPrint~year2014:indian")]
scores_popPrint <- scores_fit[,'popPrint'] 
corr_popPrint <- cor(pred_popPrint, scores_popPrint) # for interest sake
mse_popPrint <- mean((pred_popPrint - scores_popPrint)^2) # for interest sake

# vectors for news
pred_news <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                                coeff_fit[which(names(coeff_fit) == "news~year2008"):which(names(coeff_fit) == "news~year2014:indian")]
scores_news <- scores_fit[,'news'] 
corr_news <- cor(pred_news, scores_news) # for interest sake
mse_news <- mean((pred_news - scores_news)^2) # for interest sake

# vectors for social
pred_social <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                                coeff_fit[which(names(coeff_fit) == "social~year2008"):which(names(coeff_fit) == "social~year2014:indian")]
scores_social <- scores_fit[,'social'] 
corr_social <- cor(pred_social, scores_social) # for interest sake
mse_social <- mean((pred_social - scores_social)^2) # for interest sake

# vectors for afrikaans
pred_afrikaans <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                                coeff_fit[which(names(coeff_fit) == "afrikaans~year2008"):which(names(coeff_fit) == "afrikaans~year2014:indian")]
scores_afrikaans <- scores_fit[,'afrikaans'] 
corr_afrikaans <- cor(pred_afrikaans, scores_afrikaans) # for interest sake
mse_afrikaans <- mean((pred_afrikaans - scores_afrikaans)^2) # for interest sake

# vectors for african
pred_african <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                                coeff_fit[which(names(coeff_fit) == "african~year2008"):which(names(coeff_fit) == "african~year2014:indian")]
scores_african <- scores_fit[,'african'] 
corr_african <- cor(pred_african, scores_african) # for interest sake
mse_african <- mean((pred_african - scores_african)^2) # for interest sake

# vectors for soccer
pred_soccer <-  yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                               coeff_fit[which(names(coeff_fit) == "soccer~year2008"):which(names(coeff_fit) == "soccer~year2014:indian")]
scores_soccer <- scores_fit[,'soccer']
corr_soccer <- cor(pred_soccer, scores_soccer) # for interest sake
mse_soccer <- mean((pred_soccer - scores_soccer)^2) # for interest sake

# vectors for freeTV
pred_freeTV <- yhat[,which(names(data.frame(yhat)) == "year2008"):which(names(data.frame(yhat)) == "year2014.indian")] %*%
                              coeff_fit[which(names(coeff_fit) == "freeTV~year2008"):which(names(coeff_fit) == "freeTV~year2014:indian")]
scores_freeTV <- scores_fit[,'freeTV']
corr_freeTV <- cor(pred_freeTV, scores_freeTV) # for interest sake
mse_freeTV <- mean((pred_freeTV - scores_freeTV)^2) # for interest sake

# create single dataset:
set_fin <- set_min_complete %>%
        mutate(pred_popPrint = as.vector(pred_popPrint)) %>%
        mutate(scores_popPrint = as.vector(scores_popPrint)) %>%
        mutate(pred_social = as.vector(pred_social)) %>%
        mutate(scores_social = as.vector(scores_social)) %>%
        mutate(pred_afrikaans = as.vector(pred_afrikaans)) %>%
        mutate(scores_afrikaans = as.vector(scores_afrikaans)) %>%
        mutate(pred_african = as.vector(pred_african)) %>%
        mutate(scores_african = as.vector(scores_african)) %>%
        mutate(pred_news = as.vector(pred_news)) %>%
        mutate(scores_news = as.vector(scores_news)) %>%
        mutate(pred_freeTV = as.vector(pred_freeTV)) %>%
        mutate(scores_freeTV = as.vector(scores_freeTV)) %>%
        mutate(pred_soccer = as.vector(pred_soccer)) %>%
        mutate(scores_soccer = as.vector(scores_soccer))

# # library(psych)
# # lavaan.diagram(fit_sem_all, cut = 0.3, cex = 0.5, regression = TRUE)
# 
# 
# fit_sem_parameters <- as.data.frame(parameterEstimates(fit_sem_all, rsquare = TRUE, standardized = TRUE))
# fit_sem_parameters[,4:ncol(fit_sem_parameters)] <- round(fit_sem_parameters[,4:ncol(fit_sem_parameters)], 4)
# ##  try to do plots comparing fitted with actual by category and by year:
# 
# # define the model matrix (NB order of columns: freeTV, social, news, afrikaans, popPrint, soccer, african)
# mod_mat <- matrix(c(0,0,0,0,1,0,0, 
#                     0,0,0,0,1,0,0,
#                     0,0,0,1,0,0,0,
#                     0,0,0,0,1,0,0,
#                     0,0,0,0,1,0,0,
#                     0,0,0,0,0,1,0,
#                     0,0,0,0,0,0,1,
#                     0,0,0,1,0,0,0,
#                     0,0,0,0,1,0,0,
#                     0,0,0,0,0,1,0,
#                     0,0,0,0,0,0,1,
#                     0,0,0,0,1,0,0,
#                     0,0,0,0,1,0,0,
#                     0,0,0,0,1,0,0,
#                     0,0,0,1,0,0,0,
#                     0,0,0,0,1,0,0,
#                     0,1,0,0,0,0,0,
#                     0,0,0,0,0,0,1,
#                     1,0,0,0,0,0,0,
#                     1,0,0,0,0,0,0,
#                     1,0,0,0,0,0,0,
#                     1,0,0,0,0,0,0,
#                     0,1,0,0,0,0,0,
#                     0,0,1,0,0,0,0,
#                     0,1,0,0,0,0,0,
#                     0,0,1,0,0,0,0,
#                     0,1,0,0,0,0,0,
#                     0,1,0,0,0,0,0), nrow = 28, ncol = 7, byrow = TRUE)
# 
# # construct dataframe of standardised observed by the model:
# obs <-as.matrix(set_min[,22:ncol(set_min)]) %*% mod_mat
# colnames(obs) <- c("freeTV", "social", "news", "afrikaans", "popPrint","soccer", "african")
# comb <- cbind.data.frame(set_min[1:21],scale(obs))

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
                          scores_popPrint = mean(scores_popPrint),
                          pred_social = mean(pred_social),
                          scores_social = mean(scores_social),
                          pred_soccer = mean(pred_soccer),
                          scores_soccer = mean(scores_soccer),
                          pred_news = mean(pred_news),
                          scores_news = mean(scores_news),
                          pred_freeTV = mean(pred_freeTV),
                          scores_freeTV = mean(scores_freeTV),
                          pred_african = mean(pred_african),
                          scores_african = mean(scores_african),
                          pred_afrikaans = mean(pred_afrikaans),
                          scores_afrikaans = mean(scores_afrikaans))
                          
        
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
all_preds_scores <- data.frame(frame_bind_factor_all(set_fin))

# all_plots_factors <- function(data, title = "All Scores Observed") {
#         ggplot(data = data, title = title) +
#                 geom_line(aes(year, social, group = category, colour = "Social")) +
#                 geom_line(aes(year, freeTV, group = category, colour = "Free TV")) +
#                 geom_line(aes(year, afrikaans, group = category, colour = "Afrikaans")) +
#                 geom_line(aes(year, soccer, group = category, colour = "Soccer")) +
#                 geom_line(aes(year, news, group = category, colour = "News")) +
#                 geom_line(aes(year, african, group = category, colour = "African")) +
#                 geom_line(aes(year, popPrint, group = category, colour = "PopPrint")) +
#                 scale_colour_discrete(name="Factors") +
#                 facet_grid(. ~ category) +
#                 theme(axis.text.x = element_text(size = 6)) +
#                 labs(y = "aggregate scores observed", title = title)
#         
# }
# 
# all_plots_factors(factor_scores_all)

vector_row1 <- c("male", "female","15-24","25-44", "45-54","55+","black", "coloured", "indian", "white")
vector_row2 <- c("<matric", "matric",">matric", "<R2500","R2500-R6999","R7000-R11999",">=R12000", "LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10")

# now focussing on the predicted (fitted values):

# function for plotting fitted models
plot_fitted_factors <- function(data, factor) { # factor:one of: popPrint afrikaans soccer  african  social   freeTV    news
        
        if(factor == "social") {
                a <- "scores_social"
                b <- "pred_social"
                # c <- "up_f1"
                # d <- "low_f1"
                e <- "Social"
                f <- "Social with Fitted Values"
        }
        if(factor == "freeTV") {
                a <- "scores_freeTV"
                b <- "pred_freeTV"
                # c <- "up_f2"
                # d <- "low_f2"
                e <- "Free TV"
                f <- "Free TV with Fitted Values"
        }
        if(factor == "afrikaans") {
                a <- "scores_afrikaans"
                b <- "pred_afrikaans"
                # c <- "up_f3"
                # d <- "low_f3"
                e <- "Afrikaans"
                f <- "Afrikaans with Fitted Values"
        }
        if(factor == "soccer") {
                a <- "scores_soccer"
                b <- "pred_soccer"
                # c <- "up_f4"
                # d <- "low_f4"
                e <- "Soccer"
                f <- "Soccer with Fitted Values"
        }
        if(factor == "news") {
                a <- "scores_news"
                b <- "pred_news"
                # c <- "up_f5"
                # d <- "low_f5"
                e <- "News"
                f <- "News with Fitted Values"
        }
        if(factor == "african") {
                a <- "scores_african"
                b <- "pred_african"
                # c <- "up_f6"
                # d <- "low_f6"
                e <- "African"
                f <- "African with Fitted Values"
        }
        if(factor == "popPrint") {
                a <- "scores_popPrint"
                b <- "pred_popPrint"
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
pf_social_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                       factor = "social")
pf_social_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                         factor = "social")
jpeg("social_fitted.jpeg", quality = 100)
grid.arrange(pf_social_up, pf_social_down, nrow = 2)
dev.off()

## freeTV
pf_freeTV_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                             factor = "freeTV")
pf_freeTV_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                               factor = "freeTV")
jpeg("freeTV_fitted.jpeg", quality = 100)
grid.arrange(pf_freeTV_up, pf_freeTV_down, nrow = 2)
dev.off()

## afrikaans
pf_afrikaans_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                                factor = "afrikaans")
pf_afrikaans_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                                  factor = "afrikaans")
jpeg("afrikaans_fitted.jpeg", quality = 100)
grid.arrange(pf_afrikaans_up, pf_afrikaans_down, nrow = 2)
dev.off()

## soccer
pf_soccer_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                                factor = "soccer")
pf_soccer_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                                  factor = "soccer")
jpeg("soccer_fitted.jpeg", quality = 100)
grid.arrange(pf_soccer_up, pf_soccer_down, nrow = 2)
dev.off()

## african
pf_african_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                                factor = "african")
pf_african_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                                  factor = "african")
jpeg("african_fitted.jpeg", quality = 100)
grid.arrange(pf_african_up, pf_african_down, nrow = 2)
dev.off()

## popPrint
pf_popPrint_up <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row1),],
                                                factor = "popPrint")
pf_popPrint_down <- plot_fitted_factors(data = all_preds_scores[which(all_preds_scores$category %in% vector_row2),],
                                                  factor = "popPrint")
jpeg("popPrint_fitted.jpeg", quality = 100)
grid.arrange(pf_popPrint_up, pf_popPrint_down, nrow = 2)
dev.off()




# Try seperate linear regression equations to determine errors
test <- lm(set_fin$scores_freeTV ~ year + age.actual + sex + edu + hh.inc + race + as.numeric(lsm.full) +
                   year*age.actual +
                   year*sex +
                   year*edu +
                   year*hh.inc +
                   year*race +
                   year*as.numeric(lsm.full),
           data = set_min_complete)
coef(test)
preds <- predict(test, se.fit = TRUE)

# se for conf interval
se.CI <- preds$se.fit
# se for pred interval
se.PI <- sqrt(se.CI^2 + preds$residual.scale^2)

#conf/pred intervals at 90% conf
alpha <- qt((1-0.9)/2, df = preds$df)
# [1] -1.681071
CIup <- preds$fit +  alpha* se.CI# c(alpha, -alpha)
CIdown <- preds$fit -  alpha* se.CI# c(alpha, -alpha)
CI <- preds$fit + c(alpha, -alpha)* se.CI 
# [1] 87.28387 91.97880
PI <- preds$fit + c(alpha, -alpha) * se.PI
# [1]  74.46433 104.79833





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



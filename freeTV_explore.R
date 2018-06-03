library(lavaan)
library(tidyverse)

# read in lavaan sem object
fit_sem_yrCat <- readRDS("fit_sem_yrCat.rds")

# get and isolate scores for freeTV (outcome variable)
scores_sem_yrCat <- lavPredict(fit_sem_yrCat)
scores_sem_yrCat_freeTV <- scores_sem_yrCat[,'freeTV']
str(scores_sem_yrCat_freeTV)

# create single datafile for this exercise\
# read in set_min from "sem_explore.R"
set_min <- readRDS("set_min.rds")
frame_sem_yrCat_freeTV <- set_min[complete.cases(set_min),] %>%
        dplyr::select(qn, year, age.actual, sex, edu, hh.inc, race, lsm.full) %>%
        mutate(free_TV = scores_sem_yrCat_freeTV)
frame_sem_yrCat_freeTV$year <- factor(frame_sem_yrCat_freeTV$year, ordered = TRUE)
frame_sem_yrCat_freeTV$lsm.full <- as.numeric(frame_sem_yrCat_freeTV$lsm.full)
# run linear regression model with demographics as predictors
lm_sem_yrCat_freeTV <- lm(free_TV ~  . +
                                  year:age.actual +
                                  year:sex +
                                  year:edu +
                                  year:hh.inc +
                                  year:race +
                                  year:lsm.full,
                          data = frame_sem_yrCat_freeTV[,-1])
summary(lm_sem_yrCat_freeTV)

pred_lm_sem_yrCat_freeTV <- predict(lm_sem_yrCat_freeTV, se.fit = TRUE, interval = "confidence", type = "response")

#for plots:
for_plots <- cbind.data.frame(frame_sem_yrCat_freeTV,pred_lm_sem_yrCat_freeTV$fit)
for_plots$year <- as.numeric(as.character(for_plots$year))
for_plots <- for_plots %>%
        mutate(lsm = set_min$lsm[complete.cases(set_min)]) %>%
        mutate(age = set_min$age[complete.cases(set_min)])

# create means and try visuals::::
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
                summarise(scores = mean(free_TV),
                          fita = mean(fit),
                          lwra = mean(lwr),
                          upra = mean(upr),
                          me_fit = err_mean(fit),
                          me_lwr = err_mean(lwr),
                          me_upr = err_mean(upr)
                          )
        
        
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

all_trial <- data.frame(frame_bind_factor_all(for_plots))

err_mean <- function(x) {
        2*sqrt((var(x)/length(x)))
}

err_mean(freeTV_explore$fit)
all_fit <- all_trial %>%
        group_by(year) %>%
        summarise(scores = mean(scores),
                  fita = mean(fit),
                  lwra = mean(lwr),
                  upra = mean(upr),
                  me_fit = err_mean(fit),
                  me_lwr = err_mean(lwr),
                  me_upr = err_mean(upr)
                  )
        

vector_row1 <- c("male", "female","15-24","25-44", "45-54","55+","black", "coloured", "indian", "white")
vector_row2 <- c("<matric", "matric",">matric", "<R2500","R2500-R6999","R7000-R11999",">=R12000", "LSM1-2", "LSM3-4", "LSM5-6", "LSM7-8", "LSM9-10")

ggplot(all_trial[which(all_trial$category %in% vector_row2),], aes_string("year", "scores", group = "category")) +
        geom_point(color = "blue", size = 1, fill = "white", alpha = 0.5) +
        geom_line(size = 0.2) +
        geom_line(aes_string("year", "fita", group = "category"), colour = "red", size = 0.3, linetype = 2 ) +
        facet_grid(.~ category) + theme(axis.text.x = element_text(size = 6)) +
        geom_errorbar(aes_string(ymax = "upra", ymin = "lwra"), size = 0.3, width = 0.4, alpha = 0.5) +
        labs(y = "engagement", title = "freeTV")

ggplot(all_fit, aes(year, scores)) +
        geom_point(color = "blue", size = 1, fill = "white", alpha = 0.5) +
        geom_line(size = 0.2) +
        geom_line(aes(year, fit), colour = "red", size = 0.3, linetype = 2 ) +
        # facet_grid(.~ category) + theme(axis.text.x = element_text(size = 6)) +
        geom_errorbar(aes(ymax = upr, ymin = lwr), size = 0.3, width = 0.4, alpha = 0.5) +
        labs(y = "all engagement", title = "freeTV")

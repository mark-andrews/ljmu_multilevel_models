library(tidyverse)
library(lme4) # alternatively load lmerTest
library(lmerTest)
library(performance)

# Read data ---------------------------------------------------------------

sleep_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/sleepstudy.csv")

# Multilevel linear model aka linear mixed effects ------------------------

M_3 <- lmer(Reaction ~ Days + (Days|Subject), data = sleep_df )

# random intercepts only
M_4 <- lmer(Reaction ~ 1 + Days + (1 |Subject), data = sleep_df )

# * REML v ML
# * logLik, deviance, AIC, BIC
# * Wilks's theorem and log likelihood ratio tests

anova(M_4, M_3)

# Refit using MLE
M_3_mle <- lmer(Reaction ~ Days + (Days|Subject), 
                REML = FALSE, 
                data = sleep_df )

# random intercepts only
M_4_mle <- lmer(Reaction ~ 1 + Days + (1 |Subject), 
                REML = FALSE,
                data = sleep_df )

# log likelihood
logLik(M_3_mle)
logLik(M_4_mle)

# -2 * logLikeilood: aka deviance 
-2 * logLik(M_3_mle)
-2 * logLik(M_4_mle)

deviance(M_3_mle)
deviance(M_4_mle)

# AIC 
deviance(M_3_mle) + 2 * 6
AIC(M_3_mle)
deviance(M_4_mle) + 2 * 4
AIC(M_4_mle)

# BIC: -2 log like + log(n) * p
deviance(M_3_mle) + log(180) * 6
BIC(M_3_mle)

# Wilks's theorem
# under the null, this difference is a Chi sq distrib
delta_dev <- deviance(M_4_mle) - deviance(M_3_mle)

# what is the prob of a value greater than delta_dev
# in a chisq dist with 2 df?
pchisq(delta_dev, df = 2, lower.tail = F)

anova(M_3_mle, M_4_mle)


classroom_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/classroom.csv")


# * random slopes and random intercepts but NO correlation linear model


# * random slopes and random intercepts and correlation linear model
M_7 <- lmer(mathscore ~ ses + (ses|schoolid), REML = F, data = classroom_df)

# * random intercepts only linear model
M_8 <- lmer(mathscore ~ ses + (1|schoolid), REML = F, data = classroom_df)

# * random slopes only linear model
M_9 <- lmer(mathscore ~ ses + (0 + ses|schoolid), REML = F, data = classroom_df)

# * random slopes and random intercepts but NO correlation linear model
M_10 <- lmer(mathscore ~ ses + (ses||schoolid), REML = F, data = classroom_df)



ses_models <- list(rsri = M_7,
                   ri = M_8,
                   rs = M_9,
                   rsri_no_c = M_10)

# using a purrr map
map_dbl(ses_models, logLik)
map_dbl(ses_models, deviance)
map_dbl(ses_models, AIC)
map_dbl(ses_models, BIC)

anova(M_7, M_10)
anova(M_7, M_8)
confint(M_7)

anova(M_7, M_8, M_9, M_10)

# Science -----------------------------------------------------------------

science_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/science.csv")

M_11 <- lmer(like ~ sex + (sex|school), data = science_df)
summary(M_11)

M_12 <- lmer(like ~ sex + PrivPub + (sex|school), data = science_df)
summary(M_12)

M_13 <- lmer(like ~ sex + PrivPub + (sex|school) + (sex|State), data = science_df)
count(science_df, State)

M_13 <- lmer(like ~ sex + PrivPub + (sex |school) + State, data = science_df)
M_14 <- lmer(like ~ sex * State + PrivPub + (sex|school) , data = science_df)
summary(M_13)

M_15 <- lmer(like ~ sex + PrivPub + (1|school) + State, data = science_df)


# (sex + PrivPub | school)

# Anova -------------------------------------------------------------------

df_long <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/df_long.csv")
library(afex)
M_16 <- aov_car(score ~ treat + Error(id/time), data = df_long)
summary(M_16)

M_17 <- lmer(score ~ treat * time + (1|id), data = df_long)
anova(M_17)

# Crossed -----------------------------------------------------------------

blp_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/blp-short2.csv")

library(skimr)
skim(blp_df)


M_18 <- lmer(rt ~ freq + (freq|participant), data = blp_df)

blp_df <- mutate(blp_df, freq2 = scale(freq)[,])

M_19 <- lmer(rt ~ freq2 + (freq2|participant), data = blp_df)

M_20 <- lmer(rt ~ freq2 + (freq2|participant) + (1|spelling), 
             data = blp_df)

# does not work
M_21 <- lmer(rt ~ freq2 + (freq2|participant) + (freq2|spelling), 
             data = blp_df)

anova(M_19, M_20)

summary(M_20)


blp_df3 <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/blp-short3.csv")

M_22 <- glmer(accuracy ~ lex + (1|participant) + (1|spell),
              data = blp_df3,
              family = binomial())

summary(M_22)



# Power analysis using simr -----------------------------------------------

library(simr)

# using sleep study
VarCorr(M_3)
fixef(M_3)

fake_sleep_data <- expand_grid(subject = 1:10, days = 0:9)

fixed_b <- c(250, 5)

s <- 25

V <- matrix(c(25^2, 0, 0, 5^2), nrow = 2)

fake_model <- makeLmer(y ~ days + (days|subject),
                       data = fake_sleep_data,
                       fixef = fixed_b,
                       VarCorr = V,
                       sigma = s)

summary(fake_model)

powerSim(fake_model, nsim = 100)
powerSim(fake_model, test = fixed('days'), nsim = 1000)

fake_sleep_data2 <- expand_grid(subject = 1:10, days = 0:9, 
                                x = c('A', 'B'))


fixed_b <- c(250, 10, 5, 1)

s <- 25

V <- matrix(c(25^2, 0, 0, 5^2), nrow = 2)

fake_model2 <- makeLmer(y ~ x * days + (days|subject),
                        data = fake_sleep_data2,
                        fixef = fixed_b,
                        VarCorr = V,
                        sigma = s)

summary(fake_model2)
powerSim(fake_model2, test = fixed('xB:days'), nsim = 1000)
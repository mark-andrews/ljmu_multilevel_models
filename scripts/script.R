# Load packages -----------------------------------------------------------

library(tidyverse)
library(modelr)
library(lmerTest) # alternatively, use `library(lme4)`

# load some utils
source('scripts/utils.R')

# Read in sleep study data ------------------------------------------------

sleep_df <- read_csv("data/sleepstudy.csv") 


# Plot the data -----------------------------------------------------------

sleep_df %>% 
  ggplot(mapping = aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none')


# A pooled model ----------------------------------------------------------

sleep_df %>% 
  ggplot(mapping = aes(x = Days, y = Reaction)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)  




M_1 <- lm(Reaction ~ Days, data = sleep_df)

# * Problems
#   - poor fit: see residual sd, see R^2, see looic_lm
#   - pseudoreplication: we do not have 180 independent data points
#   - over-influence of one individual possibly

# poor fit:
add_predictions(sleep_df, M_1) %>% 
  ggplot(mapping = aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  geom_line(aes(y = pred), colour = 'black') + 
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none')

# poor fit: residual standard deviation
sigma(M_1)
# also, 
r <- M_1$residuals
sqrt(sum((r-mean(r))^2)/M_1$df.residual)

# poor fit: R^2
summary(M_1)$r.sq
# also: 1 - RSS/TSS
1 - var(r)/var(sleep_df$Reaction)
# also: ESS/TSS
var(predict(M_1))/var(sleep_df$Reaction)

# R^2: Def 3: var of predicted values / var of residual + var of predicted values
var(predict(M_1)) / (var(predict(M_1)) + var(residuals(M_1)))

# R^2: Def 4: 
M_0 <- lm(Reaction ~ 1, data = sleep_df)
1 - var(residuals(M_1))/var(residuals(M_0))

# poor fit: looic: -2 * elpd
looic_lm(M_1)


# No pooling model --------------------------------------------------------

# varying slope and varying intercept model
M_2 <- lm(Reaction ~ Subject*Days, data = sleep_df)

add_predictions(sleep_df, M_2) %>% 
  ggplot(mapping = aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  geom_line(aes(y = pred), colour = 'black') + 
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none')

# fit?
sigma(M_2)
summary(M_2)$r.sq
looic_lm(M_2)

# Problems:
#   - no generalization
#   - no partial pooling: 10 data points per subject only


# Partial pooling model, aka multilevel model -----------------------------

M_3 <- lmerTest::lmer(Reaction ~ Days + (Days|Subject), data = sleep_df)

summary(M_3)
fixef(M_3)
ranef(M_3)
coef(M_3)

# partial pooling
add_predictions(sleep_df, M_3) %>% 
  ggplot(mapping = aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  geom_line(aes(y = pred), colour = 'black') + 
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none')


# R^2 in mixed effects ----------------------------------------------------

# look at fixed effect predictions AND the fixed + random predictions
predict(M_3) # predictions based on fixed and random effects
predict(M_3, re.form = NA) # predictions based on fixed effects only

add_predictions(sleep_df, M_3) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RI/RS/C')

mutate(sleepstudy, pred = predict(M_3, re.form = NA)) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RI/RS/C: fixed effects only prediction')

# both
add_predictions(sleep_df, M_3) |>
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) + # adds the best fit using fixed effects and random effects
  geom_line(data = mutate(sleep_df, pred = predict(M_3, re.form = NA)), aes(y = pred)) + # add a best fit only using fixed effects
  facet_wrap(~Subject)+
  theme_minimal() +
  theme(legend.position = 'none') 


# Approx 1 of R^2 in mixed effects: Conditional R^2
var(predict(M_3)) / var(sleep_df$Reaction)

# Approx 2 of R^2 in mixed effects: Marginal R^2
var(predict(M_3, re.form = NA)) / var(sleep_df$Reaction)


performance::r2_nakagawa(M_3) 
# A general and simple method for obtaining R2 from generalized linear mixed-effects models
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x
# see also
# Quantifying Explained Variance in Multilevel Models:
#   An Integrative Framework for Defining R-Squared Measures
# Jason D. Rights and Sonya K. Sterba
# https://cdn.vanderbilt.edu/vu-my/wp-content/uploads/sites/2472/2019/09/23132255/rights_sterba_2019.pdf


# Bayes -------------------------------------------------------------------

M_4 <- brm(Reaction ~ Days + (Days|Subject), data = sleep_df)
M_5 <- brm(Reaction ~ Days , data = sleep_df)
M_6 <- brm(Reaction ~ Days * Subject, data = sleep_df)

loo(M_4, M_5, M_6)

# Alternative models ------------------------------------------------------

# random intercepts only
M_7 <- lmer(Reaction ~ 1 + Days + (1|Subject), 
            data = sleep_df)

summary(M_7)

# random slopes only
M_8 <- lmer(Reaction ~ Days + (0 + Days|Subject), 
            data = sleep_df)
summary(M_8)


# random slopes AND random intercepts AND no correlation 
M_9 <- lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), 
             data = sleep_df)

# alternatively, random slopes AND random intercepts AND no correlation 
M_10 <- lmer(Reaction ~ Days + (Days||Subject), 
             data = sleep_df)

summary(M_9)


add_predictions(sleep_df, M_3) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RI/RS/C')

add_predictions(sleep_df, M_7) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RI')

add_predictions(sleep_df, M_8) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RS')

add_predictions(sleep_df, M_9) %>% 
  ggplot(aes(x = Days, y = Reaction, colour = Subject)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~Subject) +
  theme_minimal() +
  theme(legend.position = 'none') +
  ggtitle('RI/RS no C')

# Model loglikelihood & deviances & LRT -----------------------------------


M_3_mle <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy, REML = FALSE)

M_10_mle <- lmer(Reaction ~ Days + (Days||Subject), 
                 REML = FALSE,
                 data = sleepstudy)

logLik(M_3_mle)
deviance(M_3_mle)

logLik(M_10_mle)
deviance(M_10_mle)

deviance(M_10_mle) - deviance(M_3_mle)

anova(M_10_mle, M_3_mle)

# anova automatically refits with REML = FALSE
anova(M_10, M_3) # test for correlation
anova(M_7, M_10) # test for random slopes
anova(M_8, M_10) # test for random intercepts



# Nested ------------------------------------------------------------------


classroom_df <- read_csv("data/classroom.csv")


ggplot(classroom_df, aes(x = ses, y = mathscore)) +
  geom_point() +
  stat_smooth(method = 'lm')


M_14 <- lmer(mathscore ~ ses + (ses|schoolid), data = classroom_df)
summary(M_14)
confint(M_14)


M_15 <- lmer(mathscore ~ ses + (ses|schoolid) + (ses|classid), data = classroom_df)


# remove random slopes for classes

M_18 <- lmer(mathscore ~ ses + (ses|schoolid) + (1|classid), 
             data = classroom_df)

M_19 <- lmer(mathscore ~ ses + (ses|schoolid) + (1|schoolid/classid2), 
             data = classroom_df)

M_19a <- lmer(mathscore ~ ses + (ses|schoolid) + (1|schoolid:classid2), 
             data = classroom_df)


M_20 <- lmer(mathscore ~ ses + (1|schoolid) + (1|classid), 
             data = classroom_df)
M_21 <- lmer(mathscore ~ ses + (1|schoolid/classid2), 
             data = classroom_df)
M_22 <- lmer(mathscore ~ ses + (1|schoolid) + (1|schoolid:classid2),
             data = classroom_df)



# Crossed -----------------------------------------------------------------

load('data/blp_trials.Rda')

blp_words <- filter(blp_trials, lex == 'word')

M_11 <- lmer(rt ~ freq + (freq|subject) + (1|item), data = blp_words)

# re scale freq
blp_words <- mutate(blp_words, freq = scale(log(freq))[,1])

M_11 <- lmer(rt ~ freq + (freq|subject) + (1|item), data = blp_words)



# Power analysis ----------------------------------------------------------

library(simr)

fake_sleep_data <- expand_grid(subject = 1:10, days = 0:10)

b <- c(250, 5)

V <- matrix(c(25^2, 0, 0, 5^2), nrow = 2)

s <- 25

fake_model <- makeLmer(y ~ days + (days|subject), data = fake_sleep_data, fixef = b, VarCorr = V, sigma = s)

powerSim(fake_model, nsim = 100)

M_3 <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)

powerSim(M_3, nsim = 100)



fake_sleep_data2 <- mutate(fake_sleep_data, x = rnorm(n = n()))

b2 <- c(250, 5, 1)

V <- matrix(c(25^2, 0, 0, 5^2), nrow = 2)

s <- 25

fake_model2 <- makeLmer(y ~ days + x + (days|subject), data = fake_sleep_data2, fixef = b2, VarCorr = V, sigma = s)

powerSim(fake_model2, test = fixed('x'), nsim = 100)

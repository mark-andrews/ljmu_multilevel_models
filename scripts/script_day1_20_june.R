library(tidyverse)
library(lme4) # alternatively load lmerTest
library(lmerTest)
library(performance)

# Read data ---------------------------------------------------------------

sleep_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/ljmu_multilevel_models/main/data/sleepstudy.csv")

# visualize the data ---------------------------------
ggplot(sleep_df, 
       aes(x = Days, y = Reaction, colour = Subject)
) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE) +
  facet_wrap(~Subject)


# Pooled analysis ---------------------------------------------------------

M_1 <- lm(Reaction ~ 1 + Days, data = sleep_df)

summary(M_1)


# No pooling model --------------------------------------------------------


M_2 <- lm(Reaction ~ Days * Subject, data = sleep_df)

summary(M_2)


# Multilevel linear model aka linear mixed effects ------------------------

M_3 <- lmer(Reaction ~ Days + (Days|Subject), data = sleep_df )
M_3a <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject), data = sleep_df )

summary(M_3)

coef(M_3) # coefficient for individual subjects
ranef(M_3)

confint(M_3) # confidence intervals 

# conditional and marginal R^2 
r2_nakagawa(M_3)

# random intercepts only
M_4 <- lmer(Reaction ~ 1 + Days + (1 |Subject), data = sleep_df )


# random slopes only
M_5 <- lmer(Reaction ~ 1 + Days + (0 + Days|Subject), data = sleep_df )

summary(M_5)


# random slopes and random intercepts BUT no correlation
M_6 <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0 + Days|Subject), data = sleep_df )
# same as ....
M_6a <- lmer(Reaction ~ 1 + Days + (Days||Subject) , data = sleep_df )

summary(M_6)

anova(M_3, M_6) # null hypothesis test comparing M_3 and M_6

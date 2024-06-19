looic_lm <- function(model){
  # leave-one-out cross-validation in
  # a linear model
  
  # the data used in the model is stored as `$model`
  data_df <- model$model
  
  lpd <- function(i){
    # all data except observation i
    training_df <- slice(data_df, -i)
    # observation i
    test_df <- slice(data_df, i)
    # model fit using all data except observation i
    lm_fit <- lm(as.formula(model), data = training_df)
    # predicted mean of outcome, according to `lm_fit`
    # when predictors are those of observation i
    mu <- predict(newdata = test_df, object = lm_fit)
    # predicted sd of outcome, according to `lm_fit`
    # when predictors are those of observation i
    sd <- sqrt(mean(residuals(lm_fit)^2))
    # value of the outcome on observation i
    y <- test_df[[formula.tools::lhs.vars(formula(model))]]
    # lpd: log predictive density
    # log prob of y on held-out observation
    # predicted by model trained on all but
    # one observation, and where predictors
    # have values on held-out observation:
    dnorm(x = y, mean = mu, sd = sd, log = T)
  }
  
  # -2 elpd
  -2 * sum(map_dbl(seq(nrow(data_df)), lpd))

}

looic_lmer <- function(model){
  # leave-one-out cross-validation in
  # a linear model
  
  # the data used in the model is stored as `$model`
  data_df <- model@frame
  
  lpd <- function(i){
    # all data except observation i
    training_df <- slice(data_df, -i)
    # observation i
    test_df <- slice(data_df, i)
    # model fit using all data except observation i
    suppressWarnings(
      {lmer_fit <- lme4::lmer(as.formula(model), REML = FALSE, data = training_df)}
    )
    # predicted mean of outcome, according to `lm_fit`
    # when predictors are those of observation i
    # fixed + random predictions
    mu <- predict(newdata = test_df, object = lmer_fit)
    # predicted sd of outcome, according to `lm_fit`
    # when predictors are those of observation i
    sd <- sigma(lmer_fit)
    # value of the outcome on observation i
    y <- test_df[[formula.tools::lhs.vars(formula(model))]]
    # lpd: log predictive density
    # log prob of y on held-out observation
    # predicted by model trained on all but
    # one observation, and where predictors
    # have values on held-out observation:
    dnorm(x = y, mean = mu, sd = sd, log = T)
  }
  
  # -2 elpd
  -2 * sum(map_dbl(seq(nrow(data_df)), lpd))
  
}

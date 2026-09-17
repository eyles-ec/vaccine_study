library(stats)
library(tidyverse)
library(MASS)
library(broom)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)


#function to fit CITS models, adapted from 
#https://github.com/eyles-ec/common-ambition 
#outcome_var can be changed to which outcome you like, default is hiv test count
#outcome_type can be adjusted to count or rate, which does poisson/NB in the first case or gaussian in the second
#group_var is the name of the group variable
#exposure_var is the population of the specific group selected
#groups_to_include are which comparison you're making with group_var, default is Somali v White
#reference_group allows you to set the ref group for the model of above
#event_date is the intervention date
#outbreaks = TRUE and outbreaks_dates are an optional variable to include dates of flu outbreaks, but the data are monthly
#robust_if_overdispersed allows the function to try negative binomial by default (but there's also an option for quasipoisson)
#fallback family can be set to quasipoisson if desired

#saves the model output, dispersion, family used, and predicted values for use in the counterfactual in a list
#also prints the final selected model if display = TRUE (by default)

#calculates time and period, dynamically constructs the model formula based on what needs to be included. 

fit_cits_model <- function(df, 
                           outcome_var = "Vaccinations",
                           outcome_type = c("count", "rate"),
                           group_var = "group",
                           exposure_var = "Population",
                           groups_to_include = c("Somali", "White"),
                           reference_group = "White",
                           age = "C1",
                           event_date = as.Date("2022-04-01"),
                           outbreaks = FALSE,
                           outbreaks_dates = "outbreak_dates",
                           robust_if_overdispersed = TRUE,
                           fallback_family = "neg_binomial",
                           display = TRUE) {
  

  # Filter and factorise group variable; ensure period is 0/1 integer
  df_filtered <- df %>%
    dplyr::filter(
      .data[[group_var]] %in% groups_to_include
    ) %>%
   dplyr::mutate(
      !!group_var := factor(
        .data[[group_var]],
        levels = groups_to_include
      ),
      !!group_var := stats::relevel(
        .data[[group_var]],
        ref = reference_group
      ),
      
      # months relative to intervention
      time = (
        lubridate::year(month_start) - lubridate::year(event_date)
      ) * 12 +
        (
          lubridate::month(month_start) - lubridate::month(event_date)
        ),
      
      # post-intervention indicator
      period = as.integer(time >= 0),
     )
  
  #filter out by age group
  if (!is.null(age)) {
    df_filtered <- df_filtered %>%
      filter(age_group %in% age)
  }
  
    #add in outbreaks binary if outbreaks are selected
  if (outbreaks) {
   
     df_filtered <- df_filtered %>%
      mutate(
        outbreak = as.integer(month_start %in% outbreaks_dates)
      )
    
  } else {
    
    df_filtered <- df_filtered %>%
      mutate(outbreak = 0L)
    
  }
  
  # CITS mean model + month fixed effects, formula construction
  
  #right side first
  right_side <- paste0(
    "time * ", group_var, " * period"
  )
  
  if (outbreaks) {
    right_side <- paste(right_side, "+ outbreak")
  }
  
  if (outcome_var == "Vaccinations"){
  right_side <- paste0(
    right_side,
    " + offset(log(", exposure_var, "))"
  )
  }
  formula <- stats::as.formula(
    paste(outcome_var, "~", right_side)
  )
  
  
  if (outcome_type == "count"){
    #fit initial Poisson model
    poisson_model <- stats::glm(formula,
                                family = stats::poisson(link = "log"),
                                data = df_filtered)
    
    #Pearson check for overdispersion 
    dispersion <- sum(stats::residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
    
    #Switch to NB (or quasi) if needed
    if (dispersion > 1.5 && robust_if_overdispersed) {
      warning(sprintf("Overdispersion detected (dispersion = %.2f). Switching to %s model.",
                      dispersion, fallback_family))
      
      model <- switch(fallback_family,
                      quasipoisson = stats::glm(formula,
                                                family = stats::quasipoisson(link = "log"),
                                                data = df_filtered),
                      neg_binomial = MASS::glm.nb(formula, data = df_filtered),
                      stop("Unsupported fallback_family. 'quasipoisson' or 'neg_binomial' only in function call."))
    } else {
      model <- poisson_model
    }
  
    #Add fitted values + pointwise CI for mean (link-scale SE)
    pred <- stats::predict(model, type = "link", se.fit = TRUE)
    df_filtered <- df_filtered %>%
      dplyr::mutate(
        yhat = exp(pred$fit),
        yhat_lower = exp(pred$fit - 1.96 * pred$se.fit),
        yhat_upper = exp(pred$fit + 1.96 * pred$se.fit)
      )
  }
  
  if (outcome_type == "rate"){
    
    #if it's a percentage/rate then fit a gaussian model
    model <- stats::glm(formula,
                        family = stats::gaussian(),
                        data = df_filtered)
    
        #add fitted values + pointwise CI for mean (response-scale SE)
    pred <- stats::predict(model, type = "response", se.fit = TRUE)
    df_filtered <- df_filtered %>%
        dplyr::mutate(
        yhat = pred$fit,
        yhat_lower = pred$fit - 1.96 * pred$se.fit,
        yhat_upper = pred$fit + 1.96 * pred$se.fit
      )
  }
    
  if (display) {
    cat("\n--- Final Model Summary (", stats::family(model)$family, ") ---\n", sep = "")
    print(summary(model))
  }
  
  list(
    data = df_filtered,
    model = model,
   # dispersion = dispersion,
    family_used = stats::family(model)$family
  )
}


#function to generate the counterfactual (if there was no intervention)
#needs the dataframe and model from CITS function above
#we generate this, then append it on the list output from the CITS function for the rest of the pipeline
#also possible to look at the counterfactual summary with the summary option on (by default, or off (summary =FALSE))
#because of autocorrelation we use Newey-West calculation for SE/CIs, nw_lag sets the lag for the calculation (in months)

generate_counterfactual <- function(df, model,
                                    outcome_type = c("count", "rate"),
                                    group_name  = "Somali",
                                    outcome_var = "IncidencePct",
                                    time_var    = "time",
                                    group_var   = "group",
                                    period_var  = "period",
                                    level       = 0.95,
                                    nw_lag    = 3,
                                    summary = TRUE) {
  
  #check arguments and columns since we are doing multiple outcomes/comparisons
  needed_cols <- c(outcome_var, time_var, group_var, period_var)
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("generate_counterfactual(): df is missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  outcome_type <- match.arg(outcome_type)
  
  #calculate constant for CIs 
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha/2)
  
  #evaluate offset() terms on newdata (works for glm / glm.nb)
  compute_offset <- function(model, newdata) {
    mf  <- stats::model.frame(stats::terms(model), data = newdata, na.action = stats::na.pass)
    off <- stats::model.offset(mf)
    if (is.null(off)) rep(0, nrow(newdata)) else as.numeric(off)
  }
  
  #Add row id to merge in counterfactuals later to the correct weeks
  df2 <- df |>
    dplyr::mutate(.row_id = dplyr::row_number())
  
  #Post-period rows for Somali/Weston Deprived/Etc 
  post <- df2 |>
    dplyr::filter(.data[[group_var]] == group_name,
                  .data[[period_var]] == 1)
  
  #output if there's no post rows (unlikely, but just in case)
  if (nrow(post) == 0) {
    out0 <- list(
      data = df2 |> dplyr::select(-.row_id),
      summary_table = dplyr::tibble(
        weeks_post = 0L,
        vcov_used = NA_character_
      ),
      cf_data = post
    )
    if (summary) dplyr::glimpse(out0$summary_table)
    return(out0)
  }
  
  #Offsets (same exposure, same rows)
  offset_vec <- compute_offset(model, post)
  
  #Model matrices (aligned with coefficients), in order to calculate counterfactuals
  tt <- stats::delete.response(stats::terms(model)) #term names
  X  <- stats::model.matrix(tt, post) #design matrix
  
  #Generate counterfactual matrix, e.g. 
  #What would group of interest look like post-intervention if they followed the same post change
  #as the other group + common effects, without the group specific trend
  # Keep post-period rows, but remove the treated group's post differential:
  # 1) group level change in post: group:period
  # 2) group slope change in post: time:group:period
  X_cf <- X
  
  #collect field names
  group_token <- paste0(group_var, group_name)
  col_level_diff <- paste0(group_token, ":", period_var)
  col_slope_diff <- paste0(time_var, ":", group_token, ":", period_var)
  
  #treated group's post differential
  cols_to_zero <- c(col_level_diff, col_slope_diff)
  
  #check for missing columns in case of error with names, throw error up
  missing_cols2 <- setdiff(cols_to_zero, colnames(X_cf))
  if (length(missing_cols2) > 0) {
    stop(
      "Controlled counterfactual could not find these columns in model matrix:\n  ",
      paste(missing_cols2, collapse = "\n  "),
      "\n\nAvailable columns include:\n  ",
      paste(colnames(X_cf), collapse = ", ")
    )
  }
  
  #remove group trend from counterfactual model matrix
  X_cf[, cols_to_zero] <- 0
  
  #generate variance-covariance matrix using Newey West technique,
  #which is more robust for autocorrelation present in data. (default for this is 4 weeks)
  #For some outcomes, AR1 prewhitening does not work, so we wrap it in a tryCatch, which allows
  #us to automatically fallback to not using prewhitening. 
  #define two simple functions
  vcov_try_1 <- function() {
    sandwich::NeweyWest(model, lag = nw_lag, prewhite = TRUE, adjust = TRUE)
    
  } 
  vcov_try_2 <- function() {
    sandwich::NeweyWest(model, lag = nw_lag, prewhite = FALSE, adjust = TRUE)
  }
  
  #calculate variance covariance matrix - > if prewhitening doesn't work, automatically move on
  V_full <- tryCatch(vcov_try_1(), error = function(e1) {
    message("NeweyWest (prewhite=TRUE) failed: ", conditionMessage(e1),
            " -> Trying prewhite=FALSE")
    tryCatch(vcov_try_2(), error = function(e2) {
      stop("Both prewhite=TRUE and prewhite=FALSE NeweyWest() failed:\n",
           conditionMessage(e2))
    })
  })
  
  #repeat to get labels (probably there's a better way to do this..), running the same quietly
  vcov_used <- tryCatch({
    invisible(vcov_try_1()); paste0("NW_lag_", nw_lag, "_prewhite_TRUE")
  }, error = function(e1) {
    paste0("NW_lag_", nw_lag, "_prewhite_FALSE")
  })
  
  #collect model coefficients
  beta_full <- stats::coef(model)
  
  #ensure all the columns match between the cf and the vcov matrix via intersect
  common <- intersect(colnames(X), names(beta_full))
  common <- intersect(common, rownames(V_full))
  
  #keep the common columns to ensure the matrices match 
  X    <- X[, common, drop = FALSE] #design matrix
  X_cf <- X_cf[, common, drop = FALSE] #counterfactual matrix
  beta <- beta_full[common] #coefficients from model
  V2   <- V_full[common, common, drop = FALSE] #variance-covariance matrix 
  
  #including the offset, calculate the mean and counterfacual mean
  #exponentiated due to model output
  #again, CF is what would group have looked like post-intervention if intervention hadn't happened,
  #but the common post effect, monthly dummies, and baseline trends still happened (e.g. other group)
  #offset_vec should be all 0s if it's a gaussian model so don't need to do anything fancy
  eta_hat_post <- as.numeric(X %*% beta) + offset_vec
  eta_cf_post  <- as.numeric(X_cf %*% beta) + offset_vec
  
  if (outcome_type == "count") {
    yhat_post <- exp(eta_hat_post)
    cf_post   <- exp(eta_cf_post)
  } else {
    yhat_post <- eta_hat_post
    cf_post   <- eta_cf_post
  }
  
  y_obs      <- post[[outcome_var]] #extract observed counts
  weeks_post <- nrow(post) #used for calculating annualised/weekly differences
  
  #calculate totals for summary table 
  total_observed       <- sum(y_obs, na.rm = TRUE)
  total_counterfactual <- sum(cf_post, na.rm = TRUE)
  total_difference     <- sum(y_obs - cf_post, na.rm = TRUE)
  
  total_predicted  <- sum(yhat_post, na.rm = TRUE)
  total_pred_diff  <- total_predicted - total_counterfactual
  
  #calculate percent/annualised/weekly differences
  #if/else structure avoids division by 0 and breaking the code
  percent_increase   <- if (total_counterfactual > 0) 100 * total_difference / total_counterfactual else NA_real_
  pred_perc_increase <- if (total_counterfactual > 0) 100 * total_pred_diff / total_counterfactual else NA_real_
  
  annualised_difference     <- if (weeks_post > 0) total_difference * 52 / weeks_post else NA_real_
  average_weekly_difference <- if (weeks_post > 0) total_difference / weeks_post else NA_real_
  
  #use delta method to calculate slopes with matrix crossproducts
  #for reference https://doi.org/10.1111/1475-6773.12122
  
  if (outcome_type == "count") {
    #poisson/nb
    g_cf   <- as.numeric(crossprod(X_cf, cf_post))  # d sum(cf) / d beta
    g_pred <- as.numeric(crossprod(X, yhat_post))# d sum(pred) / d beta
    
  } else {
    #gaussian/rate
    g_cf   <- colSums(X_cf)
    g_pred <- colSums(X)
    
  }
  
  #use the crossproducts to calculate standard errors 
  se_total_cf   <- sqrt(drop(t(g_cf)   %*% V2 %*% g_cf))
  se_total_pred <- sqrt(drop(t(g_pred) %*% V2 %*% g_pred))
  
  #observed total is non-random, so the se comes from modelled counterfactual
  se_total_diff <- se_total_cf
  
  #predicted difference se (using matrix crossproducts)
  g_pred_diff <- g_pred - g_cf
  se_total_pred_diff <- sqrt(drop(t(g_pred_diff) %*% V2 %*% g_pred_diff))
  
  #the se of percent increase depends on the counterfactual, using similar if/else logic to the annualised/weekly
  if (is.finite(total_counterfactual) && total_counterfactual > 0) {
    g_percent_inc <- 100 * (-total_observed / (total_counterfactual^2)) * g_cf
    se_percent_inc <- sqrt(drop(t(g_percent_inc) %*% V2 %*% g_percent_inc))
    
    g_pred_percent_inc <- 100 * ((1/total_counterfactual) * g_pred -
                                   (total_predicted/(total_counterfactual^2)) * g_cf)
    se_pred_percent_inc <- sqrt(drop(t(g_pred_percent_inc) %*% V2 %*% g_pred_percent_inc))
  } else {
    se_percent_inc <- NA_real_
    se_pred_percent_inc <- NA_real_
  }
  
  #logic as above, for se (avoids dividing by 0)
  se_annualised_diff <- if (weeks_post > 0) se_total_diff * 52 / weeks_post else NA_real_
  se_avg_weekly_diff <- if (weeks_post > 0) se_total_diff / weeks_post else NA_real_
  
  #mini helper function to calculate cis for everything, with error catching for possible induced NAs above
  ci <- function(est, se, lower_bound = -Inf) {
    if (!is.finite(est) || !is.finite(se)) return(c(NA_real_, NA_real_))
    c(max(lower_bound, est - z * se), est + z * se)
  }
  
  #generate cis for cf, differences, etc. 
  ci_total_cf        <- ci(total_counterfactual, se_total_cf, lower_bound = 0)
  ci_total_diff      <- ci(total_difference, se_total_diff)
  ci_total_pred      <- ci(total_predicted, se_total_pred, lower_bound = 0)
  ci_total_pred_diff <- ci(total_pred_diff, se_total_pred_diff)
  
  ci_percent_inc  <- ci(percent_increase, se_percent_inc)
  ci_pred_percent <- ci(pred_perc_increase, se_pred_percent_inc)
  ci_annualised   <- ci(annualised_difference, se_annualised_diff)
  ci_avg_weekly   <- ci(average_weekly_difference, se_avg_weekly_diff)
  
  #join per-week counterfactual back into the data using the row ids we made earlier
  cf_tbl <- dplyr::tibble(.row_id = post$.row_id, cf = cf_post)
  
  #add cf back into the original data, makes it easier to plot 
  out <- df2 |>
    dplyr::left_join(cf_tbl, by = ".row_id") |>
    dplyr::mutate(
      cf = dplyr::if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1, cf, NA_real_),
      cf_diff = dplyr::if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1,
                               .data[[outcome_var]] - cf, NA_real_)
    ) |>
    dplyr::select(-.row_id)
  
  #also return a tidy summary table for reporting 
  summary_df <- dplyr::tibble(
    total_observed = total_observed,
    total_counterfactual = total_counterfactual,
    total_counterfactual_lower = ci_total_cf[1],
    total_counterfactual_upper = ci_total_cf[2],
    total_difference = total_difference,
    total_difference_lower = ci_total_diff[1],
    total_difference_upper = ci_total_diff[2],
    total_predicted = total_predicted,
    total_predicted_lower = ci_total_pred[1],
    total_predicted_upper = ci_total_pred[2],
    total_pred_diff = total_pred_diff,
    total_pred_diff_lower = ci_total_pred_diff[1],
    total_pred_diff_upper = ci_total_pred_diff[2],
    percent_increase = percent_increase,
    percent_increase_lower = ci_percent_inc[1],
    percent_increase_upper = ci_percent_inc[2],
    pred_perc_increase = pred_perc_increase,
    pred_perc_increase_lower = ci_pred_percent[1],
    pred_perc_increase_upper = ci_pred_percent[2],
    weeks_post = weeks_post,
    annualised_difference = annualised_difference,
    annualised_difference_lower = ci_annualised[1],
    annualised_difference_upper = ci_annualised[2],
    average_weekly_difference = average_weekly_difference,
    average_weekly_difference_lower = ci_avg_weekly[1],
    average_weekly_difference_upper = ci_avg_weekly[2],
    vcov_used = vcov_used
  )
  
  #stick it all together in a list
  cf_data <- out |>
    dplyr::filter(.data[[group_var]] == group_name,
                  .data[[period_var]] == 1) |>
    dplyr::select(dplyr::all_of(c(time_var, group_var, period_var)), cf)
  
  result <- list(data = out, summary_table = summary_df, cf_data = cf_data)
  
  #print an easy summary of the cf table if you want 
  if (summary) dplyr::glimpse(result$summary_table)
  
  return(result)
}


#create a standardised CITS plot, including the counterfactual line for Bristol ACHC
#can save the plot if desired with own filename, overwrite protection. you also can scale it
#height and width are in inches by default 

plot_cits <- function(cits_model,
                      outcome_var = "Vaccinations",       
                      outcome_label = "Vaccinations (count)",    
                      group_var = "group",
                      groups_to_include = c("Somali", "White"),
                      main_group = "Somali",
                      save_plot = FALSE,
                      filename = "./plots/cits_plot.png",
                      overwrite = FALSE,
                      width = 10,
                      height = 6,
                      dpi = 300) {
  
  df_plot <- cits_model$data
  
  # Dynamically extract the group variable
  df_plot <- df_plot %>%
    mutate(group = .data[[group_var]]) %>%
    filter(group %in% groups_to_include)
  
  # Counterfactual line only for main group
  cf_data <- df_plot %>%
    filter(group == main_group, !is.na(cf))
  
  # Build plot
  p <- ggplot(df_plot, aes(x = time, y = .data[[outcome_var]], color = group)) +
    geom_point(alpha = 0.5)+
    geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper, fill = group),
                alpha = 0.2, color = NA) +
    geom_line(aes(y = yhat), size = 1.2) +
    geom_line(data = cf_data,
              aes(y = cf), linetype = "solid", color = "grey40", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      #title = paste("Observed vs Predicted", outcome_label, "Over Time"),
      #subtitle = "Dashed line = Interventions (Month of April 1, 2022)\nDotted line = counterfactual prediction for Bristol ACHC",
      x = "Months Since Intervention",
      y = paste0(outcome_label, " (count)"),
      color = "Group",
      fill = "Group"
    ) +
    facet_wrap(~ group, scales = "free_y") +
    theme_minimal()
  
  # Save plot if requested
  if (save_plot) {
    if (file.exists(filename) && !overwrite) {
      warning(paste("File", filename, "already exists. Set overwrite = TRUE to replace it."))
    } else {
      ggsave(filename, plot = p, width = width, height = height, dpi = dpi)
    }
  }
  
  return(p)
}

#get working directory from paths.R
setwd(wd)

#read in data, fix dates
vax <- read_excel("Maindata_13_08_26.xlsx", sheet = "BNSSG_Monthly")
vax$Month <- as.Date(as.character(vax$Month))

#reshape to long format for modelling
#change to NAs to coerce numeric
vax <- vax %>%
  mutate(
    across(
      -c(month_start, Month),
      ~ as.numeric(na_if(as.character(.), "*"))
    )
  )

#long format for cits
vax_long <- vax %>%
  
pivot_longer(cols = -c(month_start, Month),
              names_to = c(".value", "age_group", "group"),
              names_pattern = "(Vaccinations|Population|UptakePct|IncidencePct)_(C\\d+)_(.*)"
) 

rm(vax)

#add in intervention dates - this will be useful for when we make this into a loop
event_dates <- tribble(
  ~group,                ~event_date,
  "Somali",              as.Date("2023-04-01"),
  "OtherMinority",       as.Date("2023-04-01"),
  "WestonDeprived",      as.Date("2023-04-01")
)

#outbreaks dates

#how to add in months in sequence, see below to append individual dates
outbreak_dates <- c(
  seq(
    as.Date("2023-06-01"),
    as.Date("2024-03-01"),
    by = "month"
  )
)

#commented out, but append individual months as such
# outbreak_dates <- c(outbreak_dates, as.Date(c(
#   "2023-06-01",
#   "2023-07-01"
# )
#)
#)

#missingness
vax_long <- vax_long %>%
  mutate(
    Vaccinations = ifelse(
      is.na(Vaccinations),
      round(UptakePct * Population / 100),
      Vaccinations
    )
  )

#fit model. arguments are described in comments above the function
#make sure you have the right outcome type. 
#the outbreaks don't seem to add anything. i think maybe it's overparameterised? 

cits <- fit_cits_model(vax_long, outcome_var = "UptakePct", age ="C1", outcome_type = "rate", outbreaks = TRUE)

#take a peek

ggplot(cits$data,
       aes(Month, UptakePct,
           colour = group)) +
  geom_point() +
  geom_line(aes(y = yhat), linewidth = 1)

#generate cf, no need to filter by age group as if using the cits object and data it is the filtered version used for the model

cf <- generate_counterfactual(
  df = cits$data,
  model = cits$model,
  group_name  = "Somali",
  outcome_var =  "UptakePct",
  outcome_type = "rate",
  time_var    = "time",
  group_var   = "group",
  period_var  = "period",
  summary = TRUE
)

#create cf graph 

p <- plot_cits(cf,save_plot = FALSE)

p


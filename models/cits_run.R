library(stats)
library(tidyverse)
library(MASS)
library(broom)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(marginaleffects)
library(emmeans)


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
                           event_date = as.Date("2023-06-01"),
                           outbreaks = FALSE,
                           outbreaks_dates = outbreak_dates,
                           robust_if_overdispersed = FALSE,
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
      
      # months post intervention
      
      time_post = pmax(time, 0),
      
      # post-intervention indicator
      period = as.integer(time >= 0)
     )
  
  #filter out by age group
  if (!is.null(age)) {
    df_filtered <- df_filtered %>%
      filter(age_group %in% age)
  }
  
    #add in outbreaks binary if outbreaks are selected
  if (outbreaks) {
   
    if (age == "C1") {
      
      outbreak_start <- min(outbreaks_date)
      
      df_filtered <- df_filtered %>%
        mutate(
          outbreak = as.integer(month_start >= outbreak_start)
        )
      
    
    } else {
      
      df_filtered <- df_filtered %>%
        mutate(
          outbreak = as.integer(month_start %in% outbreaks_dates)
        )
    }
  
    
  } else {
    
    df_filtered <- df_filtered %>%
      mutate(outbreak = 0L)
    
  }
  
  # CITS formula construction, with separations for age and outcome type
  #right side first
  if (age == "C1") {
    #group + time + period + time_after + group x time + group x period + group x time_post
    
    right_side <- paste0(
      group_var, "* (time + period + time_post)"
    )
    
    if (outbreaks) {
      right_side <- paste0(right_side, "+ outbreak")
    }
    
  
   
  } else {
    #group + time + outbreak (pulse) + group x time + group x outbreak (pulse)
      
    right_side <- paste0(
      group_var, " * time"
    )
    
    if (outbreaks) {
      right_side <- paste0(
        group_var, " * (time + outbreak)"
      )
    }
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
  
  cat("\nFormula being fitted:\n")
  print(formula)
  
  cat("\nComparison and cohort:\n")
  print(paste0(groups_to_include, " ", age))
  
  
  if (outcome_type == "count"){
    #fit initial Poisson model
    poisson_model <- stats::glm(formula,
                                family = stats::poisson(link = "log"),
                                data = df_filtered)
    
    #Pearson check for overdispersion 
    dispersion <- sum(stats::residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
    
    cat("\nDispersion:\n")
    print(dispersion)
    
    #Switch to NB (or quasi) if needed
    if (dispersion > 2 && robust_if_overdispersed) {
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
    dispersion = dispersion,
    age_cohort = age,
    family_used = stats::family(model)$family
  )
}


#function to calculate crude vaccine event rates (not vaccination rates, since we have incident counts)
#goes by period, outbreak, or overall rates, depending on age group 

calculate_crude_rates <- function(
    df,
    age,
    groups_to_include,
    stratify_by = c("period", "outbreak", "overall"),
    outcome_var = "Vaccinations",
    exposure_var = "Population",
    age_var = "age_group",
    group_var = "group",
    multiplier = 1000
) {
  
  stratify_by <- match.arg(stratify_by)
  
  if (stratify_by == "overall") {
    
    out<- df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        events = sum(.data[[outcome_var]], na.rm = TRUE),
        population = sum(.data[[exposure_var]], na.rm = TRUE),
        crude_rate = events / population * multiplier,
        .groups = "drop"
      )
    
  } else {
    
    out <- df %>%
      dplyr::group_by(
        group,
        .data[[stratify_by]]
      ) %>%
      dplyr::summarise(
        events = sum(.data[[outcome_var]], na.rm = TRUE),
        population = sum(.data[[exposure_var]], na.rm = TRUE),
        crude_rate = events / population * multiplier,
        .groups = "drop"
      )
    
  }
  
  out
  
}

#estimated marginal means
#returns means, interactions, 
#we do not set default arguments here as they are drawn from the model output being used in the function call in the main code

generate_marginal_means <- function(
    model,
    age,
    comparison,
    response = TRUE
) {
  
  if (age == "C1") {
    
    emm <- emmeans(
      model,
      ~ group | period,
      type = "response"
    )
    
    contrasts <- contrast(
      emm,
      method = "revpairwise"
    )
    
    interaction_emm <- emmeans(
      model,
      ~ group * period,
      type = "response"
    )
    
    interaction_contrasts <- contrast(
      interaction_emm,
      interaction = "pairwise"
    )
    
  } else {
    
    emm <- emmeans(
      model,
      ~ group | outbreak,
      type = "response"
    )
    
    contrasts <- contrast(
      emm,
      method = "revpairwise"
    )
    
    interaction_emm <- emmeans(
      model,
      ~ group * outbreak,
      type = "response"
    )
    
    interaction_contrasts <- contrast(
      interaction_emm,
      interaction = "pairwise"
    )
    
  }
  
  list(
    emm = as.data.frame(emm),
    contrasts = as.data.frame(summary(contrasts)),
    interaction_contrasts = as.data.frame(summary(interaction_contrasts)),
    age = age,
    comparison = comparison
  )
  
}

#create a standardised CITS plot, including the counterfactual line for Bristol ACHC
#can save the plot if desired with own filename, overwrite protection. you also can scale it
#height and width are in inches by default 

plot_cits <- function(cits_model,
                      outcome_var = "Vaccinations",       
                      outcome_label = "Vaccinations (count)",    
                      group_var = "group",
                      groups_to_include = c("Somali", "White"),
                      group_name = "Somali",
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
    filter(group == group_name, !is.na(cf))
  
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
vax$month_start <- as.Date(vax$month_start)

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

#outbreaks dates

#how to add in months in sequence, see below to append individual dates
outbreak_dates <- c(
  seq(
    as.Date("2023-09-01"),
    as.Date("2024-04-01"),
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

#set model specifications to run the model function with purrr (pmap)
model_specs <- tibble::tribble(
  ~age, ~groups_to_include, ~reference_group, ~outbreaks, ~event_date,
  "C1", c("Somali", "White"), "White", FALSE, as.Date("2023-06-01"),
  "C2", c("Somali", "White"), "White", TRUE,  as.Date("2023-06-01"),
  "C3", c("Somali", "White"), "White", TRUE,  as.Date("2023-06-01"),
  "C1", c("WestonDeprived", "NonWestonDeprived"), "NonWestonDeprived", FALSE, as.Date("2023-06-01"),
  "C2", c("WestonDeprived", "NonWestonDeprived"), "NonWestonDeprived", TRUE,  as.Date("2023-06-01"),
  "C3", c("WestonDeprived", "NonWestonDeprived"), "NonWestonDeprived", TRUE,  as.Date("2023-06-01")
)

model_output <- pmap(
  model_specs,
  function(age,
           groups_to_include,
           reference_group,
           outbreaks,
           event_date){
    
    fit_cits_model(
      df = vax_long,
      outcome_var = "Vaccinations",
      outcome_type = "count",
      group_var = "group",
      exposure_var = "Population",
      groups_to_include = groups_to_include,
      reference_group = reference_group,
      age = age,
      event_date = event_date,
      outbreaks = outbreaks,
      outbreaks_dates = outbreak_dates,
      robust_if_overdispersed = TRUE,
      fallback_family = "neg_binomial",
      display = FALSE
    )
    
  }
)

#map outputs onto a dataframe using map2_dfr from purrr

coef_output <- purrr::map2_dfr(
  model_output,
  seq_along(model_output),
  ~{
    
    broom::tidy(
      .x$model,
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
      mutate(
        age = model_specs$age[.y],
        groups = paste(
          model_specs$groups_to_include[[.y]],
          collapse = " vs "
        ),
        event_date = model_specs$event_date[.y],
        outbreaks = model_specs$outbreaks[.y],
        dispersion = .x$dispersion,
        aic = .x$model$aic,
        model_family = .x$family_used,
        formula = as.character(paste(deparse(formula(.x$model)), collapse = ""))
      )
    
  }
)

#calculate crude rates
#these are interpreted as for ex, 24 vaccination events per 1,000 children in the population during the study period
#we can also calculate separately for before/after and outbreak/not as stratified 

#over the entire study period
crude_rates_overall <- purrr::map2_dfr(
  model_output,
  seq_along(model_output),
  ~{
    
    calculate_crude_rates(
      df = .x$data,
      stratify_by = "overall"
    ) %>%
      dplyr::mutate(
        comparison = paste(
          model_specs$groups_to_include[[.y]],
          collapse = " vs "
        ),
        age = model_specs$age[.y],
        stratification = "Overall"
      )
    
  }
)

crude_rates_stratified <- purrr::map2_dfr(
  model_output,
  seq_along(model_output),
  ~{
    
    stratified_by <- if (
      model_specs$age[.y] == "C1"
    ) {
      "period"
    } else {
      "outbreak"
    }
    
    calculate_crude_rates(
      df = .x$data,
      stratify_by = stratified_by
    ) %>%
      dplyr::mutate(
        comparison = paste(
          model_specs$groups_to_include[[.y]],
          collapse = " vs "
        ),
        age = model_specs$age[.y],
        stratification = stratified_by
      )
    
  }
)

crude_rates_all <- dplyr::bind_rows(
  crude_rates_overall,
  crude_rates_stratified
) %>%
  dplyr::relocate(
    comparison,
    age,
    stratification
  )

#generate estimated marginal means, interaction contrasts, pairwise contrasts

mm_output <- purrr::map2(
  model_output,
  seq_along(model_output),
  ~ generate_marginal_means(
    model = .x$model,
    age = model_specs$age[.y],
    comparison = paste(
      model_specs$groups_to_include[[.y]],
      collapse = " vs "
    )
  )
)

mm_means <- purrr::map_dfr(
  mm_output,
  ~{
    .x$emm %>%
      dplyr::mutate(
        comparison = .x$comparison,
        age = .x$age,
        output_type = "Marginal Mean"
      )
  }
) %>%
  dplyr::relocate(
    output_type,
    comparison,
    age
  )


mm_contrasts <- purrr::map_dfr(
  mm_output,
  ~{
    .x$contrasts %>%
      dplyr::mutate(
        comparison = .x$comparison,
        age = .x$age,
        output_type = "Pairwise Contrast"
      )
  }
) %>%
  dplyr::relocate(
    output_type,
    comparison,
    age
  )

mm_interactions <- purrr::map_dfr(
  mm_output,
  ~{
    .x$interaction_contrasts %>%
      dplyr::mutate(
        comparison = .x$comparison,
        age = .x$age,
        output_type = "Interaction Contrast"
      )
  }
) %>%
  dplyr::relocate(
    output_type,
    comparison,
    age
  )

#create results directory if it doesn't exist
if (!dir.exists("./results")) {
  dir.create("./results", recursive = TRUE)
}

#write results to directory 
write.csv(
  coef_output,
  "./results/cits_model_coefficients.csv",
  row.names = FALSE
)

write.csv(
  crude_rates_all,
  "./results/crude_event_rates.csv",
  row.names = FALSE
)

write.csv(
  mm_means, 
  "./results/mm_means.csv",
  row.names = FALSE
  )

write.csv(
  mm_contrasts,
  "./results/mm_contrasts.csv",
  row.names = FALSE
  )

write.csv(
  mm_interactions, 
  "./results/mm_interactions.csv",
    row.names = FALSE
  )
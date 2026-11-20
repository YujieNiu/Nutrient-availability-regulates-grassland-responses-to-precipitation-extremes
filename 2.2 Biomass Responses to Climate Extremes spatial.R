
# ============================================================
# Biomass response versus site precipitation
#
# Separate mixed-effects models are fitted for:
# - Extreme dry events
# - Extreme wet events
#
# Model within each extreme event type:
# LRR_biom ~ precipitation_100 * trt
#
# Estimate represents the change in LRR_biom associated with
# a 100-mm increase in site precipitation.
# ============================================================

rm(list = ls(all.names = TRUE))


# ============================================================
# Load required packages
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(purrr)
library(patchwork)


# ============================================================
# Define the data directory
# ============================================================

dir.data <- ""


if (!dir.exists(dir.data)) {
  stop(
    "The data directory does not exist: ",
    dir.data
  )
}


# ============================================================
# Read the biomass-response data
# ============================================================

data.file <- file.path(
  dir.data,
  "Biomass response result.csv"
)

if (!file.exists(data.file)) {
  stop(
    "The input file does not exist: ",
    data.file
  )
}

recovery_scores <- readr::read_csv(
  data.file,
  show_col_types = FALSE
)


# ============================================================
# Check required columns
# ============================================================

required.columns <- c(
  "site_code",
  "block",
  "plot",
  "trt",
  "extreme_type",
  "LRR_biom",
  "mean_prec_2000_2020"
)

missing.columns <- setdiff(
  required.columns,
  names(recovery_scores)
)

if (length(missing.columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing.columns, collapse = ", "),
    "\nAvailable columns: ",
    paste(names(recovery_scores), collapse = ", ")
  )
}


# ============================================================
# Clean the complete response dataset
# ============================================================

# Records with missing LRR values are retained in dat_all for
# data-availability checks, but they cannot contribute to model
# estimation.

dat_all <- recovery_scores %>%
  dplyr::mutate(
    site_code = trimws(
      as.character(site_code)
    ),
    
    block = trimws(
      as.character(block)
    ),
    
    plot = trimws(
      as.character(plot)
    ),
    
    trt = trimws(
      as.character(trt)
    ),
    
    extreme_type = trimws(
      as.character(extreme_type)
    ),
    
    LRR_biom = suppressWarnings(
      as.numeric(
        as.character(LRR_biom)
      )
    ),
    
    mean_prec_2000_2020 = suppressWarnings(
      as.numeric(
        as.character(mean_prec_2000_2020)
      )
    )
  ) %>%
  dplyr::filter(
    trt %in% c(
      "Control",
      "NPK"
    ),
    
    extreme_type %in% c(
      "Extreme dry",
      "Extreme wet"
    ),
    
    !is.na(site_code),
    site_code != "",
    
    !is.na(block),
    block != "",
    
    !is.na(plot),
    plot != "",
    
    !is.na(mean_prec_2000_2020),
    is.finite(mean_prec_2000_2020)
  )


# ============================================================
# Count all sites in the response dataset
# ============================================================

n_sites_in_file <- dat_all %>%
  dplyr::summarise(
    n_sites = dplyr::n_distinct(site_code)
  ) %>%
  dplyr::pull(n_sites)

cat(
  "Number of sites in the biomass-response dataset: ",
  n_sites_in_file,
  "\n",
  sep = ""
)


# ============================================================
# Summarize LRR availability by site
# ============================================================

site_availability <- dat_all %>%
  dplyr::group_by(site_code) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    
    n_valid_LRR = sum(
      !is.na(LRR_biom) &
        is.finite(LRR_biom)
    ),
    
    n_missing_LRR = sum(
      is.na(LRR_biom) |
        !is.finite(LRR_biom)
    ),
    
    contributes_to_model =
      n_valid_LRR > 0,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    contributes_to_model,
    site_code
  )

print(
  site_availability,
  n = Inf
)

readr::write_csv(
  site_availability,
  file.path(
    dir.data,
    "Biomass_LRR_site_availability.csv"
  )
)


# ============================================================
# Create the effective model dataset
# ============================================================

# Only records with a valid response can be used by lmer().
# Sites are not manually removed.

dat_mixed <- dat_all %>%
  dplyr::filter(
    !is.na(LRR_biom),
    is.finite(LRR_biom)
  ) %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    extreme_type = factor(
      extreme_type,
      levels = c(
        "Extreme dry",
        "Extreme wet"
      )
    ),
    
    site_code = factor(site_code),
    block = factor(block),
    plot = factor(plot),
    
    # Unique block identifier
    block_id = interaction(
      site_code,
      block,
      drop = TRUE
    ),
    
    # Unique experimental-plot identifier
    plot_id = interaction(
      site_code,
      block,
      plot,
      trt,
      drop = TRUE
    )
  ) %>%
  droplevels()

if (nrow(dat_mixed) == 0) {
  stop(
    "No valid LRR observations remain for model fitting."
  )
}


# ============================================================
# Report the effective model dataset
# ============================================================

n_sites_in_model <- dat_mixed %>%
  dplyr::summarise(
    n_sites = dplyr::n_distinct(site_code)
  ) %>%
  dplyr::pull(n_sites)

cat(
  "Number of sites contributing valid LRR observations: ",
  n_sites_in_model,
  "\n",
  sep = ""
)

cat(
  "Number of valid LRR observations: ",
  nrow(dat_mixed),
  "\n",
  sep = ""
)

cat(
  "Number of experimental plots: ",
  dplyr::n_distinct(dat_mixed$plot_id),
  "\n",
  sep = ""
)


# ============================================================
# Check site-level precipitation
# ============================================================


site_precipitation_check <- dat_mixed %>%
  dplyr::group_by(site_code) %>%
  dplyr::summarise(
    n_precipitation_values =
      dplyr::n_distinct(mean_prec_2000_2020),
    
    precipitation_min = min(
      mean_prec_2000_2020,
      na.rm = TRUE
    ),
    
    precipitation_max = max(
      mean_prec_2000_2020,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_precipitation_values > 1
  )

if (nrow(site_precipitation_check) > 0) {
  print(
    site_precipitation_check,
    n = Inf
  )
  
  stop(
    paste0(
      "Some sites contain more than one ",
      "mean_prec_2000_2020 value."
    )
  )
}


# ============================================================
# Center and scale site precipitation
# ============================================================

# A common center is retained for the two separate models.
# One unit of precipitation_100 represents 100 mm.

precipitation_mean <- mean(
  dat_mixed$mean_prec_2000_2020,
  na.rm = TRUE
)

dat_mixed <- dat_mixed %>%
  dplyr::mutate(
    precipitation_100 = (
      mean_prec_2000_2020 -
        precipitation_mean
    ) / 100
  )

cat(
  "Mean precipitation used for centering: ",
  round(precipitation_mean, 2),
  " mm\n",
  sep = ""
)

cat(
  "Observed site-precipitation range: ",
  paste(
    round(
      range(
        dat_mixed$mean_prec_2000_2020,
        na.rm = TRUE
      ),
      2
    ),
    collapse = " to "
  ),
  " mm\n",
  sep = ""
)


# ============================================================
# Calculate sample sizes
# ============================================================

sample_size_table <- dat_mixed %>%
  dplyr::group_by(
    extreme_type,
    trt
  ) %>%
  dplyr::summarise(
    n_sites = dplyr::n_distinct(
      site_code
    ),
    
    n_obs = dplyr::n(),
    
    n_blocks = dplyr::n_distinct(
      block_id
    ),
    
    n_plots = dplyr::n_distinct(
      plot_id
    ),
    
    .groups = "drop"
  )

print(
  sample_size_table,
  n = Inf
)


# ============================================================
# Confirm that both treatments are available in each event type
# ============================================================

combination_check <- tidyr::expand_grid(
  extreme_type = factor(
    c(
      "Extreme dry",
      "Extreme wet"
    ),
    levels = levels(
      dat_mixed$extreme_type
    )
  ),
  
  trt = factor(
    c(
      "Control",
      "NPK"
    ),
    levels = levels(
      dat_mixed$trt
    )
  )
) %>%
  dplyr::left_join(
    sample_size_table,
    by = c(
      "extreme_type",
      "trt"
    )
  )

if (any(is.na(combination_check$n_obs))) {
  print(
    combination_check,
    n = Inf
  )
  
  stop(
    paste0(
      "At least one extreme-event × treatment combination ",
      "has no valid observations."
    )
  )
}


# ============================================================
# Split the data into dry and wet datasets
# ============================================================

mixed_data_list <- split(
  dat_mixed,
  dat_mixed$extreme_type,
  drop = TRUE
)

print(
  names(mixed_data_list)
)


# ============================================================
# Function for fitting one event-specific model
# ============================================================

fit_event_model <- function(data) {
  
  data <- data %>%
    droplevels()
  
  if (dplyr::n_distinct(data$trt) < 2) {
    stop(
      "Both Control and NPK are required in each event model."
    )
  }
  
  # Sum-to-zero contrasts are used for the Type III tests.
  contrasts(data$trt) <- stats::contr.sum(
    nlevels(data$trt)
  )
  
  lmerTest::lmer(
    LRR_biom ~
      precipitation_100 *
      trt +
      (1 | site_code) +
      (1 | block_id) +
      (1 | plot_id),
    data = data,
    REML = TRUE,
    na.action = na.exclude,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
}


# ============================================================
# Fit separate dry and wet mixed-effects models
# ============================================================

models <- purrr::map(
  mixed_data_list,
  fit_event_model
)


# ============================================================
# Display model summaries and Type III tests
# ============================================================

for (event_name in names(models)) {
  
  current_model <- models[[event_name]]
  
  cat("\n")
  cat("====================================================\n")
  cat("Mixed-effects model: ", event_name, "\n", sep = "")
  cat("====================================================\n")
  
  print(
    summary(current_model)
  )
  
  cat("\nType III tests:\n")
  
  current_anova <- stats::anova(
    current_model,
    type = 3,
    ddf = "Satterthwaite"
  )
  
  print(current_anova)
  
  singular_flag <- lme4::isSingular(
    current_model,
    tol = 1e-4
  )
  
  cat(
    "\nSingular fit: ",
    singular_flag,
    "\n",
    sep = ""
  )
  
  if (singular_flag) {
    cat("\nRandom-effect variances:\n")
    
    print(
      lme4::VarCorr(current_model)
    )
  }
}


# ============================================================
# Export fixed-effect coefficients
# ============================================================

# ============================================================
# Define formatting functions
# ============================================================

format_three_decimals <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  dplyr::if_else(
    is.na(x),
    NA_character_,
    formatC(
      x,
      format = "f",
      digits = 3
    )
  )
}


format_p_value <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    
    x < 0.001 ~ formatC(
      x,
      format = "e",
      digits = 3
    ),
    
    TRUE ~ formatC(
      x,
      format = "f",
      digits = 3
    )
  )
}


# ============================================================
# Extract fixed-effect results
# ============================================================

fixed_effect_results_raw <- purrr::imap_dfr(
  models,
  function(model, event_name) {
    
    # Extract the data used in the model
    model_data <- stats::model.frame(
      model
    )
    
    # Number of sites included in the model
    n_sites_model <- dplyr::n_distinct(
      model_data$site_code
    )
    
    # Number of observations included in the model
    n_obs_model <- stats::nobs(
      model
    )
    
    broom.mixed::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE,
      conf.level = 0.95
    ) %>%
      dplyr::transmute(
        extreme_type = event_name,
        
        term = as.character(term),
        
        Estimate = as.numeric(
          estimate
        ),
        
        SE = as.numeric(
          std.error
        ),
        
        DF = as.numeric(
          df
        ),
        
        CI_lower = as.numeric(
          conf.low
        ),
        
        CI_upper = as.numeric(
          conf.high
        ),
        
        t_value = as.numeric(
          statistic
        ),
        
        P_value = as.numeric(
          p.value
        ),
        
        n_sites = as.integer(
          n_sites_model
        ),
        
        n_obs = as.integer(
          n_obs_model
        )
      )
  }
) %>%
  tibble::as_tibble()


# ============================================================
# Add significance information
# ============================================================

fixed_effect_results_raw <- fixed_effect_results_raw %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      is.na(P_value) ~
        "P value unavailable",
      
      P_value < 0.05 ~
        "Significant",
      
      TRUE ~
        "Not significant"
    ),
    
    sig_symbol = dplyr::case_when(
      is.na(P_value) ~ "",
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )


# ============================================================
# Create the formatted fixed-effect table
# ============================================================

fixed_effect_results <- fixed_effect_results_raw %>%
  dplyr::transmute(
    extreme_type = extreme_type,
    
    term = term,
    
    Estimate = format_three_decimals(
      Estimate
    ),
    
    SE = format_three_decimals(
      SE
    ),
    
    DF = format_three_decimals(
      DF
    ),
    
    CI_lower = format_three_decimals(
      CI_lower
    ),
    
    CI_upper = format_three_decimals(
      CI_upper
    ),
    
    t_value = format_three_decimals(
      t_value
    ),
    
    P_value = format_p_value(
      P_value
    ),
    
    significance = significance,
    
    sig_symbol = sig_symbol,
    
    n_sites = as.integer(
      n_sites
    ),
    
    n_obs = as.integer(
      n_obs
    )
  )


# ============================================================
# Display the formatted fixed-effect results
# ============================================================

print(
  fixed_effect_results,
  n = Inf,
  width = Inf
)


# ============================================================
# Export the fixed-effect results as an Excel file
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

fixed_effect_output.xlsx <- file.path(
  dir.data,
  "Biomass_LRR_separate_models_fixed_effects.xlsx"
)

writexl::write_xlsx(
  list(
    `Fixed effects` = fixed_effect_results
  ),
  path = fixed_effect_output.xlsx
)

cat(
  "\nFixed-effect results exported to: ",
  fixed_effect_output.xlsx,
  "\n",
  sep = ""
)






# ============================================================
# Extract precipitation slopes from each model
# ============================================================

extract_event_slopes <- function(
    model,
    event_name
) {
  
  slope_object <- emmeans::emtrends(
    model,
    specs = ~ trt,
    var = "precipitation_100"
  )
  
  slope_raw <- as.data.frame(
    summary(
      slope_object,
      infer = c(TRUE, TRUE),
      level = 0.95,
      adjust = "none",
      null = 0
    )
  )
  
  slope_column <- grep(
    "\\.trend$|^trend$",
    names(slope_raw),
    value = TRUE
  )
  
  lower_column <- grep(
    "^lower\\.CL$|^asymp\\.LCL$",
    names(slope_raw),
    value = TRUE
  )
  
  upper_column <- grep(
    "^upper\\.CL$|^asymp\\.UCL$",
    names(slope_raw),
    value = TRUE
  )
  
  statistic_column <- grep(
    "^t\\.ratio$|^z\\.ratio$",
    names(slope_raw),
    value = TRUE
  )
  
  if (
    length(slope_column) != 1 ||
    length(lower_column) != 1 ||
    length(upper_column) != 1 ||
    length(statistic_column) != 1 ||
    !"SE" %in% names(slope_raw) ||
    !"df" %in% names(slope_raw) ||
    !"p.value" %in% names(slope_raw)
  ) {
    stop(
      "Required emtrends columns could not be identified for ",
      event_name,
      ".\nAvailable columns: ",
      paste(
        names(slope_raw),
        collapse = ", "
      )
    )
  }
  
  slope_raw %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      extreme_type = event_name,
      
      trt = as.character(trt),
      
      Estimate = as.numeric(
        .data[[slope_column]]
      ),
      
      SE = as.numeric(SE),
      
      DF = as.numeric(df),
      
      CI_lower = as.numeric(
        .data[[lower_column]]
      ),
      
      CI_upper = as.numeric(
        .data[[upper_column]]
      ),
      
      t_value = as.numeric(
        .data[[statistic_column]]
      ),
      
      P_value = as.numeric(
        p.value
      )
    )
}


# ============================================================
# Create the final slope table
# ============================================================

slope_results <- purrr::imap_dfr(
  models,
  function(model, event_name) {
    extract_event_slopes(
      model = model,
      event_name = event_name
    )
  }
) %>%
  dplyr::left_join(
    sample_size_table %>%
      dplyr::transmute(
        extreme_type =
          as.character(extreme_type),
        
        trt =
          as.character(trt),
        
        n_sites =
          as.integer(n_sites),
        
        n_obs =
          as.integer(n_obs)
      ),
    by = c(
      "extreme_type",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    extreme_type = factor(
      extreme_type,
      levels = c(
        "Extreme dry",
        "Extreme wet"
      )
    ),
    
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    )
  ) %>%
  dplyr::arrange(
    extreme_type,
    trt
  ) %>%
  tibble::as_tibble()

cat("\n")
cat("====================================================\n")
cat("Treatment-specific precipitation slopes\n")
cat("Estimate: change in LRR_biom per 100 mm\n")
cat("====================================================\n")

print(
  slope_results,
  n = Inf
)

readr::write_csv(
  slope_results,
  file.path(
    dir.data,
    "Biomass_LRR_precipitation_slopes_separate_models.csv"
  )
)


# ============================================================
# Compare Control and NPK slopes within each event model
# ============================================================

slope_contrast_results <- purrr::imap_dfr(
  models,
  function(model, event_name) {
    
    event_trends <- emmeans::emtrends(
      model,
      specs = ~ trt,
      var = "precipitation_100"
    )
    
    as.data.frame(
      summary(
        pairs(
          event_trends,
          adjust = "none"
        ),
        infer = c(TRUE, TRUE),
        level = 0.95
      )
    ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        extreme_type = event_name,
        .before = 1
      )
  }
) %>%
  tibble::as_tibble()

cat("\n")
cat("====================================================\n")
cat("Control versus NPK slope contrasts\n")
cat("====================================================\n")

print(
  slope_contrast_results,
  n = Inf
)




# ============================================================
# Function for model-based predictions
# ============================================================

predict_event_model <- function(
    model,
    event_name,
    model_data,
    precipitation_center,
    number_of_points = 200
) {
  
  precipitation_values <- seq(
    min(
      model_data$mean_prec_2000_2020,
      na.rm = TRUE
    ),
    
    max(
      model_data$mean_prec_2000_2020,
      na.rm = TRUE
    ),
    
    length.out = number_of_points
  )
  
  precipitation_100_values <- (
    precipitation_values -
      precipitation_center
  ) / 100
  
  prediction_object <- emmeans::emmeans(
    model,
    specs =
      ~ trt *
      precipitation_100,
    at = list(
      precipitation_100 =
        precipitation_100_values
    )
  )
  
  prediction_raw <- as.data.frame(
    summary(
      prediction_object,
      infer = c(TRUE, FALSE),
      level = 0.95
    )
  )
  
  lower_column <- grep(
    "^lower\\.CL$|^asymp\\.LCL$",
    names(prediction_raw),
    value = TRUE
  )
  
  upper_column <- grep(
    "^upper\\.CL$|^asymp\\.UCL$",
    names(prediction_raw),
    value = TRUE
  )
  
  if (
    length(lower_column) != 1 ||
    length(upper_column) != 1 ||
    !"emmean" %in% names(prediction_raw)
  ) {
    stop(
      "Prediction columns could not be identified for ",
      event_name,
      ".\nAvailable columns: ",
      paste(
        names(prediction_raw),
        collapse = ", "
      )
    )
  }
  
  prediction_raw %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      extreme_type = event_name,
      
      trt = as.character(trt),
      
      precipitation_100 =
        as.numeric(precipitation_100),
      
      mean_prec_2000_2020 =
        precipitation_100 * 100 +
        precipitation_center,
      
      predicted_LRR =
        as.numeric(emmean),
      
      prediction_SE =
        as.numeric(SE),
      
      lower_CL = as.numeric(
        .data[[lower_column]]
      ),
      
      upper_CL = as.numeric(
        .data[[upper_column]]
      )
    )
}


# ============================================================
# Generate predictions from the two separate models
# ============================================================

prediction_data <- purrr::imap_dfr(
  models,
  function(model, event_name) {
    
    predict_event_model(
      model = model,
      event_name = event_name,
      model_data =
        mixed_data_list[[event_name]],
      precipitation_center =
        precipitation_mean,
      number_of_points = 200
    )
  }
) %>%
  dplyr::left_join(
    slope_results %>%
      dplyr::transmute(
        extreme_type =
          as.character(extreme_type),
        
        trt =
          as.character(trt),
        
        slope_Estimate = Estimate,
        slope_P_value = P_value
      ),
    by = c(
      "extreme_type",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    line_significance = dplyr::case_when(
      is.na(slope_P_value) ~
        "P unavailable",
      
      slope_P_value < 0.05 ~
        "P < 0.05",
      
      TRUE ~
        "P >= 0.05"
    ),
    
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    extreme_type = factor(
      extreme_type,
      levels = c(
        "Extreme dry",
        "Extreme wet"
      )
    )
  ) %>%
  dplyr::arrange(
    extreme_type,
    trt,
    mean_prec_2000_2020
  ) %>%
  tibble::as_tibble()


# ============================================================
# Export model predictions
# ============================================================

readr::write_csv(
  prediction_data,
  file.path(
    dir.data,
    "Biomass_LRR_separate_model_predictions.csv"
  )
)


# ============================================================
# Define plotting colors
# ============================================================

climate_colors <- c(
  "Extreme dry" = "#ff9f00",
  "Extreme wet" = "#1f78b4"
)


# ============================================================
# Prepare sample-size labels
# ============================================================

n_summary <- sample_size_table %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    extreme_type = factor(
      extreme_type,
      levels = c(
        "Extreme dry",
        "Extreme wet"
      )
    )
  )

print(
  n_summary,
  n = Inf
)


# ============================================================
# Define the plotting function
# ============================================================

plot_mixed_model <- function(
    treatment_name,
    model_predictions,
    sample_sizes,
    x_limits = c(0, 1100),
    y_limits = c(-2, 1),
    show_y_title = TRUE
) {
  
  treatment_title <- if (
    treatment_name == "NPK"
  ) {
    "NPK\u03bc"
  } else {
    "Control"
  }
  
  treatment_title_color <- if (
    treatment_name == "NPK"
  ) {
    "#984ea3"
  } else {
    "#33a02c"
  }
  
  prediction_subset <- model_predictions %>%
    dplyr::filter(
      as.character(trt) ==
        treatment_name
    ) %>%
    droplevels()
  
  if (nrow(prediction_subset) == 0) {
    stop(
      "No predictions were found for treatment: ",
      treatment_name
    )
  }
  
  n_label_subset <- sample_sizes %>%
    dplyr::filter(
      as.character(trt) ==
        treatment_name
    ) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        as.character(extreme_type) ==
          "Extreme dry" ~
          paste0(
            "Extreme dry: n = ",
            n_obs
          ),
        
        as.character(extreme_type) ==
          "Extreme wet" ~
          paste0(
            "Extreme wet: n = ",
            n_obs
          ),
        
        TRUE ~
          paste0(
            "n = ",
            n_obs
          )
      ),
      
      label_x =
        x_limits[2] -
        0.03 * diff(x_limits),
      
      label_y = dplyr::case_when(
        as.character(extreme_type) ==
          "Extreme dry" ~
          y_limits[1] +
          0.07 * diff(y_limits),
        
        as.character(extreme_type) ==
          "Extreme wet" ~
          y_limits[1] +
          0.14 * diff(y_limits),
        
        TRUE ~
          y_limits[1]
      )
    )
  
  ggplot2::ggplot() +
    
    # Zero-response line
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.6
    ) +
    
    # Model-based 95% confidence intervals
    ggplot2::geom_ribbon(
      data = prediction_subset,
      mapping = ggplot2::aes(
        x = mean_prec_2000_2020,
        ymin = lower_CL,
        ymax = upper_CL,
        fill = extreme_type,
        group = extreme_type
      ),
      alpha = 0.18,
      color = NA
    ) +
    
    # Fixed-effect prediction lines
    ggplot2::geom_line(
      data = prediction_subset,
      mapping = ggplot2::aes(
        x = mean_prec_2000_2020,
        y = predicted_LRR,
        color = extreme_type,
        linetype = line_significance,
        group = extreme_type
      ),
      linewidth = 1.4
    ) +
    
    # Sample-size labels
    ggplot2::geom_text(
      data = n_label_subset,
      mapping = ggplot2::aes(
        x = label_x,
        y = label_y,
        label = label,
        color = extreme_type
      ),
      hjust = 1,
      vjust = 0,
      size = 4,
      show.legend = FALSE
    ) +
    
    ggplot2::scale_color_manual(
      values = climate_colors,
      breaks = c(
        "Extreme dry",
        "Extreme wet"
      ),
      name = "Extreme event"
    ) +
    
    ggplot2::scale_fill_manual(
      values = climate_colors,
      breaks = c(
        "Extreme dry",
        "Extreme wet"
      ),
      name = "Extreme event"
    ) +
    
    # Solid lines indicate slopes significantly different from zero.
    # Dashed lines indicate nonsignificant slopes.
    ggplot2::scale_linetype_manual(
      values = c(
        "P < 0.05" = "solid",
        "P >= 0.05" = "dashed",
        "P unavailable" = "dotted"
      ),
      guide = "none"
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = seq(
        x_limits[1],
        x_limits[2],
        by = 200
      )
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = seq(
        y_limits[1],
        y_limits[2],
        by = 1
      )
    ) +
    
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    
    ggplot2::labs(
      title = treatment_title,
      x = "Site precipitation (mm)",
      y = if (show_y_title) {
        "Biomass response (LRR)"
      } else {
        NULL
      }
    ) +
    
    ggplot2::theme_classic(
      base_size = 15
    ) +
    
    ggplot2::theme(
      panel.grid =
        ggplot2::element_blank(),
      
      axis.line =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.7
        ),
      
      axis.ticks =
        ggplot2::element_line(
          color = "black",
          linewidth = 0.7
        ),
      
      axis.text =
        ggplot2::element_text(
          size = 12,
          color = "black"
        ),
      
      axis.title =
        ggplot2::element_text(
          size = 15,
          face = "plain"
        ),
      
      plot.title =
        ggplot2::element_text(
          hjust = 0,
          size = 13,
          face = "plain",
          color =
            treatment_title_color
        ),
      
      legend.position = "bottom",
      
      legend.title =
        ggplot2::element_blank(),
      
      plot.margin =
        ggplot2::margin(
          t = 5,
          r = 8,
          b = 5,
          l = 8
        )
    )
}


# ============================================================
# Create the Control panel
# ============================================================

p_control <- plot_mixed_model(
  treatment_name = "Control",
  model_predictions = prediction_data,
  sample_sizes = n_summary,
  x_limits = c(0, 1100),
  y_limits = c(-2, 1),
  show_y_title = TRUE
)

print(p_control)


# ============================================================
# Create the NPK panel
# ============================================================

p_npk <- plot_mixed_model(
  treatment_name = "NPK",
  model_predictions = prediction_data,
  sample_sizes = n_summary,
  x_limits = c(0, 1100),
  y_limits = c(-2, 1),
  show_y_title = FALSE
)

print(p_npk)


# ============================================================
# Combine the two panels
# ============================================================

final_plot <- (
  p_control +
    p_npk +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    
    axis.text.y =
      ggplot2::element_text(
        size = 12,
        color = "black"
      ),
    
    axis.ticks.y =
      ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
    
    axis.line.y =
      ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      )
  )

print(final_plot)


# ============================================================
# Save the final figure
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "LRR_biom_vs_prec_separate_models.tiff"
  ),
  plot = final_plot,
  width = 8,
  height = 6.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ============================================================
# Final summary
# ============================================================

cat(
  "\nSeparate dry and wet mixed-effects analyses completed.\n"
)

cat(
  "Sites in the response dataset: ",
  n_sites_in_file,
  "\n",
  sep = ""
)

cat(
  "Sites contributing valid LRR observations: ",
  n_sites_in_model,
  "\n",
  sep = ""
)


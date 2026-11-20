

# ============================================================
# Figure 2a-b
# Standardized biomass responses to precipitation anomalies
#
# Linear and quadratic mixed-effects models are compared for:
# 1. Overall precipitation gradient
# 2. Dry range: dif.re <= 0
# 3. Wet range: dif.re > 0
# ============================================================


# ============================================================
# Clear the workspace
# ============================================================

rm(list = ls(all.names = TRUE))


# ============================================================
# Load required packages
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(patchwork)
library(purrr)
library(tibble)


# ============================================================
# Define the data directory
# ============================================================

dir.data <- ""


if (!dir.exists(dir.data)) {
  stop("The data directory does not exist: ", dir.data)
}


# ============================================================
# Read the complete data
# ============================================================

data.file <- file.path(
  dir.data,
  "full_data_with_climate.csv"
)

if (!file.exists(data.file)) {
  stop("The input file does not exist: ", data.file)
}

df <- readr::read_csv(
  data.file,
  show_col_types = FALSE
)


# ============================================================
# Check required variables
# ============================================================

required.columns <- c(
  "site_code",
  "year",
  "block",
  "plot",
  "trt",
  "year_trt",
  "dif.re",
  "vascular_live_mass"
)

missing.columns <- setdiff(
  required.columns,
  names(df)
)

if (length(missing.columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing.columns, collapse = ", ")
  )
}


# ============================================================
# Prepare biomass data
# ============================================================

df_rel1 <- df %>%
  dplyr::mutate(
    year = suppressWarnings(
      as.integer(as.character(year))
    ),
    
    year_trt = suppressWarnings(
      as.numeric(as.character(year_trt))
    ),
    
    dif.re = suppressWarnings(
      as.numeric(as.character(dif.re))
    ),
    
    vascular_live_mass = suppressWarnings(
      as.numeric(as.character(vascular_live_mass))
    ),
    
    trt = trimws(
      as.character(trt)
    )
  ) %>%
  
  # Retain the two treatments and remove the NPK pretreatment year
  dplyr::filter(
    trt %in% c("Control", "NPK"),
    
    !(
      trt == "NPK" &
        !is.na(year_trt) &
        year_trt == 0
    )
  ) %>%
  
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c("Control", "NPK")
    ),
    
    # Log-transform biomass
    log_biomass = log10(
      vascular_live_mass + 1
    )
  ) %>%
  
  dplyr::filter(
    !is.na(site_code),
    !is.na(year),
    !is.na(block),
    !is.na(plot),
    
    !is.na(dif.re),
    is.finite(dif.re),
    
    !is.na(log_biomass),
    is.finite(log_biomass)
  ) %>%
  
  # Standardize biomass within each site and treatment
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  
  dplyr::mutate(
    mean_biomass = mean(
      log_biomass,
      na.rm = TRUE
    ),
    
    sd_biomass = stats::sd(
      log_biomass,
      na.rm = TRUE
    ),
    
    Log_Biomass_z = dplyr::if_else(
      !is.na(sd_biomass) &
        is.finite(sd_biomass) &
        sd_biomass > 0,
      
      (log_biomass - mean_biomass) /
        sd_biomass,
      
      NA_real_
    )
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::filter(
    !is.na(Log_Biomass_z),
    is.finite(Log_Biomass_z)
  ) %>%
  
  dplyr::mutate(
    site_code = factor(site_code),
    block = factor(block),
    plot = factor(plot),
    
    # Unique plot identifier
    plot_id = interaction(
      site_code,
      block,
      plot,
      trt,
      drop = TRUE
    )
  ) %>%
  
  droplevels()


# ============================================================
# Check removed NPK pretreatment observations
# ============================================================

remaining_npk_year0 <- df_rel1 %>%
  dplyr::filter(
    trt == "NPK",
    !is.na(year_trt),
    year_trt == 0
  )

if (nrow(remaining_npk_year0) > 0) {
  stop(
    "NPK observations with year_trt == 0 remain in the analysis."
  )
}

cat(
  "Remaining NPK observations with year_trt == 0:",
  nrow(remaining_npk_year0),
  "\n"
)


# ============================================================
# Check the retained data
# ============================================================

cat(
  "Number of observations:",
  nrow(df_rel1),
  "\n"
)

cat(
  "Number of sites:",
  dplyr::n_distinct(df_rel1$site_code),
  "\n"
)

cat(
  "Number of plots:",
  dplyr::n_distinct(df_rel1$plot_id),
  "\n"
)


# ============================================================
# Save the standardized biomass data
# ============================================================

readr::write_csv(
  df_rel1,
  file.path(
    dir.data,
    "Relative_Biomass_dual.csv"
  )
)


# ============================================================
# Prepare the three model subsets
# ============================================================

# The same x-axis limits are used for modeling and plotting.

df_model <- df_rel1 %>%
  dplyr::filter(
    dif.re > -2.3,
    dif.re < 3
  )

# Overall range
df_all <- df_model

# Dry range, including zero
df_neg <- df_model %>%
  dplyr::filter(
    dif.re <= 0
  )

# Wet range, excluding zero
df_pos <- df_model %>%
  dplyr::filter(
    dif.re > 0
  )


# ============================================================
# Check the three data subsets
# ============================================================

subset_summary <- dplyr::bind_rows(
  df_all %>%
    dplyr::mutate(range = "Overall"),
  
  df_neg %>%
    dplyr::mutate(range = "Dry: dif.re <= 0"),
  
  df_pos %>%
    dplyr::mutate(range = "Wet: dif.re > 0")
) %>%
  dplyr::group_by(
    range,
    trt
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_sites = dplyr::n_distinct(site_code),
    n_plots = dplyr::n_distinct(plot_id),
    n_x_values = dplyr::n_distinct(dif.re),
    x_min = min(dif.re, na.rm = TRUE),
    x_max = max(dif.re, na.rm = TRUE),
    .groups = "drop"
  )

print(subset_summary, n = Inf)


# ============================================================
# Function to compare linear and quadratic models
# ============================================================

# Model selection procedure:
#
# 1. Fit the linear and quadratic models using maximum likelihood.
# 2. Compare them using AIC and a likelihood-ratio test.
# 3. Select the quadratic model only when:
#      AIC_linear - AIC_quadratic >= 2
#    and
#      likelihood-ratio P < 0.05
# 4. Otherwise select the linear model.
# 5. Refit the selected model using REML.

compare_and_select_model <- function(
    data,
    subset_name
) {
  
  # ----------------------------------------------------------
  # Check the data
  # ----------------------------------------------------------
  
  data <- data %>%
    dplyr::filter(
      !is.na(Log_Biomass_z),
      is.finite(Log_Biomass_z),
      !is.na(dif.re),
      is.finite(dif.re),
      !is.na(trt),
      !is.na(site_code),
      !is.na(plot_id)
    ) %>%
    droplevels()
  
  if (nrow(data) < 10) {
    stop(
      "Insufficient observations for subset: ",
      subset_name
    )
  }
  
  if (dplyr::n_distinct(data$dif.re) < 3) {
    stop(
      "Fewer than three unique precipitation values for subset: ",
      subset_name
    )
  }
  
  if (dplyr::n_distinct(data$trt) < 2) {
    stop(
      "Both treatments are not available for subset: ",
      subset_name
    )
  }
  
  
  # ----------------------------------------------------------
  # Center precipitation within this subset
  # ----------------------------------------------------------
  
  x_center <- mean(
    data$dif.re,
    na.rm = TRUE
  )
  
  data <- data %>%
    dplyr::mutate(
      dif_c = dif.re - x_center
    )
  
  
  # ----------------------------------------------------------
  # Define candidate model formulas
  # ----------------------------------------------------------
  
  linear_formula <- stats::as.formula(
    paste(
      "Log_Biomass_z ~ trt * dif_c",
      "+ (1 | site_code)",
      "+ (1 | plot_id)"
    )
  )
  
  quadratic_formula <- stats::as.formula(
    paste(
      "Log_Biomass_z ~ trt *",
      "(dif_c + I(dif_c^2))",
      "+ (1 | site_code)",
      "+ (1 | plot_id)"
    )
  )
  
  
  # ----------------------------------------------------------
  # Fit candidate models using maximum likelihood
  # ----------------------------------------------------------
  
  linear_ml <- tryCatch(
    lmerTest::lmer(
      linear_formula,
      data = data,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) {
      message(
        subset_name,
        ": linear model failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(linear_ml)) {
    stop(
      "The linear model could not be fitted for: ",
      subset_name
    )
  }
  
  quadratic_ml <- tryCatch(
    lmerTest::lmer(
      quadratic_formula,
      data = data,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) {
      message(
        subset_name,
        ": quadratic model failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Select the model
  # ----------------------------------------------------------
  
  if (is.null(quadratic_ml)) {
    
    aic_linear <- stats::AIC(linear_ml)
    aic_quadratic <- NA_real_
    delta_aic <- NA_real_
    lrt_chisq <- NA_real_
    lrt_df <- NA_real_
    lrt_p <- NA_real_
    
    selected_type <- "linear"
    selection_reason <- "Quadratic model failed"
    
  } else {
    
    aic_linear <- stats::AIC(linear_ml)
    aic_quadratic <- stats::AIC(quadratic_ml)
    
    # Positive values favor the quadratic model
    delta_aic <- aic_linear - aic_quadratic
    
    model_comparison <- as.data.frame(
      stats::anova(
        linear_ml,
        quadratic_ml,
        refit = FALSE
      )
    )
    
    lrt_chisq <- model_comparison$Chisq[2]
    lrt_df <- model_comparison$`Chi Df`[2]
    lrt_p <- model_comparison$`Pr(>Chisq)`[2]
    
    if (
      is.finite(delta_aic) &&
      delta_aic >= 2 &&
      !is.na(lrt_p) &&
      lrt_p < 0.05
    ) {
      
      selected_type <- "quadratic"
      
      selection_reason <- paste0(
        "Quadratic selected: delta AIC = ",
        round(delta_aic, 2),
        "; LRT P = ",
        signif(lrt_p, 3)
      )
      
    } else {
      
      selected_type <- "linear"
      
      selection_reason <- paste0(
        "Linear selected: delta AIC = ",
        round(delta_aic, 2),
        "; LRT P = ",
        signif(lrt_p, 3)
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Refit the selected model using REML
  # ----------------------------------------------------------
  
  selected_formula <- if (
    selected_type == "quadratic"
  ) {
    quadratic_formula
  } else {
    linear_formula
  }
  
  final_model <- lmerTest::lmer(
    selected_formula,
    data = data,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  
  # ----------------------------------------------------------
  # Model-comparison result
  # ----------------------------------------------------------
  
  comparison_result <- tibble::tibble(
    subset = subset_name,
    n = nrow(data),
    n_sites = dplyr::n_distinct(data$site_code),
    n_plots = dplyr::n_distinct(data$plot_id),
    x_min = min(data$dif.re, na.rm = TRUE),
    x_max = max(data$dif.re, na.rm = TRUE),
    x_center = x_center,
    AIC_linear = aic_linear,
    AIC_quadratic = aic_quadratic,
    delta_AIC_linear_minus_quadratic = delta_aic,
    LRT_chisq = lrt_chisq,
    LRT_df = lrt_df,
    LRT_p = lrt_p,
    selected_model = selected_type,
    selection_reason = selection_reason
  )
  
  cat("\n")
  cat("====================================================\n")
  cat("Model comparison:", subset_name, "\n")
  cat("====================================================\n")
  
  print(
    comparison_result,
    n = Inf
  )
  
  cat("\nSelected model:", selected_type, "\n")
  cat(selection_reason, "\n")
  
  list(
    data = data,
    center = x_center,
    linear_ml = linear_ml,
    quadratic_ml = quadratic_ml,
    model = final_model,
    model_type = selected_type,
    comparison = comparison_result
  )
}


# ============================================================
# Compare models for the three data subsets
# ============================================================

fit_all <- compare_and_select_model(
  data = df_all,
  subset_name = "Overall precipitation gradient"
)

fit_neg <- compare_and_select_model(
  data = df_neg,
  subset_name = "Dry range: dif.re <= 0"
)

fit_pos <- compare_and_select_model(
  data = df_pos,
  subset_name = "Wet range: dif.re > 0"
)


# Replace the original subsets with the centered model data
df_all <- fit_all$data
df_neg <- fit_neg$data
df_pos <- fit_pos$data

model_all <- fit_all$model
model_neg <- fit_neg$model
model_pos <- fit_pos$model


# ============================================================
# Combine and export model-comparison results
# ============================================================

model_comparison_results <- dplyr::bind_rows(
  fit_all$comparison,
  fit_neg$comparison,
  fit_pos$comparison
)

print(
  model_comparison_results,
  n = Inf
)

readr::write_csv(
  model_comparison_results,
  file.path(
    dir.data,
    "linear_vs_quadratic_model_comparison.csv"
  )
)


# ============================================================
# Prepare model-comparison results for export
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


# Combine the model-comparison results
model_comparison_results <- dplyr::bind_rows(
  fit_all$comparison,
  fit_neg$comparison,
  fit_pos$comparison
)


# Retain only information relevant to model selection
model_comparison_export <- model_comparison_results %>%
  dplyr::transmute(
    subset = gsub(
      pattern = "dif\\.re",
      replacement = "Pz",
      x = subset
    ),
    
    n = as.integer(n),
    
    n_sites = as.integer(n_sites),
    
    n_plots = as.integer(n_plots),
    
    AIC_linear = format_three_decimals(
      AIC_linear
    ),
    
    AIC_quadratic = format_three_decimals(
      AIC_quadratic
    ),
    
    delta_AIC_linear_minus_quadratic =
      format_three_decimals(
        delta_AIC_linear_minus_quadratic
      ),
    
    LRT_chisq = format_three_decimals(
      LRT_chisq
    ),
    
    LRT_p = format_p_value(
      LRT_p
    ),
    
    selected_model = selected_model,
    
    selection_reason = gsub(
      pattern = "dif\\.re",
      replacement = "Pz",
      x = selection_reason
    )
  )


# ============================================================
# Display the model-comparison table
# ============================================================

print(
  model_comparison_export,
  n = Inf,
  width = Inf
)


# ============================================================
# Export the model-comparison table as Excel
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

model_comparison_output.xlsx <- file.path(
  dir.data,
  "linear_vs_quadratic_model_comparison.xlsx"
)

writexl::write_xlsx(
  list(
    `Model comparison` = model_comparison_export
  ),
  path = model_comparison_output.xlsx
)

cat(
  "\nModel-comparison results exported to: ",
  model_comparison_output.xlsx,
  "\n",
  sep = ""
)









# ============================================================
# Display the selected model results
# ============================================================

selected_models <- list(
  "Overall precipitation gradient" = model_all,
  "Dry range: dif.re <= 0" = model_neg,
  "Wet range: dif.re > 0" = model_pos
)

selected_model_types <- c(
  "Overall precipitation gradient" =
    fit_all$model_type,
  
  "Dry range: dif.re <= 0" =
    fit_neg$model_type,
  
  "Wet range: dif.re > 0" =
    fit_pos$model_type
)

for (model_name in names(selected_models)) {
  
  cat("\n")
  cat("====================================================\n")
  cat(model_name, "\n")
  cat("Selected form:", selected_model_types[[model_name]], "\n")
  cat("====================================================\n")
  
  print(
    summary(
      selected_models[[model_name]]
    )
  )
  
  cat("\nType III tests:\n")
  
  print(
    anova(
      selected_models[[model_name]],
      type = 3
    )
  )
}


# ============================================================
# Export fixed-effect results from selected models
# ============================================================

# Function for formatting numeric values to three decimal places
format_three_decimals <- function(x) {
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


# Function for formatting P values
format_p_value <- function(x) {
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


# Export fixed-effect results from selected models
model_results <- purrr::imap_dfr(
  selected_models,
  function(model, model_name) {
    
    model_data <- stats::model.frame(model)
    
    n_sites_model <- dplyr::n_distinct(
      model_data$site_code
    )
    
    n_obs_model <- stats::nobs(model)
    
    broom.mixed::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE,
      conf.level = 0.95
    ) %>%
      dplyr::transmute(
        model = model_name,
        
        selected_form =
          selected_model_types[[model_name]],
        
        term = term,
        Estimate = estimate,
        SE = std.error,
        DF = df,
        CI_lower = conf.low,
        CI_upper = conf.high,
        t_value = statistic,
        P_value = p.value,
        n_sites = as.integer(n_sites_model),
        n_obs = as.integer(n_obs_model)
      )
  }
) %>%
  tibble::as_tibble() %>%
  
  # Rename dif.re and the centered model variable dif_c as Pz
  dplyr::mutate(
    model = gsub(
      "dif\\.re",
      "Pz",
      model
    ),
    
    term = gsub(
      "dif\\.re",
      "Pz",
      term
    ),
    
    term = gsub(
      "dif_c",
      "Pz",
      term,
      fixed = TRUE
    ),
    
    # Format statistical estimates
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
    
    # Retain sample sizes as integers
    n_sites = as.integer(n_sites),
    n_obs = as.integer(n_obs)
  )


print(
  model_results,
  n = Inf,
  width = Inf
)


readr::write_csv(
  model_results,
  file.path(
    dir.data,
    "selected_mixed_model_results_combined.csv"
  ),
  na = ""
)

# ============================================================
# Export the formatted results as an Excel file
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

output.xlsx <- file.path(
  dir.data,
  "selected_mixed_model_results_combined.xlsx"
)

writexl::write_xlsx(
  list(
    `Mixed model results` = model_results
  ),
  path = output.xlsx
)

cat(
  "\nExcel file exported to: ",
  output.xlsx,
  "\n",
  sep = ""
)


# ============================================================
# Export Type III tests
# ============================================================

# ============================================================
# Define formatting functions
# ============================================================

# Format numeric statistics to three decimal places
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


# Format P values
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
# Extract Type III test results
# ============================================================

type3_results_raw <- purrr::imap_dfr(
  selected_models,
  function(model, model_name) {
    
    # Extract the data used in the model
    model_data <- stats::model.frame(
      model
    )
    
    # Count the number of sites
    n_sites_model <- dplyr::n_distinct(
      model_data$site_code
    )
    
    # Count the number of observations
    n_obs_model <- stats::nobs(
      model
    )
    
    # Extract the Type III ANOVA table
    result <- as.data.frame(
      stats::anova(
        model,
        type = 3
      )
    )
    
    # Convert row names to a term column
    result$term <- rownames(
      result
    )
    
    rownames(result) <- NULL
    
    
    # Check the required columns
    required_type3_columns <- c(
      "Sum Sq",
      "Mean Sq",
      "NumDF",
      "DenDF",
      "F value",
      "Pr(>F)",
      "term"
    )
    
    missing_type3_columns <- setdiff(
      required_type3_columns,
      names(result)
    )
    
    if (length(missing_type3_columns) > 0) {
      stop(
        "The following Type III ANOVA columns are missing for model '",
        model_name,
        "': ",
        paste(
          missing_type3_columns,
          collapse = ", "
        ),
        "\nAvailable columns: ",
        paste(
          names(result),
          collapse = ", "
        )
      )
    }
    
    
    # Prepare the results
    result %>%
      tibble::as_tibble() %>%
      dplyr::transmute(
        model = model_name,
        
        selected_form =
          selected_model_types[[model_name]],
        
        term = as.character(term),
        
        Sum_Sq = as.numeric(
          .data[["Sum Sq"]]
        ),
        
        Mean_Sq = as.numeric(
          .data[["Mean Sq"]]
        ),
        
        NumDF = as.numeric(
          .data[["NumDF"]]
        ),
        
        DenDF = as.numeric(
          .data[["DenDF"]]
        ),
        
        F_value = as.numeric(
          .data[["F value"]]
        ),
        
        P_value = as.numeric(
          .data[["Pr(>F)"]]
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
# Rename precipitation-anomaly terms as Pz
# ============================================================

type3_results_raw <- type3_results_raw %>%
  dplyr::mutate(
    # Replace dif.re in model names
    model = gsub(
      pattern = "dif\\.re",
      replacement = "Pz",
      x = model
    ),
    
    # Replace dif_c in model names
    model = gsub(
      pattern = "dif_c",
      replacement = "Pz",
      x = model,
      fixed = TRUE
    ),
    
    # Replace dif.re in model terms
    term = gsub(
      pattern = "dif\\.re",
      replacement = "Pz",
      x = term
    ),
    
    # Replace dif_c in model terms
    term = gsub(
      pattern = "dif_c",
      replacement = "Pz",
      x = term,
      fixed = TRUE
    )
  )


# ============================================================
# Add significance information
# ============================================================

type3_results_raw <- type3_results_raw %>%
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
# Create the formatted Type III table
# ============================================================

type3_results <- type3_results_raw %>%
  dplyr::transmute(
    model = model,
    
    selected_form = selected_form,
    
    term = term,
    
    Sum_Sq = format_three_decimals(
      Sum_Sq
    ),
    
    Mean_Sq = format_three_decimals(
      Mean_Sq
    ),
    
    # NumDF is retained as an integer
    NumDF = as.integer(
      round(NumDF)
    ),
    
    DenDF = format_three_decimals(
      DenDF
    ),
    
    F_value = format_three_decimals(
      F_value
    ),
    
    P_value = format_p_value(
      P_value
    ),
    
    significance = significance,
    
    sig_symbol = sig_symbol,
    
    # Sample sizes are retained as integers
    n_sites = as.integer(
      n_sites
    ),
    
    n_obs = as.integer(
      n_obs
    )
  )


# ============================================================
# Display the formatted Type III results
# ============================================================

print(
  type3_results,
  n = Inf,
  width = Inf
)


# ============================================================
# Export the Type III results as an Excel file
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

type3_output.xlsx <- file.path(
  dir.data,
  "selected_model_type3_tests.xlsx"
)

writexl::write_xlsx(
  list(
    `Type III tests` = type3_results
  ),
  path = type3_output.xlsx
)

cat(
  "\nType III test results exported to: ",
  type3_output.xlsx,
  "\n",
  sep = ""
)


# ============================================================
# Function for fixed-effect predictions
# ============================================================

predict_selected_model <- function(
    model,
    model_type,
    x_values,
    x_center
) {
  
  prediction_data <- expand.grid(
    dif.re = x_values,
    trt = c("Control", "NPK"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      dif_c = dif.re - x_center,
      
      trt = factor(
        trt,
        levels = c(
          "Control",
          "NPK"
        )
      )
    )
  
  if (model_type == "quadratic") {
    
    design_matrix <- stats::model.matrix(
      ~ trt * (
        dif_c + I(dif_c^2)
      ),
      data = prediction_data
    )
    
  } else {
    
    design_matrix <- stats::model.matrix(
      ~ trt * dif_c,
      data = prediction_data
    )
  }
  
  fixed_coefficients <- lme4::fixef(
    model
  )
  
  covariance_matrix <- as.matrix(
    stats::vcov(model)
  )
  
  design_matrix <- design_matrix[
    ,
    names(fixed_coefficients),
    drop = FALSE
  ]
  
  prediction_data$predicted <- as.numeric(
    design_matrix %*%
      fixed_coefficients
  )
  
  prediction_data$prediction_se <- sqrt(
    rowSums(
      (
        design_matrix %*%
          covariance_matrix
      ) *
        design_matrix
    )
  )
  
  prediction_data %>%
    dplyr::mutate(
      lower_CL =
        predicted -
        1.96 * prediction_se,
      
      upper_CL =
        predicted +
        1.96 * prediction_se,
      
      selected_model = model_type
    )
}


# ============================================================
# Generate predictions using the selected models
# ============================================================

pred_all <- predict_selected_model(
  model = model_all,
  model_type = fit_all$model_type,
  x_values = seq(
    min(df_all$dif.re, na.rm = TRUE),
    max(df_all$dif.re, na.rm = TRUE),
    length.out = 250
  ),
  x_center = fit_all$center
)

pred_neg <- predict_selected_model(
  model = model_neg,
  model_type = fit_neg$model_type,
  x_values = seq(
    min(df_neg$dif.re, na.rm = TRUE),
    max(df_neg$dif.re, na.rm = TRUE),
    length.out = 150
  ),
  x_center = fit_neg$center
)

pred_pos <- predict_selected_model(
  model = model_pos,
  model_type = fit_pos$model_type,
  x_values = seq(
    min(df_pos$dif.re, na.rm = TRUE),
    max(df_pos$dif.re, na.rm = TRUE),
    length.out = 150
  ),
  x_center = fit_pos$center
)


# ============================================================
# Define figure colors and labels
# ============================================================

treatment_colors <- c(
  "Control" = "#33a02c",
  "NPK" = "#984ea3"
)

treatment_labels <- c(
  "Control" = "Control",
  "NPK" = "NPK\u03bc"
)


# ============================================================
# Function for threshold lines and shaded regions
# ============================================================

add_climate_thresholds <- function(plot_object) {
  
  plot_object +
    
    ggplot2::annotate(
      "rect",
      xmin = -1.282,
      xmax = -1,
      ymin = -Inf,
      ymax = Inf,
      fill = "gray80",
      alpha = 0.20
    ) +
    
    ggplot2::annotate(
      "rect",
      xmin = 1,
      xmax = 1.282,
      ymin = -Inf,
      ymax = Inf,
      fill = "gray80",
      alpha = 0.20
    ) +
    
    ggplot2::geom_vline(
      xintercept = c(
        -1.282,
        -1,
        1,
        1.282
      ),
      linetype = "dashed",
      color = "gray55",
      linewidth = 0.4
    )
}


# ============================================================
# General plotting function
# ============================================================

plot_response <- function(
    raw_data,
    prediction_data,
    title,
    x_label,
    x_limits,
    show_y_title = TRUE
) {
  
  p <- ggplot2::ggplot()
  
  p <- add_climate_thresholds(p)
  
  p +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray25",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray25",
      linewidth = 0.5
    ) +
    
    # Raw observations
    ggplot2::geom_point(
      data = raw_data,
      ggplot2::aes(
        x = dif.re,
        y = Log_Biomass_z,
        fill = trt
      ),
      shape = 21,
      color = "white",
      size = 2.8,
      stroke = 0.3,
      alpha = 0.22
    ) +
    
    # Confidence intervals from the selected model
    ggplot2::geom_ribbon(
      data = prediction_data,
      ggplot2::aes(
        x = dif.re,
        ymin = lower_CL,
        ymax = upper_CL,
        fill = trt,
        group = trt
      ),
      alpha = 0.18,
      color = NA
    ) +
    
    # Fixed-effect predictions from the selected model
    ggplot2::geom_line(
      data = prediction_data,
      ggplot2::aes(
        x = dif.re,
        y = predicted,
        color = trt,
        group = trt
      ),
      linewidth = 1.3
    ) +
    
    ggplot2::scale_color_manual(
      values = treatment_colors,
      labels = treatment_labels,
      name = "Treatment"
    ) +
    
    ggplot2::scale_fill_manual(
      values = treatment_colors,
      labels = treatment_labels,
      name = "Treatment"
    ) +
    
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(-4.5, 3.5)
    ) +
    
    ggplot2::labs(
      title = title,
      x = x_label,
      y = if (show_y_title) {
        "Standardized biomass"
      } else {
        NULL
      }
    ) +
    
    ggplot2::theme_classic(
      base_size = 14
    ) +
    
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "plain",
        size = 15
      ),
      
      axis.title = ggplot2::element_text(
        size = 13,
        face = "plain"
      ),
      
      axis.text = ggplot2::element_text(
        size = 11,
        color = "black"
      ),
      
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      plot.margin = ggplot2::margin(
        t = 5,
        r = 8,
        b = 5,
        l = 8
      )
    )
}


# ============================================================
# Sample sizes for the overall panels
# ============================================================

overall_sample_sizes <- df_all %>%
  dplyr::group_by(trt) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

n_control <- overall_sample_sizes %>%
  dplyr::filter(
    trt == "Control"
  ) %>%
  dplyr::pull(n)

n_npk <- overall_sample_sizes %>%
  dplyr::filter(
    trt == "NPK"
  ) %>%
  dplyr::pull(n)

cat(
  "Control observations:",
  n_control,
  "\n"
)

cat(
  "NPK observations:",
  n_npk,
  "\n"
)


# ============================================================
# Function for an overall treatment panel
# ============================================================

plot_overall_treatment <- function(
    treatment_name,
    raw_data,
    prediction_data,
    sample_size,
    show_y_title = TRUE
) {
  
  treatment_title <- if (
    treatment_name == "NPK"
  ) {
    "NPK\u03bc"
  } else {
    "Control"
  }
  
  treatment_color <- treatment_colors[
    treatment_name
  ]
  
  raw_subset <- raw_data %>%
    dplyr::filter(
      trt == treatment_name
    )
  
  prediction_subset <- prediction_data %>%
    dplyr::filter(
      trt == treatment_name
    )
  
  p <- ggplot2::ggplot()
  
  p <- add_climate_thresholds(p)
  
  p +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray25",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray25",
      linewidth = 0.5
    ) +
    
    ggplot2::geom_ribbon(
      data = prediction_subset,
      ggplot2::aes(
        x = dif.re,
        ymin = lower_CL,
        ymax = upper_CL
      ),
      fill = treatment_color,
      alpha = 0.18,
      color = NA
    ) +
    
    ggplot2::geom_point(
      data = raw_subset,
      ggplot2::aes(
        x = dif.re,
        y = Log_Biomass_z
      ),
      shape = 21,
      color = "white",
      fill = treatment_color,
      size = 4,
      stroke = 0.3,
      alpha = 0.22
    ) +
    
    ggplot2::geom_line(
      data = prediction_subset,
      ggplot2::aes(
        x = dif.re,
        y = predicted
      ),
      color = treatment_color,
      linewidth = 1.3
    ) +
    
    ggplot2::annotate(
      geom = "text",
      x = Inf,
      y = -Inf,
      label = paste0(
        "n = ",
        sample_size
      ),
      hjust = 1.15,
      vjust = -0.8,
      size = 4,
      color = "black"
    ) +
    
    ggplot2::coord_cartesian(
      xlim = c(-2.3, 3),
      ylim = c(-4.05, 3.5),
      clip = "on"
    ) +
    
    ggplot2::labs(
      title = treatment_title,
      x = NULL,
      y = if (show_y_title) {
        "Standardized biomass"
      } else {
        NULL
      }
    ) +
    
    ggplot2::theme_classic(
      base_size = 14
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "plain",
        size = 12
      ),
      
      axis.title = ggplot2::element_text(
        face = "plain",
        size = 13
      ),
      
      axis.text = ggplot2::element_text(
        size = 11,
        color = "black"
      ),
      
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      )
    )
}


# ============================================================
# Overall Control and NPK panels
# ============================================================

p_control <- plot_overall_treatment(
  treatment_name = "Control",
  raw_data = df_all,
  prediction_data = pred_all,
  sample_size = n_control,
  show_y_title = TRUE
)

p_npk <- plot_overall_treatment(
  treatment_name = "NPK",
  raw_data = df_all,
  prediction_data = pred_all,
  sample_size = n_npk,
  show_y_title = FALSE
)


# ============================================================
# Combine the overall treatment panels
# ============================================================

combined_horizontal <- (
  p_control |
    p_npk
) +
  patchwork::plot_annotation(
    caption = "Standardized precipitation within sites",
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        vjust = 0,
        size = 15,
        face = "plain"
      )
    )
  )

print(combined_horizontal)


# ============================================================
# Dry and wet response panels
# ============================================================

dry_title <- paste0(
  "Dry response (",
  fit_neg$model_type,
  ")"
)

wet_title <- paste0(
  "Wet response (",
  fit_pos$model_type,
  ")"
)

p_dry <- plot_response(
  raw_data = df_neg,
  prediction_data = pred_neg,
  title = dry_title,
  x_label = paste0(
    "Standardized precipitation\n",
    "(Decreasing precipitation)"
  ),
  x_limits = c(-2.3, 0),
  show_y_title = TRUE
)

p_wet <- plot_response(
  raw_data = df_pos,
  prediction_data = pred_pos,
  title = wet_title,
  x_label = paste0(
    "Standardized precipitation\n",
    "(Increasing precipitation)"
  ),
  x_limits = c(0, 3),
  show_y_title = FALSE
)

print(p_dry)
print(p_wet)


# ============================================================
# Save the figures
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "stand_responses.tiff"
  ),
  plot = combined_horizontal,
  height = 7,
  width = 9.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "dry_responses.tiff"
  ),
  plot = p_dry,
  height = 7,
  width = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "wet_responses.tiff"
  ),
  plot = p_wet,
  height = 7,
  width = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)






# ============================================================
# Extract peak locations and plot predictions from the
# existing overall quadratic mixed-effects model
#
# This section uses:
#   model_all
#   fit_all
#   df_all
#   pred_all
#
# No additional model is fitted.
# ============================================================


# ============================================================
# Check required objects
# ============================================================


required_objects <- c(
  "model_all",
  "fit_all",
  "df_all",
  "pred_all"
)

object_exists <- vapply(
  required_objects,
  function(object_name) {
    exists(
      object_name,
      envir = .GlobalEnv,
      inherits = TRUE
    )
  },
  FUN.VALUE = logical(1)
)

missing_objects <- required_objects[
  !object_exists
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects do not exist: ",
    paste(
      missing_objects,
      collapse = ", "
    ),
    ". Run the overall mixed-effects-model analysis first."
  )
}

cat(
  "All required objects are available: ",
  paste(
    required_objects,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ============================================================
# Confirm that the selected overall model is quadratic
# ============================================================

if (
  is.null(fit_all$model_type) ||
  fit_all$model_type != "quadratic"
) {
  stop(
    "The selected overall model is not quadratic. ",
    "An internal biomass maximum cannot be estimated ",
    "from a linear model."
  )
}

if (
  is.null(fit_all$center) ||
  length(fit_all$center) != 1 ||
  !is.finite(fit_all$center)
) {
  stop(
    "A valid precipitation-centering value was not found ",
    "in fit_all$center."
  )
}


# ============================================================
# Check the prediction data
# ============================================================

required_prediction_columns <- c(
  "dif.re",
  "trt",
  "predicted",
  "lower_CL",
  "upper_CL"
)

missing_prediction_columns <- setdiff(
  required_prediction_columns,
  names(pred_all)
)

if (length(missing_prediction_columns) > 0) {
  stop(
    "pred_all is missing the following columns: ",
    paste(
      missing_prediction_columns,
      collapse = ", "
    )
  )
}

pred_all <- pred_all %>%
  dplyr::mutate(
    trt = factor(
      as.character(trt),
      levels = c(
        "Control",
        "NPK"
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(trt),
    is.finite(dif.re),
    is.finite(predicted),
    is.finite(lower_CL),
    is.finite(upper_CL)
  )

if (nrow(pred_all) == 0) {
  stop(
    "No valid mixed-effects-model predictions are available."
  )
}


# ============================================================
# Extract fixed-effect coefficients from model_all
# ============================================================

fixed_coefficients <- lme4::fixef(
  model_all
)

cat("\nFixed-effect coefficients used for peak calculation:\n")

print(
  fixed_coefficients
)


# ============================================================
# Function for safely extracting a coefficient
# ============================================================

get_model_coefficient <- function(
    possible_names,
    required = TRUE,
    default = 0
) {
  
  matched_names <- possible_names[
    possible_names %in%
      names(fixed_coefficients)
  ]
  
  if (length(matched_names) == 0) {
    
    if (required) {
      stop(
        "The following coefficient could not be found: ",
        paste(
          possible_names,
          collapse = " or "
        ),
        "\nAvailable coefficients are: ",
        paste(
          names(fixed_coefficients),
          collapse = ", "
        )
      )
    }
    
    return(default)
  }
  
  as.numeric(
    fixed_coefficients[
      matched_names[1]
    ]
  )
}


# ============================================================
# Extract coefficients for the Control treatment
#
# The model was fitted using:
# dif_c = dif.re - fit_all$center
# ============================================================

control_intercept <- get_model_coefficient(
  possible_names = "(Intercept)"
)

control_linear <- get_model_coefficient(
  possible_names = "dif_c"
)

control_quadratic <- get_model_coefficient(
  possible_names = "I(dif_c^2)"
)


# ============================================================
# Extract NPK main effect and interaction coefficients
# ============================================================

npk_intercept_effect <- get_model_coefficient(
  possible_names = "trtNPK",
  required = FALSE,
  default = 0
)

npk_linear_interaction <- get_model_coefficient(
  possible_names = c(
    "trtNPK:dif_c",
    "dif_c:trtNPK"
  ),
  required = TRUE
)

npk_quadratic_interaction <- get_model_coefficient(
  possible_names = c(
    "trtNPK:I(dif_c^2)",
    "I(dif_c^2):trtNPK"
  ),
  required = TRUE
)


# ============================================================
# Construct treatment-specific quadratic coefficients
#
# Control:
# y = b0 + b1*dif_c + b2*dif_c^2
#
# NPK:
# y = (b0 + treatment effect)
#   + (b1 + linear interaction)*dif_c
#   + (b2 + quadratic interaction)*dif_c^2
# ============================================================

peak_coefficients <- tibble::tibble(
  trt = c(
    "Control",
    "NPK"
  ),
  
  intercept = c(
    control_intercept,
    
    control_intercept +
      npk_intercept_effect
  ),
  
  linear_coefficient = c(
    control_linear,
    
    control_linear +
      npk_linear_interaction
  ),
  
  quadratic_coefficient = c(
    control_quadratic,
    
    control_quadratic +
      npk_quadratic_interaction
  )
)


# ============================================================
# Calculate treatment-specific peak locations
#
# On the centered precipitation scale:
#
# peak_c = -b / (2*a)
#
# Convert back to the original standardized precipitation scale:
#
# peak_x = peak_c + fit_all$center
# ============================================================

peak_data <- peak_coefficients %>%
  dplyr::mutate(
    peak_c =
      -linear_coefficient /
      (
        2 * quadratic_coefficient
      ),
    
    peak_x =
      peak_c +
      fit_all$center,
    
    peak_y =
      intercept +
      linear_coefficient * peak_c +
      quadratic_coefficient * peak_c^2
  )


# ============================================================
# Add observed precipitation ranges and sample sizes
# ============================================================

observed_ranges <- df_all %>%
  dplyr::mutate(
    trt = as.character(trt)
  ) %>%
  dplyr::filter(
    trt %in% c(
      "Control",
      "NPK"
    ),
    is.finite(dif.re)
  ) %>%
  dplyr::group_by(trt) %>%
  dplyr::summarise(
    observed_xmin = min(
      dif.re,
      na.rm = TRUE
    ),
    
    observed_xmax = max(
      dif.re,
      na.rm = TRUE
    ),
    
    n = dplyr::n(),
    
    .groups = "drop"
  )

peak_data <- peak_data %>%
  dplyr::left_join(
    observed_ranges,
    by = "trt"
  ) %>%
  dplyr::mutate(
    valid_peak =
      is.finite(quadratic_coefficient) &
      quadratic_coefficient < 0 &
      is.finite(peak_x) &
      is.finite(peak_y) &
      is.finite(observed_xmin) &
      is.finite(observed_xmax) &
      peak_x >= observed_xmin &
      peak_x <= observed_xmax &
      peak_x > -2.3 &
      peak_x < 3
  )


# ============================================================
# Check the calculated peak locations
# ============================================================

if (any(!peak_data$valid_peak)) {
  
  invalid_treatments <- peak_data %>%
    dplyr::filter(
      !valid_peak
    ) %>%
    dplyr::pull(trt)
  
  warning(
    "A valid internal maximum was not identified for: ",
    paste(
      invalid_treatments,
      collapse = ", "
    ),
    call. = FALSE
  )
}

peak_data <- peak_data %>%
  dplyr::filter(
    valid_peak
  ) %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    peak_label = dplyr::if_else(
      trt == "NPK",
      
      sprintf(
        "NPK\u03bc peak = %.2f",
        peak_x
      ),
      
      sprintf(
        "Control peak = %.2f",
        peak_x
      )
    ),
    
    # Separate the two labels slightly to reduce overlap
    peak_label_y = dplyr::if_else(
      trt == "NPK",
      peak_y + 0.10,
      peak_y + 0.055
    )
  )

if (nrow(peak_data) == 0) {
  stop(
    "No valid internal biomass peaks were identified."
  )
}


# ============================================================
# Display the mixed-effects-model peak estimates
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PEAKS FROM THE OVERALL QUADRATIC MIXED-EFFECTS MODEL\n")
cat("============================================================\n")

peak_results_display <- peak_data %>%
  dplyr::transmute(
    Treatment = as.character(trt),
    
    Intercept = round(
      intercept,
      4
    ),
    
    Linear_coefficient = round(
      linear_coefficient,
      4
    ),
    
    Quadratic_coefficient = round(
      quadratic_coefficient,
      4
    ),
    
    Peak_precipitation = round(
      peak_x,
      4
    ),
    
    Predicted_peak_biomass = round(
      peak_y,
      4
    ),
    
    Observed_xmin = round(
      observed_xmin,
      4
    ),
    
    Observed_xmax = round(
      observed_xmax,
      4
    ),
    
    n = as.integer(n)
  )

print(
  as.data.frame(
    peak_results_display
  ),
  row.names = FALSE
)


# ============================================================
# Define treatment colors and labels
# ============================================================

treatment_colors <- c(
  "Control" = "#33a02c",
  "NPK" = "#984ea3"
)

treatment_labels <- c(
  "Control" = "Control",
  "NPK" = "NPK\u03bc"
)


# ============================================================
# Prepare sample-size labels
# ============================================================

n_labels <- observed_ranges %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    label = dplyr::if_else(
      trt == "NPK",
      
      paste0(
        "NPK\u03bc: n = ",
        n
      ),
      
      paste0(
        "Control: n = ",
        n
      )
    ),
    
    x = 2.85,
    
    y = dplyr::if_else(
      trt == "Control",
      -0.72,
      -0.62
    )
  )


# ============================================================
# Create the mixed-effects-model figure
#
# The curves, confidence intervals, and peak locations all
# come from the existing overall mixed-effects model.
# ============================================================

p_quadratic_mixed <- ggplot2::ggplot() +
  
  # Horizontal reference line
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray30",
    linewidth = 0.6
  ) +
  
  # Vertical reference line
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray30",
    linewidth = 0.6
  ) +
  
  # Confidence intervals from model_all
  ggplot2::geom_ribbon(
    data = pred_all,
    mapping = ggplot2::aes(
      x = dif.re,
      ymin = lower_CL,
      ymax = upper_CL,
      fill = trt,
      group = trt
    ),
    alpha = 0.25,
    color = NA
  ) +
  
  # Fixed-effect predictions from model_all
  ggplot2::geom_line(
    data = pred_all,
    mapping = ggplot2::aes(
      x = dif.re,
      y = predicted,
      color = trt,
      group = trt
    ),
    linewidth = 1.2
  ) +
  
  # Vertical lines at the estimated peaks
  ggplot2::geom_vline(
    data = peak_data,
    mapping = ggplot2::aes(
      xintercept = peak_x,
      color = trt
    ),
    linetype = "dotdash",
    linewidth = 1,
    show.legend = FALSE
  ) +
  
  # Open circles at the estimated maxima
  ggplot2::geom_point(
    data = peak_data,
    mapping = ggplot2::aes(
      x = peak_x,
      y = peak_y,
      color = trt
    ),
    shape = 21,
    size = 3.8,
    stroke = 1.2,
    fill = "white",
    show.legend = FALSE
  ) +
  
  # Labels showing peak precipitation values
  ggplot2::geom_label(
    data = peak_data,
    mapping = ggplot2::aes(
      x = peak_x,
      y = peak_label_y,
      label = peak_label,
      color = trt
    ),
    size = 3.5,
    fill = "white",
    linewidth = 0.2,
    label.padding = grid::unit(
      0.15,
      "lines"
    ),
    show.legend = FALSE
  ) +
  
  # Sample-size labels
  ggplot2::geom_text(
    data = n_labels,
    mapping = ggplot2::aes(
      x = x,
      y = y,
      label = label,
      color = trt
    ),
    hjust = 1,
    size = 4,
    show.legend = FALSE
  ) +
  
  # Treatment colors
  ggplot2::scale_color_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment"
  ) +
  
  ggplot2::scale_fill_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment"
  ) +
  
  # Axis and legend titles
  ggplot2::labs(
    x = "Standardized growing-season precipitation within sites",
    y = "Standardized biomass",
    color = "Treatment",
    fill = "Treatment"
  ) +
  
  # Displayed ranges
  ggplot2::coord_cartesian(
    xlim = c(
      -2.3,
      3
    ),
    ylim = c(
      -0.8,
      0.5
    ),
    clip = "off"
  ) +
  
  # Figure style
  ggplot2::theme_classic(
    base_size = 15
  ) +
  
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    
    axis.title = ggplot2::element_text(
      size = 15,
      face = "plain"
    ),
    
    axis.text = ggplot2::element_text(
      size = 12,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.7
    ),
    
    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.7
    ),
    
    plot.margin = ggplot2::margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )


# ============================================================
# Display the figure
# ============================================================

print(
  p_quadratic_mixed
)


# ============================================================
# Save the figure
# ============================================================

mixed_peak_figure <- file.path(
  dir.data,
  "p_quadraticWithPeakline_mixed_model.tiff"
)

ggplot2::ggsave(
  filename = mixed_peak_figure,
  plot = p_quadratic_mixed,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat(
  "\nMixed-effects-model peak figure saved to:\n",
  mixed_peak_figure,
  "\n",
  sep = ""
)






# ============================================================
# Supplementary figure
# Site-specific relationships between precipitation and biomass
#
# GAMs are fitted separately for each site and treatment using
# all plot-level observations. A plot-level random intercept is
# included to account for repeated observations from the same
# plot.
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
library(mgcv)


# ============================================================
# Define the data directory
# ============================================================

# dir.data <- "/Volumes/U5K-Q8X/Science advances/rui_New"
dir.data <- "S:/Science advances/rui_New"


if (!dir.exists(dir.data)) {
  stop(
    "The data directory does not exist: ",
    dir.data
  )
}


# ============================================================
# Read the standardized biomass data
# ============================================================

data.file <- file.path(
  dir.data,
  "Relative_Biomass_dual.csv"
)

if (!file.exists(data.file)) {
  stop(
    "The input file does not exist: ",
    data.file
  )
}

df_raw <- readr::read_csv(
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
  "year",
  "trt",
  "dif.re",
  "Log_Biomass_z"
)

missing.columns <- setdiff(
  required.columns,
  names(df_raw)
)

if (length(missing.columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing.columns, collapse = ", "),
    "\nAvailable columns: ",
    paste(names(df_raw), collapse = ", ")
  )
}


# ============================================================
# Define the complete site list
# ============================================================

# The site list is created before filtering invalid observations.
# This allows facet_wrap(drop = FALSE) to retain every site panel.

all_site_levels <- df_raw %>%
  dplyr::filter(
    !is.na(site_code),
    trimws(as.character(site_code)) != ""
  ) %>%
  dplyr::transmute(
    site_code = trimws(
      as.character(site_code)
    )
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(site_code) %>%
  dplyr::pull(site_code)

if (length(all_site_levels) == 0) {
  stop("No valid site codes were found in the input data.")
}

cat(
  "Number of sites in the original dataset: ",
  length(all_site_levels),
  "\n",
  sep = ""
)

cat(
  "Site codes:\n",
  paste(all_site_levels, collapse = ", "),
  "\n",
  sep = ""
)


# ============================================================
# Convert variables
# ============================================================

df <- df_raw %>%
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
    
    year = suppressWarnings(
      as.integer(
        as.character(year)
      )
    ),
    
    dif.re = suppressWarnings(
      as.numeric(
        as.character(dif.re)
      )
    ),
    
    Log_Biomass_z = suppressWarnings(
      as.numeric(
        as.character(Log_Biomass_z)
      )
    ),
    
    trt = trimws(
      as.character(trt)
    )
  )


# ============================================================
# Remove the NPK pretreatment year when year_trt is available
# ============================================================

if ("year_trt" %in% names(df)) {
  
  year_trt_original <- df$year_trt
  
  df <- df %>%
    dplyr::mutate(
      year_trt = suppressWarnings(
        as.numeric(
          as.character(year_trt)
        )
      )
    )
  
  invalid_year_trt <- (
    !is.na(year_trt_original) &
      trimws(as.character(year_trt_original)) != "" &
      is.na(df$year_trt)
  )
  
  if (any(invalid_year_trt)) {
    
    invalid_values <- unique(
      year_trt_original[
        invalid_year_trt
      ]
    )
    
    warning(
      "Some year_trt values could not be converted to numeric: ",
      paste(invalid_values, collapse = ", "),
      call. = FALSE
    )
  }
  
  n_npk_year0 <- df %>%
    dplyr::filter(
      trt == "NPK",
      !is.na(year_trt),
      year_trt == 0
    ) %>%
    nrow()
  
  cat(
    "Number of NPK observations with year_trt == 0 removed: ",
    n_npk_year0,
    "\n",
    sep = ""
  )
  
  df <- df %>%
    dplyr::filter(
      !(
        trt == "NPK" &
          !is.na(year_trt) &
          year_trt == 0
      )
    )
  
} else {
  
  cat(
    paste0(
      "The year_trt column is not present. ",
      "No additional NPK pretreatment filtering was performed.\n"
    )
  )
}


# ============================================================
# Retain valid plot-level observations
# ============================================================

df <- df %>%
  dplyr::filter(
    trt %in% c(
      "Control",
      "NPK"
    ),
    
    !is.na(site_code),
    site_code != "",
    
    !is.na(block),
    block != "",
    
    !is.na(plot),
    plot != "",
    
    !is.na(year),
    
    !is.na(dif.re),
    is.finite(dif.re),
    
    !is.na(Log_Biomass_z),
    is.finite(Log_Biomass_z)
  ) %>%
  dplyr::mutate(
    site_code = factor(
      site_code,
      levels = all_site_levels
    ),
    
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    ),
    
    # Plot identifiers are unique across sites, blocks,
    # plots, and treatments.
    plot_id = interaction(
      site_code,
      block,
      plot,
      trt,
      drop = TRUE
    )
  )


# ============================================================
# Confirm removal of NPK year zero
# ============================================================

if ("year_trt" %in% names(df)) {
  
  remaining_npk_year0 <- df %>%
    dplyr::filter(
      trt == "NPK",
      !is.na(year_trt),
      year_trt == 0
    )
  
  if (nrow(remaining_npk_year0) > 0) {
    stop(
      "NPK observations with year_trt == 0 remain after filtering."
    )
  }
  
  cat(
    "Remaining NPK observations with year_trt == 0: ",
    nrow(remaining_npk_year0),
    "\n",
    sep = ""
  )
}


# ============================================================
# Report the retained plot-level data
# ============================================================

cat(
  "Number of valid plot-level observations retained: ",
  nrow(df),
  "\n",
  sep = ""
)

cat(
  "Number of sites with at least one valid observation: ",
  dplyr::n_distinct(df$site_code),
  "\n",
  sep = ""
)

cat(
  "Number of unique plots: ",
  dplyr::n_distinct(df$plot_id),
  "\n",
  sep = ""
)


# ============================================================
# Check duplicated plot-year-treatment records
# ============================================================

duplicate_plot_records <- df %>%
  dplyr::count(
    site_code,
    block,
    plot,
    year,
    trt,
    name = "n_records"
  ) %>%
  dplyr::filter(
    n_records > 1
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_records)
  )

cat(
  "Number of duplicated plot-year-treatment combinations: ",
  nrow(duplicate_plot_records),
  "\n",
  sep = ""
)

if (nrow(duplicate_plot_records) > 0) {
  
  print(
    duplicate_plot_records,
    n = 50
  )
  
  warning(
    paste0(
      "Some plot-year-treatment combinations contain more than ",
      "one record. All valid records will be displayed and used ",
      "in the GAM unless these duplicates are resolved."
    ),
    call. = FALSE
  )
}


# ============================================================
# Check precipitation values within site and year
# ============================================================

# Annual precipitation should be identical among all plots and
# treatments from the same site and year.

precipitation_check <- df %>%
  dplyr::group_by(
    site_code,
    year
  ) %>%
  dplyr::summarise(
    n_precipitation_values =
      dplyr::n_distinct(dif.re),
    
    precipitation_min = min(
      dif.re,
      na.rm = TRUE
    ),
    
    precipitation_max = max(
      dif.re,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_precipitation_values > 1
  )

if (nrow(precipitation_check) > 0) {
  
  print(
    precipitation_check,
    n = Inf
  )
  
  stop(
    paste0(
      "Some site-year combinations contain more than one ",
      "dif.re value. Resolve these inconsistencies before ",
      "fitting the site-specific GAMs."
    )
  )
}

cat(
  "Number of site-year combinations with inconsistent ",
  "precipitation values: ",
  nrow(precipitation_check),
  "\n",
  sep = ""
)


# ============================================================
# Check plot-level replication
# ============================================================

replicate_summary <- df %>%
  dplyr::group_by(
    site_code,
    year,
    trt
  ) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    
    n_plots = dplyr::n_distinct(
      interaction(
        block,
        plot,
        drop = TRUE
      )
    ),
    
    n_distinct_biomass_values =
      dplyr::n_distinct(
        Log_Biomass_z
      ),
    
    .groups = "drop"
  )





# ============================================================
# Create all site-treatment combinations
# ============================================================

complete_site_treatment <- tidyr::expand_grid(
  site_code = all_site_levels,
  trt = c(
    "Control",
    "NPK"
  )
)


# ============================================================
# Summarize plot-level data availability
# ============================================================

site_data_summary <- df %>%
  dplyr::mutate(
    site_code = as.character(site_code),
    trt = as.character(trt)
  ) %>%
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  dplyr::summarise(
    n_observations = dplyr::n(),
    
    n_years =
      dplyr::n_distinct(year),
    
    n_precipitation_values =
      dplyr::n_distinct(dif.re),
    
    n_plots =
      dplyr::n_distinct(plot_id),
    
    precipitation_min = min(
      dif.re,
      na.rm = TRUE
    ),
    
    precipitation_max = max(
      dif.re,
      na.rm = TRUE
    ),
    
    biomass_min = min(
      Log_Biomass_z,
      na.rm = TRUE
    ),
    
    biomass_max = max(
      Log_Biomass_z,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::right_join(
    complete_site_treatment,
    by = c(
      "site_code",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    n_observations = tidyr::replace_na(
      n_observations,
      0L
    ),
    
    n_years = tidyr::replace_na(
      n_years,
      0L
    ),
    
    n_precipitation_values =
      tidyr::replace_na(
        n_precipitation_values,
        0L
      ),
    
    n_plots = tidyr::replace_na(
      n_plots,
      0L
    ),
    
    site_code = factor(
      site_code,
      levels = all_site_levels
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
    trt,
    site_code
  )

print(
  site_data_summary,
  n = Inf
)

readr::write_csv(
  site_data_summary,
  file.path(
    dir.data,
    "site_treatment_plot_level_data_summary.csv"
  )
)


# ============================================================
# Identify groups for which a GAM can be estimated
# ============================================================

# These conditions are used only to determine whether a GAM
# curve can be estimated. No site or observation is removed
# from the plotted points or from the facet layout.

gam_eligible_groups <- site_data_summary %>%
  dplyr::filter(
    n_observations >= 4,
    n_precipitation_values >= 3
  ) %>%
  dplyr::transmute(
    site_code = as.character(site_code),
    trt = as.character(trt)
  )


# ============================================================
# Define an empty GAM-prediction table
# ============================================================

empty_gam_prediction <- function() {
  
  tibble::tibble(
    dif.re = double(),
    predicted = double(),
    prediction_se = double(),
    lower_CL = double(),
    upper_CL = double(),
    effective_df = double(),
    gam_p_value = double(),
    n_observations = integer(),
    n_years = integer(),
    n_plots = integer(),
    random_plot_effect = logical()
  )
}


# ============================================================
# Fit one site-treatment GAM
# ============================================================

fit_site_gam <- function(data) {
  
  data <- data %>%
    dplyr::filter(
      !is.na(Log_Biomass_z),
      is.finite(Log_Biomass_z),
      
      !is.na(dif.re),
      is.finite(dif.re),
      
      !is.na(year),
      !is.na(plot_id)
    ) %>%
    dplyr::mutate(
      plot_id = droplevels(
        factor(plot_id)
      )
    )
  
  n_observations <- nrow(data)
  
  n_years <- dplyr::n_distinct(
    data$year
  )
  
  n_unique_x <- dplyr::n_distinct(
    data$dif.re
  )
  
  n_plots <- dplyr::n_distinct(
    data$plot_id
  )
  
  # A smooth with k = 3 requires sufficient observations and
  # at least three distinct precipitation values.
  
  if (
    n_observations < 4 ||
    n_unique_x < 3
  ) {
    return(
      empty_gam_prediction()
    )
  }
  
  # Include a plot-level random intercept when at least two
  # plots are available.
  
  if (n_plots >= 2) {
    
    gam_formula <- stats::as.formula(
      paste0(
        "Log_Biomass_z ~ ",
        "s(dif.re, k = 3, bs = 'cs') + ",
        "s(plot_id, bs = 're')"
      )
    )
    
    random_plot_effect <- TRUE
    
  } else {
    
    gam_formula <- stats::as.formula(
      paste0(
        "Log_Biomass_z ~ ",
        "s(dif.re, k = 3, bs = 'cs')"
      )
    )
    
    random_plot_effect <- FALSE
  }
  
  gam_fit <- tryCatch(
    mgcv::gam(
      formula = gam_formula,
      data = data,
      method = "REML"
    ),
    error = function(e) {
      
      message(
        "GAM fitting failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  if (is.null(gam_fit)) {
    return(
      empty_gam_prediction()
    )
  }
  
  # Extract the precipitation smooth information.
  
  gam_summary <- summary(
    gam_fit
  )
  
  smooth_table <- gam_summary$s.table
  
  precipitation_row <- grep(
    "^s\\(dif\\.re\\)",
    rownames(smooth_table)
  )
  
  if (length(precipitation_row) == 1) {
    
    effective_df <- as.numeric(
      smooth_table[
        precipitation_row,
        "edf"
      ]
    )
    
    gam_p_value <- as.numeric(
      smooth_table[
        precipitation_row,
        "p-value"
      ]
    )
    
  } else {
    
    effective_df <- NA_real_
    gam_p_value <- NA_real_
  }
  
  # Generate predictions across the observed precipitation
  # range for this site and treatment.
  
  x_sequence <- seq(
    min(
      data$dif.re,
      na.rm = TRUE
    ),
    
    max(
      data$dif.re,
      na.rm = TRUE
    ),
    
    length.out = 100
  )
  
  if (random_plot_effect) {
    
    first_plot_level <- levels(
      data$plot_id
    )[1]
    
    prediction_data <- data.frame(
      dif.re = x_sequence,
      
      plot_id = factor(
        first_plot_level,
        levels = levels(
          data$plot_id
        )
      )
    )
    
    # Exclude the plot-specific random intercept from predictions.
    # The plotted curve therefore represents the site-treatment
    # precipitation response rather than one particular plot.
    
    prediction <- tryCatch(
      stats::predict(
        gam_fit,
        newdata = prediction_data,
        type = "response",
        se.fit = TRUE,
        exclude = "s(plot_id)"
      ),
      error = function(e) {
        
        message(
          "GAM prediction failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
  } else {
    
    prediction_data <- data.frame(
      dif.re = x_sequence
    )
    
    prediction <- tryCatch(
      stats::predict(
        gam_fit,
        newdata = prediction_data,
        type = "response",
        se.fit = TRUE
      ),
      error = function(e) {
        
        message(
          "GAM prediction failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
  }
  
  if (is.null(prediction)) {
    return(
      empty_gam_prediction()
    )
  }
  
  tibble::tibble(
    dif.re = x_sequence,
    
    predicted = as.numeric(
      prediction$fit
    ),
    
    prediction_se = as.numeric(
      prediction$se.fit
    ),
    
    lower_CL =
      predicted -
      1.96 * prediction_se,
    
    upper_CL =
      predicted +
      1.96 * prediction_se,
    
    effective_df = effective_df,
    
    gam_p_value = gam_p_value,
    
    n_observations = n_observations,
    
    n_years = n_years,
    
    n_plots = n_plots,
    
    random_plot_effect =
      random_plot_effect
  )
}


# ============================================================
# Prepare all eligible plot-level observations for GAM fitting
# ============================================================

gam_input_data <- df %>%
  dplyr::mutate(
    site_code = as.character(site_code),
    trt = as.character(trt)
  ) %>%
  dplyr::semi_join(
    gam_eligible_groups,
    by = c(
      "site_code",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    plot_id = factor(plot_id)
  )


# ============================================================
# Fit all site- and treatment-specific GAMs
# ============================================================

gam_prediction_data <- gam_input_data %>%
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  dplyr::group_modify(
    ~ fit_site_gam(.x)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    site_code = factor(
      site_code,
      levels = all_site_levels
    ),
    
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    )
  )


# ============================================================
# Identify groups with and without a fitted GAM
# ============================================================

groups_with_gam <- gam_prediction_data %>%
  dplyr::distinct(
    site_code,
    trt
  )

groups_without_gam <- site_data_summary %>%
  dplyr::anti_join(
    groups_with_gam,
    by = c(
      "site_code",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    gam_status = dplyr::case_when(
      n_observations == 0 ~
        "No valid plot-level observations",
      
      n_observations < 4 ~
        "Fewer than four plot-level observations",
      
      n_precipitation_values < 3 ~
        "Fewer than three distinct precipitation values",
      
      TRUE ~
        "GAM fitting or prediction failed"
    )
  )

cat(
  "Number of site-treatment combinations with a GAM curve: ",
  nrow(groups_with_gam),
  "\n",
  sep = ""
)

cat(
  "Number of site-treatment combinations without a GAM curve: ",
  nrow(groups_without_gam),
  "\n",
  sep = ""
)

if (nrow(groups_without_gam) > 0) {
  
  print(
    groups_without_gam,
    n = Inf
  )
  
  readr::write_csv(
    groups_without_gam,
    file.path(
      dir.data,
      "site_treatment_groups_without_GAM_curve.csv"
    )
  )
}


# ============================================================
# Export GAM predictions
# ============================================================

readr::write_csv(
  gam_prediction_data,
  file.path(
    dir.data,
    "site_specific_plot_level_GAM_predictions.csv"
  )
)


# ============================================================
# Export GAM model summaries
# ============================================================

gam_model_summary <- gam_prediction_data %>%
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  dplyr::summarise(
    effective_df = dplyr::first(
      effective_df
    ),
    
    gam_p_value = dplyr::first(
      gam_p_value
    ),
    
    n_observations = dplyr::first(
      n_observations
    ),
    
    n_years = dplyr::first(
      n_years
    ),
    
    n_plots = dplyr::first(
      n_plots
    ),
    
    random_plot_effect = dplyr::first(
      random_plot_effect
    ),
    
    .groups = "drop"
  )

print(
  gam_model_summary,
  n = Inf
)

readr::write_csv(
  gam_model_summary,
  file.path(
    dir.data,
    "site_specific_plot_level_GAM_summary.csv"
  )
)


# ============================================================
# Define treatment colors and labels
# ============================================================

treatment.colors <- c(
  "Control" = "#33a02c",
  "NPK" = "#984ea3"
)

treatment.labels <- c(
  "Control" = "Control",
  "NPK" = "NPK\u03bc"
)


# ============================================================
# Prepare all plot-level observations for the figure
# ============================================================

point_plot_data <- df %>%
  dplyr::mutate(
    site_code = factor(
      as.character(site_code),
      levels = all_site_levels
    ),
    
    trt = factor(
      as.character(trt),
      levels = c(
        "Control",
        "NPK"
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(site_code),
    !is.na(trt),
    
    !is.na(dif.re),
    is.finite(dif.re),
    
    !is.na(Log_Biomass_z),
    is.finite(Log_Biomass_z)
  )


# ============================================================
# Prepare GAM prediction data for plotting
# ============================================================

gam_plot_data <- gam_prediction_data %>%
  dplyr::mutate(
    site_code = factor(
      as.character(site_code),
      levels = all_site_levels
    ),
    
    trt = factor(
      as.character(trt),
      levels = c(
        "Control",
        "NPK"
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(site_code),
    !is.na(trt),
    
    !is.na(dif.re),
    is.finite(dif.re),
    
    !is.na(predicted),
    is.finite(predicted),
    
    !is.na(lower_CL),
    is.finite(lower_CL),
    
    !is.na(upper_CL),
    is.finite(upper_CL)
  )


# ============================================================
# Report plotting sample sizes
# ============================================================

cat(
  "Total number of site panels: ",
  length(all_site_levels),
  "\n",
  sep = ""
)

cat(
  "Number of plot-level observations displayed: ",
  nrow(point_plot_data),
  "\n",
  sep = ""
)

cat(
  "Number of sites with valid plot-level observations: ",
  dplyr::n_distinct(
    point_plot_data$site_code
  ),
  "\n",
  sep = ""
)


# ============================================================
# Create the combined Control–NPK figure
# ============================================================

p_all_sites <- ggplot2::ggplot() +
  
  # Horizontal zero-reference line
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.45
  ) +
  
  # Vertical zero-reference line
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.45
  ) +
  
  # GAM-based 95% confidence intervals
  # Ribbons are drawn before observations and fitted lines.
  ggplot2::geom_ribbon(
    data = gam_plot_data,
    mapping = ggplot2::aes(
      x = dif.re,
      ymin = lower_CL,
      ymax = upper_CL,
      fill = trt,
      group = trt
    ),
    alpha = 0.13,
    color = NA
  ) +
  
  # All plot-level observations
  #
  # Horizontal jitter separates observations that share the
  # same site-year precipitation value. No vertical jitter is
  # applied. Jitter is used only for visualization.
  ggplot2::geom_point(
    data = point_plot_data,
    mapping = ggplot2::aes(
      x = dif.re,
      y = Log_Biomass_z,
      fill = trt
    ),
    shape = 21,
    color = "white",
    size = 3.5,
    stroke = 0.3,
    alpha = 0.60,
    position = ggplot2::position_jitter(
      width = 0.025,
      height = 0,
      seed = 123
    )
  ) +
  
  # Site- and treatment-specific GAM curves
  # All curves are displayed as thin solid lines.
  ggplot2::geom_line(
    data = gam_plot_data,
    mapping = ggplot2::aes(
      x = dif.re,
      y = predicted,
      color = trt,
      group = trt
    ),
    linetype = "solid",
    linewidth = 0.65
  ) +
  
  # One panel for every site
  ggplot2::facet_wrap(
    facets = ~site_code,
    scales = "free_x",
    ncol = 5,
    drop = FALSE
  ) +
  
  # Treatment colors
  ggplot2::scale_color_manual(
    values = treatment.colors,
    breaks = c(
      "Control",
      "NPK"
    ),
    labels = treatment.labels,
    name = "Treatment",
    drop = FALSE
  ) +
  
  ggplot2::scale_fill_manual(
    values = treatment.colors,
    breaks = c(
      "Control",
      "NPK"
    ),
    labels = treatment.labels,
    name = "Treatment",
    drop = FALSE
  ) +
  
  # Use one treatment legend
  ggplot2::guides(
    fill = "none",
    
    color = ggplot2::guide_legend(
      override.aes = list(
        linetype = "solid",
        linewidth = 1.2
      )
    )
  ) +
  
  ggplot2::labs(
    x = "Standardized precipitation within sites",
    y = "Standardized biomass",
    color = "Treatment",
    fill = "Treatment"
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
  ) +
  
  ggplot2::theme(
    legend.position = "bottom",
    
    legend.title =
      ggplot2::element_blank(),
    
    legend.text =
      ggplot2::element_text(
        size = 11,
        color = "black"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    panel.grid.major =
      ggplot2::element_line(
        color = "gray92",
        linewidth = 0.3
      ),
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        size = 9,
        face = "plain",
        color = "black"
      ),
    
    panel.border =
      ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
    
    axis.title =
      ggplot2::element_text(
        size = 12,
        face = "plain"
      ),
    
    axis.text =
      ggplot2::element_text(
        size = 9,
        color = "black"
      ),
    
    axis.ticks =
      ggplot2::element_line(
        color = "black",
        linewidth = 0.5
      ),
    
    plot.margin =
      ggplot2::margin(
        t = 5,
        r = 8,
        b = 5,
        l = 8
      )
  )


# ============================================================
# Display the figure
# ============================================================

print(p_all_sites)


# ============================================================
# Save the combined figure
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "p_Control_NPK_allSites_siteSpecific.tiff"
  ),
  plot = p_all_sites,
  width = 10,
  height = 14,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

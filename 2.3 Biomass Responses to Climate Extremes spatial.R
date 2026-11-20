

# ============================================================
# Supplementary figure:
# Spatiotemporal relationship between precipitation and biomass
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
library(emmeans)
library(broom.mixed)


# ============================================================
# Define the data directory
# ============================================================

dir.data <- ""



if (!dir.exists(dir.data)) {
  stop("The data directory does not exist: ", dir.data)
}


# ============================================================
# Read all cleaned biomass data
# ============================================================

data.file <- file.path(
  dir.data,
  "df_clean_data.csv"
)

if (!file.exists(data.file)) {
  stop("The input file does not exist: ", data.file)
}

df <- readr::read_csv(
  data.file,
  show_col_types = FALSE
)


# ============================================================
# Check the required columns
# ============================================================

required.columns <- c(
  "site_code",
  "year",
  "block",
  "plot",
  "trt",
  "prec_sum",
  "mean_prec_2000_2020",
  "vascular_live_mass_new"
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
# Prepare all biomass observations
# ============================================================

# No filtering based on climate_display or extreme events is
# applied. All sites and years with valid biomass are retained.

df_biomass <- df %>%
  dplyr::mutate(
    year = as.integer(year),
    prec_sum = as.numeric(prec_sum),
    
    mean_prec_2000_2020 = as.numeric(
      mean_prec_2000_2020
    ),
    
    vascular_live_mass_new = as.numeric(
      vascular_live_mass_new
    ),
    
    trt = factor(
      trt,
      levels = c("Control", "NPK")
    ),
    
    site_code = factor(site_code),
    block = factor(block),
    plot = factor(plot),
    
    # Transform biomass consistently across all observations
    log_biomass = log10(
      vascular_live_mass_new + 1
    ),
    
    # Create a unique plot-treatment identifier
    plot_id = interaction(
      site_code,
      block,
      plot,
      trt,
      drop = TRUE
    )
  ) %>%
  dplyr::filter(
    trt %in% c("Control", "NPK"),
    
    !is.na(site_code),
    !is.na(year),
    !is.na(block),
    !is.na(plot),
    
    !is.na(prec_sum),
    is.finite(prec_sum),
    
    !is.na(mean_prec_2000_2020),
    is.finite(mean_prec_2000_2020),
    
    !is.na(vascular_live_mass_new),
    is.finite(vascular_live_mass_new),
    vascular_live_mass_new >= 0,
    
    !is.na(log_biomass),
    is.finite(log_biomass)
  ) %>%
  droplevels()


# ============================================================
# Check the retained observations
# ============================================================

cat(
  "Number of observations:",
  nrow(df_biomass),
  "\n"
)

cat(
  "Number of sites:",
  dplyr::n_distinct(df_biomass$site_code),
  "\n"
)

cat(
  "Number of plots:",
  dplyr::n_distinct(df_biomass$plot_id),
  "\n"
)

sample_summary <- df_biomass %>%
  dplyr::group_by(trt) %>%
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_sites = dplyr::n_distinct(site_code),
    n_plots = dplyr::n_distinct(plot_id),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    .groups = "drop"
  )

print(sample_summary)

readr::write_csv(
  sample_summary,
  file.path(
    dir.data,
    "biomass_precipitation_sample_sizes.csv"
  )
)


# ============================================================
# Calculate site-year mean biomass
# ============================================================

site_year_data <- df_biomass %>%
  dplyr::group_by(
    site_code,
    year,
    trt
  ) %>%
  dplyr::summarise(
    prec_sum = dplyr::first(prec_sum),
    
    mean_log_biomass = mean(
      log_biomass,
      na.rm = TRUE
    ),
    
    n_observations = dplyr::n(),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(prec_sum),
    is.finite(mean_log_biomass)
  )


# ============================================================
# Calculate mean positions for each site and treatment
# ============================================================

site_mean_data <- df_biomass %>%
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  dplyr::summarise(
    mean_prec_2000_2020 = dplyr::first(
      mean_prec_2000_2020
    ),
    
    mean_biomass = mean(
      vascular_live_mass_new,
      na.rm = TRUE
    ),
    
    mean_log_biomass = log10(
      mean_biomass + 1
    ),
    
    n_years = dplyr::n_distinct(year),
    n_observations = dplyr::n(),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(mean_prec_2000_2020),
    is.finite(mean_biomass),
    is.finite(mean_log_biomass)
  )


# ============================================================
# Select groups suitable for site-specific temporal lines
# ============================================================

# At least two years with two different precipitation values
# are needed to fit a site-specific line.

valid_line_groups <- site_year_data %>%
  dplyr::group_by(
    site_code,
    trt
  ) %>%
  dplyr::summarise(
    n_years = dplyr::n_distinct(year),
    n_precipitation_values = dplyr::n_distinct(prec_sum),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_years >= 2,
    n_precipitation_values >= 2
  )

site_year_lines <- site_year_data %>%
  dplyr::semi_join(
    valid_line_groups,
    by = c("site_code", "trt")
  )


# ============================================================
# Center and scale annual growing-season precipitation
# ============================================================

# One unit of precipitation_100 represents 100 mm.

precipitation_mean <- mean(
  df_biomass$prec_sum,
  na.rm = TRUE
)

df_biomass <- df_biomass %>%
  dplyr::mutate(
    precipitation_100 = (
      prec_sum - precipitation_mean
    ) / 100
  )

cat(
  "Mean annual precipitation used for centering:",
  round(precipitation_mean, 2),
  "mm\n"
)


# ============================================================
# Fit the mixed-effects model using all observations
# ============================================================

biomass_model <- lmerTest::lmer(
  log_biomass ~ precipitation_100 * trt +
    (1 | site_code) +
    (1 | plot_id),
  data = df_biomass,
  REML = TRUE,
  na.action = na.omit
)

cat("\nMixed-effects model summary:\n")
print(summary(biomass_model))

cat("\nType III tests:\n")
print(
  anova(
    biomass_model,
    type = 3
  )
)


# ============================================================
# Export fixed-effect results
# ============================================================

fixed_effect_results <- broom.mixed::tidy(
  biomass_model,
  effects = "fixed",
  conf.int = TRUE
)

print(fixed_effect_results)


# ============================================================
# Extract precipitation slopes for Control and NPK
# ============================================================

precipitation_trends <- emmeans::emtrends(
  biomass_model,
  specs = ~trt,
  var = "precipitation_100"
)

slope_results <- as.data.frame(
  summary(
    precipitation_trends,
    infer = c(TRUE, TRUE)
  )
)

# Identify the trend column generated by emtrends
slope.column <- grep(
  "\\.trend$",
  names(slope_results),
  value = TRUE
)

if (length(slope.column) != 1) {
  stop("The precipitation-slope column could not be identified.")
}

slope_results <- slope_results %>%
  dplyr::rename(
    slope_per_100mm = dplyr::all_of(
      slope.column
    )
  ) %>%
  dplyr::mutate(
    trt = as.character(trt),
    
    slope_per_1mm = slope_per_100mm / 100,
    
    line_significance = dplyr::case_when(
      is.na(p.value) ~ "P unavailable",
      p.value < 0.05 ~ "P < 0.05",
      TRUE ~ "P >= 0.05"
    ),
    
    slope_label = dplyr::case_when(
      is.na(p.value) ~
        paste0(
          "Slope = ",
          sprintf("%.3f", slope_per_100mm),
          ", P = NA"
        ),
      
      p.value < 0.001 ~
        paste0(
          "Slope = ",
          sprintf("%.3f", slope_per_100mm),
          ", P < 0.001"
        ),
      
      TRUE ~
        paste0(
          "Slope = ",
          sprintf("%.3f", slope_per_100mm),
          ", P = ",
          sprintf("%.3f", p.value)
        )
    )
  )

print(slope_results)


# ============================================================
# Compare the Control and NPK slopes
# ============================================================

slope_contrast <- as.data.frame(
  summary(
    pairs(precipitation_trends),
    infer = c(TRUE, TRUE)
  )
)

print(slope_contrast)


# ============================================================
# Generate fixed-effect model predictions
# ============================================================

precipitation_sequence <- seq(
  min(df_biomass$prec_sum, na.rm = TRUE),
  max(df_biomass$prec_sum, na.rm = TRUE),
  length.out = 200
)

prediction_data <- expand.grid(
  prec_sum = precipitation_sequence,
  trt = c("Control", "NPK"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    precipitation_100 = (
      prec_sum - precipitation_mean
    ) / 100,
    
    trt = factor(
      trt,
      levels = c("Control", "NPK")
    )
  )

# Construct the fixed-effect design matrix
design_matrix <- stats::model.matrix(
  ~ precipitation_100 * trt,
  data = prediction_data
)

fixed_coefficients <- lme4::fixef(
  biomass_model
)

covariance_matrix <- as.matrix(
  stats::vcov(biomass_model)
)

# Match the design matrix to the model coefficients
design_matrix <- design_matrix[
  ,
  names(fixed_coefficients),
  drop = FALSE
]

# Population-level fixed-effect predictions
prediction_data$predicted_log_biomass <- as.numeric(
  design_matrix %*% fixed_coefficients
)

# Standard errors of the predictions
prediction_data$prediction_se <- sqrt(
  rowSums(
    (
      design_matrix %*%
        covariance_matrix
    ) *
      design_matrix
  )
)

# Approximate 95% confidence intervals
prediction_data <- prediction_data %>%
  dplyr::mutate(
    lower_CL = predicted_log_biomass -
      1.96 * prediction_se,
    
    upper_CL = predicted_log_biomass +
      1.96 * prediction_se,
    
    trt = as.character(trt)
  ) %>%
  dplyr::left_join(
    slope_results %>%
      dplyr::select(
        trt,
        slope_per_100mm,
        p.value,
        line_significance,
        slope_label
      ),
    by = "trt"
  ) %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c("Control", "NPK")
    )
  )


# ============================================================
# Prepare slope labels
# ============================================================

# Data ranges
x_range <- range(df_biomass$prec_sum, na.rm = TRUE)
y_range <- range(df_biomass$log_biomass, na.rm = TRUE)

label_data <- slope_results %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c("Control", "NPK")
    ),
    
    # Show treatment name and P value in scientific notation
    slope_label = dplyr::case_when(
      is.na(p.value) ~ paste0(
        dplyr::if_else(trt == "NPK", "NPK\u03bc", "Control"),
        ": P = NA"
      ),
      
      TRUE ~ paste0(
        dplyr::if_else(trt == "NPK", "NPK\u03bc", "Control"),
        ": P = ",
        formatC(
          p.value,
          format = "e",
          digits = 2
        )
      )
    ),
    
    # Position the labels near the lower-right corner
    label_x = x_range[2] - 0.03 * diff(x_range),
    
    label_y = dplyr::case_when(
      trt == "Control" ~
        y_range[1] + 0.06 * diff(y_range),
      
      trt == "NPK" ~
        y_range[1] + 0.14 * diff(y_range),
      
      TRUE ~ y_range[1]
    )
  )


# ============================================================
# Define colors and labels
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
# Plot the spatiotemporal relationship
# ============================================================

Bio_Pre <- ggplot2::ggplot() +
  
  # Mixed-model 95% confidence intervals
  # Draw the confidence intervals first so that they remain
  # behind all points and lines
  ggplot2::geom_ribbon(
    data = prediction_data,
    ggplot2::aes(
      x = prec_sum,
      ymin = lower_CL,
      ymax = upper_CL,
      fill = trt,
      group = trt
    ),
    alpha = 0.18,
    color = NA
  ) +
  
  # Pale site-specific temporal relationships
  ggplot2::geom_smooth(
    data = site_year_lines,
    ggplot2::aes(
      x = prec_sum,
      y = mean_log_biomass,
      color = trt,
      group = interaction(site_code, trt)
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.4,
    alpha = 0.15,
    show.legend = FALSE
  ) +
  
  # Plot-level annual observations
  ggplot2::geom_point(
    data = df_biomass,
    ggplot2::aes(
      x = prec_sum,
      y = log_biomass,
      fill = trt
    ),
    shape = 21,
    color = "white",
    size = 2.5,
    stroke = 0.25,
    alpha = 0.18
  ) +
  
  # Site-treatment mean positions
  # Use smaller circles with black borders
  ggplot2::geom_point(
    data = site_mean_data,
    ggplot2::aes(
      x = mean_prec_2000_2020,
      y = mean_log_biomass,
      fill = trt
    ),
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.65,
    alpha = 0.6
  ) +
  
  # Draw these lines after the points so they appear on top
  ggplot2::geom_line(
    data = prediction_data,
    ggplot2::aes(
      x = prec_sum,
      y = predicted_log_biomass,
      color = trt,
      linetype = line_significance,
      group = trt
    ),
    linewidth = 1.4
  ) +
  
  # Mixed-model slopes and P values
  ggplot2::geom_text(
    data = label_data,
    ggplot2::aes(
      x = label_x,
      y = label_y,
      label = slope_label,
      color = trt
    ),
    hjust = 1,
    vjust = 0,
    size = 4,
    show.legend = FALSE
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
  
  ggplot2::scale_linetype_manual(
    values = c(
      "P < 0.05" = "solid",
      "P >= 0.05" = "dashed"
    ),
    guide = "none"
  ) +
  
  ggplot2::labs(
    x = "Annual growing-season precipitation (mm)",
    y = "Biomass (log-transformed)"
  ) +
  
  ggplot2::theme_classic(
    base_size = 15
  ) +
  
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(
      size = 14
    ),
    
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
      t = 5,
      r = 10,
      b = 5,
      l = 10
    )
  )

print(Bio_Pre)


# ============================================================
# Save the figure as TIFF
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "Bio_Pre_all_sites_spatiotemporal.tiff"
  ),
  plot = Bio_Pre,
  width = 6,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)



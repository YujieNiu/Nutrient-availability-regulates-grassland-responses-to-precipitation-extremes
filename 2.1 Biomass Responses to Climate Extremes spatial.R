
# ============================================================
# Biomass responses to extreme precipitation events
# Extreme effect = ln(Ye / Yn)
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


# ============================================================
# Define the data directory
# ============================================================

dir.data <- ""




if (!dir.exists(dir.data)) {
  stop("The data directory does not exist: ", dir.data)
}


# ============================================================
# Read the data
# ============================================================

data.file <- file.path(
  dir.data,
  "full_data_with_climate.csv"
)

if (!file.exists(data.file)) {
  stop("The input file does not exist: ", data.file)
}

data <- read.csv(
  data.file,
  header = TRUE,
  stringsAsFactors = FALSE
)


# ============================================================
# Check the required variables
# ============================================================

required.columns <- c(
  "site_code",
  "year",
  "block",
  "plot",
  "trt",
  "year_trt",
  "climate_display",
  "vascular_live_mass"
)

missing.columns <- setdiff(
  required.columns,
  names(data)
)

if (length(missing.columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing.columns, collapse = ", ")
  )
}


# ============================================================
# Clean the biomass data
# ============================================================

df_clean <- data %>%
  dplyr::mutate(
    year = as.integer(year),
    year_trt = as.numeric(year_trt),
    block = as.character(block),
    plot = as.character(plot),
    trt = as.character(trt),
    climate_display = as.character(climate_display),
    vascular_live_mass = as.numeric(vascular_live_mass),
    vascular_live_mass_new = dplyr::if_else(
      trt == "NPK" & year_trt == 0,
      NA_real_,
      vascular_live_mass
    )
  ) %>%
  dplyr::filter(
    trt %in% c("Control", "NPK")
  ) %>%
  dplyr::arrange(
    site_code,
    block,
    plot,
    trt,
    year
  )

# Save the cleaned data
readr::write_csv(
  df_clean,
  file.path(dir.data, "df_clean_data.csv")
)

cat(
  "Number of sites in the cleaned data:",
  dplyr::n_distinct(df_clean$site_code),
  "\n"
)

cat(
  "Number of plots in the cleaned data:",
  df_clean %>%
    dplyr::distinct(site_code, block, plot, trt) %>%
    nrow(),
  "\n"
)

# Check the climate categories
print(
  table(
    df_clean$climate_display,
    useNA = "ifany"
  )
)


# ============================================================
# Calculate Yn: mean biomass during normal years
# ============================================================

normal_biomass <- df_clean %>%
  dplyr::filter(
    climate_display == "Normal",
    !is.na(vascular_live_mass_new),
    is.finite(vascular_live_mass_new),
    vascular_live_mass_new > 0
  ) %>%
  dplyr::group_by(
    site_code,
    block,
    plot,
    trt
  ) %>%
  dplyr::summarise(
    Yn_biom = mean(
      vascular_live_mass_new,
      na.rm = TRUE
    ),
    n_normal_years = dplyr::n_distinct(year),
    .groups = "drop"
  )


# ============================================================
# Calculate Ye: biomass during each extreme year
# ============================================================

# If more than one observation occurs for the same plot and year,
# their mean is used as Ye.

extreme_biomass <- df_clean %>%
  dplyr::filter(
    climate_display %in% c(
      "Extreme dry",
      "Extreme wet"
    )
  ) %>%
  dplyr::group_by(
    site_code,
    block,
    plot,
    trt,
    extreme_year = year,
    extreme_type = climate_display
  ) %>%
  dplyr::summarise(
    Ye = {
      x <- vascular_live_mass_new[
        !is.na(vascular_live_mass_new) &
          is.finite(vascular_live_mass_new)
      ]
      
      if (length(x) == 0) {
        NA_real_
      } else {
        mean(x)
      }
    },
    .groups = "drop"
  )


# ============================================================
# Calculate the extreme-event biomass response
# ============================================================

recovery_scores <- extreme_biomass %>%
  dplyr::left_join(
    normal_biomass,
    by = c(
      "site_code",
      "block",
      "plot",
      "trt"
    )
  ) %>%
  dplyr::mutate(
    # Extreme-event effect: ln(Ye / Yn)
    LRR_biom = dplyr::if_else(
      !is.na(Ye) &
        is.finite(Ye) &
        Ye > 0 &
        !is.na(Yn_biom) &
        is.finite(Yn_biom) &
        Yn_biom > 0,
      log(Ye / Yn_biom),
      NA_real_
    ),
    
    # Log-transformed biomass during the extreme year
    log_biomass = dplyr::if_else(
      !is.na(Ye) &
        is.finite(Ye) &
        Ye > 0,
      log(Ye),
      NA_real_
    ),
    
    # Keep the old output-column name for compatibility
    vascular_live_mass_new = Ye
  )


# ============================================================
# Add mean precipitation information
# ============================================================

if ("mean_prec_2000_2020" %in% names(df_clean)) {
  
  site_precipitation <- df_clean %>%
    dplyr::select(
      site_code,
      mean_prec_2000_2020
    ) %>%
    dplyr::filter(
      !is.na(mean_prec_2000_2020)
    ) %>%
    dplyr::distinct(
      site_code,
      .keep_all = TRUE
    )
  
  recovery_scores <- recovery_scores %>%
    dplyr::left_join(
      site_precipitation,
      by = "site_code"
    )
  
} else {
  
  warning(
    "The variable mean_prec_2000_2020 was not found in df_clean."
  )
}


# ============================================================
# Arrange and save the biomass-response results
# ============================================================

recovery_scores <- recovery_scores %>%
  dplyr::select(
    site_code,
    block,
    plot,
    trt,
    extreme_type,
    extreme_year,
    Ye,
    vascular_live_mass_new,
    Yn_biom,
    n_normal_years,
    LRR_biom,
    log_biomass,
    dplyr::any_of("mean_prec_2000_2020")
  ) %>%
  dplyr::arrange(
    site_code,
    block,
    plot,
    trt,
    extreme_year
  )

readr::write_csv(
  recovery_scores,
  file.path(
    dir.data,
    "Biomass response result.csv"
  )
)


# ============================================================
# Check the biomass-response results
# ============================================================

cat(
  "Number of biomass-response observations:",
  nrow(recovery_scores),
  "\n"
)

cat(
  "Number of valid biomass LRR observations:",
  sum(!is.na(recovery_scores$LRR_biom)),
  "\n"
)



# ============================================================
# Summarize biomass responses
# ============================================================

plot_summary <- recovery_scores %>%
  dplyr::filter(
    trt %in% c("Control", "NPK"),
    extreme_type %in% c(
      "Extreme dry",
      "Extreme wet"
    ),
    !is.na(LRR_biom),
    is.finite(LRR_biom)
  ) %>%
  dplyr::group_by(
    trt,
    extreme_type
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_sites = dplyr::n_distinct(site_code),
    mean = mean(LRR_biom, na.rm = TRUE),
    sd = stats::sd(LRR_biom, na.rm = TRUE),
    se = sd / sqrt(n),
    
    # The confidence interval is NA when n < 2
    tval = dplyr::if_else(
      n > 1,
      stats::qt(0.975, df = n - 1),
      NA_real_
    ),
    
    ci95 = tval * se,
    ymin = mean - ci95,
    ymax = mean + ci95,
    
    # Convert the mean LRR to percentage change
    perc_change = (exp(mean) - 1) * 100,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    trt = factor(
      trt,
      levels = c("Control", "NPK")
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
    trt,
    extreme_type
  )

print(plot_summary)


# ============================================================
# Define extreme-event colors
# ============================================================

# ============================================================
# Plot biomass responses
# ============================================================

# Colors for extreme events
extreme_colors <- c(
  "Extreme dry" = "#ff9f00",
  "Extreme wet" = "#1f78b4"
)


# ============================================================
# Define a common Y-axis range
# ============================================================

# Include zero in the Y-axis range
y_min <- min(
  c(plot_summary$ymin, 0),
  na.rm = TRUE
)

y_max <- max(
  c(plot_summary$ymax, 0),
  na.rm = TRUE
)

y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- 1
}

common_y_limits <- c(
  y_min - 0.08 * y_range,
  y_max + 0.08 * y_range
)


# ============================================================
# Plotting function
# ============================================================

plot_biomass_by_treatment <- function(
    data,
    treatment,
    show_y_title = TRUE
) {
  
  plot_data <- data %>%
    dplyr::filter(
      as.character(trt) == treatment
    ) %>%
    dplyr::mutate(
      # Use the existing perc_change column
      percentage_label = dplyr::case_when(
        is.na(perc_change) ~ "",
        perc_change > 0 ~ paste0(
          "+",
          sprintf("%.1f", perc_change),
          "%"
        ),
        TRUE ~ paste0(
          sprintf("%.1f", perc_change),
          "%"
        )
      ),
      
      # Place the percentage labels toward the center
      label_hjust = dplyr::case_when(
        as.character(extreme_type) == "Extreme dry" ~ -0.35,
        as.character(extreme_type) == "Extreme wet" ~ 1.35,
        TRUE ~ 0.5
      )
    )
  
  if (nrow(plot_data) == 0) {
    stop("No data were found for treatment: ", treatment)
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = extreme_type,
      y = mean,
      color = extreme_type,
      fill = extreme_type
    )
  ) +
    
    # Zero-reference line
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray60",
      linewidth = 0.6
    ) +
    
    # 95% confidence intervals
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = ymin,
        ymax = ymax
      ),
      width = 0,
      linewidth = 1.1
    ) +
    
    # Mean LRR
    ggplot2::geom_point(
      size = 4.5,
      shape = 21,
      stroke = 0.8
    ) +
    
    # Percentage change beside each mean point
    ggplot2::geom_text(
      ggplot2::aes(
        label = percentage_label,
        hjust = label_hjust
      ),
      color = "black",
      size = 3.8,
      angle = 90,
      hjust = 0.5,
      vjust = 2,
      show.legend = FALSE
    ) +
    
    # Extreme-event colors
    ggplot2::scale_color_manual(
      values = extreme_colors
    ) +
    
    ggplot2::scale_fill_manual(
      values = extreme_colors
    ) +
    
    # Identical Y-axis limits for both plots
    ggplot2::coord_cartesian(
      ylim = common_y_limits,
      clip = "off"
    ) +
    
    # No panel title and no X-axis title
    ggplot2::labs(
      x = NULL,
      y = if (show_y_title) {
        "Biomass response (LRR)"
      } else {
        NULL
      }
    ) +
    
    ggplot2::theme_classic(
      base_size = 14
    ) +
    
    ggplot2::theme(
      # Remove the dry/wet legend
      legend.position = "none",
      
      axis.title.y = ggplot2::element_text(
        size = 14,
        face = "plain",
        margin = ggplot2::margin(r = 8)
      ),
      
      # Remove Extreme dry and Extreme wet text
      axis.text.x = ggplot2::element_blank(),
      
      axis.text.y = ggplot2::element_text(
        size = 12,
        color = "black"
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.6
      ),
      
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.6
      ),
      
      panel.grid = ggplot2::element_blank(),
      
      plot.margin = ggplot2::margin(
        t = 5,
        r = 18,
        b = 5,
        l = 18
      ),
      
      panel.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      
      legend.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      
      legend.box.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      )
    )
}


# ============================================================
# Create the Control plot
# ============================================================

(p_biom_control <- plot_biomass_by_treatment(
  data = plot_summary,
  treatment = "Control",
  show_y_title = TRUE
))


# ============================================================
# Create the NPK plot
# ============================================================

(p_biom_npk <- plot_biomass_by_treatment(
  data = plot_summary,
  treatment = "NPK",
  show_y_title = TRUE
))


# ============================================================
# Save the plots separately as TIFF files
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "Biomass_response_Control.tiff"
  ),
  plot = p_biom_control,
  width = 2.5,
  height = 3.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "transparent"
)

ggplot2::ggsave(
  filename = file.path(
    dir.data,
    "Biomass_response_NPKmu.tiff"
  ),
  plot = p_biom_npk,
  width = 2.5,
  height = 3.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "transparent"
)



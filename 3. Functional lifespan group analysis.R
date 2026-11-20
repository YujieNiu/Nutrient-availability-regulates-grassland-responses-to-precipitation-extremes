


#======================================================================
# Forest plots: LRR of cover under extreme dry / extreme wet years
# along the site precipitation gradient (Control vs NPK)
#
# Overall functional groups and lifespan
# Functional group x lifespan combinations
#
# Points: FILLED circle = significant slope
#         HOLLOW circle = non-significant slope
#======================================================================


#======================================================================
# Required packages
#======================================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(openxlsx)


#======================================================================
# Setup and data import
#======================================================================

dir.data <- ""



dir.data <- "S:/Science advances/rui_New"

if (!dir.exists(dir.data)) {
  stop(
    "The data directory does not exist: ",
    dir.data
  )
}

path.df <- file.path(
  dir.data,
  "df_clean_data_with_FG_Lifespan_cross_total.csv"
)

if (!file.exists(path.df)) {
  stop(
    "The data file does not exist: ",
    path.df
  )
}

df <- readr::read_csv(
  path.df,
  show_col_types = FALSE
)


#======================================================================
# Define variables and check the dataset
#======================================================================

keys <- c(
  "site_code",
  "block",
  "plot",
  "year",
  "trt"
)

cover_variables <- c(
  "Forb",
  "Graminoid",
  "Annual",
  "Perennial"
)

required_columns <- c(
  keys,
  "year_trt",
  "climate_display",
  "mean_prec_2000_2020",
  cover_variables
)

missing_columns <- setdiff(
  required_columns,
  names(df)
)

if (length(missing_columns) > 0) {
  stop(
    "The data file is missing the following columns: ",
    paste(
      missing_columns,
      collapse = ", "
    ),
    "\nAvailable columns: ",
    paste(
      names(df),
      collapse = ", "
    )
  )
}


#======================================================================
# Convert and check year_trt
#======================================================================

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
    trimws(
      as.character(year_trt_original)
    ) != "" &
    is.na(df$year_trt)
)

if (any(invalid_year_trt)) {
  
  invalid_values <- unique(
    year_trt_original[
      invalid_year_trt
    ]
  )
  
  stop(
    "The following year_trt values could not be converted to numeric: ",
    paste(
      invalid_values,
      collapse = ", "
    )
  )
}


#======================================================================
# Check precipitation and cover variables
#======================================================================

if (!is.numeric(df$mean_prec_2000_2020)) {
  stop(
    "The column 'mean_prec_2000_2020' must be numeric. ",
    "Its current class is: ",
    paste(
      class(df$mean_prec_2000_2020),
      collapse = ", "
    )
  )
}

non_numeric_cover <- cover_variables[
  !vapply(
    df[cover_variables],
    is.numeric,
    FUN.VALUE = logical(1)
  )
]

if (length(non_numeric_cover) > 0) {
  stop(
    "The following cover variables must be numeric: ",
    paste(
      non_numeric_cover,
      collapse = ", "
    )
  )
}

negative_cover_summary <- vapply(
  df[cover_variables],
  function(x) {
    sum(
      x < 0,
      na.rm = TRUE
    )
  },
  FUN.VALUE = numeric(1)
)



if (any(negative_cover_summary > 0)) {
  
  print(
    negative_cover_summary[
      negative_cover_summary > 0
    ]
  )
  
  stop(
    "Negative cover values were detected."
  )
}


#======================================================================
# Prepare treatment and climate variables
#======================================================================

df <- df %>%
  dplyr::mutate(
    trt = trimws(
      as.character(trt)
    ),
    
    climate_display = trimws(
      as.character(climate_display)
    )
  )

unexpected_treatments <- setdiff(
  unique(
    stats::na.omit(df$trt)
  ),
  c(
    "Control",
    "NPK"
  )
)

if (length(unexpected_treatments) > 0) {
  warning(
    "The following treatment values will be excluded: ",
    paste(
      unexpected_treatments,
      collapse = ", "
    ),
    call. = FALSE
  )
}


#======================================================================
# Check NPK observations without year_trt
#======================================================================

missing_npk_year_trt <- df %>%
  dplyr::filter(
    trt == "NPK",
    is.na(year_trt)
  )

if (nrow(missing_npk_year_trt) > 0) {
  
  print(
    as.data.frame(
      missing_npk_year_trt %>%
        dplyr::select(
          dplyr::all_of(keys),
          year_trt
        )
    ),
    row.names = FALSE
  )
  
  stop(
    "Some NPK observations have missing year_trt values. ",
    "The pretreatment year cannot be identified reliably."
  )
}


#======================================================================
# Remove NPK treatment-year-zero observations and excluded sites
#======================================================================

n_npk_year_zero <- sum(
  df$trt == "NPK" &
    !is.na(df$year_trt) &
    df$year_trt == 0,
  na.rm = TRUE
)

cat(
  "Number of NPK observations with year_trt == 0 removed: ",
  n_npk_year_zero,
  "\n",
  sep = ""
)

df <- df %>%
  dplyr::filter(
    trt %in% c(
      "Control",
      "NPK"
    ),
    
    !(
      trt == "NPK" &
        !is.na(year_trt) &
        year_trt == 0
    ),
    
    !site_code %in% c(
      "cowi.ca",
      "msla_3.us"
    )
  ) %>%
  dplyr::mutate(
    site_code = factor(site_code),
    block = factor(block),
    plot = factor(plot),
    
    trt = factor(
      trt,
      levels = c(
        "Control",
        "NPK"
      )
    )
  ) %>%
  droplevels()


#======================================================================
# Display basic dataset information
#======================================================================

cat(
  "Number of observations after filtering: ",
  nrow(df),
  "\n",
  sep = ""
)

cat(
  "Number of sites after filtering: ",
  dplyr::n_distinct(df$site_code),
  "\n",
  sep = ""
)

cat(
  "Treatments: ",
  paste(
    levels(df$trt),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Climate categories: ",
  paste(
    sort(
      unique(df$climate_display)
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Number of missing mean_prec_2000_2020 values: ",
  sum(
    is.na(df$mean_prec_2000_2020)
  ),
  "\n",
  sep = ""
)

if (all(is.na(df$mean_prec_2000_2020))) {
  stop(
    "All mean_prec_2000_2020 values are missing."
  )
}

cat(
  "Range of mean_prec_2000_2020: ",
  paste(
    range(
      df$mean_prec_2000_2020,
      na.rm = TRUE
    ),
    collapse = " to "
  ),
  "\n",
  sep = ""
)


#======================================================================
# Check duplicated observations
#======================================================================

duplicate_keys <- df %>%
  dplyr::count(
    dplyr::across(
      dplyr::all_of(keys)
    ),
    name = "n_rows"
  ) %>%
  dplyr::filter(
    n_rows > 1
  )

cat(
  "Number of duplicated site-block-plot-year-treatment records: ",
  nrow(duplicate_keys),
  "\n",
  sep = ""
)

if (nrow(duplicate_keys) > 0) {
  
  print(
    as.data.frame(duplicate_keys),
    row.names = FALSE
  )
  
  stop(
    "Duplicated site-block-plot-year-treatment records were detected."
  )
}


#======================================================================
# Check that precipitation is constant within sites
#======================================================================

site_precipitation_check <- df %>%
  dplyr::filter(
    !is.na(mean_prec_2000_2020),
    is.finite(mean_prec_2000_2020)
  ) %>%
  dplyr::group_by(site_code) %>%
  dplyr::summarise(
    n_precipitation_values = dplyr::n_distinct(
      mean_prec_2000_2020
    ),
    
    precipitation_min = min(
      mean_prec_2000_2020
    ),
    
    precipitation_max = max(
      mean_prec_2000_2020
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_precipitation_values > 1
  )

if (nrow(site_precipitation_check) > 0) {
  
  print(
    as.data.frame(site_precipitation_check),
    row.names = FALSE
  )
  
  stop(
    "Some sites have more than one mean_prec_2000_2020 value."
  )
}

if (!all(c("Control", "NPK") %in% levels(df$trt))) {
  stop(
    "Both Control and NPK treatments are required after filtering."
  )
}


#======================================================================
# Convert functional-group cover to long format
#======================================================================

long1 <- df %>%
  dplyr::select(
    dplyr::all_of(
      c(
        keys,
        "climate_display",
        "mean_prec_2000_2020"
      )
    ),
    dplyr::all_of(
      cover_variables
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      cover_variables
    ),
    names_to = "Group",
    values_to = "Cover"
  )


#======================================================================
# Calculate mean cover across normal years
#======================================================================

clean1 <- long1 %>%
  dplyr::group_by(
    site_code,
    block,
    plot,
    trt,
    Group
  ) %>%
  dplyr::mutate(
    n_normal_years = sum(
      climate_display == "Normal" &
        !is.na(Cover),
      na.rm = TRUE
    ),
    
    Yn_cover = if (
      sum(
        climate_display == "Normal" &
        !is.na(Cover),
        na.rm = TRUE
      ) > 0
    ) {
      mean(
        Cover[
          climate_display == "Normal" &
            !is.na(Cover)
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    n_normal_years > 0,
    !is.na(Yn_cover),
    is.finite(Yn_cover)
  )

cat(
  "\nDistribution of the number of normal years used:\n"
)

normal_year_summary <- clean1 %>%
  dplyr::distinct(
    site_code,
    block,
    plot,
    trt,
    Group,
    n_normal_years
  ) %>%
  dplyr::count(
    n_normal_years,
    name = "n_plot_groups"
  )

print(
  as.data.frame(normal_year_summary),
  row.names = FALSE
)


#======================================================================
# Calculate group-specific adaptive epsilon values
#======================================================================

eps1 <- clean1 %>%
  dplyr::filter(
    !is.na(Cover),
    is.finite(Cover),
    Cover > 0
  ) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    eps = max(
      1e-8,
      0.01 * median(
        Cover,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

missing_epsilon_groups <- setdiff(
  unique(clean1$Group),
  eps1$Group
)

if (length(missing_epsilon_groups) > 0) {
  stop(
    "No positive cover values were available for: ",
    paste(
      missing_epsilon_groups,
      collapse = ", "
    )
  )
}


#======================================================================
# Calculate LRR for extreme precipitation years
#======================================================================

scores1 <- clean1 %>%
  dplyr::filter(
    climate_display %in% c(
      "Extreme dry",
      "Extreme wet"
    )
  ) %>%
  dplyr::left_join(
    eps1,
    by = "Group"
  ) %>%
  dplyr::filter(
    !is.na(Cover),
    is.finite(Cover),
    
    !is.na(eps),
    is.finite(eps),
    
    !(Cover == 0 & Yn_cover == 0)
  ) %>%
  dplyr::mutate(
    extreme_type = climate_display,
    
    LRR_cover = dplyr::case_when(
      Cover > 0 & Yn_cover > 0 ~
        log(
          Cover / Yn_cover
        ),
      
      TRUE ~
        log(
          (Cover + eps) /
            (Yn_cover + eps)
        )
    )
  ) %>%
  dplyr::filter(
    !is.na(LRR_cover),
    is.finite(LRR_cover),
    
    !is.na(mean_prec_2000_2020),
    is.finite(mean_prec_2000_2020)
  )

cat(
  "\nNumber of observations available for mixed-effects models: ",
  nrow(scores1),
  "\n",
  sep = ""
)

if (nrow(scores1) == 0) {
  stop(
    "No observations are available for mixed-effects model fitting."
  )
}


#======================================================================
# Define empty result tables
#======================================================================

empty_fixed_effects <- function() {
  
  tibble::tibble(
    Group = character(),
    extreme_type = character(),
    term = character(),
    Estimate = double(),
    SE = double(),
    DF = double(),
    CI_lower = double(),
    CI_upper = double(),
    t_value = double(),
    P_value = double(),
    n_sites = integer(),
    n_obs = integer()
  )
}

empty_slope_results <- function() {
  
  tibble::tibble(
    Group = character(),
    extreme_type = character(),
    trt = character(),
    Estimate = double(),
    SE = double(),
    DF = double(),
    CI_lower = double(),
    CI_upper = double(),
    t_value = double(),
    P_value = double(),
    n_sites = integer(),
    n_obs = integer()
  )
}


#======================================================================
# Fit one mixed-effects model and extract results
#======================================================================

fit_mixed_model <- function(
    dat,
    min_n_per_treatment = 5,
    min_sites_per_treatment = 3
) {
  
  current_group <- unique(
    as.character(dat$Group)
  )
  
  current_extreme <- unique(
    as.character(dat$extreme_type)
  )
  
  if (length(current_group) != 1) {
    stop(
      "Each model dataset must contain exactly one functional group."
    )
  }
  
  if (length(current_extreme) != 1) {
    stop(
      "Each model dataset must contain exactly one extreme-event type."
    )
  }
  
  
  #--------------------------------------------------------------------
  # Prepare the model dataset
  #--------------------------------------------------------------------
  
  model_dat <- dat %>%
    dplyr::filter(
      !is.na(LRR_cover),
      is.finite(LRR_cover),
      
      !is.na(mean_prec_2000_2020),
      is.finite(mean_prec_2000_2020),
      
      !is.na(site_code),
      !is.na(block),
      !is.na(plot),
      !is.na(trt)
    ) %>%
    dplyr::mutate(
      prec_100 = mean_prec_2000_2020 / 100,
      
      site_code = factor(site_code),
      block = factor(block),
      plot = factor(plot),
      
      trt = factor(
        as.character(trt),
        levels = c(
          "Control",
          "NPK"
        )
      )
    ) %>%
    dplyr::filter(
      !is.na(trt)
    ) %>%
    droplevels()
  
  
  #--------------------------------------------------------------------
  # Check data availability by treatment
  #--------------------------------------------------------------------
  
  data_check <- model_dat %>%
    dplyr::group_by(trt) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      
      n_sites = dplyr::n_distinct(
        site_code
      ),
      
      n_precipitation_values = dplyr::n_distinct(
        prec_100
      ),
      
      .groups = "drop"
    )
  
  sufficient_data <- (
    nrow(data_check) == 2 &&
      all(
        c("Control", "NPK") %in%
          as.character(data_check$trt)
      ) &&
      all(
        data_check$n_obs >= min_n_per_treatment
      ) &&
      all(
        data_check$n_sites >= min_sites_per_treatment
      ) &&
      all(
        data_check$n_precipitation_values >= 2
      )
  )
  
  if (!sufficient_data) {
    
    warning(
      paste0(
        "Insufficient data for ",
        current_group,
        " - ",
        current_extreme,
        ". The model was skipped."
      ),
      call. = FALSE
    )
    
    return(
      list(
        fixed_effects = empty_fixed_effects(),
        slopes = empty_slope_results()
      )
    )
  }
  
  
  #--------------------------------------------------------------------
  # Fit the mixed-effects model
  #--------------------------------------------------------------------
  
  model_fit <- tryCatch(
    suppressMessages(
      lmerTest::lmer(
        LRR_cover ~ prec_100 * trt +
          (1 | site_code / block / plot),
        data = model_dat,
        REML = TRUE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 2e5
          )
        )
      )
    ),
    error = function(e) {
      
      warning(
        paste0(
          "Mixed-effects model fitting failed for ",
          current_group,
          " - ",
          current_extreme,
          ": ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (is.null(model_fit)) {
    return(
      list(
        fixed_effects = empty_fixed_effects(),
        slopes = empty_slope_results()
      )
    )
  }
  
  if (
    lme4::isSingular(
      model_fit,
      tol = 1e-4
    )
  ) {
    warning(
      paste0(
        "Singular model fit detected for ",
        current_group,
        " - ",
        current_extreme,
        "."
      ),
      call. = FALSE
    )
  }
  
  model_n_sites <- dplyr::n_distinct(
    model_dat$site_code
  )
  
  model_n_obs <- nrow(model_dat)
  
  
  #--------------------------------------------------------------------
  # Extract fixed-effects parameter estimates
  #--------------------------------------------------------------------
  
  coefficient_table <- as.data.frame(
    summary(model_fit)$coefficients
  ) %>%
    tibble::rownames_to_column(
      var = "term"
    )
  
  required_coefficient_columns <- c(
    "term",
    "Estimate",
    "Std. Error",
    "df",
    "t value",
    "Pr(>|t|)"
  )
  
  missing_coefficient_columns <- setdiff(
    required_coefficient_columns,
    names(coefficient_table)
  )
  
  if (length(missing_coefficient_columns) > 0) {
    stop(
      "The fixed-effects table is missing the following columns: ",
      paste(
        missing_coefficient_columns,
        collapse = ", "
      ),
      "\nAvailable columns: ",
      paste(
        names(coefficient_table),
        collapse = ", "
      )
    )
  }
  
  fixed_effects <- coefficient_table %>%
    dplyr::mutate(
      critical_t = stats::qt(
        p = 0.975,
        df = df
      ),
      
      CI_lower =
        Estimate -
        critical_t * `Std. Error`,
      
      CI_upper =
        Estimate +
        critical_t * `Std. Error`
    ) %>%
    dplyr::transmute(
      Group = current_group,
      
      extreme_type = current_extreme,
      
      term = as.character(
        term
      ),
      
      Estimate = as.numeric(
        Estimate
      ),
      
      SE = as.numeric(
        `Std. Error`
      ),
      
      DF = as.numeric(
        df
      ),
      
      CI_lower = as.numeric(
        CI_lower
      ),
      
      CI_upper = as.numeric(
        CI_upper
      ),
      
      t_value = as.numeric(
        `t value`
      ),
      
      P_value = as.numeric(
        `Pr(>|t|)`
      ),
      
      n_sites = as.integer(
        model_n_sites
      ),
      
      n_obs = as.integer(
        model_n_obs
      )
    )
  
  
  #--------------------------------------------------------------------
  # Estimate treatment-specific precipitation slopes
  #--------------------------------------------------------------------
  
  trend_object <- tryCatch(
    emmeans::emtrends(
      model_fit,
      specs = ~ trt,
      var = "prec_100"
    ),
    error = function(e) {
      
      warning(
        paste0(
          "emtrends() failed for ",
          current_group,
          " - ",
          current_extreme,
          ": ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (is.null(trend_object)) {
    
    slope_results <- empty_slope_results()
    
  } else {
    
    trend_table <- tryCatch(
      as.data.frame(
        summary(
          trend_object,
          infer = c(
            TRUE,
            TRUE
          ),
          level = 0.95,
          adjust = "none",
          null = 0
        )
      ),
      error = function(e) {
        
        warning(
          paste0(
            "Summarizing emtrends() failed for ",
            current_group,
            " - ",
            current_extreme,
            ": ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
        
        NULL
      }
    )
    
    if (
      is.null(trend_table) ||
      nrow(trend_table) == 0
    ) {
      
      slope_results <- empty_slope_results()
      
    } else {
      
      slope_column <- grep(
        "\\.trend$|^trend$",
        names(trend_table),
        value = TRUE
      )
      
      lower_column <- grep(
        "^lower\\.CL$|^asymp\\.LCL$",
        names(trend_table),
        value = TRUE
      )
      
      upper_column <- grep(
        "^upper\\.CL$|^asymp\\.UCL$",
        names(trend_table),
        value = TRUE
      )
      
      statistic_column <- grep(
        "^t\\.ratio$|^z\\.ratio$",
        names(trend_table),
        value = TRUE
      )
      
      if (
        length(slope_column) != 1 ||
        length(lower_column) != 1 ||
        length(upper_column) != 1 ||
        length(statistic_column) != 1 ||
        !"SE" %in% names(trend_table) ||
        !"df" %in% names(trend_table) ||
        !"p.value" %in% names(trend_table)
      ) {
        stop(
          "The required emtrends columns could not be identified.\n",
          "Available columns: ",
          paste(
            names(trend_table),
            collapse = ", "
          )
        )
      }
      
      treatment_counts <- model_dat %>%
        dplyr::group_by(trt) %>%
        dplyr::summarise(
          n_sites = dplyr::n_distinct(
            site_code
          ),
          
          n_obs = dplyr::n(),
          
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          trt = as.character(trt)
        )
      
      slope_results <- tibble::tibble(
        Group = current_group,
        
        extreme_type = current_extreme,
        
        trt = as.character(
          trend_table$trt
        ),
        
        Estimate = as.numeric(
          trend_table[[slope_column]]
        ),
        
        SE = as.numeric(
          trend_table$SE
        ),
        
        DF = as.numeric(
          trend_table$df
        ),
        
        CI_lower = as.numeric(
          trend_table[[lower_column]]
        ),
        
        CI_upper = as.numeric(
          trend_table[[upper_column]]
        ),
        
        t_value = as.numeric(
          trend_table[[statistic_column]]
        ),
        
        P_value = as.numeric(
          trend_table$p.value
        )
      ) %>%
        dplyr::left_join(
          treatment_counts,
          by = "trt"
        ) %>%
        dplyr::mutate(
          n_sites = as.integer(
            n_sites
          ),
          
          n_obs = as.integer(
            n_obs
          )
        )
    }
  }
  
  
  #--------------------------------------------------------------------
  # Return fixed effects and slopes from the same model
  #--------------------------------------------------------------------
  
  list(
    fixed_effects = fixed_effects,
    slopes = slope_results
  )
}


#======================================================================
# Split data by functional group and extreme-event type
#======================================================================

model_datasets <- scores1 %>%
  dplyr::mutate(
    Group = as.character(
      Group
    ),
    
    extreme_type = as.character(
      extreme_type
    )
  ) %>%
  dplyr::group_by(
    Group,
    extreme_type
  ) %>%
  dplyr::group_split(
    .keep = TRUE
  )


#======================================================================
# Fit all mixed-effects models
#======================================================================

model_results <- lapply(
  model_datasets,
  fit_mixed_model
)


#======================================================================
# Combine fixed-effects results
#======================================================================

fixed_effects_results <- dplyr::bind_rows(
  lapply(
    model_results,
    function(x) {
      x$fixed_effects
    }
  )
)


#======================================================================
# Combine treatment-specific slopes
#======================================================================

slope_results <- dplyr::bind_rows(
  lapply(
    model_results,
    function(x) {
      x$slopes
    }
  )
)

if (nrow(fixed_effects_results) == 0) {
  stop(
    "No fixed-effects results were generated."
  )
}

if (nrow(slope_results) == 0) {
  stop(
    "No treatment-specific slope results were generated."
  )
}


#======================================================================
# Arrange the statistical results
#======================================================================

group_order <- c(
  "Annual",
  "Forb",
  "Graminoid",
  "Perennial"
)

extreme_order <- c(
  "Extreme dry",
  "Extreme wet"
)

treatment_order <- c(
  "Control",
  "NPK"
)

fixed_term_order <- c(
  "(Intercept)",
  "prec_100",
  "trtNPK",
  "prec_100:trtNPK"
)

fixed_effects_results <- fixed_effects_results %>%
  dplyr::mutate(
    group_order_variable = factor(
      Group,
      levels = group_order
    ),
    
    extreme_order_variable = factor(
      extreme_type,
      levels = extreme_order
    ),
    
    term_order_variable = factor(
      term,
      levels = fixed_term_order
    )
  ) %>%
  dplyr::arrange(
    group_order_variable,
    extreme_order_variable,
    term_order_variable
  ) %>%
  dplyr::select(
    -group_order_variable,
    -extreme_order_variable,
    -term_order_variable
  )

slope_results <- slope_results %>%
  dplyr::mutate(
    group_order_variable = factor(
      Group,
      levels = group_order
    ),
    
    extreme_order_variable = factor(
      extreme_type,
      levels = extreme_order
    ),
    
    treatment_order_variable = factor(
      trt,
      levels = treatment_order
    )
  ) %>%
  dplyr::arrange(
    group_order_variable,
    extreme_order_variable,
    treatment_order_variable
  ) %>%
  dplyr::select(
    -group_order_variable,
    -extreme_order_variable,
    -treatment_order_variable
  )


#======================================================================
# Validate the statistical results
#======================================================================

all_p_values <- c(
  fixed_effects_results$P_value,
  slope_results$P_value
)

if (
  any(
    !is.na(all_p_values) &
    (
      all_p_values < 0 |
      all_p_values > 1
    )
  )
) {
  stop(
    "One or more P values are outside the valid range from 0 to 1."
  )
}

if (
  any(
    fixed_effects_results$CI_lower >
    fixed_effects_results$CI_upper,
    na.rm = TRUE
  )
) {
  stop(
    "One or more fixed-effect confidence intervals are invalid."
  )
}

if (
  any(
    slope_results$CI_lower >
    slope_results$CI_upper,
    na.rm = TRUE
  )
) {
  stop(
    "One or more treatment-specific slope confidence intervals are invalid."
  )
}


#======================================================================
# Check agreement between slope P values and confidence intervals
#======================================================================

slope_significance_check <- slope_results %>%
  dplyr::mutate(
    p_significant =
      !is.na(P_value) &
      P_value < 0.05,
    
    ci_significant =
      !is.na(CI_lower) &
      !is.na(CI_upper) &
      (
        CI_lower > 0 |
          CI_upper < 0
      ),
    
    p_ci_agree =
      p_significant == ci_significant
  )

if (
  any(
    !slope_significance_check$p_ci_agree,
    na.rm = TRUE
  )
) {
  
  warning(
    paste0(
      "P-value and confidence-interval significance decisions ",
      "do not agree for one or more treatment-specific slopes."
    ),
    call. = FALSE
  )
}




print(
  as.data.frame(
    fixed_effects_results
  ),
  row.names = FALSE
)


print(
  as.data.frame(
    slope_results
  ),
  row.names = FALSE
)


#======================================================================
# Define table-formatting functions
#======================================================================

format_three_decimals <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  non_missing <- !is.na(x)
  
  result[non_missing] <- formatC(
    x[non_missing],
    format = "f",
    digits = 3,
    decimal.mark = "."
  )
  
  result
}

format_p_value <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  non_missing <- !is.na(x)
  
  result[
    non_missing &
      x < 0.001
  ] <- formatC(
    x[
      non_missing &
        x < 0.001
    ],
    format = "e",
    digits = 3,
    decimal.mark = "."
  )
  
  result[
    non_missing &
      x >= 0.001
  ] <- formatC(
    x[
      non_missing &
        x >= 0.001
    ],
    format = "f",
    digits = 3,
    decimal.mark = "."
  )
  
  result
}


#======================================================================
# Prepare treatment-specific slopes for export
#======================================================================

slope_export <- slope_results %>%
  dplyr::transmute(
    Group = Group,
    
    `Extreme type` = extreme_type,
    
    Treat = trt,
    
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
    
    n_sites = as.integer(
      n_sites
    ),
    
    n_obs = as.integer(
      n_obs
    )
  )


#======================================================================
# Display formatted treatment-specific slope table
#======================================================================

cat("\n")
cat("============================================================\n")
cat("FORMATTED TREATMENT-SPECIFIC SLOPE TABLE\n")
cat("Estimates represent changes in LRR per 100 mm\n")
cat("============================================================\n")

print(
  as.data.frame(
    slope_export
  ),
  row.names = FALSE
)


#======================================================================
# Define the Excel output path
#======================================================================

slope_output <- file.path(
  dir.data,
  "Treatment_specific_slopes_functional_groups.xlsx"
)


#======================================================================
# Create the Excel workbook
#======================================================================

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  workbook,
  "Treatment slopes"
)


#======================================================================
# Write the slope table to Excel
#======================================================================

openxlsx::writeData(
  workbook,
  sheet = "Treatment slopes",
  x = slope_export,
  withFilter = TRUE
)


#======================================================================
# Define Excel styles
#======================================================================

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

significant_style <- openxlsx::createStyle(
  textDecoration = "bold"
)

center_style <- openxlsx::createStyle(
  halign = "center",
  valign = "center"
)


#======================================================================
# Format the worksheet
#======================================================================

openxlsx::addStyle(
  workbook,
  sheet = "Treatment slopes",
  style = center_style,
  rows = seq_len(
    nrow(slope_export) + 1
  ),
  cols = seq_len(
    ncol(slope_export)
  ),
  gridExpand = TRUE,
  stack = TRUE
)

openxlsx::addStyle(
  workbook,
  sheet = "Treatment slopes",
  style = header_style,
  rows = 1,
  cols = seq_len(
    ncol(slope_export)
  ),
  gridExpand = TRUE,
  stack = TRUE
)

openxlsx::setColWidths(
  workbook,
  sheet = "Treatment slopes",
  cols = seq_len(
    ncol(slope_export)
  ),
  widths = "auto"
)

openxlsx::freezePane(
  workbook,
  sheet = "Treatment slopes",
  firstRow = TRUE
)


#======================================================================
# Bold significant treatment-specific slopes and P values
#======================================================================

slope_significant_rows <- which(
  !is.na(slope_results$P_value) &
    slope_results$P_value < 0.05
) + 1

slope_significant_columns <- match(
  c(
    "Treat",
    "P_value"
  ),
  names(slope_export)
)

slope_significant_columns <- slope_significant_columns[
  !is.na(slope_significant_columns)
]

if (
  length(slope_significant_rows) > 0 &&
  length(slope_significant_columns) > 0
) {
  
  openxlsx::addStyle(
    workbook,
    sheet = "Treatment slopes",
    style = significant_style,
    rows = slope_significant_rows,
    cols = slope_significant_columns,
    gridExpand = TRUE,
    stack = TRUE
  )
}


#======================================================================
# Save the Excel workbook
#======================================================================

openxlsx::saveWorkbook(
  workbook,
  file = slope_output,
  overwrite = TRUE
)









#======================================================================
# Prepare mixed-effects-model slopes for plotting
#======================================================================

col_climate <- c(
  "Extreme dry" = "#ff9f00",
  "Extreme wet" = "#1f78b4"
)

plot_results <- slope_results %>%
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
    
    # The reversed order produces:
    # Forb, Graminoid, Annual, Perennial
    # from top to bottom in the figure.
    Group = factor(
      Group,
      levels = rev(
        c(
          "Forb",
          "Graminoid",
          "Annual",
          "Perennial"
        )
      )
    ),
    
    significant =
      !is.na(P_value) &
      P_value < 0.05,
    
    sig_symbol = dplyr::case_when(
      !is.na(P_value) &
        P_value < 0.001 ~ "***",
      
      !is.na(P_value) &
        P_value < 0.01 ~ "**",
      
      !is.na(P_value) &
        P_value < 0.05 ~ "*",
      
      TRUE ~ ""
    ),
    
    label_text = dplyr::if_else(
      significant,
      
      paste0(
        "n=",
        n_obs,
        " ",
        sig_symbol
      ),
      
      paste0(
        "n=",
        n_obs
      )
    ),
    
    point_fill = dplyr::if_else(
      significant,
      
      unname(
        col_climate[
          as.character(extreme_type)
        ]
      ),
      
      "white"
    ),
    
    label_fontface = dplyr::if_else(
      significant,
      "bold",
      "plain"
    )
  )


#======================================================================
# Check the plotting data
#======================================================================

if (nrow(plot_results) == 0) {
  stop(
    "No treatment-specific slope results are available for plotting."
  )
}

if (
  any(
    !is.finite(plot_results$Estimate) |
    !is.finite(plot_results$CI_lower) |
    !is.finite(plot_results$CI_upper)
  )
) {
  stop(
    "Some slope estimates or confidence limits are not finite."
  )
}


#======================================================================
# Calculate a common x-axis range for both treatment panels
#======================================================================

confidence_interval_values <- c(
  plot_results$CI_lower,
  plot_results$CI_upper
)

confidence_interval_values <- confidence_interval_values[
  is.finite(confidence_interval_values)
]

if (length(confidence_interval_values) == 0) {
  stop(
    "No finite confidence-interval values are available for plotting."
  )
}

confidence_interval_range <- range(
  confidence_interval_values,
  na.rm = TRUE
)

confidence_interval_width <- diff(
  confidence_interval_range
)

if (
  !is.finite(confidence_interval_width) ||
  confidence_interval_width == 0
) {
  
  confidence_interval_width <- max(
    abs(confidence_interval_range),
    1
  ) * 0.2
}

x_axis_limits <- c(
  confidence_interval_range[1] -
    0.08 * confidence_interval_width,
  
  confidence_interval_range[2] +
    0.65 * confidence_interval_width
)

label_x_position <-
  confidence_interval_range[2] +
  0.08 * confidence_interval_width


#======================================================================
# Define the forest-plot function
#======================================================================

plot_forest_panel <- function(
    data,
    treatment_name,
    show_y_axis = TRUE
) {
  
  panel_data <- data %>%
    dplyr::filter(
      trt == treatment_name
    )
  
  if (nrow(panel_data) == 0) {
    stop(
      "No plotting results are available for treatment: ",
      treatment_name
    )
  }
  
  treatment_title <- if (
    treatment_name == "NPK"
  ) {
    "NPK\u03bc"
  } else {
    "Control"
  }
  
  dodge_position <- ggplot2::position_dodge(
    width = 0.60
  )
  
  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      y = Group,
      x = Estimate,
      color = extreme_type,
      group = extreme_type
    )
  ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.6
    ) +
    
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      orientation = "y",
      width = 0,
      linewidth = 1,
      alpha = 0.85,
      position = dodge_position
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        fill = point_fill
      ),
      shape = 21,
      size = 4.5,
      stroke = 1.3,
      position = dodge_position
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        x = label_x_position,
        label = label_text,
        fontface = label_fontface
      ),
      hjust = 0,
      vjust = 0.5,
      size = 3.8,
      position = dodge_position,
      show.legend = FALSE
    ) +
    
    ggplot2::scale_color_manual(
      values = col_climate,
      drop = FALSE
    ) +
    
    ggplot2::scale_fill_identity() +
    
    ggplot2::coord_cartesian(
      xlim = x_axis_limits,
      clip = "off"
    ) +
    
    ggplot2::labs(
      x = "Slope (change in LRR per 100 mm)",
      y = NULL,
      title = treatment_title,
      color = NULL
    ) +
    
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          linewidth = 1.2,
          size = 4
        )
      )
    ) +
    
    ggplot2::theme_bw(
      base_size = 13
    ) +
    
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      
      legend.text = ggplot2::element_text(
        size = 12
      ),
      
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.6
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 11,
        color = "black"
      ),
      
      axis.text.y = if (show_y_axis) {
        ggplot2::element_text(
          size = 12,
          color = "black"
        )
      } else {
        ggplot2::element_blank()
      },
      
      axis.ticks.y = if (show_y_axis) {
        ggplot2::element_line(
          color = "black"
        )
      } else {
        ggplot2::element_blank()
      },
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 13,
        face = "plain"
      ),
      
      plot.margin = ggplot2::margin(
        t = 5.5,
        r = 15,
        b = 5.5,
        l = 5.5
      )
    )
}


#======================================================================
# Create the Control and NPKμ panels
#======================================================================

control_panel <- plot_forest_panel(
  data = plot_results,
  treatment_name = "Control",
  show_y_axis = TRUE
)

npk_panel <- plot_forest_panel(
  data = plot_results,
  treatment_name = "NPK",
  show_y_axis = FALSE
)


#======================================================================
# Combine the two panels
#======================================================================

fig1 <- (
  control_panel +
    npk_panel
) +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

print(fig1)


#======================================================================
# Save the forest plot
#======================================================================

figure_output_tiff <- file.path(
  dir.data,
  "combined_plot_overall_func.tiff"
)

ggplot2::ggsave(
  filename = figure_output_tiff,
  plot = fig1,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)









#======================================================================
# Functional group × lifespan combinations
# Mixed-effects models and figure only
#======================================================================


#----------------------------------------------------------------------
# Check required packages
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Check the required objects
#----------------------------------------------------------------------

if (!exists("df")) {
  stop(
    "The object 'df' does not exist."
  )
}

if (!exists("dir.data")) {
  stop(
    "The object 'dir.data' does not exist."
  )
}


#----------------------------------------------------------------------
# Define colors and category order
#----------------------------------------------------------------------

col_climate <- c(
  "Extreme dry" = "#ff9f00",
  "Extreme wet" = "#1f78b4"
)

order2 <- c(
  "Annual Forb",
  "Annual Graminoid",
  "Perennial Forb",
  "Perennial Graminoid"
)


#----------------------------------------------------------------------
# Convert functional group × lifespan variables to long format
#----------------------------------------------------------------------

long2 <- df %>%
  dplyr::select(
    dplyr::all_of(
      c(
        keys,
        "climate_display",
        "mean_prec_2000_2020"
      )
    ),
    Graminoid_Annual,
    Forb_Annual,
    Graminoid_Perennial,
    Forb_Perennial
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      Graminoid_Annual,
      Forb_Annual,
      Graminoid_Perennial,
      Forb_Perennial
    ),
    names_to = c(
      "Functional_group",
      "Lifespan"
    ),
    names_pattern =
      "(Graminoid|Forb)_(Annual|Perennial)",
    values_to = "Cover"
  ) %>%
  dplyr::mutate(
    Group_Label = paste(
      Lifespan,
      Functional_group,
      sep = "_"
    )
  )


#----------------------------------------------------------------------
# Calculate normal-year reference cover
#----------------------------------------------------------------------

clean2 <- long2 %>%
  dplyr::group_by(
    site_code,
    block,
    plot,
    trt,
    Group_Label
  ) %>%
  dplyr::mutate(
    n_normal_years = sum(
      climate_display == "Normal" &
        !is.na(Cover),
      na.rm = TRUE
    ),
    
    Yn_cover = if (
      sum(
        climate_display == "Normal" &
        !is.na(Cover),
        na.rm = TRUE
      ) > 0
    ) {
      mean(
        Cover[
          climate_display == "Normal" &
            !is.na(Cover)
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    n_normal_years > 0,
    !is.na(Yn_cover),
    is.finite(Yn_cover)
  )

if (nrow(clean2) == 0) {
  stop(
    "No observations with valid normal-year reference cover are available."
  )
}


#----------------------------------------------------------------------
# Calculate adaptive epsilon
#----------------------------------------------------------------------

eps2 <- clean2 %>%
  dplyr::filter(
    !is.na(Cover),
    is.finite(Cover),
    Cover > 0
  ) %>%
  dplyr::group_by(Group_Label) %>%
  dplyr::summarise(
    eps = max(
      1e-8,
      0.01 * median(
        Cover,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )


#----------------------------------------------------------------------
# Calculate LRR for extreme dry and extreme wet years
#----------------------------------------------------------------------

scores2 <- clean2 %>%
  dplyr::filter(
    climate_display %in% c(
      "Extreme dry",
      "Extreme wet"
    )
  ) %>%
  dplyr::left_join(
    eps2,
    by = "Group_Label"
  ) %>%
  dplyr::filter(
    !is.na(Cover),
    is.finite(Cover),
    
    !is.na(Yn_cover),
    is.finite(Yn_cover),
    
    !is.na(eps),
    is.finite(eps),
    
    !(Cover == 0 & Yn_cover == 0),
    
    !is.na(mean_prec_2000_2020),
    is.finite(mean_prec_2000_2020)
  ) %>%
  dplyr::mutate(
    extreme_type = climate_display,
    
    LRR_cover = dplyr::case_when(
      Cover > 0 & Yn_cover > 0 ~
        log(
          Cover / Yn_cover
        ),
      
      TRUE ~
        log(
          (Cover + eps) /
            (Yn_cover + eps)
        )
    )
  ) %>%
  dplyr::filter(
    !is.na(LRR_cover),
    is.finite(LRR_cover)
  )

if (nrow(scores2) == 0) {
  stop(
    "No valid extreme-event observations are available for modeling."
  )
}


#----------------------------------------------------------------------
# Fit mixed-effects model and extract treatment-specific slopes
#
# Slopes represent changes in LRR per 100-mm increase in mean
# growing-season precipitation.
#----------------------------------------------------------------------

fit_slopes <- function(
    dat,
    min_n_per_treatment = 5,
    min_sites_per_treatment = 3
) {
  
  empty_result <- tibble::tibble(
    trt = character(),
    slope = double(),
    SE = double(),
    DF = double(),
    lower = double(),
    upper = double(),
    t_value = double(),
    p_val = double(),
    n_sites = integer(),
    n = integer()
  )
  
  model_dat <- dat %>%
    dplyr::filter(
      !is.na(LRR_cover),
      is.finite(LRR_cover),
      
      !is.na(mean_prec_2000_2020),
      is.finite(mean_prec_2000_2020),
      
      !is.na(site_code),
      !is.na(block),
      !is.na(plot),
      !is.na(trt)
    ) %>%
    dplyr::mutate(
      prec_100 =
        mean_prec_2000_2020 / 100,
      
      site_code = factor(site_code),
      block = factor(block),
      plot = factor(plot),
      
      trt = factor(
        as.character(trt),
        levels = c(
          "Control",
          "NPK"
        )
      )
    ) %>%
    dplyr::filter(
      !is.na(trt)
    ) %>%
    droplevels()
  
  
  # Check the available data in each treatment
  data_check <- model_dat %>%
    dplyr::group_by(trt) %>%
    dplyr::summarise(
      n = dplyr::n(),
      
      n_sites =
        dplyr::n_distinct(site_code),
      
      n_precipitation_values =
        dplyr::n_distinct(prec_100),
      
      .groups = "drop"
    )
  
  sufficient_data <- (
    nrow(data_check) == 2 &&
      all(
        c(
          "Control",
          "NPK"
        ) %in%
          as.character(data_check$trt)
      ) &&
      all(
        data_check$n >=
          min_n_per_treatment
      ) &&
      all(
        data_check$n_sites >=
          min_sites_per_treatment
      ) &&
      all(
        data_check$n_precipitation_values >= 2
      )
  )
  
  if (!sufficient_data) {
    return(
      empty_result
    )
  }
  
  
  # Fit the mixed-effects model
  model_fit <- tryCatch(
    suppressMessages(
      lmerTest::lmer(
        LRR_cover ~ prec_100 * trt +
          (1 | site_code / block / plot),
        data = model_dat,
        REML = TRUE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 2e5
          )
        )
      )
    ),
    error = function(e) {
      warning(
        "Mixed-effects model fitting failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (is.null(model_fit)) {
    return(
      empty_result
    )
  }
  
  
  # Calculate precipitation slopes by treatment
  trend_object <- tryCatch(
    emmeans::emtrends(
      model_fit,
      specs = ~ trt,
      var = "prec_100"
    ),
    error = function(e) {
      warning(
        "emtrends() failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (is.null(trend_object)) {
    return(
      empty_result
    )
  }
  
  
  # Extract estimates, confidence intervals, and P values
  trend_table <- tryCatch(
    as.data.frame(
      summary(
        trend_object,
        infer = c(
          TRUE,
          TRUE
        ),
        level = 0.95,
        adjust = "none",
        null = 0
      )
    ),
    error = function(e) {
      warning(
        "Summarizing emtrends() failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (
    is.null(trend_table) ||
    nrow(trend_table) == 0
  ) {
    return(
      empty_result
    )
  }
  
  
  # Identify the output columns
  slope_column <- grep(
    "\\.trend$|^trend$",
    names(trend_table),
    value = TRUE
  )
  
  lower_column <- grep(
    "^lower\\.CL$|^asymp\\.LCL$",
    names(trend_table),
    value = TRUE
  )
  
  upper_column <- grep(
    "^upper\\.CL$|^asymp\\.UCL$",
    names(trend_table),
    value = TRUE
  )
  
  statistic_column <- grep(
    "^t\\.ratio$|^z\\.ratio$",
    names(trend_table),
    value = TRUE
  )
  
  if (
    length(slope_column) != 1 ||
    length(lower_column) != 1 ||
    length(upper_column) != 1 ||
    length(statistic_column) != 1 ||
    !"SE" %in% names(trend_table) ||
    !"df" %in% names(trend_table) ||
    !"p.value" %in% names(trend_table)
  ) {
    stop(
      "The required emtrends columns could not be identified.\n",
      "Available columns: ",
      paste(
        names(trend_table),
        collapse = ", "
      )
    )
  }
  
  
  # Count observations and sites within each treatment
  treatment_counts <- model_dat %>%
    dplyr::group_by(trt) %>%
    dplyr::summarise(
      n_sites =
        dplyr::n_distinct(site_code),
      
      n = dplyr::n(),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      trt = as.character(trt)
    )
  
  
  # Return slope results required for plotting
  tibble::tibble(
    trt = as.character(
      trend_table$trt
    ),
    
    slope = as.numeric(
      trend_table[[slope_column]]
    ),
    
    SE = as.numeric(
      trend_table$SE
    ),
    
    DF = as.numeric(
      trend_table$df
    ),
    
    lower = as.numeric(
      trend_table[[lower_column]]
    ),
    
    upper = as.numeric(
      trend_table[[upper_column]]
    ),
    
    t_value = as.numeric(
      trend_table[[statistic_column]]
    ),
    
    p_val = as.numeric(
      trend_table$p.value
    )
  ) %>%
    dplyr::left_join(
      treatment_counts,
      by = "trt"
    ) %>%
    dplyr::mutate(
      n_sites = as.integer(
        n_sites
      ),
      
      n = as.integer(
        n
      )
    )
}


#----------------------------------------------------------------------
# Add information required for plotting
#----------------------------------------------------------------------

add_plot_information <- function(x) {
  
  x %>%
    dplyr::mutate(
      sig_flag =
        !is.na(p_val) &
        p_val < 0.05,
      
      sig_symbol = dplyr::case_when(
        !is.na(p_val) &
          p_val < 0.001 ~ "***",
        
        !is.na(p_val) &
          p_val < 0.01 ~ "**",
        
        !is.na(p_val) &
          p_val < 0.05 ~ "*",
        
        TRUE ~ ""
      ),
      
      label_text = dplyr::if_else(
        sig_flag,
        
        paste0(
          "n=",
          n,
          " ",
          sig_symbol
        ),
        
        paste0(
          "n=",
          n
        )
      ),
      
      fill_col = dplyr::if_else(
        sig_flag,
        
        unname(
          col_climate[
            as.character(extreme_type)
          ]
        ),
        
        "white"
      ),
      
      label_fontface = dplyr::if_else(
        sig_flag,
        "bold",
        "plain"
      )
    )
}


#----------------------------------------------------------------------
# Fit models for each functional group × lifespan combination
#----------------------------------------------------------------------

sens2 <- scores2 %>%
  dplyr::group_by(
    Group_Label,
    extreme_type
  ) %>%
  dplyr::group_modify(
    ~ fit_slopes(
      dat = .x,
      min_n_per_treatment = 5,
      min_sites_per_treatment = 3
    )
  ) %>%
  dplyr::ungroup() %>%
  add_plot_information() %>%
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
    
    Group_Label = factor(
      gsub(
        "_",
        " ",
        Group_Label
      ),
      levels = rev(
        order2
      )
    )
  )

if (nrow(sens2) == 0) {
  stop(
    "No mixed-effects-model results are available for plotting."
  )
}


#----------------------------------------------------------------------
# Remove rows with invalid estimates
#----------------------------------------------------------------------

sens2 <- sens2 %>%
  dplyr::filter(
    is.finite(slope),
    is.finite(lower),
    is.finite(upper)
  )

if (nrow(sens2) == 0) {
  stop(
    "No finite slope estimates are available for plotting."
  )
}


#----------------------------------------------------------------------
# Define the two-panel forest-plot function
#----------------------------------------------------------------------

make_figure <- function(
    sens,
    y_var
) {
  
  confidence_values <- c(
    sens$lower,
    sens$upper
  )
  
  confidence_values <- confidence_values[
    is.finite(confidence_values)
  ]
  
  result_range <- range(
    confidence_values,
    na.rm = TRUE
  )
  
  range_width <- diff(
    result_range
  )
  
  if (
    !is.finite(range_width) ||
    range_width == 0
  ) {
    range_width <- 1
  }
  
  x_limits <- c(
    result_range[1] -
      0.05 * range_width,
    
    result_range[2] +
      0.70 * range_width
  )
  
  label_position <-
    result_range[2] +
    0.07 * range_width
  
  
  create_panel <- function(
    treatment_name,
    show_y_axis
  ) {
    
    panel_data <- sens %>%
      dplyr::filter(
        trt == treatment_name
      )
    
    if (nrow(panel_data) == 0) {
      stop(
        "No results are available for treatment: ",
        treatment_name
      )
    }
    
    panel_title <- if (
      treatment_name == "NPK"
    ) {
      "NPK\u03bc"
    } else {
      "Control"
    }
    
    dodge_position <- ggplot2::position_dodge(
      width = 0.6
    )
    
    ggplot2::ggplot(
      data = panel_data,
      mapping = ggplot2::aes(
        y = .data[[y_var]],
        x = slope,
        color = extreme_type,
        group = extreme_type
      )
    ) +
      
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray40",
        linewidth = 0.6
      ) +
      
      ggplot2::geom_errorbar(
        mapping = ggplot2::aes(
          xmin = lower,
          xmax = upper
        ),
        orientation = "y",
        width = 0,
        linewidth = 1,
        alpha = 0.85,
        position = dodge_position
      ) +
      
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          fill = fill_col
        ),
        shape = 21,
        size = 4.5,
        stroke = 1.3,
        position = dodge_position
      ) +
      
      ggplot2::geom_text(
        mapping = ggplot2::aes(
          x = label_position,
          label = label_text,
          fontface = label_fontface
        ),
        hjust = 0,
        vjust = 0.5,
        size = 3.8,
        position = dodge_position,
        show.legend = FALSE
      ) +
      
      ggplot2::scale_color_manual(
        values = col_climate,
        drop = FALSE
      ) +
      
      ggplot2::scale_fill_identity() +
      
      ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(
            shape = 21,
            fill = unname(
              col_climate
            ),
            size = 4,
            linetype = 0
          )
        )
      ) +
      
      ggplot2::coord_cartesian(
        xlim = x_limits,
        clip = "off"
      ) +
      
      ggplot2::labs(
        x = "Slope (change in LRR per 100 mm)",
        y = NULL,
        title = panel_title,
        color = NULL
      ) +
      
      ggplot2::theme_bw(
        base_size = 13
      ) +
      
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        
        legend.text = ggplot2::element_text(
          size = 12
        ),
        
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        
        axis.line = ggplot2::element_line(
          color = "black",
          linewidth = 0.6
        ),
        
        axis.text.x = ggplot2::element_text(
          size = 11,
          color = "black"
        ),
        
        axis.text.y = if (show_y_axis) {
          ggplot2::element_text(
            size = 12,
            color = "black"
          )
        } else {
          ggplot2::element_blank()
        },
        
        axis.ticks.y = if (show_y_axis) {
          ggplot2::element_line(
            color = "black"
          )
        } else {
          ggplot2::element_blank()
        },
        
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = 13,
          face = "plain"
        ),
        
        plot.margin = ggplot2::margin(
          t = 5.5,
          r = 15,
          b = 5.5,
          l = 5.5
        )
      )
  }
  
  
  control_plot <- create_panel(
    treatment_name = "Control",
    show_y_axis = TRUE
  )
  
  npk_plot <- create_panel(
    treatment_name = "NPK",
    show_y_axis = FALSE
  )
  
  (
    control_plot +
      npk_plot +
      patchwork::plot_layout(
        ncol = 2,
        guides = "collect"
      )
  ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}


#----------------------------------------------------------------------
# Create the figure
#----------------------------------------------------------------------

fig2 <- make_figure(
  sens = sens2,
  y_var = "Group_Label"
)

print(fig2)


#----------------------------------------------------------------------
# Save the figure
#----------------------------------------------------------------------

figure_output <- file.path(
  dir.data,
  "p_final_func_mix.tiff"
)

ggplot2::ggsave(
  filename = figure_output,
  plot = fig2,
  height = 6,
  width = 8.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


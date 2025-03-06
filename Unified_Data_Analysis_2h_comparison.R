# This script orchestrates the analysis of 2-hour exposure experiments across different conditions
# It performs comparative statistical analysis to identify significant differences in thermal tolerance
# between experimental treatments. The script focuses specifically on comparing the effects of 
# pressure on thermal tolerance in marine organisms using a 2-hour exposure protocol.

# Import utility functions and analysis modules
source("plot_utilities.R")        # For plotting and chart management functions
source("non_linear_regression.R") # For fitting sigmoidal survival curves to raw data
source("linear_regression.R")     # For thermal death time analysis and statistical comparisons

# Configure graphics based on operating system
# This ensures consistent rendering across platforms (particularly important for Windows)
check_os_func()

# Main title for charts and analysis output - defines the experiment being visualized
# This title is used in chart headers and file naming
main_title <- "2h tubeworm oxic pressure-nopressure"
#main_title <- "2h popmpeii worm oxic pressure-nopressure"
#main_title <- "2h popmpeii worm anoxic pressure-nopressure"

# Directories where to scan for CSV data - each directory represents a different experimental condition
# For pressure comparison experiments, typically two conditions are compared:
# 1. With high hydrostatic pressure (simulating deep-sea conditions)
# 2. Without additional pressure (control/reference condition)
dir_paths <- list(
  "data/cruise/t_test/tubeworms_mussels/oxic"
  #"data/cruise/t_test/pompeii_worms/oxic"
  #"data/cruise/t_test/pompeii_worms/anoxic"
)

# Initial parameter values for the non-linear regression model
# These starting values are critical for successful convergence of the model fitting
# Parameters represent:
# 1. Maximum survival percentage (typically 100%)
# 2. Temperature at 50% survival (initial guess, typically around 30°C)
# 3. Slope parameter controlling steepness of the survival curve (typically 2-5)
# Multiple sets of parameters are provided for different organism types
nls_param_list <- list(
  list(100, 30, 4)
)

k <- 0
# Iterate over all directories to process each experimental condition
for (dir_path in dir_paths) {
  k <- k + 1
  # Generate a descriptive subtitle based on the directory structure
  chart_subtitle <- chart_subtitle_func(dir_path)
  
  # Perform non-linear regression analysis on all files in the current directory
  # This fits survival curves and extracts LD50 values for each exposure time
  non_linear_regression_result <- non_linear_regression_dir_func(
    dir_path = dir_path,
    nls_param_list = nls_param_list,
    k = k,
    chart_subtitle = chart_subtitle,
    main_title = main_title
  )
}

# ANOVA model analysis test to compare regression coefficients between conditions
# This statistical test determines if there are significant differences in thermal tolerance
# between the experimental conditions (e.g., with vs. without pressure)

# Extract model parameters for statistical comparison
# p1 Estimate from first model (typically parameter for condition 1)
beta1 <- non_linear_regression_result$nls_coefficients[[1]]
# p1 Standard Error from first model
se1 <- non_linear_regression_result$nls_coefficients[[4]]
# p1 Estimate from second model (typically parameter for condition 2)
beta2 <- non_linear_regression_result$nls_coefficients[[13]]
# p1 Standard Error from second model
se2 <- non_linear_regression_result$nls_coefficients[[16]]

# Calculate t-statistic for the difference between parameters
# Formula: t = (β₁ - β₂) / sqrt(SE₁² + SE₂²)
t_statistics <- (beta1 - beta2) / sqrt(se1^2 + se2^2)

# Calculate p-value from t-distribution (two-tailed test)
# Uses minimum degrees of freedom to provide conservative estimate
p_value <- 2 * (1 - pt(abs(t_statistics), non_linear_regression_result$min_df))

# Output statistical results
cat("########## ANOVA ANALYSIS ##########\n\n")
cat("P Value: ", p_value, "\n")
# P-value interpretation:
# - p < 0.05: Significant difference between conditions (reject null hypothesis)
# - p >= 0.05: No significant difference detected (fail to reject null hypothesis) - defines the experiment being visualized
# This title is used in chart headers and file naming
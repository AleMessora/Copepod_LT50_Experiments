if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

source("plot_utilities.R")

# Function for processing a single CSV file using non-linear regression
# Fits a 3-parameter sigmoidal model to the survival data: p1 / (1 + (x / p2)^p3)
# Where:
#   - p1 represents maximum survival (typically 100%)
#   - p2 represents the temperature at which survival is 50% (LD50)
#   - p3 controls the steepness/slope of the curve
#
# Parameters:
#   - dir_path: Path to the directory containing the CSV file
#   - csv_file: Name of the CSV file to process
#   - nls_param_list: List of initial parameter values for the model
#   - time_list: List of exposure times extracted from file names
#   - lines_color: Color palette for visualization
#   - k: Index for accessing the correct initial parameters from nls_param_list
#   - i: Index for the current file in the sequence
#
# Returns a list containing:
#   - LD50 value and its standard error
#   - Model fit parameters
#   - Degrees of freedom for statistical analysis
#   - Data for visualization (points, lines, error bars)
non_linear_regression_file_func <- function(dir_path, csv_file, nls_param_list, time_list, lines_color, k, i) {
  csv_path <- sprintf("%s/%s", dir_path, csv_file)
  cat("########## FILE: ", csv_path, " ###########\n\n")
  data <- read.csv(csv_path, row.names = NULL)
  data$Survival <- data$Alive / (data$Alive + data$Dead)
  
  p1_start <- nls_param_list[[k]][[1]]
  p2_start <- nls_param_list[[k]][[2]]
  p3_start <- nls_param_list[[k]][[3]]
  xdata <- data[, 1]
  ydata <- data[, 4]
  
  # Fit non-linear model to data using the 3-parameter sigmoidal function
  # Maximum iterations increased to ensure convergence for difficult datasets
  fit <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
             start = list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
             control = nls.control(maxiter = 500))
  
  new <- data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
  fit_params <- summary(fit)
  print(summary(fit))
  cat("\n")
  
  ld50 <- c(fit_params$parameters[2])
  
  # Calculate R-squared by comparing residual sum of squares to total sum of squares
  # This measures how well the model explains the variation in the data
  rss <- sum(residuals(fit)^2)
  tss <- sum((ydata - mean(ydata))^2)
  r_squared <- 1 - (rss / tss)
  cat("R squared for fit model: ", r_squared, "\n\n")
  
  model_dose <- fit_params$parameters[2]
  model_dose_y <- predict(fit, newdata = data.frame(xdata = model_dose))
  model_dose_stderr <- fit_params$parameters[2, 2]
  
  # Extract the left part of the CSV filename to check if data needs special treatment
  # This handles the case where multiple experimental conditions are in the same directory
  # with filename format like "Condition_Data01.csv"
  left_part_csv_file <- ""
  if (grepl("_Data", csv_file)) {
    left_part_csv_file <- sub("(.*)_Data.*", "\\1", csv_file)
    left_part_csv_file <- gsub("_", " ", left_part_csv_file)
  }
  
  # Return a comprehensive list of results for further analysis and visualization
  return(list(
    ld50 = ld50,
    fit_params = fit_params,
    df = df.residual(fit),
    new_data = new,
    fit_pred = predict(fit, newdata = new),
    time = ifelse(
      left_part_csv_file == "",
      paste(sprintf("%02d", time_list[i]), "h exposure"),
      paste(left_part_csv_file, "-", sprintf("%02d", time_list[i]), "h exposure")
    ),
    jitter_xdata = jitter(xdata, factor = 0.5),
    ydata = ydata,
    model_dose = model_dose,
    model_dose_y = model_dose_y,
    model_dose_stderr = model_dose_stderr
  ))
}

# Main function for processing all files in a directory through non-linear regression
# Aggregates results from individual files and creates a combined visualization
#
# Parameters:
#   - dir_path: Path to the directory containing CSV files
#   - nls_param_list: List of initial parameter values for the model
#   - k: Index for accessing the correct parameters from nls_param_list
#   - main_title: Main title for the generated charts
#   - chart_subtitle: Subtitle for the generated charts
#   - csv_file_list: Optional list of specific CSV files to process (if NULL, processes all)
#
# Returns a list containing:
#   - LD50 values for all processed files
#   - Time values from file names
#   - Model coefficients for further statistical analysis
#   - Minimum degrees of freedom across all models
#   - Survival function parameters for 3D visualization
non_linear_regression_dir_func <- function(
    dir_path,
    nls_param_list,
    k,
    main_title,
    chart_subtitle,
    csv_file_list = NULL
) {
  cat("########## NON LINEAR REGRESSION ##########\n\n")
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")
  
  # Get list of CSV files to process - either all in directory or from specified list
  if (is.null(csv_file_list)) {
    csv_files <- list.files(path = dir_path)
  } else {
    csv_files <- csv_file_list
  }
  
  # Extract exposure times from filenames using regular expressions
  # Assumes filenames contain numeric characters representing exposure time
  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  time_list <- as.numeric(time_character)
  ld50 <- c()
  nls_coefficients <- list()
  min_df <- .Machine$double.xmax
  survival_params <- c(0, 0, 0, 0)
  
  # Color palette for distinguishing different exposure times in plots
  lines_color <- c("green", "red", "blue", "orange", "purple")
  i <- 0
  
  # Container for accumulated data from all files for combined plotting
  data_list <- list()
  
  # First pass: Process each file to extract model parameters and create curve data
  for (csv_file in csv_files) {
    i <- i + 1
    non_linear_regression_file_result <- non_linear_regression_file_func(
      dir_path = dir_path,
      csv_file = csv_file,
      nls_param_list = nls_param_list,
      time_list = time_list,
      k = k,
      i = i
    )
    ld50 <- c(ld50, non_linear_regression_file_result$ld50)
    nls_coefficients <- append(nls_coefficients, coef(non_linear_regression_file_result$fit_params))
    
    # Track minimum degrees of freedom for later statistical tests
    if (non_linear_regression_file_result$df < min_df) {
      min_df <- non_linear_regression_file_result$df
    }
    
    # Store parameters from first file for 3D visualization
    # This assumes shortest exposure time (typically first file) as reference
    if (i == 1) {
      survival_params <- c(
        non_linear_regression_file_result$fit_params$parameters[1],
        non_linear_regression_file_result$fit_params$parameters[2],
        non_linear_regression_file_result$fit_params$parameters[3],
        time_list[1]
      )
    }
    
    # Store fitted curve data for visualization
    data_list[[i]] <- data.frame(
      xdata = non_linear_regression_file_result$new_data$xdata,
      ydata = non_linear_regression_file_result$fit_pred,
      time = non_linear_regression_file_result$time
    )
  }
  
  # Combine all curve data into a single dataframe for plotting
  plot_data <- do.call(rbind, data_list)
  
  # Create the main plot with fitted curves for all exposure times
  non_linear_regression_plot <- ggplot(plot_data, aes(x = xdata, y = ydata, color = as.factor(time))) +
    geom_line() +
    labs(title = sprintf("LT50 Survival Curve\n%s", chart_subtitle),
         x = "Temperature (°C)", y = "Proportional Survival") +
    scale_color_manual(values = lines_color, name = "Exposure Time (h)")
  
  # Second pass: Add original data points, LD50 markers, and error bars to the plot
  for (i in 1:length(data_list)) {
    non_linear_regression_file_result <- non_linear_regression_file_func(
      dir_path = dir_path,
      csv_file = csv_files[i],
      nls_param_list = nls_param_list,
      time_list = time_list,
      k = k,
      i = i
    )
    
    # Add original sample data points (jittered slightly for better visibility)
    non_linear_regression_plot <- non_linear_regression_plot + geom_point(data = data.frame(
      x = non_linear_regression_file_result$jitter_xdata,
      y = non_linear_regression_file_result$ydata,
      time = non_linear_regression_file_result$time
    ), aes(x = x, y = y, color = as.factor(time)), size = 0.5)
    
    # Add LD50 points (intersection with 50% survival) highlighted with diamond markers
    non_linear_regression_plot <- non_linear_regression_plot + geom_point(data = data.frame(
      x = non_linear_regression_file_result$model_dose,
      y = non_linear_regression_file_result$model_dose_y,
      time = non_linear_regression_file_result$time
    ), aes(x = x, y = y, color = as.factor(time)), shape = 18, size = 4,
    inherit.aes = FALSE)
    
    # Add horizontal error bars to show standard error of LD50 estimates
    non_linear_regression_plot <- non_linear_regression_plot + geom_errorbarh(data = data.frame(
      xmin = non_linear_regression_file_result$model_dose - non_linear_regression_file_result$model_dose_stderr,
      xmax = non_linear_regression_file_result$model_dose + non_linear_regression_file_result$model_dose_stderr,
      y = non_linear_regression_file_result$model_dose_y
    ), aes(xmin = xmin, xmax = xmax, y = y), height = 0.02,
    inherit.aes = FALSE)
    
    # Add temperature labels for each LD50 point
    non_linear_regression_plot <- non_linear_regression_plot + geom_text(data = data.frame(
      x = non_linear_regression_file_result$model_dose,
      y = non_linear_regression_file_result$model_dose_y,
      label = sprintf("%.1f", non_linear_regression_file_result$model_dose)
    ), aes(x = x, y = y, label = label),
    position = position_nudge(y = -0.05),
    size = 2,
    inherit.aes = FALSE)
  }
  
  # Display the plot in the R graphics device
  print(non_linear_regression_plot)
  
  # Save the plot to disk with standardized naming and formatting
  save_plot_func(
    plot = non_linear_regression_plot,
    path = paste("charts/", gsub("\\\\", "/", gsub(" ", "", tools::toTitleCase(trimws(main_title)))), sep = ""),
    filename = paste(gsub(" ", "", trimws(chart_subtitle)), "_lt50", sep = ""),
    width = 1920,
    height = 1080
  )
  
  # Return comprehensive results for further analysis
  return(list(
    ld50 = ld50,
    time_list = time_list,
    nls_coefficients = nls_coefficients,
    min_df = min_df,
    survival_params = survival_params
  ))
}
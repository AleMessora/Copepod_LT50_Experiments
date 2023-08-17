if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

# Initialize dataframes for ANOVA analysis
# Dataframe for LD50
anova_data <- data.frame(
  ld50 = numeric(0),
  time_list = numeric(0),
  dir = character((0))
)
# Dataframe for linear regression params
anova_slopes <- data.frame(
  slope = numeric(0),
  intercept = numeric(0),
  dir = character(0)
)
# Dataframe for T test
t_test_data <- data.frame(
  slope = numeric(0),
  slope_std_err = numeric(0),
  intercept = numeric(0),
  intercept_std_err = numeric(0),
  data_points_length = numeric(0),
  dir = character(0)
)

tdt_line_func <- function(slope, intercept, z) {
  return(intercept + slope * z)
}

linear_regression_func <- function(dir_path, ld50, time_list, chart_subtitle) {
  # Create the Thermal death time curve (TDT) plot

  cat("########## LINEAR REGRESSION ##########\n\n")
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")

  # Create a data frame with LD50 and LT50 to construct TDT curve
  lt50 <- data.frame(ld50, time_list)
  lt50_temp <- lt50
  lt50_temp$dir <- chart_subtitle
  anova_data <- rbind(anova_data, lt50_temp)
  cat("LD50 values and corresponding LT50 exposure values\n")
  print(lt50)
  cat("\n")

  # Run a linear regression with time in hours
  tdt_results_hours <- lm(lt50$ld50 ~ log10(lt50$time_list), data = lt50)
  tdt_slope <- tdt_results_hours$coefficients["log10(lt50$time_list)"]
  tdt_intercept_hours <- tdt_results_hours$coefficients["(Intercept)"]

  new_anova_element <- data.frame(
    slope = tdt_slope,
    intercept = tdt_intercept_hours,
    dir = chart_subtitle
  )
  anova_slopes <- rbind(anova_slopes, new_anova_element)

  new_t_test_data <- data.frame(
    slope = summary(tdt_results_hours)$coefficients[2, 1],
    slope_std_err = summary(tdt_results_hours)$coefficients[2, 2],
    intercept = summary(tdt_results_hours)$coefficients[1, 1],
    intercept_std_err = summary(tdt_results_hours)$coefficients[1, 2],
    data_points_length = length(lt50_temp$time_list),
    dir = dir_path
  )
  t_test_data <- rbind(t_test_data, new_t_test_data)

  # Create plot with LD50/LT50 values in log-scale
  plot(log10(lt50$time_list), lt50$ld50,
    xlim = c(0, 2), ylim = c(0, 55),
    xlab = expression("log"[10] * "Time (LT50, hours)"),
    ylab = "Temperature (LD50, °C)",
    main = paste("Thermal death time (TDT) Curve\n", chart_subtitle)
  )

  # Create a vector of LD50 values and use it as input to plot the curve
  z <- seq(0, 2, 0.01)
  lines(z, tdt_line_func(tdt_slope, tdt_intercept_hours, z), new = TRUE)

  # Calculate R squared goodness of fit for the TDT curve
  cat("Summary statistics for TDT curve\n")
  print(summary(tdt_results_hours))
  cat("R squared value for TDT curve: ", summary(tdt_results_hours)$r.squared, "\n\n")

  # Run a linear regression with time in minutes
  lt50_minutes <- lt50
  lt50_minutes$time_list <- lt50_minutes$time_list * 60
  tdt_results_minutes <- lm(lt50_minutes$ld50 ~ log10(lt50_minutes$time_list), data = lt50_minutes)
  tdt_intercept_minutes <- tdt_results_minutes$coefficients["(Intercept)"]
  cat("TDT intercept in minutes: ", tdt_intercept_minutes, "\n\n")
  return(list(
    anova_data = anova_data,
    anova_slopes = anova_slopes,
    t_test_data = t_test_data,
    survival_param_z = tdt_slope
  ))
}

anova_analysis <- function(anova_data, anova_slopes) {
  # Plot ANOVA analysis
  anova_data$time_list <- log10(anova_data$time_list)
  anova_plot <- ggplot(
    anova_data,
    aes(
      x = time_list,
      y = ld50,
      color = factor(dir)
    )
  ) +
    geom_point() +
    labs(
      x = "LOG10 Time (LT50, hours)",
      y = "Temperature (LD50, °C)",
      title = "Thermal death time (TDT) Curve"
    ) +
    scale_color_discrete(name = "Sample")
  for (i in 1:nrow(anova_slopes)) { # nolint: seq_linter.
    z <- seq(0, 2, 0.01)
    anova_slope_data <- tdt_line_func(anova_slopes[i, "slope"], anova_slopes[i, "intercept"], z)
    anova_slope_df <- data.frame(ld50 = anova_slope_data, time_list = z, dir = anova_slopes[i, "dir"])
    anova_plot <- anova_plot +
      geom_line(data = anova_slope_df, aes(group = dir))
  }
  print(anova_plot)
  anova_analysis <- anova(lm(ld50 ~ log(time_list) * dir, anova_data))
  cat("########## ANOVA ANALYSIS ##########\n\n")
  print(anova_analysis)
  cat("\n")
}

t_test_func <- function(t_test_data) {
  # T test slope
  t_value_slope <-
    (t_test_data$slope[1] - t_test_data$slope[2]) /
    sqrt(t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)
  df_value_slope <-
    (t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)^2 /
    (t_test_data$slope_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
     t_test_data$slope_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
  p_value_slope <- 2 * (1 - pt(abs(t_value_slope), df_value_slope))

  t_value_intercept <-
    (t_test_data$intercept[1] - t_test_data$intercept[2]) /
    sqrt(t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)
  df_value_intercept <-
    (t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)^2 /
    (t_test_data$intercept_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
     t_test_data$intercept_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
  p_value_intercept <- 2 * (1 - pt(abs(t_value_intercept), df_value_intercept))
  cat("########## T TEST ##########\n\n")
  cat("P Value slopes: ", p_value_slope, "\n")
  cat("P Value intercept: ", p_value_intercept, "\n\n")
}
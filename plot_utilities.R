if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

# Function to detect operating system and configure appropriate graphics driver
# Windows requires Cairo for proper rendering of complex plots
# This ensures consistent output across different operating systems
check_os_func <- function() {
  sys_info <- Sys.info()
  os_type <- sys_info["sysname"]
  if (os_type == "Windows") {
    if (!require("Cairo")) {
      install.packages("Cairo")
      library(Cairo)
    }
    CairoWin()
  }
}

# Function to generate a standardized subtitle for charts based on the directory path
# Extracts the last two parts of the path, formats them for display by:
# - Replacing slashes with hyphens
# - Replacing underscores with spaces
# - Converting to title case for readability
chart_subtitle_func <- function(dir_path) {
  chart_subtitle_parts <- unlist(strsplit(dir_path, "/"))
  chart_subtitle <- paste(
    chart_subtitle_parts[(length(chart_subtitle_parts) - 1):length(chart_subtitle_parts)],
    collapse = "/"
  )
  chart_subtitle <- gsub("/", " - ", chart_subtitle)
  chart_subtitle <- gsub("_", " ", chart_subtitle)
  chart_subtitle <- tools::toTitleCase(chart_subtitle)
  return(chart_subtitle)
}

# Function to save plots with consistent settings across the project
# Parameters:
# - plot: The ggplot object to save
# - path: The relative path where the plot should be saved
# - filename: Base filename without extension
# - width/height: Dimensions in pixels
#
# Creates directory if it doesn't exist
# Saves as high-resolution (300 DPI) PNG with Cairo renderer for best quality
save_plot_func <- function(plot, path, filename, width, height) {
  dpi <- 300
  width_in_inches <- width / dpi
  height_in_inches <- height / dpi
  # Determine current script directory to use as base path
  script_path <- normalizePath(sys.frame(1)$ofile)
  script_dir <- dirname(script_path)
  plot_dir <- file.path(script_dir, path)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  ggsave(
    plot = plot,
    width = width_in_inches,
    height = height_in_inches,
    dpi = dpi,
    path = plot_dir,
    filename = paste(filename, ".png", sep = ""),
    type = "cairo-png"
  )
}
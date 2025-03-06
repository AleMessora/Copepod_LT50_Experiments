# Thermal Tolerance Analysis Framework

This repository contains a comprehensive R framework for analyzing thermal tolerance in marine organisms through LD50/LT50 experiments. The framework provides tools for data analysis, statistical modeling, and visualization of thermal tolerance curves.

## Overview

The framework processes experimental data where organisms are exposed to different temperatures for various durations, analyzing survival rates to determine lethal temperatures (LD50) and lethal times (LT50). The analysis produces three types of visualizations:

1. **LT50 Survival Curves**: Shows the relationship between temperature and survival probability for different exposure durations
2. **Thermal Death Time (TDT) Curves**: Plots LD50 temperatures against log-transformed exposure times
3. **3D Thermal Tolerance Landscapes**: Visualizes the survival rates across both temperature and time dimensions

## Key Components

### 1. plot_utilities.R

This module provides utility functions for creating and saving plots:

- `check_os_func()`: Detects the operating system and loads appropriate graphics drivers (Cairo for Windows)
- `chart_subtitle_func(dir_path)`: Generates descriptive subtitles for charts based on directory paths
- `save_plot_func(plot, path, filename, width, height)`: Saves plots with consistent settings

### 2. non_linear_regression.R

This module handles non-linear regression modeling for thermal tolerance data:

#### Key Functions:

- `non_linear_regression_file_func(dir_path, csv_file, nls_param_list, time_list, lines_color, k, i)`:
  - Processes a single CSV file containing survival data
  - Applies non-linear regression using the model: `survival = p1 / (1 + (temperature / p2)^p3)`
  - Returns model parameters, LD50 values, and visualization data

- `non_linear_regression_dir_func(dir_path, nls_param_list, k, main_title, chart_subtitle, csv_file_list)`:
  - Iterates through all CSV files in a directory (or specified files)
  - Aggregates results from individual file analyses
  - Creates combined LT50 survival curves with error bars and statistics
  - Returns parameters needed for further analysis (3D visualization, linear regression)

The non-linear regression uses a specialized sigmoid model where:
- `p1` represents the maximum survival rate
- `p2` is the temperature at which survival is reduced by half (LD50)
- `p3` controls the slope of the curve (steepness of the decline)

### 3. linear_regression.R

This module analyzes thermal death time patterns using linear regression on LD50 values:

#### Key Functions:

- `tdt_line_func(slope, intercept, z)`: Helper function to calculate points along the TDT line

- `linear_regression_func(dir_path, ld50, time_list, main_title, chart_subtitle)`:
  - Takes LD50 values from the non-linear regression step
  - Performs linear regression of LD50 against log-transformed exposure times
  - Generates TDT curves
  - Returns parameters for statistical comparisons

- `anova_analysis_func(anova_data, anova_slopes, main_title)`:
  - Performs ANOVA analysis to compare slopes between different experimental conditions
  - Visualizes multiple TDT curves for comparison

- `t_test_func(t_test_data, main_title)`:
  - Performs t-tests to compare regression slopes and intercepts between experiments
  - Calculates p-values to determine statistical significance of differences

### 4. Unified_Data_Analysis_2h_comparison.R

This script orchestrates the comparative analysis of 2-hour exposure experiments across different conditions:

#### Key Steps:

1. Sources necessary utility and analysis modules
2. Checks the operating system for proper graphics setup
3. Specifies the main title and data directories to analyze
4. Sets initial parameter values for non-linear regression models
5. Iterates through directories, performing for each:
   - Non-linear regression to determine LD50 values
   - Statistical tests to compare results between conditions
6. Performs ANOVA analysis to test for significant differences between experimental conditions
7. Calculates and reports p-values for statistical comparisons

## Data Requirements

- CSV files must follow a specific naming convention: `Dataxx.csv` where xx is a number (01, 02, 03... 10)
- Each CSV file should contain columns for Temperature, Alive, and Dead organisms
- Files should be organized in directories according to experimental conditions

## Usage

To perform a thermal tolerance analysis:

1. Organize your data in appropriate directories following the required format
2. Adjust parameters in the main script for your specific experiment:
   ```R
   main_title <- "Your experiment title"
   dir_paths <- list("path/to/your/data/directory")
   nls_param_list <- list(list(100, 30, 4))  # Initial parameter guesses
   ```
3. Run the script using Ctrl+Shift+S (Source) to see statistics in the console
4. Generated plots will be saved in the `charts/` directory

## Statistical Analysis

The framework provides comprehensive statistical analysis capabilities:

- R-squared values for goodness-of-fit of non-linear and linear models
- ANOVA tests to compare regression lines between different experimental conditions
- T-tests to specifically compare slopes and intercepts of TDT curves
- Standard error calculations and visualization for LD50 estimates

## Visualization Outputs

All visualizations are automatically saved in high resolution (300 DPI) PNG format with consistent styling:

- **LT50 Curves**: Show proportional survival vs. temperature for each exposure time
- **TDT Curves**: Plot LD50 temperatures against log-transformed exposure time
- **Comparative plots**: Overlay results from different experimental conditions for direct comparison

## Example Applications

The repository includes several application scripts:
- `Unified_Data_Analysis_FR.R`: Full-range analysis with 3D visualization
- `Unified_Data_Analysis_FR_anova.R`: ANOVA comparisons across multiple conditions
- `Unified_Data_Analysis_FR_single-csv.R`: Analysis of a single data file
- `Unified_Data_Analysis_2h_comparison.R`: Specific comparison of 2-hour exposure experiments

## Mathematical Models

The framework employs two primary mathematical models:

1. **Non-linear sigmoid model for survival rate**:
   ```
   Survival = p1 / (1 + ((Temperature - (z * log10(Time / t0))) / p2)^p3)
   ```
   
2. **Linear model for Thermal Death Time**:
   ```
   LD50 = intercept + slope * log10(Time)
   ```

Where:
- `p1`, `p2`, `p3` are parameters of the sigmoid curve
- `z` is the temperature coefficient from TDT analysis
- `t0` is the reference time point

## License

This project is licensed under CC0 1.0 Universal, placing it in the public domain.

 # ----------------------------------------------------------------------------------------------
 library(pacman)     # for silent loading of libraries
 INSTALL = FALSE
 # Load required libraries:
 # Data manipulation:
 p_load(DBI, install = INSTALL)          # R Database Interface
 p_load(lubridate, install = INSTALL)    # Make Dealing with Dates a Little Easier
 p_load(magrittr, install = INSTALL)     # A Forward-Pipe Operator for R
 p_load(tidyverse, install = INSTALL)    # Easily Install and Load the 'Tidyverse'
 p_load(tibbletime, install = INSTALL)   # Time Aware Tibbles
 p_load(openxlsx, install = INSTALL)     # Read, Write and Edit XLSX Files
 p_load(readxl, install = INSTALL)       # Read Excel Files
 p_load(reshape, install = INSTALL)      # Flexibly Reshape Data
 p_load(stringi, install = INSTALL)      # Character String Processing Facilities
 # Plots and visuals:
 p_load(ggplot2, install = INSTALL)      # Create Elegant Data Visualisations Using the Grammar of Graphics
 p_load(gridExtra, install = INSTALL)    # Miscellaneous Functions for "Grid" Graphics
 p_load(ggplotify, install = INSTALL)    # Convert Plot to 'grob' or 'ggplot' Object
 p_load(scales, install = INSTALL)       # Scale Functions for Visualization
 p_load(ztable, install = INSTALL)       # Zebra-Striped Tables in LaTeX and HTML Formats
 options(ztable.type = "viewer")         # use "html" or "latex" for extraction of the codes
 # Time series models:
 p_load(dlm, install = INSTALL)          # Bayesian and Likelihood Analysis of Dynamic Linear Models
 p_load(glmnet, install = INSTALL)       # Lasso and Elastic-Net Regularized Generalized Linear Models
 p_load(forecast, install = INSTALL)     # Forecasting Functions for Time Series and Linear Models
 p_load(imputeTS, install = INSTALL)     # Time Series Missing Value Imputation
 p_load(fable, install = INSTALL)        # Forecasting Models for Tidy Time Series
 p_load(mFilter, install = INSTALL)      # Miscellaneous Time Series Filters
 p_load(seasonal, install = INSTALL)     # R Interface to X-13-ARIMA-SEATS
 p_load(tempdisagg, install = INSTALL)   # Methods for Temporal Disaggregation and Interpolation of Time Series
 p_load(xts, install = INSTALL)          # eXtensible Time Series (extended zoo)
 p_load(WaveletComp, install = INSTALL)  # Computational Wavelet Analysis
 p_load(zoo, install = INSTALL)          # S3 Infrastructure for Regular and Irregular Time Series
 p_load(pls, install = INSTALL)          # Partial Least Squares and Principal Component Regression
 # Anomalies detection:
 p_load(corrplot, install = INSTALL)     # Visualization of a Correlation Matrix
 p_load(tsoutliers, install = INSTALL)   # Detection of Outliers in Time Series
 # Transformations
 p_load(TTR, install = INSTALL)          # Technical Trading Rules
 p_load(coda, install = INSTALL)         # Output Analysis and Diagnostics for MCMC
 p_load(urca, install = INSTALL)         # Unit Root and Cointegration Tests for Time Series Data
 library(MASS, exclude = "select")       # Support Functions and Datasets for Venables and Ripley's MASS
 p_load(dplyr, install = INSTALL)        # A Grammar of Data Manipulation
 p_load(dfoptim, install = INSTALL)      # Derivative-Free optimization algorithms
 p_load(mFilter, install = INSTALL)      # Miscellaneous Time Series Filters
 # Time series analysis:
 p_load(BCDating, install = INSTALL)     # Business Cycle Dating and Plotting Tools
 
 # ----------------------------------------------------------------------------------------------
 # Load required functions:
 setwd(paste(path, "/functions", sep = ""))
 accuracy             <- dget("accuracy.R")
 add.lags             <- dget("add.lags.R")
 arima_analysis       <- dget("arima_analysis.R")
 apply_BBQ            <- dget("apply_BBQ.R")
 bivariate_analysis   <- dget("bivariate_analysis.R")
 capwords             <- dget("capwords.R")
 create.dir           <- dget("create.dir.R")
 create_lagged_vars   <- dget("create_lagged_vars.R")
 cv_pca               <- dget("cv_pca.R")
 cv_relaxed           <- dget("cv_relaxed.R")
 cv_pls               <- dget("cv_pls.R")
 dlmMLE_m             <- dget("dlmMLE_m.R")
 find_model           <- dget("find_model.R")
 forecast_oos         <- dget("forecast_oos.R")            
 get_country_data     <- dget("get_country_data.R")
 hp.filter            <- dget("hp.filter.R")
 hp.q                 <- dget("HP.Q.R")
 hp_smoother          <- dget("hp_smoother.R")
 lm_analysis          <- dget("lm_analysis.R")
 max_column           <- dget("max_column.R")
 numeric_to_yqrt      <- dget("numeric_to_yqrt.R")
 prescreen_features   <- dget("prescreen_features.R")
 reverse_forecast     <- dget("reverse_forecast.R")
 reverse_ts           <- dget("reverse_ts.R")
 right_na             <- dget("right_na.R")
 run_arima_with_lasso <- dget("run_arima_with_lasso.R")
 run_lasso            <- dget("run_lasso.R")
 sc.length            <- dget("sc.length.R")
 scale_over_fullcycle <- dget("scale_over_fullcycle.R")
 select_model         <- dget("select_model.R")
 simulate.cts         <- dget("simulating.cyclical.time.series.r")
 stochastic.cycle     <- dget("stochastic.cycle.R")
 stochastic.cycle.nc  <- dget("stochastic.cycle.nc.R")
 tcd_eda              <- dget("tcd_eda.R")
 tcd.one              <- dget("tcd.one.R")
 theme_custom         <- dget("theme_custom.R")
 trim_and_smooth      <- dget("trim_and_smooth.R")
 uc.ecb               <- dget("ucm.ecb.R")

 
 source("tcs.R") # Matlab routines with auxiliary matrix operations functions translated into R
 # ----------------------------------------------------------------------------------------------
 # Some custom colours:
 col.custom <- c("#63003d", "#f06b94", "#6eaad2")
 col.lg.rev     <- colorRampPalette(c("#f06b94", "white","#6eaad2"))
 col.lg.countries <- c(EE = "#99C6DC", LV = "#DB8993", LT = "#A8B596", BB = "#36102D")
 theme_set(theme_custom()) # theme_bw()
 
 setwd(path)
 

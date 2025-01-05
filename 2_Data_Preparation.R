
 # Description:   - Prepares Z input data
 #                - Prepares macro data:
 #                  1) keeps only data sufficiently long
 #                  2) extracts growth cycle out of this data
 #                - Stores some exploratory data analysis

 
 #  TCD decision log:
 #  Keep 3 methods in the mini suite of models, use UCM for computation of cyclical properties:
 #    1) CF filter full band (6:80) (sub-bands for driving cycle weights analysis)
 #    2) HP filter with optimal lambda cycle length dependent
 #    3) Prior consistent TCS filter Mohr (2004) with updated rho and cycle length 
 #    4) 1-3 combination
 #  Heuristics:
 #   if the number of diffs (ND) = 2, then CF and HP could have larger leakage effects
 #   if max(abs(range(cycle_final)*100)) < [0.5, 1.5] threshold, then conclude 
 #   that data has no severe cyclical component relative importance factor 
 #   based on variance explained also could be used
 
 # Continue from Step 1 (may require re-run)
 t3 <- Sys.time()
 
 # Re-load required libraries, install missing ones:
 source("@CCFunctions.R")
 
 # Make a bulk analysis function start with LT
 CN <- 'LV'
 CN_long <- 'Latvia' # for plots only
 PL <- 6   # lower limit in quarters
 PU <- 80  # upper limit in quarters up to 20 years cycles (~ ECB WGEM paper (2018))
 REVIEW_SCALING  <- FALSE     # flag to re-run BBQ full cycle periods analysis
 scaling_folder  <- "2024 Q3" # provide scaling folder in Results, from where to take periods
 SCALED_ANALYSIS <- FALSE     # plot and print prepared macro data analysis for documentation
 
 data <- cc_data[[CN]]
 data_pessimistic <- cc_data_pessimistic[[CN]]
 data_optimistic  <- cc_data_optimistic[[CN]]
 odf_z <- odf_z_data[[CN]]
 
 z_data_le <- odf_z$le
 z_data_le$Date[1] %>% {c(year(.), quarter(.))} -> z_start_date_le
 z_data_le$Z %>% ts(start = z_start_date_le, frequency = 4) -> z_le
 
 z_data_pi <- odf_z$pi
 z_data_pi$Date[1] %>% {c(year(.), quarter(.))} -> z_start_date_pi
 z_data_pi$Z %>% ts(start = z_start_date_pi, frequency = 4) -> z_pi
 
 # Find time series starting points (sufficiently long if starts before 2007 Q1):
 data %>% 
 lapply(function(x) x %>% na.trim() %>% {tsp(.)[1]}) %>% 
 unlist() %>% 
 yearqtr()
 
 # Define feature candidates for trend cycle decomposition G1 and G2
 TCD_Variables <- c("HIC", "HPI", "WAG", "WAR", "YER", 
 "CFR", "CGR", "CHR", "GCR", "ITR", "JVN", "MTR",
 "OMX", "PCR", "WIR", "XTR", 
 "d4HIC", "d4HPI", "d4WAG", "d4WAR")
 TCD_Variables %>%
 lapply(function(x) tcd.one(var_nm = x,
 show_plots = FALSE, # no plots for faster computation 
 use_all_tcs = FALSE
 )) -> CC_macro_cycles
 
 # See "Heuristics" in description of the script:
 CC_macro_cycles %>% lapply(function(x) x$manual) %>% unlist() %>% which() -> problems
 CC_macro_cycles %>% lapply(function(x) x$sig_cycle) %>% unlist() %>% which() -> sig_cycles
 
 if(sum(problems %in% sig_cycles) > 0L){                                                     
 problems_to_check <- intersect(problems, sig_cycles)
 warning("Variables: ", paste(TCD_Variables[problems_to_check], collapse = " "), 
 " require manual judgement before approving the results!", 
 immediate. = TRUE)
 # Example of single analysis, say 5 -- YER is problematic:
 problem_index <- 5
 var_nm <- TCD_Variables[problem_index]
 y <- na.omit(data[, var_nm])
 tcd_one <- tcd.one(var_nm, show_plots = TRUE)
 test.data <- data.frame(date = as.yearmon(time(y)), y = y)
 # use 2^j, j = 1, 2, \dots for upper lower limits; 
 # dj controls for more pixels detected (more makes slower, but smoother detection of features)
 # loess.span shows how local the trend approximation should be
 # dt shows one unit resolution, if 24hrs 1/24 setting will turn a unit to day
 my.w <- analyze.wavelet(test.data, "y", loess.span = 0.8, dt = 1, dj = 1/22, 
 lowerPeriod = PL, upperPeriod = PU, make.pval = TRUE, n.sim = 100)
 wt.image(my.w, color.key = "quantile", n.levels = 100, color.palette = "col.lg.rev(n.levels)", 
 legend.params = list(lab = "wavelet power levels", mar = 2), show.date = TRUE,
 periodlab  = "periods (quarters)")
 wt.avg(my.w)
 reconstruct(my.w, plot.waves = FALSE, lwd = c(1, 2))
 }
 
 # Extract cycles into separate data frame:
 CC_macro_cycles %>% lapply(function(x)x$cycles[, ncol(x$cycles)]) -> cc.list
 CC_macro_cycles %>% lapply(function(x)x$trends[, ncol(x$trends)]) -> te.list
 names(cc.list) <- names(te.list) <- TCD_Variables
 CC_macro_suite   <- NULL
 TE_macro_suite   <- NULL #errors and trends
 for(i in 1:length(te.list)) {
 CC_macro_suite <- ts.union(CC_macro_suite, cc.list[[i]] / 100) 
 TE_macro_suite <- ts.union(TE_macro_suite, te.list[[i]])
 }
 colnames(CC_macro_suite) <- colnames(TE_macro_suite) <- TCD_Variables
 
 # Define G3 group of cyclical variables:
 Scale_only_Variables <- c("E6M", "URX", "CCI", "ESI", "FCI", "ICI", "RCI", "SCI")
 
 # Remove noise below one year frequency applying Hodrick-Prescott filter:
 data[, Scale_only_Variables] %>%
 lapply(function(x){x <- na.omit(x);return(hp_smoother(x))}) -> smoothed_ts
 smoothed_data <- smoothed_ts[[1]]
 for(index in 2:length(smoothed_ts)){
 smoothed_data %<>% {ts.union(., smoothed_ts[[index]])}
 }
 removed_noise <- data[, Scale_only_Variables] - smoothed_data
 colnames(removed_noise) <- colnames(smoothed_data) <- Scale_only_Variables
 
 cyclical_data <-  ts.union(z_le, z_pi, CC_macro_suite, smoothed_data)
 colnames(cyclical_data) <- c("Z_LE", "Z_PI", colnames(CC_macro_suite), colnames(smoothed_data))
 
 ### Scaling over full cycle period only
 # Apply peak-to-peak / trough-to-trough full cycle period scaling
 # Review scaling periods at most once per year else use saved intervals for each country
 scaling_periods_path <- paste0("../4_Results/", scaling_folder, "/ScalingPeriods_", CN,".csv")
 if(REVIEW_SCALING){
 BBQ_results <- apply_BBQ(data = cyclical_data, max_date = 2024) 
 # reminder, 2024 == "2024 Q1", Q step by 0.25
 # Prepare the data frame for CSV export
 scaling_periods <- data.frame()
 
 # Populate the data frame with only scaling period information
 for (i in seq_along(BBQ_results)) {
 if (!is.na(BBQ_results[[i]]) && !is.null(BBQ_results[[i]]$scaling_period)) {
 scaling_periods <- rbind(scaling_periods, c(colnames(cyclical_data)[i], 
 BBQ_results[[i]]$all_data))
 }
 }
 
 scaling_periods$Overriden <- 0
 colnames(scaling_periods) <- c("Variable", "Min", "Max", 
 "Min_Peak", "Max_Peak", "Min_Trough",
 "Max_Trough", "Start", "End", "Overridden")
 
 # Write the data frame to a CSV file 
 write.csv(scaling_periods, scaling_periods_path, row.names = FALSE)
 # Open csv override starting or end periods applying expert judgement
 # Base your judgement on visual inspection of BBQ plots
 # Mark any overridden periods with 1 in "Overridden" column
 # Save and close csv and rerun scaling by loading the file (standard scaling approach)
 }
 # load periods data for country CN
 scaling_periods_path %>%
 read.csv(stringsAsFactors = FALSE) %>%
 select(Variable, Start, End) -> scaling_periods
 
 list_all_scaled <- scale_over_fullcycle(cyclical_data, scaling_periods)
 
 all_scaled_data <- list_all_scaled$data
 scaled_center   <- list_all_scaled$scaled_center
 scaled_scale    <- list_all_scaled$scaled_scale
 
 # Unsupervised learning analysis of the data 
 if(SCALED_ANALYSIS){
 # All in one plot
 autoplot(all_scaled_data, color = "#6eaad2", alpha = 0.5, xlab = "", ylab = "") +
 ggtitle(paste(CN_long, "all scaled macroeconomic data cycles", sep = ",")) +
 theme_custom()
 
 common_scaled_data <- na.omit(all_scaled_data)
 # Correlation analysis, choose small subgroups for that
 set.seed(2024)
 M <- cor(common_scaled_data, use = "pairwise")
 corrplot(M, method= "ellipse", order = "hclust", addrect = 4, # 2 or 3 clusters
 rect.col = "darkblue", col = col.lg.rev(22), hclust.method = "ward.D2",
 insig = "blank", tl.col  = "black", tl.cex = 0.5)
 
 group_pca <- prcomp(common_scaled_data,
 center = TRUE,
 scale. = TRUE)
 # how many PC to display?
 variances <- data.frame(x = 1:10, y = head((group_pca$sdev)^2, 10))  
 ggplot(variances, aes(x = x, y = y)) +
 geom_point() +
 geom_line() +
 labs(x = "", y = "") + 
 scale_x_continuous(breaks = seq(0, 10, by = 1))
 summary(group_pca)
 biplot(group_pca, cex = 0.75)
 group_pca$x %>% ts(start = tsp(common_scaled_data)[1], frequency = 4) -> group_pca_ts
 # The Economic composite indicator could be defined as the first PC
 eci <- -ts(scale(group_pca_ts[,1]), start = tsp(common_scaled_data)[1], frequency = 4)
 autoplot(common_scaled_data, color = "#6eaad2", alpha = 0.5, xlab = "", ylab = "") +
 autolayer(eci, color = 1, lwd = 1.2)
 }
 
 #extended group of cycles G1
 G1 <- c("E6M", "URX", "HIC", "HPI", "WAG", "WAR", "YER", "d4HPI", "d4HIC", "d4WAR", "d4WAG") 
 exclusions <- ts.union(TE_macro_suite, removed_noise)
 colnames(exclusions) <- c(colnames(TE_macro_suite), colnames(removed_noise))
 exclusions_selected <- exclusions[, G1]
 
 scaled_m  <- scaled_center[G1]
 scaled_sd <- scaled_scale[G1]
 
 baseline <- all_scaled_data[, G1]
 
 pessimistic <- cc_data_pessimistic[[CN]][, G1] - exclusions_selected
 pessimistic <- (pessimistic - matrix(scaled_m, nrow = nrow(pessimistic), ncol = length(G1), byrow = TRUE)) / 
 matrix(scaled_sd, nrow = nrow(pessimistic), ncol = length(G1), byrow = TRUE)
 colnames(pessimistic) <- G1
 
 optimistic <- cc_data_optimistic[[CN]][, G1] - exclusions_selected
 optimistic <- (optimistic - matrix(scaled_m, nrow = nrow(optimistic), ncol = length(G1), byrow = TRUE)) / 
 matrix(scaled_sd, nrow = nrow(optimistic), ncol = length(G1), byrow = TRUE)
 colnames(optimistic) <- G1
 
 #Use to check the logic for correlation with scenarios optimistic (red), pessimistic (blue)
 # Scenarios interpretation:
 #  pessimistic = cost of living scenario (high externally driven (energy, commodity etc) inflation, low growth)
 #  baseline = standstill (low growth, low inflation)
 #  optimistic = economic growth accelerates (moderate growth, moderate inflation)
 #  due to this HIC baseline is lowest of 3 scenarios!
 
 for(var_nm in colnames(baseline)){
 df <- data.frame(dates = time(baseline[,var_nm]), baseline = baseline[,var_nm],
 pessimistic = pessimistic[, var_nm], optimistic = optimistic[, var_nm])
 df <- na.omit(df)
 ggplot(data = df, mapping = aes(x = dates, y = optimistic, color = "Optimistic")) +
 geom_line(lwd = 1, lty = 2) + xlab("Year") + ylab(var_nm) +
 geom_line(aes(y = pessimistic, color = "Pessimistic"), lwd = 1, lty = 2) +
 geom_line(aes(y = baseline, color = "Baseline"), lwd = 1) + 
 scale_color_manual(values = col.custom) +
 geom_hline(yintercept = 0, lty = 1, lwd = 0.5) -> pic
 print(pic)
 }
 
 # If no mistakes, save raw data image and go to the next step on trend-cycle decomposition
 cat("Saving results! \n")
 
 objects_to_save <- c("CN", "path", "results_path", "data", "odf_z", "baseline",
 "optimistic", "pessimistic", 
 "data_pessimistic", "data_optimistic",
 "CC_macro_cycles", "problems", "sig_cycles", 
 "CC_macro_suite", "smoothed_data",
 "all_scaled_data")
 setwd(results_path) 
 output_file_name <- paste("CC_macro_cycles_", CN,".Rdata", sep = "")
 save(list = objects_to_save, 
 file = output_file_name)
 setwd(path)
 
 
 cat("Data preparation completed, proceed to modelling step...\n")
 round(Sys.time() - t3, 1) %>% print()


 # Description:   Combine a single country dataset

 # Attributes:
 #   country         -- "EE", "LV" or "LT"
 #   scenario_name   -- "baseline", "optimistic" or "pessimistic"
 
 # Combine data together
 get_country_data <- function(country = 'LT', scenario_name = "baseline"){
 macro_scenarios %>%
 filter(Region   == "EU",
 MacroVar == "EURIBOR_6M",
 Scenario == scenario_name) %>%
 select(ProjectedValue) %>% 
 ts(frequency = 4, start = start_year_quarter) -> euribor_forecast
 
 ts(c(euribor_history, window(euribor_forecast, tsp(euribor_history)[2] + 0.25)), 
 start = tsp(euribor_history)[1], frequency = FREQ) -> euribor
 
 
 # Select baseline scenarios data for 6 variables
 macro_scenarios %>%
 filter(Region == country, Scenario == scenario_name) %>%
 cast(Period ~ MacroVar, value = "ProjectedValue") %>%
 mutate(Period = as.integer(substring(Period, 2))) %>%
 arrange(Period) %>%
 select(-Period) %>%
 log() %>%
 mutate(UR = exp(UR)) %>%
 ts(start = start_year_quarter, frequency = FREQ) -> macro_projected
 names_projected <- c("YER0", "HIC", "HPI", "URX0", "WAG")
 colnames(macro_projected) <- names_projected
 
 # Extract quarterly data:
 macro_history_selected %>%
 filter(Region == country, Frequency == "Q") %>%
 {min(.$DateKey)} -> min_date
 
 MACRO_START <- c(year(min_date), quarter(min_date))
 
 macro_history_selected %>% 
 filter(Region == country, Frequency == "Q") %>%
 cast(DateKey ~ MacroVar, value = "ObservedValue") %>%
 select(-DateKey) %>%
 log() %>%
 mutate(UR_1574 = exp(UR_1574),
 JVAC    = exp(JVAC)) %>%
 ts(start = MACRO_START, frequency = FREQ) -> macro_q
 colnames(macro_q) <- c("YER0", "YEN0", "CGN", "CHN", "HPI", "JVN", "CFN", "URX0", "WAG")
 
 # Correct wage definition via proportional splicing for Q1 2019 wage definition change
 if(country == 'LT'){
 wage <- macro_q[, 'WAG']
 d <- c(window(wage, 2019, 2019)) - c(window(wage, 2018.75, 2018.75)) - 0.022  # value to add
 wage_adj <- c(window(wage, end = 2018.75) + d, window(wage, 2019)) %>%
 ts(start = tsp(wage)[1], frequency = tsp(wage)[3])
 macro_q[,'WAG'] <- wage_adj
 }
 
 # Extract monthly data
 macro_history_selected %>% 
 filter(Region == country, MacroVar == "HICP_2015") %>%
 select(ObservedValue) %>% 
 ts(frequency = 12, start = c(1996, 1)) %>%
 log() %>%
 ta(conversion = "average", to = "quarterly") -> hicp
 
 # Take Eurostat data:
 eurostat_data <- macro_data[[country]]
 
 names <- c(colnames(macro_q), colnames(eurostat_data), "HIC", "E6M")
 
 macro_q <- ts.union(macro_q, eurostat_data, hicp, euribor)
 colnames(macro_q) <- names
 
 # Add baseline projections for G1 variables:
 names_projected %>%
 lapply(function(x){
 y <- na.trim(macro_q[,x], sides = "right")
 ts(c(y, window(macro_projected[,x], tsp(y)[2] + 0.25)), 
 start = tsp(y)[1], frequency = FREQ)
 }) -> extended_data
 
 # Make real variables:
 
 macro_q %<>%
 as_tibble()  %>%
 mutate(YER0 = extended_data[[1]],
 URX0 = extended_data[[4]] / 100,
 WAG = extended_data[[5]],
 CFR = CFN - YED,
 CGR = CGN - YED,
 CHR = CHN - YED,
 WIR = WIN - YED,
 HIC = extended_data[[2]] - log(100),  # fit to units in non-log scale
 HPI = extended_data[[3]] - log(100),
 URX = URX / 100,
 E6M = E6M / 100) %>%
 mutate(WAR = WAG - HIC) %>%  # Real wage, HICP deflated
 mutate(d4HIC = c(rep(NA, FREQ), diff(HIC, 4, 1)),  # add acceleration cycles
 d4HPI = c(rep(NA, FREQ), diff(HPI, 4, 1)),
 d4WAG = c(rep(NA, FREQ), diff(WAG, 4, 1)),
 d4WAR = c(rep(NA, FREQ), diff(WAR, 4, 1))) %>%
 ts(frequency = 4, start = c(1994, 1))
 
 # YER seasonally adjusted should grow as YER0 projected after seasonal adjustment
 
 macro_q[,"YER0"] -> y0
 x <- stl(trim_and_smooth(y0), s.window = 7)
 y_sa <- x$time.series[, "trend"] + x$time.series[, "remainder"]
 macro_q[,"YER"] -> y
 y_adj <- trim_and_smooth(y)
 
 z0 <- window(y_sa/stats::lag(y_sa, -1), start = tsp(y_adj)[2] + 0.25)
 extension <- window(y, start = end(y_adj))
 
 for(index in 1:right_na(y)){
 extension[1 + index] <- extension[index] * z0[index]
 }
 y_final <- ts(c(y_adj, extension[-1]), start = start(y_adj), frequency = frequency(y_adj))
 
 macro_q[,"URX0"] -> u0
 x <- stl(trim_and_smooth(u0), s.window = 7)
 u_sa <- x$time.series[, "trend"] + x$time.series[, "remainder"]
 macro_q[,"URX"] -> u
 u_adj <- trim_and_smooth(u)
 
 z0 <- window(u_sa/stats::lag(u_sa, -1), start = tsp(u_adj)[2] + 0.25)
 extension <- window(u, start = end(u_adj))
 for(index in 1:right_na(u)){
 extension[1 + index] <- extension[index] * z0[index]
 }
 u_final <- ts(c(u_adj, extension[-1]), start = start(u_adj), frequency = frequency(u_adj))
 
 extended_data <- ts.union(y, y_final, u_final)
 
 macro_q %<>%
 as_tibble() %>%
 mutate(YER = extended_data[, 2],
 URX = extended_data[, 3]) %>%
 ts(frequency = 4, start = c(1994, 1))
 
 return(macro_q)
 }

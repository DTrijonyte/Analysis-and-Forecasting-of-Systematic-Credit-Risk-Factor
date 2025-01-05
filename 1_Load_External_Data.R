# Description:      Load external data and merge for the use with EE, LT, LV


# Clear everything out so we're starting clean
rm(list = ls())

# Get required libraries
library(pacman)
p_load(yfR, forecast, readxl, lubridate, eurostat,
       dplyr, magrittr, reshape, tempdisagg, tidyr, 
       zoo, openxlsx, readr, install = TRUE)

# Process OMX data from Yahoo Finance (year 2020 - ongoing)
omx_yahoo <- yf_get(tickers = c("^OMXVGI", "^OMXRGI", "^OMXTGI"), 
                    first_date = as.Date("2020-01-01")) %>%
  mutate(date = as.Date(ref_date),
         values = price_adjusted, 
         geo = case_when(
           ticker == "^OMXVGI" ~ "LT",
           ticker == "^OMXRGI" ~ "LV",
           ticker == "^OMXTGI" ~ "EE"
         )) %>%
  complete(date = seq(min(date), max(date), by = "day"), nesting(geo)) %>%
  group_by(geo) %>%
  mutate(values = na.interp(values)) %>%
  mutate(time = as.yearqtr(date)) %>%
  group_by(geo, time) %>%
  summarise(values = mean(log(values))) %>%
  ungroup()

# Process OMX data from Nasdaq (year 2000-2020)
omx_nasdaq <- read_excel("charts_20000101_20200101.xlsx")

omx_nasdaq <- omx_nasdaq %>%
  rename_with(.cols = 1, ~"geo") %>%
  mutate(date = as.Date(Date),
         values = Value,
         geo = case_when(geo == "OMXVGI" ~ "LT",
                         geo == "OMXRGI" ~ "LV",
                         geo == "OMXTGI" ~ "EE")) %>% 
  complete(date = seq(min(date), max(date), by = "day"), nesting(geo)) %>%
  group_by(geo) %>%
  mutate(values = na.interp(values)) %>%
  mutate(time = as.yearqtr(date)) %>%
  group_by(geo, time) %>%
  summarise(values = mean(log(values))) %>%
  ungroup()

# Combine the datasets
omx_data <- bind_rows(omx_nasdaq, omx_yahoo) %>% 
  group_by(geo, time) %>%
  summarize(values = mean(values), .groups = 'drop')


# Get all Eurostat tables and investigate which to use: 
# Search for needed table names: 
search_eurostat("A*")  %>% 
  filter(grepl("Q", data.start)) %>%
  distinct(.keep_all = TRUE) -> sectoral_tables 
search_eurostat("GDP") %>% 
  filter(grepl("Q", data.start)) %>%
  distinct(.keep_all = TRUE) -> gdp_tables 
search_eurostat("Unemplo") %>% 
  filter(grepl("Q", data.start)) %>% 
  distinct(.keep_all = TRUE) -> un_tables 
search_eurostat("Sentim") %>% 
  distinct(.keep_all = TRUE) -> ci_tables 

# Set filtering parameters: 
GEO      <- c('EE', 'LT', 'LV')
REAL_NA_ITEMS <- c( 'B1GQ',    # Gross domestic product 
                    'P3_S13',  # General Government consumption 
                    'P31_S14', # Households consumption
                    'P51G',    # Gross Fixed Capital Formation 
                    'P6',      # Exports of goods and services 
                    'P7'       # Imports  of goods and services 
) 
NOMINAL_NA_ITEMS <- c( 
  'B1GQ',    # Gross domestic product
  'D1'       # Compensation of employees with social contributions 
) 

gdp_data <- get_eurostat('namq_10_gdp') # Bulk load always works 
gdp_data %>%
  filter(unit  == "CP_MEUR",
         s_adj == "SCA",   # seasonally and calendar adjusted
         na_item %in% NOMINAL_NA_ITEMS, 
         geo %in% GEO) %>% 
  transmute(geo, na_item, time = TIME_PERIOD, values) %>% 
  arrange(geo, na_item, time) %>%
  mutate(values = log(values)) -> nominal_data
# Nominal wages should be deflated with HICP, nominal GDP is used for GDP deflator definition 

gdp_data %<>% 
  filter(unit  == "CLV15_MEUR", 
         s_adj == "SCA",   # seasonally and calendar adjusted 
         na_item %in% REAL_NA_ITEMS, 
         geo   %in% GEO) %>% 
  transmute(geo, na_item, time = TIME_PERIOD, values) %>% 
  arrange(geo, na_item, time) %>% 
  mutate(values = log(values)) 

# Add unemployment data: 
un_data <- get_eurostat('une_rt_q', time = "date")

un_data %>% 
  filter(age   == "Y15-74",
         s_adj == "SA", 
         unit  == "PC_ACT", 
         sex   == "T",
         geo %in% GEO) %>% 
  transmute(geo, time = TIME_PERIOD, values) %>% 
  arrange(geo, time) -> un_data_current 

un_data_h <- get_eurostat('une_rt_q_h', time = "date") 

un_data_h %>% 
  filter(age   == "Y15-74",
         s_adj == "SA",
         unit  == "PC_ACT", 
         sex   == "T", 
         geo %in% GEO) %>% 
  transmute(geo, time = TIME_PERIOD, values) %>% 
  arrange(geo, time) -> un_data_historical

un_data_historical %>% 
  filter(time < "2009-01-01") %>% 
  {rbind(., un_data_current)} %>% 
  arrange(geo, time) -> un_data_merged 

# Add confidence indicators: 
conf_data <- get_eurostat('ei_bssi_m_r2', time = "date")

keep_middle <- function(x, split = "-"){strsplit(x, "-") %>% unlist() %>%.[2]} 

conf_data %<>% 
  filter(s_adj == "SA",
         geo %in% GEO) %>% 
  transmute(geo, indic = as.factor(indic), time = TIME_PERIOD, values) %>%
  arrange(geo, indic, time) 
levels(conf_data$indic) %<>% lapply(keep_middle) %>% unlist() 
conf_data %<>% mutate(indic = as.character(indic)) 

# Description (we keep only middle): 
# BS-CCI-BAL construction confidence indicator 
# BS-ESI-I economic sentiment indicator (average of other indicators) 
# BS-ICI-BAL industrial confidence indicator 
# BS-RCI-BAL retail confidence indicator
# BS-CSMSCI-BAL consumer confidence indicator (CSM to distinguish with C) 
# BS-SCI-BAL services confidence indicator 


# Data merging function: 
get_country_data <- function(country = "EE"){ 
  
  # GDP by expenditure approach with STL 
  gdp_data %>%
    filter(geo == country) %>% 
    cast(time ~ na_item, value = "values") %>% 
    select(-time) %>% 
    ts(start = c(1995, 1), frequency = 4) -> real_q 
  # nominal GDP, wages 
  nominal_data %>% 
    filter(geo == country) %>% 
    cast(time ~ na_item, value = "values") %>%
    select(-time) %>% 
    ts(start = c(1995, 1), frequency = 4) -> nominal_q 
  
  # unemployment rate 
  un_data_merged %>% 
    filter(geo == country) %>% 
    select(values) %>% 
    ts(start = c(1997, 1), frequency = 4) -> unempl 
  
  # confidence indicators
  conf_data %>%
    filter(geo == country, time > "1994-12-01") %>%
    cast(time ~ indic, value = "values") %>%
    select(-time) %>% ts(start = c(1995, 1), frequency = 12) %>%
    ta(conversion = "average", to = "quarterly") -> conf_ind
  
  # # OMX data
  omx_data %>%
    filter(geo == country) %>%
    select(values) %>%
    ts(start = c(2000, 1), frequency = 4) -> omx
  
  yed <- nominal_q[,"B1GQ"] - real_q[,"B1GQ"] # since after log transformation 
  
  # Combine data
  data_names <- c("YER", "GCR", "PCR", "ITR", "XTR", "MTR", "YEN", "WIN", "URX",
                  "FCI", "CCI", "ESI", "ICI", "RCI", "SCI", "OMX", "YED")
  macro_q <- ts.union(real_q, nominal_q, unempl, conf_ind, omx, yed)
  colnames(macro_q) <- data_names
  
  return(macro_q)
} 

# Store results in a list of country data tables
GEO %>% sapply(get_country_data, USE.NAMES = TRUE, simplify = FALSE) -> macro_data 

DO_PLOTS = T
if(DO_PLOTS){ 
  macro_data[["LT"]][,1:10] %>% plot(las = 1, cex.lab =.7, main = "Lithuania 1") 
  macro_data[["LT"]][,11:17] %>% plot(las = 1, cex.lab =.7, main = "Lithuania 2")
  macro_data[["LV"]][,1:10] %>% plot(las = 1, cex.lab =.7, main = "Latvia 1") 
  macro_data[["LV"]][,11:17] %>% plot(las = 1, cex.lab =.7, main = "Latvia 2") 
  macro_data[["EE"]][,1:10] %>% plot(las = 1, cex.lab =.7, main = "Estonia 1") 
  macro_data[["EE"]][,11:17] %>% plot(las = 1, cex.lab =.7, main = "Estonia 2") 
}

save(macro_data, file = "macro_data.Rdata")

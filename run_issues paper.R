
## load library
source("./load library.R")

## Step 1. Import raw dataset
sim_results <- readxl::read_excel(path = "./rawdata - csv/SEM paper_Table 3_4_5.xlsx", sheet = "Sheet1", col_names = TRUE)

## Step 2. Preprocess

## -- remove space in column names
names(sim_results) <- names(sim_results) %>%
  stringr::str_replace_all(pattern = "\\s+", replacement = "\\_")

## Step 3. create analysis dataset
## Nonconvergence cases
nonconvergence <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%  
  dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_nonconvergence)

## Heywood cases
heywood <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>% 
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_heywood)

## Quality of Recovering the True factor scores
fIndeterminacy <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%   
  dplyr::summarise(avg_fIndeterminacy = round(x = mean(Table_5_Factor_Indeterminacy), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_fIndeterminacy)

## -- Save analysis datasets
readr::write_csv(x = nonconvergence, path = "./analysis datasets/nonconvergence_avg.csv")
readr::write_csv(x = heywood, path = "./analysis datasets/heywood_avg.csv")
readr::write_csv(x = fIndeterminacy, path = "./analysis datasets/fIndeterminacy_avg.csv")

## Step 4. Plot simulation results





## load library
source("./load library.R")

##################################
## Step 1. Import raw dataset.  ##
##################################
sim_results <- readxl::read_excel(path = "./rawdata - csv/SEM paper_Table 3_4_5.xlsx", sheet = "Sheet1", col_names = TRUE)

##################################
## Step 2. Preprocess.          ##
##################################
## remove space in column names
names(sim_results) <- names(sim_results) %>%
  stringr::str_replace_all(pattern = "\\s+", replacement = "\\_")

#######################################
## Step 3. create analysis dataset.  ##
#######################################
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

#######################################
## Step 4. Plot simulation results.  ##
#######################################

## Read analysis dataset
nonconvergence <- readr::read_csv(file = "./analysis datasets/nonconvergence_avg.csv", col_names = T)
heywood <- readr::read_csv(file = "./analysis datasets/heywood_avg.csv", col_names = T)
fIndeterminacy <- readr::read_csv(file = "./analysis datasets/fIndeterminacy_avg.csv", col_names = T)

## ---------------------------- ##
## Plot Nonconvergence cases.   ##
## ---------------------------- ##
## adjust variable type
nonconvergence$Type_of_Indicators <- as.factor(nonconvergence$Type_of_Indicators)
nonconvergence$Number_of_Indicators <- as.factor(nonconvergence$Number_of_Indicators)

plt <- ggplot2::ggplot(data = nonconvergence, aes(x = Type_of_Indicators, y = avg_nonconvergence, group = Number_of_Indicators,
                                           color = Number_of_Indicators,
                                           shape = Number_of_Indicators))

plt_interaction <- plt +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Nonconvergence (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.93, 0.88),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

plt_interaction

## Save plot
ggplot2::ggsave(filename = "./plots/plt_nonconvergence.png", plot = last_plot(),
                width = 30, height = 20, units = "cm", dpi = 300)


## ------------------------- ##
## Plot Heywood cases.       ##
## ------------------------- ##
## adjust variable type
heywood$Type_of_Indicators <- as.factor(heywood$Type_of_Indicators)
heywood$Number_of_Indicators <- as.factor(heywood$Number_of_Indicators)

plt <- ggplot2::ggplot(data = heywood, aes(x = Type_of_Indicators, y = avg_heywood, group = Number_of_Indicators,
                                           color = Number_of_Indicators,
                                           shape = Number_of_Indicators))

plt_interaction <- plt +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Heywood (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.93, 0.88),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

plt_interaction

## Save plot
ggplot2::ggsave(filename = "./plots/plt_heywood.png", plot = last_plot(),
                width = 30, height = 20, units = "cm", dpi = 300)


## ----------------------------- ##
## Plot Factor Indeterminacy.    ##
## ----------------------------- ##
## adjust variable type
fIndeterminacy$Type_of_Indicators <- as.factor(fIndeterminacy$Type_of_Indicators)
fIndeterminacy$Number_of_Indicators <- as.factor(fIndeterminacy$Number_of_Indicators)

plt <- ggplot2::ggplot(data = fIndeterminacy, aes(x = Type_of_Indicators, y = avg_fIndeterminacy, group = Number_of_Indicators,
                                                  color = Number_of_Indicators,
                                                  shape = Number_of_Indicators))

plt_interaction <- plt +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Quality of Recovering Factor Scores") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.93, 0.23),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

plt_interaction

## Save plot
ggplot2::ggsave(filename = "./plots/plt_fIndeterminacy.png", plot = last_plot(),
                width = 30, height = 20, units = "cm", dpi = 300)





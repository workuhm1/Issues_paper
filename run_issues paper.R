
## load library
source("./load library.R")

##################################
## Step 1. Import raw dataset.  ##
##################################
# sim_results <- readxl::read_excel(path = "./rawdata - csv/SEM paper_Table 3_4_5.xlsx", sheet = "Sheet1", col_names = TRUE)
sim_results <- readxl::read_excel(path = "./rawdata - csv/SEM paper_Table 3_4_5_new.xlsx", sheet = "Sheet1", col_names = TRUE)

##################################
## Step 2. Preprocess.          ##
##################################
## remove space in column names
names(sim_results) <- names(sim_results) %>%
  stringr::str_replace_all(pattern = "\\s+", replacement = "\\_")

## convert type of the design factors
sim_results$Type_of_Indicators <- factor(sim_results$Type_of_Indicators, levels = c("BLR", "BMR", "Interval"))
sim_results$Communality <- factor(sim_results$Communality, levels = c("Weak", "Moderate", "Strong"))
sim_results$Correlation_between_Factors <- factor(sim_results$Correlation_between_Factors, levels = c("Independence", "Moderate", "Strong"))
sim_results$Sample_size <- factor(sim_results$Sample_size, levels = c("50", "100", "300", "3000"))
sim_results$Number_of_Indicators <- factor(sim_results$Number_of_Indicators, levels = c("6", "10", "16"))


#######################################
## Step 3. create analysis dataset.  ##
#######################################

## ---------------------------- ##
## -- Nonconvergence cases      ##
## ---------------------------- ##
# nonconvergence <- sim_results %>% 
#   dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%  
#   dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
#   dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_nonconvergence)

## Nonconvergence: interaction between Type of Indicators and Number of Indicators
nonconvg_ab <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%  
  dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_nonconvergence)

## Nonconvergence: interaction between Type of Indicators and Factor structure
nonconvg_ac <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Communality) %>%  
  dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Communality, avg_nonconvergence)

## Nonconvergence: interaction between Type of Indicators and Sample size
nonconvg_ae <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Sample_size) %>%  
  dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Sample_size, avg_nonconvergence)

## Nonconvergence: interaction between Type of Indicators and Sample size
nonconvg_be <- sim_results %>% 
  dplyr::group_by(Number_of_Indicators, Sample_size) %>%  
  dplyr::summarise(avg_nonconvergence = round(x = mean(Table_3_Nonconvergence), digits = 2)) %>% 
  dplyr::select(Number_of_Indicators, Sample_size, avg_nonconvergence)

## -- combine the tables to get max/min for plotting later
nonconvg_all <- list(nonconvg_ab, nonconvg_ac, nonconvg_ae, nonconvg_be)
max_nonconvg <- max(unlist(nonconvg_all))
min_nonconvg <- min(unlist(nonconvg_all))

## ------------------------- ##
## Plot nonconvergence.      ##
## ------------------------- ##
# tikz(file = "./plots/plt_nonconvergence_largeEffects.tex", width = 5, height = 5)

p1 <- ggplot2::ggplot(data = nonconvg_ab, 
                      aes(x = Type_of_Indicators, y = avg_nonconvergence, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_nonconvg, max_nonconvg)) + 
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Nonconvergence (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    # legend.position = c(0.93, 0.88),
    legend.position = c(0.86, 0.78),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p1

p2 <- ggplot2::ggplot(data = nonconvg_ac, 
                      aes(x = Type_of_Indicators, y = avg_nonconvergence, 
                          group = Communality, 
                          color = Communality, shape = Communality)) +
  geom_line(aes(linetype = Communality), size = 0.75) +
  geom_point(aes(shape = Communality), size = 2) +
  scale_y_continuous(limits = c(min_nonconvg, max_nonconvg)) + 
  guides(linetype = guide_legend("Factor structure"),
         colour = guide_legend("Factor structure"),
         shape = guide_legend("Factor structure")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Nonconvergence (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.86, 0.78),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p2

p3 <- ggplot2::ggplot(data = nonconvg_ae, 
                      aes(x = Type_of_Indicators, y = avg_nonconvergence, 
                          group = Sample_size, 
                          color = Sample_size, shape = Sample_size)) +
  geom_line(aes(linetype = Sample_size), size = 0.75) +
  geom_point(aes(shape = Sample_size), size = 2) +
  scale_y_continuous(limits = c(min_nonconvg, max_nonconvg)) + 
  guides(linetype = guide_legend("Sample size"),
         colour = guide_legend("Sample size"),
         shape = guide_legend("Sample size")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Nonconvergence (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.86, 0.78),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p3

p4 <- ggplot2::ggplot(data = nonconvg_be, 
                      aes(x = Sample_size, y = avg_nonconvergence, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_nonconvg, max_nonconvg)) + 
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Sample size",
       y = "Prevalence of Nonconvergence (%)") +
  theme_classic() +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.86, 0.78),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p4

## combine plots
g_nonconvg <- gridExtra::arrangeGrob(p1, p2, p3, p4)
# gridExtra::arrangeGrob(p1, p2, p3, p4)

## Save plot
# ggplot2::ggsave(filename = "./plots/plt_nonconvergence_all interactions.png", 
#                 plot = last_plot(),
#                 width = 30, height = 20, units = "cm", dpi = 300)

ggplot2::ggsave(filename = "./plots/plt_nonconvergence_largeEffects.png",
# ggplot2::ggsave(filename = "./plots/plt_nonconvergence_largeEffects.eps",
                plot = g_nonconvg,
                # width = 30, height = 20, units = "cm", 
                width = 20, height = 10, units = "cm",
                dpi = 300
                )

# print(g_nonconvg)

dev.off()


## ------------------- ##
## Heywood cases.      ##
## ------------------- ##
# heywood <- sim_results %>% 
#   dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>% 
#   dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
#   dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_heywood)

## Heywood: interaction between Type of Indicators and Number of Indicators
heywood_ab <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_heywood)

## Heywood: interaction between Type of Indicators and Factor structure
heywood_ac <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Communality) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Communality, avg_heywood)

## Heywood: interaction between Type of Indicators and Correlation btn factors
heywood_ad <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Correlation_between_Factors) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Correlation_between_Factors, avg_heywood)

## Heywood: interaction between Type of Indicators and Sample size
heywood_ae <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Sample_size) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Sample_size, avg_heywood)

## Heywood: interaction between Number of Indicators and Factor structure
heywood_bc <- sim_results %>% 
  dplyr::group_by(Number_of_Indicators, Communality) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Number_of_Indicators, Communality, avg_heywood)

## Heywood: interaction between Number of Indicators and Sample size
heywood_be <- sim_results %>% 
  dplyr::group_by(Number_of_Indicators, Sample_size) %>%  
  dplyr::summarise(avg_heywood = round(x = mean(Table_4_Heywood), digits = 2)) %>% 
  dplyr::select(Number_of_Indicators, Sample_size, avg_heywood)

## -- combine the tables to get max/min for plotting later
heywood_all <- list(heywood_ab, heywood_ac, heywood_ad, heywood_ae,
                    heywood_bc, heywood_be)
max_heywood <- max(unlist(heywood_all))
min_heywood <- min(unlist(heywood_all))

## ------------------------- ##
## Plot nonconvergence.      ##
## ------------------------- ##
p1 <- ggplot2::ggplot(data = heywood_ab, 
                      aes(x = Type_of_Indicators, y = avg_heywood, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Heywood (%)") +
  theme_classic(base_size = 7) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    # legend.position = c(0.86, 0.70),
    legend.position = c(0.92, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p1

p2 <- ggplot2::ggplot(data = heywood_ac, 
                      aes(x = Type_of_Indicators, y = avg_heywood, 
                          group = Communality, 
                          color = Communality, shape = Communality)) +
  geom_line(aes(linetype = Communality), size = 0.75) +
  geom_point(aes(shape = Communality), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("Factor structure"),
         colour = guide_legend("Factor structure"),
         shape = guide_legend("Factor structure")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Heywood (%)") +
  theme_classic(base_size = 7) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.91, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p2

p3 <- ggplot2::ggplot(data = heywood_ad, 
                      aes(x = Type_of_Indicators, y = avg_heywood, 
                          group = Correlation_between_Factors, 
                          color = Correlation_between_Factors, shape = Correlation_between_Factors)) +
  geom_line(aes(linetype = Correlation_between_Factors), size = 0.75) +
  geom_point(aes(shape = Correlation_between_Factors), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("Corr. Factors"),
         colour = guide_legend("Corr. Factors"),
         shape = guide_legend("Corr. Factors")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Heywood (%)") +
  theme_classic(base_size = 7) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.91, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p3

p4 <- ggplot2::ggplot(data = heywood_ae, 
                      aes(x = Type_of_Indicators, y = avg_heywood, 
                          group = Sample_size, 
                          color = Sample_size, shape = Sample_size)) +
  geom_line(aes(linetype = Sample_size), size = 0.75) +
  geom_point(aes(shape = Sample_size), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("Sample size"),
         colour = guide_legend("Sample size"),
         shape = guide_legend("Sample size")) +
  labs(x = "Type of Indicators",
       y = "Prevalence of Heywood (%)") +
  # theme_classic() +                                ## to put graph in a box
  theme_classic(base_size = 7) +                     ## to resize legend text size
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.92, 0.70),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p4

p5 <- ggplot2::ggplot(data = heywood_bc, 
                      aes(x = Communality, y = avg_heywood, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Factor structure",
       y = "Prevalence of Heywood (%)") +
  theme_classic(base_size = 7) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.92, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p5

p6 <- ggplot2::ggplot(data = heywood_be, 
                      aes(x = Sample_size, y = avg_heywood, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_heywood, max_heywood)) + 
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Sample size",
       y = "Prevalence of Heywood (%)") +
  theme_classic(base_size = 7) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.92, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# p6

## combine plots
g_heywood <- gridExtra::arrangeGrob(p1, p2, p3, p4, p5, p6)

## Save plot
# ggplot2::ggsave(filename = "./plots/plt_nonconvergence_all interactions.png", 
#                 plot = last_plot(),
#                 width = 30, height = 20, units = "cm", dpi = 300)

ggplot2::ggsave(filename = "./plots/plt_heywood_largeEffects.png", 
                plot = g_heywood,
                width = 30, height = 20, units = "cm", dpi = 300)

dev.off()

## -------------------------- ##
## Factor Indeterminacy.      ##
## -------------------------- ##
# fIndeterminacy <- sim_results %>% 
#   dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%   
#   dplyr::summarise(avg_fIndeterminacy = round(x = mean(Table_5_Factor_Indeterminacy), digits = 2)) %>% 
#   dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_fIndeterminacy)

## Factor Indeterminacy: interaction between Type of Indicators and Number of Indicators
fIndeterminacy_ab <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Number_of_Indicators) %>%  
  dplyr::summarise(avg_fIndeterminacy = round(x = mean(Table_5_Factor_Indeterminacy), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Number_of_Indicators, avg_fIndeterminacy)

## Factor Indeterminacy: interaction between Type of Indicators and Factor structure
fIndeterminacy_ac <- sim_results %>% 
  dplyr::group_by(Type_of_Indicators, Communality) %>%  
  dplyr::summarise(avg_fIndeterminacy = round(x = mean(Table_5_Factor_Indeterminacy), digits = 2)) %>% 
  dplyr::select(Type_of_Indicators, Communality, avg_fIndeterminacy)


## -- combine the tables to get max/min for plotting later
fIndeterminacy_all <- list(fIndeterminacy_ab, fIndeterminacy_ac)

# max_fIndeterminacy <- max(unlist(fIndeterminacy_all))
# min_fIndeterminacy <- min(unlist(fIndeterminacy_all))

max_fIndeterminacy <- max(max(fIndeterminacy_ab$avg_fIndeterminacy),
                          max(fIndeterminacy_ac$avg_fIndeterminacy))

## min_fIndeterminacy = 0.42
min_fIndeterminacy <- min(min(fIndeterminacy_ab$avg_fIndeterminacy),
                          min(fIndeterminacy_ac$avg_fIndeterminacy))


## ------------------------- ##
## Plot nonconvergence.      ##
## ------------------------- ##
p1 <- ggplot2::ggplot(data = fIndeterminacy_ab, 
                      aes(x = Type_of_Indicators, y = avg_fIndeterminacy, 
                          group = Number_of_Indicators, 
                          color = Number_of_Indicators, shape = Number_of_Indicators)) +
  geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
  geom_point(aes(shape = Number_of_Indicators), size = 2) +
  scale_y_continuous(limits = c(min_fIndeterminacy, max_fIndeterminacy)) +
  guides(linetype = guide_legend("No. Indicators"),
         colour = guide_legend("No. Indicators"),
         shape = guide_legend("No. Indicators")) +
  labs(x = "Type of Indicators",
       y = "Quality of Recovering Factor Scores") +
  theme_classic(base_size = 8) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.91, 0.25),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p1

p2 <- ggplot2::ggplot(data = fIndeterminacy_ac, 
                      aes(x = Type_of_Indicators, y = avg_fIndeterminacy, 
                          group = Communality, 
                          color = Communality, shape = Communality)) +
  geom_line(aes(linetype = Communality), size = 0.75) +
  geom_point(aes(shape = Communality), size = 2) +
  scale_y_continuous(limits = c(min_fIndeterminacy, max_fIndeterminacy)) +
  guides(linetype = guide_legend("Factor structure"),
         colour = guide_legend("Factor structure"),
         shape = guide_legend("Factor structure")) +
  labs(x = "Type of Indicators",
       y = "Quality of Recovering Factor Scores") +
  theme_classic(base_size = 8) +                                ## to put graph in a box
  theme(
    panel.background = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.position = c(0.91, 0.25),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )
# p2

## combine plots
g_fIndeterminacy <- gridExtra::arrangeGrob(p1, p2, ncol = 2)

## Save plot
ggplot2::ggsave(filename = "./plots/plt_QualityRecooveryFS_largeEffects.png", 
                plot = g_fIndeterminacy,
                width = 30, height = 20, units = "cm", dpi = 300)



# ## -- Save analysis datasets
# readr::write_csv(x = nonconvergence, path = "./analysis datasets/nonconvergence_avg.csv")
# readr::write_csv(x = heywood, path = "./analysis datasets/heywood_avg.csv")
# readr::write_csv(x = fIndeterminacy, path = "./analysis datasets/fIndeterminacy_avg.csv")
# 
# #######################################
# ## Step 4. Plot simulation results.  ##
# #######################################
# 
# ## Read analysis dataset
# nonconvergence <- readr::read_csv(file = "./analysis datasets/nonconvergence_avg.csv", col_names = T)
# heywood <- readr::read_csv(file = "./analysis datasets/heywood_avg.csv", col_names = T)
# fIndeterminacy <- readr::read_csv(file = "./analysis datasets/fIndeterminacy_avg.csv", col_names = T)
# 
# ## ---------------------------- ##
# ## Plot Nonconvergence cases.   ##
# ## ---------------------------- ##
# ## adjust variable type
# nonconvergence$Type_of_Indicators <- as.factor(nonconvergence$Type_of_Indicators)
# nonconvergence$Number_of_Indicators <- as.factor(nonconvergence$Number_of_Indicators)
# 
# plt <- ggplot2::ggplot(data = nonconvergence, aes(x = Type_of_Indicators, y = avg_nonconvergence, group = Number_of_Indicators,
#                                            color = Number_of_Indicators,
#                                            shape = Number_of_Indicators))
# 
# plt_interaction <- plt +
#   geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
#   geom_point(aes(shape = Number_of_Indicators), size = 2) +
#   guides(linetype = guide_legend("No. Indicators"),
#          colour = guide_legend("No. Indicators"),
#          shape = guide_legend("No. Indicators")) +
#   labs(x = "Type of Indicators",
#        y = "Prevalence of Nonconvergence (%)") +
#   theme_classic() +                                ## to put graph in a box
#   theme(
#     panel.background = element_rect(colour = "black", fill = NA),
#     axis.text = element_text(colour = 1, size = 12),
#     legend.position = c(0.93, 0.88),
#     legend.background = element_blank(),
#     legend.box.background = element_rect(colour = "black")
#   )
# 
# plt_interaction
# 
# ## Save plot
# ggplot2::ggsave(filename = "./plots/plt_nonconvergence.png", plot = last_plot(),
#                 width = 30, height = 20, units = "cm", dpi = 300)
# 
# 
# ## ------------------------- ##
# ## Plot Heywood cases.       ##
# ## ------------------------- ##
# ## adjust variable type
# heywood$Type_of_Indicators <- as.factor(heywood$Type_of_Indicators)
# heywood$Number_of_Indicators <- as.factor(heywood$Number_of_Indicators)
# 
# plt <- ggplot2::ggplot(data = heywood, aes(x = Type_of_Indicators, y = avg_heywood, group = Number_of_Indicators,
#                                            color = Number_of_Indicators,
#                                            shape = Number_of_Indicators))
# 
# plt_interaction <- plt +
#   geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
#   geom_point(aes(shape = Number_of_Indicators), size = 2) +
#   guides(linetype = guide_legend("No. Indicators"),
#          colour = guide_legend("No. Indicators"),
#          shape = guide_legend("No. Indicators")) +
#   labs(x = "Type of Indicators",
#        y = "Prevalence of Heywood (%)") +
#   theme_classic() +                                ## to put graph in a box
#   theme(
#     panel.background = element_rect(colour = "black", fill = NA),
#     axis.text = element_text(colour = 1, size = 12),
#     legend.position = c(0.93, 0.88),
#     legend.background = element_blank(),
#     legend.box.background = element_rect(colour = "black")
#   )
# 
# plt_interaction
# 
# ## Save plot
# ggplot2::ggsave(filename = "./plots/plt_heywood.png", plot = last_plot(),
#                 width = 30, height = 20, units = "cm", dpi = 300)
# 
# 
# ## ----------------------------- ##
# ## Plot Factor Indeterminacy.    ##
# ## ----------------------------- ##
# ## adjust variable type
# fIndeterminacy$Type_of_Indicators <- as.factor(fIndeterminacy$Type_of_Indicators)
# fIndeterminacy$Number_of_Indicators <- as.factor(fIndeterminacy$Number_of_Indicators)
# 
# plt <- ggplot2::ggplot(data = fIndeterminacy, aes(x = Type_of_Indicators, y = avg_fIndeterminacy, group = Number_of_Indicators,
#                                                   color = Number_of_Indicators,
#                                                   shape = Number_of_Indicators))
# 
# plt_interaction <- plt +
#   geom_line(aes(linetype = Number_of_Indicators), size = 0.75) +
#   geom_point(aes(shape = Number_of_Indicators), size = 2) +
#   guides(linetype = guide_legend("No. Indicators"),
#          colour = guide_legend("No. Indicators"),
#          shape = guide_legend("No. Indicators")) +
#   labs(x = "Type of Indicators",
#        y = "Quality of Recovering Factor Scores") +
#   theme_classic() +                                ## to put graph in a box
#   theme(
#     panel.background = element_rect(colour = "black", fill = NA),
#     axis.text = element_text(colour = 1, size = 12),
#     legend.position = c(0.93, 0.23),
#     legend.background = element_blank(),
#     legend.box.background = element_rect(colour = "black")
#   )
# 
# plt_interaction
# 
# ## Save plot
# ggplot2::ggsave(filename = "./plots/plt_fIndeterminacy.png", plot = last_plot(),
#                 width = 30, height = 20, units = "cm", dpi = 300)





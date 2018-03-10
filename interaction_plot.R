
interaction_plot <- function() {
  
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
}
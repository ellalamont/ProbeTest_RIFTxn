# Correlations between samples (scatterplots)
# E. Lamont
# 7/8/25

source("Import_data.R") # to get All_tpm




# Log10 transform the data
All_tpm_Log10 <- All_tpm %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values



# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

###################################################################
########################### MAKE AVERAGES #########################

# Add columns for averages
All_tpm_Log10 <- All_tpm_Log10 %>% 
  mutate(
    AVERAGE_Broth_Untreated = rowMeans(select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE),
    AVERAGE_Broth_RIF = rowMeans(select(., c(Ra_broth_RIF_1, Ra_broth_RIF_2, Ra_broth_RIF_3, Ra_broth_RIF_4)), na.rm = TRUE),
  )
# NOT FINISHED!!!

###################################################################
####################### AVERAGES GGCORRPLOT #######################

Sample1 <- "AVERAGE_THP1Spiked" # THP1 spiked Captured
Sample2 <- "AVERAGE_BrothNotCaptured" # Broth Not Captured
ScatterCorr <- my_tpm_subset_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 1e6 Ra THP1 spiked captured VS Broth Not captured (Not scaled)",
       x = paste0("Log10(TPM+1) Ra1e6 THP1 samples averaged"), y = paste0("Log10(TPM+1) NOT captured Broth samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_vs_BrothNOTCaptured_Averages.pdf"),
       path = "Correlation_Figures",
       width = 7, height = 5, units = "in")
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_vs_BrothNOTCaptured_Averages.png"),
       path = "Correlation_Figures",
       width = 7, height = 5, units = "in")
ggsave(ScatterCorr,
       file = paste0("THP1Spiked1e6_vs_BrothNOTCaptured_Averages_v2.png"),
       path = "Correlation_Figures",
       width = 5, height = 5, units = "in")




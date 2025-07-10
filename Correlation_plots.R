# Correlations between samples (scatterplots)
# E. Lamont
# 7/8/25

source("Import_data.R") # to get All_tpm


options(scipen = 0)

# Log10 transform the data
All_tpm_Log10 <- All_tpm %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values
All_tpm_Log10_numeric <- All_tpm_Log10 %>% select(-Gene)


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


###########################################################
############### RARITY ALL GRAPHS TOGETHER ################

# Need to remove the gene column or it won't work
All_tpm_numeric <- All_tpm %>% select(-Gene)

# pairs(All_tpm_numeric)

# https://borisleroy.com/en/2013/06/09/correlation-plots-in-r/
# install.packages("Rarity")
library(Rarity)

# Pearson
# pdf("ggcorrplot_Figures/rarity_PearsonLog10_v1.pdf", width = 10, height = 10)
corPlot(All_tpm_Log10_numeric, method = "pearson",
        title = "Pearson Correlation Log10(TPM+1)")
# dev.off()

# Just the THP1 spiked RIF samples because they don't look so similar on the PCA
pdf("Figures/Correlations_Scatterplots/rarity_THP1Spiked.RIF_v1.pdf", width = 10, height = 10)
All_tpm_Log10_numeric_SpikedRIF <- All_tpm_Log10_numeric %>% select(THP1_spiked_Ra_RIF_1, THP1_spiked_Ra_RIF_2, THP1_spiked_Ra_RIF_3, THP1_spiked_Ra_RIF_4, THP1_spiked_Ra_RIF_5)
corPlot(All_tpm_Log10_numeric_SpikedRIF, method = "pearson",
        title = "THP1 spiked RIF Pearson Correlation Log10(TPM+1)")
dev.off()

# Just the THP1 spiked samples to compare
All_tpm_Log10_numeric_SpikedUntreated <- All_tpm_Log10_numeric %>% select(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)
pdf("Figures/Correlations_Scatterplots/rarity_THP1Spiked.Untreated_v1.pdf", width = 10, height = 10)
corPlot(All_tpm_Log10_numeric_SpikedUntreated, method = "pearson",
        title = "THP1 spiked Untreated Pearson Correlation Log10(TPM+1)")
dev.off()

# Broth RIF
All_tpm_Log10_numeric_BrothRIF <- All_tpm_Log10_numeric %>% select(Ra_broth_RIF_1, Ra_broth_RIF_2, Ra_broth_RIF_3, Ra_broth_RIF_4)
pdf("Figures/Correlations_Scatterplots/rarity_Broth.RIF_v1.pdf", width = 10, height = 10)
corPlot(All_tpm_Log10_numeric_BrothRIF, method = "pearson",
        title = "Broth RIF Pearson Correlation Log10(TPM+1)")
dev.off()

# Broth Untreated
All_tpm_Log10_numeric_BrothUntreated <- All_tpm_Log10_numeric %>% select(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)
pdf("Figures/Correlations_Scatterplots/rarity_Broth.Untreated_v1.pdf", width = 10, height = 10)
corPlot(All_tpm_Log10_numeric_BrothRIF, method = "pearson",
        title = "Broth Untreated Pearson Correlation Log10(TPM+1)")
dev.off()


###################################################################
########################### MAKE AVERAGES #########################

# Add columns for averages
All_tpm_Log10_averages <- All_tpm_Log10 %>% 
  mutate(
    AVERAGE_Broth_Untreated = rowMeans(select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE),
    AVERAGE_Broth_RIF = rowMeans(select(., c(Ra_broth_RIF_1, Ra_broth_RIF_2, Ra_broth_RIF_3, Ra_broth_RIF_4)), na.rm = TRUE),
    AVERAGE_Spiked_Untreated = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    AVERAGE_Spiked_RIF = rowMeans(select(., c(THP1_spiked_Ra_RIF_1, THP1_spiked_Ra_RIF_2, THP1_spiked_Ra_RIF_3, THP1_spiked_Ra_RIF_4, THP1_spiked_Ra_RIF_5)))
  )

###################################################################
####################### AVERAGES GGCORRPLOT #######################

Sample1 <- "AVERAGE_Spiked_RIF" # THP1 spiked Captured
Sample2 <- "AVERAGE_Broth_RIF" # Broth Not Captured
ScatterCorr <- All_tpm_Log10_averages %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; All samples with RIF",
       x = paste0("Log10(TPM+1) Ra1e6 THP1 samples averaged"), y = paste0("Log10(TPM+1) Broth samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("RIF_THP1Spiked1e6_vs_Broth_Averages.pdf"),
       path = "Figures/Correlations_Scatterplots",
       width = 7, height = 5, units = "in")



Sample1 <- "AVERAGE_Broth_Untreated" # THP1 spiked Captured
Sample2 <- "AVERAGE_Broth_RIF" # Broth Not Captured
ScatterCorr <- All_tpm_Log10_averages %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; All samples Broth",
       x = paste0("Log10(TPM+1) Broth untreated samples averaged"), y = paste0("Log10(TPM+1) Broth RIF samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("Broth.Untreated_vs_Broth.RIF.pdf"),
       path = "Figures/Correlations_Scatterplots",
       width = 7, height = 5, units = "in")


Sample1 <- "AVERAGE_Spiked_Untreated" # THP1 spiked Captured
Sample2 <- "AVERAGE_Spiked_RIF" # Broth Not Captured
ScatterCorr <- All_tpm_Log10_averages %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; All samples THP1 spiked with 1e6 Ra",
       x = paste0("Log10(TPM+1) THP1 spiked untreated samples averaged"), y = paste0("Log10(TPM+1) THp1 spiked RIF samples averaged")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("Spiked.Untreated_vs_Spiked.RIF.pdf"),
       path = "Figures/Correlations_Scatterplots",
       width = 7, height = 5, units = "in")

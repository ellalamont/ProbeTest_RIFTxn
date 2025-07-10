# Comparing the log2fold change of RIF and untreated
# E. Lamont
# 7/10/25

source("Import_data.R")


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


################################################################
###### Spiked.UNTREATED.vs.RIF to Broth.UNTREATED.vs.RIF #######

# list_dfs_2$THP1_Untreat.ComparedTo.THP1_RIF
# list_dfs_2$Broth_Untreat.ComparedTo.Broth_RIF

df_spike <- list_dfs_2$THP1_Untreat.ComparedTo.THP1_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_THP1Spiked = LOG2FOLD,
         DE_THP1Spiked = DE) 
df_broth <- list_dfs_2$Broth_Untreat.ComparedTo.Broth_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_Broth = LOG2FOLD,
         DE_Broth = DE) 
df_DEG_Untreat.vs.RIF <- merge(df_spike, df_broth)


Sample1 <- "LOG2FOLD_THP1Spiked" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_Broth" # Broth Not Captured
my_plot <- df_DEG_Untreat.vs.RIF %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID), alpha = 0.7, size = 2, color = "black") +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples: ", Sample1, " vs ", Sample2),
       subtitle = "DEG between untreated and RIF",
       x = paste0("Log2Fold Ra1e6 THP1 samples"), y = paste0("Log2Fold Broth")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot


############################################################################
###### Broth.UNTREATED.vs.Spiked.RIF to Broth.UNTREATED.vs.Broth.RIF #######

# list_dfs_2$Broth_Untreat.ComparedTo.THP1_RIF
# list_dfs_2$Broth_Untreat.ComparedTo.Broth_RIF

df_spike <- list_dfs_2$Broth_Untreat.ComparedTo.THP1_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_THP1Spiked = LOG2FOLD,
         DE_THP1Spiked = DE) 
df_broth <- list_dfs_2$Broth_Untreat.ComparedTo.Broth_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_Broth = LOG2FOLD,
         DE_Broth = DE) 
df_combined <- merge(df_spike, df_broth)

df_combined <- df_combined %>%
  mutate(Combined_Sig = case_when(
    DE_THP1Spiked == "not significant" & DE_Broth == "not significant" ~ "both not significant",
    DE_THP1Spiked == "significant up" & DE_Broth == "significant down" ~ "significant in opposite directions",
    DE_THP1Spiked == "significant down" & DE_Broth == "significant up" ~ "significant in opposite directions",
    DE_THP1Spiked == "significant up" & DE_Broth == "significant up" ~ "both significant the same",
    DE_THP1Spiked == "significant down" & DE_Broth == "significant down" ~ "both significant the same",
    DE_THP1Spiked == "significant down" & DE_Broth == "not significant" ~ "only one significant",
    DE_THP1Spiked == "not significant" & DE_Broth == "significant down" ~ "only one significant",
    DE_THP1Spiked == "significant up" & DE_Broth == "not significant" ~ "only one significant",
    DE_THP1Spiked == "not significant" & DE_Broth == "significant up" ~ "only one significant"
  ))

Sample1 <- "LOG2FOLD_THP1Spiked" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_Broth" # Broth Not Captured
my_plot <- df_combined %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID, color = Combined_Sig), alpha = 0.7, size = 2) +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("black", "blue4", "#999990")) + 
  labs(title = "DEG both RIF compared to untreated broth",
       subtitle = "Pearson correlation",
       x = paste0("Log2Fold Broth.UNTREATED.vs.Spiked.RIF"), y = paste0("Log2Fold Broth.UNTREATED.vs.Broth.RIF")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot
ggsave(my_plot,
       file = paste0("RIFs_ComparedTo_Broth.Untreated.pdf"),
       path = "Figures/DEG_Comparisons",
       width = 7, height = 5, units = "in")


############################################################################
###### Spiked.UNTREATED.vs.Spiked.RIF to Spiked.UNTREATED.vs.Broth.RIF #####

# list_dfs_2$THP1_Untreat.ComparedTo.THP1_RIF
# list_dfs_2$THP1_Untreat.ComparedTo.Broth_RIF

df_spike <- list_dfs_2$THP1_Untreat.ComparedTo.THP1_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_THP1Spiked = LOG2FOLD,
         DE_THP1Spiked = DE) 
df_broth <- list_dfs_2$THP1_Untreat.ComparedTo.Broth_RIF %>% 
  select(GENE_ID, LOG2FOLD, DE) %>% 
  rename(LOG2FOLD_Broth = LOG2FOLD,
         DE_Broth = DE) 
df_combined <- merge(df_spike, df_broth)

df_combined <- df_combined %>%
  mutate(Combined_Sig = case_when(
    DE_THP1Spiked == "not significant" & DE_Broth == "not significant" ~ "both not significant",
    DE_THP1Spiked == "significant up" & DE_Broth == "significant down" ~ "significant in opposite directions",
    DE_THP1Spiked == "significant down" & DE_Broth == "significant up" ~ "significant in opposite directions",
    DE_THP1Spiked == "significant up" & DE_Broth == "significant up" ~ "both significant the same",
    DE_THP1Spiked == "significant down" & DE_Broth == "significant down" ~ "both significant the same",
    DE_THP1Spiked == "significant down" & DE_Broth == "not significant" ~ "only one significant",
    DE_THP1Spiked == "not significant" & DE_Broth == "significant down" ~ "only one significant",
    DE_THP1Spiked == "significant up" & DE_Broth == "not significant" ~ "only one significant",
    DE_THP1Spiked == "not significant" & DE_Broth == "significant up" ~ "only one significant"
  ))

Sample1 <- "LOG2FOLD_THP1Spiked" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_Broth" # Broth Not Captured
my_plot <- df_combined %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID, color = Combined_Sig), alpha = 0.7, size = 2) +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("black", "blue4", "#999990")) + 
  labs(title = "DEG both RIF compared to untreated THP1 spiked with 1e6 Ra",
       subtitle = "Pearson correlation",
       x = paste0("Log2Fold THP1Spiked.UNTREATED.vs.Spiked.RIF"), y = paste0("Log2Fold THP1Spiked.UNTREATED.vs.Broth.RIF")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot
ggsave(my_plot,
       file = paste0("RIFs_ComparedTo_Broth.Untreated.pdf"),
       path = "Figures/DEG_Comparisons",
       width = 7, height = 5, units = "in")



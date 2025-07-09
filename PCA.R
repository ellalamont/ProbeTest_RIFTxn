# PCA 
# E. Lamont
# 7/8/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

# source("Import_data.R") # to get All_tpm

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )



##########################################################
###################### ALL SAMPLES #######################

my_metadata <- All_pipeSummary
rownames(my_metadata) <- my_metadata[,1] # add the rownames
# my_metadata <- my_metadata[,-1] # Remove the old column of rownames
my_metadata <- my_metadata %>% 
  mutate(SampleType_Drug = paste0(Sample_Type, "_", Drug))

# Start with All_tpm

# Need to remove the gene column or it won't work
All_tpm_numeric <- All_tpm %>% select(-Gene)

# Transform the data
All_tpm_t <- as.data.frame(t(All_tpm_numeric))

# Remove columns that are all zero so the scale works for prcomp
All_tpm_t2 <- All_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(All_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 32.4% of variance
summary_PCA[2,1] # PC2 explains 6.2% of variance
summary_PCA[3,1] # PC3 explains 12.4% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID", )

fig_PC1vsPC2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = SampleType_Drug, shape = SampleType_Drug), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`Broth_Untreated` = "palegreen2", `Broth_RIF` = "green4", `THP1_Untreated`= "#FDBF6F", `THP1_RIF` = "#FF7F00")) +  
  scale_shape_manual(values=c(`Broth_Untreated` = 21, `Broth_RIF` = 22, `THP1_Untreated`= 21, `THP1_RIF` = 22)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA All samples",
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2
ggsave(fig_PC1vsPC2,
       file = "PCA_AllSamples.pdf",
       path = "Figures/PCA",
       width = 8, height = 5, units = "in")


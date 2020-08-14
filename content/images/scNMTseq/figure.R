
load('content/images/scNMTseq/liger-scnmtseq.RData')
ggsave(plot = liger_plots_promoter[[1]] + labs(col = 'Modality'), "figures/liger-promoter-modality.pdf")

# Create colour palettes for stages:
stages <- c("E4.5", "E5.5", "E6.5", "E7.5")
# col_palette <- dput(viridisLite::viridis(n = length(stages)))
col_palette <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
names(col_palette) <- stages

ggsave(plot = liger_plots_promoter[[2]] + labs(col = 'Stage') + scale_color_manual(values = col_palette), "figures/liger-promoter-stage.pdf")

ggsave(plot = liger_plots_genebody[[1]] + labs(col = 'Modality'), "figures/liger-genebody-modality.pdf")

ggsave(plot = liger_plots_genebody[[2]] + labs(col = 'Stage') + scale_color_manual(values = col_palette), "figures/liger-genebody-stage.pdf")
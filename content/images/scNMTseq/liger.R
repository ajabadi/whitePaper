library(liger)
library(MultiAssayExperiment)
library(BIRSBIO2020.scNMTseq.LIGER)

gastru.mae <- readRDS(url('https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=scnmtseq_gastrulation_mae_AllCells.rds'))
drop_lineages = c('Primitive_endoderm','Visceral_endoderm', 'ExE_ectoderm')
metadata = read.table("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=sample-metadata_matching-cells.txt",header=T,sep=",")
gastru.mae <- gastru.mae[,,grep("rna|met",names(assays(gastru.mae)))]
gastru.mae <- gastru.mae[, !(gastru.mae$lineage %in% drop_lineages),]
gastru.mae <- gastru.mae[,gastru.mae$pass_metQC==TRUE & gastru.mae$pass_rnaQC==TRUE]

proms <- experiments(gastru.mae)[["met_promoter"]]
genes <- experiments(gastru.mae)[["met_genebody"]]
rna <- counts(experiments(gastru.mae)[["rna"]])
proms[is.na(proms)] = 0
genes[is.na(genes)] = 0
all(colnames(proms)==colnames(rna))
colnames(proms)=paste0(colnames(proms),"_met")
colnames(genes)=paste0(colnames(genes),"_met")

rna_liger = createLiger(list(rna=rna,met=proms))
rna_liger = liger::normalize(rna_liger)
rna_liger = liger::selectGenes(rna_liger,var.thresh = 1,datasets.use=1)
rna_liger = liger::scaleNotCenter(rna_liger)

rna_liger@norm.data[["met"]] = rna_liger@raw.data[["met"]] = as.matrix(100-proms)
rna_liger@var.genes = rna_liger@var.genes[rna_liger@var.genes %in% rownames(rna_liger@raw.data[["met"]])]
rna_liger@raw.data[["met"]] = Matrix(rna_liger@raw.data[["met"]],sparse=T)
rna_liger@norm.data[["met"]] = Matrix(rna_liger@norm.data[["met"]],sparse=T)
rna_liger@scale.data[["met"]] = t(as.matrix(rna_liger@norm.data[["met"]][rna_liger@var.genes,]))

rna_liger = optimizeALS(rna_liger,k=20)
rna_liger = quantile_norm(rna_liger)
rna_liger = runTSNE(rna_liger)
metadata$id_met=paste0(metadata$id_met,"_met")
stage = c(as.character(metadata$stage),as.character(metadata$stage))
stage = as.factor(stage)
names(stage)=c(as.character(metadata$sample),as.character(metadata$id_met))
liger_plots_promoter <- plotByDatasetAndCluster(rna_liger,clusters=stage,pt.size=1,text.size = 0, return.plots = TRUE)

rna_liger = createLiger(list(rna=rna,met=genes))
rna_liger = liger::normalize(rna_liger)
rna_liger = liger::selectGenes(rna_liger,var.thresh = 1,datasets.use=1)
rna_liger = liger::scaleNotCenter(rna_liger)

rna_liger@norm.data[["met"]] = rna_liger@raw.data[["met"]] = genes
rna_liger@var.genes = rna_liger@var.genes[rna_liger@var.genes %in% rownames(rna_liger@raw.data[["met"]])]
rna_liger@raw.data[["met"]] = Matrix(rna_liger@raw.data[["met"]],sparse=T)
rna_liger@norm.data[["met"]] = Matrix(rna_liger@norm.data[["met"]],sparse=T)
rna_liger@scale.data[["met"]] = t(as.matrix(rna_liger@norm.data[["met"]][rna_liger@var.genes,]))

rna_liger = optimizeALS(rna_liger,k=20)
rna_liger = quantile_norm(rna_liger)
rna_liger = runTSNE(rna_liger)
liger_plots_genebody <- plotByDatasetAndCluster(rna_liger,clusters=stage,pt.size=1,text.size = 0, return.plots = TRUE)
rm(gastru.mae)
save.image(file = 'content/images/scNMTseq/liger-scnmtseq.RData')

# load('content/images/scNMTseq/liger-scnmtseq.RData')
# ggsave(plot = liger_plots_promoter[[1]] + labs(col = 'Modality'), "figures/liger-promoter-modality.pdf")
# 
# # Create colour palettes for stages:
# stages <- c("E4.5", "E5.5", "E6.5", "E7.5")
# # col_palette <- dput(viridisLite::viridis(n = length(stages)))
# col_palette <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
# names(col_palette) <- stages
# 
# ggsave(plot = liger_plots_promoter[[2]] + labs(col = 'Stage') + scale_color_manual(values = col_palette), "figures/liger-promoter-stage.pdf")
# 
# ggsave(plot = liger_plots_genebody[[1]] + labs(col = 'Modality'), "figures/liger-genebody-modality.pdf")
# 
# ggsave(plot = liger_plots_genebody[[2]] + labs(col = 'Stage') + scale_color_manual(values = col_palette), "figures/liger-genebody-stage.pdf")
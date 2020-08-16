library(BIRSBIO2020.scNMTseq.PLS)
library(MultiAssayExperiment)
library(scater)
library(scran)
library(mixOmics)
library(ggplot2)
library(magrittr)
library(reshape2)
library(uwot)

params <- list(
    on_my_mac = c(user = TRUE),
    save_output = FALSE,
    local_data = FALSE,
    mini_run = FALSE,
    matching_rna_for_umap = FALSE,
    drop_lineages = c("Primitive_endoderm", "Visceral_endoderm",
                      "ExE_ectoderm"),
    umap_params = c(
        run.seed = 42,
        n_neighbors = 15,
        n_components = 2,
        min_dist = 0.55
    )
)

cat('loading data from Cloudstor ...\n')
gastru.mae_path <- url('https://cloudstor.aarnet.edu.au/plus/s/jsW7nh4YwThw8Q5/download')

gastru.mae <- readRDS(gastru.mae_path)

## ------------- filter
## subset RNA expression and DNA methylation modalities
keep_assays <- grep("rna|met",names(assays(gastru.mae)))
gastru.mae <- gastru.mae[,,keep_assays]
## remove putative extraembryonic cells
cat(sprintf('dropping lineages %s, plus the unassigned lineages.\n', paste(params$drop_lineages, collapse = ', ')))
gastru.mae <- gastru.mae[,!(gastru.mae$lineage %in% params$drop_lineages),]
## keep cells which pass RNA QC
gastru.mae <- gastru.mae[,gastru.mae$pass_rnaQC==TRUE,]

## keep cells that also pass QC for DNA methylation
gastru.mae <- gastru.mae[,gastru.mae$pass_metQC==TRUE,]
## rna SCE for cells passing met and rna QC
rna.sce.matching <- gastru.mae@ExperimentList$rna
## ------------- RNA
# Create colour palettes for stages:
stages <- c("E4.5", "E5.5", "E6.5", "E7.5")
# col_palette <- dput(viridisLite::viridis(n = length(stages)))
col_palette <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
names(col_palette) <- stages

## helper function to create UMAP plots coloured by stage
plot_reducedDims <- function(sce, reducedDim = 'UMAP', 
                             dims = c(1,2), col_palette) {
    comps.prefix <- ifelse(reducedDim == 'UMAP', 'UMAP_', 'PC_')
    df <- data.frame(reducedDim(sce, reducedDim))[,dims]
    axes <- colnames(df) <- paste0(comps.prefix, dims)
    df$stage <- sce$stage
    ggplot(df, aes_string(axes[1], axes[2])) + geom_point(aes(col=stage)) +
        theme_classic()+ scale_color_manual(values = col_palette) + labs(col = 'Stage')
}

## ------------------------------ Matching Cells ------------------------------- ##
decomp <- modelGeneVar(rna.sce.matching)
## filter by mean expression and significance of biological variation signal
hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 0.01]
length(hvgs)
## ------------- UMAP
npcs <- 15
## PCA first: retrieve npcs PCs
rna.sce <- runPCA(rna.sce,  ncomponents = npcs, subset_row=hvgs, name = 'PCA')

## run UMAP
npcs <- 15
## PCA first: retrieve npcs PCs
rna.sce.matching <- runPCA(rna.sce.matching, ncomponents = npcs, subset_row=hvgs, name = 'PCA')

## UMAP parameters used:
cat(sprintf('Running UMAP with parameters %s\n', paste(names(params$umap_params), ':',params$umap_params, collapse = ', ')))

## run UMAP
set.seed(params$umap_params['run.seed'])
rna.sce.matching <- runUMAP(rna.sce.matching, dimred="PCA", 
                            ncomponents = params$umap_params['n_components'], 
                            n_neighbors = params$umap_params['n_neighbors'], 
                            min_dist = 0.35)
saveRDS(rna.sce, 'content/images/scNMTseq/rna.sce.matching.umap.rds')
## ------------------------------ All Cells ------------------------------- ##
decomp <- modelGeneVar(rna.sce)
## filter by mean expression and significance of biological variation signal
hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 0.01]
length(hvgs)
## ------------- UMAP
npcs <- 15
## PCA first: retrieve npcs PCs
rna.sce <- runPCA(rna.sce,  ncomponents = npcs, subset_row=hvgs, name = 'PCA')

## run UMAP
npcs <- 15
## PCA first: retrieve npcs PCs
rna.sce <- runPCA(rna.sce, ncomponents = npcs, subset_row=hvgs, name = 'PCA')

## UMAP parameters used:
cat(sprintf('Running UMAP with parameters %s\n', paste(names(params$umap_params), ':',params$umap_params, collapse = ', ')))

## run UMAP
set.seed(params$umap_params['run.seed'])
rna.sce <- runUMAP(rna.sce, dimred="PCA", 
                            ncomponents = params$umap_params['n_components'], 
                            n_neighbors = params$umap_params['n_neighbors'], 
                            min_dist = 0.35)

saveRDS(rna.sce, 'content/images/scNMTseq/rna.sce.umap.rds')

'/Users/alabadi/Projects/analysis/R/_fork/multiOmics/BIRSBIO2020/whitePaper/content/images/scNMTseq/rna.sce.umap.rds'

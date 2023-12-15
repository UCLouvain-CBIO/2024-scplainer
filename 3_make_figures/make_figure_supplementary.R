#### SUPPLEMENTARY FIGURES


## This script generates the supplementary figures about miscellaneous
## topics to support the claims in the body of the work.

####---- Loading libraries and data ----####

## libraries
library("scp")
library("ggplot2")
library("patchwork")
library("limma")
library("MsCoreUtils")
library("scater")
library("dplyr")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"

## define colors
populationColors <- c(
    Monocyte = "coral", 
    Melanoma = "skyblue3",
    "Melanoma_main" = "skyblue3",
    "Undefined Melanoma" = "gray",
    "Resistant Melanoma" = "skyblue",
    "Melanoma_sub" = "skyblue",
    PDAC = "khaki"
)

####---- Workflow: order matters ----####

## Batch correction and imputation is not the same as imputation and 
## batch correction
sce <- readRDS("scripts/data/leduc2022_pSCoPE_modelled.rds")
# sce <- filterNA(sce, pNA = 0.9)
Y <- scpModelInput(sce)
annot <- colData(sce)
## Batch correct then impute
Ybi <- removeBatchEffect(Y, batch = annot$Set,
                         design = model.matrix(~ SampleType, annot))
Ybi <- impute_knn(Ybi, k = 10, rowmax = 1, colmax = 1)
YbiPca <- svd(Ybi, nv = 2, nu = 2)
## Impute then batch correct
Yib <- impute_knn(Y, k = 10, rowmax = 1, colmax = 1)
Yib <- removeBatchEffect(Yib, batch = annot$Set,
                         design = model.matrix(~ SampleType, annot))
YibPca <- svd(Yib, nv = 2, nu = 2)
## Plot comparison
m <- as.vector(Yib - Ybi)
a <- as.vector(Yib + Ybi) / 2
set.seed(1234)
sub <- sample(1:length(m), 1E5) ## avoid too many points
panel1 <- ggplot(data.frame(M = m[sub], A = a[sub])) +
    aes(x = M,
        y = A,) +
    labs(x = "Method1 - Method2",
         y = "(Method1 + Method2) / 2") +
    geom_point(alpha = 0.1, size = 0.5) +
    theme_minimal()
pcLabs <- scp:::.createAxisLabelsForPCA(YbiPca$d / sum(YbiPca$d), comp = 1:2)
panel2 <- ggplot(data.frame(PC = YbiPca$v, colData(sce))) +
    aes(x = PC.1,
        y = PC.2,
        colour = SampleType) +
    labs(x = pcLabs$x,
         y = pcLabs$y,
         title = "Method1\nBatch correct and impute") +
    geom_point(alpha = 0.3) +
    theme_minimal()
pcLabs <- scp:::.createAxisLabelsForPCA(YibPca$d / sum(YibPca$d), comp = 1:2)
panel3 <- ggplot(data.frame(PC = YibPca$v, colData(sce))) +
    aes(x = PC.1,
        y = PC.2,
        colour = SampleType) +
    labs(x = pcLabs$x,
         y = pcLabs$y,
         title = "Method2\nImpute and batch correct") +
    geom_point(alpha = 0.3) +
    theme_minimal()
(fig <- panel2 +
        panel3 +
        panel1 +
        plot_layout(design = "12
                              33",
                    guides = "collect") +
        plot_annotation(tag_levels = "a") &
        scale_colour_manual(values = populationColors))
ggsave("scripts/figs/supp_workflow_order_matters.pdf", fig, 
       width = 9, height = 7)

####---- Workflow: steps inside or outside modelling ----####

modelData <- function(object, normType, batchCorrectionType, subset) {
    object <- object[subset, ]
    fString <- "~ 1 + SampleType"
    if (normType == "before") {
        assay(object) <- sweep(assay(object), 2, object$MedianIntensity, "-")
    } else if (normType == "model") {
        fString <- paste(fString, "+ MedianIntensity")
    } else {
        stopifnot(normType == "none")
    }
    if (batchCorrectionType == "before") {
        object <- batchCorrect(object)
    } else if (batchCorrectionType == "model") {
        fString <- paste(fString, "+ Channel + Set")
    } else {
        stopifnot(batchCorrectionType == "none")
    }
    scpModelWorkflow(object, as.formula(fString))
}

batchCorrect <- function(object) {
    assay(object) <- removeBatchEffect(
        assay(object), batch = object$Set, batch2 = object$Channel, 
        design = model.matrix(~ SampleType, colData(object))
    )
    object
}

strategies <- do.call(rbind, list(
    c(1, "none", "none"),
    c(2, "before", "none"),
    c(3, "model", "none"),
    c(4, "none", "before"),
    c(5, "none", "model"),
    c(6, "before", "before"),
    c(7, "before", "model"),
    c(8, "model", "before"),
    c(9, "model", "model")
))
colnames(strategies) <- c("strategy", "normalisation", "batch_correction")
strategies <- data.frame(strategies)

sce <- readRDS("scripts/data/leduc2022_pSCoPE_modelled.rds")
source("scripts/3_make_figures/utils.R")

results <- list()
system.time({
    for (i in 1:nrow(strategies)) {
        normType <- strategies$normalisation[[i]]
        batchCorrectionType <- strategies$batch_correction[[i]]
        print(strategies[i, ])
        bc <- modelData(
            sce, normType, batchCorrectionType,
            subset = ascatmp:::scpModelFeatureNames(sce)
        )
        results[[i]] <- benchmarkData(bc, bioVar = "SampleType", techVar = "Set")
        results[[i]]$strategy <- strategies[i, ]
    }
})
## Takes a bit less than 1h
# save(results, file = "package/data/results.rda")
results <- lapply(results, function(r) {
    r$clustering <- clusterData(r$dimred, sce$SampleType)
    set.seed(1234)
    r$tsne <- calculateTSNE(t(r$dimred))
    r
})
## Plot metrics
(figMetrics <- lapply(results, function(r) cbind(r$metrics, r$strategy)) |> 
        do.call(what = rbind) |> 
        mutate(label =  paste("Strategy", strategy), 
               strategy = paste(
                   "Strategy", strategy,
                   "\nNorm:", normalisation, 
                   "\nBatch correction: ", batch_correction, "\n"
               )
        ) |> 
        ggplot() +
        aes(x = biological,
            y = 1 - technical,
            colour = strategy) +
        geom_point() +
        geom_text_repel(aes(label = label), size = 3) +
        facet_wrap(~ metric, scales = "free", ncol = 1) +
        theme_minimal())

## Plot PCA
(fig <- wrap_plots(lapply(results, function(r) {
    strategy <- paste(
        "Norm:", r$strategy$normalisation, 
        "\nBatch correction: ", r$strategy$batch_correction
    )
    ggplot(data.frame(r$dimred, cluster = r$clustering)) +
        aes(x = PC1,
            y = PC2,
            colour = sce$SampleType) +
        geom_point(size = 1, alpha = 0.5) +
        ggtitle(strategy) +
        theme_minimal()
}), guides = "collect"))

## Plot t-SNE
(figTSNE <- wrap_plots(lapply(results, function(r) {
    strategy <- paste("Strategy", r$strategy)
    # strategy <- paste(
    #     "Norm:", r$strategy$normalisation, 
    #     "\nBatch correction: ", r$strategy$batch_correction
    # )
    ggplot(data.frame(r$tsne, cluster = r$clustering)) +
        aes(x = TSNE1,
            y = TSNE2,
            colour = sce$SampleType) +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = populationColors, name = "Cell Type") +
        ggtitle(strategy) +
        theme_minimal()
})))

(fig <- figMetrics +
        figTSNE +
        plot_annotation(tag_levels = "a") +
        plot_layout(guides = "collect", 
                    design = "1222"))
ggsave("scripts/figs/supp_workflow_steps.pdf", fig, 
       width = 12, height = 8)

####---- Workflow: n/p filter ----####

sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))
pNA <- nNA(sce)$nNArows$pNA
npRatio <- scpModelFilterNPRatio(sce)
(fig <- ggplot(data.frame(pNA, npRatio)) +
        aes(x = pNA,
            y = npRatio) +
        geom_point(size = 0.5, alpha = 0.2) +
        geom_vline(xintercept = 0.95, colour = "firebrick", linewidth = 1.5) +
        geom_hline(yintercept = 1, colour = "orange2", linewidth = 1.5) +
        theme_minimal())
ggsave("scripts/figs/supp_workflow_filter.pdf", fig, 
       width = 4, height = 4)

####---- Variance analysis ----####

## VIM is not influenced by cell size
sce <- readRDS("scripts/data/leduc2022_pSCoPE_modelled.rds")
vaRes <- scpVarianceAnalysis(sce) |> 
    scpAnnotateResults(
        rowData(sce)[, c("Sequence", "gene")], 
        by = "feature", by2 = "Sequence"
    )
## Best peptide
target <- "CORO1A"
i <- data.frame(vaRes$SampleType) |> 
    filter(!is.na(gene) & gene == target) |> 
    slice_max(SS) |> 
    pull(feature)
## Create annotation with cell sizes
cellenoneData <- read.csv("~/PhD/.localData/SCP/leduc2022/pSCoPE/cellenoneData.csv", row.names = 1)
cellenoneData$Channel <- as.numeric(sub("^RI(\\d*)$", "\\1", cellenoneData$Label))
tmts <- c("126", paste0(rep(127:134, each = 2), rep(c("N", "C"), 8)), "135N")
cellenoneData$Channel <- paste0("TMT", tmts[cellenoneData$Channel])
annot <- left_join(data.frame(colData(sce)),
                   cellenoneData[, c("Set", "Channel", "Diameter", "Elongation")],
                   by = c("Set", "Channel"))
(fig <- ggplot(data.frame(logIntensity = assay(scpModelEffects(sce)[["MedianIntensity"]])[i, ],
                          annot)) +
        ggtitle("Normalization effect") +
        ggplot(data.frame(logIntensity = assay(scpKeepEffect(sce, c("SampleType")))[i, ],
                          annot)) +
        ggtitle("Cell type effect") +
        plot_layout(guides = "collect") +
        plot_annotation(title = paste(
            "Effect of cell diameter on", i, "(", target, ") intensity"
        ), tag_levels = "a") &
        aes(x = Diameter,
            y = logIntensity,
            colour = SampleType) &
        geom_point() &
        scale_colour_manual(values = populationColors) &
        labs(x = "Cell diameter (µm)",
             y = "log2(Intensity)",
             colour = "Cell type") &
        theme_minimal())
ggsave(
    "scripts/figs/supp_variance_diameter.pdf", fig,
    width = 8, height = 4
)

####---- Differential analysis ----####

## none

####---- Component analysis ----####

sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))
## Adapt annotations
sce$cell <- colnames(sce)
sce$Population <- ifelse(
    !is.na(sce$MelanomaSubCluster) & sce$MelanomaSubCluster == "B",
    "Resistant Melanoma", sce$SampleType
)
sce$Population <- ifelse(
    is.na(sce$MelanomaSubCluster) & sce$Population == "Melanoma",
    "Undefined Melanoma", sce$Population
)
caRes <- scpComponentAnalysis(
    sce,
    ncomp = 20,
    maxiter = 50,
    method = "APCA",
    effect = "SampleType",
    residuals = FALSE
)

## Variance explained by t-SNE
panels <- list()
for (what in c("unmodelled", "APCA_SampleType")) {
    if (what == "unmodelled") {
        pvar <- metadata(caRes$bySample$unmodelled)$proportionVariance
        title <- "Raw (unmodelled) data"
    } else {
        pvar <- metadata(caRes$bySample$APCA_SampleType)$proportionVariance
        title <- "Batch corrected data"
    }
    df <- data.frame(pvar = pvar, comp = seq_along(pvar)
    )
    pl <- ggplot(df) +
        aes(y = pvar * 100, x = comp) +
        geom_bar(stat = "identity") +
        labs(
            x = "Component",
            y = "Percentage variance explained",
            title = paste0(
                "Variance explained by top 20 PCs\n",
                "for the ", title),
            subtitle = paste0(
                round(sum(df$pvar) * 100), 
                "% cumulated variance explained"
            )) +
        theme_minimal()
    panels <- c(panels, list(pl))
}
fig <- wrap_plots(panels) +
    plot_annotation(tag_levels = "a")
ggsave("scripts/figs/supp_component_explained_var.pdf", fig, height = 4, width = 7)

## t-SNE for different variables
pcaUnmodelled <- caRes$bySample$unmodelled
pcaUnmodelled <- pcaUnmodelled[, grepl("^PC", colnames(pcaUnmodelled))]
set.seed(1234)
tsneUnmodelled <-calculateTSNE(t(as.matrix(pcaUnmodelled)))
apcaSampleType <- caRes$bySample$APCA_SampleType
apcaSampleType <- apcaSampleType[, grepl("^PC", colnames(apcaSampleType))]
set.seed(1234)
tsneSampleType <-calculateTSNE(t(as.matrix(apcaSampleType)))
supp <- list()
for (what in c("unmodelled", "APCA_SampleType")) {
    if (what == "unmodelled") {
        x <- tsneUnmodelled
        title <- "Raw (unmodelled) data"
    } else {
        x <- tsneSampleType
        title <- "Batch corrected data"
    }
    for (effect in scpModelEffectNames(sce)) {
        pl <- data.frame(x, colData(sce)) |> 
            ggplot() +
            aes(x = TSNE1,
                y = TSNE2,
                colour = .data[[effect]]) +
            geom_point(alpha = 0.6) +
            labs(title = title, 
                 subtitle = paste0("Coloured by ", effect)) +
            theme_minimal() +
            theme(legend.position = "none")
        supp <- c(supp, list(pl))
    }
}
(fig <- wrap_plots(supp) +
        plot_layout(nrow = 2) +
        plot_annotation(tag_levels = "a"))
ggsave("scripts/figs/supp_component_batch_qc.pdf", fig, height = 8, width = 12)

## Cell cluster analysis on raw data
clusters <- list()
set.seed(11) ## set seed to match clustering of the main component fig
clusters$unmodelled <- as.factor(kmeans(pcaUnmodelled, 3)$cluster)
set.seed(11) ## set seed to match clustering of the main component fig
clusters$batchCorrected <- as.factor(kmeans(apcaSampleType, 3)$cluster)
tsnes <- list(
    unmodelled = tsneUnmodelled, 
    batchCorrected = tsneSampleType
)
panels <- list()
for (cl in names(clusters)) {
    for (dat in names(tsnes)) {
        df <- data.frame(
            tsnes[[dat]], colData(sce), ClusterID = clusters[[cl]]
        )
        pl <- ggplot(df) +
            aes(x = TSNE1,
                y = TSNE2,
                colour = ClusterID,
                shape = lcbatch) +
            geom_point() +
            labs(title = paste("t-SNE on", dat, "data"),
                 subtitle = paste("Clustering on", cl, "data")) +
            theme_minimal()
        panels <- c(panels, list(pl))
    }
}
(fig <- wrap_plots(panels) +
        plot_layout(guide = "collect"))
ggsave("scripts/figs/supp_component_clustering.pdf", fig, height = 7, width = 7)

## Peptide cluster analysis for SampleType

if (FALSE) { ## in case I would need this
    ## Cluster analysis
    caResProts <- lapply(caResProts, function(x) {
        pci <- grep("^PC\\d", colnames(x))
        x$cluster <- as.factor(kmeans(x[, pci], 150)$cluster)
        x
    })
    ## GO enrichment analysis on sample type
    compBiol <- caResProts$APCA_SampleType
    enrichResBiol <- lapply(split(compBiol, compBiol$cluster), function(x) {
        res <- enrichGO(
            x$gene, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
            pAdjustMethod = "BH", qvalueCutoff = 0.05, 
            universe = compBiol$gene, minGSSize = 10, maxGSSize = 50
        )
        pathway <- head(res[, "Description"], 1)
        if (!length(pathway)) pathway <- NA
        data.frame(
            pathway = pathway, cluster = unique(x$cluster)
        )
    })
    enrichResBiol <- do.call(rbind, enrichResBiol)
    ## Plot
    caResProtsSampleType <- scpAnnotateResults(
        caResProts["APCA_SampleType"], enrichResBiol, by = "cluster"
    )
    (panel3 <- scpComponentPlot(
        caResProtsSampleType, pointAes = aes(colour = pathway)
    )[[1]] +
            geom_label_repel(data = data.frame(caResProtsSampleType[[1]]) |> 
                                 group_by(pathway) |> 
                                 summarise(PC1 = median(PC1),
                                           PC2 = median(PC2)),
                             aes(x = PC1, 
                                 y = PC2,
                                 colour = pathway,
                                 label = pathway),
                             size = 3) +
            theme(legend.position = "none"))
    
    ## Peptide cluster analysis for residuals
    
    ## GO enrichment analysis on residuals
    compResid <- caResProts$residuals
    enrichResResid <- lapply(split(compResid, compResid$cluster), function(x) {
        res <- enrichGO(
            x$gene, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
            pAdjustMethod = "BH", qvalueCutoff = 0.05, 
            universe = compResid$gene, minGSSize = 10, maxGSSize = 50
        )
        pathway <- head(res[, "Description"], 1)
        if (!length(pathway)) pathway <- NA
        data.frame(
            pathway = pathway, cluster = unique(x$cluster)
        )
    })
    enrichResResid <- do.call(rbind, enrichResResid)
    ## Plot
    caResProtsResid <- scpAnnotateResults(
        caResProts["residuals"], enrichRes, by = "cluster"
    )
    (panel4 <- scpComponentPlot(
        caResProtsResid["residuals"], pointAes = aes(colour = pathway)
    )[[1]] +
            geom_label_repel(data = data.frame(caResProtsResid[[1]]) |> 
                                 group_by(pathway) |> 
                                 summarise(PC1 = median(PC1),
                                           PC2 = median(PC2)),
                             aes(x = PC1, 
                                 y = PC2,
                                 colour = pathway,
                                 label = pathway),
                             size = 3) +
            theme(legend.position = "none"))
    
    ## Emert's markers
    
    ## The marker table was retrieved from the slavov lab in the repo:
    ## https://drive.google.com/drive/folders/116XGOyPpn6g0qWz8i32ZU215kjSZm-KJ
    markers <- read.csv("~/PhD/.localData/SCP/leduc2022/pSCoPE/misc/Melanoma_diff_exp_raj.csv")
    markers <- markers[!is.na(markers$padj) & markers$padj < 0.1, ]
    up <- markers$gene_id[markers$log2FoldChange > 0]
    down <- markers$gene_id[markers$log2FoldChange < 0]
    proteinIds <- unique(rowData(sce)$Leading.razor.protein.id)
    proteinConversionDf <- transcripts(
        EnsDb.Hsapiens.v86, 
        columns = "gene_id",
        return.type = "data.frame",
        filter = UniprotFilter(proteinIds)
    )
    annot <- proteinConversionDf[, -2]
    annot$ResistanceMarker <- annot$gene_id %in% markers$gene_id
    annot$MarkerDirection <- "none"
    annot$MarkerDirection[annot$gene_id %in% up] <- "up"
    annot$MarkerDirection[annot$gene_id %in% down] <- "down"
    caResPeps <- scpAnnotateResults(
        caResPeps, annot, by = "Leading.razor.protein.id", by2 = "uniprot_id"
    )
    plEmert <- scpComponentPlot(
        .scaleComponentsToUnity(caResPeps, paste0("PC", 1:2)),
        pointAes = aes(alpha = ResistanceMarker,
                       colour = MarkerDirection))
    (panel3 <- plEmert$APCA_SampleType +
            scale_colour_manual(values = c(down = "dodgerblue", 
                                           up = "coral", 
                                           none = "gray30")) +
            scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.05),
                               na.value = 0))
}

####---- Integration analysis ----####

## Compare with previous analysis

leduc <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))
integration <- readRDS(paste0(dataDir, "integration_modelled.rds"))
rowData(integration)$Peptide <- rownames(integration)
dfLeduc <- scpDifferentialAnalysis(
    leduc, contrasts = list(c("SampleType", "Melanoma", "Monocyte"))
) |>
    scpAnnotateResults(
        rowData(leduc)[, c("Sequence", "Leading.razor.protein.id", "gene")],
        by = "feature", by2 = "Sequence"
    ) |>
    getListElement(1) |> 
    data.frame() |> 
    dplyr::select(padj, Estimate, feature) |> 
    dplyr::rename(padj_leduc = padj, logfc_leduc = Estimate)
dfIntegration <- scpDifferentialAnalysis(
    integration, contrasts = list(c("Celltype", "Melanoma", "Monocyte"))
) |>
    scpAnnotateResults(
        rowData(integration)[, c("Peptide", "ProteinId", "Gene")],
        by = "feature", by2 = "Peptide"
    ) |>
    getListElement(1) |> 
    data.frame() |> 
    dplyr::select(padj, Estimate, feature) |> 
    dplyr::filter(!is.na(Estimate)) |> 
    dplyr::rename(padj_integr = padj, logfc_integr = Estimate)
dfAll <- inner_join(dfIntegration, dfLeduc, by = "feature")
nrow(dfAll)
table(sign(dfAll$logfc_integr), sign(dfAll$logfc_leduc))
#     -1   1
# -1 874 542
# 1  342 803
table(sign(dfAll$logfc_integr), sign(dfAll$logfc_leduc)) / nrow(dfAll)
#           -1         1
# -1 0.3412729 0.2116361
# 1  0.1335416 0.3135494
table(dfAll$padj_integr < 0.05, dfAll$padj_leduc < 0.05)
#       FALSE TRUE
# FALSE   703  396
# TRUE    519  531
table(dfAll$padj_integr < 0.05, dfAll$padj_leduc < 0.05) / nrow(dfAll)
#           FALSE      TRUE
# FALSE 0.3121670 0.1758437
# TRUE  0.2304618 0.2357904

## The two analyses provide contradictory results for some peptides
i <- "SFVLNLGK" 
prot <- rowData(integration)[i, "Gene"]
(fig <- ggplot(data.frame(
    logIntensity = assay(integration)[i, ], 
    colData(integration)
)) +
        aes(x = dataset,
            y = logIntensity,
            fill = Celltype) +
        labs(y = "Unmodelled log2(Intensity)", 
             x = "Data set") +
        ggtitle("Integration data set\nBefore modelling") + 
        ggplot(data.frame(
            logIntensity = assay(leduc)[i, ], 
            colData(leduc)
        )) +
        aes(x = lcbatch,
            y = logIntensity,
            fill = SampleType) +
        labs(y = "Unmodelled log2(Intensity)", 
             x = "Chromatographic batch") +
        ggtitle("leduc2022 data set\nBefore modelling") + 
        ggplot(data.frame(
            logIntensity = assay(scpKeepEffect(integration, "Celltype"))[i, ], 
            colData(integration)
        )) +
        aes(x = dataset,
            y = logIntensity,
            fill = Celltype) +
        labs(y = "Batch-corrected log2(Intensity)", 
             x = "Data set") +
        ggtitle("Integration data set\nAfter modelling") + 
        ggplot(data.frame(
            logIntensity = assay(scpKeepEffect(leduc, "SampleType"))[i, ], 
            colData(leduc)
        )) +
        aes(x = lcbatch,
            y = logIntensity,
            fill = SampleType) +
        labs(y = "Batch-corrected log2(Intensity)", 
             x = "Chromatographic batch") +
        ggtitle("leduc2022 data set\nAfter modelling") + 
        plot_layout(nrow = 2) +
        plot_annotation(paste(i, "(", prot, ")"), tag_levels = "a") &
        geom_boxplot() &
        scale_fill_manual(values = populationColors) &
        theme_minimal() &
        theme(axis.text.x = element_text(angle = 30, hjust = 1)))
ggsave("scripts/figs/supp_integration.pdf", fig, height = 9, width = 9)

####---- Re-explore leduc2022 with clusters ----####

sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_re_modelled.rds"))

## Variance analysis
vaRes <- scpVarianceAnalysis(sce, "model_with_cluster")
vaRes <- scpAnnotateResults(
    vaRes, rowData(sce), by = "feature", by2 = "Sequence"
)
(panel1 <- scpVariancePlot(vaRes))
## Component analysis
caRes <- scpComponentAnalysis(
    sce, name = "model_with_cluster",
    ncomp = 2,
    maxiter = 50,
    method = "APCA",
    effect = "Cluster",
    residuals = FALSE,
    unmodelled = FALSE
)
caResCells <- caRes$bySample
sce$cell <- colnames(sce)
caResCells <- scpAnnotateResults(caResCells, colData(sce), by = "cell")
(panel2 <- 
        scpComponentPlot(
            caResCells, pointParams = list(aes(colour = Cluster))
        )[[1]] +
        scale_colour_manual(values = populationColors, 
                            labels = c("Melanoma main", "Melanoma subpopulation", "Monocyte"))
)
## Differential analysis on clusters
daResCluster <- scpDifferentialAnalysis(
    sce, name = "model_with_cluster",
    contrasts = list(
        c("Cluster", "Melanoma_main", "Monocyte"),
        c("Cluster", "Melanoma_main", "Melanoma_sub")
    )
)
daResCluster <- scpAnnotateResults(
    daResCluster, rowData(sce), by = "feature", by2 = "Sequence"
)
## PSEA on differences between main melanoma and sub melanoma
library("msigdbr")
library("fgsea")
library("BiocParallel")
goBiologicalProcess <- msigdbr(species = "Homo sapiens", category = "C5", 
                               subcategory = "BP")
goBiologicalProcess <- split(goBiologicalProcess$gene_symbol, 
                             goBiologicalProcess$gs_name)
daResProt <- scpDifferentialAggregate(daResCluster, "gene")
rank <- daResProt$Cluster_Melanoma_main_vs_Melanoma_sub$Estimate
names(rank) <-  daResProt$Cluster_Melanoma_main_vs_Melanoma_sub$gene
rank <- sort(rank, decreasing = FALSE)
set.seed(1234)
pseaRes <- out <- fgsea(goBiologicalProcess, 
                         rank,
                         minSize = 20,
                         maxSize = 200,
                         nPermSimple = 2000,
                         BPPARAM = MulticoreParam(workers = 4))
pseaRes <- pseaRes[pseaRes$padj < 0.05, ]
(pseaRes <- pseaRes[order(pseaRes$padj), ])
## Save PSEA results
out$leadingEdge <- sapply(out$leadingEdge, paste, collapse = ", ")
write.csv(
    x = out, file = "data/PSEA_results.csv",
    row.names = FALSE
)
## Create annotation table with 2 pathways
genes <- rowData(sce)$gene
p <- "GOBP_PYRUVATE_METABOLIC_PROCESS"
i <- which(pseaRes$pathway == p)
pathway <- ifelse(genes %in% pseaRes$leadingEdge[[i]], p, "other")
p <- "GOBP_OXIDATIVE_PHOSPHORYLATION"
i <- which(pseaRes$pathway == p)
pathway <- ifelse(genes %in% pseaRes$leadingEdge[[i]], p, pathway)
pathway <- gsub("GOBP_", "", pathway)
pathway <- gsub("_", " ", pathway)
pathway <- tolower(pathway)
daResCluster <- scpAnnotateResults(
    daResCluster, data.frame(gene = genes, pathway = pathway), by = "gene"
)
## Plot volcano with PSEA
volcanos <- scpVolcanoPlot(
    daResCluster, textBy = "gene", 
    labelParams = list(max.overlaps = 30),
    pointParams = list(
        aes(colour = pathway, 
            alpha = pathway == "other")
    )
)
volcanos <- lapply(names(volcanos), function(i) {
    title <- if (grepl("Mono", i)) {
        "Monocyte <---> Melanoma cell"
    } else {
        "Melanoma\nSubpopulation <---> Main"
    }
    volcanos[[i]] +
        ggtitle(title)
})
(panel3 <- volcanos[[2]] +
        scale_alpha_manual(values = c(1, 0.2)) +
        scale_color_manual(values = c("grey50", "green4", "orange2")) +
        guides(alpha = "none") +
        theme_minimal())

## Create figure
(fig <- panel1 +
        wrap_elements(full = panel2) +
        panel3 +
        plot_layout(
            design = "1222
                      3333") +
        plot_annotation(tag_levels = "a"))
ggsave("figs/supp_remodel.pdf", fig, height = 6, width = 7)

## Plot normalised enrichment scores
(fig <- pseaRes |> 
        filter(padj < 0.05) |> 
        mutate(
            pathway = tolower(pathway),
            pathway = sub("gobp_", "", pathway),
            pathway = gsub("_", " ", pathway),
            pathway = factor(pathway, levels = pathway[order(NES)]),
            type = ifelse(NES > 0, "Melanoma main", "Melanoma subpopulation")
        ) |> 
        ggplot() +
        aes(y = pathway, 
            x = NES, 
            fill = type) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
            values = c(
                "Melanoma main" =  "skyblue3", "Melanoma subpopulation" = "skyblue"
            )
        ) +
        labs(x = "Normalised enrichment score",
             y = "",
             fill = "") +
        theme_minimal() +
        theme(legend.position = "none"))
ggsave("figs/supp_NES.pdf", fig, height = 4, width = 6)

## Compare the differential analysis between modelling sample type vs
## cluster
## Differential analysis on sampletype
daResSt <- scpDifferentialAnalysis(
    sce, name = "model",
    contrasts = list(c("SampleType", "Melanoma", "Monocyte"))
)
## Compare sample type vs cluster analysis
common <- intersect(daResSt[[1]]$feature, daResCluster[[1]]$feature)
data.frame(
    SampleTypeAnalysis = daResSt[[1]][match(common, daResSt[[1]]$feature), "Estimate"],
    ClusterAnalysis = daResCluster[[1]][match(common, daResCluster[[1]]$feature), "Estimate"]
) |> 
    ggplot() +
    aes(x = SampleTypeAnalysis,
        y = ClusterAnalysis) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0)


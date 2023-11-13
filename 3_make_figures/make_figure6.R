#### INTEGRATION STUDY

## This script generates the figure where we demonstrate the results of
## dataset integration

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("ggplot2")
library("patchwork")
library("RColorBrewer")
library("scater")
library("dplyr")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"
sce <- readRDS(paste0(dataDir, "integration_modelled.rds"))

cellTypeCol <- c(Monocyte = "coral", Melanoma = "skyblue3", PDAC = "khaki")

####---- Panel: Variance analysis ----####
vaRes <- scpVarianceAnalysis(sce)
cols <- brewer.pal(n = 8, name = "Set2")
cols <- c("white", cols[c(4, 6, 2, 5)])
names(cols) <- names(vaRes)
(panel1 <- scpVariancePlot(vaRes) +
        scale_fill_manual(
            values = cols,
            labels = c(
                "Residuals", "Normalization Factor", "Batch", 
                "Label", "Biology"
            )
        )
)

####---- Panel: Volcano plot ----####

## Note that VIM is no longer present because not selected for 
## modelling
daRes <- scpDifferentialAnalysis(
    sce, contrasts = list(
        c("Celltype", "Melanoma", "Monocyte"),
        c("Celltype", "Melanoma", "PDAC"),
        c("Celltype", "Monocyte", "PDAC")
    )
)
rowData(sce)$Peptide <- rownames(sce)
daRes <- scpAnnotateResults(
    daRes, rowData(sce)[, c("Peptide", "ProteinId", "Gene")],
    by = "feature", by2 = "Peptide"
)
vplots <- scpVolcanoPlot(
    daRes, textBy = "Gene", top = 200,
    labelParams = list(
        max.overlaps = 20,
        size = 2
    )
)
(panel2 <- vplots$Celltype_Melanoma_vs_Monocyte +
        labs(x = "log (Fold change)", 
             title = "Monocyte <---> Melanoma cell") +
        geom_point(data = data.frame(daRes[[1]]) |> 
                       dplyr::filter(!is.na(Gene) & Gene == "VIM"),
                   colour = "red3") +
        geom_point(data = data.frame(daRes[[1]]) |> 
                       dplyr::filter(!is.na(Gene) & Gene == "CTTN"),
                   colour = "orange2") +
        geom_point(data = data.frame(daRes[[1]]) |> 
                       dplyr::filter(!is.na(Gene) & Gene == "ARHGDIB"),
                   colour = "dodgerblue") +
        geom_point(data = data.frame(daRes[[1]]) |> 
                       dplyr::filter(!is.na(Gene) & Gene == "LCP1"),
                   colour = "purple2") +
        theme_minimal())

####---- Panel: Component analysis ----####

caRes <- scpComponentAnalysis(
    sce, method = "APCA", effects = "Celltype",
    ncomp = 20, maxiter = 50, residuals = FALSE, 
)

## Data before correction
pcaUnmodelled <- caRes$bySample$unmodelled
pcaUnmodelled <- pcaUnmodelled[, grepl("^PC", colnames(pcaUnmodelled))]
set.seed(1111)
tsneUnmodelled <- calculateTSNE(t(as.matrix(pcaUnmodelled)))
(panel3 <- ggplot(data.frame(tsneUnmodelled, colData(sce))) +
        aes(x = TSNE1,
            y = TSNE2,
            fill = Celltype,
            shape = dataset) +
        geom_point(alpha = 0.6, size = 2) +
        ggtitle("Raw (unmodelled)") +
        scale_shape_manual(values = 21:24) +
        scale_fill_manual(values = cellTypeCol) + 
        guides(fill = guide_legend(override.aes=list(shape=21))) +
        theme_minimal())

## Data after batch correction
pcaBc <- caRes$bySample$APCA_Celltype
pcaBc <- pcaBc[, grepl("^PC", colnames(pcaBc))]
set.seed(1111)
tsneBc <- calculateTSNE(t(as.matrix(pcaBc)))
(panel4 <- ggplot(data.frame(tsneBc, colData(sce))) +
        aes(x = TSNE1, 
            y = TSNE2, 
            fill = Celltype,
            shape = dataset) +
        geom_point(alpha = 0.6, size = 2) +
        ggtitle("Batch corrected") +
        scale_shape_manual(values = 21:24) +
        scale_fill_manual(values = cellTypeCol) + 
        guides(fill = guide_legend(override.aes=list(shape=21))) +
        theme_minimal())

####---- Create figure ----####

(fig <- wrap_elements(full = (
    panel3 +
        panel4 +
        plot_layout(guides = "collect") &
        labs(shape = "Data set")
)) +
    (panel1 + 
         labs(fill = "Model variable")) +
    (panel2 +
         guides(shape = "none")) +
    plot_layout(design = "1111
                           2333") +
    plot_annotation(tag_levels = "a"))
ggsave("scripts/figs/integration.pdf", fig, height = 6, width = 7)

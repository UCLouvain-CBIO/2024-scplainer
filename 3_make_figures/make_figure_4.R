#### COMPONENT ANALYSIS

## This script generates figure... where we demonstrate the results of
## differential abundance analysis on the modelled data

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("ggplot2")
library("patchwork")
library("dplyr")
library("tidyr")
library("nipals")
library("scater")

dataDir <- "data/"
figDir <- "figs/"

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

####---- Component analysis ----####

caRes <- scpComponentAnalysis(
    sce,
    ncomp = 20,
    maxiter = 50,
    method = "APCA",
    effect = "SampleType",
    residuals = FALSE
)

####---- t-SNE for unmodelled data vs batch corrected data ----####

pcaUnmodelled <- caRes$bySample$unmodelled
pcaUnmodelled <- pcaUnmodelled[, grepl("^PC", colnames(pcaUnmodelled))]
set.seed(1234)
tsneUnmodelled <-calculateTSNE(t(as.matrix(pcaUnmodelled)))
(panel1 <- data.frame(tsneUnmodelled, colData(sce)) |>
        ggplot() +
        aes(x = TSNE1,
            y = TSNE2,
            colour = Population,
            shape = lcbatch) +
        geom_point(alpha = 0.6) +
        theme_minimal())

apcaSampleType <- caRes$bySample$APCA_SampleType
apcaSampleType <- apcaSampleType[, grepl("^PC", colnames(apcaSampleType))]
set.seed(1234)
tsneSampleType <-calculateTSNE(t(as.matrix(apcaSampleType)))
(panel2 <- data.frame(tsneSampleType, colData(sce)) |>
        ggplot() +
        aes(x = TSNE1,
            y = TSNE2,
            colour = Population,
            shape = lcbatch) +
        geom_point(alpha = 0.6) +
        theme_minimal())

####---- Cell cluster analysis on batch correction ----####

set.seed(2222)
sce$Cluster <- as.factor(kmeans(apcaSampleType, 3)$cluster)
(panel3 <- data.frame(tsneSampleType, colData(sce)) |>
        ggplot() +
        aes(x = TSNE1,
            y = TSNE2,
            colour = Cluster,
            shape = lcbatch) +
        geom_point() +
        theme_minimal())

####---- APCA for sample type ----####

(panel4 <- scpComponentBiplot(
    caRes$bySample |>
        scpAnnotateResults(colData(sce), by = "cell"),
    caRes$byFeature |>
        scpAnnotateResults(rowData(sce), by = "feature", by2 = "Sequence"),
    pointParams = list(
        aes(colour = Population, shape = lcbatch),
        alpha = 0.6
    ),
    textBy = "gene", top = 40
)$APCA_SampleType)

####---- Create figure ----####

populationColorScale <- scale_color_manual(
    values = c(Monocyte = "coral",
               Melanoma = "skyblue3",
               "Undefined Melanoma" = "gray",
               "Resistant Melanoma" = "skyblue"),
    name = "Cell type"
)

(fig <- panel1 +
        ggtitle("Raw (unmodelled) data") +
        guides(shape = FALSE) +
        populationColorScale +
        panel2 +
        ggtitle("Batch corrected data") +
        guides(shape = FALSE) +
        populationColorScale +
        panel3 +
        ggtitle("Cluster analysis") +
        panel4 +
        ggtitle("APCA+ for cell type effects") +
        populationColorScale +
        plot_layout(design = "123
                              444
                              444
                              444",
                    guides = "collect") +
        plot_annotation(tag_levels = "a") &
        labs(shape = "Chromatographic batch"))
ggsave(paste0(figDir, "component.pdf"), fig, height = 10, width = 10)

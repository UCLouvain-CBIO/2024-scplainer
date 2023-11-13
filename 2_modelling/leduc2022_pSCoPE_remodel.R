#### SCP DATA MODELLING

## This script models the leduc2022_pSCoPE dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"
sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))

## Component analysis
caRes <- scpComponentAnalysis(
    sce,
    ncomp = 20,
    maxiter = 50,
    method = "APCA",
    effect = "SampleType",
    residuals = FALSE, 
    unmodelled = FALSE
)
set.seed(11)
apcaRes <- caRes$bySample$APCA_SampleType
apcaRes <- as.matrix(apcaRes[, grep("PC", colnames(apcaRes))])
sce$Cluster <- as.factor(kmeans(apcaRes, 3)$cluster)
sce$Cluster <- dplyr::recode(
    sce$Cluster, "1" = "Melanoma_main", "2" = "Melanoma_sub", "3" = "Monocyte"
)
ggplot(data.frame(tsneSampleType, colData(sce))) +
    aes(x = TSNE1, y = TSNE2, colour = Cluster) +
    geom_point() +
    ggplot(data.frame(tsneSampleType, colData(sce))) +
    aes(x = TSNE1, y = TSNE2, colour = SampleType) +
    geom_point()

####---- Model the data ----####

sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalization
        MedianIntensity +
        ## batch effects
        Channel + Set + 
        ## biological variability
        Cluster,
    name = "model_with_cluster"
)
scpModelFilterPlot(sce, "model_with_cluster")
scpModelFilterThreshold(sce, "model_with_cluster") <- 3
saveRDS(sce, paste0(dataDir, "leduc2022_pSCoPE_re_modelled.rds"))

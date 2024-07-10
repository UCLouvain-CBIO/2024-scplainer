#### SCP DATA MODELLING

## This script models the leduc2022_pSCoPE dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("ggplot2")
library("dplyr")

## data
dataDir <- "data/"
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
apcaRes <- caRes$bySample$APCA_SampleType
apcaRes <- as.matrix(apcaRes[, grep("PC", colnames(apcaRes))])
set.seed(11)
sce$Cluster <- as.factor(kmeans(apcaRes, 3)$cluster)

## Recode cluster with known labels
ggplot(data.frame(apcaRes, colData(sce))) +
    aes(x = PC1, 
        y = PC2, 
        colour = Cluster,
        shape = SampleType) +
    geom_point() 
newLabel <- recode(sce$MelanomaSubCluster, A = "main", B = "sub", .missing = "")
newLabel <- sub("[_]$", "", paste0(sce$SampleType, "_", newLabel))
labelMap <- table(sce$Cluster, newLabel)
clusterLabels <- colnames(labelMap)[apply(labelMap, 1, which.max)]
names(clusterLabels) <- rownames(labelMap)
sce$Cluster <- clusterLabels[sce$Cluster]

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

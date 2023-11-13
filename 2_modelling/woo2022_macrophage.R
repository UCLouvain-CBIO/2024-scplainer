#### SCP DATA MODELLING

## This script models the woo2022_macrophage dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"
woo <- readRDS(paste0(dataDir, "woo2022_macrophage_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(woo, "peptides_log")
sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalization
        MedianIntensity +
        ## batch effects
        Chip +
        ## biological variability
        Treatment)
pl <- scpModelFilterPlot(sce)
scpModelFilterThreshold(sce) <- 2
saveRDS(sce, paste0(dataDir, "woo2022_macrophage_modelled.rds"))

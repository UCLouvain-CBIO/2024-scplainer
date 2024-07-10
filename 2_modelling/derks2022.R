#### SCP DATA MODELLING

## This script models the derks2022 dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "data/"
derks <- readRDS(paste0(dataDir, "derks2022_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(derks, "peptides_log")
sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalization
        MedianIntensity +
        ## batch effects
        Set + Label +
        ## biological variability
        Celltype)
scpModelFilterPlot(sce)
saveRDS(sce, paste0(dataDir, "derks2022_modelled.rds"))

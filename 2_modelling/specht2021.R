#### SCP DATA MODELLING

## This script models the specht2021 dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "data/"
specht <- readRDS(paste0(dataDir, "specht2021_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(specht, "peptides_log")
sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalization
        MedianIntensity +
        ## batch effects
        Channel + Set +
        ## biological variability
        SampleType
)
scpModelFilterPlot(sce)
scpModelFilterThreshold(sce) <- 2
saveRDS(sce, paste0(dataDir, "specht2021_modelled.rds"))

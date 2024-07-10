#### SCP DATA MODELLING

## This script models the schoof2021 dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "data/"
schoof <- readRDS(paste0(dataDir, "schoof2021_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(schoof, "peptides_log")
sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalisation
        MedianIntensity +
        ## batch effects
        Channel + File.ID +
        ## biological variability
        Population
)
scpModelFilterPlot(sce)
scpModelFilterThreshold(sce) <- 3
saveRDS(sce, paste0(dataDir, "schoof2021_modelled.rds"))

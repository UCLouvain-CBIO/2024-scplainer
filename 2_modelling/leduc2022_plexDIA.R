#### SCP DATA MODELLING

## This script models the leduc2022_plexDIA dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"
leduc <- readRDS(paste0(dataDir, "leduc2022_plexDIA_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(leduc, "peptides_log")
sce <- scpModelWorkflow(
    sce,
    formula = ~ 1 + ## intercept
        ## normalization
        MedianIntensity +
        ## batch effects
        Run + Label)
scpModelFilterPlot(sce)
scpModelFilterThreshold(sce) <- 1.25
saveRDS(sce, paste0(dataDir, "leduc2022_plexDIA_modelled.rds"))

#### SCP DATA MODELLING

## This script models the leduc2022_pSCoPE dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")

## data
dataDir <- "~/PhD/asca-scp/scripts/data/"
leduc <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_processed.rds"))

####---- Model the data ----####

sce <- getWithColData(leduc, "peptides_log")
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
scpModelFilterThreshold(sce) <- 3
saveRDS(sce, paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))


#### MINIMAL PROCESSING

## This script performs the minimal processing on the
## leduc2022_plexDIA dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
library("dplyr")
library("camprotR")

## data
leduc <- leduc2022_plexDIA()

## keep only required data
leduc <- removeAssay(leduc, c("peptides", "proteins"))
requiredRowData <- c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id"
)
leduc <- selectRowData(leduc, requiredRowData)

## format missing values
leduc <- zeroIsNA(leduc, i = names(leduc))

####---- Feature quality control ----####

## Remove contaminants
## inspired from https://cambridgecentreforproteomics.github.io/camprotR/articles/crap.html
contaminants <- get_ccp_crap()
contaminants <- paste(contaminants, collapse = "|")
leduc <- filterFeatures(leduc, ~ !grepl(contaminants, Protein.Ids))

####---- Sample quality control ----####

## Number of detected peptides per sample
leduc <- countUniqueFeatures(
    leduc, i = "Ms1Extracted", groupBy = "Stripped.Sequence",
    colDataName = "NumberPeptides"
)
## Median intensity per sample
logInt <- log(assay(leduc[["Ms1Extracted"]]))
MedianIntensity <- colMedians(logInt, na.rm = TRUE)
names(MedianIntensity) <- colnames(leduc)[["Ms1Extracted"]]
colData(leduc)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
## Median CV per sample
leduc <- medianCVperCell(
    leduc, i = "Ms1Extracted", groupBy = "Protein.Group", nobs = 3,
    na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)
## Plot
data.frame(colData(leduc)) %>%
    ggplot() +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = SampleType
    ) +
    geom_point(size = 2) +
    scale_color_continuous(type = "viridis", limits = c(NA, 0.8))
## Filter
leduc$passQC <- leduc$MedianIntensity > 9.5 &
    leduc$SampleType == "Melanoma"
leduc <- subsetByColData(leduc, leduc$passQC)

####---- Building the peptide matrix ----####

## Aggregate PSMs to peptides
leduc <- aggregateFeatures(
    leduc, i = "Ms1Extracted", fcol = "Stripped.Sequence",
    name = "peptides", fun = colMedians, na.rm = TRUE
)

####---- Log-transformation ----####

leduc <- logTransform(leduc, i = "peptides", name = "peptides_log")

####---- Save results ----####

saveRDS(leduc, "data/leduc2022_plexDIA_processed.rds")

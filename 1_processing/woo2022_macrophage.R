
#### MINIMAL PROCESSING

## This script performs the minimal processing on the 
## woo2022_macrophage dataset. 

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
library("patchwork")
library("dplyr")
library("ensembldb")
library("AnnotationHub")

## data
woo <- woo2022_macrophage()

## keep only required data
woo <- woo[, , "peptides_intensity"]
requiredRowData <- c(
    "Sequence", "Leading.razor.protein", "Reverse", 
    "Potential.contaminant"
)
woo <- selectRowData(woo, requiredRowData)

## format missing values
woo <- zeroIsNA(woo, i = 1)

if (FALSE) { ## does not work because of annotation mismatch...
    ## convert protein Uniprot ID to gene names
    proteinIds <- rowData(woo)[[1]]$Leading.razor.protein
    ## mus musculs EnsemblDb annotations (GRCm39) from AnnotationHub, 
    ## latest version at time of writing was v109
    mmus <- query(AnnotationHub(), c("EnsDb", "Mus musculus", "GRCm39", 109))[[1]]
    proteinConversionDf <- transcripts(
        mmus, columns = "gene_name", return.type = "data.frame",
        filter = UniprotFilter(proteinIds, condition = "startsWith")
    )
    matchedIndex <- match(proteinIds, proteinConversionDf$uniprot_id)
    geneName <- proteinConversionDf$gene_name[matchedIndex]
    rowData(woo)[[1]]$gene <- geneName
} 
####---- Feature quality control ----####

## Contaminant plot
df <- data.frame(rowData(woo)[[1]])
df$ContaminantOrReverse <- !(df$Reverse != "+" & 
                                 df$Potential.contaminant != "+" & 
                                 !grepl("REV|CON", df$Leading.razor.protein))
ggplot(df) + 
    aes(x = ContaminantOrReverse) +
    geom_bar()
## filter
woo <- filterFeatures(
    woo, ~ Reverse != "+" & 
        Potential.contaminant != "+" & 
        !grepl("REV|CON", Leading.razor.protein)
)

####---- Sample quality control ----####

## Number of detected peptides per sample
woo <- countUniqueFeatures(
    woo, i = 1, groupBy = "Sequence", colDataName = "NumberPeptides"
)
## Median intensity per sample
woo$MedianIntensity <- colMedians(assay(woo[[1]]), na.rm = TRUE)
## Median CV per sample
woo <- medianCVperCell(
    woo, i = "peptides_intensity", groupBy = "Leading.razor.protein",
    nobs = 5, norm = "SCoPE2", colDataName = "MedianCV"
)
## Plot
data.frame(colData(woo)) %>%
    ggplot() +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = Treatment
    ) +
    geom_point(size = 2) +
    scale_color_continuous(type = "viridis")
## Filter
woo$passQC <- woo$NumberPeptides > 400
woo <- subsetByColData(woo, woo$passQC)

####---- Log-transformation ----####

woo <- logTransform(woo, i = 1, name = "peptides_log")

####---- Save results ----####

saveRDS(woo, "../data/woo2022_macrophage_processed.rds")

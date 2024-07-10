#### MINIMAL PROCESSING

## This script performs the minimal processing on the specht2021
## data set.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
library("patchwork")
library("ensembldb")
library("EnsDb.Hsapiens.v86")
library("dplyr")

## data
specht <- specht2019v3()

## keep only required data
specht <- removeAssay(specht, c("peptides", "proteins"))
requiredRowData <- c(
    "Sequence", "Reverse", "Potential.contaminant",
    "protein", "PIF", "dart_qval"
)
specht <- selectRowData(specht, requiredRowData)

## format missing values
specht <- zeroIsNA(specht, i = names(specht))

####---- Feature quality control ----####

specht <- computeSCR(
    specht, names(specht), colvar = "SampleType",
    samplePattern = "Mel|Macro", carrierPattern = "Carrier",
    sampleFUN = "mean", rowDataName = "MeanSCR"
)
## Contaminant plot
data.frame(rbindRowData(specht, names(specht))) |>
    mutate(ContaminantOrReverse = (Reverse == "+" |
                                       Potential.contaminant == "+" |
                                       grepl("REV|CON", protein))) |>
    ggplot() +
    aes(x = ContaminantOrReverse) +
    geom_bar() +
    ## PIF plot
    ggplot(df) +
    aes(x = PIF) +
    geom_histogram() +
    ## q-value plot
    ggplot(df) +
    aes(x = log10(dart_qval)) +
    geom_histogram() +
    ## mean SCR plot
    ggplot(df) +
    aes(x = log10(MeanSCR)) +
    geom_histogram()
specht <- filterFeatures(
    specht, ~ Reverse != "+" &
        Potential.contaminant != "+" &
        !grepl("REV|CON", protein) &
        !is.na(PIF) & PIF > 0.6 &
        dart_qval < 0.01 &
        !is.na(MeanSCR) & MeanSCR < 0.05
)

####---- Sample quality control ----####

## Number of detected peptides per sample
specht <- countUniqueFeatures(
    specht, i = names(specht), groupBy = "Sequence",
    colDataName = "NumberPeptides"
)
## Median intensity per sample
MedianIntensity <- lapply(experiments(specht), function(x) {
    out <- colMedians(log(assay(x)), na.rm = TRUE)
    names(out) <- colnames(x)
    out
})
names(MedianIntensity) <- NULL
MedianIntensity <- unlist(MedianIntensity)
colData(specht)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
## Median CV per sample
specht <- medianCVperCell(
    specht, i = names(specht), groupBy = "protein",
    nobs = 3, na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)

df <- data.frame(colData(specht))
ggplot(df) +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = SampleType
    ) +
    geom_point(size = 2) +
    # scale_color_continuous(type = "viridis") +
    ggplot(df) +
    aes(
        x = MedianCV,
        fill = SampleType
    ) +
    geom_histogram(bins = 50) +
    plot_layout(ncol = 1, heights = c(0.75, 0.25))

specht$passQC <- !is.na(specht$MedianCV) & specht$MedianCV < 0.32 &
    specht$MedianIntensity > 6 &
    specht$NumberPeptides > 600 &
    grepl("Mono|Macro", specht$SampleType)
specht <- subsetByColData(specht, specht$passQC)
specht <- dropEmptyAssays(specht)

####---- Building the peptide matrix ----####

## Aggregate PSMs to peptides
peptideAssays <- paste0("peptides_", names(specht))
specht <- aggregateFeatures(specht,
                            i = names(specht),
                            fcol = "Sequence",
                            name = peptideAssays,
                            fun = colMedians,
                            na.rm = TRUE)
## Apply majority vote for peptide to protein mapping
ppMap <- rbindRowData(specht, i = grep("^pep", names(specht))) |>
    data.frame() |>
    group_by(Sequence) |>
    ## The majority vote happens here
    mutate(protein = names(
        sort(table(protein), decreasing = TRUE))[1]
    ) |>
    dplyr::select(Sequence, protein) %>%
    dplyr::filter(!duplicated(Sequence, protein))
consensus <- endoapply(rowData(specht)[peptideAssays], function(x) {
    ind <- match(x$Sequence, ppMap$Sequence)
    DataFrame(protein = ppMap$protein[ind])
})
rowData(specht) <- consensus
## Join peptide assays
specht <- joinAssays(specht, i = peptideAssays,
                     name = "peptides")

## Add gene name information
proteinIds <- rowData(specht)[["peptides"]]$protein
proteinConversionDf <- transcripts(
    EnsDb.Hsapiens.v86,
    columns = "gene_name",
    return.type = "data.frame",
    filter = UniprotFilter(proteinIds)
)
matchedIndex <- match(proteinIds, proteinConversionDf$uniprot_id)
geneName <- proteinConversionDf$gene_name[matchedIndex]
rowData(specht)[["peptides"]]$Gene <- geneName

####---- Log-transformation ----####

specht <- logTransform(specht, i = "peptides", name = "peptides_log")

####---- Save results ----####

saveRDS(specht, "data/specht2021_processed.rds")

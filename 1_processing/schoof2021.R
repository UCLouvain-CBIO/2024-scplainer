#### MINIMAL PROCESSING

## This script performs the minimal processing on the schoof2021
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
schoof <- schoof2021()

## keep only required data
schoof <- removeAssay(schoof, c("proteins", "logNormProteins"))
requiredRowData <- c(
    "Annotated.Sequence", "isContaminant",
    "Master.Protein.Accessions", "Isolation.Interference.in.Percent",
    "Percolator.q.Value"
)
schoof <- selectRowData(schoof, requiredRowData)

## format missing values
schoof <- zeroIsNA(schoof, i = names(schoof))

####---- Feature quality control ----####

schoof <- computeSCR(
    schoof, names(schoof), colvar = "SampleType", 
    samplePattern = "sc", carrierPattern = "^boost",
    sampleFUN = "mean", rowDataName = "MeanSCR"
)
## Contaminant plot
df <- data.frame(rbindRowData(schoof, names(schoof)))
ggplot(df) +
    aes(x = isContaminant) +
    geom_bar() +
    ## PIF plot
    ggplot(df) + 
    aes(x = Isolation.Interference.in.Percent) +
    geom_histogram() +
    ## q-value plot
    ggplot(df) + 
    aes(x = log10(Percolator.q.Value)) +
    geom_histogram() +
    ## mean SCR plot
    ggplot(df) + 
    aes(x = log10(MeanSCR)) +
    geom_histogram()

schoof <- filterFeatures(
    schoof, ~ Percolator.q.Value < 0.01 &
        !isContaminant & Master.Protein.Accessions != "" &
        Isolation.Interference.in.Percent < 30 &
        !is.na(MeanSCR) & MeanSCR < 0.2
)

####---- Sample quality control ----####

## Number of detected peptides per sample
schoof <- countUniqueFeatures(
    schoof, i = names(schoof), groupBy = "Annotated.Sequence",
    colDataName = "NumberPeptides"
)
## Median intensity per sample
MedianIntensity <- lapply(experiments(schoof), function(x) {
    out <- colMedians(log(assay(x)), na.rm = TRUE)
    names(out) <- colnames(x)
    out
})
names(MedianIntensity) <- NULL
MedianIntensity <- unlist(MedianIntensity)
colData(schoof)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
## Median CV per sample
schoof <- medianCVperCell(
    schoof, i = names(schoof), groupBy = "Master.Protein.Accessions",
    nobs = 3, na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)

df <- data.frame(colData(schoof))
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

schoof$passQC <- !is.na(schoof$MedianCV) & schoof$MedianCV < 0.4 &
    schoof$MedianIntensity < 3 &
    schoof$NumberPeptides > 1200 &
    schoof$SampleType == "sc"
schoof <- subsetByColData(schoof, schoof$passQC)
schoof <- dropEmptyAssays(schoof)

####---- Building the peptide matrix ----####

## Check PSM to peptide mapping
rd <- rbindRowData(schoof, grep("^F", names(schoof)))
psmsPerPeptide <- sapply(split(rd, rd$assay), function(x) {
    counts <- table(table(x$Annotated.Sequence))
    out <- rep(NA, 10)
    out[1:length(counts)] <- counts
    out
})
psmsPerPeptide <- rowSums(psmsPerPeptide, na.rm = TRUE)
psmsPerPeptide <- psmsPerPeptide[psmsPerPeptide != 0]
ggplot(data.frame(
    nPsms = seq_along(psmsPerPeptide),
    counts = psmsPerPeptide
)) +
    aes(x = nPsms,
        y = counts) +
    geom_bar(stat = "identity") +
    labs(x = "Number PSMs per peptide", 
         title = "schoof2021")

## Aggregate PSMs to peptides
peptideAssays <- paste0("peptides_", names(schoof))
schoof <- aggregateFeatures(schoof,
                            i = names(schoof),
                            fcol = "Annotated.Sequence",
                            name = peptideAssays,
                            fun = colMedians,
                            na.rm = TRUE)
## Join peptide assays 
schoof <- joinAssays(schoof, i = peptideAssays, 
                     name = "peptides")

## Add gene name information
allProts <- rowData(schoof)[["peptides"]][, "Master.Protein.Accessions"]
allProts <- unlist(sapply(allProts, function(x) {
    strsplit(x, "; ")[[1]]
}, USE.NAMES = FALSE))
allProts <- sub("[-]\\d*$", "", allProts)
convert <- transcripts(
    EnsDb.Hsapiens.v86, 
    columns = "gene_name",
    return.type = "data.frame",
    filter = UniprotFilter(allProts)
)
rowData(schoof)[["peptides"]]$Gene <- 
    sapply(rowData(schoof)[["peptides"]]$Master.Protein.Accessions, function(x) {
        out <- strsplit(x, "; ")[[1]]
        out <- sapply(out, function(xx) {
            gene <- convert$gene_name[convert$uniprot_id == sub("[-]\\d*$", "", xx)]
            if (!length(gene)) return(NA)
            unique(gene)
        })
        paste(out, collapse = "; ")
    })
head(rowData(schoof)[["peptides"]][, "Master.Protein.Accessions"])
head(rowData(schoof)[["peptides"]][, "Gene"])

####---- Log-transformation ----####

schoof <- logTransform(schoof, i = "peptides", name = "peptides_log")

####---- Save results ----####

saveRDS(schoof, "../data/schoof2021_processed.rds")

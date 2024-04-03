#### MINIMAL PROCESSING

## This script performs the minimal processing on the
## derks2022 dataset.

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
## BiocManager::install("CambridgeCentreForProteomics/camprotR")
library("camprotR")

## data
derks <- derks2022()
## keep only required data
assaysToRemove <- c(
    grep("119[456]", names(derks), value = TRUE),
    "bulk_prec_extracted", "proteins"
)
derks <- removeAssay(derks, assaysToRemove)
requiredRowData <- c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id"
)
derks <- selectRowData(derks, requiredRowData)
## reformat annotations
derks$Celltype <- sub("_t", "", derks$Celltype)

## format missing values
derks <- zeroIsNA(derks, i = names(derks))

####---- Feature quality control ----####

## Remove contaminants
## inspird from https://cambridgecentreforproteomics.github.io/camprotR/articles/crap.html
contaminants <- get_ccp_crap()
contaminants <- paste(contaminants, collapse = "|")
## remove contaminants (including all keratins)
derks <- filterFeatures(derks, ~ !grepl(contaminants, Protein.Ids) &
                            !grepl("KRT|^K\\d.\\d", Protein.Names))

####---- Sample quality control ----####

sel <- grep("extracted", names(derks), value = TRUE)
## Number of detected peptides per sample
derks <- countUniqueFeatures(
    derks, i = sel,
    groupBy = "Stripped.Sequence", colDataName = "NumberPeptides"
)
## Median intensity per sample
MedianIntensity <- lapply(experiments(derks)[sel], function(x) {
    out <- colMedians(log(assay(x)), na.rm = TRUE)
    names(out) <- colnames(x)
    out
})
names(MedianIntensity) <- NULL
MedianIntensity <- unlist(MedianIntensity)
colData(derks)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
## Median CV per sample
derks <- medianCVperCell(
    derks, i = sel, groupBy = "Protein.Group", nobs = 3,
    na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)
## Q-Exactive
cd <- data.frame(colData(derks))
cd <- dplyr::filter(cd, Instrument == "Q-Exactive")
## Plot
ggplot(cd) +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = Celltype
    ) +
    geom_point(size = 2) +
    scale_color_continuous(type = "viridis")
## Filter
derks$passQC <- derks$Instrument != "Q-Exactive" |
    (derks$MedianIntensity > 9.4 &
         derks$NumberPeptides > 750 &
         derks$Celltype != "Neg")
derks <- subsetByColData(derks, derks$passQC)

## timsTOF-SCP
cd <- data.frame(colData(derks))
cd <- dplyr::filter(cd, Instrument == "timsTOFSCP")
## Plot
ggplot(cd) +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = Celltype
    ) +
    geom_point(size = 2) +
    scale_color_continuous(type = "viridis")
## Filter
derks$passQC <- derks$Instrument != "timsTOFSCP" |
    (derks$MedianIntensity > 6.5 &
         !grepl("DB$", derks$Celltype))
derks <- subsetByColData(derks, derks$passQC)

####---- Building the peptide matrix ----####

## Check precursor to peptide mapping
rd <- rbindRowData(derks, sel)
psmsPerPeptide <- sapply(split(rd, rd$assay), function(x) {
    counts <- table(table(x$Stripped.Sequence))
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
         title = "derks2022")

## Aggregate precursors to peptides
derks <- aggregateFeatures(
    derks, i = sel, fcol = "Stripped.Sequence",
    name = paste0("peptides_", sel),
    fun = colMedians, na.rm = TRUE
)
## Join peptide assays
derks <- joinAssays(
    derks, i = paste0("peptides_", sel), name = "peptides"
)

####---- Log-transformation ----####

derks <- logTransform(derks, i = "peptides", name = "peptides_log")

####---- Save results ----####

saveRDS(derks, "../data/derks2022_processed.rds")

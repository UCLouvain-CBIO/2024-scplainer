#### SCP DATA MODELLING

## This script models the leduc2022_pSCoPE, leduc2022_plexDIA and 
## derks2022 dataset together. This demonstrates that the modelling
## approach can be used to model different types of data (DIA and DDA)
## from different sets of experiments/authors (same lab though).

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("dplyr")
library("ensembldb")
library("EnsDb.Hsapiens.v86")

## data
dataDir <- "data/"
leduc_plexDIA <- readRDS(paste0(dataDir, "leduc2022_plexDIA_processed.rds"))
derks <- readRDS(paste0(dataDir, "derks2022_processed.rds"))

## keep only log peptide assays
leduc_plexDIA <- leduc_plexDIA[, , "peptides_log"]
derks <- derks[, , "peptides_log"]

####---- Combine datasets ----####

## Unify annotation names
colnames(colData(leduc_plexDIA))[[4]] <- "Celltype"
colnames(colData(leduc_plexDIA))[[1]] <- "Set"
leduc_plexDIA$dataset <- "leduc_plexDIA"

derks$dataset <- paste0("derks_", derks$Instrument)
derks$Celltype <- ifelse(derks$Celltype == "U-937", "Monocyte", derks$Celltype)

## rename assays
names(leduc_plexDIA) <- "peptides_leduc_plexDIA"
names(derks) <- "peptides_derks"

## combine data
all <- c(leduc_plexDIA, derks)
all <- joinAssays(all, i = names(all), name = "peptides_log")

## RowData is lost, so need to restore the protein names
pepProtMap <- list(
    rowData(all[["peptides_leduc_plexDIA"]])[, c("Stripped.Sequence", "Protein.Ids")],
    rowData(all[["peptides_derks"]])[, c("Stripped.Sequence", "Protein.Ids")]
)
pepProtMap <- do.call(
    rbind, lapply(pepProtMap, function(x) {
        colnames(x) <- c("Peptide", "Protein")
        x
    })
)
combineProteinNames <- function(prot) {
    prot <- lapply(prot, function(x) strsplit(x, ";")[[1]])
    prot <- sub("^(.*)-.*$", "\\1", unlist(prot))
    paste(sort(unique(prot)), collapse = ";")
}
pepProtMap <- group_by(data.frame(pepProtMap), Peptide) |> 
    summarise(Protein = combineProteinNames(Protein))
ind <- match(rownames(all)[["peptides_log"]], pepProtMap$Peptide)
prots <- pepProtMap$Protein[ind]
rowData(all[["peptides_log"]])$ProteinId <- prots
## Add gene name information
proteinConversionDf <- transcripts(
    EnsDb.Hsapiens.v86, 
    columns = "gene_name",
    return.type = "data.frame",
    filter = UniprotFilter(prots)
)
ind <- match(prots, proteinConversionDf$uniprot_id)
geneName <- proteinConversionDf$gene_name[ind]
rowData(all[["peptides_log"]])$Gene <- geneName

####---- Data modelling ----####

sce <- getWithColData(all, "peptides_log")
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
scpModelFilterThreshold(sce) <- 1.2
saveRDS(sce, paste0(dataDir, "integration_modelled.rds"))

#### VARIANCE ANALYSIS

## This script generates figure... where we demonstrate the results of
## vairance analysis on the modelled data

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("ggplot2")
library("patchwork")
library("dplyr")
library("RColorBrewer")

dataDir <- "data/"
figDir <- "figs/"


####---- Compare sources of variance ----####

## Apply analysis of variance on all modelled data sets
dataFiles <- list.files(
    dataDir, pattern = "modelled.rds", full.names = TRUE
)
dataFiles <- dataFiles[!grepl("integration", dataFiles) &
                           !grepl("re_model", dataFiles)]
vaResAll <- lapply(dataFiles, function(dataFile) {
    dataset <- sub("^.*data//(.*)_mod.*$", "\\1", dataFile)
    cat(dataset)
    x <- scpVarianceAnalysis(readRDS(dataFile))
    x <- scp:::.gatherVarianceData(x)
    x <- group_by(as.data.frame(x), effectName)
    x <- summarise(x, percentExplainedVar = mean(percentExplainedVar))
    x$dataset <- dataset
    gc()
    x
})
vaResCombined <- do.call(rbind, vaResAll)
vaResCombined$effectName <- recode(
    vaResCombined$effectName,
    Channel = "Label",
    Set = "Batch", Run = "Batch", Chip = "Batch",  File.ID = "Batch",
    SampleType = "Biology", Population = "Biology",
    Celltype = "Biology", Treatment = "Biology",
    MedianIntensity = "Normalization Factor"
)
vaResCombined$effectName <- factor(
    vaResCombined$effectName,
    levels = c("Residuals", "Normalization Factor", "Label", "Batch", "Biology")
)
cols <- brewer.pal(n = 8, name = "Set2")
cols <- c("white", cols[c(2, 4, 5, 6)])
names(cols) <- c(
    "Residuals", "Label", "Normalization Factor", "Biology", "Batch"
)
(panel1 <- ggplot(vaResCombined) +
        aes(x = dataset,
            y = percentExplainedVar,
            fill = effectName) +
        geom_bar(stat = "identity", colour = "black") +
        #facet_wrap(~ dataset, nrow = 1) +
        labs(x = "",
             y = "Explained variance") +
        scale_x_discrete(position = "top") +
        scale_fill_manual(values = cols) +
        theme_classic())

####---- Proteins affected by model variables ----####

sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))
vaRes <- scpVarianceAnalysis(sce)
vaRes <- scpAnnotateResults(
    vaRes, rowData(sce)[, c("Sequence", "gene")],
    by = "feature", by2 = "Sequence"
)
## Top 20 peptides affected by sample type
panel2 <- scpVariancePlot(
    vaRes, effect = "SampleType", top = 20, fcol = "gene",
    combined = FALSE
)
## Adapt the labels with more generic names
cols <- brewer.pal(n = 8, name = "Set2")
cols <- c("white", cols[c(2, 4, 5, 6)])
names(cols) <- c(
    "Residuals", "Channel", "MedianIntensity", "SampleType", "Set"
)
(panel2 <- panel2 +
        scale_fill_manual(
            values = cols,
            labels = c("Residuals", "Normalization Factor", "Label", "Batch", "Biology"),
            name = "Descriptor"))

## Top 20 peptides affected by batch effects
(panel3 <- scpVariancePlot(
    vaRes, effect = "Set", top = 20, fcol = "gene", combined = FALSE
))

####---- Peptide with strong biology and batch effect ----####

# IGGSTDTGK is an excellent example of how crucial it is to correctly
# model batch effects!
# Other examples: APNVVVTR, ISGLIYEETR, AEADKNDK, IGGSTDTGK
i <- "APNVVVTR"
prot <- vaRes$Residuals$gene[vaRes$Residuals$feature == i]
df <- data.frame(logIntensity = scpModelInput(sce)[i, ], colData(sce))
(panel4 <- ggplot(df) +
        aes(x = Set,
            y = logIntensity,
            colour = SampleType) +
        geom_boxplot() +
        scale_color_manual(values = c(
            Monocyte = "coral", Melanoma = "skyblue3"
        )) +
        labs(x = "MS acquisition batch",
             y = "Raw log2(Intensity)",
             title = paste0(i, " (", prot, ")"),
             colour = "Cell type") +
        theme_minimal())

(panel5 <- scpVariancePlot(
    endoapply(vaRes, function(x) x[x$feature == i, ]),
    effect = "Set", combined = FALSE
))
panel5$layers[[1]]$geom_params$width <- 0.2

####---- Create figure ----####

(fig <-
     wrap_elements(
         full = panel1 +
             theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90, hjust = 0),
                   axis.line.x = element_blank(),
                   axis.ticks.x = element_blank())
     ) +
     panel2 +
     scale_x_discrete(position = "top") +
     theme(axis.ticks = element_blank(),
           axis.text.x = element_text(angle = 90, hjust = 0),
           strip.text = element_text(angle = 90)) +
     panel4 +
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.line.x = element_blank(),
           panel.grid.major.x = element_blank()) +
     panel5 +
     theme(legend.position = "none",
           axis.text.x = element_text(angle = 0, hjust = 0.5)) +
     plot_layout(design = "
                11112222
                33333334", guides = "collect") +
     plot_annotation(tag_levels = "a"))
ggsave(paste0(figDir, "variance.pdf"), fig, height = 7, width = 10)

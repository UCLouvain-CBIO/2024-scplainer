#### DIFFERENTIAL ABUNDANCE ANALYSIS

## This script generates figure... where we demonstrate the results of
## differential abundance analysis on the modelled data

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("ggplot2")
library("patchwork")
library("dplyr")
library("tidyr")

dataDir <- "data/"
figDir <- "figs/"

sce <- readRDS(paste0(dataDir, "leduc2022_pSCoPE_modelled.rds"))

####---- Differential analysis ----####


daResLeduc <- scpDifferentialAnalysis(
    sce,
    contrasts = list(c("SampleType", "Melanoma", "Monocyte"))
)
daResLeduc <- scpAnnotateResults(
    daResLeduc, rowData(sce)[, c("Sequence", "Leading.razor.protein.id", "gene")],
    by = "feature", by2 = "Sequence"
)
daResLeducProts <- scpDifferentialAggregate(
    daResLeduc, fcol = "Leading.razor.protein.id"
)

####---- Panel: significant vs non-significant ----####

padjPep <- daResLeduc$SampleType_Melanoma_vs_Monocyte$padj
table(padjPep < 0.05)
padjProt <- daResLeducProts$SampleType_Melanoma_vs_Monocyte$padj
table(padjProt < 0.05)
(panel1 <- ggplot(rbind(data.frame(padj = padjPep, type = "peptide"),
                        data.frame(padj = padjProt, type = "protein"))) +
        aes(x = padj < 0.05) +
        geom_bar(stat = "count") +
        facet_grid(~ type ) +
        labs(x = "Significant at 5% FDR") +
        theme_classic())

####---- Panel: volcano plot ----####

## at peptide level
vpl <- scpVolcanoPlot(daResLeduc, textBy = "gene", top = 200)
vpl$SampleType_Melanoma_vs_Monocyte
(panel2 <- vpl$SampleType_Melanoma_vs_Monocyte +
        geom_point(data = data.frame(daResLeduc$SampleType_Melanoma_vs_Monocyte) |>
                       dplyr::filter(!is.na(gene) & gene == "VIM"),
                   colour = "red3") +
        geom_point(data = data.frame(daResLeduc$SampleType_Melanoma_vs_Monocyte) |>
                       dplyr::filter(!is.na(gene) & gene == "CTTN"),
                   colour = "orange2") +
        geom_point(data = data.frame(daResLeduc$SampleType_Melanoma_vs_Monocyte) |>
                       dplyr::filter(!is.na(gene) & gene == "ARHGDIB"),
                   colour = "dodgerblue") +
        geom_point(data = data.frame(daResLeduc$SampleType_Melanoma_vs_Monocyte) |>
                       dplyr::filter(!is.na(gene) & gene == "LCP1"),
                   colour = "purple2") +
        labs(y = "-log10(Adjusted p-value)",
             x = "log2(Fold change)",
             title = "Monocyte <---> Melanoma cell") +
        theme_minimal())

## at protein level
vplProts <- scpVolcanoPlot(daResLeducProts, textBy = "gene", top = 200)
vplProts$SampleType_Melanoma_vs_Monocyte

####---- LFC directions are consistent but influenced by baseline----####

targetProt <- "VIM"
scebc <- scpKeepEffect(sce, "SampleType")
sel <- !is.na(rowData(scebc)$gene) & rowData(scebc)$gene == targetProt
i <- rownames(scebc)[sel]
baseline <- scp:::scpModelIntercept(sce)[i]
lfc <- daResLeduc$SampleType_Melanoma_vs_Monocyte
lfc <- lfc[lfc$feature %in% i, "Estimate"]

logInt <- sweep(assay(scebc)[i, ], 1, baseline, "+")
df <- data.frame(logInt)
df$feature <- factor(
    rownames(df), levels = rownames(df)[order(baseline)]
)
df <- pivot_longer(df, cols = -feature)
df <- cbind(df, colData(scebc)[df$name, ])
(panel3 <- ggplot(df) +
        aes(x = feature,
            y = value) +
        geom_jitter(aes(colour = SampleType), size = 0.25, height = 0) +
        geom_line(data = data.frame(
            feature = i,
            value = baseline
        ), aes(group = 1)) +
        geom_segment(data = data.frame(
            feature = i,
            baseline = baseline,
            lfc = lfc
        ), aes(x = feature,
               xend = feature,
               y = baseline - lfc/2,
               yend = baseline + lfc/2)) +
        scale_color_manual(
            values = c(Monocyte = "coral", Melanoma = "skyblue3")
        ) +
        labs(x = "", y = "log2(Intensity)",
             title = paste("Batch corrected data for", targetProt, "peptides")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)))

####---- Create figure ----####

(fig <- panel1 +
     geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.25) +
     ylim(0, 5250) + ## leave space for geom_text
     panel2 +
     panel3 +
     plot_layout(design = "1222
                           3333") +
     plot_annotation(tag_levels = "a"))
ggsave(paste0(figDir, "differential.pdf"), fig, height = 8, width = 10)

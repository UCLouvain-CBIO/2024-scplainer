#### BENCHMARKING

## This script generates figure... where we compare our approach to 
## concurrent approaches

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
library("patchwork")
library("scater")
library("ggrepel")
## benchmarking utils
source("scripts/3_make_figures/utils.R")

populationColors <- c(
    Monocyte = "coral", 
    Melanoma = "skyblue3",
    "Undefined Melanoma" = "gray",
    "Resistant Melanoma" = "skyblue",
    PDAC = "khaki"
)

####---- Utility functions ----####

benchmarkDataSet <- function(dataSets, bioVars, techVars) {
    stopifnot(!is.null(names(dataSets)))
    out <- mapply(function(d, bioVar, techVar) {
        cat("Benchmarking", d, ":")
        dataSet <- dataSets[[d]]
        res <- benchmarkDataSet(
            dataSet, bioVar, techVar
        )
        cat("\n")
        res[["dataset"]] <- d
        res
    }, names(dataSets), bioVars, techVars)
    do.call(rbind, out)
}

benchmarkMethods <- function(dataSet, bioVar, techVar, normVar,
                             methods = c("None", "Leduc", "HarmonizR_ComBat", "limma", "scplainer")) {
    out <- lapply(methods, function(m) {
        cat(m, ": batch correction...")
        batchCorrectFun <- get(paste0("batchCorrect", m))
        t0 <- Sys.time()
        bc <- batchCorrectFun(dataSet, bioVar, techVar, normVar)
        timing <- Sys.time() - t0
        cat("benchmarking... ")
        res <- benchmarkData(bc, bioVar, techVar)
        res$timing <- timing
        cat("OK\n")
        res
    })
    names(out) <- methods
    out
}

batchCorrectNone <- function(x, bioVar, techVar, normVar) {
    mat <- scpModelInput(x)
    mat <- sweep(mat, 2, x[[normVar]], "-")
    mat2SCE(mat, x)
}

batchCorrectLeduc <- function(x, bioVar, techVar, normVar) {
    suppressMessages(leduc <- leduc2022_pSCoPE())
    sce <- getWithColData(leduc, "proteins_processed")
    sce$Population <- adaptMelanomaAnnotations(sce$SampleType, sce$MelanomaSubCluster)
    sce
}

batchCorrectscplainer <- function(x, bioVar, techVar, normVar) {
    x <- scpModelWorkflow(x, ~ 1 + MedianIntensity + Set + Channel + Population) ## refit model for accurate timing
    scpModelFilterThreshold(x) <- 3
    scpKeepEffect(x, bioVar)
}

batchCorrectHarmonizR_ComBat <- function(x, bioVar, techVar, normVar) {
    require("HarmonizR")
    mat <- scpModelInput(x)
    mat <- sweep(mat, 2, x[[normVar]], "-")
    mat <- as.data.frame(mat)
    design <- data.frame(
        ID = colnames(x), sample = 1:ncol(x), 
        batch = as.factor(x[[techVar]])
    )
    mat <- harmonizR(
        mat, design, algorithm = "ComBat", ComBat_mode = 1, 
        cores = 2, verbosity = 1)
    mat2SCE(as.matrix(mat), x)
}

batchCorrectRfComBat <- function(x, bioVar, techVar, normVar) {
    require("MsCoreUtils")
    mat <- scpModelInput(x)
    mat <- sweep(mat, 2, x[[normVar]], "-")
    mat <- impute_RF(mat, ntree = 100, verbose = TRUE)
    mm <- model.matrix(as.formula(paste("~ 1 +", bioVar)), colData(x))
    mat <- ComBat(mat, batch = x[[techVar]], mod = mm)
    mat2SCE(mat, x)
}

batchCorrectlimma <- function(x, bioVar, techVar, normVar) {
    require("limma")
    mat <- scpModelInput(x)
    techVars <- grep(
        bioVar, all.vars(scpModelFormula(x)), invert = TRUE, 
        value = TRUE
    )
    isFactor <- sapply(techVars, function(i) is.factor(x[[i]]))
    contrs <- lapply(techVars[isFactor], function(i) "contr.sum")
    names(contrs) <- techVars[isFactor]
    techmm <- model.matrix(
        as.formula(paste("~ 1 +", paste(techVars, collapse = " +"))), 
        colData(x), contrasts.arg = contrs
    )
    biomm <- model.matrix(as.formula(paste("~ 1 +", bioVar)), colData(x))
    mat <- removeBatchEffect(mat, covariates = techmm, design = biomm)
    mat2SCE(mat, x)
}

mat2SCE <- function(mat, sce) {
    SingleCellExperiment(
        assays = list(mat),
        colData = colData(sce)[colnames(mat), ],
        rowData = rowData(sce)[rownames(mat), ]
    )
}

adaptMelanomaAnnotations <- function(mainLabel, subLabel) {
    new <- ifelse(
        !is.na(subLabel) & subLabel == "B",
        "Resistant Melanoma", mainLabel
    )
    new <- ifelse(
        is.na(subLabel) & new == "Melanoma",
        "Undefined Melanoma", new
    )
    new
}

####---- Running the benchmark ----####

sce <- readRDS("~/PhD/asca-scp/scripts/data/leduc2022_pSCoPE_modelled.rds")
sce$Population <- adaptMelanomaAnnotations(sce$SampleType, sce$MelanomaSubCluster)
sce <- sce[, sce$Population != "Undefined Melanoma"]
benchmarking <- benchmarkMethods(
    sce, bioVar = "Population", techVar = "Set", 
    normVar = "MedianIntensity"
)
save(benchmarking, file = "scripts/data/benchmarking.rda")

####---- Plot clustering metrics ----####

(metricsFig <- lapply(names(benchmarking), function(n) {
    df <- benchmarking[[n]]$metrics
    df$method <- n
    df
}) |> 
    do.call(what = rbind) |> 
    ggplot() +
    aes(x = biological,
        y = 1 - technical,
        colour = method) +
    geom_point() +
    geom_text_repel(aes(label = method), size = 3) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    theme_minimal() +
    theme(legend.position = "none"))

####---- Plot dimension reductions ----####

wrap_plots(lapply(names(benchmarking), function(n) {
    pcs <- benchmarking[[n]]$pca
    pvar <- round(attr(pcs, "pvar") * 100, 1)
    ggplot(data.frame(pcs, benchmarking[[n]]$annotations)) +
        aes(x = PC1,
            y = PC2,
            colour = Population) +
        geom_point() +
        scale_colour_manual(values = populationColors) +
        labs(x = paste("PC1 (", pvar[[1]], "%)"),
             y = paste("PC2 (", pvar[[2]], "%)"),
             title = n) +
        theme_minimal()
})) +
    plot_layout(guides = "collect")

(dimredFig <- wrap_plots(lapply(names(benchmarking), function(n) {
    tsnes <- benchmarking[[n]]$tsne
    ggplot(data.frame(tsnes, benchmarking[[n]]$annotations)) +
        aes(x = TSNE1,
            y = TSNE2,
            colour = as.factor(Population)) +
        geom_point() +
        labs(x = "t-SNE1", y = "t-SNE2") +
        scale_colour_manual(values = populationColors, name = "Cell Type") +
        ggtitle(n) +
        theme_minimal()
})) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom"))

####---- Create figure ----####

(fig <- metricsFig +
     dimredFig +
     plot_annotation(tag_levels = "a") +
     plot_layout(heights = c(0.25, 0.75)))
ggsave("scripts/figs/benchmark.pdf", fig, height = 8, width = 8)

####---- Supplementary ----####

panels <- list()
for (descriptor in c("Set", "Channel", "lcbatch")) {
    for (n in names(benchmarking)) {
        tsnes <- benchmarking[[n]]$tsne
        pl <- ggplot(data.frame(tsnes, benchmarking[[n]]$annotations)) +
            aes(x = TSNE1,
                y = TSNE2,
                colour = .data[[descriptor]]) +
            geom_point() +
            labs(x = "t-SNE1", y = "t-SNE2") +
            ggtitle(paste("Corrected with", n), 
                    subtitle = paste("Coloured by", descriptor)) +
            theme_minimal()
        panels <- c(panels, list(pl))
    }
}
(suppFig <- wrap_plots(panels) +
    plot_layout(ncol = 3, byrow = FALSE) &
    theme(legend.position = "none"))
ggsave("scripts/figs/supp_benchmark.pdf", suppFig, height = 16, width = 10)

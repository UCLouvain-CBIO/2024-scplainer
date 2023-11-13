library("lisi")
library("aricode")
library("cluster")
library("nipals")
library("scater")

####---- Benchmarking functions ----####

## These function perform the actual benchmarking. The workflow is:
## 1. Dimension reduction
## 2. Clustering: using the K-means clustering
## 3. Assess clustering performance
## 
## Clustering performance is assessed using clustering based on known
## labels. We assess biological labels, for which clustering 
## performance should be maximized, and the technical labels, for 
##  which clustering performance should be minimized.

benchmarkData <- function(dataSet, bioVar, techVar,
                          metrics = c("ASW", "ARI", "NMI", "PS")) {
    pca <- calculateNIPALS(dataSet, bioVar)
    set.seed(1234)
    tsne <- calculateTSNE(t(pca))
    rownames(tsne) <- colnames(dataSet)
    k <- length(unique(dataSet[[bioVar]]))
    dataSet[["cluster"]] <- clusterData(pca, k = k)
    bioMetrics <- computeMetrics(
        pca, dataSet[[bioVar]], dataSet[["cluster"]], metrics
    )
    techMetrics <- computeMetrics(
        pca, dataSet[[techVar]], dataSet[["cluster"]], metrics
    )
    list(
        pca = pca,
        tsne = tsne,
        annotations = colData(dataSet),
        metrics = data.frame(
            metric = metrics, biological = bioMetrics, 
            technical = techMetrics
        )
    )
}

calculateNIPALS <- function(dataSet, bioVar, ncomp = 20, maxiter = 50) {
    if (!is.null(metadata(dataSet)$model)) {
        out <- scpComponentAnalysis(
            dataSet, ncomp = ncomp, maxiter = maxiter, method = "APCA",
            effect = bioVar, residuals = FALSE, unmodelled = FALSE
        )
        out <- out$bySample[[paste0("APCA_", bioVar)]]
        attr(out, "pvar") <- metadata(out)$proportionVariance
        out <- as.matrix(out[, grepl("^PC", colnames(out))])
    } else {
        res <- nipals(
            t(assay(dataSet)), ncomp = ncomp, center = TRUE, 
            maxiter = maxiter, startcol = function(x) sum(!is.na(x)), 
            scale = TRUE
        )
        out <- res$scores %*% diag(res$eig)
        attr(out, "pvar") <- res$R2
        dimnames(out) <- list(colnames(dataSet), paste0("PC", 1:ncomp))
    }
    out
}

computeMetrics <- function(mat, label, cluster, metrics) {
    sapply(metrics, function(m) {
        fun <- get(paste0("compute", m))
        fun(mat, label, cluster)
    })
}

clusterData <- function(mat, k) {
    set.seed(1234)
    clust <- kmeans(mat, centers = k, nstart = 10)
    clust$cluster
}

computeASW <- function(mat, label, cluster) {
    d <- dist(mat)
    sil <- silhouette(as.numeric(as.factor(label)), dist = d)
    mean(sil[, "sil_width"])
}

computeARI <- function(mat, label, cluster) {
    ARI(label, cluster)
}

computeNMI <- function(mat, label, cluster) {
    NMI(label, cluster)
}

computePS <- function(mat, label, cluster) {
    n <- length(label)
    cTable <- table(label, cluster)
    maxInCluster <- colMaxs(cTable)
    sum(maxInCluster) / n
}

computeLISI <- function(mat, label, cluster) {
    lisi <- compute_lisi(mat, data.frame(label = label), "label")
    median(lisi[, 1])
}


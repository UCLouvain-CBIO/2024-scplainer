
## BIOCONDUCTOR
## Pull the development version of the Bioconductor Docker image
FROM bioconductor/bioconductor_docker:devel

## Change home directory
WORKDIR /home/rstudio/

## Install QFeatures, scpdata, and dependencies (2 passes)
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
## Install latest scp version (containing the scplainer implementation) and dependencies (2 passes)
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', dependencies = TRUE)"
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', dependencies = TRUE)"
## Install further dependencies of the project(2 passes)
RUN R -e "BiocManager::install(c('CambridgeCentreForProteomics/camprotR', 'patchwork', 'ensembldb', 'AnnotationHub', 'EnsDb.Hsapiens.v86', 'limma', 'MSCoreUtils', 'scater', 'RColorBrewer', 'tidyr', 'dplyr', 'nipals', 'ggplot2', 'ggrepel', 'immunogenomics/lisi', 'aricode', 'cluster', 'HarmonizR', 'msigdbr', 'fgsea', 'BiocParallel'), dependencies = TRUE)"
RUN R -e "BiocManager::install(c('CambridgeCentreForProteomics/camprotR', 'patchwork', 'ensembldb', 'AnnotationHub', 'EnsDb.Hsapiens.v86', 'limma', 'MSCoreUtils', 'scater', 'RColorBrewer', 'tidyr', 'dplyr', 'nipals', 'ggplot2', 'ggrepel', 'immunogenomics/lisi', 'aricode', 'cluster', 'HarmonizR', 'msigdbr', 'fgsea', 'BiocParallel'), dependencies = TRUE)"

## Download all scpdata data sets used in the project
RUN R -e "scpdata::derks2022()"
RUN R -e "scpdata::leduc2022_pSCoPE()"
RUN R -e "scpdata::leduc2022_plexDIA()"
RUN R -e "scpdata::specht2019v3()"
RUN R -e "scpdata::woo2022_macrophage()"
RUN R -e "scpdata::schoof2021()"

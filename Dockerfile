
## BIOCONDUCTOR
## Pull the development version of the Bioconductor Docker image
FROM bioconductor/bioconductor_docker:devel

## Change home directory
WORKDIR /home/rstudio/

## Clone the scplainer paper repo
RUN git clone https://github.com/UCLouvain-CBIO/2023-scplainer.git .

## Install QFeatures, scpdata, and dependencies (2 passes)
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
## Install latest scp version (containing the scplainer implementation) and dependencies (2 passes)
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', ref = 'scplainer', dependencies = TRUE)"
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', ref = 'scplainer', dependencies = TRUE)"

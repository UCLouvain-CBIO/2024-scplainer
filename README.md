# scplainer: using linear models to understand mass spectrometry-based single-cell proteomics data (source code)

The repository contains the scripts to reproduce the results and
figures of the *scplainer* manuscript.


**scplainer: using linear models to understand mass spectrometry-based single-cell proteomics data**
Christophe Vanderaa, Laurent Gatto *bioRxiv* 2023.12.14.571792;
doi: [https://doi.org/10.1101/2023.12.14.571792](https://doi.org/10.1101/2023.12.14.571792).

## Abstract

Analysing mass spectrometry (MS)-based single-cell proteomics (SCP)
data is challenging. The data analysis must address numerous problems
that are inherent to both MS-based proteomics technologies and
single-cell experiments. This has led to the development of complex
and divergent data processing workflows within the field. In this
work, we present scplainer, a principled and standardised approach
for extracting meaningful insights from SCP data. The approach relies
on minimal data processing combined with linear modelling. The
approach is a simple yet powerful approach for exploring and
interpreting various types of SCP data. scplainer performs variance
analysis, differential abundance analysis and component analysis while
streamlining the visualization of the results. This thorough
exploration enhances our capacity to gain a deeper understanding of
the biological processes hidden in the data. Finally, we demonstrate
that scplainer corrects for technical variability, and even enables
the integration of data sets from different SCP experiments. The
approach effectively generates high-quality data that is amenable to
perform downstream analyses. In conclusion, this work reshapes the
analysis of SCP data by moving efforts from dealing with the technical
aspects of data analysis to focusing on answering biologically
relevant questions.

## Installation instructions

You should install a recent version of R (version >= 4.4.0)
and run the following code chunk to install the necessary
dependencies:

```r
install.packages("BiocManager")
pkgs <- c("QFeatures", "SingleCellExperiment", "scp", "scpdata",
          "ggplot2", "dplyr", "patchwork", "scater")
BiocManager::install(pkgs)
```

For the second part, also run the following command to install the
*scplainer* modelling code:

```r
BiocManager::install("scp")
```

## Useful links

The scplainer approach is implemented in
[scp](https://github.com/UCLouvain-CBIO/scp)

<img
src="https://raw.githubusercontent.com/UCLouvain-CBIO/scp/master/sticker/sticker.png"
height="150">

> Vanderaa, Christophe, and Laurent Gatto. 2021. “Replication of
> Single-Cell Proteomics Data Reveals Important Computational
> Challenges.” Expert Review of Proteomics 18 (10): 835–43.

Data were retrieved from
[scpdata](https://github.com/UCLouvain-CBIO/scpdata)

<img
src="https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/master/sticker/sticker.png"
height="150">

> Vanderaa, Christophe, and Laurent Gatto. 2023. “The Current State of
> Single-Cell Proteomics Data Analysis.” Current Protocols 3 (1):
> e658.

## Repository organisation

- `1_processing` contains R scripts to perform scplainer's minimal
  data processing approach, one script per SCP data set.
- `2_modelling` contains R scripts to perform splainer's modelling
 approach, one script per SCP data set.
- `3_make_figures` contains R scripts to generate the figures from the
  article.

## Reproducing the analysis

You can reproduce the analysis presented in the scplainer paper on
your local machine using `Docker`. You must first install Docker.
Then, pull the image from the
[DockerHub repository](https://hub.docker.com/repository/docker/cvanderaa/scplainer_paper_docker).

```bash
docker pull cvanderaa/scplainer_paper_docker:latest
```

Then, clone the scplainer GitHub repo


## Reproducing analysis through a Rstudio interface

Move your working directory into the cloned github repository, e.g. 
using `cd 2024-scplainer`. 

You can start a Rstudio session within a Docker container using the
following command through your computer terminal:

```bash
docker run -e PASSWORD=bioc -p 8787:8787 -v `pwd`:/home/rstudio/2024-scplainer/ cvanderaa/scplainer_paper_docker:latest
```

Note you should use `%CD%` instead of `pwd` when using Windows.

Open your browser and go to http://localhost:8787. The USER is
`rstudio` and the password is `bioc`. See the [DockerHub
repository](https://hub.docker.com/repository/docker/cvanderaa/scplainer_paper_docker)
for more detailed information on getting started with `Docker`.

Once the Rstudio/R session is opened, you have access to all the
scripts. Any file saved during the session will also be saved on your
local computer.

## Reproducing analysis through the command line

If you don't want to run the analysis in Rstudio through the browser,
you can access the environment from the command line using:

```bash
docker run  -v `pwd`:/home/rstudio/2024-scplainer/ -it cvanderaa/scplainer_paper_docker:latest bash
```

## Note to future self

### Build the docker image

If new dependencies are required, update the `Dockerfile` accordingly.
Build the image (make sure to `cd` in the `2024-scplainer/Docker`
folder):

```bash
docker build -t cvanderaa/scplainer_paper_docker:latest .
```

### Publish the docker image

When complete, push the new image to DockerHub:

```bash
docker push cvanderaa/scplainer_paper_docker:latest
```

### Run complete analysis

Run the following command to fully reproduce the analysis:

```bash
docker run  -v `pwd`:/home/rstudio/2024-scplainer/ -it cvanderaa/scplainer_paper_docker:latest bash

cd 2024-scplainer/
make all
make clean
```

This will automatically update all figures in `figs/` and intermediate
data is stored in `data/`. `make clean` will remove any unnecessary
intermediated files.

## Licence

The source code in this repository is provided under a permissive
[Artistic 2.0 license](https://opensource.org/licenses/Artistic-2.0).
The documentation, including the manual pages and the vignettes, are
distributed under a
[CC BY-SA license](https://creativecommons.org/licenses/by-sa/2.0/).

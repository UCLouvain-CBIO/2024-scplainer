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

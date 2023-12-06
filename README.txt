DimensionReduction folder consists of R code and tables for running PCA, T-SNE, and UMAP

The following files are expression datasets
transcriptome.txt
filtered_transcriptome.txt
proteome.txt

The following files are python scripts:
datafilter.py - filter expression datasets with gene expression values = 0; output filtered_transcriptome.txt
GENIE3.py - GENIE3 function
Run_GENIE3.py - create network using GENIE3 from the filtered_transcriptome.txt or proteome.txt

The following files are graphs to be input into Cytoscape:
TranscriptomicGraph.txt
ProteomicGraph.txt
TwoOmicsGraph_v2.txt

Folders (MixOmics, Proteome, and Transcriptome) contain images and corresponding Cytoscape file

This file (metadata.txt) is just label file to be used in R script

The following files are R scripts that should be able to run using the R project (FinalProject):
DifferentialExpressionTest.R - differential expression analysis
TwoOmics.R - sPLS data integration

Other txt files are output from R scripts, for example:
UvsP_proteome = differential analysis on URA3 vs Prototroph using the proteomic dataset
UvsP_upreg_proteome = log2(fold change) values above 0.2 from the differential analysis

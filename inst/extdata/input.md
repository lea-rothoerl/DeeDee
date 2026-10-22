---
title: "input.md"
author: "Lea Schwarz"
date: "7/12/2021"
output: html_document
---

<span style="color:darkblue"><font face="courier"> <font size="4">***i***</font></font></span> &nbsp;&nbsp;&nbsp;
Upload a DeeDeeExperiment object saved as an .RDS file (created with `saveRDS()`). 
If you started the app with `deedee_app(your_dde)`, the data is already loaded and the upload field is collapsed.

The object needs Ensembl gene IDs as row names and gene symbols in the `SYMBOL` column of its `rowData`. 
Human and mouse data are supported.

**Overview:** The table lists every contrast in your object with the tool that produced it (e.g. DESeq2), the number of genes, and the number of differentially expressed genes in total, up-regulated and down-regulated. The bar plot shows the up- and down-regulated genes per contrast. A gene counts as differentially expressed if its p-value is below the threshold set above the table. The overview always shows the object as submitted.

**Settings:**

* **Show gene symbols:** Switches between gene symbols and Ensembl IDs in plots and tables.
* **Species:** Detected automatically from the Ensembl IDs (`ENSG...` = human, `ENSMUSG...` = mouse). 
It is used for the gene links and the over-representation analysis. 
If the IDs are neither human nor mouse, the scatter plot and the overlap list cannot be created.
* **Contrasts to Analyze** (right side): Untick the contrasts you do not want to use. They are removed from the analysis in all other tabs.

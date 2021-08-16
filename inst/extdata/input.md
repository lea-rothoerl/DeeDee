---
title: "input.md"
author: "Lea Rothörl"
date: "7/12/2021"
output: html_document
---

<span style="color:darkblue"><font face="courier"> <font size="4">***i***</font></font></span> &nbsp;&nbsp;&nbsp;
Upload one or more files. They can hold DeeDee objects or raw DEA result data from DESeq2, limma or edgeR. The raw results need to be supplied in .RDS format, while single DeeDee tables can be supplied as .RDS, .txt or .xlsx files. Lists of DeeDee tables can be provided as .RDS or .xlsx files. 

The .txt files have to be of the following form: The first row is made up of the column names '.logFC' and '.pval', denoting the second and third column in the file. The first column, unnamed in the .txt, holds the gene identifiers as row names. The next two columns comprise the logFC and p-values, respectively. In the .txt, the columns are separated by single spaces. The filename is used as the contrast name in DeeDee.

.xlsx files with one DeeDee table have to consist of one sheet with one column called 'rowname', holding the gene identifiers, and two columns called 'logFC' and `pval', containing the respective values for each gene. In this case, the filename is used as the contrast name. An .xlsx file that comprises a list of multiple DeeDee tables needs to store each table, formatted like above, in a seperate sheet with the sheet names denoting the contrast names.


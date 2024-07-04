# Overview
This folder contains collections of codes free to use within the PMC_Kool group.
It is a continuous Work In Progress, therefore, feedback is highly appreciated :)


Data_fetching: Contains pipeline to fetch bulk RNAsequencing runs from the bioinformatic server of the PMC, it is semi automated, the only thing to be changed are the PMCID of the given samples you want to retrieve in the R script ID_conversion.R (code line 10) and the date and requester in the Data_fetching.sh job script (code line 14 to 16).

DESeq2: Contains pipeline to fetch the bulk RNA sequencing fetched with the 'Data_fetching' pipeline and process it to perform DESeq analysis. The only thing to be changed are the PMCID of the given samples you want to analise in the DESeq2.R script (code line 10) - even though would be nice to have a total automated pipeline, the nature of DESeq analysis varies based on the hypothesis you want to test (different designs), therefore, several minor adjustment might be needed when this analysis is required.
Ideally, one would first run Data_fetching to retrieve the bulk raw counts from the server and then DESeq2 to perform the actual analysis.

Overview: IMPORTANT, overview stores a spreadsheet(best for internal usage for other members of the group) with all the corresponding ID for a sample, don't move it, or else the Data_fetching pipeline won't work. In case of update, replace the spreadsheet with the new one.

remapping: Contains all the necessary elements to perform alignment of single cell sequencing fastq files.

...

Main contributors:
Francesco Valzano
...

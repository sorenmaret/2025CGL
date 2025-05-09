01_Human_Nextflow.Rmd
================
SJM
2025-04-25

# Goal: This is a file to outline everything nessisary to run the RNA-seq on human dox data.

##### We use the RNA-seq pipeline that can be found at the URL <https://nf-co.re/rnaseq/3.18.0/>. Importantly version 3.18.0 was used in our analysis. Instillation instructions can be found at the above URL for the curious. Lets move onto the nessisary files required to run the anaysis.

# Nextflow requires 3 inputs.

##### 1. Raw Sequencing reads in a sample sheet. They need to be in a Fastq format. All reads must be laid out in a sample sheet which organizes these files by number of replicates and strandedness. Our samplesheet can be found at the path below and is called samplesheet.csv:

``` bash
cd /scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_DOX
```

##### 2. Next we need a config file which will configure the analysis into our specific computer and ensures that all resources required for the analysis are availible before and during the analysis. this file is labelled nextflow.config and can be found at the directory below:

``` bash
cd /scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_DOX
```

##### 3. finally we need a shell script which calls up version softwares, performes different opperations in different directories and ensures that nextflow loads with the right version of javascript and outputs files to the desired directories. once again this file GFP.sh can be found at the following directory:

``` bash
cd /scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_DOX
```

##### Conclusion: Given these files we get a large series of analyses and outputs but importantly the rwo files 1. salmon.merged.gene_counts.tsv and 2. salmon.merged.gene_tpm.tsv eable the rest of our analyses and can be found at the directory:

``` bash
cd /scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/HUMAN
```

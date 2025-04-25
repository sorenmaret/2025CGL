
``` r
#BiocManager::install("biomaRt")
library(biomaRt)
```

``` r
# Use the Ensembl mirror
#human_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www")
#mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www")

# Retrieve orthologs
#orthologs <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name"),
                    #filters = "ensembl_gene_id",
                    #values = NULL, # Specify values if needed
                    #mart = human_mart,
                    #attributesL = c("ensembl_gene_id", "external_gene_name"),
                    #martL = mouse_mart)
###!!Error: biomaRt has encountered an unexpected server error.!!###
```

``` r
# download orthologs info from https://www.ensembl.org/biomart/martview
# about how to download, refer https://www.ensembl.org/info/data/biomart/how_to_use_biomart.html'

# load to R
orthologs <- read.csv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/mart_export.txt", header = TRUE, stringsAsFactors = FALSE, sep = ",")
synteny <- read.delim("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/human_mouse_syntenic_lncRNA.txt",
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")
gencode_synteny <- read.delim("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/GENCODEvGENCODE_orthologs_v3.txt",
                      header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# View the first few rows
head(orthologs)
```

    ##    Gene.stable.ID Gene.stable.ID.version Transcript.stable.ID
    ## 1 ENSG00000210049      ENSG00000210049.1      ENST00000387314
    ## 2 ENSG00000211459      ENSG00000211459.2      ENST00000389680
    ## 3 ENSG00000210077      ENSG00000210077.1      ENST00000387342
    ## 4 ENSG00000210082      ENSG00000210082.2      ENST00000387347
    ## 5 ENSG00000209082      ENSG00000209082.1      ENST00000386347
    ## 6 ENSG00000198888      ENSG00000198888.2      ENST00000361390
    ##   Transcript.stable.ID.version Gene.name
    ## 1            ENST00000387314.1     MT-TF
    ## 2            ENST00000389680.2   MT-RNR1
    ## 3            ENST00000387342.1     MT-TV
    ## 4            ENST00000387347.2   MT-RNR2
    ## 5            ENST00000386347.1    MT-TL1
    ## 6            ENST00000361390.2    MT-ND1
    ##                                                                                           Gene.description
    ## 1                              mitochondrially encoded tRNA-Phe (UUU/C) [Source:HGNC Symbol;Acc:HGNC:7481]
    ## 2                                      mitochondrially encoded 12S rRNA [Source:HGNC Symbol;Acc:HGNC:7470]
    ## 3                                mitochondrially encoded tRNA-Val (GUN) [Source:HGNC Symbol;Acc:HGNC:7500]
    ## 4                                      mitochondrially encoded 16S rRNA [Source:HGNC Symbol;Acc:HGNC:7471]
    ## 5                            mitochondrially encoded tRNA-Leu (UUA/G) 1 [Source:HGNC Symbol;Acc:HGNC:7490]
    ## 6 mitochondrially encoded NADH:ubiquinone oxidoreductase core subunit 1 [Source:HGNC Symbol;Acc:HGNC:7455]
    ##   Mouse.gene.stable.ID Mouse.gene.name Mouse.homology.type
    ## 1                                                         
    ## 2                                                         
    ## 3                                                         
    ## 4                                                         
    ## 5                                                         
    ## 6   ENSMUSG00000064341          mt-Nd1    ortholog_one2one
    ##   X.id..target.Mouse.gene.identical.to.query.gene
    ## 1                                              NA
    ## 2                                              NA
    ## 3                                              NA
    ## 4                                              NA
    ## 5                                              NA
    ## 6                                          77.044
    ##   X.id..query.gene.identical.to.target.Mouse.gene
    ## 1                                              NA
    ## 2                                              NA
    ## 3                                              NA
    ## 4                                              NA
    ## 5                                              NA
    ## 6                                          77.044
    ##   Mouse.Gene.order.conservation.score Mouse.Whole.genome.alignment.coverage
    ## 1                                  NA                                    NA
    ## 2                                  NA                                    NA
    ## 3                                  NA                                    NA
    ## 4                                  NA                                    NA
    ## 5                                  NA                                    NA
    ## 6                                  50                                   100
    ##   Mouse.orthology.confidence..0.low..1.high.
    ## 1                                         NA
    ## 2                                         NA
    ## 3                                         NA
    ## 4                                         NA
    ## 5                                         NA
    ## 6                                          1

``` r
colnames(orthologs)
```

    ##  [1] "Gene.stable.ID"                                 
    ##  [2] "Gene.stable.ID.version"                         
    ##  [3] "Transcript.stable.ID"                           
    ##  [4] "Transcript.stable.ID.version"                   
    ##  [5] "Gene.name"                                      
    ##  [6] "Gene.description"                               
    ##  [7] "Mouse.gene.stable.ID"                           
    ##  [8] "Mouse.gene.name"                                
    ##  [9] "Mouse.homology.type"                            
    ## [10] "X.id..target.Mouse.gene.identical.to.query.gene"
    ## [11] "X.id..query.gene.identical.to.target.Mouse.gene"
    ## [12] "Mouse.Gene.order.conservation.score"            
    ## [13] "Mouse.Whole.genome.alignment.coverage"          
    ## [14] "Mouse.orthology.confidence..0.low..1.high."

``` r
head(synteny)
```

    ##        human.gene         mouse.gene Proportion.of.common.genomic.anchors
    ## 1 ENSG00000182873 ENSMUSG00000085063                            0.4821429
    ## 2 ENSG00000182873 ENSMUSG00000098729                            0.6243568
    ## 3 ENSG00000215014 ENSMUSG00000073679                            0.5983936
    ## 4 ENSG00000215014 ENSMUSG00000086549                            0.4779116
    ## 5 ENSG00000215014 ENSMUSG00000086682                            0.5857418
    ## 6 ENSG00000215014 ENSMUSG00000087381                            0.5020080
    ##   Proportion.of.common.protein.coding.genes
    ## 1                                 0.7049180
    ## 2                                 0.8166667
    ## 3                                 0.6833333
    ## 4                                 0.4166667
    ## 5                                 0.7666667
    ## 6                                 0.4333333

``` r
orthologs_filtered <- orthologs[orthologs$Mouse.homology.type == "ortholog_one2one", ]

save(orthologs, orthologs_filtered, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs.RData")

#now, we are ready to find genes are changed in both human and mouse stem cells
```

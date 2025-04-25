02_Finding_conserved_genes.rmd
================
Soren Maret
03.29.25

GOAL: Identify genes that changed in humans and mice stem cells in
response to dox.To do this we want to compare our lists of sig genes
from mouse and humans to a list of orthologues. This is in effect the
same thing as taking the union between the 2 lists.

``` r
# loading realivent dataframes:
# Loading orthologues list
load ("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs.RData")
# Loading human sig genes
load ("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/filtered_res_df_2")
#loading mouse sig genes
load ("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/soma3841/MASTER_CLASS/lessons/06_Differential_expression_analyses/results/DESEQ_results.rdata")
```

First things first we need to make our data compatable. right now if we
try and filter our dataframe we wont get any genes back because of the
decimal values that are present in our sig gene dataframes will block
our downstream analysis.

``` r
# lets make a new dataframe with gene_id and gene_name
mouse_gene_df <- filtered_res_df[, c("gene_id", "gene_name")]

head(mouse_gene_df)
```

    ##                gene_id gene_name
    ## 1 ENSMUSG00000056904.3    Gm5620
    ## 2 ENSMUSG00000058625.6   Gm17383
    ## 3 ENSMUSG00000059751.7    Rps3a3
    ## 4 ENSMUSG00000060198.7   Gm11353
    ## 5 ENSMUSG00000064357.1   mt-Atp6
    ## 6 ENSMUSG00000066270.3   Gm10157

``` r
#we now need to remove the decimal values from our gene_id collumn for orthology analysis
mouse_gene_df <- mouse_gene_df %>%
  mutate(gene_id = sub("\\..*", "", gene_id))
# View result
head(mouse_gene_df)
```

    ##              gene_id gene_name
    ## 1 ENSMUSG00000056904    Gm5620
    ## 2 ENSMUSG00000058625   Gm17383
    ## 3 ENSMUSG00000059751    Rps3a3
    ## 4 ENSMUSG00000060198   Gm11353
    ## 5 ENSMUSG00000064357   mt-Atp6
    ## 6 ENSMUSG00000066270   Gm10157

``` r
#great!! that worked. now lets do the same for our Human d
```

``` r
human_gene_df <- filtered_res_df_2[, c("gene_id", "gene_name")]
#lets check
head(human_gene_df)
```

    ##              gene_id gene_name
    ## 1  ENSG00000001561.7     ENPP4
    ## 2 ENSG00000011201.12     ANOS1
    ## 3 ENSG00000012124.17      CD22
    ## 4 ENSG00000023445.16     BIRC3
    ## 5 ENSG00000041353.10    RAB27B
    ## 6 ENSG00000043355.12      ZIC2

``` r
#we now need to remove the decimal values from our gene_id collumn for orthology analysis
human_gene_df <- human_gene_df %>%
  mutate(gene_id = sub("\\..*", "", gene_id))
# view result
head(human_gene_df)
```

    ##           gene_id gene_name
    ## 1 ENSG00000001561     ENPP4
    ## 2 ENSG00000011201     ANOS1
    ## 3 ENSG00000012124      CD22
    ## 4 ENSG00000023445     BIRC3
    ## 5 ENSG00000041353    RAB27B
    ## 6 ENSG00000043355      ZIC2

``` r
# Filter the orthologue table
conserved_changing_orthologues <- orthologs_filtered %>%
  filter(Gene.stable.ID %in% human_gene_df$gene_id,
         Mouse.gene.stable.ID %in% mouse_gene_df$gene_id)

# View the result
head(conserved_changing_orthologues)
```

    ##    Gene.stable.ID Gene.stable.ID.version Transcript.stable.ID
    ## 1 ENSG00000159231      ENSG00000159231.6      ENST00000290354
    ## 2 ENSG00000151468     ENSG00000151468.11      ENST00000378825
    ## 3 ENSG00000151468     ENSG00000151468.11      ENST00000378839
    ## 4 ENSG00000187569      ENSG00000187569.3      ENST00000345088
    ## 5 ENSG00000167077     ENSG00000167077.13      ENST00000401548
    ## 6 ENSG00000167077     ENSG00000167077.13      ENST00000540833
    ##   Transcript.stable.ID.version Gene.name
    ## 1            ENST00000290354.6      CBR3
    ## 2            ENST00000378825.5     CCDC3
    ## 3            ENST00000378839.1     CCDC3
    ## 4            ENST00000345088.3     DPPA3
    ## 5            ENST00000401548.8      MEI1
    ## 6            ENST00000540833.1      MEI1
    ##                                                                        Gene.description
    ## 1                               carbonyl reductase 3 [Source:HGNC Symbol;Acc:HGNC:1549]
    ## 2                   coiled-coil domain containing 3 [Source:HGNC Symbol;Acc:HGNC:23813]
    ## 3                   coiled-coil domain containing 3 [Source:HGNC Symbol;Acc:HGNC:23813]
    ## 4           developmental pluripotency associated 3 [Source:HGNC Symbol;Acc:HGNC:19199]
    ## 5 meiotic double-stranded break formation protein 1 [Source:HGNC Symbol;Acc:HGNC:28613]
    ## 6 meiotic double-stranded break formation protein 1 [Source:HGNC Symbol;Acc:HGNC:28613]
    ##   Mouse.gene.stable.ID Mouse.gene.name Mouse.homology.type
    ## 1   ENSMUSG00000022947            Cbr3    ortholog_one2one
    ## 2   ENSMUSG00000026676           Ccdc3    ortholog_one2one
    ## 3   ENSMUSG00000026676           Ccdc3    ortholog_one2one
    ## 4   ENSMUSG00000046323           Dppa3    ortholog_one2one
    ## 5   ENSMUSG00000068117            Mei1    ortholog_one2one
    ## 6   ENSMUSG00000068117            Mei1    ortholog_one2one
    ##   X.id..target.Mouse.gene.identical.to.query.gene
    ## 1                                         85.1986
    ## 2                                         90.3704
    ## 3                                         90.3704
    ## 4                                         32.0755
    ## 5                                         80.0628
    ## 6                                         80.0628
    ##   X.id..query.gene.identical.to.target.Mouse.gene
    ## 1                                         85.1986
    ## 2                                         89.3773
    ## 3                                         89.3773
    ## 4                                         34.0000
    ## 5                                         77.0393
    ## 6                                         77.0393
    ##   Mouse.Gene.order.conservation.score Mouse.Whole.genome.alignment.coverage
    ## 1                                 100                                100.00
    ## 2                                 100                                100.00
    ## 3                                 100                                100.00
    ## 4                                 100                                  1.65
    ## 5                                  75                                100.00
    ## 6                                  75                                100.00
    ##   Mouse.orthology.confidence..0.low..1.high.
    ## 1                                          1
    ## 2                                          1
    ## 3                                          1
    ## 4                                          0
    ## 5                                          1
    ## 6                                          1

``` r
#note how our dataframe is small than either of our 2 input dataframes suggesting that our analysis worked. 
```

\#RESULT:

!! We have extracted 276 genes that both signifigantly changed in our
human and mouse samples !!

!! We have several duplicated genes names !!

``` r
#lets go and get rid of our duplicates based on the gene name
# Remove duplicate rows based on "Gene.name"
conserved_changing_orthologues_dereplicated <- conserved_changing_orthologues %>%
  distinct(Gene.name, .keep_all = TRUE)
# WHOA that was a nearly 10x reduction. lets see what these genes are!
cat(conserved_changing_orthologues_dereplicated$Gene.name, sep = "\n")
```

    ## CBR3
    ## CCDC3
    ## DPPA3
    ## MEI1
    ## RBM24
    ## PSTPIP2
    ## DLL1
    ## BCHE
    ## CPM
    ## OTX2
    ## SLC40A1
    ## BIK
    ## PTPRN
    ## GCK
    ## IER3
    ## PPP1R3C
    ## ELAVL4
    ## OAS3
    ## ABCC2
    ## LHX5
    ## PADI3
    ## HTR5A
    ## PSD4
    ## GADL1
    ## FGFR3
    ## CYP1A1
    ## CHAC1
    ## ETV4
    ## RHCG
    ## SP5
    ## ROBO4
    ## ITGA10
    ## LY6H
    ## GAP43
    ## IGF2
    ## ATP1A4

``` r
# copy and paste into https://maayanlab.cloud/Enrichr/ 
```

Lets save our files!

``` r
save(human_gene_df, mouse_gene_df, conserved_changing_orthologues, conserved_changing_orthologues_dereplicated, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs_results.RData")
#checking... 
load("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/MOUSEVSHUMAN/orthologs_results.RData")
```

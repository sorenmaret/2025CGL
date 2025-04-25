Human DESeq
================
ReAnna Owens
3/17/2025

\#Objective: Run DESEQ2 on human data

\#Raw data file path:
/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_DOX

``` r
knitr::opts_chunk$set(echo = TRUE)
library(IRanges)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(readr)
library(ggplot2)
library(purrr)
```

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     reduce

``` r
library(magrittr)
```

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

``` r
library(pheatmap)
library(textshape)
```

    ## 
    ## Attaching package: 'textshape'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(Rcpp)
library(DESeq2)
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'GenomicRanges'

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     subtract

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(tibble)
```

    ## 
    ## Attaching package: 'tibble'

    ## The following object is masked from 'package:textshape':
    ## 
    ##     column_to_rownames

``` r
# We need to make the function write to an output object. To do so we use the 
# <- assisgnement operator, here we are calling the salmon merged gene counts "counts_matrix"
# Then we will call the read.table() function that first requires a file path to the counts in quotes
# Note that in R functions follow with () where the parameter inputs are placed.
# Then we will call up different parameters with a comma, in this case header = TRUE, row.names =1
hDOX_counts <- read.table("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/HUMAN/salmon.merged.gene_counts.tsv", header=TRUE, row.names=1)

# Nice, now we have an object called counts_matrix in our Environment Panel (top right)

# There are time points here that don't exist for mouse data, so these should be removed
hDOX_counts <- hDOX_counts %>%
  select(gene_name, gfp_0_1, gfp_0_2, gfp_0_3, gfp_12_1, gfp_12_2, gfp_12_3, gfp_24_1, gfp_24_2, gfp_24_3, gfp_48_1, gfp_48_2, gfp_48_3, gfp_96_1, gfp_96_2, gfp_96_3)
```

# 2) Making a dataframe of gene_id and gene_name

Now let’s make a table of gene_id and gene_name from the first two
columns of counts_matrix The reason being that the gene_id is a unique
identifying number and the gene_symbol is more intuitive name. We will
use this a lot throughout our analyses.

Note: that the gene_id is currently in a “meta” format so we need to
move it into the table. To do this we will use the rownames function

``` r
# Creating a dataframe with gene_id and gene_symbol
# note creating columns with name = data
g2s <- data.frame(gene_id = rownames(hDOX_counts), gene_symbol = hDOX_counts[ , 1])
#first value is rows, second value is commas in brackets
```

``` r
# The first two columns are gene_id and gene_name
# Those will be very handy let's make an object of these

g2s <- data.frame(
  gene_id = rownames(hDOX_counts),
  gene_name = hDOX_counts[, 1]
)

# Note that rows or columns in black are "meta" meaning not really part of the matrix - just a label of sorts
# Note that there is a column with characters for gene name - we want a pure number matrix, but can have meta labels.
# Let's remove the gene name column for now and we will bring it back later.

# removing gene names by indexing rows and columns via object[ rows, columns ]
hDOX_counts <- hDOX_counts %>% 
  select(-gene_name)

# turning into a matrix using as.matrix() function
hDOX_counts <- as.matrix(hDOX_counts) 

# Rounding the numbers with the round() function
hDOX_counts_rounded <-round(hDOX_counts)

# Note meta nature of cols and rows
```

``` r
# Filtering so Deseq is not getting too many 0's for normalization
hDOX_counts_filtered <- hDOX_counts_rounded[rowSums(hDOX_counts_rounded) > 1, ]

## Check out values of counts_matrix and counts_filtered in Environment panel (top right)
# How many genes were filtered out by this indexing? Hint compare counts_matrix_rounded with counts_filtered
```

``` r
# First we will create a dataframe from counts matrix that includes all the column names
# To do this we call the data.frame() function and inside make col called sample_id
# that is equal to the names of the columns in count matrix (e.g., WT_0_1)
human_deseq_samples <- data.frame(
  sample_id = colnames(hDOX_counts))

# There is good information in each sample_id and we can use code to separate out.
# First we will make an object called "split_values" to split sample_id by any underscore
# Note we are indexing into deseq_samples using '$' to denote index $column_name
# Note we use the function strsplit() to go into sample_id col and split at " _ "
split_values <- strsplit(human_deseq_samples$sample_id, "_")
```

# We just created a new dataframe to input our samples into DESEQ2

Note we used the sample_id to keep track of the original naming This
helps with reproducibility for tracking names across code. Let’s keep
using this practice and get more sample info from sample_id

# 2) Learning Rfuncitons & indexing to retreive name info

``` r
# Now we are going to learn another key function "sapply()"
# This is essentially a built in for loop to go over all rows and perform the same function.
# So here we will go through each row of split_values and run a "generic function(x)" 
# We will then retain the second item which is the time point value in sample_id
time_values <- sapply(split_values, function(x) x[[2]])

# Similar to above we are using sapply to grab the third fragment in split_values
# after two breaks at each "_" there are three parts - here we are grabbing the 3rd 
# NOTE indexing the function of x, which number from the split to a specific number with [[ ]]
replicate_values <- sapply(split_values, function(x) x[[3]])

# Ok now we can do more indexing and make more columns in our deseq_samples dataframe!
# The $ can allow you to index into a specific column or even create a new one.
# so we are making a new column by adding $time_point onto deseq_sample. 
# Then we fill this column by using the assignment operator to place time points in
human_deseq_samples$time_point <- time_values

# Now let's add another column for replicate
human_deseq_samples$replicate <- replicate_values
```

# 3) Factoring columns

``` r
human_deseq_samples$time_point <- factor(human_deseq_samples$time_point)
human_deseq_samples$replicate <- factor(human_deseq_samples$replicate)
```

# 4) Save the sample sheet

``` r
save(human_deseq_samples, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_deseq_samples.RData")
```

# (1) first run DESeq2 by creating a dds (DESEQ Data Structure)

We use the function “DESeqDataSetFromMatrix” With parameters of
countData, colData, design

dds \<- DESeqDataSetFromMatrix(countData = hDOX_counts_filtered, \# this
is our counts data we read in colData = human_deseq_samples, \# telling
DeSeq what env variable to use for sample sheet design = ~ time_point)
\# perhaps most important is “condition” is a factor in deseq_samples \#
we will be using time_point

# Ok let’s make a dds !

``` r
dds <- DESeqDataSetFromMatrix(countData = hDOX_counts_filtered,
                              colData = human_deseq_samples,
                              design = ~ time_point)
```

    ## converting counts to integer mode

# Now run DESeq function that does all the magic !

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# resultsNames function retreives names of results
resultsNames(dds)
```

    ## [1] "Intercept"          "time_point_12_vs_0" "time_point_24_vs_0"
    ## [4] "time_point_48_vs_0" "time_point_96_vs_0"

``` r
# Nice we see this output:
# "Intercept" "time_point_12_vs_0" "time_point_24_vs_0" "time_point_48_vs_0" "time_point_96_vs_0"
# DESEQ found genes differential expression values for each time point relative to 0 time point.
```

# Normalize counts (rlog function)

This basically is rank counts normalized to std error in replicates.

``` r
rlog_counts <- rlog(dds, blind = TRUE)
```

# now we retrieve the values using the “assay” function that converts to rlog_counts

``` r
rlog_counts_matrix <- assay(rlog_counts)
```

# Finally, let’s save this rlog matrix to use in our future analyses and plotting.

We will save this as R Data Structure .rds – this will keep the object
stored properly to be loaded into the environment to use for downstream
analyses.

``` r
write_rds(rlog_counts_matrix, "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/rlog_counts_hDOX.rds")
```

``` r
hDOX_counts_integer <- round(hDOX_counts)
```

``` r
# Check ordering

stopifnot(all(colnames(hDOX_counts_integer) == rownames(human_deseq_samples$sample_id)))

# Nice our columns in the counts are the same as rows in sample_id
```

# 5 Compile all results 0-vs-12, 0-vs-24, 0-vs-48, 0-vs-96

First we are going to make a datframe structure to store the results \#
Note this is a common strategy before a for-loop : make file to populate

``` r
# Let's find the names of the results we want
result_names <- resultsNames(dds)

# Let's get rid of intercept
results_names <- result_names[-1]

# Now the empty data frame with values we want to grab from results

res_df <- data.frame("gene_id" = character(), 
                     "baseMean" = numeric(), 
                     "log2FoldChange" = numeric(), 
                     "lfcSE" = numeric(),
                     "stat" = numeric(),
                     "pvalue" = numeric(),
                     "padj" = numeric(),
                     "gene_name" = character(),
                     "result_name" = character())
```

# 6 FOR-LOOP: for-loop to put the values for each time point analysis

``` r
# Let's figure out each step - for loops are common practice so let's get to know them

for(i in 1:length(results_names)) {
  # grabbing the result name for i in for loop - will repeat for all result names
  results_name <- results_names[i]
  # grabbing the results using DESEQ2 results function for time point i
  res <- results(dds, name = results_name)
  # creating a temporary results data frame in the for loop memory
  tmp_res_df <- res %>% 
    # converting this to dataframe
    as.data.frame() %>%
    # Moving the gene_id into a column now called gene_id
    rownames_to_column("gene_id") %>%
    # Merge in g2s (by gene_id)
    merge(g2s) %>%
    # Add a column to track result name for i
    mutate(result_name = results_name)
  # This will keep adding new results (as rows) for each additional i in for loop
  res_df <- bind_rows(res_df, tmp_res_df)
  
}

# Let's take a look at what we got !
```

# Result: we now have all our deseq2 results in res_df object !

``` r
# Checking NAs
  sum(is.na(res_df$padj))
```

    ## [1] 58683

``` r
  # If we do want to clear them out:
  NoNA <- complete.cases(res_df)
  # WHAAATTTT WHY NO WORK
```

# saving results in res_df

``` r
save(res_df, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/time_point_res_df.RData")
# let's also save our DDS
save(dds, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/dds_time_point_human.RData")
```

# First let’s filter to find genes that significantly changed across time points

For this we will use the filter function in baseR

``` r
filtered_res_df <- res_df %>%
  filter(padj < 0.05)

# Wow we can see a major reduction in our dataframe
```

# Now let’s see what happens when we filter for fold change across time as well

``` r
filtered_res_df_2 <- filtered_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
```

# How many unique genes does this represent?

``` r
human_sig_genes <- unique(filtered_res_df_2$gene_id)

# cool so we have 669 genes that are significantly changing in at least one time point
```

# Result (again keeping track of good results)

# Result: There are X unique genes that change across time (padj \< 0.05 log2 \> 1)

# Lets do gene pathway analysis on the significant genes affected by dox across time

``` r
# printing out the gene_symbol column 
human_genes <- cat(paste(filtered_res_df_2$gene_name, collapse = "\n"))
```

    ## ENPP4
    ## ANOS1
    ## CD22
    ## BIRC3
    ## RAB27B
    ## ZIC2
    ## LMO3
    ## CREB3L3
    ## DLX3
    ## NGFR
    ## SNAP91
    ## FRMPD1
    ## PTGS2
    ## PAG1
    ## JADE1
    ## EDN1
    ## ZCWPW1
    ## MEF2C
    ## ABCB1
    ## WNT11
    ## PHACTR3
    ## GP6
    ## SIRT4
    ## RPH3A
    ## CMTM1
    ## SMPX
    ## NLRP1
    ## SEMA6A
    ## CYP26A1
    ## MISP
    ## BIK
    ## CABP7
    ## MYLK2
    ## HAS3
    ## OCA2
    ## SLC30A4
    ## EBI3
    ## NOD1
    ## PPP1R17
    ## CUBN
    ## VSIR
    ## TLX1
    ## CCL2
    ## UCP1
    ## ZBTB16
    ## FOLR1
    ## IL10RA
    ## NANOG
    ## RBM24
    ## CDH6
    ## CNTN3
    ## HTRA2
    ## DNAH6
    ## IL1R1
    ## CHD5
    ## ACTL8
    ## ALDH8A1
    ## ELOVL3
    ## KCNJ8
    ## PTGFR
    ## RAMP3
    ## DNAI1
    ## DBH
    ## KCNJ2
    ## GPR83
    ## GNMT
    ## MT1G
    ## HS3ST3B1
    ## GRPR
    ## WNK4
    ## HIVEP3
    ## PTPRB
    ## F2RL3
    ## RASL11B
    ## KDR
    ## VPREB3
    ## PODXL
    ## PDE11A
    ## DLL4
    ## KCNA5
    ## HHLA1
    ## PTPRE
    ## ATP1A4
    ## HAPLN2
    ## SYT4
    ## MICAL2
    ## SYT6
    ## NGF
    ## MYCN
    ## CYP2J2
    ## FAM189A2
    ## STX11
    ## CPM
    ## RGS8
    ## RNASEL
    ## CHRND
    ## THSD1
    ## EDNRB
    ## GALNT5
    ## BLK
    ## MYC
    ## IER3
    ## TMPRSS4
    ## TRIM29
    ## ARHGAP20
    ## SQOR
    ## GCOM1
    ## LOXL4
    ## STAT4
    ## CDK15
    ## SLC40A1
    ## PCDH10
    ## GPAT3
    ## DUSP6
    ## SLC46A3
    ## ACVRL1
    ## ZIC5
    ## RHCG
    ## ADAMTS18
    ## IRF8
    ## PRDM16
    ## PADI3
    ## ITGA10
    ## GPA33
    ## VASH2
    ## LYST
    ## ACKR2
    ## ECE2
    ## BHMT
    ## SLC25A48
    ## RASGEF1C
    ## HMGCLL1
    ## FGD2
    ## NLGN4X
    ## CCNB3
    ## PRDM14
    ## ZEB1
    ## HTR7
    ## NPFFR1
    ## ADM
    ## KIAA1755
    ## LYPD1
    ## FOXO1
    ## MAT1A
    ## CCDC3
    ## TMEM45B
    ## PSTPIP2
    ## ZFP36L2
    ## DYNLT5
    ## NR4A2
    ## SCN3A
    ## GRAP
    ## FAM167A
    ## SORBS2
    ## PDE1C
    ## SLC24A2
    ## CFAP70
    ## B3GNT7
    ## NRG1
    ## ERG
    ## SLC30A2
    ## LRRC43
    ## NRG2
    ## NPM2
    ## DMTN
    ## WNT9B
    ## CLIC6
    ## ADGRG5
    ## FCN2
    ## BSND
    ## IL23R
    ## NFIA
    ## ARHGAP25
    ## GABRG1
    ## SPTA1
    ## KIAA0895
    ## GPR85
    ## DLC1
    ## RNF183
    ## LRFN5
    ## RPL10L
    ## EML5
    ## VSTM4
    ## STOX1
    ## SMPD1
    ## CCDC68
    ## MESP1
    ## ANPEP
    ## NOD2
    ## SCARA5
    ## ADAM9
    ## KCTD19
    ## SH3TC2
    ## ADRB2
    ## MT1E
    ## ONECUT1
    ## SERPINA6
    ## DCLK2
    ## TMC7
    ## SLC26A5
    ## SDR16C5
    ## PLAC1
    ## PRKCE
    ## CDK5R2
    ## RXFP1
    ## NMUR1
    ## SPSB1
    ## P2RY6
    ## CYP7B1
    ## ABLIM3
    ## PPP1R3B
    ## TRIB1
    ## CD34
    ## MSRB3
    ## ZCCHC12
    ## PABPC5
    ## FAM241A
    ## CD164L2
    ## VWA3A
    ## P2RY2
    ## PLEKHD1
    ## SYNE3
    ## KCNA2
    ## TGIF1
    ## NAP1L5
    ## RP4-597J3.1
    ## CD19
    ## ERICH5
    ## JUN
    ## ZNRF3-AS1
    ## AMER3
    ## SHISA3
    ## HTR1A
    ## DNAJC22
    ## KLHDC7A
    ## TRIML2
    ## TSHZ1
    ## SERTM1
    ## SKIDA1
    ## SHISA2
    ## HLA-V
    ## PLGLB1
    ## BCOR
    ## TENT5C
    ## LINC00518
    ## B3GALT5
    ## ZNF703
    ## DOCK8-AS1
    ## ADAP2
    ## CLDN5
    ## NCMAP
    ## DHRS7C
    ## DUSP8
    ## TRARG1
    ## TMEM121
    ## ADAMTSL5
    ## CALHM1
    ## FAM86MP
    ## SPDYE17
    ## MITF
    ## SAMD11
    ## SPRY4
    ## C19orf67
    ## TMEM215
    ## CERKL
    ## HMX2
    ## ASTL
    ## OR7E12P
    ## OTUD6A
    ## NCR1
    ## PAX5
    ## SPOCK3
    ## FAT4
    ## AJAP1
    ## HRH1
    ## CTD-2192J16.17
    ## POTEI
    ## CRACDL
    ## FAM163B
    ## ZNF398
    ## RP11-573D15.8
    ## SLC6A17
    ## ADAM32
    ## HMGA2-AS1
    ## AVPR1B
    ## RTP2
    ## RD3
    ## SMOC1
    ## GJC2
    ## SOWAHA
    ## TGM2
    ## SNORA71C
    ## C1orf53
    ## SAMD5
    ## MAFB
    ## OR2H2
    ## MIR1915HG
    ## SPDYC
    ## FAM221B
    ## LINC00654
    ## CCDC85C
    ## RP11-497E19.1
    ## CYS1
    ## HCP5
    ## OR2A42
    ## CSH2
    ## MEF2B
    ## SPDYE3
    ## PLIN5
    ## POLR3DP1
    ## ARHGEF28
    ## RP11-43F13.3
    ## SCRT2
    ## RP11-474L11.5
    ## CELA3B
    ## RP11-3B12.3
    ## LINGO3
    ## TMEM14EP
    ## OR2A1
    ## ROR1-AS1
    ## AC008264.4
    ## RP11-30B1.1
    ## LINC00240
    ## LARGE-AS1
    ## RP11-405O10.2
    ## CT62
    ## AGKP1
    ## LINC01108
    ## LINC01535
    ## RP4-756H11.3
    ## RHEBP2
    ## RP11-365D9.1
    ## LINC01806
    ## FGF7P6
    ## RP11-381K7.1
    ## CFL1P3
    ## RP11-91I20.3
    ## RLIMP1
    ## NKX1-2
    ## MIR205HG
    ## RP11-38L15.3
    ## FAM225A
    ## LINC02595
    ## SHQ1P1
    ## LNCPRESS1
    ## ITGA6-AS1
    ## CYP1B1-AS1
    ## FAM53B-AS1
    ## AP000688.29
    ## AP001442.2
    ## LINC01426
    ## CLDN34
    ## LRRTM4-AS1
    ## LINC00458
    ## RP3-497J21.1
    ## RP11-109I13.2
    ## DNM1P51
    ## SSR4P1
    ## TRDN-AS1
    ## AC000089.3
    ## CPB2-AS1
    ## CLEC2L
    ## RPL5P9
    ## TRHDE-AS1
    ## AC012513.4
    ## ZNF593OS
    ## RPS3AP39
    ## RP4-680D5.2
    ## SHISA9
    ## AP003025.2
    ## LINC01370
    ## RP11-771F20.1
    ## PRR20G
    ## L1TD1
    ## RPL30P12
    ## XACT
    ## RP11-435F17.3
    ## TCP10L
    ## HLA-DMB
    ## MTND5P16
    ## RN7SL189P
    ## OR2A7
    ## ENO1P3
    ## DUXAP10
    ## RP11-260E18.1
    ## RP11-446J8.1
    ## RP11-713M15.1
    ## LINC02275
    ## APELA
    ## SEMA6A-AS1
    ## PCP4L1
    ## LY75-CD302
    ## RP11-131K17.1
    ## LNCPRESS2
    ## RP11-428L21.2
    ## PGAM1P9
    ## CTD-2593A12.3
    ## APOBEC3B-AS1
    ## SERBP1P5
    ## AC008697.1
    ## RP11-114H21.2
    ## CTB-114C7.4
    ## RP11-1391J7.1
    ## RP11-119H12.3
    ## RP11-215P8.4
    ## CTD-2297D10.2
    ## TMEM30BP1
    ## RP11-10L12.2
    ## NAIPP4
    ## CTC-756D1.2
    ## RP11-531A24.5
    ## RP11-30L15.4
    ## LINC01484
    ## PDXDC2P
    ## NANOGP8
    ## DDX18P5
    ## RP5-901A4.1
    ## LINC02700
    ## ABCC6P1
    ## DND1
    ## RP11-234B24.4
    ## PGAM1P5
    ## LINC-ROR
    ## RP11-356K23.1
    ## FOXN3-AS2
    ## NDUFC2-KCTD14
    ## CAP2P1
    ## RP11-209K10.2
    ## AC100830.3
    ## RP11-261B23.1
    ## RP11-923I11.3
    ## RP13-395E19.2
    ## FRRS1L
    ## RP11-297M9.1
    ## RP11-375I20.6
    ## RP11-1006G14.1
    ## SERTM2
    ## RP11-615I2.3
    ## AC009120.5
    ## DLGAP1-AS5
    ## ABCB10P3
    ## LY6L
    ## DLGAP1-AS2
    ## LINC02175
    ## RP1-59D14.1
    ## SPDYE9
    ## U91319.1
    ## RP11-166P13.4
    ## MYZAP
    ## CTD-2541J13.2
    ## PRKCA-AS1
    ## RN7SL648P
    ## BAHCC1
    ## RP11-769O8.1
    ## RP11-214O1.2
    ## RP11-686D22.7
    ## CTB-50L17.14
    ## HOMER3-AS1
    ## MIR222HG
    ## GAS5-AS1
    ## CAHM
    ## QTRT1P1
    ## RP11-348N5.7
    ## SMIM32
    ## RP11-506O24.2
    ## RP4-736L20.3
    ## RP11-506K6.4
    ## MAGI1-IT1
    ## RP11-66B24.7
    ## AGBL1
    ## ABCB10P1
    ## RP11-92F20.1
    ## U2
    ## CTD-2281E23.1
    ## LLNLR-304A6.2
    ## CH507-24F1.2
    ## RIMBP3
    ## RP11-269C23.5
    ## PPP1R26P4
    ## RP11-55J15.2
    ## CTD-2373N4.3
    ## RP11-595O22.1
    ## AC009303.2
    ## RP11-391L3.3
    ## RP11-181K12.1
    ## RP11-573G6.4
    ## RP11-288A5.2
    ## RP11-701H24.10
    ## TMEM75
    ## RP11-562F9.3
    ## RP1-131F15.3
    ## RP11-1110F20.1
    ## RP11-9H20.2
    ## LLNLR-307A6.1
    ## JAKMIP2-AS1
    ## BLACAT1
    ## RP4-590F24.2
    ## LLNLF-176F2.1
    ## RP11-477N12.5
    ## CENPVL2
    ## HSFX4
    ## LMLN2
    ## ATP6V0CP4
    ## PERCC1
    ## RP11-556O5.7
    ## CTC-564N23.6
    ## RP1-166D18.1
    ## RP11-319F12.3
    ## RP11-217L21.1
    ## RP11-495F22.1
    ## RP11-699F21.1
    ## RP11-305M3.6
    ## RP11-71K7.1
    ## RP11-112E16.2
    ## RP11-497D24.1
    ## RP11-547J14.1
    ## CTC-756D1.4
    ## CTC-89C10.1
    ## RP11-6G22.1
    ## RP11-1046B16.4
    ## RP11-416A14.2
    ## RP11-1150C11.2
    ## AC005753.3
    ## RP11-575A19.3
    ## CTD-2008P7.11
    ## RP11-297M9.3
    ## RP13-726E6.5
    ## RP11-325J6.2
    ## RP11-1020A11.3
    ## RP11-358L8.1
    ## RP11-177J6.2
    ## ENPP4
    ## PRSS22
    ## TTC22
    ## ANOS1
    ## YAF2
    ## DPEP1
    ## SNAI2
    ## PLEKHB1
    ## BIRC3
    ## ABCC2
    ## CDH10
    ## RAB27B
    ## TNC
    ## ZIC2
    ## EPHA3
    ## LMO3
    ## METTL24
    ## HHAT
    ## MCOLN3
    ## SOAT1
    ## CREB3L3
    ## SEZ6
    ## DLX3
    ## NGFR
    ## BORCS8-MEF2B
    ## ANKRD44
    ## SNAP91
    ## NAV3
    ## PDZD4
    ## PITX1
    ## NEDD4
    ## SPTB
    ## FGF10
    ## MGAT4A
    ## TRHDE
    ## PTGS2
    ## SEMA3C
    ## ATP12A
    ## KIFAP3
    ## FMO4
    ## PAG1
    ## JADE1
    ## LAMP3
    ## ADCY2
    ## MOXD1
    ## MEF2C
    ## ABCB1
    ## WNT11
    ## FOLH1
    ## GP6
    ## DNMT3B
    ## LHX5
    ## SIRT4
    ## RPH3A
    ## SMPX
    ## NLRP1
    ## SEMA6A
    ## CYP26A1
    ## NRP1
    ## PCDH11Y
    ## SEZ6L
    ## LGALS1
    ## BIK
    ## HRH3
    ## MYLK2
    ## HCK
    ## SLC32A1
    ## RUBCNL
    ## KLF5
    ## TRADD
    ## HAS3
    ## OCA2
    ## SLC30A4
    ## CPQ
    ## NDRG1
    ## SPAG1
    ## RSPH6A
    ## EBI3
    ## CD33
    ## SIGLEC6
    ## UPK1A
    ## RUNDC3B
    ## PON3
    ## NOD1
    ## PPP1R17
    ## OGN
    ## PTGDS
    ## CUBN
    ## VSIR
    ## TLX1
    ## KAZALD1
    ## DNAJC12
    ## DLX4
    ## GABRA4
    ## UCP1
    ## TBC1D9
    ## GLRB
    ## ZBTB16
    ## JHY
    ## FOLR1
    ## IL10RA
    ## MYF6
    ## TRPV4
    ## OAS3
    ## NANOG
    ## HCFC2
    ## PHC1
    ## RBM24
    ## CLIC5
    ## POLR3G
    ## CDH6
    ## STC2
    ## EHHADH
    ## CNTN3
    ## UPK1B
    ## TFCP2L1
    ## CLIP4
    ## HTRA2
    ## DNAH6
    ## EFHD1
    ## IL1R1
    ## IL18R1
    ## SLC5A7
    ## PRRX1
    ## WLS
    ## TFAP2E
    ## NT5C1A
    ## PROX1
    ## NRP2
    ## ALDH8A1
    ## SLC16A7
    ## ABCG2
    ## CCDC92
    ## ELOVL3
    ## PPP1R3C
    ## EGR1
    ## TNFRSF8
    ## KCNJ8
    ## TMEM54
    ## XPNPEP2
    ## PTGFR
    ## INHBA
    ## RAMP3
    ## GLIPR2
    ## NEUROG3
    ## BICC1
    ## PDE1B
    ## RAB9B
    ## CTCFL
    ## IQSEC2
    ## F13A1
    ## GNMT
    ## CDKN1A
    ## CPNE5
    ## NRN1
    ## RAB17
    ## LRRC29
    ## MT1G
    ## TMEM255A
    ## BMP4
    ## MYH2
    ## HS3ST3B1
    ## PSD4
    ## NKX2-4
    ## MMP24
    ## ID1
    ## GRPR
    ## F10
    ## FLRT1
    ## WNK4
    ## OMD
    ## HIVEP3
    ## PTPRB
    ## KLF2
    ## TNFRSF19
    ## STEAP4
    ## KDR
    ## VGF
    ## PODXL
    ## FOXP2
    ## FEZF1
    ## PDE11A
    ## CHAC1
    ## AIPL1
    ## CDH15
    ## UNC13A
    ## GDF15
    ## COX4I2
    ## RHOXF2
    ## PPARG
    ## HHLA1
    ## RGS22
    ## RIN2
    ## VSTM2L
    ## SYT4
    ## MTUS2
    ## ZDHHC8P1
    ## LRRIQ1
    ## SYT6
    ## NGF
    ## MYCN
    ## RERG
    ## ADAMTS8
    ## FAM189A2
    ## B4GALNT1
    ## CPM
    ## RGS8
    ## RNASEL
    ## CHRND
    ## WNT10A
    ## SLC41A2
    ## THSD1
    ## EDNRB
    ## GALNT5
    ## SCN7A
    ## BLK
    ## SKIL
    ## LMX1B
    ## MYC
    ## SLC22A23
    ## IER3
    ## TRIM29
    ## ARHGAP20
    ## SQOR
    ## LOXL4
    ## EXOC6
    ## STAT4
    ## CDK15
    ## SLC40A1
    ## ARHGAP24
    ## HERC6
    ## PCDH10
    ## GPAT3
    ## ENPEP
    ## INHBE
    ## SLC46A3
    ## ACVRL1
    ## KCNH5
    ## DUOXA1
    ## CYP1A1
    ## RHCG
    ## IRF8
    ## SLC13A5
    ## ZNF750
    ## PNMT
    ## SLC6A3
    ## PRDM16
    ## PADI3
    ## KCNH1
    ## VASH2
    ## DUSP10
    ## ACTA1
    ## LYST
    ## FBLN7
    ## MARCHF4
    ## ALDH1L1
    ## TM4SF19
    ## ECE2
    ## ADAMTS16
    ## UGT3A1
    ## BHMT
    ## CXCL14
    ## SLC25A48
    ## GFRA3
    ## HMGCLL1
    ## PRSS35
    ## SLC22A3
    ## NLGN4X
    ## CCNB3
    ## HTR2C
    ## SNTG1
    ## PRDM14
    ## TERF1
    ## COMMD3
    ## ZEB1
    ## LRMDA
    ## NPFFR1
    ## KIRREL3
    ## TMEM25
    ## KIAA1755
    ## VEGFC
    ## PPP1R1C
    ## IL18
    ## FOXO1
    ## SLC7A11
    ## CCDC3
    ## DRD3
    ## AKR1C2
    ## VENTX
    ## TMEM45B
    ## FBXO4
    ## GFRA1
    ## PSTPIP2
    ## ZFP36L2
    ## SPEF2
    ## DYNLT5
    ## HHEX
    ## NR4A2
    ## GRAP
    ## NRGN
    ## ANGPT1
    ## C4orf19
    ## FAM167A
    ## GPR26
    ## SORBS2
    ## PDE1C
    ## WNT7A
    ## SLC24A2
    ## CFAP70
    ## GNA14
    ## GDF6
    ## PPP2R2B
    ## B3GNT7
    ## TIMP4
    ## NRG1
    ## FAM81A
    ## PWWP3B
    ## ACAN
    ## SLC30A2
    ## GRHL3
    ## LRRC43
    ## RIBC1
    ## NPM2
    ## VWA5B1
    ## DMTN
    ## WNT9B
    ## CLIC6
    ## CLDN14
    ## UROC1
    ## TPPP3
    ## AGPAT3
    ## S100B
    ## FCN2
    ## KASH5
    ## IP6K3
    ## HS3ST6
    ## FGF19
    ## ELAVL4
    ## BSND
    ## TMEM82
    ## IL23R
    ## KLHDC8A
    ## LRRTM1
    ## BBS5
    ## ANKRD23
    ## ARHGAP25
    ## GABRG1
    ## CCDC141
    ## SPTA1
    ## ADAMTS9
    ## TM4SF18
    ## ENPP6
    ## TXLNB
    ## IL31RA
    ## KIAA0895
    ## GPR85
    ## FNDC1
    ## PGAM2
    ## DLC1
    ## TMC1
    ## ALDH1A1
    ## TRPV6
    ## CFAP47
    ## RNF183
    ## PHYHIPL
    ## GJB2
    ## RPL10L
    ## VSTM4
    ## SPACA9
    ## STOX1
    ## PPP1R36
    ## CMTM5
    ## SMPD1
    ## WDR72
    ## CCDC68
    ## GALR1
    ## CFAP52
    ## HTR3A
    ## SLFN5
    ## MESP1
    ## ANPEP
    ## EVA1C
    ## MEI1
    ## NOD2
    ## IGF2
    ## GPD1
    ## CYP2S1
    ## TMC4
    ## KRT80
    ## TMEM88
    ## DNASE1L2
    ## SCARA5
    ## ADAM9
    ## KCTD19
    ## CXXC4
    ## KLK7
    ## SH3TC2
    ## ADRB2
    ## RASSF6
    ## MT1E
    ## ONECUT1
    ## P2RY1
    ## FRMPD4
    ## SERPINA6
    ## DCLK2
    ## KISS1
    ## SLC26A5
    ## SDR16C5
    ## PLAC1
    ## KCNG3
    ## PRKCE
    ## NPTX1
    ## CDK5R2
    ## RXFP1
    ## SPSB1
    ## P2RY6
    ## UTF1
    ## MIR31HG
    ## CD8B
    ## SLFN12
    ## MAB21L4
    ## KLHL6
    ## SLFN11
    ## CYP7B1
    ## SH3RF3
    ## DELEC1
    ## ABLIM3
    ## TRIB1
    ## TAS1R1
    ## HES3
    ## CMB9-22P13.1
    ## CD34
    ## MSRB3
    ## C11orf45
    ## ZCCHC12
    ## FGFBP3
    ## PABPC5
    ## BTC
    ## GPR171
    ## VWA3A
    ## P2RY2
    ## FOSL1
    ## RP11-276H1.3
    ## PLEKHD1
    ## SYNE3
    ## GCNT4
    ## LY6H
    ## KCNA2
    ## NAP1L5
    ## RP4-597J3.1
    ## CD19
    ## ERICH5
    ## IRX3
    ## JUN
    ## ZNRF3-AS1
    ## GRAMD1C
    ## SLITRK1
    ## SHISA3
    ## DNAJC22
    ## ERBB4
    ## KLHDC7A
    ## TRIML2
    ## GPC5
    ## ALOX12B
    ## TSHZ1
    ## ITPRID1
    ## SERTM1
    ## CHRM4
    ## SHISA2
    ## SAGE1
    ## ZNF467
    ## OR52A1
    ## GABRG3
    ## NLRP10
    ## AP1S2
    ## IGIP
    ## PCP4
    ## TMEM119
    ## PLGLB1
    ## SYNDIG1L
    ## TENT5C
    ## COL18A1-AS1
    ## LINC00518
    ## DOCK8-AS1
    ## ADAP2
    ## CLDN5
    ## PCDH9
    ## GDF3
    ## SLIT3
    ## NCMAP
    ## DHRS7C
    ## DUSP8
    ## TRARG1
    ## PURA
    ## METTL7A
    ## AHNAK2
    ## ACTBP7
    ## POU3F1
    ## ADAMTSL5
    ## RTN4RL1
    ## CALHM1
    ## KLHL32
    ## FAM86MP
    ## CLCNKA
    ## DTX2P1
    ## MITF
    ## MT1X
    ## MAGI2
    ## C11orf96
    ## SAMD11
    ## SPRY4
    ## HNRNPA1P61
    ## MEIS3P2
    ## C19orf67
    ## HES4
    ## SELL
    ## CERKL
    ## SERPINA5
    ## FBLL1
    ## NOXA1
    ## ENTPD8
    ## LRRK2
    ## AADACL3
    ## DPY19L2P1
    ## C15orf61
    ## OR7E12P
    ## OTUD6A
    ## IL1RAP
    ## PAX5
    ## SPOCK3
    ## TDRD7
    ## AKR1C3
    ## FAT4
    ## KIF19
    ## ZNF781
    ## S100A5
    ## CRYBA4
    ## AJAP1
    ## HRH1
    ## POU3F4
    ## CRACDL
    ## SCOC-AS1
    ## ZNF398
    ## SLC6A17
    ## ADAM32
    ## HMGA2-AS1
    ## CYP2B6
    ## DNAH10
    ## MAP1LC3C
    ## MB
    ## RTP2
    ## RD3
    ## ARC
    ## DLL1
    ## ANKRD13B
    ## SMOC1
    ## STYXL2
    ## CES1
    ## SOWAHA
    ## TGM2
    ## RORB
    ## SNORA71C
    ## SERTAD4-AS1
    ## C1orf53
    ## SAMD5
    ## ECEL1P1
    ## C9orf129
    ## CDSN
    ## RP11-304F15.3
    ## MIR1915HG
    ## SPDYC
    ## C9orf135
    ## LINC02397
    ## LINC00654
    ## KRT17P5
    ## ADGRG1
    ## MT1H
    ## CCDC85C
    ## RP11-497E19.1
    ## CYS1
    ## OR2A42
    ## CRIP1
    ## CCDC183
    ## CSH2
    ## PRDX1P1
    ## CTAGE14P
    ## LRRC37A11P
    ## GOT2P6
    ## COL28A1
    ## TCEA1P3
    ## RP11-43F13.3
    ## SCRT2
    ## LINC01356
    ## RPS16P5
    ## CELA3B
    ## RP11-3B12.3
    ## LINGO3
    ## HMSD
    ## APOL6
    ## OR2A1
    ## RPL17P17
    ## RP11-216N14.7
    ## UHRF2P1
    ## ROR1-AS1
    ## LINC01767
    ## LINC00240
    ## LARGE-AS1
    ## OR52K3P
    ## AC073957.15
    ## RP5-997D24.3
    ## LINC01771
    ## FAM225B
    ## RP3-495K2.3
    ## RP11-375H17.1
    ## VN2R10P
    ## TERF1P5
    ## RP3-462C17.1
    ## AGKP1
    ## LINC01108
    ## RP11-365D9.1
    ## LINC01806
    ## RP1-213J1P__B.1
    ## RP11-381K7.1
    ## FAM66A
    ## CFL1P3
    ## OR2M3
    ## MIR34AHG
    ## RP11-108M9.3
    ## RP11-91I20.3
    ## VIM-AS1
    ## LINC02069
    ## RLIMP1
    ## NKX1-2
    ## ZRANB2-AS2
    ## RP1-20B11.2
    ## LINC00102
    ## RP4-753M9.1
    ## MIR205HG
    ## RP11-38L15.3
    ## LINC02884
    ## PCED1CP
    ## RPL13AP17
    ## RP1-225E12.2
    ## HLA-DPA1
    ## RP1-45C12.1
    ## FAM225A
    ## LINC02595
    ## SHQ1P1
    ## LUARIS
    ## AP002856.5
    ## IMPDH1P10
    ## LNCPRESS1
    ## AC021876.4
    ## AC002480.3
    ## RPL23AP50
    ## AP000688.29
    ## AC005083.1
    ## TAS2R62P
    ## RP11-115A15.2
    ## RHBDF1P1
    ## TMEM229A
    ## RP11-146I2.1
    ## LINC01426
    ## RP11-561O23.5
    ## LINC00458
    ## RP11-119H12.4
    ## RPS3AP28
    ## RP11-374P20.4
    ## DNM1P51
    ## SSR4P1
    ## AC009299.2
    ## MSH2-OT1
    ## AC000089.3
    ## CPB2-AS1
    ## CHCHD4P3
    ## TRHDE-AS1
    ## RP11-525G13.2
    ## RP11-69I8.2
    ## PRMT5-AS1
    ## RPS3AP39
    ## RP4-680D5.2
    ## FOXD2-AS1
    ## SHISA9
    ## AP003025.2
    ## CT75
    ## DCTN1-AS1
    ## LINC01370
    ## RP11-400K9.4
    ## PLCH1-AS1
    ## RPS14P8
    ## PRR20G
    ## L1TD1
    ## RPL30P12
    ## SEPTIN7P6
    ## RPL17P19
    ## MCRIP2P1
    ## XACT
    ## ANOS2P
    ## HLA-DMB
    ## RPS3AP34
    ## RPS3AP43
    ## LINC00882
    ## AC005062.2
    ## MTND5P16
    ## RN7SL189P
    ## RPS3AP42
    ## DUXAP10
    ## RAD21L1
    ## MALINC1
    ## LINC02466
    ## RP11-629G13.1
    ## CARS1-AS1
    ## RP11-446J8.1
    ## ANK2-AS1
    ## RP11-713M15.1
    ## LINC02275
    ## GAPDHP61
    ## SEMA6A-AS1
    ## PCP4L1
    ## RP11-131K17.1
    ## LINC02432
    ## RP11-304F15.4
    ## LNCPRESS2
    ## RP11-428L21.2
    ## PGAM1P9
    ## CTD-2593A12.3
    ## APOBEC3B-AS1
    ## RP11-893F2.5
    ## PVT1
    ## TMEM158
    ## SLC7A11-AS1
    ## DNAH10OS
    ## RHOQ-AS1
    ## CTB-114C7.4
    ## CTD-2213F21.2
    ## RP11-1391J7.1
    ## RP11-631M6.2
    ## RP11-119H12.3
    ## RP11-215P8.4
    ## CTD-2297D10.2
    ## SELENOP
    ## TMEM30BP1
    ## TMEM108-AS1
    ## LINC02714
    ## DDX18P4
    ## LINC02482
    ## LINC02754
    ## BRCC3P1
    ## CTC-575N7.1
    ## MIR124-1HG
    ## RP11-697N18.3
    ## CTC-756D1.2
    ## CTD-2501M5.1
    ## RP11-531A24.5
    ## LINC01484
    ## MAP2K1P1
    ## NPIPB11
    ## LRRC24
    ## RP11-522I20.3
    ## RP11-660L16.2
    ## EGLN1P1
    ## STX16-NPEPL1
    ## RP11-326C3.2
    ## RP11-567M21.3
    ## PDXDC2P
    ## NANOGP8
    ## NOX5
    ## CTC-342M10.2
    ## LINC02700
    ## ABCC6P1
    ## DND1
    ## LINC02361
    ## RP11-17G12.2
    ## RP11-234B24.4
    ## PGAM1P5
    ## LINC01490
    ## RP11-478C19.2
    ## RP11-240B13.2
    ## LINC-ROR
    ## RP11-356K23.1
    ## FOXN3-AS2
    ## NDUFC2-KCTD14
    ## CAP2P1
    ## CTD-2054N24.2
    ## RP11-37C7.2
    ## AC100830.3
    ## DNM1P47
    ## RP11-227D13.1
    ## THSD4-AS1
    ## RP11-261B23.1
    ## RP13-395E19.2
    ## PPIAP51
    ## FRRS1L
    ## RP11-817O13.8
    ## RP11-297M9.1
    ## CTD-2542L18.1
    ## SERTM2
    ## AC009120.5
    ## KB-1836B5.1
    ## RP11-6O2.4
    ## RP4-545K15.5
    ## RP11-104N10.1
    ## LINC01003
    ## DLGAP1-AS5
    ## ABCB10P3
    ## AC005592.3
    ## DLGAP1-AS2
    ## U91319.1
    ## RP11-166P13.4
    ## MYZAP
    ## MIR3648-2
    ## PRKCA-AS1
    ## RN7SL648P
    ## BAHCC1
    ## RP11-214O1.3
    ## RP11-214O1.2
    ## RP11-686D22.7
    ## PIN1-DT
    ## RP11-47L3.1
    ## CTB-50L17.14
    ## LINC01901
    ## CTD-2537I9.12
    ## RP11-687F6.4
    ## CTC-559E9.2
    ## RP4-806M20.3
    ## DMRTC1
    ## CTC-339O9.1
    ## RP11-463C8.7
    ## RP11-932O9.8
    ## GAS5-AS1
    ## CAHM
    ## RP11-420L9.5
    ## QTRT1P1
    ## YPEL5P2
    ## RP11-506O24.2
    ## BCAN-AS1
    ## RP11-506K6.4
    ## TTC23L-AS1
    ## CD24
    ## RP11-440L14.4
    ## RP3-508I15.21
    ## RP11-674N23.4
    ## SNURF
    ## RP11-38M8.1
    ## RP11-138A9.2
    ## AGBL1
    ## U2
    ## RP11-7F17.8
    ## RP11-670E13.6
    ## LINC02334
    ## U2
    ## RP11-92F20.1
    ## U2
    ## RP11-395B7.2
    ## CTD-2281E23.1
    ## CH507-24F1.2
    ## RP5-836N17.4
    ## RIMBP3
    ## DACH1
    ## CTD-2024F21.1
    ## KB-68A7.2
    ## PPP1R26P4
    ## RDM1
    ## U2
    ## RP11-55J15.2
    ## RP13-895J2.4
    ## AC005753.1
    ## RP11-264M12.4
    ## RP11-395N21.2
    ## RP11-595O22.1
    ## RP11-861E23.2
    ## RP3-326L13.2
    ## RP11-826F13.1
    ## RP11-391L3.3
    ## RP11-181K12.1
    ## RP11-368N21.5
    ## RP11-528M18.2
    ## RP11-1334A24.5
    ## RP4-671O14.5
    ## TMEM75
    ## RP1-131F15.3
    ## RP11-463O12.5
    ## CTD-2012K14.5
    ## RP11-81A22.4
    ## RP11-1110F20.1
    ## RP11-20G6.1
    ## RP11-20I23.10
    ## LLNLR-307A6.1
    ## JAKMIP2-AS1
    ## RP11-485F13.1
    ## RP11-259O18.5
    ## RP3-467D16.3
    ## RP4-590F24.2
    ## RP11-477N12.5
    ## C19orf85
    ## RP11-434E6.5
    ## AC116366.7
    ## CTD-2010I16.2
    ## RP11-606E8.2
    ## RP11-511P7.5
    ## RP11-674E16.5
    ## CTC-564N23.6
    ## RP11-36B6.2
    ## RP11-384L19.1
    ## TOP1P1
    ## RP11-325K19.7
    ## RP1-166D18.1
    ## RP11-388A16.2
    ## RP5-914M6.2
    ## RP11-401F2.5
    ## RP11-217L21.1
    ## COPG2IT1
    ## RP11-314O13.2
    ## RP11-359N11.3
    ## RP11-342I1.3
    ## RP11-268J15.6
    ## RP11-1150C11.1
    ## CTB-88F18.6
    ## RP11-699F21.1
    ## RP11-112E16.2
    ## RP11-567P19.3
    ## RP11-497D24.1
    ## RP11-547J14.1
    ## CTC-756D1.4
    ## CTC-89C10.1
    ## RP11-6G22.1
    ## RP13-726E6.4
    ## RP11-506H21.6
    ## CH17-262O2.5
    ## RP11-1150C11.2
    ## AC005753.3
    ## RP11-575A19.3
    ## RP11-32P3.1
    ## RP11-372J12.1
    ## RP11-155D18.17
    ## RP13-726E6.5
    ## RP11-791I24.5
    ## RP11-113H14.3
    ## RP11-1020A11.3
    ## RP5-1055C14.9
    ## RP11-177J6.2
    ## TNMD
    ## MATK
    ## MLXIPL
    ## ETV7
    ## SLC7A14
    ## CLDN11
    ## GPRC5A
    ## YAF2
    ## RGPD5
    ## IL20RA
    ## VSIG2
    ## RUNX3
    ## PLEKHB1
    ## ABCC2
    ## VIM
    ## CD44
    ## AGPAT4
    ## DAPK2
    ## TLL1
    ## CDH10
    ## TNC
    ## FAM184B
    ## LMO3
    ## METTL24
    ## HHAT
    ## MCOLN3
    ## SOAT1
    ## GNA15
    ## CREB3L3
    ## COL11A1
    ## WNT8A
    ## SEZ6
    ## DLX3
    ## BORCS8-MEF2B
    ## SNAP91
    ## NAV3
    ## PDZD4
    ## FGFR3
    ## PITX1
    ## GAL
    ## SPTB
    ## FRMPD1
    ## TRHDE
    ## PTGS2
    ## SEMA3C
    ## FGF4
    ## PAG1
    ## MBNL3
    ## SPAG6
    ## MAP2
    ## ADCY2
    ## AFP
    ## MEF2C
    ## SMARCD3
    ## MAPRE3
    ## BCORL1
    ## FOLH1
    ## SEPHS1
    ## LHX5
    ## SIRT4
    ## RPH3A
    ## CMTM1
    ## KCNH4
    ## CHRD
    ## ABCC6
    ## SMPX
    ## NLRP1
    ## ZFHX4
    ## CCDC80
    ## DAZL
    ## SEMA6A
    ## TBC1D2
    ## CYP26A1
    ## PCDH11Y
    ## SEZ6L
    ## LGALS1
    ## CRYBB1
    ## BIK
    ## RASD2
    ## TRIM9
    ## CHGA
    ## SLC8A3
    ## DHRS2
    ## SLA2
    ## HRH3
    ## SLC17A9
    ## NKAIN4
    ## MYLK2
    ## PDYN
    ## MAP1LC3A
    ## ZNF516
    ## CDH20
    ## RNF125
    ## NOL4
    ## CACNA1F
    ## SYP
    ## TSNAXIP1
    ## CBLN1
    ## NECAB2
    ## CORO7-PAM16
    ## OCA2
    ## SLC30A4
    ## NEFM
    ## KCNN4
    ## FCER2
    ## RSPH6A
    ## TLE6
    ## CCNP
    ## EBI3
    ## SIGLEC6
    ## MAG
    ## RASA4
    ## PON3
    ## VIPR2
    ## NOD1
    ## PPP1R17
    ## MOGAT3
    ## RARRES2
    ## MYL7
    ## GCK
    ## DOCK8
    ## ELAVL2
    ## CUBN
    ## VSIR
    ## TLX1
    ## KAZALD1
    ## RAI1
    ## ICAM2
    ## DLX4
    ## HLF
    ## GLRB
    ## ZBTB16
    ## JHY
    ## FOLR1
    ## IL10RA
    ## SLC1A2
    ## ASIC1
    ## IL23A
    ## MYF6
    ## TRPV4
    ## OAS3
    ## NANOG
    ## HCFC2
    ## PHC1
    ## RBM24
    ## KHDRBS2
    ## EYA4
    ## CLIC5
    ## NME5
    ## LOX
    ## PCDHB7
    ## POLR3G
    ## STC2
    ## EHHADH
    ## CNTN3
    ## FGF12
    ## UPK1B
    ## OTOF
    ## CLIP4
    ## HTRA2
    ## GALNT3
    ## FN1
    ## DNAH6
    ## IGFBP5
    ## OTX1
    ## IL18R1
    ## MLPH
    ## KISS1R
    ## ATP6V1B1
    ## PARD3B
    ## PRRX1
    ## OPRD1
    ## WLS
    ## SIPA1L2
    ## CD48
    ## PROX1
    ## ESYT2
    ## B4GALT6
    ## FILIP1
    ## CCN2
    ## ABCG2
    ## SPP1
    ## PCDH17
    ## ECRG4
    ## ONECUT2
    ## PROX2
    ## ESRRB
    ## CNRIP1
    ## PYROXD2
    ## PRLHR
    ## MYCT1
    ## ENOX1
    ## ADRA1A
    ## TNFRSF8
    ## TAS2R10
    ## HS3ST2
    ## NXPH1
    ## INHBA
    ## FKBP9
    ## GLIPR2
    ## NEUROG3
    ## EGR2
    ## BHLHE41
    ## ARHGAP9
    ## PDE1B
    ## NFE2
    ## SEMG1
    ## IQSEC2
    ## XG
    ## CDKN1A
    ## SOX4
    ## NRN1
    ## ATXN1
    ## RAB17
    ## AHNAK
    ## SH3TC1
    ## LRRC29
    ## MT1G
    ## MT2A
    ## SOX21
    ## TMEM255A
    ## BMP4
    ## PTGER2
    ## MYH2
    ## HS3ST3B1
    ## OPRL1
    ## PSD4
    ## FOXA2
    ## NKX2-4
    ## FLRT3
    ## GRPR
    ## F10
    ## NR1D1
    ## WNK4
    ## IFI6
    ## HIVEP3
    ## MASP1
    ## PTPRB
    ## TNFRSF19
    ## SPINK2
    ## KDR
    ## RFPL3
    ## FOXP2
    ## FEZF1
    ## OPN1SW
    ## PDE11A
    ## CHAC1
    ## AIPL1
    ## CDKN1C
    ## CDH15
    ## SHC2
    ## THEM6
    ## GADD45G
    ## BST2
    ## MAS1
    ## EPO
    ## UNC13A
    ## GDF15
    ## IQCN
    ## CALY
    ## PNCK
    ## NUTM2F
    ## COX4I2
    ## GJA9
    ## ABHD12B
    ## MATN3
    ## HHLA1
    ## IQCA1
    ## RGS22
    ## RIN2
    ## VSTM2L
    ## SYT4
    ## MTUS2
    ## DCLK1
    ## CCNA1
    ## IRS4
    ## NTS
    ## LRRIQ1
    ## MICAL2
    ## MRO
    ## MEIS2
    ## CRB1
    ## SOX5
    ## SOX3
    ## PHC2
    ## DSC3
    ## FAM189A2
    ## AGAP2
    ## B4GALNT1
    ## HEY2
    ## CPM
    ## RGS8
    ## RNASEL
    ## CHRND
    ## WNT10A
    ## ALDH1L2
    ## THSD1
    ## EDNRB
    ## DGKB
    ## GALNT5
    ## SCN7A
    ## BLK
    ## GATA4
    ## SKIL
    ## SLC31A2
    ## GABBR2
    ## MYC
    ## IER3
    ## MYO7A
    ## SULF1
    ## TRIM29
    ## SQOR
    ## TMEM62
    ## SEMA6D
    ## LOXL4
    ## CUZD1
    ## EXOC6
    ## STAT4
    ## CDK15
    ## SLC40A1
    ## PCDH10
    ## LEF1
    ## PRICKLE1
    ## COL2A1
    ## INHBE
    ## LGR5
    ## TMEM132B
    ## FAM222A
    ## SLC46A3
    ## ACVRL1
    ## SSTR1
    ## KCNH5
    ## DUOXA1
    ## DUOX2
    ## NTRK3
    ## SLC5A2
    ## MYLK3
    ## ADAMTS18
    ## UNC45B
    ## RHBDL3
    ## GREB1L
    ## SLC13A5
    ## PNMT
    ## STAC2
    ## NLRP12
    ## PADI3
    ## FNDC7
    ## FCRL5
    ## VASH2
    ## SMYD2
    ## DUSP10
    ## ACTA1
    ## LYST
    ## CNIH3
    ## SYT2
    ## ANKRD53
    ## ST6GAL2
    ## THNSL2
    ## FBLN7
    ## CDCA7
    ## FAM171B
    ## MARCHF4
    ## GADL1
    ## ACKR2
    ## STAC
    ## IGSF11
    ## FAM86B2
    ## ECE2
    ## PDGFC
    ## ADAMTS16
    ## UGT3A1
    ## HAPLN1
    ## IQGAP2
    ## CXCL14
    ## SLC25A48
    ## GFRA3
    ## DCDC2
    ## HMGCLL1
    ## SLC22A3
    ## CCNB3
    ## HTR2C
    ## MAGEA4
    ## DOK2
    ## SNTG1
    ## RSPO2
    ## SHC3
    ## COMMD3
    ## ZEB1
    ## NPFFR1
    ## SLC5A12
    ## PAMR1
    ## SLC43A1
    ## NCAM1
    ## TMEM25
    ## TAGLN
    ## KIAA1755
    ## GPM6A
    ## IL18
    ## DRD3
    ## AKR1C2
    ## VENTX
    ## TMEM45B
    ## ZNF385D
    ## FBXO4
    ## OBI1
    ## RIT2
    ## GRIA4
    ## SPEF2
    ## SPARCL1
    ## CAPSL
    ## RASGRP3
    ## DYNLT5
    ## MR1
    ## ACOXL
    ## PTPRR
    ## NR4A2
    ## PLA2R1
    ## SCN3A
    ## C16orf89
    ## GRAP
    ## AK5
    ## NRGN
    ## ANGPT1
    ## C4orf19
    ## FAM167A
    ## C10orf90
    ## SORBS2
    ## NCAM2
    ## PDE1C
    ## SLFN13
    ## PIEZO2
    ## CA10
    ## SAMSN1
    ## ADCY8
    ## CFAP70
    ## GNA14
    ## SORCS3
    ## FGF18
    ## PPP2R2B
    ## LRFN2
    ## B3GNT7
    ## NRG1
    ## HTR5A
    ## KCNJ6
    ## ERG
    ## SLC34A2
    ## GRHL3
    ## RHBDL2
    ## AUTS2
    ## RIBC1
    ## HTR6
    ## NPM2
    ## VWA5B1
    ## DMTN
    ## CLIC6
    ## CBR3
    ## TPPP3
    ## AGPAT3
    ## FCN2
    ## IKZF3
    ## KASH5
    ## JOSD2
    ## MEIOB
    ## FGF19
    ## BSND
    ## TMEM82
    ## MXRA8
    ## IL23R
    ## NBPF20
    ## CAPN13
    ## LRRTM1
    ## KCNF1
    ## VSNL1
    ## CCDC74A
    ## BBS5
    ## ARHGAP25
    ## ALPG
    ## PROK2
    ## CCDC141
    ## DPPA2
    ## SPTA1
    ## ADAMTS9
    ## HESX1
    ## PRRT3
    ## CXCL5
    ## TM4SF18
    ## LIPH
    ## ZDHHC19
    ## HPGD
    ## TXLNB
    ## IL31RA
    ## KIAA0895
    ## GPR85
    ## SHH
    ## FNDC1
    ## PGAM2
    ## SLC30A8
    ## SYK
    ## ABCA1
    ## C8orf34
    ## ALDH1A1
    ## TRPV6
    ## CFAP47
    ## FAT3
    ## LRFN5
    ## PKNOX2
    ## RPL10L
    ## DEPP1
    ## OTX2
    ## DACT1
    ## UCMA
    ## VSTM4
    ## STOX1
    ## CMTM5
    ## FBN1
    ## AVPR1A
    ## ASPG
    ## GABRB3
    ## SMPD1
    ## PRTG
    ## A2ML1
    ## GALR1
    ## ANPEP
    ## MAP1A
    ## IGF2
    ## CA4
    ## NIBAN3
    ## RHEBL1
    ## TMC4
    ## IGFBP6
    ## TBX10
    ## TMEM88
    ## PLAAT5
    ## PTGDR
    ## GNG4
    ## BMP1
    ## ADAM9
    ## GDNF
    ## LDLRAD4
    ## KCTD19
    ## CA7
    ## CXXC4
    ## PPIC
    ## PXDC1
    ## KLK7
    ## RSPO1
    ## SHE
    ## ZEB2
    ## CLIC3
    ## MT1E
    ## LDB2
    ## ONECUT1
    ## P2RY1
    ## AHSP
    ## PYDC1
    ## FRMPD4
    ## FAM153A
    ## SERPINA6
    ## FOS
    ## DCLK2
    ## KISS1
    ## SLC26A5
    ## RASA4B
    ## SDR16C5
    ## TEX13B
    ## PLAC1
    ## ALK
    ## KCNG3
    ## PRKCE
    ## NPTX1
    ## LPAR3
    ## DSCAM
    ## UTF1
    ## FBXL14
    ## MIR31HG
    ## GAP43
    ## MOB3A
    ## SLFN12
    ## PDE3A
    ## EFEMP2
    ## MOS
    ## SLFN11
    ## CYP7B1
    ## CES3
    ## ABLIM3
    ## INSM1
    ## TAS1R1
    ## GPR160
    ## MSRB3
    ## ZHX3
    ## MIR1-1HG-AS1
    ## TRHR
    ## ZCCHC12
    ## LEP
    ## PABPC5
    ## BTC
    ## RTP1
    ## DES
    ## INHBC
    ## VWA3A
    ## PCSK1
    ## LPL
    ## FOSL1
    ## RP11-276H1.3
    ## ETV4
    ## CREG2
    ## DOK7
    ## SYNE3
    ## LRRC37A3
    ## TYMSOS
    ## LY6H
    ## FIBIN
    ## ZDHHC22
    ## OLIG3
    ## IRX3
    ## ST8SIA3
    ## ZNF491
    ## CTD-2530N21.5
    ## FAM87B
    ## LCN15
    ## ODF3B
    ## GRAMD1C
    ## RNF212
    ## TMEM151B
    ## SLITRK1
    ## HTR1A
    ## DNAJC22
    ## CPNE7
    ## KLHDC7A
    ## EXOC3L1
    ## HES7
    ## EGR3
    ## LBX2
    ## ZBTB42
    ## PCED1B
    ## LRRC3B
    ## PHC1P1
    ## NAIPP2
    ## TSHZ1
    ## C3orf80
    ## OR3A1
    ## FAR2P1
    ## MYLPF
    ## LINC00304
    ## SHISA2
    ## SAGE1
    ## ZNF467
    ## TIGIT
    ## OR52A1
    ## UNC5C
    ## ASB18
    ## NLRP10
    ## TSPAN10
    ## IGIP
    ## RGS6
    ## SLC35D3
    ## GALNT9
    ## GLUD2
    ## SOX1
    ## PCP4
    ## GAS6
    ## CSMD1
    ## TMEM119
    ## CALN1
    ## C2CD4C
    ## PLGLB1
    ## SYNDIG1L
    ## LINC00518
    ## EFHC2
    ## NOG
    ## TMPRSS2
    ## GOLGA6L4
    ## PCDH9
    ## GDF3
    ## TMEM255B
    ## DHRS7C
    ## DUSP8
    ## PDE4B
    ## KCNH7
    ## RALYL
    ## TRARG1
    ## SRY
    ## EFCAB10
    ## MANEAL
    ## PURA
    ## LINC00313
    ## SOCS1
    ## HS6ST3
    ## METTL7A
    ## DLK1
    ## INKA1
    ## POU3F1
    ## KCNQ5
    ## ADAMTSL5
    ## IFITM1
    ## CCIN
    ## FAM86MP
    ## GABRA5
    ## MST1P2
    ## TRDN
    ## MIR22HG
    ## MITF
    ## AKR1C1
    ## VSTM2B
    ## FOXD3
    ## MT1X
    ## DPPA3
    ## SPRY4
    ## C2orf78
    ## FYB2
    ## LDLRAD2
    ## C19orf67
    ## WNT7B
    ## HES4
    ## SELL
    ## SERPINA5
    ## FBLL1
    ## VWC2
    ## OPTC
    ## ADRB3
    ## ENTPD8
    ## AADACL3
    ## KCTD21
    ## RIMKLBP1
    ## IL1RAPL2
    ## ANKRD34B
    ## PLAC9
    ## PCDH18
    ## DPY19L2P1
    ## OR7E12P
    ## SYCP2
    ## SPOCK3
    ## TDRD7
    ## AKR1C3
    ## SPATS2L
    ## KIF19
    ## ACADSB
    ## ZNF781
    ## PRTN3
    ## CRYBA4
    ## FGF16
    ## SPTSSB
    ## AJAP1
    ## HRH1
    ## SLC30A10
    ## CTD-2192J16.17
    ## CRACDL
    ## LAMB3
    ## TMEM26
    ## ZNF398
    ## ADAM32
    ## PIWIL2
    ## ZNF860
    ## CYP2B6
    ## PDGFA
    ## SSPOP
    ## MYH6
    ## DNAH10
    ## CR1L
    ## CTC-301O7.4
    ## HES5
    ## KEL
    ## MB
    ## MT1F
    ## RTP2
    ## ARMH1
    ## C2CD4A
    ## RD3
    ## LPA
    ## DLL1
    ## ANKRD13B
    ## ALPK2
    ## GRM3
    ## GJC2
    ## CES1
    ## TGM2
    ## RORB
    ## SNORA71C
    ## C1orf53
    ## SAMD5
    ## RIPPLY2
    ## TDRG1
    ## MAFB
    ## ECEL1P1
    ## HLA-DMA
    ## PRRT1
    ## C9orf129
    ## ACOXL-AS1
    ## LINC01123
    ## FOXB2
    ## C9orf135
    ## SMPD5
    ## PKHD1L1
    ## LINC02397
    ## KRT17P5
    ## MT1H
    ## CCDC85C
    ## RP11-497E19.1
    ## RP9P
    ## CYS1
    ## TTC23L
    ## LINC02656
    ## OR2A42
    ## FAM24B
    ## CSH2
    ## PTPRCAP
    ## KRT222
    ## SLC25A14P1
    ## RP5-849H19.2
    ## FAM209B
    ## RPL3P6
    ## PRDX1P1
    ## RPL18P10
    ## KANSL1-AS1
    ## PLIN5
    ## LRRC37A11P
    ## COL28A1
    ## TCEA1P3
    ## BASP1-AS1
    ## RP11-43F13.3
    ## SCRT2
    ## ZC3H11B
    ## ATAD3C
    ## RPL22P12
    ## CEBPZOS
    ## CELA3B
    ## RP11-3B12.3
    ## LINGO3
    ## HTR5A-AS1
    ## HMSD
    ## OR2A1
    ## AC005481.5
    ## AC004471.9
    ## RP11-545E17.3
    ## LINP1
    ## CDYLP1
    ## RPL36P4
    ## RORB-AS1
    ## LINC00863
    ## RP11-82L18.2
    ## LARGE-AS1
    ## XXbac-BPG55C20.7
    ## SLC9A3-AS1
    ## ZNRF2P2
    ## RP11-415J8.3
    ## AC104134.2
    ## LINC01107
    ## AC007679.3
    ## KCNMA1-AS3
    ## DBH-AS1
    ## FGF13-AS1
    ## TERF1P5
    ## RP11-526D8.7
    ## LINC01108
    ## RPL15P14
    ## RP11-365D9.1
    ## LINC01806
    ## RPS3AP37
    ## SOX21-AS1
    ## INKA2-AS1
    ## RP11-381K7.1
    ## FAM66A
    ## AC007383.3
    ## OR2M3
    ## MIR34AHG
    ## RP5-968D22.3
    ## LINC00954
    ## RP11-91I20.3
    ## RP11-404F10.2
    ## OR2A9P
    ## VIM-AS1
    ## RP11-168O16.1
    ## CICP3
    ## RP11-274B21.12
    ## LINC02069
    ## NKX1-2
    ## ZRANB2-AS2
    ## LRRC37A15P
    ## MTCO1P42
    ## KARS1P2
    ## AC159540.1
    ## AC112229.1
    ## ANKRD63
    ## FOXD3-AS1
    ## SRGAP2-AS1
    ## RPS6P15
    ## MIR205HG
    ## RP11-38L15.3
    ## LINC02884
    ## PCED1CP
    ## BEND3P1
    ## HLA-DPA1
    ## RP1-45C12.1
    ## PTCHD3P2
    ## AP002856.5
    ## IMPDH1P10
    ## KIF26B-AS1
    ## LNCPRESS1
    ## EMSLR
    ## SNX18P2
    ## CYP1B1-AS1
    ## RASA4DP
    ## AP000688.29
    ## RP11-115A15.2
    ## TMEM229A
    ## LINC01426
    ## HMGN2P5
    ## DNM1P51
    ## TRDN-AS1
    ## AC009299.2
    ## LINC00607
    ## AC000089.3
    ## BHLHE40-AS1
    ## RMDN2-AS1
    ## CPB2-AS1
    ## H2BW3P
    ## FBXL12P1
    ## CHCHD4P3
    ## TRHDE-AS1
    ## LINC02653
    ## VLDLR-AS1
    ## RP11-89N17.4
    ## RP11-69I8.2
    ## RP11-390F4.6
    ## AC068134.6
    ## RPS3AP39
    ## TMEM26-AS1
    ## SHISA9
    ## LINC01275
    ## OXCT2P1
    ## LINC01370
    ## RP6-24A23.3
    ## LINC01410
    ## RP11-384F7.2
    ## RN7SL674P
    ## RPL30P12
    ## ARHGEF25
    ## RP11-372E1.1
    ## SEPTIN7P6
    ## RPL17P19
    ## MCRIP2P1
    ## AF011889.5
    ## RP11-435F17.3
    ## ARHGDIG
    ## HLA-DMB
    ## KRT8P35
    ## AC005062.2
    ## CFB
    ## MTND5P16
    ## RN7SL189P
    ## ENO1P3
    ## ARHGEF35-AS1
    ## DUXAP10
    ## OR2A1-AS1
    ## RAD21L1
    ## MALINC1
    ## LINC00461
    ## RP11-629G13.1
    ## RP11-446J8.1
    ## CICP16
    ## LINC02275
    ## SEMA6A-AS1
    ## GDNF-AS1
    ## OXCT1-AS1
    ## RP11-432B6.3
    ## C4orf54
    ## R3HDM2P1
    ## C8orf34-AS1
    ## RP11-304F15.4
    ## LNCPRESS2
    ## SEMA6A-AS2
    ## PGAM1P9
    ## CTD-2593A12.3
    ## APOBEC3B-AS1
    ## MTCO1P24
    ## MIR302CHG
    ## TMEM158
    ## SLC7A11-AS1
    ## DNAH10OS
    ## RP11-9G1.3
    ## CTB-114C7.4
    ## RP11-305O6.3
    ## TRMT9B
    ## RP11-215P8.4
    ## CTD-2297D10.2
    ## TMEM30BP1
    ## LINC01843
    ## SEMA5A-AS1
    ## RP11-586D19.2
    ## NAIPP4
    ## RP11-48B3.5
    ## CTD-2501M5.1
    ## LINC01484
    ## ETV3L
    ## RP11-48B3.3
    ## RPL7P20
    ## RP11-267L5.1
    ## CTB-180C19.1
    ## LRRC24
    ## AP000662.4
    ## RP11-660L16.2
    ## CKLF-CMTM1
    ## INMT-MINDY4
    ## RP11-326C3.2
    ## FAM66D
    ## RP11-567M21.3
    ## GLYATL1P1
    ## NANOGP8
    ## CTC-342M10.2
    ## TBX5-AS1
    ## LINC02700
    ## DND1
    ## RP11-73M18.2
    ## AP003419.11
    ## RP13-895J2.3
    ## RP11-114G22.1
    ## GPR142
    ## LINC02293
    ## RP11-328C8.4
    ## ZNF578
    ## RP11-240B13.2
    ## LINC-ROR
    ## LINC00638
    ## AC007040.11
    ## RP11-371E8.4
    ## NDUFC2-KCTD14
    ## RP11-209K10.2
    ## RP11-82L7.4
    ## RP11-3D4.2
    ## CTD-2054N24.2
    ## GH1
    ## RP11-37C7.2
    ## DNM1P47
    ## CSPG4P11
    ## HERC2P10
    ## RP11-261B23.1
    ## RP13-395E19.2
    ## FRRS1L
    ## RP11-297M9.1
    ## CTD-2541J13.1
    ## RP11-375I20.6
    ## CTA-29F11.1
    ## SERTM2
    ## LINC00622
    ## SPECC1-DT
    ## RP11-6O2.4
    ## LINC02473
    ## RBFADN
    ## LINC02188
    ## FIGNL2
    ## TMEM249
    ## CTC-523E23.1
    ## TCF24
    ## AC144831.1
    ## DKFZP434A062
    ## TMEM256-PLSCR3
    ## U91319.1
    ## RP11-621L6.2
    ## ZNF516-AS1
    ## BAHCC1
    ## RP11-640I15.1
    ## RP11-214O1.2
    ## RP11-799D4.4
    ## CTBP2P7
    ## RP11-686D22.7
    ## CTD-2162K18.3
    ## RP11-47L3.1
    ## RP11-552F3.12
    ## RP11-686D22.4
    ## RP11-322E11.5
    ## RP11-577H5.5
    ## RP11-686D22.3
    ## CTC-559E9.2
    ## CTB-50L17.9
    ## CTD-2619J13.9
    ## RP11-369G6.2
    ## RP11-32B5.8
    ## RP11-316O14.1
    ## CTD-2619J13.17
    ## CTC-513N18.7
    ## ERVV-1
    ## CTB-60E11.9
    ## RP11-463C8.7
    ## RP11-497H16.9
    ## RP11-526I2.5
    ## RP11-793H13.11
    ## BNIP3P37
    ## RP11-420L9.5
    ## RP11-51F16.9
    ## MARK2P16
    ## RP4-657E11.10
    ## ITGB8-AS1
    ## CTD-2589H19.6
    ## SMIM32
    ## JARID2-DT
    ## RP4-555D20.4
    ## RP11-465B22.8
    ## RP11-506K6.4
    ## RP11-44N22.3
    ## RP11-801F7.1
    ## CTA-390C10.10
    ## RP11-508N22.12
    ## RP11-38M8.1
    ## AGBL1
    ## RP11-7F17.8
    ## CTA-315H11.2
    ## RP13-714J12.1
    ## DGKK
    ## RP11-92F20.1
    ## RP11-457D13.4
    ## U2
    ## RP11-395B7.2
    ## H4C13
    ## HNF1B
    ## CNTNAP3P2
    ## DACH1
    ## RP13-895J2.6
    ## CICP25
    ## RP11-325L12.6
    ## PPP1R26P4
    ## RP5-1056H1.2
    ## CTC-439O9.1
    ## AC008985.1
    ## CTC-490E21.11
    ## RP11-271K21.12
    ## AC083843.4
    ## CTB-96E2.6
    ## RP11-286N22.6
    ## RP11-185E8.2
    ## DNM1P41
    ## RP11-562F9.3
    ## RP1-131F15.3
    ## RP11-290O12.2
    ## RP1-72A23.4
    ## RP11-20G6.1
    ## LLNLR-307A6.1
    ## PCAT14
    ## AC009133.23
    ## LINC01338
    ## CICP14
    ## RP11-259O18.5
    ## RP3-467D16.3
    ## RP4-590F24.2
    ## LLNLF-173C4.1
    ## RP11-477N12.5
    ## C19orf85
    ## EXOC3L2
    ## RP11-511P7.5
    ## CTC-564N23.6
    ## LINC02795
    ## RP1-37E16.13
    ## RP11-457M11.8
    ## RP11-187K14.1
    ## CTD-2544D21.2
    ## TOP1P1
    ## RP1-166D18.1
    ## RP5-914M6.2
    ## RP11-401F2.5
    ## RP11-553D4.3
    ## RP11-319F12.3
    ## AC025442.4
    ## RP13-111A12.1
    ## RP11-196B3.5
    ## RP11-314O13.2
    ## RP11-660K21.1
    ## RP4-769N13.9
    ## RP11-460N20.8
    ## RP11-158H18.2
    ## RP11-249L21.6
    ## RP11-486M3.5
    ## RP11-497D24.1
    ## RP11-1145F21.2
    ## RP11-547J14.1
    ## RP11-262H14.14
    ## RP11-6G22.1
    ## RP11-379L12.1
    ## RENO1
    ## RP11-1150C11.2
    ## RP11-941F15.2
    ## RP11-72O9.5
    ## RP11-506N21.1
    ## RP1-237I15.2
    ## RP11-281P23.3
    ## RP13-726E6.5
    ## RP11-416B23.2
    ## RP11-1020A11.3
    ## RP5-1055C14.9
    ## PDCD6-AHRR
    ## TNMD
    ## CALCR
    ## SCIN
    ## MATK
    ## CAMK1G
    ## DLEC1
    ## ETV7
    ## SLC7A14
    ## CLDN11
    ## YAF2
    ## DPEP1
    ## NPC1L1
    ## IL20RA
    ## CNTN1
    ## RUNX3
    ## PLEKHB1
    ## GABRA1
    ## VIM
    ## FAS
    ## CD44
    ## AGPAT4
    ## DAPK2
    ## TLL1
    ## CDH10
    ## TNC
    ## EPHA3
    ## LMO3
    ## PTPRN
    ## HHAT
    ## SOAT1
    ## PSD
    ## GNA15
    ## CREB3L3
    ## WNT8A
    ## SEZ6
    ## DLX3
    ## ANKRD44
    ## SNAP91
    ## NAV3
    ## PITX1
    ## GAL
    ## MYO3B
    ## TRHDE
    ## SLC12A1
    ## SEMA3C
    ## FGF4
    ## PAG1
    ## MBNL3
    ## MAP2
    ## NEBL
    ## ADCY2
    ## RDH8
    ## MEF2C
    ## CACNA1S
    ## CADPS2
    ## SMARCD3
    ## P2RX5
    ## BCORL1
    ## FOLH1
    ## ACHE
    ## GP6
    ## FER1L4
    ## SLC4A11
    ## LHX5
    ## RPH3A
    ## CMTM1
    ## LAG3
    ## CHRD
    ## SMPX
    ## NLRP1
    ## ZFHX4
    ## DAZL
    ## TBX15
    ## TGFB2
    ## TBC1D2
    ## CYP26A1
    ## NRP1
    ## RASD2
    ## CHADL
    ## DHRS2
    ## PROCR
    ## SLA2
    ## HRH3
    ## MYLK2
    ## PDYN
    ## RNF125
    ## MCTS2P
    ## SYP
    ## NALCN
    ## RGCC
    ## NECAB2
    ## SLC7A5
    ## OCA2
    ## SLC30A4
    ## KCNN4
    ## TULP2
    ## FCER2
    ## RSPH6A
    ## TLE6
    ## EBI3
    ## CD33
    ## MAG
    ## RASA4
    ## CAV1
    ## VIPR2
    ## NOD1
    ## RARRES2
    ## DOCK8
    ## ELAVL2
    ## PDLIM1
    ## CUBN
    ## VSIR
    ## TLX1
    ## KAZALD1
    ## RAI1
    ## DLX4
    ## HLF
    ## CPZ
    ## GLRB
    ## CRYAB
    ## JHY
    ## FOLR1
    ## IL10RA
    ## SLC1A2
    ## SLC15A3
    ## MYF6
    ## OAS3
    ## NANOG
    ## RIPOR2
    ## RBM24
    ## KHDRBS2
    ## CRYBG1
    ## SMOC2
    ## NME5
    ## LOX
    ## CDH9
    ## CDH6
    ## STC2
    ## CNTN3
    ## BCHE
    ## FGF12
    ## OTOF
    ## SLC30A3
    ## CLIP4
    ## HTRA2
    ## GALNT3
    ## FN1
    ## DNAH6
    ## IGFBP5
    ## IL18R1
    ## MLPH
    ## KISS1R
    ## OPRD1
    ## WLS
    ## SIPA1L2
    ## CD48
    ## PLPPR4
    ## PROX1
    ## ESYT2
    ## B4GALT6
    ## CCN2
    ## ABCG2
    ## SPP1
    ## KLF9
    ## ECRG4
    ## CCDC92
    ## PROX2
    ## CNRIP1
    ## PRLHR
    ## MYCT1
    ## SLC25A2
    ## ENOX1
    ## SOHLH2
    ## EGR1
    ## ADRA1A
    ## TNFRSF8
    ## XPNPEP2
    ## HS3ST2
    ## NXPH1
    ## FKBP9
    ## EGR2
    ## BHLHE41
    ## KCNJ2
    ## MATN4
    ## SEMG1
    ## IQSEC2
    ## XG
    ## CDKN1A
    ## NRN1
    ## RAB17
    ## SH3TC1
    ## LRRC29
    ## MT1G
    ## MT2A
    ## SOX21
    ## TMEM255A
    ## BMP4
    ## MYH2
    ## HS3ST3B1
    ## OPRL1
    ## PSD4
    ## CD70
    ## FOXA2
    ## FAM182A
    ## NKX2-4
    ## FLRT3
    ## MMP24
    ## GRPR
    ## F10
    ## WNK4
    ## HIVEP3
    ## MASP1
    ## F2RL3
    ## FBXL16
    ## TNFRSF19
    ## SPINK2
    ## KDR
    ## FOXP2
    ## LRRC17
    ## OPN1SW
    ## CHAC1
    ## TPH1
    ## FGF13
    ## LBP
    ## THEM6
    ## BST2
    ## MAS1
    ## GDF15
    ## RBBP8NL
    ## NUTM2F
    ## COX4I2
    ## IDO1
    ## GALNT15
    ## TMEM204
    ## TNS4
    ## ZNF132
    ## ABHD12B
    ## PPARG
    ## HHLA1
    ## IQCA1
    ## RAMP1
    ## RIN2
    ## KANK4
    ## SYT4
    ## MTUS2
    ## DCLK1
    ## CCNA1
    ## IRS4
    ## MYO18B
    ## NTS
    ## LRRIQ1
    ## MRO
    ## MEIS2
    ## CRB1
    ## SOX3
    ## PHC2
    ## DSC3
    ## APLNR
    ## ANXA1
    ## FAM189A2
    ## AGAP2
    ## B4GALNT1
    ## SLC26A10
    ## NIBAN1
    ## CHRND
    ## WNT10A
    ## ALDH1L2
    ## THSD1
    ## EDNRB
    ## DGKB
    ## GALNT5
    ## BLK
    ## SKIL
    ## GABBR2
    ## MYC
    ## IER3
    ## SLCO2B1
    ## TRIM29
    ## TMPRSS13
    ## CYP19A1
    ## SEMA6D
    ## LOXL4
    ## EXOC6
    ## STAT4
    ## CDK15
    ## SLC40A1
    ## ARHGAP24
    ## PCDH10
    ## INHBE
    ## TCHP
    ## ACVRL1
    ## LPAR6
    ## SSTR1
    ## KCNH5
    ## SLC24A4
    ## CYP1A1
    ## CCDC33
    ## ITGAX
    ## NLRC5
    ## CDH11
    ## ADAD2
    ## DPEP3
    ## GREB1L
    ## SLC13A5
    ## PNMT
    ## NLRP12
    ## PRDM16
    ## SMYD2
    ## LYST
    ## SYT2
    ## FBLN7
    ## FAM171B
    ## MARCHF4
    ## GADL1
    ## STAC
    ## RP11-977G19.10
    ## IGSF11
    ## ECE2
    ## SLC10A4
    ## ADAMTS16
    ## UGT3A1
    ## CXCL14
    ## GFRA3
    ## DCDC2
    ## RASGEF1C
    ## HMGCLL1
    ## SLC22A3
    ## CCNB3
    ## HTR2C
    ## MAGEA4
    ## DOK2
    ## SNTG1
    ## TERF1
    ## RSPO2
    ## VLDLR
    ## STOM
    ## COMMD3
    ## ZEB1
    ## NPFFR1
    ## ADAM12
    ## SLC43A1
    ## NCAM1
    ## KIRREL3
    ## TMEM25
    ## KIAA1755
    ## CACNA2D4
    ## AKAP6
    ## AKR1C2
    ## VENTX
    ## TMEM45B
    ## FBXO4
    ## OBI1
    ## GRIA4
    ## SPEF2
    ## MR1
    ## ACOXL
    ## BMP6
    ## NR4A2
    ## SCN3A
    ## GRAP
    ## AK5
    ## ROBO4
    ## NRGN
    ## ANGPT1
    ## C4orf19
    ## FAM167A
    ## CNTNAP3B
    ## SORBS2
    ## NCAM2
    ## PDE1C
    ## SLFN13
    ## APCDD1
    ## CA10
    ## SAMSN1
    ## TTN
    ## SPAG17
    ## ADCY8
    ## CFAP70
    ## GNA14
    ## WIF1
    ## SORCS3
    ## FUT6
    ## PPP2R2B
    ## LRFN2
    ## NRG1
    ## KCNJ6
    ## SLC34A2
    ## GRHL3
    ## AUTS2
    ## NRG2
    ## HTR6
    ## NPM2
    ## DMTN
    ## CLIC6
    ## CBR3
    ## GJD2
    ## FCN2
    ## IKZF3
    ## KASH5
    ## PRR35
    ## FGF19
    ## BSND
    ## IL23R
    ## ADGRL4
    ## ATF3
    ## NBPF20
    ## CAPN13
    ## KCNJ3
    ## ALMS1P1
    ## VSNL1
    ## BBS5
    ## ARHGAP25
    ## PROK2
    ## CCDC141
    ## EOMES
    ## DPPA2
    ## SPTA1
    ## ADAMTS9
    ## HESX1
    ## LIPH
    ## HPGD
    ## SERINC5
    ## ACSL6
    ## SAMD3
    ## KIAA0895
    ## GPR85
    ## PGAM2
    ## MAMDC2
    ## C8orf34
    ## ALDH1A1
    ## TRPV6
    ## FBP1
    ## RNF183
    ## LRFN5
    ## RPL10L
    ## NGB
    ## STOX1
    ## IFI27
    ## GPR176
    ## FBN1
    ## CRABP1
    ## PRTG
    ## PKD1L2
    ## A2ML1
    ## GALR1
    ## MESP1
    ## ANPEP
    ## IGF2
    ## CA4
    ## NIBAN3
    ## GPD1
    ## KLK6
    ## IGFBP6
    ## TMEM88
    ## PLAAT5
    ## SCARA5
    ## GNG4
    ## FAM107A
    ## BMP1
    ## GBX2
    ## ADAM9
    ## LDLRAD4
    ## KCTD19
    ## CXXC4
    ## NSG1
    ## FSTL5
    ## ATOH8
    ## AFAP1L2
    ## SH3TC2
    ## SHE
    ## ZNF280A
    ## ZEB2
    ## BNC1
    ## MT1E
    ## UGP2
    ## ONECUT1
    ## PYDC1
    ## FRMPD4
    ## FAM153A
    ## SERPINA6
    ## DCLK2
    ## SDR16C5
    ## TEX13B
    ## PLAC1
    ## KCNG3
    ## PRKCE
    ## NPTX1
    ## LPAR3
    ## DSCAM
    ## NHLH1
    ## UTF1
    ## PCDHB1
    ## FBXL14
    ## MIR31HG
    ## FOXB1
    ## DRC3
    ## GAP43
    ## SLFN12
    ## PDE3A
    ## SLFN11
    ## CYP7B1
    ## DELEC1
    ## ABLIM3
    ## KCNK7
    ## TAS1R1
    ## RP11-97O12.7
    ## PIFO
    ## MSRB3
    ## ZCCHC12
    ## LEP
    ## RTP1
    ## DES
    ## VWA3A
    ## PCSK1
    ## LPL
    ## LONRF3
    ## FOSL1
    ## CH17-12M21.1
    ## SYNE3
    ## PLAAT3
    ## LRRC37A3
    ## GCNT4
    ## KCNA2
    ## IRX3
    ## GRAMD1C
    ## RNF212
    ## TMEM151B
    ## SLITRK1
    ## DNAJC22
    ## CPNE7
    ## EXOC3L1
    ## HTR1D
    ## CDH5
    ## LRRC3B
    ## NAIPP2
    ## TSHZ1
    ## C3orf80
    ## OR3A1
    ## RRH
    ## SHISA2
    ## ZNF467
    ## OR52A1
    ## ASB18
    ## NLRP10
    ## IGIP
    ## SLC35D3
    ## SOX1
    ## PCP4
    ## CSMD1
    ## ANKRD20A2P
    ## TMEM119
    ## CTNNA3
    ## PLGLB1
    ## MARCHF11
    ## LINC00518
    ## OPCML
    ## ACP7
    ## ZNF703
    ## FAM162B
    ## CCR4
    ## LRRC55
    ## NPW
    ## OAF
    ## TACSTD2
    ## GDF3
    ## SLIT3
    ## DHRS7C
    ## PDE4B
    ## KCNH7
    ## RALYL
    ## TRARG1
    ## SRY
    ## MANEAL
    ## PURA
    ## LINC00313
    ## KCNQ5
    ## IFITM1
    ## RTN4RL1
    ## TEX38
    ## FAM86MP
    ## GABRA5
    ## TRDN
    ## HPDL
    ## TRABD2A
    ## FOXD3
    ## MT1X
    ## MAGI2
    ## DPPA3
    ## PEAR1
    ## C2orf78
    ## DNER
    ## SELL
    ## CERKL
    ## SERPINA5
    ## FBLL1
    ## VWC2
    ## OPTC
    ## HMX2
    ## ASTL
    ## AADACL3
    ## ANKRD34B
    ## PCDH18
    ## DPY19L2P1
    ## OR7E12P
    ## PAX5
    ## SPOCK3
    ## AKR1C3
    ## SPATS2L
    ## KIF19
    ## STK31
    ## ZNF781
    ## NOXO1
    ## CRYBA4
    ## FGF16
    ## SPTSSB
    ## AJAP1
    ## HRH1
    ## CRACDL
    ## TMEM26
    ## SCOC-AS1
    ## ADAM32
    ## SULT1A2
    ## PIWIL2
    ## UAP1L1
    ## ZNF860
    ## CYP2B6
    ## AC005276.1
    ## SSPOP
    ## MYH6
    ## DNAH10
    ## CR1L
    ## CTC-301O7.4
    ## MT1F
    ## RTP2
    ## C2CD4A
    ## RD3
    ## RYR2
    ## DLL1
    ## ANKRD13B
    ## COLGALT2
    ## GRM3
    ## GJC2
    ## CES1
    ## CCDC167
    ## TGM2
    ## RORB
    ## SNORA71C
    ## C1orf53
    ## SAMD5
    ## RIPPLY2
    ## DPPA5
    ## MAFB
    ## SP5
    ## C9orf129
    ## PRAMEF7
    ## ACOXL-AS1
    ## FOXB2
    ## MOG
    ## C9orf135
    ## KRT17P5
    ## MT1H
    ## CCDC85C
    ## RP9P
    ## C21orf62
    ## CFAP44
    ## RP11-316F12.1
    ## COLQ
    ## XKR4
    ## TRBC2
    ## YWHAZP4
    ## KRT222
    ## GPSM3
    ## PRDX1P1
    ## NOTO
    ## NPIPA7
    ## TCEA1P3
    ## SCRT2
    ## ZNF663P
    ## LINC01356
    ## ATAD3C
    ## RPL22P12
    ## DNAJA1P4
    ## CEBPZOS
    ## CELA3B
    ## RP11-3B12.3
    ## LINGO3
    ## HMSD
    ## RP11-545E17.3
    ## AC027612.4
    ## TSSC2
    ## HRAT92
    ## MIR181A2HG
    ## CDYLP1
    ## DNAJC27-AS1
    ## HSPA8P7
    ## RP11-82L18.2
    ## SLC9A3-AS1
    ## RP11-415J8.3
    ## AC104134.2
    ## LINC01107
    ## FGF13-AS1
    ## NOP56P2
    ## LINC01108
    ## RPL15P14
    ## RP11-432J24.5
    ## AC096574.5
    ## AC092535.3
    ## LINC01806
    ## NCAM1-AS1
    ## RPS3AP37
    ## INKA2-AS1
    ## RP11-381K7.1
    ## FAM66A
    ## OR2M3
    ## MIR34AHG
    ## RP11-91I20.3
    ## VIM-AS1
    ## LINC01981
    ## CICP3
    ## MAST4-AS1
    ## AC018462.2
    ## ZRANB2-AS2
    ## LINC00102
    ## LINC00102
    ## AC092162.1
    ## AC159540.1
    ## AC112229.1
    ## AC092171.4
    ## ANKRD63
    ## FOXD3-AS1
    ## MIR205HG
    ## RP11-38L15.3
    ## LINC02884
    ## PCED1CP
    ## CFAP97D1
    ## BEND3P1
    ## HLA-DPA1
    ## RP1-45C12.1
    ## OR55B1P
    ## DIRC3
    ## TSPAN19
    ## RFX3-AS1
    ## IMPDH1P10
    ## AC073257.1
    ## LINC00865
    ## AC002480.3
    ## SNX18P2
    ## DIRC3-AS1
    ## AP000688.29
    ## LINC01783
    ## CICP27
    ## RP11-115A15.2
    ## RHBDF1P1
    ## RP11-146I2.1
    ## LINC01426
    ## RP11-561O23.5
    ## LINC00458
    ## DNM1P51
    ## SSR4P1
    ## TRDN-AS1
    ## AF127936.9
    ## AC009299.2
    ## PHACTR2-AS1
    ## AC000089.3
    ## CPB2-AS1
    ## APOA1-AS
    ## H2BW3P
    ## TRHDE-AS1
    ## VLDLR-AS1
    ## RP11-69I8.2
    ## RP11-73E6.2
    ## TMEM26-AS1
    ## RP4-680D5.2
    ## RPL21P121
    ## SHISA9
    ## RP11-383J24.2
    ## OXCT2P1
    ## CT75
    ## AC079922.3
    ## LINC01370
    ## RP11-216M21.1
    ## RP6-24A23.3
    ## PSME2P1
    ## LINC01410
    ## RPL30P12
    ## ARHGEF25
    ## AF011889.5
    ## LINC02029
    ## XACT
    ## RP11-435F17.3
    ## HLA-DMB
    ## LINC00882
    ## AC005062.2
    ## MTND5P16
    ## RN7SL189P
    ## DUXAP10
    ## OR2A1-AS1
    ## MALINC1
    ## RP11-334E6.10
    ## LINC00461
    ## CRNDE
    ## RP11-629G13.1
    ## RP11-446J8.1
    ## LINC02275
    ## R3HDM2P1
    ## C8orf34-AS1
    ## RP11-419C19.1
    ## LNCPRESS2
    ## SEMA6A-AS2
    ## PGAM1P9
    ## CTD-2593A12.3
    ## APOBEC3B-AS1
    ## MTCO1P24
    ## CARMN
    ## TMEM158
    ## SLC7A11-AS1
    ## DNAH10OS
    ## CTB-114C7.4
    ## TRMT9B
    ## TUNAR
    ## RP11-631M6.2
    ## CTD-2297D10.2
    ## LINC02381
    ## TMEM30BP1
    ## RP11-115D19.1
    ## CTBP2P4
    ## DDX18P4
    ## NAIPP4
    ## OC90
    ## RP11-659E9.4
    ## CTD-2501M5.1
    ## SLC10A5
    ## ETV3L
    ## MAP2K1P1
    ## RPL7P20
    ## RP11-267L5.1
    ## LRRC24
    ## AP000662.4
    ## CKLF-CMTM1
    ## FAM66D
    ## RP11-567M21.3
    ## CTC-342M10.2
    ## TBX5-AS1
    ## LINC02700
    ## SALL3
    ## AP000438.2
    ## LINC02293
    ## LINC01490
    ## RP11-887P2.5
    ## RP11-240B13.2
    ## RP11-255M2.1
    ## LINC-ROR
    ## RP11-688G15.3
    ## RP11-718G2.5
    ## NDUFC2-KCTD14
    ## RP11-82L7.4
    ## CTD-2054N24.2
    ## TPM1-AS
    ## RBM17P4
    ## DNM1P47
    ## CSPG4P11
    ## RP11-261B23.1
    ## RP13-395E19.2
    ## FRRS1L
    ## RP11-297M9.1
    ## RP11-375I20.6
    ## RP11-184E9.2
    ## RP5-882C2.2
    ## SERTM2
    ## RP11-256I9.2
    ## FBXL19-AS1
    ## RP11-274H2.5
    ## RP11-6O2.4
    ## LINC02473
    ## RBFADN
    ## FIGNL2
    ## CH17-260O16.1
    ## TMEM249
    ## SCRT1
    ## CTC-523E23.1
    ## TCF24
    ## DKFZP434A062
    ## SPDYE9
    ## U91319.1
    ## DLGAP1-AS3
    ## TWSG1-DT
    ## ZNF516-AS1
    ## BAHCC1
    ## LINC00683
    ## RP11-214O1.2
    ## RP11-799D4.4
    ## CTBP2P7
    ## RP11-686D22.7
    ## RP11-795F19.5
    ## RP11-47L3.1
    ## RP11-13K12.5
    ## RP11-686D22.4
    ## STK25P1
    ## RP11-686D22.3
    ## CTC-559E9.2
    ## LINC01801
    ## RP11-171I2.1
    ## RP11-369G6.2
    ## RP11-316O14.1
    ## BNIP3P26
    ## DMRTC1
    ## NBPF9
    ## AF003625.3
    ## SCARNA2
    ## RP11-362K14.6
    ## RP11-526I2.5
    ## RP11-793H13.11
    ## NBPF8
    ## BNIP3P37
    ## RP11-420L9.5
    ## NAV2-AS6
    ## RP4-657E11.10
    ## AC015849.19
    ## RAP2CP1
    ## LRRC37A9P
    ## ITGB8-AS1
    ## CTD-2589H19.6
    ## RP4-736L20.3
    ## RP4-555D20.4
    ## RP3-325F22.5
    ## DOC2B
    ## RP11-66B24.7
    ## RP11-484K9.4
    ## CTA-390C10.10
    ## RP11-38M8.1
    ## RP11-248G5.9
    ## AGBL1
    ## LINC02334
    ## RP11-92F20.1
    ## U2
    ## RP11-15E1.5
    ## H4C13
    ## RP13-753N3.3
    ## HNF1B
    ## CH507-24F1.2
    ## ZNF516-DT
    ## RIMBP3
    ## CNTNAP3P2
    ## DACH1
    ## RP13-895J2.6
    ## RP11-622C24.2
    ## RP4-794I6.4
    ## RP11-325L12.6
    ## RP11-129M16.4
    ## RP11-679B19.1
    ## PPP1R26P4
    ## RP11-4B16.4
    ## CTC-439O9.1
    ## RP11-399B17.1
    ## CTC-490E21.11
    ## AC083843.4
    ## RP11-133K1.8
    ## RP11-4L24.4
    ## RP11-185E8.2
    ## DNM1P41
    ## RP1-131F15.3
    ## RP11-463O12.5
    ## RP11-20G6.1
    ## LLNLR-307A6.1
    ## RP11-35J1.2
    ## AC009133.23
    ## CICP14
    ## RP11-259O18.5
    ## RP3-467D16.3
    ## RP4-590F24.2
    ## RP11-477N12.5
    ## LSP1P5
    ## CTC-591M7.1
    ## C19orf85
    ## RP11-511P7.5
    ## CTC-564N23.6
    ## LINC02795
    ## RP11-457M11.8
    ## RP11-187K14.1
    ## TOP1P1
    ## RP11-401F2.5
    ## RP11-553D4.3
    ## AC025442.4
    ## FAS-AS1
    ## RP11-400O10.1
    ## RP11-314O13.2
    ## RP4-769N13.9
    ## RP11-249L21.6
    ## RP11-1150C11.1
    ## CTB-88F18.6
    ## RP11-399E6.5
    ## CTD-2335I23.3
    ## CTC-756D1.4
    ## RP11-262H14.14
    ## CH17-353B19.2
    ## RP11-6G22.1
    ## RENO1
    ## RP11-1150C11.2
    ## RP5-1055C4.1
    ## RP11-941F15.2
    ## RP11-791G22.3
    ## RP11-72O9.5
    ## RP11-506N21.1
    ## AC061961.2
    ## RP13-726E6.5
    ## RP11-416B23.2
    ## TMX2-CTNND1
    ## RP11-1020A11.3
    ## RP5-1055C14.9
    ## PDCD6-AHRR
    ## RP11-177J6.2

``` r
# as data frame:
human_genes <- as.data.frame(filtered_res_df_2$gene_name, collapse = "\n")
# Let's make sure to only get the unique names
human_genes <- unique(human_genes)

# Now let's write this out and do gene enrichment analysis
write.table(human_genes["filtered_res_df_2$gene_name"], row.names = FALSE, col.names = FALSE, "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/human_gene_names.csv")

save(filtered_res_df_2, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/HUMAN/filtered_res_df_2")
# copy and paste into EnrichR 
#https://maayanlab.cloud/Enrichr/
```

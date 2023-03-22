Class\_exercise
================
Megan Hupka
3/22/2023

# Load the libraries you need

# Load functions you need “my\_class\_functions”

# load in your peak files for each replicate of each protein

# Here I am starting to analyze my data for my proteins of interest:

# proteinX, Y, Z …..

# First I will read in each replicate file

``` r
# Loading in the peak files into the variable peak_list
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# Creating a talbe of the number of peaks in each file
peak_num <- sapply(peak_list, length) %>% as.data.frame(row.names = T)
```

    ## Warning in as.data.frame.integer(., row.names = T): 'row.names' is not a
    ## character vector of length 18 -- omitting it. Will be an error!

``` r
# label column
names(peak_num) <- c("num_peaks")

# make dbp name a col.
peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")
  # peak_num <- separate(peak_num, col = dbp,  into = c('dbp', 'replicate'), sep = "_")

# let's save this for our records 
write_csv(peak_num, "results/num_peaks_df.csv")


# printing out a table of the number of peaks in each file:
peak_num
```

    ##        dbp replicate num_peaks
    ## 1   ARID3A        R1     47962
    ## 2   ARID3A        R2     30842
    ## 3     ATF3        R1      3848
    ## 4     ATF3        R2     35018
    ## 5     ATF3        R3     49209
    ## 6     ATF3        R4     60518
    ## 7  BHLHE40        R1     18183
    ## 8  BHLHE40        R2      2618
    ## 9  BHLHE40        R3     10119
    ## 10 BHLHE40        R4      8074
    ## 11   BRCA1        R1      2212
    ## 12   BRCA1        R2      4073
    ## 13   BRCA1        R3     44978
    ## 14   BRCA1        R4     44173
    ## 15   CEBPB        R1     47888
    ## 16   CEBPB        R2     26625
    ## 17   CEBPB        R3     53316
    ## 18   CEBPB        R4    187109

# Now I am going to create consensus peaks for each protein

``` r
# List of unique dbps
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

# running the function consensus_from_reduced
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

# export consensus peaks to results folder

# setting file path to export
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/mehu6123"
consensus_path <- "CLASS_2023/CLASSES/class_exercise/results/"
exportpath <- file.path(basepath, consensus_path)

# exporting as .bed files
for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus_peaks.bed") )}
```

# Now I am going to make my consensus peaks compatable with UCSC genome browser

``` r
consensusFilePath <- "/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results"
exportFilePath <- "/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/UCSC/"

ucsc_formating(consensusFilePath = consensusFilePath, export_path = exportFilePath)
```

    ## [1] 10

    ## [1] "done?"

``` r
# print out consensus peak files in a results/UCSC directory       (??????)
```

# I am curious if my proteins are transcription factors so I will use the annotations

# in a cell paper I found and see

``` r
# creating a data frame with our consensus peaks to use
num_peaks_df <- data.frame("dbp" = names(consensus_list),
                           "num_peaks" = sapply(consensus_list, length))


# Downloading the information on transcription factors
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"
destination_for_url <- "results/TF_annotations.xlsx"
# to download we can use download.file
download.file(url, destination_for_url)

#redx1::read_excel to import
human_tfs <- readxl::read_excel("results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)
```

    ## Warning: Expecting logical in M1006 / R1006C13: got 'Contains a SANT and
    ## multiple DNA-binding C2H2 domains. Motif is 99% AA ID from mouse (Transfac).'

    ## Warning: Expecting logical in M1021 / R1021C13: got 'Close ortholog (PP1RA)
    ## binds to mRNA; single-stranded DNA (ssDNA); poly(A) and poly(G) homopolymers
    ## (Uniprot)'

    ## Warning: Expecting logical in M1542 / R1542C13: got 'Contains 1 SANT domain'

    ## Warning: Expecting logical in M1543 / R1543C13: got 'Contains 2 Myb DBDs.
    ## Sources of Hocomoco/Transfac motifs are unclear. However these sequences look
    ## similar to in vitro sites selected by SELEX (PMID:11082045)'

    ## Warning: Expecting logical in M1544 / R1544C13: got 'Although CHD2 has weak
    ## similarity to a Myb domain (PMID:9326634), it's more closely related to the
    ## non-DNA-binding SANT domain based on our alignment analysis. The data showing
    ## that show that CHD2 binding histone H3.3 (PMID:22569126) further support the
    ## conclusion that the Myb domain is probably a SANT domain facilitating the
    ## histone interaction'

    ## Warning: Expecting logical in M1545 / R1545C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1546 / R1546C13: got 'Contains 2 Myb DBDs'

    ## Warning: Expecting logical in M1547 / R1547C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1548 / R1548C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1549 / R1549C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1550 / R1550C13: got 'Domain is truncated, and
    ## there is nothing known about this gene'

    ## Warning: Expecting logical in M1551 / R1551C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1552 / R1552C13: got 'MIER2's Myb domain is more
    ## similar to the non-DNA-binding SANT domain'

    ## Warning: Expecting logical in M1553 / R1553C13: got 'MIER3's Myb domain is more
    ## similar to the non-DNA-binding SANT domain'

    ## Warning: Expecting logical in M1554 / R1554C13: got 'Contains 1 SANT domain,
    ## and a SANTA domain'

    ## Warning: Expecting logical in M1555 / R1555C13: got 'Contains a single Myb-like
    ## domain with an insertion in the middle. It is ambiguous whether Myb-like
    ## domains are DNA or protein binding. Since it has a single domain it's likely
    ## non-specific, but future experiments should be performed to assay it's
    ## specificity'

    ## Warning: Expecting logical in M1556 / R1556C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1557 / R1557C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1558 / R1558C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1559 / R1559C13: got 'Contains a single Myb-like
    ## domain. Mouse ortholog has motif'

    ## Warning: Expecting logical in M1560 / R1560C13: got 'MYSM1 has been shown to
    ## bind DNA ? interaction with DNA requires the MYSM1 Myb but not the SWIRM domain
    ## (PMID:17428495). Domain sequence alignment places it near DNA-binding Myb
    ## domains but scores slightly higher as a SANT rather than Myb domain based on
    ## Prosite patterns. Given that most Myb proteins that bind DNA sequence
    ## specifically have multiple Myb domains in an array this protein could bind DNA
    ## sequence non-specifically with it?s single Myb domain. Future experiments
    ## should assay MYSM1?s specificity'

    ## Warning: Expecting logical in M1561 / R1561C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1562 / R1562C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1564 / R1564C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1565 / R1565C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1566 / R1566C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains. RCOR3 SANT domains are known to
    ## facilitate PPIs'

    ## Warning: Expecting logical in M1567 / R1567C13: got 'SMARCA1 contains a
    ## truncated Myb-like and SANT domain. Given the presence of the Myb-like domain,
    ## and other domains known to associated with DNA (DEAD box helicase) it likely
    ## associates with DNA non-sequence-specifically'

    ## Warning: Expecting logical in M1568 / R1568C13: got 'Contains a SANT, and
    ## Myb-like domain'

    ## Warning: Expecting logical in M1569 / R1569C13: got 'Contains 1 SANT domain,
    ## and no other putative DNA-binding domains. Motif logos look like bZIP dimeric
    ## binding sites, and are thus likely specificifities of SMARCC1 interactors'

    ## Warning: Expecting logical in M1570 / R1570C13: got 'Contains 1 SANT domain,
    ## and no other putative DNA-binding domains. Motif logos ares likely
    ## specificifities of SMARCC2 interactors'

    ## Warning: Expecting logical in M1571 / R1571C13: got 'Contains only Myb DBDs'

    ## Warning: Expecting logical in M1572 / R1572C13: got 'Contains 1 SANT domain'

    ## Warning: Expecting logical in M1573 / R1573C13: got 'TADA2B contains a single
    ## SANT domain and is thus unlikely to bind DNA'

    ## Warning: Expecting logical in M1574 / R1574C13: got 'Contains a single Myb
    ## domain (with slightly less simialrity to a SANT domain.) This domain has been
    ## shown to be involved in PPIs but this may not be mutually exclusive with
    ## DNA-binding. The sequence-specificity of CCDC79 should be investigated in the
    ## future'

    ## Warning: Expecting logical in M1575 / R1575C13: got 'Contains 1 Myb domain, and
    ## has structural evidence of DNA-binding'

    ## Warning: Expecting logical in M1576 / R1576C13: got 'Motif is inferred from
    ## mouse (92% DBD AA ID)'

    ## Warning: Expecting logical in M1577 / R1577C13: got 'TERF2IP contains a single
    ## Myb-like domain. While it's unclear if TERF2IP (Human Rap1) contacts DNA
    ## directly it has been shown to affect the DNA binding activity of TRF2'

    ## Warning: Expecting logical in M1578 / R1578C13: got 'This protein contains Myb,
    ## and Myb-like domains and is annotated as a Pol1 terminator. TTF1 DNA-binding
    ## has been demonstrated in vitro (PMID: 7597036), but it's specificity has not
    ## been determined'

    ## Warning: Expecting logical in M1579 / R1579C13: got 'Contains 1 Myb DBD'

    ## Warning: Expecting logical in M1580 / R1580C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). The GATA ZnF domain of MTA1 is required for it's
    ## interaction with RBBP4 and RBBP7 (PMID:18067919). Full-length protein has been
    ## tried in HT-SELEX and did not yield a motif'

    ## Warning: Expecting logical in M1581 / R1581C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Full-length protein has been tried in HT-SELEX,
    ## and DBD has been tried on PBM - neither yielded motifs'

    ## Warning: Expecting logical in M1582 / R1582C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Hasn't been tried in any in vitro assays'

    ## Warning: Expecting logical in M1583 / R1583C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Has been tried as a DBD in HT-SELEX but did not
    ## yield a motif'

    ## Warning: Expecting logical in M1791 / R1791C13: got 'CNOT3 is a part of the
    ## CCR4-NOT complex involved in mRNA decay'

    ## Warning: Expecting logical in M1932 / R1932C13: got '"Prosite identifies a
    ## low-confidence Myb-like domain (e.g. can?t decide between Myb and SANT) so it?s
    ## probably not a TF"'

    ## New names:
    ## • `` -> `...4`

``` r
# let's rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"
# now let's intersect gene names that are in our ChIP data and has TF identity.
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 5

``` r
#
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]
# adding new column names
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

# Merging
num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)

# Let's check how many NAs -- we should have some missing values.
dim(num_peaks_df[is.na(num_peaks_df$tf),])
```

    ## [1] 0 5

``` r
# Adding a few more features to num_peaks
num_peaks_df$total_peak_length <- sapply(consensus_list, function(x) sum(width(x)))

# Saving the results
write_csv(num_peaks_df, "results/num_peaks_df.csv")

num_peaks_df
```

    ##       dbp num_peaks      ensembl_id         dbd  tf total_peak_length
    ## 1  ARID3A     23669 ENSG00000116017 ARID/BRIGHT Yes           9316873
    ## 2    ATF3      2876 ENSG00000162772        bZIP Yes           3049787
    ## 3 BHLHE40      1234 ENSG00000134107        bHLH Yes            481839
    ## 4   BRCA1      1160 ENSG00000012048     Unknown  No           1496055
    ## 5   CEBPB     12519 ENSG00000172216        bZIP Yes           4109415

# Now I want to compare a protein with a previous analysis

``` r
# goto UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses


knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/UCSC/ATF3.jpg")
```

<img src="results/UCSC/ATF3.jpg" width="1894" />

``` r
knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/UCSC/CEBPB.jpg")
```

<img src="results/UCSC/CEBPB.jpg" width="1902" />

``` r
knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/UCSC/BRCA1.jpg")
```

<img src="results/UCSC/BRCA1.jpg" width="1888" />

``` r
knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/UCSC/BHLHE40.jpg")
```

<img src="results/UCSC/BHLHE40.jpg" width="1896" />

# Now I am going to determine how my peaks for each protein overlap annotations of the genome

# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
lncrna_mrna_promoters <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/05_R_analyses/01_peak_features/results/gene_annotations/lncrna_mrna_promoters.gtf")

promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, type = "counts")
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

write_csv(num_peaks_df, "results/num_peaks_df.csv")

num_peaks_df
```

    ##       dbp num_peaks      ensembl_id         dbd  tf total_peak_length
    ## 1  ARID3A     23669 ENSG00000116017 ARID/BRIGHT Yes           9316873
    ## 2    ATF3      2876 ENSG00000162772        bZIP Yes           3049787
    ## 3 BHLHE40      1234 ENSG00000134107        bHLH Yes            481839
    ## 4   BRCA1      1160 ENSG00000012048     Unknown  No           1496055
    ## 5   CEBPB     12519 ENSG00000172216        bZIP Yes           4109415
    ##   peaks_overlapping_promoters
    ## 1                        4601
    ## 2                        2659
    ## 3                         628
    ## 4                        1481
    ## 5                        1742

## results:

\#1) What can you determine from these overlaps?

# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
#lncRNA
lncrna_gene_ids <- lncrna_mrna_promoters$gene_id[lncrna_mrna_promoters$gene_type == "lncRNA"]

num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mRNA
mrna_gene_ids <- lncrna_mrna_promoters$gene_id[lncrna_mrna_promoters$gene_type == "protein_coding"]

num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

num_peaks_df
```

    ##       dbp num_peaks      ensembl_id         dbd  tf total_peak_length
    ## 1  ARID3A     23669 ENSG00000116017 ARID/BRIGHT Yes           9316873
    ## 2    ATF3      2876 ENSG00000162772        bZIP Yes           3049787
    ## 3 BHLHE40      1234 ENSG00000134107        bHLH Yes            481839
    ## 4   BRCA1      1160 ENSG00000012048     Unknown  No           1496055
    ## 5   CEBPB     12519 ENSG00000172216        bZIP Yes           4109415
    ##   peaks_overlapping_promoters peaks_overlapping_lncrna_promoters
    ## 1                        4601                               1259
    ## 2                        2659                                532
    ## 3                         628                                159
    ## 4                        1481                                286
    ## 5                        1742                                529
    ##   peaks_overlapping_mrna_promoters
    ## 1                             3342
    ## 2                             2127
    ## 3                              469
    ## 4                             1195
    ## 5                             1213

## results:

# 1) What is the difference in overlaps between mRNA and lncRNA promoters

Most of the peaks seem to overlap with mRNA promoters instead of lncRNA
promoters.

# Now I am going to test if there is more binding over gene bodies than promoters

# I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
genebody_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, 
                                                consensus_list, 
                                                type = "counts")

num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)


# lncRNA
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA
num_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])

num_peaks_df
```

    ##       dbp num_peaks      ensembl_id         dbd  tf total_peak_length
    ## 1  ARID3A     23669 ENSG00000116017 ARID/BRIGHT Yes           9316873
    ## 2    ATF3      2876 ENSG00000162772        bZIP Yes           3049787
    ## 3 BHLHE40      1234 ENSG00000134107        bHLH Yes            481839
    ## 4   BRCA1      1160 ENSG00000012048     Unknown  No           1496055
    ## 5   CEBPB     12519 ENSG00000172216        bZIP Yes           4109415
    ##   peaks_overlapping_promoters peaks_overlapping_lncrna_promoters
    ## 1                        4601                               1259
    ## 2                        2659                                532
    ## 3                         628                                159
    ## 4                        1481                                286
    ## 5                        1742                                529
    ##   peaks_overlapping_mrna_promoters peaks_overlapping_genebody
    ## 1                             3342                       4601
    ## 2                             2127                       2659
    ## 3                              469                        628
    ## 4                             1195                       1481
    ## 5                             1213                       1742
    ##   peaks_overlapping_lncrna_genebody peaks_overlapping_mrna_genebody
    ## 1                              1259                            3342
    ## 2                               532                            2127
    ## 3                               159                             469
    ## 4                               286                            1195
    ## 5                               529                            1213

## results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

# It is nice and all to find overlaps, but I am interested in how many proteins

# bind a specific promoter. I will use my handy “occurence” parameter in

# " count peaks per feature"

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, 
                                               type = "occurrence")

# Let's double check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

# Great we will use this quite a bit moving forward so let's write it out! 
write.table(promoter_peak_occurence, "results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Creating a dataframe
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))

write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")
```

## results: I find the max number of proteins on a promoter to be X

# Now I want to start plotting my results

# First I will see if there is a realtionship between peak number and total DNA covered

``` r
num_peaks_df <- read_csv('/scratch/Shares/rinnclass/CLASS_2023/mehu6123/CLASS_2023/CLASSES/class_exercise/results/num_peaks_df.csv')
```

    ## Rows: 5 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): dbp, ensembl_id, dbd, tf
    ## dbl (3): num_peaks, total_peak_length, peaks_overlapping_promoters
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
  geom_point() + 

  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Now I want to color my plot by wether the protein is a TF or not.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                 color = tf == "Yes")) +
  geom_point() +
  
  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggplot
```

    ## function (data = NULL, mapping = aes(), ..., environment = parent.frame()) 
    ## {
    ##     UseMethod("ggplot")
    ## }
    ## <bytecode: 0x3ea9580>
    ## <environment: namespace:ggplot2>

# I want to make a histogram of the number of peaks for each of my proteins

``` r
ggplot(num_peaks_df, aes(x = num_peaks)) +
  geom_histogram(bins = 30) + 
  
  ylab("Number of Proteins") +
  xlab("Number of peaks") +
  ggtitle("Number of Peaks for Each Protein")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
hist
```

    ## function (x, ...) 
    ## UseMethod("hist")
    ## <bytecode: 0x10a7ab8>
    ## <environment: namespace:graphics>

# Now I want to facet this by the type of DNA binding domain my protein has.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, fill = dbd)) +
  geom_histogram(bins = 30) + 
  
    
  ylab("Number of Proteins") +
  xlab("Number of peaks") +
  ggtitle("Number of Peaks for Each Protein")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Cool now I am ready to send my result to my collaborator as a

# Knitted document

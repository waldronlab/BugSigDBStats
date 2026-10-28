# 

## Setup

``` r

library(bugsigdbr)
library(BugSigDBStats)
library(ComplexHeatmap)
library(ggpubr)
```

## Reading data

Get bulk export from bugsigdb.org:

``` r

full.dat <- bugsigdbr::importBugSigDB(version = "devel", cache = FALSE)
dim(full.dat)
```

    ## [1] 16226    53

``` r

colnames(full.dat)
```

    ##  [1] "BSDB ID"                    "Study"                     
    ##  [3] "Study design"               "PMID"                      
    ##  [5] "DOI"                        "URL"                       
    ##  [7] "Authors list"               "Title"                     
    ##  [9] "Journal"                    "Year"                      
    ## [11] "Keywords"                   "Review state"              
    ## [13] "Quality control tag"        "Experiment"                
    ## [15] "Location of subjects"       "Host species"              
    ## [17] "Body site"                  "UBERON ID"                 
    ## [19] "Condition"                  "EFO ID"                    
    ## [21] "Group 0 name"               "Group 1 name"              
    ## [23] "Group 1 definition"         "Group 0 sample size"       
    ## [25] "Group 1 sample size"        "Antibiotics exclusion"     
    ## [27] "Sequencing type"            "16S variable region"       
    ## [29] "Sequencing platform"        "Data transformation"       
    ## [31] "Statistical test"           "Significance threshold"    
    ## [33] "MHT correction"             "LDA Score above"           
    ## [35] "Matched on"                 "Confounders controlled for"
    ## [37] "Pielou"                     "Shannon"                   
    ## [39] "Chao1"                      "Simpson"                   
    ## [41] "Inverse Simpson"            "Richness"                  
    ## [43] "Signature page name"        "Source"                    
    ## [45] "Curated date"               "Curator"                   
    ## [47] "Revision editor"            "Description"               
    ## [49] "Abundance in Group 1"       "MetaPhlAn taxon names"     
    ## [51] "NCBI Taxonomy IDs"          "State"                     
    ## [53] "Reviewer"

Stripping illformed entries:

``` r

is.study <- !is.na(full.dat[["Study"]])
is.exp <- !is.na(full.dat[["Experiment"]])
full.dat <- full.dat[is.study & is.exp, ]
```

## Curation output

Number of papers and signatures curated:

``` r

pmids <- unique(full.dat[,"PMID"])
length(pmids)
```

    ## [1] 2292

``` r

nrow(full.dat)
```

    ## [1] 16226

### Publication date of the curated papers:

``` r

pmids <- pmids[!is.na(pmids)]
pubyear <- pmid2pubyear(pmids)
head(cbind(pmids, pubyear))
```

``` r

tab <- table(pubyear)
tab <- tab[order(as.integer(names(tab)))]
df <- data.frame(year = names(tab), papers = as.integer(tab))
ggbarplot(df, x = "year", y = "papers", 
          label = TRUE, fill = "steelblue",
          ggtheme = theme_bw())
```

Stripping empty signatures:

``` r

ind1 <- lengths(full.dat[["MetaPhlAn taxon names"]]) > 0
ind2 <- lengths(full.dat[["NCBI Taxonomy IDs"]]) > 0
dat <- full.dat[ind1 & ind2,]
nrow(dat)
```

    ## [1] 16226

Papers containing only empty UP and DOWN signatures (under curation?):

``` r

setdiff(pmids, unique(dat[,"PMID"]))
```

    ## character(0)

### Progress over time:

``` r

dat[,"Curated date"] <- as.character(lubridate::dmy(dat[,"Curated date"]))
plotProgressOverTime(dat)
```

![](BugSigDBStats_files/figure-html/unnamed-chunk-4-1.png)

``` r

plotProgressOverTime(dat, diff = TRUE)
```

![](BugSigDBStats_files/figure-html/unnamed-chunk-4-2.png)

Stratified by curator:

``` r

npc <- stratifyByCurator(dat)
plotCuratorStats(dat, npc)
```

![](BugSigDBStats_files/figure-html/curatorOutput-1.png)

Number of complete and revised signatures: Turned off because it’s way
too long these days

``` r

table(dat[["State"]])
table(dat[,"Revision editor"])
```

## Study stats

### Study design

``` r

spl <- split(dat[["Study"]], dat[["Study design"]])
sds <- lapply(spl, unique)
sort(lengths(sds), decreasing = FALSE)
```

    ##              cross-sectional observational, not case-control,randomized controlled trial 
    ##                                                                                        1 
    ##                                        laboratory experiment,randomized controlled trial 
    ##                                                                                        1 
    ##                                                         meta-analysis,prospective cohort 
    ##                                                                                        1 
    ##                                                 prospective cohort,laboratory experiment 
    ##                                                                                        1 
    ##                                           prospective cohort,randomized controlled trial 
    ##                                                                                        1 
    ##                                                 randomized controlled trial,case-control 
    ##                                                                                        1 
    ## time series / longitudinal observational,cross-sectional observational, not case-control 
    ##                                                                                        1 
    ##                    cross-sectional observational, not case-control,laboratory experiment 
    ##                                                                                        2 
    ##                                                      meta-analysis,laboratory experiment 
    ##                                                                                        2 
    ##                                    time series / longitudinal observational,case-control 
    ##                                                                                        2 
    ##                     time series / longitudinal observational,randomized controlled trial 
    ##                                                                                        2 
    ##                                                          case-control,prospective cohort 
    ##                                                                                        3 
    ##                    laboratory experiment,cross-sectional observational, not case-control 
    ##                                                                                        3 
    ##                                                 laboratory experiment,prospective cohort 
    ##                                                                                        3 
    ##                     randomized controlled trial,time series / longitudinal observational 
    ##                                                                                        3 
    ##                                                       laboratory experiment,case-control 
    ##                                                                                        4 
    ##                           laboratory experiment,time series / longitudinal observational 
    ##                                                                                        4 
    ##                                                               meta-analysis,case-control 
    ##                                                                                        4 
    ##                       prospective cohort,cross-sectional observational, not case-control 
    ##                                                                                        4 
    ##                                        randomized controlled trial,laboratory experiment 
    ##                                                                                        4 
    ##                           time series / longitudinal observational,laboratory experiment 
    ##                                                                                        4 
    ##                              time series / longitudinal observational,prospective cohort 
    ##                                                                                        4 
    ##                                                               case-control,meta-analysis 
    ##                                                                                        5 
    ##                                    case-control,time series / longitudinal observational 
    ##                                                                                        7 
    ##                                                       case-control,laboratory experiment 
    ##                                                                                       10 
    ##                                                                            meta-analysis 
    ##                                                                                       37 
    ##                                                              randomized controlled trial 
    ##                                                                                      102 
    ##                                                 time series / longitudinal observational 
    ##                                                                                      171 
    ##                                                                       prospective cohort 
    ##                                                                                      191 
    ##                                                                    laboratory experiment 
    ##                                                                                      238 
    ##                                          cross-sectional observational, not case-control 
    ##                                                                                      579 
    ##                                                                             case-control 
    ##                                                                                      901

## Experiment stats

Columns of the full dataset that describe experiments:

``` r

# Experiment ID
exp.cols <- c("Study", "Experiment")

# Subjects
sub.cols <- c("Host species",    
              "Location of subjects", 
              "Body site",
              "Condition", 
              "Antibiotics exclusion",
              "Group 0 sample size",
              "Group 1 sample size")

# Lab analysis              
lab.cols <-  c("Sequencing type",
              "16S variable region",
              "Sequencing platform")

# Statistical analysis
stat.cols <-  c("Statistical test",
              "MHT correction",
              "Significance threshold")

# Alpha diversity
div.cols <- c("Pielou",
              "Shannon",
              "Chao1",
              "Simpson", 
              "Inverse Simpson",
              "Richness")
```

Restrict dataset to experiment information:

``` r

exps <- dat[,c(exp.cols, sub.cols, lab.cols, stat.cols, div.cols)]
exps <- unique(exps)
```

### Subjects

Number of experiments for the top 10 categories for each subjects
column:

``` r

sub.tab <- lapply(sub.cols[1:5], tabCol, df = exps, n = 10)
names(sub.tab) <- sub.cols[1:5]
sub.tab
```

    ## $`Host species`
    ## 
    ##           Homo sapiens           Mus musculus      Rattus norvegicus 
    ##                   7979                   1000                    245 
    ##  Sus scrofa domesticus Canis lupus familiaris          Not specified 
    ##                    142                    135                     44 
    ##         Macaca mulatta             Bos taurus             Ovis aries 
    ##                     33                     30                     28 
    ##    Macaca fascicularis 
    ##                     24 
    ## 
    ## $`Location of subjects`
    ## 
    ##                    China United States of America                    Japan 
    ##                     3755                     1421                      362 
    ##              South Korea                  Germany                    Spain 
    ##                      231                      226                      215 
    ##                    Italy                   Canada                  Denmark 
    ##                      212                      193                      193 
    ##                Australia 
    ##                      162 
    ## 
    ## $`Body site`
    ## 
    ##                     Feces                    Saliva                    Vagina 
    ##                      6319                       460                       285 
    ## Subgingival dental plaque                    Caecum               Oral cavity 
    ##                       149                        99                        88 
    ##                     Colon               Nasopharynx                    Rectum 
    ##                        87                        81                        68 
    ##                     Mouth 
    ##                        67 
    ## 
    ## $Condition
    ## 
    ##         Colorectal cancer                      Diet       Parkinson's disease 
    ##                       374                       297                       235 
    ##                   Obesity                  COVID-19          Response to diet 
    ##                       176                       140                       137 
    ##    Response to transplant          Diet measurement Polycystic ovary syndrome 
    ##                       126                       120                       120 
    ## Type II diabetes mellitus 
    ##                        97 
    ## 
    ## $`Antibiotics exclusion`
    ## 
    ##                  3 months                   1 month                  2 months 
    ##                      1344                      1121                       450 
    ##                  6 months                   2 weeks                    1 week 
    ##                       412                       298                        87 
    ## Recent use of antibiotics                   30 days  Currently on antibiotics 
    ##                        55                        49                        48 
    ##                   4 weeks 
    ##                        40

Proportions instead:

``` r

sub.tab <- lapply(sub.cols[1:5], tabCol, df = exps, n = 10, perc = TRUE)
names(sub.tab) <- sub.cols[1:5]
sub.tab
```

    ## $`Host species`
    ## 
    ##           Homo sapiens           Mus musculus      Rattus norvegicus 
    ##                0.79900                0.10000                0.02450 
    ##  Sus scrofa domesticus Canis lupus familiaris          Not specified 
    ##                0.01420                0.01350                0.00441 
    ##         Macaca mulatta             Bos taurus             Ovis aries 
    ##                0.00331                0.00301                0.00281 
    ##    Macaca fascicularis 
    ##                0.00240 
    ## 
    ## $`Location of subjects`
    ## 
    ##                    China United States of America                    Japan 
    ##                   0.3760                   0.1420                   0.0363 
    ##              South Korea                  Germany                    Spain 
    ##                   0.0232                   0.0227                   0.0216 
    ##                    Italy                   Canada                  Denmark 
    ##                   0.0213                   0.0194                   0.0194 
    ##                Australia 
    ##                   0.0162 
    ## 
    ## $`Body site`
    ## 
    ##                     Feces                    Saliva                    Vagina 
    ##                   0.63600                   0.04630                   0.02870 
    ## Subgingival dental plaque                    Caecum               Oral cavity 
    ##                   0.01500                   0.00996                   0.00886 
    ##                     Colon               Nasopharynx                    Rectum 
    ##                   0.00876                   0.00815                   0.00684 
    ##                     Mouth 
    ##                   0.00674 
    ## 
    ## $Condition
    ## 
    ##         Colorectal cancer                      Diet       Parkinson's disease 
    ##                   0.03830                   0.03040                   0.02410 
    ##                   Obesity                  COVID-19          Response to diet 
    ##                   0.01800                   0.01430                   0.01400 
    ##    Response to transplant          Diet measurement Polycystic ovary syndrome 
    ##                   0.01290                   0.01230                   0.01230 
    ## Type II diabetes mellitus 
    ##                   0.00993 
    ## 
    ## $`Antibiotics exclusion`
    ## 
    ##                  3 months                   1 month                  2 months 
    ##                   0.30500                   0.25400                   0.10200 
    ##                  6 months                   2 weeks                    1 week 
    ##                   0.09350                   0.06760                   0.01970 
    ## Recent use of antibiotics                   30 days  Currently on antibiotics 
    ##                   0.01250                   0.01110                   0.01090 
    ##                   4 weeks 
    ##                   0.00908

Sample size:

``` r

ssize <- apply(exps[,sub.cols[6:7]], 2, summary)
ssize
```

    ##         Group 0 sample size Group 1 sample size
    ## Min.                  0.000             1.00000
    ## 1st Qu.              12.000            10.00000
    ## Median               25.000            22.00000
    ## Mean                363.109            67.05831
    ## 3rd Qu.              52.000            43.00000
    ## Max.             308633.000         10413.00000
    ## NA's               1977.000          1971.00000

### Lab analysis

Number of experiments for the top 10 categories for each lab analysis
column:

``` r

lab.tab <- lapply(lab.cols, tabCol, df = exps, n = 10)
names(lab.tab) <- lab.cols
lab.tab
```

    ## $`Sequencing type`
    ## 
    ##        16S        WMS        PCR ITS / ITS2    RNA-Seq        18S 
    ##       7731       1780         88         77          9          5 
    ## 
    ## $`16S variable region`
    ## 
    ##        34         4        12       123        45 123456789       345         3 
    ##      3635      1990       393       255       204       173       154        86 
    ##        56  23456789 
    ##        63        53 
    ## 
    ## $`Sequencing platform`
    ## 
    ##                              Illumina                              Roche454 
    ##                                  8188                                   362 
    ##                           Ion Torrent                               RT-qPCR 
    ##                                   350                                   145 
    ## PacBio Vega (VS)/Revio (RS)/Sequel II                              Nanopore 
    ##                                    88                                    83 
    ##                           MGISEQ-2000                 BGISEQ-500 Sequencing 
    ##                                    61                                    37 
    ##                             DNBSEQ-T7           Human Intestinal Tract Chip 
    ##                                    31                                    31

Proportions instead:

``` r

lab.tab <- lapply(lab.cols, tabCol, df = exps, n = 10, perc = TRUE)
names(lab.tab) <- lab.cols
lab.tab
```

    ## $`Sequencing type`
    ## 
    ##        16S        WMS        PCR ITS / ITS2    RNA-Seq        18S 
    ##   0.798000   0.184000   0.009080   0.007950   0.000929   0.000516 
    ## 
    ## $`16S variable region`
    ## 
    ##        34         4        12       123        45 123456789       345         3 
    ##   0.50500   0.27600   0.05460   0.03540   0.02830   0.02400   0.02140   0.01190 
    ##        56  23456789 
    ##   0.00875   0.00736 
    ## 
    ## $`Sequencing platform`
    ## 
    ##                              Illumina                              Roche454 
    ##                               0.86400                               0.03820 
    ##                           Ion Torrent                               RT-qPCR 
    ##                               0.03690                               0.01530 
    ## PacBio Vega (VS)/Revio (RS)/Sequel II                              Nanopore 
    ##                               0.00928                               0.00876 
    ##                           MGISEQ-2000                 BGISEQ-500 Sequencing 
    ##                               0.00644                               0.00390 
    ##                             DNBSEQ-T7           Human Intestinal Tract Chip 
    ##                               0.00327                               0.00327

### Statistical analysis

Number of experiments for the top 10 categories for each statistical
analysis column:

``` r

# Define the columns to analyze
stat.cols <- c("Statistical test", "MHT correction", "Significance threshold")
# Get top 10 counts
stat.tab.count <- lapply(stat.cols, tabCol, df = exps, n = 10)
names(stat.tab.count) <- stat.cols
# Get top 10 proportions
stat.tab.prop <- lapply(stat.cols, tabCol, df = exps, n = 10, perc = TRUE)
names(stat.tab.prop) <- stat.cols
```

#### Top 10 Counts

``` r

count_df <- data.frame(
  Name = names(stat.tab.count[["Statistical test"]]), 
  Count = as.vector(stat.tab.count[["Statistical test"]])
)
DT::datatable(
  count_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top categories for: Statistical test"
)
```

``` r

count_df <- data.frame(
  Name = names(stat.tab.count[["MHT correction"]]), 
  Count = as.vector(stat.tab.count[["MHT correction"]])
)
DT::datatable(
  count_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top categories for: MHT correction"
)
```

``` r

count_df <- data.frame(
  Value = names(stat.tab.count[["Significance threshold"]]), 
  Count = as.vector(stat.tab.count[["Significance threshold"]])
)
DT::datatable(
  count_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top categories for: Significance threshold"
)
```

#### Top 10 Proportions

``` r

prop_df <- data.frame(
  Name = names(stat.tab.prop[["Statistical test"]]), 
  Proportion = as.vector(stat.tab.prop[["Statistical test"]])
)
DT::datatable(
  prop_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top proportions for: Statistical test"
)
```

``` r

prop_df <- data.frame(
  Name = names(stat.tab.prop[["MHT correction"]]), 
  Proportion = as.vector(stat.tab.prop[["MHT correction"]])
)
DT::datatable(
  prop_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top proportions for: MHT correction"
)
```

``` r

prop_df <- data.frame(
  Value = names(stat.tab.prop[["Significance threshold"]]), 
  Proportion = as.vector(stat.tab.prop[["Significance threshold"]])
)
DT::datatable(
  prop_df,
  rownames = FALSE,
  options = list(pageLength = 10),
  caption = "Top proportions for: Significance threshold"
)
```

### Alpha diversity

Overall distribution:

``` r

apply(exps[,div.cols], 2, table)
```

    ##           Pielou Shannon Chao1 Simpson Inverse Simpson Richness
    ## decreased     98    1093   657     364              87      695
    ## increased     77     809   481     235              72      513
    ## unchanged    389    3083  1551    1255             309     1556

Correspondence of Shannon diversity and Richness:

``` r

table(exps$Shannon, exps$Richness)
```

    ##            
    ##             decreased increased unchanged
    ##   decreased       382        22        84
    ##   increased        15       246        78
    ##   unchanged       158       130      1250

Conditions with consistently increased or decreased alpha diversity:

``` r

tabDiv(exps, "Shannon", "Condition")
```

    ##                                                           increased decreased
    ## Pulmonary tuberculosis                                            4        31
    ## Diet                                                             19        44
    ## Type I diabetes mellitus                                          0        17
    ## Gastric cancer                                                    6        22
    ## Polycystic ovary syndrome                                         6        21
    ## Periodontitis                                                    17         3
    ## Non-alcoholic fatty liver disease                                 5        18
    ## Ulcerative colitis                                                1        14
    ## COVID-19                                                         14        26
    ## Human papilloma virus infection                                  13         1
    ## Lung cancer                                                       2        13
    ## Crohn's disease                                                   0        10
    ## Multiple sclerosis                                                2        12
    ## Obesity                                                           7        17
    ## Systemic inflammatory response syndrome                           5        15
    ## Alzheimer's disease                                               2        11
    ## Chronic constipation                                              9         0
    ## Clostridium difficile infection                                  10         1
    ## Colorectal cancer                                                26        35
    ## HIV infection                                                     3        12
    ## Gingivitis                                                        1         9
    ## Pancreatic carcinoma                                              9         1
    ## Smoking behaviour measurement                                     8         0
    ## Cervical cancer                                                   8         1
    ## Dry eye syndrome                                                  1         8
    ## Heart failure                                                     0         7
    ## Inflammatory bowel disease                                        0         7
    ## Mastitis                                                          0         7
    ## Mycobacterium tuberculosis                                        7         0
    ## Atopic eczema                                                     5        11
    ## Autism spectrum disorder                                          7         1
    ## Cesarean section                                                  7         1
    ## Graves disease                                                    1         7
    ## Parkinson's disease                                              20        14
    ## Response to allogeneic hematopoietic stem cell transplant         0         6
    ## Spontaneous preterm birth                                        13         7
    ## Acute pancreatitis                                                0         5
    ## Age                                                               8        13
    ## Aging                                                             2         7
    ## Bronchiectasis                                                    0         5
    ## Cystic fibrosis                                                   0         5
    ## Epilepsy                                                          5         0
    ## Helminthiasis                                                     5         0
    ## Non-alcoholic steatohepatitis                                     1         6
    ## Oxygen                                                            5         0
    ## Rheumatoid arthritis                                             10         5
    ## Urinary tract infection                                           1         6
    ## Acute lymphoblastic leukemia                                      0         4
    ## Body mass index                                                   0         4
    ## Chronic kidney disease                                            2         6
    ## Chronic periodontitis                                             4         0
    ## Colitis                                                           4         0
    ## Constipation                                                      7         3
    ## Depressive disorder                                               0         4
    ## Endometrial cancer                                                6         2
    ## Ethnic group                                                      4         8
    ## Female infertility                                                4         0
    ## Food allergy                                                      6         2
    ## Gestational diabetes                                              4         0
    ## Human immunodeficiency virus                                      0         4
    ## Insecticide                                                       0         4
    ## Type II diabetes mellitus                                        10         6
    ## Alcohol drinking                                                  3         0
    ## Antimicrobial agent                                               7        10
    ## Atopic asthma                                                     4         1
    ## Bacterial vaginosis                                               4         1
    ## Birth measurement                                                 3         0
    ## Disease recurrence                                                0         3
    ## Extraction protocol                                              23        26
    ## Hematopoietic stem cell                                           0         3
    ## Hypertension                                                      7         4
    ## Hypothyroidism                                                    0         3
    ## Lifestyle measurement                                             5         2
    ## Male homosexuality                                                3         0
    ## Maternal milk                                                     3         6
    ## Mood disorder                                                     0         3
    ## Postmenopausal osteoporosis                                       0         3
    ## Pregnancy                                                         5         2
    ## Response to antibiotic                                            4         7
    ## Response to antiviral drug                                        2         5
    ## Response to immune checkpoint inhibitor                           3         0
    ## SARS-CoV-2-related disease                                        0         3
    ## Schizophrenia                                                     1         4
    ## Traditional Chinese medicine type                                 2         5
    ## Tuberculosis                                                      2         5
    ## Acne                                                              0         2
    ## Adenomatous colon polyp                                           2         0
    ## Age at assessment                                                 3         1
    ## Alcohol use disorder measurement                                  0         2
    ## Breastfeeding duration                                            3         5
    ## Breed                                                             0         2
    ## Cervical glandular intraepithelial neoplasia                      2         0
    ## Cognitive impairment                                              1         3
    ## Dental caries                                                     2         0
    ## Diffuse large B-cell lymphoma                                     0         2
    ## Eczema                                                            0         2
    ## Environmental exposure measurement                                3         5
    ## Environmental factor                                              7         5
    ## Esophageal adenocarcinoma                                         0         2
    ## Hepatic steatosis                                                 0         2
    ## Iron biomarker measurement                                        1         3
    ## Irritable bowel syndrome                                          5         7
    ## Milk allergic reaction                                            2         0
    ## Pancreatic ductal adenocarcinoma                                  0         2
    ## Papillary thyroid carcinoma                                       2         0
    ## Phenotype                                                         2         0
    ## Phenylketonuria                                                   1         3
    ## Population                                                        5         7
    ## Reproductive behaviour measurement                                2         0
    ## Response to anti-tuberculosis drug                                8        10
    ## Response to diet                                                 10         8
    ## Response to immunochemotherapy                                    3         1
    ## Response to opioid                                                3         1
    ## Sampling site                                                     3         1
    ## Simian immunodeficiency virus infection                           0         2
    ## Sleep duration                                                    3         1
    ## Small for gestational age                                         0         2
    ## Smoking status measurement                                        2         0
    ## Species design                                                    9        11
    ## Squamous cell carcinoma                                           2         0
    ## Streptococcus pneumoniae                                          0         2
    ## Stroke                                                            2         0
    ## Systemic lupus erythematosus                                      0         2
    ## Treatment                                                         5         3
    ## Acute respiratory failure                                         6         5
    ## Age-related macular degeneration                                  0         1
    ## Air pollution                                                     7         6
    ## Anxiety disorder                                                  0         1
    ## Attention deficit-hyperactivity disorder                          2         1
    ## Breast cancer                                                     4         5
    ## Chronic fatigue syndrome                                          0         1
    ## Chronic hepatitis B virus infection                               0         1
    ## Chronic obstructive pulmonary disease                             3         2
    ## Clinical treatment                                                1         2
    ## Coccidiosis                                                       0         1
    ## Delivery method                                                   3         4
    ## Diabetes mellitus                                                 0         1
    ## Diarrhea                                                          3         4
    ## Endometriosis                                                     3         4
    ## Esophageal cancer                                                 1         2
    ## Esophageal squamous cell carcinoma                                1         0
    ## Estradiol measurement                                             1         2
    ## Fasting                                                           0         1
    ## Fungal infectious disease                                         2         1
    ## Glaucoma                                                          4         5
    ## Hepatocellular carcinoma                                          4         5
    ## HIV mother to child transmission                                  1         0
    ## Hypertrophy                                                       1         0
    ## Infant                                                            1         0
    ## Liver enzyme measurement                                          0         1
    ## Major depressive disorder                                         0         1
    ## Metastatic colorectal cancer                                      1         2
    ## Myocardial infarction                                             1         2
    ## Obstructive sleep apnea                                           0         1
    ## Oral cavity carcinoma                                             0         1
    ## Oral squamous cell carcinoma                                      3         2
    ## Osteoporosis                                                      1         0
    ## Ovarian cancer                                                    4         3
    ## Overweight body mass index status                                 0         1
    ## Peri-Implantitis                                                  2         1
    ## Pneumonia                                                         1         2
    ## Prediabetes syndrome                                              0         1
    ## Prostate cancer                                                   0         1
    ## Psoriasis                                                         1         0
    ## Respiratory Syncytial Virus Infection                             0         1
    ## Respiratory tract infectious disease                              0         1
    ## Response to exercise                                              1         0
    ## Response to stress                                                1         0
    ## Response to transplant                                           10         9
    ## Response to vaccine                                               1         2
    ## Rosacea                                                           0         1
    ## Sample treatment protocol                                         1         0
    ## Sampling time                                                     4         3
    ## Sleep Disorder                                                    1         0
    ## Social interaction measurement                                    2         1
    ## Socioeconomic status                                              3         4
    ## Stimulus or stress design                                         0         1
    ## Timepoint                                                         1         0
    ## Treatment outcome measurement                                     1         2
    ## Vesicle membrane                                                  3         2
    ## Vitiligo                                                          0         1
    ## Abnormal stool composition                                        0         0
    ## Acute myeloid leukemia                                            1         1
    ## Arthritis                                                         0         0
    ## Asthma                                                            3         3
    ## Biological sex                                                    1         1
    ## Bipolar disorder                                                  0         0
    ## Celiac disease                                                    0         0
    ## Chlamydia trachomatis                                             2         2
    ## Colon carcinoma                                                   1         1
    ## Colorectal adenoma                                                2         2
    ## Contraception                                                     0         0
    ## COVID-19 symptoms measurement                                     0         0
    ## Diet measurement                                                  0         0
    ## Disease progression measurement                                   0         0
    ## Exposure                                                          1         1
    ## Gastric adenocarcinoma                                            0         0
    ## Glioma                                                            0         0
    ## Head and neck squamous cell carcinoma                             0         0
    ## Health study participation                                        3         3
    ## Hepatic fibrosis                                                  0         0
    ## High fat diet                                                     1         1
    ## Hydroxyproline measurement                                        1         1
    ## Hypoxia                                                           4         4
    ## Ischemic stroke                                                   0         0
    ## Lactose intolerance                                               0         0
    ## Lung transplantation                                              2         2
    ## Metabolic process                                                 0         0
    ## Neurodevelopmental delay                                          0         0
    ## Obsessive-compulsive disorder                                     0         0
    ## Oral lichen planus                                                3         3
    ## Oral mucositis                                                    2         2
    ## Oxalate measurement                                               1         1
    ## Pathological complete response                                    0         0
    ## Preeclampsia                                                      0         0
    ## Premature birth                                                   2         2
    ## Psoriasis vulgaris                                                0         0
    ## Response to dietary antigen                                       0         0
    ## Response to ketogenic diet                                        2         2
    ## Sample collection protocol                                        0         0
    ## SARS coronavirus                                                  0         0
    ## Smoking behavior                                                 10        10
    ## Social deprivation,Psychological measurement                      0         0
    ## Transplant outcome measurement                                    0         0
    ## Viral load                                                        0         0
    ## Waist circumference                                               0         0
    ##                                                           unchanged
    ## Pulmonary tuberculosis                                           20
    ## Diet                                                             74
    ## Type I diabetes mellitus                                         11
    ## Gastric cancer                                                   28
    ## Polycystic ovary syndrome                                        38
    ## Periodontitis                                                    16
    ## Non-alcoholic fatty liver disease                                25
    ## Ulcerative colitis                                                5
    ## COVID-19                                                         50
    ## Human papilloma virus infection                                  32
    ## Lung cancer                                                      10
    ## Crohn's disease                                                   8
    ## Multiple sclerosis                                               32
    ## Obesity                                                          76
    ## Systemic inflammatory response syndrome                           4
    ## Alzheimer's disease                                              30
    ## Chronic constipation                                             12
    ## Clostridium difficile infection                                   1
    ## Colorectal cancer                                               105
    ## HIV infection                                                    28
    ## Gingivitis                                                       13
    ## Pancreatic carcinoma                                              6
    ## Smoking behaviour measurement                                     9
    ## Cervical cancer                                                  13
    ## Dry eye syndrome                                                 14
    ## Heart failure                                                     2
    ## Inflammatory bowel disease                                        1
    ## Mastitis                                                          0
    ## Mycobacterium tuberculosis                                        5
    ## Atopic eczema                                                    72
    ## Autism spectrum disorder                                         12
    ## Cesarean section                                                 22
    ## Graves disease                                                    7
    ## Parkinson's disease                                              95
    ## Response to allogeneic hematopoietic stem cell transplant         0
    ## Spontaneous preterm birth                                         8
    ## Acute pancreatitis                                                3
    ## Age                                                              15
    ## Aging                                                             4
    ## Bronchiectasis                                                    0
    ## Cystic fibrosis                                                   3
    ## Epilepsy                                                          5
    ## Helminthiasis                                                     8
    ## Non-alcoholic steatohepatitis                                     4
    ## Oxygen                                                            0
    ## Rheumatoid arthritis                                             15
    ## Urinary tract infection                                          11
    ## Acute lymphoblastic leukemia                                     10
    ## Body mass index                                                   5
    ## Chronic kidney disease                                           11
    ## Chronic periodontitis                                             8
    ## Colitis                                                           2
    ## Constipation                                                     10
    ## Depressive disorder                                               7
    ## Endometrial cancer                                                3
    ## Ethnic group                                                     11
    ## Female infertility                                                3
    ## Food allergy                                                     19
    ## Gestational diabetes                                             41
    ## Human immunodeficiency virus                                      6
    ## Insecticide                                                       4
    ## Type II diabetes mellitus                                        43
    ## Alcohol drinking                                                  2
    ## Antimicrobial agent                                              25
    ## Atopic asthma                                                     7
    ## Bacterial vaginosis                                               0
    ## Birth measurement                                                 4
    ## Disease recurrence                                                7
    ## Extraction protocol                                              21
    ## Hematopoietic stem cell                                           3
    ## Hypertension                                                     13
    ## Hypothyroidism                                                    4
    ## Lifestyle measurement                                            26
    ## Male homosexuality                                                6
    ## Maternal milk                                                     6
    ## Mood disorder                                                     2
    ## Postmenopausal osteoporosis                                       2
    ## Pregnancy                                                        15
    ## Response to antibiotic                                           24
    ## Response to antiviral drug                                       17
    ## Response to immune checkpoint inhibitor                           2
    ## SARS-CoV-2-related disease                                        5
    ## Schizophrenia                                                    24
    ## Traditional Chinese medicine type                                 8
    ## Tuberculosis                                                      0
    ## Acne                                                              3
    ## Adenomatous colon polyp                                           3
    ## Age at assessment                                                 1
    ## Alcohol use disorder measurement                                  5
    ## Breastfeeding duration                                           10
    ## Breed                                                             7
    ## Cervical glandular intraepithelial neoplasia                      9
    ## Cognitive impairment                                             14
    ## Dental caries                                                     4
    ## Diffuse large B-cell lymphoma                                     9
    ## Eczema                                                           10
    ## Environmental exposure measurement                                4
    ## Environmental factor                                             22
    ## Esophageal adenocarcinoma                                         4
    ## Hepatic steatosis                                                 5
    ## Iron biomarker measurement                                        2
    ## Irritable bowel syndrome                                         27
    ## Milk allergic reaction                                            5
    ## Pancreatic ductal adenocarcinoma                                  4
    ## Papillary thyroid carcinoma                                      10
    ## Phenotype                                                        19
    ## Phenylketonuria                                                   4
    ## Population                                                       28
    ## Reproductive behaviour measurement                                3
    ## Response to anti-tuberculosis drug                               13
    ## Response to diet                                                 49
    ## Response to immunochemotherapy                                    3
    ## Response to opioid                                                2
    ## Sampling site                                                     7
    ## Simian immunodeficiency virus infection                           9
    ## Sleep duration                                                    2
    ## Small for gestational age                                         6
    ## Smoking status measurement                                        4
    ## Species design                                                   13
    ## Squamous cell carcinoma                                           4
    ## Streptococcus pneumoniae                                          4
    ## Stroke                                                           16
    ## Systemic lupus erythematosus                                      5
    ## Treatment                                                        10
    ## Acute respiratory failure                                         0
    ## Age-related macular degeneration                                  4
    ## Air pollution                                                     3
    ## Anxiety disorder                                                  7
    ## Attention deficit-hyperactivity disorder                          6
    ## Breast cancer                                                    31
    ## Chronic fatigue syndrome                                          4
    ## Chronic hepatitis B virus infection                              10
    ## Chronic obstructive pulmonary disease                             5
    ## Clinical treatment                                               12
    ## Coccidiosis                                                       4
    ## Delivery method                                                  12
    ## Diabetes mellitus                                                 9
    ## Diarrhea                                                         15
    ## Endometriosis                                                    14
    ## Esophageal cancer                                                 2
    ## Esophageal squamous cell carcinoma                                5
    ## Estradiol measurement                                             2
    ## Fasting                                                           4
    ## Fungal infectious disease                                         4
    ## Glaucoma                                                          6
    ## Hepatocellular carcinoma                                         10
    ## HIV mother to child transmission                                  9
    ## Hypertrophy                                                       4
    ## Infant                                                            4
    ## Liver enzyme measurement                                          5
    ## Major depressive disorder                                         6
    ## Metastatic colorectal cancer                                      3
    ## Myocardial infarction                                             8
    ## Obstructive sleep apnea                                          11
    ## Oral cavity carcinoma                                             7
    ## Oral squamous cell carcinoma                                      3
    ## Osteoporosis                                                      7
    ## Ovarian cancer                                                   35
    ## Overweight body mass index status                                10
    ## Peri-Implantitis                                                  3
    ## Pneumonia                                                         3
    ## Prediabetes syndrome                                              6
    ## Prostate cancer                                                   9
    ## Psoriasis                                                        14
    ## Respiratory Syncytial Virus Infection                             5
    ## Respiratory tract infectious disease                              5
    ## Response to exercise                                              4
    ## Response to stress                                                9
    ## Response to transplant                                           34
    ## Response to vaccine                                               5
    ## Rosacea                                                           5
    ## Sample treatment protocol                                         4
    ## Sampling time                                                     7
    ## Sleep Disorder                                                    4
    ## Social interaction measurement                                    6
    ## Socioeconomic status                                              8
    ## Stimulus or stress design                                         4
    ## Timepoint                                                         5
    ## Treatment outcome measurement                                     3
    ## Vesicle membrane                                                  1
    ## Vitiligo                                                          4
    ## Abnormal stool composition                                        6
    ## Acute myeloid leukemia                                            7
    ## Arthritis                                                         6
    ## Asthma                                                           18
    ## Biological sex                                                   14
    ## Bipolar disorder                                                  5
    ## Celiac disease                                                    6
    ## Chlamydia trachomatis                                             2
    ## Colon carcinoma                                                  10
    ## Colorectal adenoma                                               23
    ## Contraception                                                     5
    ## COVID-19 symptoms measurement                                     5
    ## Diet measurement                                                 60
    ## Disease progression measurement                                   5
    ## Exposure                                                          5
    ## Gastric adenocarcinoma                                            8
    ## Glioma                                                            5
    ## Head and neck squamous cell carcinoma                             8
    ## Health study participation                                       18
    ## Hepatic fibrosis                                                  5
    ## High fat diet                                                    18
    ## Hydroxyproline measurement                                        3
    ## Hypoxia                                                           0
    ## Ischemic stroke                                                   5
    ## Lactose intolerance                                               5
    ## Lung transplantation                                              2
    ## Metabolic process                                                 7
    ## Neurodevelopmental delay                                          6
    ## Obsessive-compulsive disorder                                     5
    ## Oral lichen planus                                                4
    ## Oral mucositis                                                    3
    ## Oxalate measurement                                               8
    ## Pathological complete response                                    7
    ## Preeclampsia                                                      6
    ## Premature birth                                                   6
    ## Psoriasis vulgaris                                               14
    ## Response to dietary antigen                                       6
    ## Response to ketogenic diet                                        3
    ## Sample collection protocol                                        9
    ## SARS coronavirus                                                  6
    ## Smoking behavior                                                 26
    ## Social deprivation,Psychological measurement                     12
    ## Transplant outcome measurement                                   11
    ## Viral load                                                        6
    ## Waist circumference                                               5

``` r

tabDiv(exps, "Shannon", "Condition", perc = TRUE)
```

    ##                                                           increased decreased
    ## Pulmonary tuberculosis                                        0.073     0.560
    ## Diet                                                          0.140     0.320
    ## Type I diabetes mellitus                                      0.000     0.610
    ## Gastric cancer                                                0.110     0.390
    ## Polycystic ovary syndrome                                     0.092     0.320
    ## Periodontitis                                                 0.470     0.083
    ## Non-alcoholic fatty liver disease                             0.100     0.380
    ## Ulcerative colitis                                            0.050     0.700
    ## COVID-19                                                      0.160     0.290
    ## Human papilloma virus infection                               0.280     0.022
    ## Lung cancer                                                   0.080     0.520
    ## Crohn's disease                                               0.000     0.560
    ## Multiple sclerosis                                            0.043     0.260
    ## Obesity                                                       0.070     0.170
    ## Systemic inflammatory response syndrome                       0.210     0.620
    ## Alzheimer's disease                                           0.047     0.260
    ## Chronic constipation                                          0.430     0.000
    ## Clostridium difficile infection                               0.830     0.083
    ## Colorectal cancer                                             0.160     0.210
    ## HIV infection                                                 0.070     0.280
    ## Gingivitis                                                    0.043     0.390
    ## Pancreatic carcinoma                                          0.560     0.062
    ## Smoking behaviour measurement                                 0.470     0.000
    ## Cervical cancer                                               0.360     0.045
    ## Dry eye syndrome                                              0.043     0.350
    ## Heart failure                                                 0.000     0.780
    ## Inflammatory bowel disease                                    0.000     0.880
    ## Mastitis                                                      0.000     1.000
    ## Mycobacterium tuberculosis                                    0.580     0.000
    ## Atopic eczema                                                 0.057     0.120
    ## Autism spectrum disorder                                      0.350     0.050
    ## Cesarean section                                              0.230     0.033
    ## Graves disease                                                0.067     0.470
    ## Parkinson's disease                                           0.160     0.110
    ## Response to allogeneic hematopoietic stem cell transplant     0.000     1.000
    ## Spontaneous preterm birth                                     0.460     0.250
    ## Acute pancreatitis                                            0.000     0.620
    ## Age                                                           0.220     0.360
    ## Aging                                                         0.150     0.540
    ## Bronchiectasis                                                0.000     1.000
    ## Cystic fibrosis                                               0.000     0.620
    ## Epilepsy                                                      0.500     0.000
    ## Helminthiasis                                                 0.380     0.000
    ## Non-alcoholic steatohepatitis                                 0.091     0.550
    ## Oxygen                                                        1.000     0.000
    ## Rheumatoid arthritis                                          0.330     0.170
    ## Urinary tract infection                                       0.056     0.330
    ## Acute lymphoblastic leukemia                                  0.000     0.290
    ## Body mass index                                               0.000     0.440
    ## Chronic kidney disease                                        0.110     0.320
    ## Chronic periodontitis                                         0.330     0.000
    ## Colitis                                                       0.670     0.000
    ## Constipation                                                  0.350     0.150
    ## Depressive disorder                                           0.000     0.360
    ## Endometrial cancer                                            0.550     0.180
    ## Ethnic group                                                  0.170     0.350
    ## Female infertility                                            0.570     0.000
    ## Food allergy                                                  0.220     0.074
    ## Gestational diabetes                                          0.089     0.000
    ## Human immunodeficiency virus                                  0.000     0.400
    ## Insecticide                                                   0.000     0.500
    ## Type II diabetes mellitus                                     0.170     0.100
    ## Alcohol drinking                                              0.600     0.000
    ## Antimicrobial agent                                           0.170     0.240
    ## Atopic asthma                                                 0.330     0.083
    ## Bacterial vaginosis                                           0.800     0.200
    ## Birth measurement                                             0.430     0.000
    ## Disease recurrence                                            0.000     0.300
    ## Extraction protocol                                           0.330     0.370
    ## Hematopoietic stem cell                                       0.000     0.500
    ## Hypertension                                                  0.290     0.170
    ## Hypothyroidism                                                0.000     0.430
    ## Lifestyle measurement                                         0.150     0.061
    ## Male homosexuality                                            0.330     0.000
    ## Maternal milk                                                 0.200     0.400
    ## Mood disorder                                                 0.000     0.600
    ## Postmenopausal osteoporosis                                   0.000     0.600
    ## Pregnancy                                                     0.230     0.091
    ## Response to antibiotic                                        0.110     0.200
    ## Response to antiviral drug                                    0.083     0.210
    ## Response to immune checkpoint inhibitor                       0.600     0.000
    ## SARS-CoV-2-related disease                                    0.000     0.380
    ## Schizophrenia                                                 0.034     0.140
    ## Traditional Chinese medicine type                             0.130     0.330
    ## Tuberculosis                                                  0.290     0.710
    ## Acne                                                          0.000     0.400
    ## Adenomatous colon polyp                                       0.400     0.000
    ## Age at assessment                                             0.600     0.200
    ## Alcohol use disorder measurement                              0.000     0.290
    ## Breastfeeding duration                                        0.170     0.280
    ## Breed                                                         0.000     0.220
    ## Cervical glandular intraepithelial neoplasia                  0.180     0.000
    ## Cognitive impairment                                          0.056     0.170
    ## Dental caries                                                 0.330     0.000
    ## Diffuse large B-cell lymphoma                                 0.000     0.180
    ## Eczema                                                        0.000     0.170
    ## Environmental exposure measurement                            0.250     0.420
    ## Environmental factor                                          0.210     0.150
    ## Esophageal adenocarcinoma                                     0.000     0.330
    ## Hepatic steatosis                                             0.000     0.290
    ## Iron biomarker measurement                                    0.170     0.500
    ## Irritable bowel syndrome                                      0.130     0.180
    ## Milk allergic reaction                                        0.290     0.000
    ## Pancreatic ductal adenocarcinoma                              0.000     0.330
    ## Papillary thyroid carcinoma                                   0.170     0.000
    ## Phenotype                                                     0.095     0.000
    ## Phenylketonuria                                               0.120     0.380
    ## Population                                                    0.120     0.180
    ## Reproductive behaviour measurement                            0.400     0.000
    ## Response to anti-tuberculosis drug                            0.260     0.320
    ## Response to diet                                              0.150     0.120
    ## Response to immunochemotherapy                                0.430     0.140
    ## Response to opioid                                            0.500     0.170
    ## Sampling site                                                 0.270     0.091
    ## Simian immunodeficiency virus infection                       0.000     0.180
    ## Sleep duration                                                0.500     0.170
    ## Small for gestational age                                     0.000     0.250
    ## Smoking status measurement                                    0.330     0.000
    ## Species design                                                0.270     0.330
    ## Squamous cell carcinoma                                       0.330     0.000
    ## Streptococcus pneumoniae                                      0.000     0.330
    ## Stroke                                                        0.110     0.000
    ## Systemic lupus erythematosus                                  0.000     0.290
    ## Treatment                                                     0.280     0.170
    ## Acute respiratory failure                                     0.550     0.450
    ## Age-related macular degeneration                              0.000     0.200
    ## Air pollution                                                 0.440     0.380
    ## Anxiety disorder                                              0.000     0.120
    ## Attention deficit-hyperactivity disorder                      0.220     0.110
    ## Breast cancer                                                 0.100     0.120
    ## Chronic fatigue syndrome                                      0.000     0.200
    ## Chronic hepatitis B virus infection                           0.000     0.091
    ## Chronic obstructive pulmonary disease                         0.300     0.200
    ## Clinical treatment                                            0.067     0.130
    ## Coccidiosis                                                   0.000     0.200
    ## Delivery method                                               0.160     0.210
    ## Diabetes mellitus                                             0.000     0.100
    ## Diarrhea                                                      0.140     0.180
    ## Endometriosis                                                 0.140     0.190
    ## Esophageal cancer                                             0.200     0.400
    ## Esophageal squamous cell carcinoma                            0.170     0.000
    ## Estradiol measurement                                         0.200     0.400
    ## Fasting                                                       0.000     0.200
    ## Fungal infectious disease                                     0.290     0.140
    ## Glaucoma                                                      0.270     0.330
    ## Hepatocellular carcinoma                                      0.210     0.260
    ## HIV mother to child transmission                              0.100     0.000
    ## Hypertrophy                                                   0.200     0.000
    ## Infant                                                        0.200     0.000
    ## Liver enzyme measurement                                      0.000     0.170
    ## Major depressive disorder                                     0.000     0.140
    ## Metastatic colorectal cancer                                  0.170     0.330
    ## Myocardial infarction                                         0.091     0.180
    ## Obstructive sleep apnea                                       0.000     0.083
    ## Oral cavity carcinoma                                         0.000     0.120
    ## Oral squamous cell carcinoma                                  0.380     0.250
    ## Osteoporosis                                                  0.120     0.000
    ## Ovarian cancer                                                0.095     0.071
    ## Overweight body mass index status                             0.000     0.091
    ## Peri-Implantitis                                              0.330     0.170
    ## Pneumonia                                                     0.170     0.330
    ## Prediabetes syndrome                                          0.000     0.140
    ## Prostate cancer                                               0.000     0.100
    ## Psoriasis                                                     0.067     0.000
    ## Respiratory Syncytial Virus Infection                         0.000     0.170
    ## Respiratory tract infectious disease                          0.000     0.170
    ## Response to exercise                                          0.200     0.000
    ## Response to stress                                            0.100     0.000
    ## Response to transplant                                        0.190     0.170
    ## Response to vaccine                                           0.120     0.250
    ## Rosacea                                                       0.000     0.170
    ## Sample treatment protocol                                     0.200     0.000
    ## Sampling time                                                 0.290     0.210
    ## Sleep Disorder                                                0.200     0.000
    ## Social interaction measurement                                0.220     0.110
    ## Socioeconomic status                                          0.200     0.270
    ## Stimulus or stress design                                     0.000     0.200
    ## Timepoint                                                     0.170     0.000
    ## Treatment outcome measurement                                 0.170     0.330
    ## Vesicle membrane                                              0.500     0.330
    ## Vitiligo                                                      0.000     0.200
    ## Abnormal stool composition                                    0.000     0.000
    ## Acute myeloid leukemia                                        0.110     0.110
    ## Arthritis                                                     0.000     0.000
    ## Asthma                                                        0.120     0.120
    ## Biological sex                                                0.062     0.062
    ## Bipolar disorder                                              0.000     0.000
    ## Celiac disease                                                0.000     0.000
    ## Chlamydia trachomatis                                         0.330     0.330
    ## Colon carcinoma                                               0.083     0.083
    ## Colorectal adenoma                                            0.074     0.074
    ## Contraception                                                 0.000     0.000
    ## COVID-19 symptoms measurement                                 0.000     0.000
    ## Diet measurement                                              0.000     0.000
    ## Disease progression measurement                               0.000     0.000
    ## Exposure                                                      0.140     0.140
    ## Gastric adenocarcinoma                                        0.000     0.000
    ## Glioma                                                        0.000     0.000
    ## Head and neck squamous cell carcinoma                         0.000     0.000
    ## Health study participation                                    0.120     0.120
    ## Hepatic fibrosis                                              0.000     0.000
    ## High fat diet                                                 0.050     0.050
    ## Hydroxyproline measurement                                    0.200     0.200
    ## Hypoxia                                                       0.500     0.500
    ## Ischemic stroke                                               0.000     0.000
    ## Lactose intolerance                                           0.000     0.000
    ## Lung transplantation                                          0.330     0.330
    ## Metabolic process                                             0.000     0.000
    ## Neurodevelopmental delay                                      0.000     0.000
    ## Obsessive-compulsive disorder                                 0.000     0.000
    ## Oral lichen planus                                            0.300     0.300
    ## Oral mucositis                                                0.290     0.290
    ## Oxalate measurement                                           0.100     0.100
    ## Pathological complete response                                0.000     0.000
    ## Preeclampsia                                                  0.000     0.000
    ## Premature birth                                               0.200     0.200
    ## Psoriasis vulgaris                                            0.000     0.000
    ## Response to dietary antigen                                   0.000     0.000
    ## Response to ketogenic diet                                    0.290     0.290
    ## Sample collection protocol                                    0.000     0.000
    ## SARS coronavirus                                              0.000     0.000
    ## Smoking behavior                                              0.220     0.220
    ## Social deprivation,Psychological measurement                  0.000     0.000
    ## Transplant outcome measurement                                0.000     0.000
    ## Viral load                                                    0.000     0.000
    ## Waist circumference                                           0.000     0.000
    ##                                                           unchanged
    ## Pulmonary tuberculosis                                        0.360
    ## Diet                                                          0.540
    ## Type I diabetes mellitus                                      0.390
    ## Gastric cancer                                                0.500
    ## Polycystic ovary syndrome                                     0.580
    ## Periodontitis                                                 0.440
    ## Non-alcoholic fatty liver disease                             0.520
    ## Ulcerative colitis                                            0.250
    ## COVID-19                                                      0.560
    ## Human papilloma virus infection                               0.700
    ## Lung cancer                                                   0.400
    ## Crohn's disease                                               0.440
    ## Multiple sclerosis                                            0.700
    ## Obesity                                                       0.760
    ## Systemic inflammatory response syndrome                       0.170
    ## Alzheimer's disease                                           0.700
    ## Chronic constipation                                          0.570
    ## Clostridium difficile infection                               0.083
    ## Colorectal cancer                                             0.630
    ## HIV infection                                                 0.650
    ## Gingivitis                                                    0.570
    ## Pancreatic carcinoma                                          0.380
    ## Smoking behaviour measurement                                 0.530
    ## Cervical cancer                                               0.590
    ## Dry eye syndrome                                              0.610
    ## Heart failure                                                 0.220
    ## Inflammatory bowel disease                                    0.120
    ## Mastitis                                                      0.000
    ## Mycobacterium tuberculosis                                    0.420
    ## Atopic eczema                                                 0.820
    ## Autism spectrum disorder                                      0.600
    ## Cesarean section                                              0.730
    ## Graves disease                                                0.470
    ## Parkinson's disease                                           0.740
    ## Response to allogeneic hematopoietic stem cell transplant     0.000
    ## Spontaneous preterm birth                                     0.290
    ## Acute pancreatitis                                            0.380
    ## Age                                                           0.420
    ## Aging                                                         0.310
    ## Bronchiectasis                                                0.000
    ## Cystic fibrosis                                               0.380
    ## Epilepsy                                                      0.500
    ## Helminthiasis                                                 0.620
    ## Non-alcoholic steatohepatitis                                 0.360
    ## Oxygen                                                        0.000
    ## Rheumatoid arthritis                                          0.500
    ## Urinary tract infection                                       0.610
    ## Acute lymphoblastic leukemia                                  0.710
    ## Body mass index                                               0.560
    ## Chronic kidney disease                                        0.580
    ## Chronic periodontitis                                         0.670
    ## Colitis                                                       0.330
    ## Constipation                                                  0.500
    ## Depressive disorder                                           0.640
    ## Endometrial cancer                                            0.270
    ## Ethnic group                                                  0.480
    ## Female infertility                                            0.430
    ## Food allergy                                                  0.700
    ## Gestational diabetes                                          0.910
    ## Human immunodeficiency virus                                  0.600
    ## Insecticide                                                   0.500
    ## Type II diabetes mellitus                                     0.730
    ## Alcohol drinking                                              0.400
    ## Antimicrobial agent                                           0.600
    ## Atopic asthma                                                 0.580
    ## Bacterial vaginosis                                           0.000
    ## Birth measurement                                             0.570
    ## Disease recurrence                                            0.700
    ## Extraction protocol                                           0.300
    ## Hematopoietic stem cell                                       0.500
    ## Hypertension                                                  0.540
    ## Hypothyroidism                                                0.570
    ## Lifestyle measurement                                         0.790
    ## Male homosexuality                                            0.670
    ## Maternal milk                                                 0.400
    ## Mood disorder                                                 0.400
    ## Postmenopausal osteoporosis                                   0.400
    ## Pregnancy                                                     0.680
    ## Response to antibiotic                                        0.690
    ## Response to antiviral drug                                    0.710
    ## Response to immune checkpoint inhibitor                       0.400
    ## SARS-CoV-2-related disease                                    0.620
    ## Schizophrenia                                                 0.830
    ## Traditional Chinese medicine type                             0.530
    ## Tuberculosis                                                  0.000
    ## Acne                                                          0.600
    ## Adenomatous colon polyp                                       0.600
    ## Age at assessment                                             0.200
    ## Alcohol use disorder measurement                              0.710
    ## Breastfeeding duration                                        0.560
    ## Breed                                                         0.780
    ## Cervical glandular intraepithelial neoplasia                  0.820
    ## Cognitive impairment                                          0.780
    ## Dental caries                                                 0.670
    ## Diffuse large B-cell lymphoma                                 0.820
    ## Eczema                                                        0.830
    ## Environmental exposure measurement                            0.330
    ## Environmental factor                                          0.650
    ## Esophageal adenocarcinoma                                     0.670
    ## Hepatic steatosis                                             0.710
    ## Iron biomarker measurement                                    0.330
    ## Irritable bowel syndrome                                      0.690
    ## Milk allergic reaction                                        0.710
    ## Pancreatic ductal adenocarcinoma                              0.670
    ## Papillary thyroid carcinoma                                   0.830
    ## Phenotype                                                     0.900
    ## Phenylketonuria                                               0.500
    ## Population                                                    0.700
    ## Reproductive behaviour measurement                            0.600
    ## Response to anti-tuberculosis drug                            0.420
    ## Response to diet                                              0.730
    ## Response to immunochemotherapy                                0.430
    ## Response to opioid                                            0.330
    ## Sampling site                                                 0.640
    ## Simian immunodeficiency virus infection                       0.820
    ## Sleep duration                                                0.330
    ## Small for gestational age                                     0.750
    ## Smoking status measurement                                    0.670
    ## Species design                                                0.390
    ## Squamous cell carcinoma                                       0.670
    ## Streptococcus pneumoniae                                      0.670
    ## Stroke                                                        0.890
    ## Systemic lupus erythematosus                                  0.710
    ## Treatment                                                     0.560
    ## Acute respiratory failure                                     0.000
    ## Age-related macular degeneration                              0.800
    ## Air pollution                                                 0.190
    ## Anxiety disorder                                              0.880
    ## Attention deficit-hyperactivity disorder                      0.670
    ## Breast cancer                                                 0.780
    ## Chronic fatigue syndrome                                      0.800
    ## Chronic hepatitis B virus infection                           0.910
    ## Chronic obstructive pulmonary disease                         0.500
    ## Clinical treatment                                            0.800
    ## Coccidiosis                                                   0.800
    ## Delivery method                                               0.630
    ## Diabetes mellitus                                             0.900
    ## Diarrhea                                                      0.680
    ## Endometriosis                                                 0.670
    ## Esophageal cancer                                             0.400
    ## Esophageal squamous cell carcinoma                            0.830
    ## Estradiol measurement                                         0.400
    ## Fasting                                                       0.800
    ## Fungal infectious disease                                     0.570
    ## Glaucoma                                                      0.400
    ## Hepatocellular carcinoma                                      0.530
    ## HIV mother to child transmission                              0.900
    ## Hypertrophy                                                   0.800
    ## Infant                                                        0.800
    ## Liver enzyme measurement                                      0.830
    ## Major depressive disorder                                     0.860
    ## Metastatic colorectal cancer                                  0.500
    ## Myocardial infarction                                         0.730
    ## Obstructive sleep apnea                                       0.920
    ## Oral cavity carcinoma                                         0.880
    ## Oral squamous cell carcinoma                                  0.380
    ## Osteoporosis                                                  0.880
    ## Ovarian cancer                                                0.830
    ## Overweight body mass index status                             0.910
    ## Peri-Implantitis                                              0.500
    ## Pneumonia                                                     0.500
    ## Prediabetes syndrome                                          0.860
    ## Prostate cancer                                               0.900
    ## Psoriasis                                                     0.930
    ## Respiratory Syncytial Virus Infection                         0.830
    ## Respiratory tract infectious disease                          0.830
    ## Response to exercise                                          0.800
    ## Response to stress                                            0.900
    ## Response to transplant                                        0.640
    ## Response to vaccine                                           0.620
    ## Rosacea                                                       0.830
    ## Sample treatment protocol                                     0.800
    ## Sampling time                                                 0.500
    ## Sleep Disorder                                                0.800
    ## Social interaction measurement                                0.670
    ## Socioeconomic status                                          0.530
    ## Stimulus or stress design                                     0.800
    ## Timepoint                                                     0.830
    ## Treatment outcome measurement                                 0.500
    ## Vesicle membrane                                              0.170
    ## Vitiligo                                                      0.800
    ## Abnormal stool composition                                    1.000
    ## Acute myeloid leukemia                                        0.780
    ## Arthritis                                                     1.000
    ## Asthma                                                        0.750
    ## Biological sex                                                0.880
    ## Bipolar disorder                                              1.000
    ## Celiac disease                                                1.000
    ## Chlamydia trachomatis                                         0.330
    ## Colon carcinoma                                               0.830
    ## Colorectal adenoma                                            0.850
    ## Contraception                                                 1.000
    ## COVID-19 symptoms measurement                                 1.000
    ## Diet measurement                                              1.000
    ## Disease progression measurement                               1.000
    ## Exposure                                                      0.710
    ## Gastric adenocarcinoma                                        1.000
    ## Glioma                                                        1.000
    ## Head and neck squamous cell carcinoma                         1.000
    ## Health study participation                                    0.750
    ## Hepatic fibrosis                                              1.000
    ## High fat diet                                                 0.900
    ## Hydroxyproline measurement                                    0.600
    ## Hypoxia                                                       0.000
    ## Ischemic stroke                                               1.000
    ## Lactose intolerance                                           1.000
    ## Lung transplantation                                          0.330
    ## Metabolic process                                             1.000
    ## Neurodevelopmental delay                                      1.000
    ## Obsessive-compulsive disorder                                 1.000
    ## Oral lichen planus                                            0.400
    ## Oral mucositis                                                0.430
    ## Oxalate measurement                                           0.800
    ## Pathological complete response                                1.000
    ## Preeclampsia                                                  1.000
    ## Premature birth                                               0.600
    ## Psoriasis vulgaris                                            1.000
    ## Response to dietary antigen                                   1.000
    ## Response to ketogenic diet                                    0.430
    ## Sample collection protocol                                    1.000
    ## SARS coronavirus                                              1.000
    ## Smoking behavior                                              0.570
    ## Social deprivation,Psychological measurement                  1.000
    ## Transplant outcome measurement                                1.000
    ## Viral load                                                    1.000
    ## Waist circumference                                           1.000

``` r

tabDiv(exps, "Richness", "Condition")
```

    ##                                                           increased decreased
    ## Pulmonary tuberculosis                                            6        22
    ## Diet                                                              6        20
    ## Helminthiasis                                                    13         0
    ## Periodontitis                                                    16         3
    ## HIV infection                                                     3        15
    ## Multiple sclerosis                                                2        14
    ## Polycystic ovary syndrome                                         9        18
    ## Aging                                                             0         8
    ## Parkinson's disease                                              19        27
    ## Phenotype                                                         9         1
    ## Small for gestational age                                         0         8
    ## Chronic constipation                                              7         0
    ## Gastric cancer                                                    5        12
    ## Gingivitis                                                        1         8
    ## High fat diet                                                     0         7
    ## Opioid dependence,Opioid use disorder                             7         0
    ## Chronic kidney disease                                            0         6
    ## Colorectal cancer                                                22        28
    ## Increased intestinal transit time                                 6         0
    ## Mastitis                                                          0         6
    ## Response to allogeneic hematopoietic stem cell transplant         0         6
    ## Schizophrenia                                                     1         7
    ## Alcohol drinking                                                  5         0
    ## Antimicrobial agent                                               2         7
    ## COVID-19                                                         15        20
    ## Environmental factor                                              0         5
    ## Human immunodeficiency virus                                      1         6
    ## Human papilloma virus infection                                   7         2
    ## Species design                                                    7        12
    ## Systemic lupus erythematosus                                      5         0
    ## Type II diabetes mellitus                                         8         3
    ## Acute lymphoblastic leukemia                                      5         1
    ## Air pollution                                                     9         5
    ## Cervical glandular intraepithelial neoplasia                      4         0
    ## Colon carcinoma                                                   0         4
    ## Crohn's disease                                                   2         6
    ## Diarrhea                                                          5         1
    ## Dry eye syndrome                                                  0         4
    ## Epilepsy                                                          4         0
    ## Fungal infectious disease                                         0         4
    ## Inflammatory bowel disease                                        2         6
    ## Lung cancer                                                       0         4
    ## Major depressive disorder                                         4         0
    ## Maternal milk                                                     2         6
    ## Non-alcoholic steatohepatitis                                     1         5
    ## Rheumatoid arthritis                                              8         4
    ## Smoking behavior                                                  6        10
    ## Stress-related disorder                                           0         4
    ## Ulcerative colitis                                                0         4
    ## Vesicle membrane                                                  5         1
    ## Alzheimer's disease                                               8         5
    ## Atopic asthma                                                     4         1
    ## Delivery method                                                   4         1
    ## Depressive disorder                                               0         3
    ## Endometriosis                                                     4         1
    ## Food allergy                                                      0         3
    ## Graves disease                                                    2         5
    ## Hematopoietic stem cell                                           0         3
    ## Hypertension                                                      1         4
    ## Hypertrophy                                                       3         0
    ## Iron biomarker measurement                                        1         4
    ## Oral squamous cell carcinoma                                      1         4
    ## Pneumonia                                                         0         3
    ## Pregnancy                                                         1         4
    ## Response to diet                                                  6         9
    ## Tuberculosis                                                      2         5
    ## Acute myeloid leukemia                                            0         2
    ## Age                                                               7         5
    ## Asthma                                                            2         4
    ## Bone mineral content measurement                                  0         2
    ## Breast cancer                                                     2         0
    ## Cognitive impairment                                              0         2
    ## Esophageal adenocarcinoma                                         0         2
    ## Exercise                                                          3         1
    ## Gestational diabetes                                              4         6
    ## Health study participation                                        2         0
    ## Hypothyroidism                                                    0         2
    ## Lifestyle measurement                                             4         2
    ## Metabolic process                                                 0         2
    ## Non-alcoholic fatty liver disease                                 0         2
    ## Obesity                                                          10         8
    ## Osteoporosis                                                      2         0
    ## Overweight body mass index status                                 2         4
    ## Phenylketonuria                                                   1         3
    ## Population                                                        2         0
    ## Response to antiviral drug                                        0         2
    ## Simian immunodeficiency virus infection                           0         2
    ## Smoking behaviour measurement                                     2         0
    ## Smoking status measurement                                        2         0
    ## Streptococcus pneumoniae                                          0         2
    ## Traditional Chinese medicine type                                 1         3
    ## Transplant outcome measurement                                    0         2
    ## Treatment                                                         5         7
    ## Atopic eczema                                                     2         1
    ## Autism spectrum disorder                                          5         6
    ## Biological sex                                                    0         1
    ## Bone density                                                      0         1
    ## Breastfeeding duration                                            2         3
    ## Cervical cancer                                                   1         2
    ## Cesarean section                                                  3         2
    ## Chronic hepatitis B virus infection                               0         1
    ## Chronic periodontitis                                             1         0
    ## Coccidiosis                                                       0         1
    ## Colorectal adenoma                                                1         2
    ## Endometrial cancer                                                1         2
    ## Ethnic group                                                      3         2
    ## Glioma                                                            1         2
    ## Heart failure                                                     2         3
    ## Hepatic steatosis                                                 0         1
    ## HIV mother to child transmission                                  1         0
    ## Infant                                                            1         0
    ## Ischemic stroke                                                   2         1
    ## Neurodevelopmental delay                                          1         0
    ## Obsessive-compulsive disorder                                     0         1
    ## Ovarian cancer                                                    1         0
    ## Psoriasis                                                         0         1
    ## Reproductive behaviour measurement                                1         0
    ## Response to antibiotic                                            5         4
    ## Response to stress                                                1         0
    ## Response to transplant                                            9         8
    ## Rosacea                                                           0         1
    ## Sampling site                                                     1         2
    ## Socioeconomic status                                              2         1
    ## Transport                                                         1         2
    ## Treatment outcome measurement                                     2         3
    ## Type I diabetes mellitus                                          0         1
    ## Abnormal stool composition                                        0         0
    ## Attention deficit-hyperactivity disorder                          0         0
    ## Bacterial vaginosis                                               1         1
    ## Chlamydia trachomatis                                             1         1
    ## Constipation                                                      7         7
    ## Diet measurement                                                  0         0
    ## Head and neck squamous cell carcinoma                             0         0
    ## Hepatocellular carcinoma                                          3         3
    ## Irritable bowel syndrome                                          5         5
    ## Male homosexuality                                                0         0
    ## Myocardial infarction                                             0         0
    ## Papillary thyroid carcinoma                                       0         0
    ## Physical activity                                                 2         2
    ## Premature birth                                                   2         2
    ## Prostate cancer                                                   0         0
    ## Psoriasis vulgaris                                                0         0
    ## Sample collection protocol                                        0         0
    ## Sleep Disorder                                                    2         2
    ## Stroke                                                            2         2
    ## Urinary tract infection                                           1         1
    ## Viral load                                                        0         0
    ##                                                           unchanged
    ## Pulmonary tuberculosis                                            9
    ## Diet                                                             34
    ## Helminthiasis                                                     0
    ## Periodontitis                                                    23
    ## HIV infection                                                    10
    ## Multiple sclerosis                                               29
    ## Polycystic ovary syndrome                                        19
    ## Aging                                                             2
    ## Parkinson's disease                                              33
    ## Phenotype                                                        11
    ## Small for gestational age                                         0
    ## Chronic constipation                                              7
    ## Gastric cancer                                                   15
    ## Gingivitis                                                        8
    ## High fat diet                                                     8
    ## Opioid dependence,Opioid use disorder                             0
    ## Chronic kidney disease                                            7
    ## Colorectal cancer                                                52
    ## Increased intestinal transit time                                 0
    ## Mastitis                                                          1
    ## Response to allogeneic hematopoietic stem cell transplant         0
    ## Schizophrenia                                                    14
    ## Alcohol drinking                                                  0
    ## Antimicrobial agent                                              10
    ## COVID-19                                                         26
    ## Environmental factor                                             17
    ## Human immunodeficiency virus                                      2
    ## Human papilloma virus infection                                  13
    ## Species design                                                   14
    ## Systemic lupus erythematosus                                      5
    ## Type II diabetes mellitus                                        25
    ## Acute lymphoblastic leukemia                                      1
    ## Air pollution                                                     6
    ## Cervical glandular intraepithelial neoplasia                      2
    ## Colon carcinoma                                                  10
    ## Crohn's disease                                                   4
    ## Diarrhea                                                          3
    ## Dry eye syndrome                                                  4
    ## Epilepsy                                                          3
    ## Fungal infectious disease                                         2
    ## Inflammatory bowel disease                                        1
    ## Lung cancer                                                      10
    ## Major depressive disorder                                         3
    ## Maternal milk                                                     1
    ## Non-alcoholic steatohepatitis                                     1
    ## Rheumatoid arthritis                                              3
    ## Smoking behavior                                                  8
    ## Stress-related disorder                                           1
    ## Ulcerative colitis                                                2
    ## Vesicle membrane                                                  0
    ## Alzheimer's disease                                              27
    ## Atopic asthma                                                     7
    ## Delivery method                                                   8
    ## Depressive disorder                                               4
    ## Endometriosis                                                     9
    ## Food allergy                                                      9
    ## Graves disease                                                    6
    ## Hematopoietic stem cell                                           3
    ## Hypertension                                                     12
    ## Hypertrophy                                                       2
    ## Iron biomarker measurement                                        1
    ## Oral squamous cell carcinoma                                      0
    ## Pneumonia                                                         2
    ## Pregnancy                                                        10
    ## Response to diet                                                 20
    ## Tuberculosis                                                      0
    ## Acute myeloid leukemia                                            4
    ## Age                                                               3
    ## Asthma                                                           13
    ## Bone mineral content measurement                                  8
    ## Breast cancer                                                    17
    ## Cognitive impairment                                              7
    ## Esophageal adenocarcinoma                                         4
    ## Exercise                                                          2
    ## Gestational diabetes                                             31
    ## Health study participation                                       11
    ## Hypothyroidism                                                    4
    ## Lifestyle measurement                                             3
    ## Metabolic process                                                 5
    ## Non-alcoholic fatty liver disease                                14
    ## Obesity                                                          29
    ## Osteoporosis                                                      4
    ## Overweight body mass index status                                 3
    ## Phenylketonuria                                                   4
    ## Population                                                        4
    ## Response to antiviral drug                                       13
    ## Simian immunodeficiency virus infection                           4
    ## Smoking behaviour measurement                                     5
    ## Smoking status measurement                                        5
    ## Streptococcus pneumoniae                                          3
    ## Traditional Chinese medicine type                                 4
    ## Transplant outcome measurement                                    5
    ## Treatment                                                        11
    ## Atopic eczema                                                     6
    ## Autism spectrum disorder                                          0
    ## Biological sex                                                   11
    ## Bone density                                                      5
    ## Breastfeeding duration                                            8
    ## Cervical cancer                                                   4
    ## Cesarean section                                                 12
    ## Chronic hepatitis B virus infection                               7
    ## Chronic periodontitis                                             6
    ## Coccidiosis                                                       4
    ## Colorectal adenoma                                               11
    ## Endometrial cancer                                                3
    ## Ethnic group                                                      1
    ## Glioma                                                            2
    ## Heart failure                                                     4
    ## Hepatic steatosis                                                 6
    ## HIV mother to child transmission                                  9
    ## Infant                                                            4
    ## Ischemic stroke                                                   3
    ## Neurodevelopmental delay                                          5
    ## Obsessive-compulsive disorder                                     4
    ## Ovarian cancer                                                   39
    ## Psoriasis                                                         8
    ## Reproductive behaviour measurement                                4
    ## Response to antibiotic                                            9
    ## Response to stress                                                9
    ## Response to transplant                                           14
    ## Rosacea                                                           4
    ## Sampling site                                                     2
    ## Socioeconomic status                                              2
    ## Transport                                                         3
    ## Treatment outcome measurement                                     1
    ## Type I diabetes mellitus                                          4
    ## Abnormal stool composition                                        6
    ## Attention deficit-hyperactivity disorder                          5
    ## Bacterial vaginosis                                               3
    ## Chlamydia trachomatis                                             3
    ## Constipation                                                     12
    ## Diet measurement                                                  7
    ## Head and neck squamous cell carcinoma                            11
    ## Hepatocellular carcinoma                                          0
    ## Irritable bowel syndrome                                         16
    ## Male homosexuality                                                9
    ## Myocardial infarction                                             6
    ## Papillary thyroid carcinoma                                      12
    ## Physical activity                                                 2
    ## Premature birth                                                   1
    ## Prostate cancer                                                   7
    ## Psoriasis vulgaris                                               14
    ## Sample collection protocol                                        9
    ## Sleep Disorder                                                    3
    ## Stroke                                                           17
    ## Urinary tract infection                                           9
    ## Viral load                                                        5

``` r

tabDiv(exps, "Richness", "Condition", perc = TRUE)
```

    ##                                                           increased decreased
    ## Pulmonary tuberculosis                                        0.160     0.590
    ## Diet                                                          0.100     0.330
    ## Helminthiasis                                                 1.000     0.000
    ## Periodontitis                                                 0.380     0.071
    ## HIV infection                                                 0.110     0.540
    ## Multiple sclerosis                                            0.044     0.310
    ## Polycystic ovary syndrome                                     0.200     0.390
    ## Aging                                                         0.000     0.800
    ## Parkinson's disease                                           0.240     0.340
    ## Phenotype                                                     0.430     0.048
    ## Small for gestational age                                     0.000     1.000
    ## Chronic constipation                                          0.500     0.000
    ## Gastric cancer                                                0.160     0.380
    ## Gingivitis                                                    0.059     0.470
    ## High fat diet                                                 0.000     0.470
    ## Opioid dependence,Opioid use disorder                         1.000     0.000
    ## Chronic kidney disease                                        0.000     0.460
    ## Colorectal cancer                                             0.220     0.270
    ## Increased intestinal transit time                             1.000     0.000
    ## Mastitis                                                      0.000     0.860
    ## Response to allogeneic hematopoietic stem cell transplant     0.000     1.000
    ## Schizophrenia                                                 0.045     0.320
    ## Alcohol drinking                                              1.000     0.000
    ## Antimicrobial agent                                           0.110     0.370
    ## COVID-19                                                      0.250     0.330
    ## Environmental factor                                          0.000     0.230
    ## Human immunodeficiency virus                                  0.110     0.670
    ## Human papilloma virus infection                               0.320     0.091
    ## Species design                                                0.210     0.360
    ## Systemic lupus erythematosus                                  0.500     0.000
    ## Type II diabetes mellitus                                     0.220     0.083
    ## Acute lymphoblastic leukemia                                  0.710     0.140
    ## Air pollution                                                 0.450     0.250
    ## Cervical glandular intraepithelial neoplasia                  0.670     0.000
    ## Colon carcinoma                                               0.000     0.290
    ## Crohn's disease                                               0.170     0.500
    ## Diarrhea                                                      0.560     0.110
    ## Dry eye syndrome                                              0.000     0.500
    ## Epilepsy                                                      0.570     0.000
    ## Fungal infectious disease                                     0.000     0.670
    ## Inflammatory bowel disease                                    0.220     0.670
    ## Lung cancer                                                   0.000     0.290
    ## Major depressive disorder                                     0.570     0.000
    ## Maternal milk                                                 0.220     0.670
    ## Non-alcoholic steatohepatitis                                 0.140     0.710
    ## Rheumatoid arthritis                                          0.530     0.270
    ## Smoking behavior                                              0.250     0.420
    ## Stress-related disorder                                       0.000     0.800
    ## Ulcerative colitis                                            0.000     0.670
    ## Vesicle membrane                                              0.830     0.170
    ## Alzheimer's disease                                           0.200     0.120
    ## Atopic asthma                                                 0.330     0.083
    ## Delivery method                                               0.310     0.077
    ## Depressive disorder                                           0.000     0.430
    ## Endometriosis                                                 0.290     0.071
    ## Food allergy                                                  0.000     0.250
    ## Graves disease                                                0.150     0.380
    ## Hematopoietic stem cell                                       0.000     0.500
    ## Hypertension                                                  0.059     0.240
    ## Hypertrophy                                                   0.600     0.000
    ## Iron biomarker measurement                                    0.170     0.670
    ## Oral squamous cell carcinoma                                  0.200     0.800
    ## Pneumonia                                                     0.000     0.600
    ## Pregnancy                                                     0.067     0.270
    ## Response to diet                                              0.170     0.260
    ## Tuberculosis                                                  0.290     0.710
    ## Acute myeloid leukemia                                        0.000     0.330
    ## Age                                                           0.470     0.330
    ## Asthma                                                        0.110     0.210
    ## Bone mineral content measurement                              0.000     0.200
    ## Breast cancer                                                 0.110     0.000
    ## Cognitive impairment                                          0.000     0.220
    ## Esophageal adenocarcinoma                                     0.000     0.330
    ## Exercise                                                      0.500     0.170
    ## Gestational diabetes                                          0.098     0.150
    ## Health study participation                                    0.150     0.000
    ## Hypothyroidism                                                0.000     0.330
    ## Lifestyle measurement                                         0.440     0.220
    ## Metabolic process                                             0.000     0.290
    ## Non-alcoholic fatty liver disease                             0.000     0.120
    ## Obesity                                                       0.210     0.170
    ## Osteoporosis                                                  0.330     0.000
    ## Overweight body mass index status                             0.220     0.440
    ## Phenylketonuria                                               0.120     0.380
    ## Population                                                    0.330     0.000
    ## Response to antiviral drug                                    0.000     0.130
    ## Simian immunodeficiency virus infection                       0.000     0.330
    ## Smoking behaviour measurement                                 0.290     0.000
    ## Smoking status measurement                                    0.290     0.000
    ## Streptococcus pneumoniae                                      0.000     0.400
    ## Traditional Chinese medicine type                             0.120     0.380
    ## Transplant outcome measurement                                0.000     0.290
    ## Treatment                                                     0.220     0.300
    ## Atopic eczema                                                 0.220     0.110
    ## Autism spectrum disorder                                      0.450     0.550
    ## Biological sex                                                0.000     0.083
    ## Bone density                                                  0.000     0.170
    ## Breastfeeding duration                                        0.150     0.230
    ## Cervical cancer                                               0.140     0.290
    ## Cesarean section                                              0.180     0.120
    ## Chronic hepatitis B virus infection                           0.000     0.120
    ## Chronic periodontitis                                         0.140     0.000
    ## Coccidiosis                                                   0.000     0.200
    ## Colorectal adenoma                                            0.071     0.140
    ## Endometrial cancer                                            0.170     0.330
    ## Ethnic group                                                  0.500     0.330
    ## Glioma                                                        0.200     0.400
    ## Heart failure                                                 0.220     0.330
    ## Hepatic steatosis                                             0.000     0.140
    ## HIV mother to child transmission                              0.100     0.000
    ## Infant                                                        0.200     0.000
    ## Ischemic stroke                                               0.330     0.170
    ## Neurodevelopmental delay                                      0.170     0.000
    ## Obsessive-compulsive disorder                                 0.000     0.200
    ## Ovarian cancer                                                0.025     0.000
    ## Psoriasis                                                     0.000     0.110
    ## Reproductive behaviour measurement                            0.200     0.000
    ## Response to antibiotic                                        0.280     0.220
    ## Response to stress                                            0.100     0.000
    ## Response to transplant                                        0.290     0.260
    ## Rosacea                                                       0.000     0.200
    ## Sampling site                                                 0.200     0.400
    ## Socioeconomic status                                          0.400     0.200
    ## Transport                                                     0.170     0.330
    ## Treatment outcome measurement                                 0.330     0.500
    ## Type I diabetes mellitus                                      0.000     0.200
    ## Abnormal stool composition                                    0.000     0.000
    ## Attention deficit-hyperactivity disorder                      0.000     0.000
    ## Bacterial vaginosis                                           0.200     0.200
    ## Chlamydia trachomatis                                         0.200     0.200
    ## Constipation                                                  0.270     0.270
    ## Diet measurement                                              0.000     0.000
    ## Head and neck squamous cell carcinoma                         0.000     0.000
    ## Hepatocellular carcinoma                                      0.500     0.500
    ## Irritable bowel syndrome                                      0.190     0.190
    ## Male homosexuality                                            0.000     0.000
    ## Myocardial infarction                                         0.000     0.000
    ## Papillary thyroid carcinoma                                   0.000     0.000
    ## Physical activity                                             0.330     0.330
    ## Premature birth                                               0.400     0.400
    ## Prostate cancer                                               0.000     0.000
    ## Psoriasis vulgaris                                            0.000     0.000
    ## Sample collection protocol                                    0.000     0.000
    ## Sleep Disorder                                                0.290     0.290
    ## Stroke                                                        0.095     0.095
    ## Urinary tract infection                                       0.091     0.091
    ## Viral load                                                    0.000     0.000
    ##                                                           unchanged
    ## Pulmonary tuberculosis                                         0.24
    ## Diet                                                           0.57
    ## Helminthiasis                                                  0.00
    ## Periodontitis                                                  0.55
    ## HIV infection                                                  0.36
    ## Multiple sclerosis                                             0.64
    ## Polycystic ovary syndrome                                      0.41
    ## Aging                                                          0.20
    ## Parkinson's disease                                            0.42
    ## Phenotype                                                      0.52
    ## Small for gestational age                                      0.00
    ## Chronic constipation                                           0.50
    ## Gastric cancer                                                 0.47
    ## Gingivitis                                                     0.47
    ## High fat diet                                                  0.53
    ## Opioid dependence,Opioid use disorder                          0.00
    ## Chronic kidney disease                                         0.54
    ## Colorectal cancer                                              0.51
    ## Increased intestinal transit time                              0.00
    ## Mastitis                                                       0.14
    ## Response to allogeneic hematopoietic stem cell transplant      0.00
    ## Schizophrenia                                                  0.64
    ## Alcohol drinking                                               0.00
    ## Antimicrobial agent                                            0.53
    ## COVID-19                                                       0.43
    ## Environmental factor                                           0.77
    ## Human immunodeficiency virus                                   0.22
    ## Human papilloma virus infection                                0.59
    ## Species design                                                 0.42
    ## Systemic lupus erythematosus                                   0.50
    ## Type II diabetes mellitus                                      0.69
    ## Acute lymphoblastic leukemia                                   0.14
    ## Air pollution                                                  0.30
    ## Cervical glandular intraepithelial neoplasia                   0.33
    ## Colon carcinoma                                                0.71
    ## Crohn's disease                                                0.33
    ## Diarrhea                                                       0.33
    ## Dry eye syndrome                                               0.50
    ## Epilepsy                                                       0.43
    ## Fungal infectious disease                                      0.33
    ## Inflammatory bowel disease                                     0.11
    ## Lung cancer                                                    0.71
    ## Major depressive disorder                                      0.43
    ## Maternal milk                                                  0.11
    ## Non-alcoholic steatohepatitis                                  0.14
    ## Rheumatoid arthritis                                           0.20
    ## Smoking behavior                                               0.33
    ## Stress-related disorder                                        0.20
    ## Ulcerative colitis                                             0.33
    ## Vesicle membrane                                               0.00
    ## Alzheimer's disease                                            0.68
    ## Atopic asthma                                                  0.58
    ## Delivery method                                                0.62
    ## Depressive disorder                                            0.57
    ## Endometriosis                                                  0.64
    ## Food allergy                                                   0.75
    ## Graves disease                                                 0.46
    ## Hematopoietic stem cell                                        0.50
    ## Hypertension                                                   0.71
    ## Hypertrophy                                                    0.40
    ## Iron biomarker measurement                                     0.17
    ## Oral squamous cell carcinoma                                   0.00
    ## Pneumonia                                                      0.40
    ## Pregnancy                                                      0.67
    ## Response to diet                                               0.57
    ## Tuberculosis                                                   0.00
    ## Acute myeloid leukemia                                         0.67
    ## Age                                                            0.20
    ## Asthma                                                         0.68
    ## Bone mineral content measurement                               0.80
    ## Breast cancer                                                  0.89
    ## Cognitive impairment                                           0.78
    ## Esophageal adenocarcinoma                                      0.67
    ## Exercise                                                       0.33
    ## Gestational diabetes                                           0.76
    ## Health study participation                                     0.85
    ## Hypothyroidism                                                 0.67
    ## Lifestyle measurement                                          0.33
    ## Metabolic process                                              0.71
    ## Non-alcoholic fatty liver disease                              0.88
    ## Obesity                                                        0.62
    ## Osteoporosis                                                   0.67
    ## Overweight body mass index status                              0.33
    ## Phenylketonuria                                                0.50
    ## Population                                                     0.67
    ## Response to antiviral drug                                     0.87
    ## Simian immunodeficiency virus infection                        0.67
    ## Smoking behaviour measurement                                  0.71
    ## Smoking status measurement                                     0.71
    ## Streptococcus pneumoniae                                       0.60
    ## Traditional Chinese medicine type                              0.50
    ## Transplant outcome measurement                                 0.71
    ## Treatment                                                      0.48
    ## Atopic eczema                                                  0.67
    ## Autism spectrum disorder                                       0.00
    ## Biological sex                                                 0.92
    ## Bone density                                                   0.83
    ## Breastfeeding duration                                         0.62
    ## Cervical cancer                                                0.57
    ## Cesarean section                                               0.71
    ## Chronic hepatitis B virus infection                            0.88
    ## Chronic periodontitis                                          0.86
    ## Coccidiosis                                                    0.80
    ## Colorectal adenoma                                             0.79
    ## Endometrial cancer                                             0.50
    ## Ethnic group                                                   0.17
    ## Glioma                                                         0.40
    ## Heart failure                                                  0.44
    ## Hepatic steatosis                                              0.86
    ## HIV mother to child transmission                               0.90
    ## Infant                                                         0.80
    ## Ischemic stroke                                                0.50
    ## Neurodevelopmental delay                                       0.83
    ## Obsessive-compulsive disorder                                  0.80
    ## Ovarian cancer                                                 0.98
    ## Psoriasis                                                      0.89
    ## Reproductive behaviour measurement                             0.80
    ## Response to antibiotic                                         0.50
    ## Response to stress                                             0.90
    ## Response to transplant                                         0.45
    ## Rosacea                                                        0.80
    ## Sampling site                                                  0.40
    ## Socioeconomic status                                           0.40
    ## Transport                                                      0.50
    ## Treatment outcome measurement                                  0.17
    ## Type I diabetes mellitus                                       0.80
    ## Abnormal stool composition                                     1.00
    ## Attention deficit-hyperactivity disorder                       1.00
    ## Bacterial vaginosis                                            0.60
    ## Chlamydia trachomatis                                          0.60
    ## Constipation                                                   0.46
    ## Diet measurement                                               1.00
    ## Head and neck squamous cell carcinoma                          1.00
    ## Hepatocellular carcinoma                                       0.00
    ## Irritable bowel syndrome                                       0.62
    ## Male homosexuality                                             1.00
    ## Myocardial infarction                                          1.00
    ## Papillary thyroid carcinoma                                    1.00
    ## Physical activity                                              0.33
    ## Premature birth                                                0.20
    ## Prostate cancer                                                1.00
    ## Psoriasis vulgaris                                             1.00
    ## Sample collection protocol                                     1.00
    ## Sleep Disorder                                                 0.43
    ## Stroke                                                         0.81
    ## Urinary tract infection                                        0.82
    ## Viral load                                                     1.00

Body sites with consistently increased or decreased alpha diversity:

``` r

tabDiv(exps, "Shannon", "Body site")
```

    ##                                                 increased decreased unchanged
    ## Feces                                                 404       684      1873
    ## Vagina                                                 41        12        57
    ## Sputum                                                  7        23        17
    ## Posterior fornix of vagina                             13         0        10
    ## Gastrointestinal system mucosa                          0        11         0
    ## Blood serum                                             9         0         0
    ## Buccal mucosa                                          12         3         8
    ## Oral cavity                                            15         6        30
    ## Stomach                                                 5        14         5
    ## Subgingival dental plaque                              17         8        44
    ## Saliva                                                 42        50       186
    ## Uterine cervix                                          9         1        21
    ## Uterine cervix,Vaginal fluid                            9         1         0
    ## Ileum                                                   3        10        14
    ## Skin of body                                            8        15         8
    ## Supragingival dental plaque                             3        10        14
    ## Meconium                                                8         2        18
    ## Lung                                                    4         9        14
    ## Pharynx                                                 6         1         2
    ## Semen                                                   0         5         3
    ## Space surrounding organism                              2         7        13
    ## Tongue                                                  0         5        16
    ## Ascending colon                                         4         0         7
    ## Axilla skin                                             5         1        12
    ## Caecum                                                  4         8        35
    ## Cervical cavity                                         4         0         3
    ## Colorectal mucosa                                       0         4        17
    ## Colorectum                                              4         0         1
    ## Throat                                                  0         4        11
    ## Vaginal fluid                                           4         0        11
    ## Ascending colon,Colorectal mucosa,Sigmoid colon         3         0         2
    ## Duodenum                                                0         3         6
    ## feces                                                   1         4        15
    ## Intestine                                               1         4        20
    ## Rectum                                                  3         6        22
    ## Skin of forearm                                         3         0         3
    ## Sputum,Feces                                            1         4         1
    ## Urine                                                   4         1        20
    ## Bile                                                    2         0         3
    ## Brachialis muscle                                       0         2         3
    ## Cecum mucosa                                            2         4         6
    ## Colon                                                   6         8        29
    ## Conjunctiva                                             1         3        12
    ## Conjunctival sac                                        1         3         3
    ## Esophagus                                               0         2         4
    ## Forelimb skin                                           2         0         4
    ## Milk                                                    3         5        24
    ## Mouth                                                   8         6        28
    ## Oropharynx                                              1         3         5
    ## Pancreas                                                0         2         4
    ## Rumen                                                   2         0         6
    ## Thyroid gland                                           2         0        10
    ## Uterus                                                  3         1        11
    ## Blood                                                   1         2        14
    ## Body proper,Insect leg                                  0         1         5
    ## Breast                                                  3         4        10
    ## Breast,Milk                                             1         0         4
    ## Bulbar conjunctiva                                      3         2         5
    ## Dental plaque                                           3         4         6
    ## Digestive tract                                         2         1         2
    ## Gingival groove                                         2         1         6
    ## Intestinal mucosa                                       1         2        11
    ## Liver                                                   4         3         2
    ## Mouth mucosa                                            1         2         2
    ## Mucosa of stomach                                       1         0         4
    ## Small intestine                                         3         4         1
    ## Vagina,Uterine cervix                                   3         2         7
    ## Bronchus                                                0         0         6
    ## Colonic mucosa                                          0         0         5
    ## Duodenal mucosa                                         0         0         6
    ## Endothelium of trachea                                  3         3         0
    ## Jejunum                                                 1         1         8
    ## Nasal cavity                                            1         1         5
    ## Nasopharynx                                             6         6        40
    ## Ovary                                                   0         0         7
    ## Peritoneal fluid                                        0         0         6
    ## Posterior wall of oropharynx                            2         2         1
    ## Skin epidermis                                          2         2         4
    ## Skin of abdomen                                         0         0         5
    ## Skin of face                                            1         1         3
    ## Surface of tongue                                       2         2         3
    ## Ventral side of post-anal tail                          0         0         6

``` r

tabDiv(exps, "Shannon", "Body site", perc = TRUE)
```

    ##                                                 increased decreased unchanged
    ## Feces                                               0.140     0.230      0.63
    ## Vagina                                              0.370     0.110      0.52
    ## Sputum                                              0.150     0.490      0.36
    ## Posterior fornix of vagina                          0.570     0.000      0.43
    ## Gastrointestinal system mucosa                      0.000     1.000      0.00
    ## Blood serum                                         1.000     0.000      0.00
    ## Buccal mucosa                                       0.520     0.130      0.35
    ## Oral cavity                                         0.290     0.120      0.59
    ## Stomach                                             0.210     0.580      0.21
    ## Subgingival dental plaque                           0.250     0.120      0.64
    ## Saliva                                              0.150     0.180      0.67
    ## Uterine cervix                                      0.290     0.032      0.68
    ## Uterine cervix,Vaginal fluid                        0.900     0.100      0.00
    ## Ileum                                               0.110     0.370      0.52
    ## Skin of body                                        0.260     0.480      0.26
    ## Supragingival dental plaque                         0.110     0.370      0.52
    ## Meconium                                            0.290     0.071      0.64
    ## Lung                                                0.150     0.330      0.52
    ## Pharynx                                             0.670     0.110      0.22
    ## Semen                                               0.000     0.620      0.38
    ## Space surrounding organism                          0.091     0.320      0.59
    ## Tongue                                              0.000     0.240      0.76
    ## Ascending colon                                     0.360     0.000      0.64
    ## Axilla skin                                         0.280     0.056      0.67
    ## Caecum                                              0.085     0.170      0.74
    ## Cervical cavity                                     0.570     0.000      0.43
    ## Colorectal mucosa                                   0.000     0.190      0.81
    ## Colorectum                                          0.800     0.000      0.20
    ## Throat                                              0.000     0.270      0.73
    ## Vaginal fluid                                       0.270     0.000      0.73
    ## Ascending colon,Colorectal mucosa,Sigmoid colon     0.600     0.000      0.40
    ## Duodenum                                            0.000     0.330      0.67
    ## feces                                               0.050     0.200      0.75
    ## Intestine                                           0.040     0.160      0.80
    ## Rectum                                              0.097     0.190      0.71
    ## Skin of forearm                                     0.500     0.000      0.50
    ## Sputum,Feces                                        0.170     0.670      0.17
    ## Urine                                               0.160     0.040      0.80
    ## Bile                                                0.400     0.000      0.60
    ## Brachialis muscle                                   0.000     0.400      0.60
    ## Cecum mucosa                                        0.170     0.330      0.50
    ## Colon                                               0.140     0.190      0.67
    ## Conjunctiva                                         0.062     0.190      0.75
    ## Conjunctival sac                                    0.140     0.430      0.43
    ## Esophagus                                           0.000     0.330      0.67
    ## Forelimb skin                                       0.330     0.000      0.67
    ## Milk                                                0.094     0.160      0.75
    ## Mouth                                               0.190     0.140      0.67
    ## Oropharynx                                          0.110     0.330      0.56
    ## Pancreas                                            0.000     0.330      0.67
    ## Rumen                                               0.250     0.000      0.75
    ## Thyroid gland                                       0.170     0.000      0.83
    ## Uterus                                              0.200     0.067      0.73
    ## Blood                                               0.059     0.120      0.82
    ## Body proper,Insect leg                              0.000     0.170      0.83
    ## Breast                                              0.180     0.240      0.59
    ## Breast,Milk                                         0.200     0.000      0.80
    ## Bulbar conjunctiva                                  0.300     0.200      0.50
    ## Dental plaque                                       0.230     0.310      0.46
    ## Digestive tract                                     0.400     0.200      0.40
    ## Gingival groove                                     0.220     0.110      0.67
    ## Intestinal mucosa                                   0.071     0.140      0.79
    ## Liver                                               0.440     0.330      0.22
    ## Mouth mucosa                                        0.200     0.400      0.40
    ## Mucosa of stomach                                   0.200     0.000      0.80
    ## Small intestine                                     0.380     0.500      0.12
    ## Vagina,Uterine cervix                               0.250     0.170      0.58
    ## Bronchus                                            0.000     0.000      1.00
    ## Colonic mucosa                                      0.000     0.000      1.00
    ## Duodenal mucosa                                     0.000     0.000      1.00
    ## Endothelium of trachea                              0.500     0.500      0.00
    ## Jejunum                                             0.100     0.100      0.80
    ## Nasal cavity                                        0.140     0.140      0.71
    ## Nasopharynx                                         0.120     0.120      0.77
    ## Ovary                                               0.000     0.000      1.00
    ## Peritoneal fluid                                    0.000     0.000      1.00
    ## Posterior wall of oropharynx                        0.400     0.400      0.20
    ## Skin epidermis                                      0.250     0.250      0.50
    ## Skin of abdomen                                     0.000     0.000      1.00
    ## Skin of face                                        0.200     0.200      0.60
    ## Surface of tongue                                   0.290     0.290      0.43
    ## Ventral side of post-anal tail                      0.000     0.000      1.00

``` r

tabDiv(exps, "Richness", "Body site")
```

    ##                              increased decreased unchanged
    ## Feces                              273       430       951
    ## Oral cavity                         16         4        22
    ## Sputum                               0        12         5
    ## Vagina                              18         8        28
    ## Posterior fornix of vagina          10         1         2
    ## Supragingival dental plaque          0         9         9
    ## Ileum                                3        11        10
    ## Stomach                              4        12         3
    ## Mouth                               10         3         9
    ## Oropharynx                           0         6         6
    ## Rectum                               1         7         8
    ## Subgingival dental plaque           11         5        32
    ## Uterine cervix                       8         2        11
    ## Skin epidermis                       0         5         2
    ## Skin of body                         4         9         6
    ## Uterine cervix,Vaginal fluid         7         2         1
    ## Conjunctival sac                     0         4         1
    ## Throat                               1         5         5
    ## Ascending colon                      3         0         6
    ## Colon                               10         7        17
    ## Lung                                 0         3         4
    ## Saliva                              27        30       102
    ## Small intestine                      1         4         0
    ## Vaginal fluid                        3         0         2
    ## Cecum mucosa                         2         4         1
    ## Colorectal mucosa                    2         0         9
    ## Dental plaque                        1         3         6
    ## Ear                                  2         0         3
    ## Esophagus                            0         2         4
    ## Surface of tongue                    4         2         1
    ## Tongue                               2         4         7
    ## Urine                                4         2        16
    ## Blood                                1         2         2
    ## Buccal mucosa                        3         2         6
    ## Caecum                               5         6        12
    ## Intestine                            2         1        16
    ## Meconium                             3         4        11
    ## Nasal cavity                         1         2         6
    ## Nasopharynx                          8         9        26
    ## Semen                                2         3         3
    ## Vagina,Uterine cervix                1         0        11
    ## Breast                               1         1         7
    ## Bronchus                             0         0         6
    ## Conjunctiva                          1         1         5
    ## feces                                1         1         4
    ## Liver                                3         3         0
    ## Milk                                 4         4         5
    ## Ovary                                0         0         7
    ## Peritoneal fluid                     0         0         6
    ## Rumen                                0         0         6
    ## Thyroid gland                        0         0        12

``` r

tabDiv(exps, "Richness", "Body site", perc = TRUE)
```

    ##                              increased decreased unchanged
    ## Feces                            0.170     0.260      0.57
    ## Oral cavity                      0.380     0.095      0.52
    ## Sputum                           0.000     0.710      0.29
    ## Vagina                           0.330     0.150      0.52
    ## Posterior fornix of vagina       0.770     0.077      0.15
    ## Supragingival dental plaque      0.000     0.500      0.50
    ## Ileum                            0.120     0.460      0.42
    ## Stomach                          0.210     0.630      0.16
    ## Mouth                            0.450     0.140      0.41
    ## Oropharynx                       0.000     0.500      0.50
    ## Rectum                           0.062     0.440      0.50
    ## Subgingival dental plaque        0.230     0.100      0.67
    ## Uterine cervix                   0.380     0.095      0.52
    ## Skin epidermis                   0.000     0.710      0.29
    ## Skin of body                     0.210     0.470      0.32
    ## Uterine cervix,Vaginal fluid     0.700     0.200      0.10
    ## Conjunctival sac                 0.000     0.800      0.20
    ## Throat                           0.091     0.450      0.45
    ## Ascending colon                  0.330     0.000      0.67
    ## Colon                            0.290     0.210      0.50
    ## Lung                             0.000     0.430      0.57
    ## Saliva                           0.170     0.190      0.64
    ## Small intestine                  0.200     0.800      0.00
    ## Vaginal fluid                    0.600     0.000      0.40
    ## Cecum mucosa                     0.290     0.570      0.14
    ## Colorectal mucosa                0.180     0.000      0.82
    ## Dental plaque                    0.100     0.300      0.60
    ## Ear                              0.400     0.000      0.60
    ## Esophagus                        0.000     0.330      0.67
    ## Surface of tongue                0.570     0.290      0.14
    ## Tongue                           0.150     0.310      0.54
    ## Urine                            0.180     0.091      0.73
    ## Blood                            0.200     0.400      0.40
    ## Buccal mucosa                    0.270     0.180      0.55
    ## Caecum                           0.220     0.260      0.52
    ## Intestine                        0.110     0.053      0.84
    ## Meconium                         0.170     0.220      0.61
    ## Nasal cavity                     0.110     0.220      0.67
    ## Nasopharynx                      0.190     0.210      0.60
    ## Semen                            0.250     0.380      0.38
    ## Vagina,Uterine cervix            0.083     0.000      0.92
    ## Breast                           0.110     0.110      0.78
    ## Bronchus                         0.000     0.000      1.00
    ## Conjunctiva                      0.140     0.140      0.71
    ## feces                            0.170     0.170      0.67
    ## Liver                            0.500     0.500      0.00
    ## Milk                             0.310     0.310      0.38
    ## Ovary                            0.000     0.000      1.00
    ## Peritoneal fluid                 0.000     0.000      1.00
    ## Rumen                            0.000     0.000      1.00
    ## Thyroid gland                    0.000     0.000      1.00

## Signature stats

``` r

sigs <- bugsigdbr::getSignatures(dat, tax.id.type = "metaphlan")
```

### Unique microbes

Number unique microbes contained in the signatures:

``` r

(nuniq <- length(unique(unlist(sigs))))
```

    ## [1] 8356

Development of unique microbes captured over time:

``` r

plotUniqueMicrobesOverTime(dat)
```

![](BugSigDBStats_files/figure-html/unnamed-chunk-26-1.png)

### Microbe set size distribution

``` r

summary(lengths(sigs))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    1.00    2.00    4.00    7.98    9.00  467.00

``` r

gghistogram(lengths(sigs), bins = 30, ylab = "number of signatures",
    xlab = "signature size", fill = "#00AFBB", ggtheme = theme_bw())
```

![](BugSigDBStats_files/figure-html/setSize-1.png)

``` r

sum(lengths(sigs) > 4)
```

    ## [1] 7437

### Microbe co-occurrence

``` r

dat.feces <- subset(dat, `Body site` == "Feces")
cooc.mat <- microbeHeatmap(dat.feces, tax.level = "genus", anno = "genus")
```

    ## Loading required namespace: safe

![](BugSigDBStats_files/figure-html/coocc-1.png)

``` r

antag.mat <- microbeHeatmap(dat.feces, tax.level = "genus", anno = "genus", antagonistic = TRUE)
```

![](BugSigDBStats_files/figure-html/coocc-2.png)

Get the top 20 genera most frequently reported as differentially
abundant:

``` r

sigs.feces <- getSignatures(dat.feces, tax.id.type = "taxname", 
                            tax.level = "genus", exact.tax.level = FALSE) 
top20 <- sort(table(unlist(sigs.feces)), decreasing = TRUE)[1:20]
top20
```

    ## 
    ##      Bacteroides  Bifidobacterium Faecalibacterium          Blautia 
    ##              963              698              690              644 
    ##    Lactobacillus    Streptococcus        Roseburia     Ruminococcus 
    ##              597              592              561              554 
    ##      Clostridium       Prevotella  Parabacteroides        Alistipes 
    ##              550              545              510              469 
    ##      Akkermansia            Dorea      Coprococcus      Veillonella 
    ##              409              399              374              348 
    ##     Anaerostipes      Lachnospira     Enterococcus        Dialister 
    ##              344              337              334              331

Subset heatmaps to the top 20 genera most frequently reported as
differentially abundant:

``` r

all(names(top20) %in% rownames(cooc.mat))
```

    ## [1] TRUE

``` r

cooc.mat <- cooc.mat[names(top20), names(top20)]
all(names(top20) %in% rownames(antag.mat))
```

    ## [1] TRUE

``` r

antag.mat <- antag.mat[names(top20), names(top20)]
```

Distinguish by direction of abundance change (increased / decreased):

``` r

# increased
sub.dat.feces <- subset(dat.feces, `Abundance in Group 1` == "increased")
sigs.feces.up <- getSignatures(sub.dat.feces, tax.id.type = "taxname", 
                               tax.level = "genus", exact.tax.level = FALSE) 
top20.up <- table(unlist(sigs.feces.up))[names(top20)]
top20.up
```

    ## 
    ##      Bacteroides  Bifidobacterium Faecalibacterium          Blautia 
    ##              485              372              242              301 
    ##    Lactobacillus    Streptococcus        Roseburia     Ruminococcus 
    ##              377              369              193              248 
    ##      Clostridium       Prevotella  Parabacteroides        Alistipes 
    ##              279              265              296              213 
    ##      Akkermansia            Dorea      Coprococcus      Veillonella 
    ##              265              180              156              209 
    ##     Anaerostipes      Lachnospira     Enterococcus        Dialister 
    ##              162              125              234              148

``` r

# decreased
sub.dat.feces <- subset(dat.feces, `Abundance in Group 1` == "decreased")
sigs.feces.down <- getSignatures(sub.dat.feces, tax.id.type = "taxname", 
                                 tax.level = "genus", exact.tax.level = FALSE) 
top20.down <- table(unlist(sigs.feces.down))[names(top20)]
top20.down
```

    ## 
    ##      Bacteroides  Bifidobacterium Faecalibacterium          Blautia 
    ##              473              321              444              339 
    ##    Lactobacillus    Streptococcus        Roseburia     Ruminococcus 
    ##              218              218              364              302 
    ##      Clostridium       Prevotella  Parabacteroides        Alistipes 
    ##              266              277              210              252 
    ##      Akkermansia            Dorea      Coprococcus      Veillonella 
    ##              141              215              214              136 
    ##     Anaerostipes      Lachnospira     Enterococcus        Dialister 
    ##              178              208               98              181

Plot the heatmap

``` r

# annotation
mat <- matrix(nc = 2, cbind(top20.up, top20.down))
bp <- ComplexHeatmap::anno_barplot(mat, gp = gpar(fill = c("#D55E00", "#0072B2"),
                                                  col = c("#D55E00", "#0072B2")),
                                   height = unit(2, "cm"))
banno <- ComplexHeatmap::HeatmapAnnotation(`Abundance in Group 1` = bp)

lgd_list <- list(
    Legend(labels = c("increased", "decreased"), 
           title = "Abundance in Group 1", 
           type = "grid",
           legend_gp = gpar(col = c("#D55E00", "#0072B2"), fill = c("#D55E00", "#0072B2"))))
                                            
# same direction
# lcm <- sweep(cooc.mat, 2, matrixStats::colMaxs(cooc.mat), FUN = "/")
# we need to dampen the maximum here a bit down,
# otherwise 100% self co-occurrence takes up a large fraction of the colorscale,
sec <- apply(cooc.mat, 2, function(x) sort(x, decreasing = TRUE)[2])
cooc.mat2 <- cooc.mat
for(i in 1:ncol(cooc.mat2)) cooc.mat2[i,i] <- min(cooc.mat2[i,i], 1.4 * sec[i])
lcm <- sweep(cooc.mat2, 2, matrixStats::colMaxs(cooc.mat2), FUN = "/")
col <- circlize::colorRamp2(c(0,1), c("#EEEEEE", "red"))
ht1 <- ComplexHeatmap::Heatmap(lcm,
                               col = col, 
                               name = "Relative frequency (top)",
                               cluster_columns = FALSE, 
                               row_km = 3, 
                               row_title = "same direction", 
                               column_names_rot = 45,
                               row_names_gp = gpar(fontsize = 8),
                               column_names_gp = gpar(fontsize = 8))

# opposite direction
acm <- sweep(antag.mat, 2, matrixStats::colMaxs(antag.mat), FUN = "/")
col <- circlize::colorRamp2(c(0,1), c("#EEEEEE", "blue"))
ht2 <- ComplexHeatmap::Heatmap(acm,
                               col = col, 
                               name = "Relative frequency (bottom)",
                               cluster_columns = FALSE, 
                               row_title = "opposite direction", 
                               row_km = 3, 
                               column_names_rot = 45,
                               row_names_gp = gpar(fontsize = 8),
                               column_names_gp = gpar(fontsize = 8))

# phylum
sfp <- bugsigdbr::getSignatures(dat.feces, tax.id.type = "metaphlan", 
                                tax.level = "genus", exact.tax.level = FALSE) 
sfp20 <- sort(table(unlist(sfp)), decreasing = TRUE)[1:20]
uanno <- bugsigdbr::extractTaxLevel(names(sfp20),
                                    tax.id.type = "taxname",
                                    tax.level = "mixed",
                                    exact.tax.level = FALSE) 
phyla.grid <- seq_along(unique(uanno))
panno <- ComplexHeatmap::HeatmapAnnotation(phylum = uanno)

uanno <- matrix(uanno, nrow = 1)
colnames(uanno) <- names(top20)
pcols <- c("#CC79A7", "#F0E442", "#009E73", "#56B4E9", "#E69F00")
uanno <- ComplexHeatmap::Heatmap(uanno, name = "Phylum",
                                 col = pcols[phyla.grid],
                                 cluster_columns = FALSE,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize = 8))
```

    ## There are 20 unique colors in the vector `col` and 20 unique values in
    ## `matrix`. `Heatmap()` will treat it as an exact discrete one-to-one
    ## mapping. If this is not what you want, slightly change the number of
    ## colors, e.g. by adding one more color or removing a color.

``` r

# put everything together
ht_list <- ht1 %v% banno %v% ht2 %v% uanno
ComplexHeatmap::draw(ht_list, annotation_legend_list = lgd_list, merge_legend = TRUE)

decorate_annotation("Abundance in Group 1", {
    grid.text("# signatures", x = unit(-1, "cm"), rot = 90, just = "bottom", gp = gpar(fontsize = 8))
    grid.text("*", x = unit(2.45, "cm"), y = unit(1.2, "cm"))
    grid.text("*", x = unit(5.18, "cm"), y = unit(1, "cm"))
    grid.text("*", x = unit(6.55, "cm"), y = unit(0.95, "cm"))
    grid.text("*", x = unit(8.6, "cm"), y = unit(0.85, "cm"))
    grid.text("*", x = unit(10, "cm"), y = unit(0.7, "cm"))
    grid.text("*", x = unit(10.7, "cm"), y = unit(0.7, "cm"))
})
```

![](BugSigDBStats_files/figure-html/mheat-1.png)

### Signature similarity

#### Jaccard index

Inspect signature similarity for signatures from stomach samples based
on Jaccard index:

``` r

stomachsub <- subset(dat, `Body site` == "Stomach")
sigsub <- bugsigdbr::getSignatures(stomachsub)
pair.jsim <- calcJaccardSimilarity(sigsub)
```

Create a dendrogram of Jaccard dissimilarities (1.0 has no overlap, 0.0
are identical signatures).

``` r

jdist <- as.dist(1 - pair.jsim)
plot(hclust(jdist))
```

![](BugSigDBStats_files/figure-html/unnamed-chunk-31-1.png)

# Test replicability of microbiome changes across studies

Function that implements a statistical test for replicability of
microbiome changes for experimental conditions in BugSigDB

## Usage

``` r
testReplicability(
  df,
  onto = NULL,
  multi.contrast = c("all", "first", "largest", "merge"),
  min.studies = 2,
  min.taxa = 5
)
```

## Arguments

- df:

  `data.frame` storing BugSigDB data. Typically obtained via
  [`importBugSigDB`](https://rdrr.io/pkg/bugsigdbr/man/importBugSigDB.html).

- onto:

  ontology. An object of class `ontology_index` storing the NCBI
  Taxonomy. Typically obtained via
  [`getNcbiTaxonomyObo`](http://waldronlab.io/BugSigDBStats/reference/getNcbiTaxonomyObo.md).

- multi.contrast:

  character. How to treat studies that report multiple contrasts for a
  study to avoid detection of duplication within studies as replication?
  Select one of

  - `"all"` incorporates all contrasts reported by a study. This is
    equivalent to do nothing, ie computes similarity between signatures
    taking all, and thus potential similar / duplicated, signatures for
    each study into account.

  - `"first"` take only the first contrast reported by a study, ie the
    first signature with increased abundance and the first signature
    with decreased abundance in the study group.

  - `"largest"` take only the largest signature with increased and
    decreased abundance in the study group, respectively.

  - `"merge"` merge signatures with the same direction of abundance
    change (increased or decreased) within studies.

- min.studies:

  integer. Minimum number of studies for a condition to be tested.
  Defaults to 2, which will then only test a condition investigated by
  at least two studies.

- min.taxa:

  integer. Minimum size for a signature to be included. Defaults to 5,
  which will then only include signatures containing at least 5 taxa.

## Value

A data.frame reporting semantic similarity and re-sampling p-value for
each condition under investigation. Results are stratified by direction
of abundance change in the study group.

## Examples

``` r
 dat <- bugsigdbr::importBugSigDB(version = "10.5281/zenodo.5904281")
#> Using cached version from 2026-07-10 06:24:11
 dat.feces <- subset(dat, `Body site` == "feces")
 res <- testReplicability(dat.feces)
#> Using cached version from 2026-07-10 06:23:35
#> Retrieveing NCBI taxonomy ontology from cache.
```

# Synchronize signatures with NCBI Taxonomy

Helper function for synchronizing signatures with NCBI Taxonomy

## Usage

``` r
syncWithNCBI(sigs, onto)
```

## Arguments

- sigs:

  signatures. A list of character vectors. Typically obtained via
  [`getSignatures`](https://rdrr.io/pkg/bugsigdbr/man/getSignatures.html).

- onto:

  ontology. An object of class `ontology_index` storing the NCBI
  Taxonomy. Typically obtained via
  [`getNcbiTaxonomyObo`](http://waldronlab.io/BugSigDBStats/reference/getNcbiTaxonomyObo.md).

## Value

The input signatures synchronized with the NCBI Taxonomy

## Examples

``` r
 library(bugsigdbr)
#> Note: After Feb. 16, 2025 PubMed ID replaced Study ID in BugSigDb. See https://github.com/waldronlab/BugSigDB/issues/263.
 df <- importBugSigDB()
#> Using cached version from 2026-06-24 06:22:26
 sigs <- getSignatures(df)
 onto <- getNcbiTaxonomyObo()
#> Using cached version from 2026-06-24 06:23:21
#> Retrieveing NCBI taxonomy ontology from cache.
 sigs <- syncWithNCBI(sigs, onto)
 
```

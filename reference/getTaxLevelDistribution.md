# Get distribution of taxonomic levels

Get distribution of taxonomic levels

## Usage

``` r
getTaxLevelDistribution(sigs)
```

## Arguments

- sigs:

  a `list` of BugSigDB signatures in metaphlan format. Typically
  obtained via `getSignatures`.

## Value

A sorted table listing absolute frequencies of unique taxa for each
taxonomic level.

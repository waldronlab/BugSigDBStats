# Taxon stats

Taxon stats

## Usage

``` r
getTaxonStats(
  dat,
  taxon,
  column = "Condition",
  direction = c("BOTH", "UP", "DOWN")
)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data

- taxon:

  character. A taxonomic name.

- column:

  character. A column of `dat` that should be tabled.

- direction:

  character. Indicates direction of abundance change for signatures to
  be included. Use `"UP"` to restrict computation to signatures with
  increased abundance in the exposed group. Use `"DOWN"` to restrict to
  signatures with decreased abundance in the exposed group. Defaults to
  `"BOTH"` which will not filter signatures by direction of abundance
  change.

## Value

A sorted table with taxon stats for the chosen metadata column

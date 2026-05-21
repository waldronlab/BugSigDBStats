# Get most frequent taxa

Get most frequent taxa

## Usage

``` r
getMostFrequentTaxa(dat, n = 10, direction = c("BOTH", "UP", "DOWN"), ...)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data

- n:

  integer. Number of taxa to show. Defaults to `10` which will then show
  the 10 most frequent taxa.

- direction:

  character. Indicates direction of abundance change for signatures to
  be included. Use `"UP"` to restrict computation to signatures with
  increased abundance in the exposed group. Use `"DOWN"` to restrict to
  signatures with decreased abundance in the exposed group. Defaults to
  `"BOTH"` which will not filter signatures by direction of abundance
  change.

- ...:

  Additional arguments passed on to
  [`bugsigdbr::getSignatures`](https://rdrr.io/pkg/bugsigdbr/man/getSignatures.html)

## Value

A sorted table listing absolute frequencies for the most frequent taxa

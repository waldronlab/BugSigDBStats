# Remove entries that do not contain signatures

Remove entries that do not contain signatures

## Usage

``` r
stripEmptySignatures(dat, col = "NCBI Taxonomy IDs")
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data

- col:

  character. A column of `dat` that contain signatures as a comma
  separated list

## Value

A `data.frame` with empty signature entries removed.

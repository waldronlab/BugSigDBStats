# Calculate pairwise Jaccard similarity between all signatures

Calculate pairwise Jaccard similarity between all signatures

## Usage

``` r
calcJaccardSimilarity(sets)
```

## Arguments

- sets:

  a named list of signatures

## Value

a Jaccard similarity matrix

## Examples

``` r
testlist <- list(a = 1:3, b = 3, c = 3:4)
jsim <- calcJaccardSimilarity(testlist)
```

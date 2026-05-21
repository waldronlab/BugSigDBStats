# Plot curation output as a function of time

Plot curation output as a function of time

## Usage

``` r
plotProgressOverTime(dat, col = "Curated date", diff = FALSE)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data.

- col:

  character. A column of `dat` that contain the curation date in dmy
  format.

- diff:

  logical. Display only the difference between months? Defaults to
  `FALSE` which will then display total cumulative numbers.

## Value

None. Plots to a graphics device.

# Plot microbe co-occurrence in a heatmap

Plot microbe co-occurrence in a heatmap

## Usage

``` r
microbeHeatmap(
  dat,
  sig.type = c("both", "increased", "decreased"),
  tax.level = "mixed",
  exact.tax.level = FALSE,
  antagonistic = FALSE,
  anno = "phylum",
  top = 100,
  fontsize = 6,
  ...
)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data.

- sig.type:

  character. Signature type. Use either `"increased"` or `"decreased"`
  to subset to signatures with either increased or decreased abundance
  in the exposed group, respectively. Default is `"both"` which will not
  subset by the direction of abundance change.

- tax.level:

  character. Either `"mixed"` or any subset of
  `c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")`.
  This full vector is equivalent to `"mixed"`.

- exact.tax.level:

  logical. Should only the exact taxonomic level specified by
  `tax.level` be returned? Defaults to `TRUE`. If `FALSE`, a more
  general `tax.level` is extracted for microbes given at a more specific
  taxonomic level.

- antagonistic:

  logical. Antagonistic co-occurrence, ie occurring together but one
  microbe up the other one is down? Defaults to `FALSE`.

- anno:

  character. Taxonomic level that should be displayed as an annotation
  bar. Defaults to `"phylum"`.

- top:

  integer. The number of microbes to display. Defaults to `100`.

- fontsize:

  integer. Fontsize for text. Defaults to `6`.

- ...:

  additional arguments to
  [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

Plots to a graphics device. Returns the co-oocurence matrix invisibly.

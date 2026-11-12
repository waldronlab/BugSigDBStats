# Plot curation output for each curator

Plot curation output for each curator

## Usage

``` r
plotCuratorStats(dat, npc)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data. Used for the plot title (total
  number of papers and signatures).

- npc:

  a numeric `matrix` with rows `"signatures"` and `"papers"` and one
  column per curator, as returned by
  [`stratifyByCurator`](http://waldronlab.io/BugSigDBStats/reference/stratifyByCurator.md).

## Value

None. Plots to a graphics device. Graphical parameters are restored on
exit.

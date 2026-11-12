test_that("stratifyByCurator splits multi-person entries and drops WikiWorks", {
    dat <- data.frame(PMID = c(1, 1, 2, 3, 4),
                      `Revision editor` = c("Ann Lee,WikiWorks743", "Ann Lee, Bob",
                                            "Bob", "WikiWorks", NA),
                      check.names = FALSE)
    npc <- stratifyByCurator(dat, col = "Revision editor")

    expect_equal(rownames(npc), c("signatures", "papers"))
    expect_setequal(colnames(npc), c("Ann", "Bob"))
    expect_equal(npc["signatures", "Ann"], 2)
    expect_equal(npc["papers", "Ann"], 1)
    expect_equal(npc["signatures", "Bob"], 2)
    expect_equal(npc["papers", "Bob"], 2)

    # keep everyone and full names
    all <- stratifyByCurator(dat, col = "Revision editor",
                             exclude = NULL, short.names = FALSE)
    expect_setequal(colnames(all),
                    c("Ann Lee", "Bob", "WikiWorks", "WikiWorks743"))
})

test_that("stratifyByCurator keeps full names when first names collide", {
    dat <- data.frame(PMID = 1:3,
                      Curator = c("Ann Lee", "Ann Smith", "Bob Jones"))
    npc <- stratifyByCurator(dat)
    expect_setequal(colnames(npc), c("Ann Lee", "Ann Smith", "Bob"))
})

.collectorDat <- function()
{
    dat <- data.frame(PMID = c(1, 1, 2, NA),
                      Study = paste("Study", c(1, 1, 2, 3)),
                      `Curated date` = c("10 March 2021", "10 January 2021",
                                         "2021-02-01", ""),
                      check.names = FALSE)
    dat[["MetaPhlAn taxon names"]] <- list(c("a", "b"), "b", c("b", "c"), "d")
    dat
}

test_that("plotCollectorsCurve accumulates unique taxa in chronological order", {
    p <- plotCollectorsCurve(.collectorDat(), by = "signatures", n.perm = 0)
    expect_s3_class(p, "ggplot")
    # order: Jan {b}, Feb {b,c}, Mar {a,b}, undated {d}
    expect_equal(p$data$observed, c(1, 2, 3, 4))
    expect_equal(p$data$effort, 1:4)
})

test_that("plotCollectorsCurve groups by paper and computes permutation band", {
    p <- plotCollectorsCurve(.collectorDat(), by = "papers", n.perm = 20)
    # papers: PMID 1 {a,b}, PMID 2 {b,c}, Study 3 {d}; PMID 1 first by date
    expect_equal(p$data$observed, c(2, 3, 4))
    expect_true(all(p$data$lower <= p$data$expected))
    expect_true(all(p$data$expected <= p$data$upper))
    expect_equal(tail(p$data$expected, 1), 4)
})

test_that("stratifyByCurator warns and returns zero columns when nobody is left", {
    dat <- data.frame(PMID = 1:3,
                      `Revision editor` = c("WikiWorks", "WikiWorks743", NA),
                      check.names = FALSE)
    expect_warning(npc <- stratifyByCurator(dat, col = "Revision editor"),
                   "No people left")
    expect_equal(dim(npc), c(2L, 0L))
    expect_error(plotCuratorStats(dat, npc), "no curators to plot")
})

test_that("plotCuratorStats plots and restores graphical parameters", {
    dat <- data.frame(PMID = c(1, 1, 2),
                      Curator = c("A Very Long Curator Name", "Bob", "Bob"))
    npc <- stratifyByCurator(dat, short.names = FALSE)
    pdf(NULL)
    on.exit(dev.off())
    before <- par(c("las", "mar", "cex"))
    expect_silent(plotCuratorStats(dat, npc))
    expect_equal(par(c("las", "mar", "cex")), before)
})

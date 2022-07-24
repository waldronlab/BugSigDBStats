context("Replicability")

dat <- bugsigdbr::importBugSigDB(version = "10.5281/zenodo.5904281")
dat.feces <- subset(dat, `Body site` == "feces")
dat.saliva <- subset(dat, `Body site` == "saliva")
res.feces <- testReplicability(dat.feces)
res.saliva <- testReplicability(dat.saliva)

checkResult <- function(res)
{
    expect_is(res, "data.frame")
    expect_true(nrow(res) > 1)
    expect_true(ncol(res) == 5)
    rel.cols <- c("CONDITION", "SEMSIM", "PVAL", "NR.STUDIES", "DIRECTION")
    expect_true(all(colnames(res) == rel.cols))
    exp.classes <- c("character", "numeric", "numeric", "integer", "character")
    obs.classes <- sapply(res, class)
    expect_true(all(obs.classes == exp.classes))
}

test_that("valid input df", {
    checkResult(res.feces)
    checkResult(res.saliva) 
    expect_error(testReplicability(dat), "more than one body site")
})

test_that("multi.contrast modes", {
    mc.modes <- c("first", "largest", "merge")
    for(m in mc.modes)
    {
        res <- testReplicability(dat.feces, multi.contrast = m)
        checkResult(res)
    }
})

test_that("min.studies", {
    expect_true(all(res.feces$NR.STUDIES > 1))
    expect_true(all(res.saliva$NR.STUDIES > 1))
    res.feces <- testReplicability(dat.feces, min.studies = 3)
    checkResult(res.feces)
    expect_true(all(res.feces$NR.STUDIES > 2))
})

test_that("min.taxa", {
    res.feces <- testReplicability(dat.feces, min.taxa = 3)
    checkResult(res.feces)
})

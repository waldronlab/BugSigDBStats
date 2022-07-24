context("Replicability")

dat <- bugsigdbr::importBugSigDB(version = "10.5281/zenodo.5904281")
dat.feces <- subset(dat, `Body site` == "feces")
dat.saliva <- subset(dat, `Body site` == "saliva")
onto <- getNcbiTaxonomyObo()
res <- testReplicability(dat, onto = onto)
res.feces <- testReplicability(dat.feces, onto = onto)
res.saliva <- testReplicability(dat.saliva, onto = onto)

checkResult <- function(res, multi.bs = FALSE)
{
    expect_is(res, "data.frame")
    expect_true(nrow(res) > 1)
    expect_true(ncol(res) == ifelse(multi.bs, 6, 5))
    rel.cols <- c("CONDITION", "SEMSIM", "PVAL", "NR.STUDIES", "DIRECTION")
    if(multi.bs) rel.cols <- c(rel.cols[1], "BODY.SITE", rel.cols[2:5])
    expect_true(all(colnames(res) == rel.cols))
    exp.classes <- c("character", "numeric", "numeric", "integer", "character")
    if(multi.bs) exp.classes <- c(exp.classes[1], "character", exp.classes[2:5])
    obs.classes <- sapply(res, class)
    expect_true(all(obs.classes == exp.classes))
}

test_that("valid input df", {
    checkResult(res, multi.bs = TRUE)
    checkResult(res.feces)
    checkResult(res.saliva) 
})

test_that("multi.contrast modes", {
    mc.modes <- c("first", "largest", "merge")
    for(m in mc.modes)
    {
        res <- testReplicability(dat.feces, onto = onto, multi.contrast = m)
        checkResult(res)
    }
})

test_that("min.studies", {
    expect_true(all(res.feces$NR.STUDIES > 1))
    expect_true(all(res.saliva$NR.STUDIES > 1))
    res.feces <- testReplicability(dat.feces, onto = onto, min.studies = 3)
    checkResult(res.feces)
    expect_true(all(res.feces$NR.STUDIES > 2))
})

test_that("min.taxa", {
    res.feces <- testReplicability(dat.feces, onto = onto, min.taxa = 3)
    checkResult(res.feces)
})

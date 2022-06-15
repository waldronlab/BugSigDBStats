############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-01-15 17:59:25
# 
# descr: curation content analysis
# 
############################################################

#' Signature stats
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @return A sorted table with signature stats for the chosen metadata column
#' @export
sigStat <- function(dat, col) sort(table(dat[,col]), decreasing=TRUE)

#' Paper stats
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @return A sorted table with paper stats for the chosen metadata column
#' @export
paperStat <- function(dat, col)
{
    l <- split(dat[,col], dat[,"PMID"])
    ul <- vapply(l, function(x) x[1], dat[1,col])
    sort(table(ul), decreasing=TRUE)
}

#' Taxon stats
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param taxon character. A taxonomic name.
#' @param column character. A column of \code{dat} that should be tabled.
#' @param direction character. Indicates direction of abundance change for signatures
#' to be included. Use \code{"UP"} to restrict
#' computation to signatures with increased abundance in the exposed group. Use 
#' \code{"DOWN"} to restrict to signatures with decreased abundance in the exposed
#' group. Defaults to \code{"BOTH"} which will not filter signatures by direction
#' of abundance change.
#' @return A sorted table with taxon stats for the chosen metadata column
#' @export
getTaxonStats <- function(dat, taxon, 
                          column = "Condition",
                          direction = c("BOTH", "UP", "DOWN"))
{
    direction <- match.arg(direction)
    if(direction %in% c("UP", "DOWN"))
    {
        ind <- !is.na(dat[["Abundance in Group 1"]])
        dat <- dat[ind,]
        direction <- ifelse(direction == "UP",
                           "increased",
                           "decreased")
        dat <- subset(dat, `Abundance in Group 1` == direction)
    }
    tnames <- lapply(dat[["MetaPhlAn taxon names"]], 
                     bugsigdbr::extractTaxLevel,
                     tax.id.type = "taxname")
    ind <- vapply(tnames, 
                  function(x) taxon %in% x,  
                  logical(1))
    tab <- table(dat[ind,column])
    tab <- sort(tab, decreasing = TRUE)
    return(tab)
}

#' Test association of a taxon with a certain category of a column
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param taxon character. A taxonomic name.
#' @param column character. A column of \code{dat} that contains the category 
#' thats hould be tested for association.
#' @param category character. One or more categories that should be tested for
#' association.
#' @param direction character. Indicates direction of abundance change for signatures
#' to be included. Use \code{"UP"} to restrict
#' computation to signatures with increased abundance in the exposed group. Use 
#' \code{"DOWN"} to restrict to signatures with decreased abundance in the exposed
#' group. Defaults to \code{"BOTH"} which will not filter signatures by direction
#' of abundance change.
#' @return A sorted table with taxon stats for the chosen metadata column
#' @importFrom stats prop.test
#' @export
testAssociation <- function(dat, taxon, column = "Condition", 
                            category = "", 
                            direction = c("BOTH", "UP", "DOWN"))
{
    stats <- getTaxonStats(dat, taxon, column, direction)
    stopifnot(all(category %in% names(stats)))

    fg <- sum(stats[category])
    fg.total <- sum(stats)    

    direction <- match.arg(direction)
    if(direction %in% c("UP", "DOWN"))
    {
        ind <- !is.na(dat[["Abundance in Group 1"]])
        dat <- dat[ind,]
        direction <- ifelse(direction == "UP",
                           "increased",
                           "decreased")
        dat <- subset(dat, `Abundance in Group 1` == direction)
    }

    bg <- sum(dat[,column] %in% category)
    bg.total <- nrow(dat)

    prop.test(c(fg, bg), c(fg.total, bg.total))
}

#' Table column
#'
#' @param df a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @param n integer. Number of categories to show. Defaults to \code{Inf}
#' which will then show all categories in the chosen column.
#' @param perc logical. Return absolute frequencies or relative frequencies?
#' Defaults to \code{FALSE} which will return absolute frequencies. 
#' @return A sorted table with table stats for the chosen metadata column
#' @importFrom utils head
#' @export
tabCol <- function(df, col, n = Inf, perc = FALSE)
{
    tab <- table(df[[col]])
    tab <- sort(tab, decreasing = TRUE)
    n <- ifelse(n < Inf, min(n, length(tab)), length(tab)) 
    htab <- head(tab, n)
    if(perc) htab <- signif(htab / sum(tab), digits = 3)
    return(htab)
}

#' Table diversity column
#'
#' @param exps a \code{data.frame} storing BugSigDB experiment data
#' @param div.col character. A column of \code{exps} that should be tabled. 
#' @param spl.col character. A column of \code{exps} that should be used to
#' split \code{div.col} into groups.
#' @param min.exps integer. Minimum number of experiments for a category in
#' \code{spl.col} to be included in the result. Defaults to \code{5}.
#' @param perc logical. Return absolute frequencies or relative frequencies?
#' Defaults to \code{FALSE} which will return absolute frequencies. 
#' @return A sorted table with table stats for the chosen diversity column
#' @export
tabDiv <- function(exps, div.col, spl.col, min.exps = 5, perc = FALSE)
{
    if(div.col == "any") exps <- .summarizeDiv(exps)
    spl <- split(exps[[div.col]], exps[[spl.col]])
    spl <- lapply(spl, function(x) x[!is.na(x)])
    spl <- spl[lengths(spl) >= min.exps]
    tab <- lapply(spl, table)
    tab <- vapply(tab, function(x) x[c("increased", "decreased", "unchanged")], integer(3))
    tab[is.na(tab)] <- 0
    rownames(tab) <- c("increased", "decreased", "unchanged")
    tab <- t(tab)
    abs.diff <- abs(tab[,"increased"] - tab[,"decreased"])
    tab <- tab[order(abs.diff, decreasing = TRUE),]
    if(perc) tab <- t(apply(tab, 1, function(x) signif(x / sum(x), digits = 2)))
    return(tab)
}

#' Get most frequent taxa
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param n integer. Number of taxa to show. Defaults to \code{10}
#' which will then show the 10 most frequent taxa.
#' @param direction character. Indicates direction of abundance change for signatures
#' to be included. Use \code{"UP"} to restrict
#' computation to signatures with increased abundance in the exposed group. Use 
#' \code{"DOWN"} to restrict to signatures with decreased abundance in the exposed
#' group. Defaults to \code{"BOTH"} which will not filter signatures by direction
#' of abundance change.
#' @param ... Additional arguments passed on to \code{bugsigdbr::getSignatures} 
#' @return A sorted table listing absolute frequencies for the most frequent taxa
#' @export
getMostFrequentTaxa <- function(dat,
                                n = 10,
                                direction = c("BOTH", "UP", "DOWN"),
                                ...)
{
    direction <- match.arg(direction)
    if(direction %in% c("UP", "DOWN"))
    {
        ind <- !is.na(dat[["Abundance in Group 1"]])
        dat <- dat[ind,]
        direction <- ifelse(direction == "UP",
                           "increased",
                           "decreased")
        dat <- subset(dat, `Abundance in Group 1` == direction)
    }
    msc <- bugsigdbr::getSignatures(dat, ...)
    msc.tab <- sort(table(unlist(msc)), decreasing=TRUE)
    names(msc.tab) <- vapply(names(msc.tab), bugsigdbr:::.getTip,
                                character(1), USE.NAMES = FALSE)
    head(msc.tab, n=n)
}

.summarizeDiv <- function(exps)
{
    div.cols <- c("Pielou", "Shannon", "Chao1", 
                  "Simpson", "Inverse Simpson", "Richness")
    .hasConflict <- function(x) all(c("increased", "decreased") %in% x)
    ind <- apply(exps[,div.cols], 1, .hasConflict)
    exps <- exps[!ind,]
    .decideDiv <- function(x)
    {
        y <- "unchanged"
        if("increased" %in% x) y <- "increased"
        else if("decreased" %in% x) y <- "decreased"  
        return(y)
    }
    any.col <- apply(exps[,div.cols], 1, .decideDiv)
    exps <- cbind(exps, any = any.col)
    return(exps)
} 

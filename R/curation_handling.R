############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-01-15 17:59:18
# 
# descr: curation data handling
# 
############################################################

TAX.LEVELS <- c("kingdom", "phylum", "class", "order",
                "family", "genus", "species", "strain")
MPA.TAX.LEVELS <- c(substring(TAX.LEVELS[1:7], 1, 1), "t")
names(MPA.TAX.LEVELS) <- TAX.LEVELS

MPA.REGEXP <- "^[kpcofgst]__"

#' Get the publication year for a PMID
#'
#' @param pmid character. A PubMed ID. 
#' @return A character containg the publication year for the given \code{pmid}. 
#' @examples 
#'  pmid2pubyear("32026945")
#' @export
pmid2pubyear <- function(pmid)
{
    x <- annotate::pubmed(pmid)
    a <- XML::xmlRoot(x)
    .f <- function(i)
    {
        art <- annotate::buildPubMedAbst(a[[i]])
        pd <- annotate::pubDate(art)
        unlist(strsplit(pd, " "))[[2]]
    }
    vapply(seq_along(pmid), .f, character(1))
}

#' Remove entries that do not contain signatures
#
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that contain signatures
#' as a comma separated list 
#' @return A \code{data.frame} with empty signature entries removed. 
#' @export
stripEmptySignatures <- function(dat, col = "NCBI Taxonomy IDs")
{
    dat[is.na(dat[,col]), col] <- ""
    non.empty.sigs <- dat[,col] != ""
    dat[non.empty.sigs,]
}


#' Stratify curation output by curator 
#' @param dat a \code{data.frame} storing BugSigDB data
#' @return a numeric \code{matrix} storing number of papers and number of 
#' signatures curated for each curator 
#' @export
stratifyByCurator <- function(dat)
{
    pc <- split(dat[,"PMID"], dat[,"Curator"])
    spc <- sort(lengths(pc))
    ppc <- sort(lengths(lapply(pc, unique)))

    npc <- rbind(spc[names(ppc)], ppc)
    colnames(npc) <- vapply(colnames(npc), function(n) unlist(strsplit(n, " "))[1], character(1))
    rownames(npc) <- c("signatures", "papers")
    npc[,!is.na(colnames(npc))]
}


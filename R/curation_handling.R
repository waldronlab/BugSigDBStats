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
pmid2pubyear <- function(pmid) {
  chunk_size <- 250
  pmid_chunks <- split(pmid, ceiling(seq_along(pmid) / chunk_size))
  pub_years <- character(0)
  
  for (chunk in pmid_chunks) {
    x <- annotate::pubmed(chunk)
    a <- XML::xmlRoot(x)
    .f <- function(i) {
      art <- annotate::buildPubMedAbst(a[[i]])
      pd <- annotate::pubDate(art)
      unlist(strsplit(pd, " "))[[2]]
    }
    
    pub_years_chunk <- vapply(seq_along(chunk), .f, character(1))
    pub_years <- c(pub_years, pub_years_chunk)
    Sys.sleep(2)
  }
  
  return(pub_years)
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
#'
#' Counts the number of signatures and papers each person contributed to.
#' Columns such as \code{"Revision editor"} can list several people per
#' signature (e.g. \code{"Claregrieve1,Chloe,WikiWorks743"}); these entries
#' are split so that each person is credited individually, regardless of
#' their position in the list.
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. The column of \code{dat} listing the people to
#' count. Defaults to \code{"Curator"}; use \code{"Revision editor"} to count
#' revision contributions.
#' @param sep character. Separator between multiple people in \code{col}.
#' @param exclude character. A regular expression; people whose name matches
#' it are dropped. Defaults to dropping the \code{WikiWorks} maintenance
#' accounts. Use \code{NULL} to keep everyone.
#' @param short.names logical. Shorten names to the first word for display
#' (e.g. \code{"Aleru Divine"} becomes \code{"Aleru"})? Names whose first word
#' is shared by more than one person are kept in full.
#' @return a numeric \code{matrix} with rows \code{"signatures"} and
#' \code{"papers"} and one column per person, ordered by increasing number
#' of papers. If no people remain after removing empty and excluded entries,
#' a warning is issued and a matrix with zero columns is returned.
#' @examples
#' dat <- data.frame(PMID = c(1, 1, 2, 3),
#'                   `Revision editor` = c("Ann,WikiWorks743", "Ann,Bob",
#'                                         "Bob", "WikiWorks"),
#'                   check.names = FALSE)
#' stratifyByCurator(dat, col = "Revision editor")
#' @export
stratifyByCurator <- function(dat, col = "Curator", sep = ",",
                              exclude = "^WikiWorks", short.names = TRUE)
{
  stopifnot(all(c("PMID", col) %in% colnames(dat)))

  # one row per (signature, person)
  people <- as.character(dat[[col]])
  people[is.na(people)] <- ""
  people <- lapply(strsplit(people, sep, fixed = TRUE),
                   function(p) unique(trimws(p)))
  long <- data.frame(PMID = rep(dat[["PMID"]], lengths(people)),
                     person = unlist(people, use.names = FALSE))

  keep <- !is.na(long$person) & long$person != ""
  if(length(exclude)) keep <- keep & !grepl(exclude, long$person)
  long <- long[keep,]
  if(!nrow(long))
    warning("No people left in column '", col,
            "' after removing empty and excluded entries")

  pc <- split(long$PMID, long$person)
  ppc <- sort(lengths(lapply(pc, unique)))
  spc <- lengths(pc)[names(ppc)]

  npc <- rbind(signatures = spc, papers = ppc)
  if(short.names) colnames(npc) <- .shortenNames(colnames(npc))
  npc
}

# first word of each name, unless that would merge different people
.shortenNames <- function(x)
{
  first <- sub(" .*", "", x)
  ambiguous <- first %in% first[duplicated(first)]
  first[ambiguous] <- x[ambiguous]
  first
}
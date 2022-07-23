#' Calculate semantic similarity
#'
#' \code{calcSemanticSimilarity} calculates the pairwise semantic similarity
#' between two sets of signature files exported from BugSigDB.org.
#'
#' @param file1 Signature exported from BugSigDB.org.
#' @param file2 Signature exported from BugSigDB.org.
#' @export
#'
#' @examples 
#'
#' \dontrun{
#'
#' library(tidyr)
#' library(tibble)
#'
#' ## Download files froom BugSigDB's drilldown ##
#'
#' file1_url <- "https://bugsigdb.org/w/index.php?title=Special:Ask&x=-5B-5BCategory%3ASignatures-5D-5D-20-5B-5BModification-20date%3A%3A%2B-5D-5D-20-5B-5BBase-20page.Location-20of-20subjects%3A%3AAustralia-5D-5D%2F-3FOriginal-20page-20name%3DSignature-20page-20name%2F-3FRelated-20experiment%3DExperiment%2F-3FRelated-20study%3DStudy%2F-3FSource-20data%3DSource%2F-3FCurated-20date%2F-3FCurator%2F-3FRevision-20editor%2F-3FDescription%2F-3FAbundance-20in-20Group-201%2F-3FNCBI-20export%3DMetaPhlAn-20taxon-20names%2F-3FNCBI-20export-20ids-20sc%3DNCBI-20Taxonomy-20IDs%2F-3FState%2F-3FReviewer&mainlabel=-&limit=5000&offset=0&format=csv&searchlabel=%3Cdiv%20class%3D%22mw-ui-button%20mw-ui-quiet%20mw-ui-progressive%20rounded-0%22%3ESignatures%3C%2Fdiv%3E&filename=signatures-filtered-australia.csv"
#' file1 <- tempfile()
#' download.file(url = file1_url, destfile = file1)
#' 
#' file2_url <- "https://bugsigdb.org/w/index.php?title=Special:Ask&x=-5B-5BCategory%3ASignatures-5D-5D-20-5B-5BModification-20date%3A%3A%2B-5D-5D-20-5B-5BBase-20page.Location-20of-20subjects%3A%3AUnited-20States-20of-20America-5D-5D%2F-3FOriginal-20page-20name%3DSignature-20page-20name%2F-3FRelated-20experiment%3DExperiment%2F-3FRelated-20study%3DStudy%2F-3FSource-20data%3DSource%2F-3FCurated-20date%2F-3FCurator%2F-3FRevision-20editor%2F-3FDescription%2F-3FAbundance-20in-20Group-201%2F-3FNCBI-20export%3DMetaPhlAn-20taxon-20names%2F-3FNCBI-20export-20ids-20sc%3DNCBI-20Taxonomy-20IDs%2F-3FState%2F-3FReviewer&mainlabel=-&limit=5000&offset=0&format=csv&searchlabel=%3Cdiv%20class%3D%22mw-ui-button%20mw-ui-quiet%20mw-ui-progressive%20rounded-0%22%3ESignatures%3C%2Fdiv%3E&filename=signatures-filtered-united-states-of-america.csv"
#' file2 <- tempfile()
#' download.file(url = file2_url, destfile = file2)
#'
#' ## Calculate semantic similarity between pairs of signatures ##
#' 
#' output <- calcSemanticSimilarity(file1, file2)
#' head(output)
#'
#' ## Conver output to a matrix ##
#' 
#' mat <- pivot_wider(
#'     output, names_from = sig2, values_from = semantic_similarity
#'     ) |> 
#'     tibble::column_to_rownames(var = "sig1") |> 
#'     as.data.frame() |> 
#'     as.matrix()
#' 
#' dim(mat)
#' mat[1:5, 1:5]
#'
#' }
#' 
calcSemanticSimilarity <- function(file1, file2) {
    
    on.exit(gc())
    
    sig1 <- .bugSigDBExportToSig(file1)
    sig2 <- .bugSigDBExportToSig(file2)
    ncbi_onto <- getNcbiTaxonomyObo() 
    
    utax <- unique(unlist(sig1))
    nt <- utax[!(utax %in% ncbi_onto$id)]
    sig1 <- lapply(sig1, function(s) setdiff(s, nt))
    
    utax <- unique(unlist(sig2))
    nt <- utax[!(utax %in% ncbi_onto$id)]
    sig2 <- lapply(sig2, function(s) setdiff(s, nt))
    
    mat <- ontologySimilarity::get_sim_grid(
        ontology = ncbi_onto, term_sets = sig1, term_sets2 = sig2
    )
    
    output <- as.data.frame(as.table(mat))
    names(output) <- c("sig1", "sig2", "semantic_similarity")
    output
}

#' Get NCBI Taxonomy in OBO format
#'
#' Gets the NCBI taxonomy in OBO format.
#'
#' @return NCBI taxonomy in OBO format
#' @importFrom utils download.file
#' @export
#' 
getNcbiTaxonomyObo <- function() {

    onto <- .getResourceFromCache("ncbi.onto")

    if (is.null(onto)) {

        url <- paste0("https://sandbox.zenodo.org/record/954943/files/",
                      "ncbitaxon.rds?download=1")
        temp_file <- tempfile()
        message("Getting NCBI taxonomy ontology.")
        utils::download.file(url = url, destfile = temp_file)
        onto <- readRDS(temp_file)
        .cacheResource(onto, "ncbi.onto")
        message("NCBI taxonomy ontology has been saved in cache.")
        return(onto)

    } else {
        
        message("Retrieveing NCBI taxonomy ontology from cache.")
        return(onto)

    }
}

.bugSigDBExportToSig <- function(fname) {
    df <- utils::read.csv(
        fname, header = TRUE, comment.char = "#", check.names = FALSE
    )
    
    id_cols <- c("Signature page name", "Experiment", "Study")
    ids <- lapply(df[, id_cols], function(x) sub("^.+ ", "", x))
    sig_names <- vector("character", length(ids[[1]]))
    for (i in seq_along(sig_names)) {
        sig_names[i] <- 
            paste0("bsdb:", ids[[1]][i], "/", ids[[2]][i], "/", ids[[3]][i])
    }
    sig <- lapply(df[["NCBI Taxonomy IDs"]], function(x) {
        sig <- strsplit(x, ";")[[1]] |>
            vapply(\(y) sub("^.+\\|", "", y), character(1))
        names(sig) <- NULL
        paste0("NCBITaxon:", sig)
    })
    names(sig) <- sig_names
    sig
}

#' Weighted semantic similarity
#'
#' Incorporation of weights into the computation of semantic similarity
#'
#' @param o ontology. An object of class \code{ontologyIndex}.
#' @param ic information content. Typically obtained via 
#' \code{ontologySimilarity::descendants_IC(o)}. 
#' @param set1 character. First term set.
#' @param set2 character. Secound term set.
#' @param weights1 numeric. Weights for each term of the first term set, ie. must
#' be parallel to \code{set1}. 
#' @param weights2 numeric. Weights for each term of the second term set, ie. must
#' be parallel to \code{set2}. 
#' @return a numeric value expressing semantic similarity between the two
#' input term sets, weighted by the given individual term weights. 
#' @export
weightedBMA <- function(o, ic, set1, set2, 
                        weights1 = rep(1, length(set1)), 
                        weights2 = rep(1, length(set2))) 
{
    m <- ontologySimilarity::get_term_sim_mat(ontology = o,
                                              information_content = ic,
                                              method = "lin",
                                              row_terms = set1,
                                              col_terms = set2
    )
    rmax <- matrixStats::rowMaxs(m, na.rm = TRUE)
    cmax <- matrixStats::colMaxs(m, na.rm = TRUE)
    mean1 <- sum(weights1 * rmax) / sum(weights1)
    mean2 <- sum(weights2 * cmax) / sum(weights2)
    m <- mean(c(mean1, mean2))
    return(m)
}

#' Synchronize signatures with NCBI Taxonomy
#'
#' Helper function for synchronizing signatures with NCBI Taxonomy
#'
#' @param sigs signatures. A list of character vectors. Typically obtained 
#' via \code{\link{getSignatures}}. 
#' @param onto ontology. An object of class \code{ontologyIndex} storing the 
#' NCBI Taxonomy. Typically obtained via \code{\link{getNcbiTaxonomyObo}}.
#' @return The input signatures synchronized with the NCBI Taxonomy
#' @examples
#'  library(bugsigdbr)
#'  df <- importBugSigDB()
#'  sigs <- getSignatures(df)
#'  onto <- getNcbiTaxonomyObo()
#'  sigs <- syncWithNCBI(sigs, onto)
#'  
#' @export
syncWithNCBI <- function(sigs, onto)
{
    sigs <- lapply(sigs, function(s) paste0("NCBITaxon:", s))
    utax <- unique(unlist(sigs))
    nt <- utax[!(utax %in% onto$id)]
    sigs <- lapply(sigs, function(s) setdiff(s, nt))
    return(sigs)
}

#' Test replicability of microbiome changes across studies
#'
#' Function that implements a statistical test for replicability of microbiome
#' changes for experimental conditions in BugSigDB
#'
#' @param df \code{data.frame} storing BugSigDB data. Typically obtained via
#' \code{\link{importBugSigDB}}.
#' @param multi.contrast character. How to treat studies that report multiple
#' contrasts for a study to avoid detection of duplication within studies as
#' replication? Select one of \itemize{
#' \item \code{"all"} incorporates all contrasts reported by a study. This is 
#' equivalent to do nothing, ie computes similarity between signatures taking all,
#' and thus potential similar / duplicated, signatures for each study into account.  
#' \item \code{"first"} take only the first contrast reported by a study, ie the
#' first signature with increased abundance and the first signature with decreased
#' abundance in the study group.
#' \item \code{"largest"} take only the largest signature with increased and
#' decreased abundance in the study group, respectively. 
#' \item \code{"merge"} merge signatures with the same direction of abundance
#' change (increased or decreased) within studies. 
#' }
#' @param min.studies integer. Minimum number of studies for a condition to be
#' tested. Defaults to 2, which will then only test a condition investigated by
#' at least two studies. 
#' @param min.taxa integer. Minimum size for a signature to be included. Defaults to 5, which
#' will then only include signatures containing at least 5 taxa.
#' @return A data.frame reporting semantic similarity and re-sampling p-value for
#' each condition under investigation. Results are stratified by direction of abundance
#' change in the study group.
#' @examples
#'  dat <- bugsigdbr::importBugSigDB(version = "10.5281/zenodo.5904281")
#'  dat.feces <- subset(dat, `Body site` == "feces")
#'  res <- testReplicability(dat.feces)
#' @export
testReplicability <- function(df, 
                              multi.contrast = c("all", "first", "largest", "merge"),
                              min.studies = 2,
                              min.taxa = 5)
{
    # sanity check on column names
    stopifnot("Condition" %in% colnames(df))
    stopifnot("Body site" %in% colnames(df))
    stopifnot("NCBI Taxonomy IDs" %in% colnames(df))

    if(length(unique(df[["Body site"]])) > 1) 
        stop("Found more than one body site.\n", 
             "Replicability testing is supported for one body site at a time only.\n",
             "Subset by body site of interest first.")

    # clean pubmed column
    ind <- !is.na(df$PMID)
    df <- df[ind,]

    # clean condition column
    column <- "Condition"
    ind <- !is.na(df[[column]]) & !grepl(",", df[[column]])
    df <- df[ind,]

    # treat multi contrast studies
    df.up <- .restrictByDirection(df, direction = "UP")
    df.down <- .restrictByDirection(df, direction = "DOWN")
    df.up <- .treatMultiContrastStudies(df.up, multi.contrast)
    df.down <- .treatMultiContrastStudies(df.down, multi.contrast)
    df <- rbind(df.up, df.down)

    # include only signatures with defined min number of taxa
    ind <- lengths(df[["NCBI Taxonomy IDs"]]) >= min.taxa
    df <- df[ind,]

    # calculate semantic similarity
    onto <- getNcbiTaxonomyObo()
    sigs <- bugsigdbr::getSignatures(df, tax.id.type = "ncbi")
    sigs <- syncWithNCBI(sigs, onto) 
    sim.mat <- ontologySimilarity::get_sim_grid(ontology = onto, term_sets = sigs)

    # include only categories with defined min number of studies
    conds.to.test <- .getCategoriesToTest(df, column, min.studies)
    
    # test conditions
    ps.up <- .testConditions(names(conds.to.test), df, sim.mat, "increased")
    ps.up$NR.STUDIES <- unname(conds.to.test[rownames(ps.up)])
    ps.down <- .testConditions(names(conds.to.test), df, sim.mat, "decreased") 
    ps.down$NR.STUDIES <- unname(conds.to.test[rownames(ps.down)])

    # combine into result table
    res <- .createResultTable(ps.up, ps.down) 
    return(res)
}

# take the results by direction of abundance change (increased / decreased)
# and combine them into on result table
.createResultTable <- function(ps.up, ps.down)
{
    isect <- intersect(rownames(ps.up), rownames(ps.down))
    res <- cbind(ps.up[isect,], ps.down[isect,])
    uranks <- .getRanks(ps.up)
    dranks <- .getRanks(ps.down)
    mranks <- (uranks[isect] + dranks[isect]) / 2 
    ind <- order(mranks[rownames(res)])
    res <- res[ind,]   
    colnames(res)[c(2:4,6:8)] <- paste(colnames(res)[c(2:4,6:8)], 
                                      rep(c("UP", "DOWN"), each = 2), 
                                      sep = ".") 
    res <- data.frame(CONDITION = rep(rownames(res), 2),
                      SEMSIM =  c(res$SEMSIM.UP, res$SEMSIM.DOWN), 
                      PVAL = c(res$PVAL.UP, res$PVAL.DOWN), 
                      NR.STUDIES = c(res$NR.STUDIES.UP, res$NR.STUDIES.DOWN),   
                      DIRECTION = rep(c("UP", "DOWN"), each = nrow(res)))
    return(res)
} 

# get ranks for a result data frame eg sorted by p-value
.getRanks <- function (res, rank.fun = c("comp.ranks", "rel.ranks", "abs.ranks"), 
                       rank.col = "PVAL", name.col = "CONDITION", decreasing = FALSE) 
{
    if (is.function(rank.fun)) 
        ranks <- rank.fun(res)
    else {
        rank.fun <- match.arg(rank.fun)
        rcol <- res[, rank.col]
        if (decreasing) 
            rcol <- -rcol
        names(rcol) <- res[, name.col]
        if (rank.fun == "comp.ranks") 
            ranks <- vapply(rcol, function(p) mean(rcol <= p) * 
                100, numeric(1))
        else {
            ucats <- unique(rcol)
            ranks <- match(rcol, ucats)
            if (rank.fun == "rel.ranks") 
                ranks <- ranks/length(ucats) * 100
        }
        return(ranks)
    }
    return(ranks)
}


# tests a character vector of conditions for semantic similarity of the
# signatures for one condition at a time
.testConditions <- function(conditions, df, sim.mat,
                            direction = c("increased", "decreased"))
{
    direction <- match.arg(direction)
    ps <- vapply(conditions, .testCondition, numeric(2), 
                 df = df, sim.mat = sim.mat, direction = direction)
    ps <- t(ps)
    ps <- ps[!is.na(ps[,"p"]),]
    ps <- ps[order(ps[,"p"]),]
    res <- data.frame(CONDITION = rownames(ps), 
                     SEMSIM = ps[,"sim"],
                     PVAL = ps[,"p"])
    return(res)
}

# test signatures of a specified condition for semantic similarity
.testCondition <- function(condition, df, sim.mat, direction)
{
    cond <- df$Condition
    dir <- df[["Abundance in Group 1"]]
    ind <- which(cond == condition & dir == direction)
    if(length(ind) < 2) res <- c(sim = NA, p = NA)
    else res <- c(sim = ontologySimilarity::get_sim(sim.mat, group = ind),
                  p = ontologySimilarity::get_sim_p(sim.mat, group = ind))
    return(res)
}

# solutions for dealing with multiple contrasts for a study to avoid detection
# of duplication within studies as replication
# l <- rle(df$PMID)
# cs <- cumsum(l$lengths)
# ind <- c(1, cs[-length(cs)] + 1)
.first <- function(s) s[1,]

.largest <- function(s)
{
    l <- lengths(s[["NCBI Taxonomy IDs"]])
    s[which.max(l),] 
}

.merge <- function(s)
{
    s[["NCBI Taxonomy IDs"]][[1]] <- Reduce(union, s[["NCBI Taxonomy IDs"]])
    s[1,"Group 1 name"] <- paste(s[,"Group 1 name"], collapse = ";") 
    s[1,"Group 0 name"] <- paste(s[,"Group 0 name"], collapse = ";") 
    .first(s)
}

.treatMultiContrastStudies <- function(df, 
                                       multi.contrast = c("all", "first",
                                                          "largest", "merge"))
{
    multi.contrast <- match.arg(multi.contrast)
    if(multi.contrast == "all") return(df)
    spl <- split(df, df$PMID)
    .f <- switch(multi.contrast,
                 first = .first,
                 largest = .largest,
                 merge = .merge) 
    spl <- lapply(spl, .f)
    df <- do.call(rbind, spl)
    return(df)
}

# obtain categories of a column with defined minimum number of studies
# investigating this categories
# can eg. be used to obtain the conditions that are studies by at least 
# two studies for replicability testing
.getCategoriesToTest <- function(df, column, min.studies)
{
    spl <- split(df$PMID, df[[column]])
    spl <- lapply(spl, unique)
    lens <- lengths(spl)
    names(lens) <- names(spl)
    lens <- sort(lens, decreasing = TRUE)
    incl <- lens[lens >= min.studies]
    return(incl)
}

# restrict a BugSigDB data frame by direction of abundance change
# UP: increased abundance in the study group
# DOWN: decreased abundance in the study group
.restrictByDirection <- function(df, direction = c("BOTH", "UP", "DOWN"))
{
    direction <- match.arg(direction)
    if(direction != "BOTH")
    {
        direction <- ifelse(direction == "UP", "increased", "decreased")
        df <- subset(df, `Abundance in Group 1` == direction)
    }
    return(df)
}

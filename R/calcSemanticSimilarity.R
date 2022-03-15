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
#' @param method Measure for semantic similarity. Default's to \code{"lin"}
#' which will then use Lin's measure.
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


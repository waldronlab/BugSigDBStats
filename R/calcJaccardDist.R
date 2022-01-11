#' Calculate pairwise Jaccard/ Overlap distance between all signatures
#'
#' @param sets
#' a named list of signatures
#' @param use
#' Either "overlap" or "jaccard"
#'
#' @return a Dist object of (1 - jaccard) or (1 - overlap)
#' @export
#' @importFrom stats as.dist
#'
#' @examples
#' testlist <- list(a = 1:3, b = 3, c = 3:4)
#' jdist <- calcJaccardDist(testlist, "jaccard")
#' 
#' @details 
#' See https://blog.jdblischak.com/posts/pairwise-overlaps/ page for the definition
#' of overlap and jaccard distance by John Blischak
calcJaccardDist <- function(sets, use = "jaccard") {
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x)
    length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }
  
  if(!identical(use %in% c("jaccard", "overlap"), TRUE))
    stop("overlap must be either 'jaccard' or 'overlap'")
  
  n_sets <- length(sets)
  iseq <- seq_len(n_sets - 1)
  
  set_names <- names(sets)
  overlaps_index <- 1
  
  vec_name1 <- character()
  vec_name2 <- character()
  vec_overlap <- numeric()
  vec_jaccard <- numeric()
  
  jsim <- matrix(0.00, length(set_names), length(set_names), dimnames=list(set_names, set_names))

  for (i in iseq) {
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      set2 <- sets[[j]]
      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <-
        unique(
          c(set1, set2),
          incomparables = FALSE,
          fromLast = FALSE,
          nmax = NA
        )
      length_intersection <- length(set_intersect)
      if (length_intersection > 0.00){
        length_union <- length(set_union)
        length_set1 <- length(set1)
        length_set2 <- length(set2)
        overlap <- (length_intersection / min(length_set1, length_set2))
        jaccard <- (length_intersection / length(set_union))
        name1 <- set_names[i]
        name2 <- set_names[j]
        
        vec_name1[overlaps_index] <- name1
        vec_name2[overlaps_index] <- name2
        if(use == "jaccard"){
          jsim[name1, name2] <- jaccard
          jsim[name2, name1] <- jaccard
        }
        if(use == "overlap"){
          jsim[name1, name2] <- overlap
          jsim[name2, name1] <- overlap
        }
        
        overlaps_index <- overlaps_index + 1
      }
    }
  }
  
  jdiff <- as.dist(1 - jsim)
  return(jdiff)
}

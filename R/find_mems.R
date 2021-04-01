#' @title Find MEMs that maximize the variance explained on a set of variables
#' @description Moran's Eigenvetor Maps are a synthetic orthogonal basis that
#'   decompose the spatial structure among a set of sites according to a spatial
#'   weighting matrix (SWM). This functions find the MEMs that explains the
#'   highest proportion of variability in a set of variables.
#' @param df dataframe, with coordinates and variables for which MEMs are
#'   calculated.
#' @param lon_name string, the column name of longitude coordinates in
#'   \code{df}.
#' @param lat_name string, the column name of latitude coordinates in \code{df}.
#' @param scannf boolean, if to display the PCA plot to choose the number of
#'   principal components.
#' @param nf integer, number of principal components, only if \code{scannf =
#'   FALSE}.
#' @param weights string vector, the weights of the SWM. Possible values are:
#'   - "bin": binary
#'   - "flin": linear
#'   - "fdown": concave down
#'   - "fup": concave up
#' @param nb string vector, how the SMW is build. Possible values are:
#'   - "del": Delaunay triangulation
#'   - "gab": Gabriel's graph
#'   - "rel": relative neighbourgood graph
#'   - "mst": minimum spanning tree
#'   - "pcnm": distance-based SWM based on PCNM
#'   - "dnear": distance-based
#'   
#' @param d1 integer, minimum distance for the "dnear" method
#' @param d2 integer, maximum distance for the "dnear" method
find_mems <- function(
  df,
  lon_name = "lon",
  lat_name = "lat",
  scannf = TRUE,
  nf = NULL,
  weights = "fup",
  nb = c("pcnm", "gab", "dnear"),
  d1 = NULL,
  d2 = NULL
) {
  "%out%" <- Negate("%in%")
  if ("data.frame" %out% class(df)) {
    stop("df must be a data.frame or a tibble")
  }
  if ("tbl" %in% class(df)) {
    df <- tibble::as_data_frame(df)
  }
  # split coordinates and variables and remove non-numeric ones.
  coords <- df[, colnames(df) %in% c("lon", "lat")]
  variables <- df[, colnames(df) %out% c("lon", "lat")]
  type <- sapply(colnames(variables), function(x) class(df[, x]))
  if (any(type != "numeric")) {
    message("Removing non-numeric variables from the data.frame")
    variables <- variables[, which(type == "numeric")]
  }
  if (!scannf & is.null(nf)) {
    stop("Specify number of component to keep in the PCA")
  }
  # PCA on variables
  pca <- ade4::dudi.pca(variables, scale = TRUE, scannf = scannf, nf = nf)
  # create all possible combination of selected methods to compute MEMs
  if (is.null(d1) & is.null(d2)) {
    cand_lw <- adespatial::listw.candidates(coords,
                                            nb = nb,
                                            weights = weights)
  } else if (is.null(d1)) {
    cand_lw <- adespatial::listw.candidates(coords,
                                            nb = nb,
                                            weights = weights,
                                            d1 = d1)
  } else if (is.null(d2)) {
    cand_lw <- adespatial::listw.candidates(coords,
                                            nb = nb,
                                            weights = weights,
                                            d2 = d2)
  } else {
    cand_lw <- adespatial::listw.candidates(coords,
                                            nb = nb,
                                            weights = weights,
                                            d1 = d1,
                                            d2 = d2)
  }
  # select most explicative method
  sel_lw <- listw.select(pca$tab,
                         candidates = cand_lw,
                         nperm = 99)
  lw_best <- cand_lw[[sel_lw$best.id]]
  ans <- list(
    best.method.name = sel_lw$best.id,
    best.method = sel_lw$best,
    lw.best = lw_best
  )
  return(ans)
}

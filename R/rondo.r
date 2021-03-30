#' rondo: Analytical calculation of groundwater flow in a single aquifer under circular shaped polders.
#'
#' Based on the publication "Wegzijging en kwel; de grondwaterstroming van hogere naar lagere gebieden".
#' L.F. Ernst. I.C.W.-rapport nr. 7 (1983).
#' \href{https://edepot.wur.nl/192343}{Rapport "Wegzijging en kwel" L.F. Ernst (1983).}
#'
#' This package exports the following functions:
#'
#' * \code{\link{rd_init}}: Initialise a list (y) with parameters used for the calculations of the groundwater heads and fluxes.
#'
#' * \code{\link{rd_phi}}: Calculate the head at radius x.
#'
#' * \code{\link{rd_q}}: Calculate the lateral discharge at radius x.
#'
#' * \code{\link{rd_seep}}: Calculate the seepage intensity at radius x.
#'
#' * \code{\link{av_seep}}: Calculate the average seepage intensity between radius r1 and r2.
#'
#' @docType package
#' @name rondo
#'
#' @importFrom magrittr %>%
#'
#' @importFrom stats approx
#'
NULL

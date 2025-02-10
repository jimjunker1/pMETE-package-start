#' Fish community data
#' 
#' An example data set to estimate the metabolic scaling exponent \eqn{\beta} at the community, \eqn{\beta_c}, and population, \eqn{\beta_p}, levels.
#'
#' @format ## `size_data_com`
#' A data frame with 366 rows and 7 columns:
#' \describe{
#'  \item{ope_id}{sampling event ID}
#'  \item{species_name}{species scientific name}
#'  \item{mass}{individual mass in grams, determined from length}
#'  \item{length}{individual length in mm}
#'  \item{comments_length}{notes for length measurements}
#'  \item{comments_mass}{notes for mass measurements}
#'  \item{outlier}{integer 0/1 if observation is an outlier}
#'  }
#'   
#'
"size_data_com"
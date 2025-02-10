#' "eqn_beta" estimate $beta$ for the lagrange multipliers
#' 
#' @param x 
#' @param N_0 
#' @param S_0 
eqn_beta <- function(x,N_0,S_0){
  sum(exp(-x*seq(1,N_0)))/sum((exp(-x*seq(1,N_0)))/seq(1,N_0))-N_0/S_0
}
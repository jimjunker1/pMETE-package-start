#' negLL_com_pmete
#' @description 
#' Deriving the log-likelihood function for the individual size distribution to solve \eqn{\beta_{c}} and \eqn{\lambda_1} and \eqn{\lambda_2} 
#' @param par 
#' @param x 
#' @param xmin 
#' @param beta 
#' @param N_0 
#' @param S_0 
negLL_com_pmete <- function(par,x,xmin,beta,N_0,S_0){
  bc = par[1]  # community-level metabolic scaling exponent
  L2 = S_0/(sum((x/xmin)^bc)-N_0) # lambda_2 by eqn 7.26 in Harte (2011)
  L1 = beta-L2  # lambda_1
  negLL_exact = 
    -N_0*log(bc*L2)+(1-bc)*sum(log(x/xmin))+sum(L1+L2*(x/xmin)^bc)+sum(log((1-exp(-L1-L2*(x/xmin)^bc))^2))-sum(log(1-(N_0+1)*exp(-N_0*(L1+L2*(x/xmin)^bc))+N_0*exp(-(N_0+1)*(L1+L2*(x/xmin)^bc))))+N_0*log((exp(-beta)-exp(-beta*(N_0+1)))/(1-exp(-beta))-(exp(-(L1+L2*sum((x/xmin)^bc)))-exp(-(L1+L2*sum((x/xmin)^bc))*(N_0+1)))/(1-exp(-(L1+L2*sum((x/xmin)^bc)))))
  return(negLL_exact)
}

# left-truncated exponential
#' Title
#'
#' @param par 
#' @param x 
#' @param xmin 
#' @param n 
#'
#' @export
#'
negLL_exp = function(par,x,xmin,n){
  lambda = par[1]  # metabolic scaling exponent
  negLL = -n*log(lambda)-n*lambda+lambda*sum(x/xmin)
  return(negLL)
}

# left-truncated power-law
#' Title
#'
#' @param par 
#' @param x 
#' @param xmin 
#' @param n 
#'
#' @export
#'
negLL_pow = function(par,x,xmin,n){
  a = par[1]  # power-law distribution parameter
  negLL = -n*log(a-1)+a*sum(log(x/xmin))
  return(negLL)
}

# left-truncated Weibull
#' Title
#'
#' @param par 
#' @param x 
#' @param n 
#' @param xmin 
#'
#' @export
#'
negLL_wei = function(par,x,n,xmin){
  a = par[1]  # scale parameter
  b = par[2]  # shape parameter
  negll = -n*log(b)+n*b*log(a)+(1-b)*sum(log(x/xmin))+sum(((x/xmin)/a)^b)-n*(1/a)^b
  return(negll)
}

# left-truncated quasi-Weibull
#' Title
#'
#' @param par 
#' @param x 
#' @param xmin 
#' @param n 
#'
#' @export
#'
negLL_qwei = function(par,x,xmin,n){
  qal = par[1]  # quasi-Weibull distribution parameter alpha
  qbe = par[2]  # quasi-Weibull distribution parameter beta
  qga = par[3]  # quasi-Weibull distribution parameter gamma
  negLL = -n*log(qal+qga)+qal*n*log(qbe)+n*log(Igamma(qal/(qal+qga),(1/qbe)^(qal+qga),lower=F))-(qal-1)*sum(log(x/xmin))+sum(((x/xmin)/qbe)^(qal+qga))
  return(negLL)
}

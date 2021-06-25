.make_ao <- function(n) {
  ao <- rep(0, n)
  ao[sample(1:n,1)] <- 1
  return(ao)
}

.make_ls <- function(n) {
  ls <- rep(0, n)
  ls[sample(2:n,1):n] <- 1
  return(ls)
}

.make_tc <- function(n, delta=0.7) {
  tc <- rep(0, n)
  s <- sample(2:n,1):n
  tc[s] <- 1 * delta^(1:length(s))
  return(tc)
}

.make_random <- function(n, type=c("AO", "LS", "TC")) {  # Use randomly one of the make_yy functions
  f <- match.fun(paste0(".make_", sample(tolower(type),1)))
  return(f(n))
}

#' Simulate an outlier
#' 
#' Simulate an outlier
#' @param n Time series length
#' @param k Number of outliers
#' @param freq Frequency of the time series
#' @param type Type of outlier
#' @param effect_size Mean size of outlier 
#' @param start Start date of output time series
#' @param multiplicative Boolean. Is multiplicative time series model assumed?
#' @return The function returns k time series of class \code{xts} containing the k outlier effects
#' @examples plot(sim_outlier(60, 4, type=c("AO", "LS")))
#' @author Daniel Ollech
#' @details Three types of outliers are implemented: AO=Additive outlier, LS=Level shift, TC=Temporary Change. The effect size is stochastic as it is drawn from a normal distribution with mean equal to the specified effect_size and a standard deviation of 1/4*effect_size. This is multiplied randomly with -1 or 1 to get negative shocks as well.
#' If multiplicative is true, the effect size is measured in percentage. If is not true, the effect size is unit less and thus adopts the unit of the time series the outliers are added to.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export
#' 
sim_outlier <- function(n, k, freq=12, type=c("AO", "LS", "TC"), effect_size=10, start=c(2020,1), multiplicative=TRUE) {
  pure_outliers <- sapply(1:k, function(x) .make_random(n=n, type=type))
  outlier_effect <- stats::rnorm(k, effect_size, effect_size/4) * (stats::rbinom(k, 1, 0.5)*2-1)
  outlier <- stats::ts(t(outlier_effect * t(pure_outliers)), freq=freq, start=start)
  if (multiplicative) outlier <- 1+outlier/100
  
  return(outlier)
}

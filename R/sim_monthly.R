#' Simulate a monthly seasonal series
#' 
#' Simulate a monthly seasonal series
#' @param N Length in years
#' @param sd Standard deviation for all seasonal factors
#' @param change_sd Standard deviation of shock to seasonal factor
#' @param beta_1 Persistance wrt to previous period of the seasonal change 
#' @param beta_tau Persistence wrt to one year/cycle of the seasonal change
#' @param moving Is the seasonal pattern allowed to change over time
#' @param extra_smooth Boolean. Should the seasonal factors be smooth on a period-by-period basis
#' @param start Start date of output time series
#' @param multiplicative Boolean. Should multiplicative seasonal factors be simulated
#' @param model Model for non-seasonal time series. A list.
#' @return Multiple simulated monthly time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series without seasonal effects}
#' \item{sfac}{The seasonal effect}
#' }
#' @examples x=sim_monthly(5, multiplicative=TRUE)
#' ts.plot(x[,1])
#' @author Daniel Ollech
#' @details Standard deviation of the seasonal factor is in percent if a multiplicative time series model is assumed. Otherwise it is in unitless.
#' Using a non-seasonal ARIMA model for the initialization of the seasonal factor does not impact the seasonality of the time series. It can just make it easier for human eyes to grasp the seasonal nature of the series. The definition of the ar and ma parameter needs to be inline with the chosen model.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export

sim_monthly <- function(N, sd=5, change_sd = sd / 10, beta_1= 0.6, beta_tau = 0.4, moving = TRUE, model=list(order=c(3,1,1), ma=0.5, ar=c(0.2, -0.4, 0.1)), start=c(2010,1), multiplicative=TRUE, extra_smooth=FALSE) {
  # Basic series
  series <- stats::ts(stats::arima.sim(n=12*N-1, model=model), start=start, frequency=12)
  series <- series - min(series) + 100
  
  # Add seasonal effect
  sfac <- sim_sfac(length(series), freq=12, sd=sd, change_sd = change_sd, moving = moving, beta_1= beta_1, beta_tau = beta_tau, ar=0.99, ma=0.99, start=start, burnin=3, multiplicative=multiplicative, extra_smooth=extra_smooth)
  
  series <- tsbox::ts_xts(series)
  seas_adj <- series

  if (multiplicative) {
    sfac <- (sfac-1) / stats::sd(sfac) * sd/100 
    sfac <- 1+sfac 
    series <- series * as.numeric(sfac) } else {
      sfac <- sfac / stats::sd(sfac) * sd
      series <- series + as.numeric(sfac)   
    }
 
  # Post adjustment
  out <- xts::merge.xts(series, seas_adj, sfac)
  colnames(out) <- c("original", "seas_adj", "sfac")
  
  return(out)
  
}





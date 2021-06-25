#' Simulate a seasonal factor
#' 
#' Simulate a seasonal factor
#' @param n Number of observations
#' @param freq Frequency of the time series
#' @param sd Standard deviation of the seasonal factor
#' @param change_sd Standard deviation of simulation change to seasonal factor
#' @param beta_1 Persistance wrt to previous period of the seasonal change 
#' @param beta_tau Persistance wrt to one year/cycle before of the seasonal change
#' @param start Start date of output time series
#' @param multiplicative Boolean. Should multiplicative seasonal factors be simulated
#' @param ar AR parameter
#' @param ma MA parameter
#' @param model Model for initial seasonal factor
#' @param sc_model Model for the seasonal change
#' @param smooth Boolean. Should initial seasonal factor be smoothed
#' @param burnin (burnin*n-n) is the burn-in period
#' @param extra_smooth Boolean. Should the seasonal factor be smooth on a period-by-period basis
#' @return The function returns a time series of class \code{ts} containing a seasonal or periodic effect.
#' @examples ts.plot(sim_sfac(60))
#' @author Daniel Ollech
#' @details Standard deviation of the seasonal factor is in percent if a multiplicative time series model is assumed. Otherwise it is in unitless.
#' Using a non-seasonal ARIMA model does not impact the seasonality of the time series. It can just make it easier for human eyes to grasp the seasonal nature of the series. The definition of the ar and ma parameter needs to be inline with the chosen model.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export

sim_sfac <- function(n, freq=12, sd = 1, change_sd = 0.02, beta_1=0.9, beta_tau=0, start=c(2020,1), multiplicative=TRUE, ar=NULL, ma=NULL, model=c(1,1,1), sc_model=list(order=c(1,1,1), ar=0.65, ma=0.25), smooth=TRUE, burnin=3, extra_smooth=FALSE) {
  
if (sd < 0) {stop("Obviously, the standard deviation of the seasonal factor cannot be less than 0")}
  
if (sd==0) {
  z_t <- stats::ts(rep(0, n) + ifelse(multiplicative, 1, 0), start=start, freq=freq)
  return(z_t)
}

# Initial seasonal factor
if (is.null(ar) | is.null(ma)) {
  e_t <- stats::rnorm(freq , 0, sd)
  e_t <- e_t - mean(e_t) 
} else {
  e_t <- utils::tail(stats::arima.sim(n=freq*3,model=list(order=model, ar=ar, ma=ma)), freq)
  e_t <- e_t - mean(e_t)
  e_t <- (e_t / stats::sd(e_t) * sd) 
}

if(smooth){
  e_t <- stats::smooth.spline(c(e_t, e_t), spar=0.75 - 2*log(freq)/freq)$y[1:length(e_t)]
}

e_t <- e_t[1:freq]

y_t <- c(e_t, rep(NA, n+n*burnin))


# Creating a moving seasonal pattern
if (change_sd > 0) {
seasonal_change <- utils::tail(stats::arima.sim(n=(freq-1)*3,model=sc_model), freq)
seasonal_change <- seasonal_change / sd(seasonal_change) * sd * change_sd
seasonal_change <- seasonal_change - mean(seasonal_change) # Normalization of change
counter <- 0

for (j in freq+1:length(y_t)) {
 if (counter >=freq) {counter <- 1} else {counter <- counter + 1}
  index <- ifelse(counter>1, counter-1, freq)
  seasonal_change[counter] <- beta_1 * seasonal_change[counter] + beta_tau * seasonal_change[index] + stats::rnorm(1,0,sd*change_sd) 
  
  y_t[j] <- y_t[j-freq] + seasonal_change[counter]
  
  
  
  if (extra_smooth & counter==freq) {
    y_t[(j-freq+1):(j)] <- stats::smooth.spline(y_t[(j-freq+1):(j)], spar=0.75 - 2*log(freq)/freq)$y
    }
  
  
}
} else {
  y_t <- rep(e_t, length(y_t)/freq)
}

# Moving seasonal pattern from start to random point (so that the highest point does not necessarily tend to be the last observation in a period)
if (burnin > 0) {y_t <- y_t[1:(length(y_t)-round(stats::runif(1, -0.4999, freq-0.5001)))]}

# Ensuring that the (moving) annual mean is close to 0 
mafilterFreq <- rep(1/freq,freq)
mafilter2 <- rep(1/2, 2)

x1 <- stats::filter(y_t, mafilterFreq, sides=1)
x2 <- stats::filter(x1, mafilter2, sides=1)
z_t <- y_t - x2  
z_t <- utils::tail(z_t, n) - mean(utils::tail(z_t, n)) 

z_t <- stats::ts(z_t, start=start, frequency = freq)

# Ensuring that the standard deviation is as specified
z_t <- (z_t / stats::sd(z_t) * sd)

# Ensuring that mean of all complete years is zero
complete_years <- length(z_t) %/% freq
z_t <- z_t - mean(z_t[1:(complete_years*freq)])

# Changing additive seasonal factors to multiplicative
if (multiplicative) z_t <- 1 + z_t / 100



return(z_t)
}

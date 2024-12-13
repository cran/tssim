#' Simulate a daily time series based on the HS model 
#' 
#' This function simulates a daily time series with a Monte Carlo simulation based 
#' on an STS model based on Harvey and Shephard (1993) (HS model). 
#' The daily data consists of a trend, weekly seasonal, annual seasonal and 
#' irregular component. The components are each simulated by a transition process 
#' with daily random shocks. At the end of the simulation the components are combined and 
#' normalized to form the complete time series.
#' @param N Length of the simulated time series in years.
#' @param multiplicative If TRUE, a multiplicative model is simulated, an additive model if FALSE.
#' @param sizeWeeklySeas Size and stability of the weekly seasonal factor.
#' @param sizeAnnualSeas Size and stability of the annual seasonal factor.
#' @param sizeTrend Size of the trend component.
#' @param sizeDrift Size of the drift of the trend component.
#' @param varIrregularity Variance of the random irregular component.
#' @param sizeWeeklySeasAux Size of the auxiliary variable for the weekly seasonal factor.
#' @param sizeAnnualSeasAux size of the auxiliary variable for the annual seasonal factor.
#' @param start The initial date or year.
#' @param sizeBurnIn Size of burn-in sample in days.
#' @param shockLevel Variance of the shock to the level (trend).
#' @param shockDrift Variance of the shock to the drift (trend).
#' @param shockWeeklySeas Variance of the shock to the weekly seasonal.
#' @param shockAnnualSeas Variance of the shock to the annual seasonal.
#' @param index A value to which the mean of the base year (first effective year) of the time series is normalized.
#' @return Multiple simulated daily time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series seasonal effects}
#' \item{sfac7}{The day-of-the-week effect}
#' \item{sfac365}{The day-of-the-year effect}
#' }
#' @examples x <- sim_daily_hs(4)
#' ts.plot(x[,1])
#' @author Nikolas Fritz , Daniel Ollech, based on code provided by Ángel Cuevas and Enrique M Quilis
#' @details The size of the components and the variance of the irregular component are defaulted to 100 each 
#' and the variances of the shocks are defaulted to 1.
#' @details The first effective year serves as base year for the time series
#' @details The impact of a seasonal factor on the time series depends on its ratio to 
#' the other components. To increase (decrease) a factor's impact, the value of 
#' the size needs to be increased (decreased) while the other components need to be 
#' kept constant. Therefore, the stability of the seasonal factor also grows as 
#' the shocks on the given component have less impact. 
#' In order to increase the impact without increasing the stability, the variance 
#' of the shock also needs to be raised accordingly.
#' @references Cuevas, Ángel and Quilis, Enrique M., Seasonal Adjustment Methods for Daily Time Series. A Comparison by a Monte Carlo Experiment (December 20, 2023). Available at SSRN: https://ssrn.com/abstract=4670922 or http://dx.doi.org/10.2139/ssrn.4670922
#' @references Structural Time Series (STS) Monte Carlo simulation Z = trend + seasonal_weekly + seasonal_annual + irregular, according to Harvey and Shephard (1993): "Structural Time Series Models", in Maddala, G.S., Rao, C.R. and Vinod, H.D. (Eds.) Handbook of Statistics, vol. 11, Elsevier Science Publishers.
#' @export

sim_daily_hs <- function(N,
                         multiplicative = TRUE,
                         sizeWeeklySeas = 100,
                         sizeAnnualSeas = 100,
                         sizeTrend = 100,
                         sizeDrift = 100,
                         varIrregularity = 100,
                         sizeWeeklySeasAux = 100,
                         sizeAnnualSeasAux = 100,
                         start = 2020,
                         sizeBurnIn = 730,
                         shockLevel = 1,
                         shockDrift = 1,
                         shockWeeklySeas = 1,
                         shockAnnualSeas = 1,
                         index = 100) {
  # Adjustment of parameters and the variances of the shock
  sizeWeeklySeas <- sizeWeeklySeas / 10
  sizeAnnualSeas <- sizeAnnualSeas / 2
  sizeTrend <- sizeTrend * 10
  sizeDrift <- sizeDrift / 50000
  sizeWeeklySeasAux <- sizeWeeklySeasAux / 10
  sizeAnnualSeasAux <- sizeAnnualSeasAux / 2
  varIrregularity <- varIrregularity / 5
  shockLevel <- shockLevel / 1000
  shockDrift <- shockDrift / 100000
  shockWeeklySeas <- shockWeeklySeas / 100
  shockAnnualSeas <- shockAnnualSeas / 10
  
  # Combining components for the initial condition
  initialConditions = c(
    sizeTrend,
    sizeDrift,
    sizeWeeklySeas,
    sizeWeeklySeasAux,
    sizeAnnualSeas,
    sizeAnnualSeasAux
  )
  
  # Size of effective sample
  sizeEffective <-
    N * 366  # The overestimated days will be deleted later
  
  # Effective number of observations
  sizeTotal <- sizeBurnIn + sizeEffective
  
  # Dimension of the state vector
  ks <- length(initialConditions)
  
  # VCV matrix of transition equation shock
  vcvMatrix <-
    diag(c(
      shockLevel,
      shockDrift,
      rep(shockWeeklySeas, 2),
      rep(shockAnnualSeas, 2)
    ))
  
  # Implicit frequencies (in rads.)
  frequencyW <- (2 / 7) * pi
  frequencyA <- (2 / 365) * pi
  
  # Transition matrix: trend
  tmTrend <- matrix(1,
                    nrow = 2,
                    ncol = 2)
  tmTrend[2, 1] <- 0
  
  # Transition matrix: weekly seasonality
  tmWeekly <- diag(rep(cos(frequencyW), 2))
  tmWeekly[1, 2] <- sin(frequencyW)
  tmWeekly[2, 1] <- -tmWeekly[1, 2]
  
  # Transition matrix: annual seasonality
  tmAnnual <- diag(rep(cos(frequencyA), 2))
  tmAnnual[1, 2] <- sin(frequencyA)
  tmAnnual[2, 1] <- -tmAnnual[1, 2]
  
  # Transition matrix: 2x2 zero matrix
  tmZero <- matrix(0, ncol = 2, nrow = 2)
  
  # Combining block matrices to create the transition matrix and an additional transition matrix for leap years
  transitionMatrix <- rbind(
    cbind(tmTrend, tmZero, tmZero),
    cbind(tmZero, tmWeekly, tmZero),
    cbind(tmZero, tmZero, tmAnnual)
  )
  
  # Delete auxiliary matrices
  rm(tmTrend, tmWeekly, tmAnnual, tmZero)
  
  # State (measurement) equation
  selectComponents <- as.matrix(rep(c(1, 0), 3))
  
  # Monte Carlo simulation
  
  # Pre-allocation of state vector
  stateVector <- matrix(NaN,
                        nrow = sizeTotal + 1,
                        ncol = ks)
  # Initial condition
  stateVector[1,] <- initialConditions
  
  # Shock simulation
  shock <- mvtnorm::rmvnorm(n = sizeTotal + 1,
                            mean = rep(0, ks),
                            sigma = vcvMatrix)
  # Burn In sample
  for (t in seq(2, sizeBurnIn)) {
    stateVector[t,] <-
      transitionMatrix %*% stateVector[t - 1,] + shock[t,]
  }
  
  # Time loop with consideration of leap days
  c <- 0
  j <- sizeBurnIn + 1
  while (c < N) {
    if ((start[1] + c) %% 4 == 0) {
      k <- j + 365
    }
    else {
      k <- j + 364
    }
    
    for (t in seq(j, k)) {
      stateVector[t,] <-
        transitionMatrix %*% stateVector[t - 1,] + shock[t,]
    }
    
    # Adding the annual factor of a leap day
    if ((start[1] + c) %% 4 == 0) {
      leapDayFactor <-
        1 / 2 * (stateVector[j + 58, 5] + stateVector[j + 59, 5])  # Annual factor of the leap day as the mean of 02-28 and 03-01
      stateVector[(j + 60):(k), 5] <-
        stateVector[(j + 59):(k - 1), 5]  # Downshift of the days after the leap day
      stateVector[j + 59, 5] <-
        leapDayFactor  # Set leap day factor
      stateVector[k, 6] <-
        stateVector[k - 1, 6]  # Adjust auxiliary variable of the annual seasonal factor
    }
    
    c <- c + 1
    j <- k + 1
  }
  
  # Deleting initial condition
  stateVector <- stateVector[-1,]
  sizeBurnIn <- sizeBurnIn - 1
  
  # Deleting unused days
  for (t in k:sizeTotal) {
    stateVector <- stateVector[-k,]
    shock <- shock[-k,]
  }
  
  # Pre-allocations
  components <- matrix(NaN,
                       nrow = k - 1,
                       ncol = 4)  # Components: trend, 2 seasonals, irregularity
  series <- matrix(NaN,
                   nrow = k - 1,
                   ncol = 1)  # Aggregate series = trend + seasonals + irregularity
  
  # Measurement error (irregularity) simulation
  irreg <- stats::rnorm(k - 1,
                        sd = varIrregularity)
  
  # State equation
  series <- stateVector %*% selectComponents + irreg
  
  # Selecting components from Y and e
  components <- cbind(stateVector[, c(1, 3, 5)], irreg)
  
  # End of Monte Carlo simulation
  
  # Deleting burn-in sample
  components <- components[-(1:sizeBurnIn), ]
  series <- series[-(1:sizeBurnIn)]
  
  # Normalizing series by setting the first effective year as base year and setting its mean to the index (100)
  if (start[1] %% 4 == 0) {
    meanBaseYear <- mean(series[1:366])
  }
  else {
    meanBaseYear <- mean(series[1:365])
  }
  normSeries <- series * index / meanBaseYear
  
  # Adjust components
  adjComponents <- (normSeries / series) * components
  
  # Implement normalized series and adjusted components
  series <- normSeries
  components <- adjComponents
  
  # Transition into multiplicative model
  if (multiplicative == TRUE) {
    #    components[,1] <- components[,1] * 100 / mean(components[,1])
    components[, (2:4)] <- (components[, (2:4)] + 100) / 100
    
    # Combining trend component and irregularity
    components[, 1] <- components[, 1] * components[, 4]
    
    # Build time series
    series <- components[, 1] * components[, 2] * components[, 3]
  }
  else {
    # Combining trend component and irregularity
    components[, 1] <- components[, 1] + components[, 4]
  }
  
  # Remove vector of the irregular component
  components <- components[, -4]
  
  # Building the complete time series including its components
  data <-
    stats::ts(cbind(series, components),
              start = start,
              frequency = 365.2425)
  colnames(data) <- c('original', 'seas_adj', 'sfac7', 'sfac365')
  
  # Convert into xts format
  out <- tsbox::ts_xts(data)
  
  return(out)
}

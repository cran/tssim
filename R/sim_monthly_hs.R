#' Simulate a monthly time series based on the HS model
#' 
#' This function simulates a monthly time series with a Monte Carlo simulation based 
#' on an STS model based on Harvey and Shephard (1993) (HS model). 
#' The monthly data consists of a trend, annual seasonal and 
#' irregular component. The components are each simulated by a transition process 
#' with monthly random shocks and then combined at the end of the simulation to form the complete time series.
#' @param N Length of the simulated time series in years.
#' @param multiplicative If true, a multiplicative model is simulated, an additive model if FALSE.
#' @param sizeSeasonality Size and stability of the annual seasonal factor.
#' @param sizeTrend Size and stability of the trend component.
#' @param sizeDrift Size and stability of the drift of the trend component.
#' @param sizeSeasonalityAux Size of the auxiliary variable for the annual seasonal factor.
#' @param varIrregularity Variance of the random irregular component.
#' @param start The initial date or year.
#' @param sizeBurnIn Size of burn-in sample in months.
#' @param shockLevel Variance of the shock to the level (trend).
#' @param shockDrift Variance of the shock to the drift (trend).
#' @param shockSeasonality Variance of the shock to the annual seasonal.
#' @param index A value to which the mean of the base year (first effective year) of the time series is normalized.
#' @return Multiple simulated monthly time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series without seasonal effects}
#' \item{sfac}{The seasonal effect}
#' }
#' @examples x <- sim_monthly_hs(4)
#' ts.plot(x[,1])
#' @author Nikolas Fritz, Daniel Ollech, based on code provided by Ángel Cuevas and Enrique M Quilis
#' @details The impact of a component on the time series depends on its ratio to 
#' the other components. To increase (decrease) a component's impact, the value of 
#' the size needs to be increased (decreased) while the other components need to be 
#' kept constant. Therefore, the stability of the component (e.g. the shape of a seasonal component) 
#' also grows as the shocks on the given component have less impact. 
#' In order to increase the impact without increasing the stability, the variance 
#' of the shock also needs to be raised accordingly. 
#' The size of the components are defaulted to 100 each and the variances of the shocks are defaulted to 1.
#' @details The first effective year serves as base year for the time series
#' @references Cuevas, Ángel and Quilis, Enrique M., Seasonal Adjustment Methods for Daily Time Series. A Comparison by a Monte Carlo Experiment (December 20, 2023). Available at SSRN: https://ssrn.com/abstract=4670922 or http://dx.doi.org/10.2139/ssrn.4670922
#' @references Structural Time Series (STS) Monte Carlo simulation Z = trend + seasonal_weekly + seasonal_annual + irregular, according to Harvey and Shephard (1993): "Structural Time Series Models",in Maddala, G.S., Rao, C.R. and Vinod, H.D. (Eds.) Handbook of Statistics, vol. 11, Elsevier Science Publishers.
#' @export

sim_monthly_hs <- function(N,
                           multiplicative = TRUE,
                           sizeSeasonality = 100,
                           sizeTrend = 100,
                           sizeDrift = 100,
                           sizeSeasonalityAux = 100,
                           varIrregularity = 1,
                           start = 2020,
                           sizeBurnIn = 24,
                           shockLevel = 1,
                           shockDrift = 1,
                           shockSeasonality = 1,
                           index = 100) {
  
  # Adjustment of parameters and the variances of the shock
  sizeSeasonality <- sizeSeasonality/10
  sizeTrend <- sizeTrend*5
  sizeDrift <- sizeDrift/100000
  sizeSeasonalityAux <- sizeSeasonalityAux/10
  varIrregularity <- varIrregularity*15
  shockLevel <- shockLevel/100
  shockDrift <- shockDrift/1000
  shockSeasonality <- shockSeasonality*4
  
  # Combining components for the initial condition
  initialConditions = c(sizeTrend, sizeDrift, sizeSeasonality, sizeSeasonalityAux)
  
  # Size of effective sample 
  sizeEffective <- N * 12
  
  # Total number of observations
  sizeTotal <- sizeBurnIn + sizeEffective
  
  # Dimension of the state vector
  ks <- length(initialConditions)
  
  # VCV matrix of transition equation shock
  vcvMatrix <- diag(c(shockLevel, shockDrift, rep(shockSeasonality,2)))
  
  # Implicit frequencies (in rads.)
  frequencyA <- (2/12)*pi
  
  # Transition matrix: trend
  tmTrend <- matrix(1, 
                    nrow = 2, 
                    ncol = 2)
  tmTrend[2,1] <- 0
  
  # Transition matrix: annual seasonality
  tmSeas <- diag(rep(cos(frequencyA), 2))
  tmSeas[1,2] <- sin(frequencyA)
  tmSeas[2,1] <- -tmSeas[1,2]
  
  # Transition matrix: 2x2 zero matrix
  tmZero <- matrix(0, ncol = 2, nrow = 2)
  
  # Combining block matrices to create the transition matrix
  transitionMatrix <- rbind(cbind(tmTrend, tmZero),
                            cbind(tmZero, tmSeas))
  
  # Delete auxiliary matrices
  rm(tmTrend, tmSeas, tmZero)
  
  # State (measurement) equation
  selectComponents <- as.matrix(rep(c(1,0), 2))
  
  # Monte Carlo simulation
  
  # Pre-allocation of state vector
  stateVector <- matrix(NaN, 
                        nrow = sizeTotal + 1, 
                        ncol = ks)
  # Initial condition
  stateVector[1, ] <- initialConditions
  
  # Shock simulation
  shock <- mvtnorm::rmvnorm(n = sizeTotal+1,
                            mean = rep(0, ks),
                            sigma = vcvMatrix)
  
  # Simulation of components
  for (t in seq(2,sizeTotal+1)) {
    stateVector[t, ] <- transitionMatrix %*% stateVector[t-1, ] + shock[t, ]
  }
  
  # Deleting initial condition
  stateVector <- stateVector[-1, ]
  
  # Pre-allocations
  components <- matrix(NaN,
                       nrow = sizeTotal,
                       ncol = 3)  # Components: trend, annual seasonal, irregularity
  series <- matrix(NaN, 
                   nrow = sizeTotal, 
                   ncol = 1)  # Aggregate series = trend + annual seasonal + irregularity
  
  # Measurement error (irregularity) simulation 
  irreg <- stats::rnorm(sizeTotal, 
                        sd = varIrregularity)
  
  # State equation
  series <- stateVector %*% selectComponents + irreg
  
  # Selecting components from Y and e
  components <- cbind(stateVector[,c(1,3)], irreg)
  
  # End of Monte Carlo simulation
  
  # Deleting burn-in sample
  components <- components[-(1:sizeBurnIn),]
  series <- series[-(1:sizeBurnIn)]
  
  # Normalizing series by setting the first effective year as base year and setting its mean to the index (100)
  meanBaseYear <- mean(series[1:12])
  norSeries <- series * index / meanBaseYear
  
  # Adjust components
  adjComponents <- (norSeries / series) * components
  
  # Implement normalized series and adjusted components
  series <- norSeries
  components <- adjComponents
  
  # Transition into a multiplicative model
  if (multiplicative == TRUE) {
#    components[,1] <- components[,1] * 100 / mean(components[,1])
    components[,2] <- (components[,2] + 100) / 100
    components[,3] <- (components[,3] + 100) / 100
    series <- components[,1] * components[,2] * components[,3] 
    
    # Combining trend component and irregularity
    components[,1] <- components[,1] * components[,3]
  }
  else {
    # Combining trend component and irregularity
    components[,1] <- components[,1] + components[,3]
  }
  
  # Remove vector of the irregular component
  components <- components[,-3]
  
  # Building the complete time series including its components
  data <- stats::ts(cbind(series, components), start=start, frequency=12)
  colnames(data) <- c('original', 'seas_adj', 'sfac')
  
  # Convert into xts format
  out <- tsbox::ts_xts(data)
  
  return(out)
}

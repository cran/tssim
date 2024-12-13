#' Daily time series simulation for the MSTL-algorithm
#'
#' This function simulates a daily time series according to the simulation model of
#' Bandara, Hyndman and Bergmeir (2021) about the MSTL-algorithm for seasonal-trend decomposition.
#' The simulated time series consists of a trend, weekly, annual and irregular component which are
#' each simulated independently from each other. After the simulation process they are normalized and then combined
#' to form the complete time series. As in the paper, this simulation function has the option to distinguish between
#' a deterministic and a stochastic data generation process.
#' @param N length in years
#' @param multiplicative If TRUE, a multiplicative model is simulated, if FALSE, the model is additive
#' @param start Start year or start date of the simulation.
#' @param sizeAnnualSeas Size of the annual seasonal factor, defaulted to 100.
#' @param sizeWeeklySeas Size of the weekly seasonal factor, defaulted to 100.
#' @param sizeIrregularity Size of the irregular component, defaulted to 100.
#' @param shockAnnualSeas Shock to the annual seasonal coefficient, defaulted to 1.
#' @param shockWeeklySeas Shock to the weekly seasonal coefficient, defaulted to 1.
#' @param deterministic If TRUE, the seasonal coefficients are deterministic, meaning they do not change after a seasonal cycle. If FALSE, the coefficients are stochastic, meaning they change randomly after a seasonal cycle.
#' @return Multiple simulated daily time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series without seasonal effects}
#' \item{sfac7}{The day-of-the-week effect}
#' \item{sfac365}{The day-of-the-year effect}
#' }
#' @examples x <- sim_daily_mstl(4)
#' ts.plot(x[,1])
#' @author Nikolas Fritz, Daniel Ollech
#' @references Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021). MSTL: A seasonal-trend decomposition
#' algorithm for time series with multiple seasonal patterns. arXiv preprint arXiv:2107.13462.
#' @export

sim_daily_mstl <- function(N,
                           multiplicative = TRUE,
                           start = 2020,
                           sizeAnnualSeas = 100,
                           sizeWeeklySeas = 100,
                           sizeIrregularity = 100,
                           shockAnnualSeas = 1,
                           shockWeeklySeas = 1,
                           deterministic = FALSE) {
  # Adjustment of the parameters
  sizeAnnualSeas <- sizeAnnualSeas / 33
  sizeWeeklySeas <- sizeWeeklySeas / 100
  sizeIrregularity <- sizeIrregularity / 150
  shockAnnualSeas <- shockAnnualSeas / 10
  shockWeeklySeas <- shockWeeklySeas / 20
  
  # Size of sample
  sizeSample <- N * 366  # Overestimated days will be deleted later
  
  # Trend component
  
  if (deterministic == TRUE) {
    # Random variables
    randomVar1 <- stats::rnorm(1)
    randomVar2 <- stats::rnorm(1)
    
    # Pre-allocation
    t <- 1:sizeSample
    trend <- matrix(NaN,
                    nrow = sizeSample,
                    ncol = 1)
    
    # Deterministic trend simulation
    trend <-
      randomVar1 * (t + sizeSample / 2 * (randomVar2 - 1)) ^ 2
  }
  
  else {
    # Stochastic trend simulation
    trend <-
      stats::arima.sim(model = list(order = c(0, 2, 0)), sizeSample)
    trend <- trend[3:(sizeSample + 2)]
  }
  
  # Seasonal components
  
  # Fourier terms pre-allocation
  fourierAnnual <- stats::ts(rep(0, sizeSample), frequency = 365)
  fourierWeekly <- stats::ts(rep(0, sizeSample), frequency = 7)
  kA <- 5
  kW <- 2
  
  # Generate Fourier terms
  fourierAnnual <- forecast::fourier(fourierAnnual, K = kA)
  fourierWeekly <- forecast::fourier(fourierWeekly, K = kW)
  
  # Seasonal coefficients
  coefAnnual <- stats::rnorm(2 * kA)
  coefWeekly <- stats::rnorm(2 * kW)
  
  # Pre-allocations seasonal components
  annualComp <- rep(0, sizeSample)
  weeklyComp <- rep(0, sizeSample)
  
  # Helpers
  c <- 0
  j <- 1
  leap <- 0
  
  # Annual component
  while (c < N) {
    k <- j + 364
    jj <-
      j - leap  # Correction after a leap year to use the correct values of the fourier terms
    kk <- k - leap
    
    # Annual Component in the given year
    annualComp[j:k] <- fourierAnnual[jj:kk,] %*% coefAnnual
    
    # Calculation of leap days
    if ((start[1] + c) %% 4 == 0) {
      k <- k + 1
      leap <- leap + 1
      leapDayFactor <-
        1 / 2 * (annualComp[j + 58] + annualComp[j + 59])  # Annual factor of the leap day as the mean of 02-28 and 03-01
      annualComp[(j + 60):(k)] <-
        annualComp[(j + 59):(k - 1)]  # Downshift of the days after the leap day
      annualComp[j + 59] <- leapDayFactor  # Set leap day factor
    }
    
    # Stochastic annual component
    if (deterministic == FALSE) {
      # Shock on the annual coefficient
      shock <- stats::rnorm(2 * kA, sd = shockAnnualSeas)
      coefAnnual <- coefAnnual + shock
    }
    
    j <- k + 1
    c <- c + 1
  }
  
  # Weekly component
  
  if (deterministic == TRUE) {
    # Deterministic weekly component
    weeklyComp <- fourierWeekly %*% coefWeekly
  }
  
  else {
    # Stochastic weekly component
    c <- 1
    j <- 1
    cor <- sizeSample %% 7
    weeks <- (sizeSample - cor) / 7
    while (c < weeks) {
      k <- j + 6
      weeklyComp[j:k] <- fourierWeekly[j:k,] %*% coefWeekly
      
      # Shock on the weekly coefficient
      shock <- stats::rnorm(2 * kW, sd = shockWeeklySeas)
      coefWeekly <- coefWeekly + shock
      
      j <- k + 1
      c <- c + 1
    }
    
    # Add the days of the last week of the sample
    weeklyComp[j:sizeSample] <-
      fourierWeekly[j:sizeSample,] %*% coefWeekly
  }
  
  # Delete overestimated days
  overEst <- k + 1
  trend <- trend[-(overEst:sizeSample)]
  annualComp <- annualComp[-(overEst:sizeSample)]
  weeklyComp <- weeklyComp[-(overEst:sizeSample)]
  
  # Normalize trend and seasonal components
  trend <- scale(trend)
  annualComp <- scale(annualComp)
  weeklyComp <- scale(weeklyComp)
  
  # Simulate Irregularity
  irregular <- matrix(stats::rnorm(k),
                      ncol = 1)
  
  # Components and series
  if (multiplicative == TRUE) {
    irregular <- sizeIrregularity * irregular + 100
    seasAdj <- (trend + 100) * irregular / 100
    sfac365 <- (sizeAnnualSeas * annualComp + 100) / 100
    sfac7 <- (sizeWeeklySeas * weeklyComp + 100) / 100
    
    series <- seasAdj * sfac365 * sfac7
  }
  else {
    seasAdj <- trend + sizeIrregularity * irregular + 100
    sfac365 <- sizeAnnualSeas * annualComp
    sfac7 <- sizeWeeklySeas * weeklyComp
    
    series <- seasAdj + sfac7 + sfac365
  }
  
  # Build data
  data <-
    stats::ts(cbind(series, seasAdj, sfac7, sfac365),
              start = start,
              frequency = 365.2425)
  colnames(data) <- c('original', 'seas_adj', 'sfac7', "sfac365")
  
  # Outcome
  out <- tsbox::ts_xts(data)
  
  return(out)
}

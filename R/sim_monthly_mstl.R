#' Monthly time series simulation for the MSTL-algorithm
#' 
#' This function simulates a monthly time series according to the simulation model of 
#' Bandara, Hyndman and Bergmeir (2021) about the MSTL-algorithm for seasonal-trend decomposition.
#' The simulated time series consists of a trend, annual seasonal and irregular component which are 
#' each simulated independently from each other. After the simulation process they are normalized and then combined
#' to form the complete time series. As in the paper, this simulation function has the option to distinguish between 
#' a deterministic and a stochastic data generation process.
#' @param N length in years
#' @param multiplicative If TRUE, a multiplicative model is simulated, if FALSE, the model is additive
#' @param start Start year or start date of the simulation.
#' @param sizeSeasonality Size of the annual seasonal factor.
#' @param sizeIrregularity Size of the irregular component.
#' @param sizeTrend Size of trend component.
#' @param shockSeasonality Variance of the shock to the annual seasonal coefficient, defaulted to 1.
#' @param deterministic If TRUE, the seasonal coefficients are deterministic, meaning they do not change after a seasonal cycle. If FALSE, the coefficients are stochastic, meaning they change by random shocks after a seasonal cycle.
#' @return Multiple simulated monthly time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series without seasonal effects}
#' \item{sfac}{The seasonal effect}
#' }
#' @examples x <- sim_monthly_mstl(4)
#' ts.plot(x[,1])
#' @author Nikolas Fritz, Daniel Ollech
#' @references Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021). MSTL: A seasonal-trend decomposition 
#' algorithm for time series with multiple seasonal patterns. arXiv preprint arXiv:2107.13462.
#' @export


sim_monthly_mstl <- function(N,
                           multiplicative = TRUE,
                           start = 2020,
                           sizeSeasonality = 100,
                           sizeIrregularity = 100,
                           sizeTrend = 100,
                           shockSeasonality = 1,
                           deterministic = FALSE){
  
  # Adjustment of the parameters
  sizeSeasonality <- sizeSeasonality / 80
  sizeIrregularity <- sizeIrregularity / 100
  sizeTrend <- sizeTrend / 200
  shockSeasonality <- shockSeasonality / 5
  
  # Size of sample
  sizeSample <- N * 12
  
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
    trend <- randomVar1 * (t + sizeSample/2 * (randomVar2 - 1))^2
  }
  
  else {
    
    # Stochastic trend simulation
    trend <- stats::arima.sim(model=list(order=c(0,2,0)), sizeSample)
    trend <- trend[3:(sizeSample+2)]
  }
  
  # Seasonal components
  
  # Fourier terms pre-allocation
  fourierAnnual <- stats::ts(rep(0, sizeSample), frequency = 12)
  kA <- 5
  
  # Generate Fourier terms
  fourierAnnual <- forecast::fourier(fourierAnnual, K = kA)

  # Seasonal coefficients
  coefAnnual <- stats::rnorm(2*kA)

  # Pre-allocations seasonal component
  annualComp <- rep(0, sizeSample)
  
  # Helpers
  c <- 0
  j <- 1
  
  # Seasonal component
  while (c < N) {
    k <- j + 11
    
    # Seasonal Component in the current year
    annualComp[j:k] <- fourierAnnual[j:k,] %*% coefAnnual  
    
    # Stochastic seasonal component
    if (deterministic == FALSE) {
      # Shock on the seasonal coefficient
      shock <- stats::rnorm(2*kA, sd = shockSeasonality)
      coefAnnual <- coefAnnual + shock
    }
    
    j <- k + 1
    c <- c + 1
  }
  
  # Normalize trend and seasonal component
  trend <- scale(trend)
  annualComp <- scale(annualComp)
  
  # Simulate Irregularity
  irregular <- matrix(stats::rnorm(k),
                      ncol = 1)
  
  # Components and series
  if (multiplicative == TRUE) {
    irregular <- sizeIrregularity * irregular + 100
    seasAdj <- (trend + 100) * irregular / 100
    sfac12 <- (sizeSeasonality * annualComp + 100) / 100

    series <- seasAdj * sfac12
  }
  else {
    seasAdj <- trend + sizeIrregularity * irregular + 100
    sfac12 <- sizeSeasonality * annualComp

    series <- seasAdj + sfac12
  }


  # Build data
  data <- stats::ts(cbind(series, seasAdj, sfac12), start=start, frequency=12)
  colnames(data) <- c('original', 'seas_adj', 'sfac')
  
  # Outcome
  out <- tsbox::ts_xts(data)
  
  return(out)
}

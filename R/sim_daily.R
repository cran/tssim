#' Simulate a daily seasonal series
#' 
#' Simulate a daily seasonal series as described in Ollech (2021).
#' @param N length in years
#' @param sd Standard deviation for all seasonal factors
#' @param change_sd Standard deviation of simulated change for all seasonal factors
#' @param week_sd Standard deviation of the seasonal factor for day-of-the-week
#' @param month_sd Standard deviation of the seasonal factor for day-of-the-month
#' @param year_sd Standard deviation of the seasonal factor for day-of-the-year
#' @param week_change_sd Standard deviation of simulated change to seasonal factor for day-of-the-week
#' @param month_change_sd Standard deviation of simulated change to seasonal factor for month-of-the-week
#' @param year_change_sd Standard deviation of simulated change to seasonal factor for year-of-the-week
#' @param innovations_sd Standard deviation of the innovations used in the non-seasonal regarima model
#' @param sa_sd Standard deviation of the non-seasonal time series
#' @param extra_smooth Boolean. Should the seasonal factors be smooth on a period-by-period basis
#' @param start Start date of output time series
#' @param multiplicative Boolean. Should multiplicative seasonal factors be simulated
#' @param model Model for non-seasonal time series. A list.
#' @param beta_1 Persistance wrt to previous period of the seasonal change 
#' @param beta_tau Persistance wrt to one year/cycle before of the seasonal change
#' @param calendar Parameters for calendar effect, a list, see sim_calendar
#' @param outlier Parameters for outlier effect, a list, see sim_outlier
#' @param timewarping Should timewarping be used to obtain the day-of-the-month factors
#' @param as_index Shall series be made to look like an index (i.e. shall values be relative to reference year = second year)
#' @return Multiple simulated daily time series of class xts including:
#' \describe{
#' \item{original}{The original series}
#' \item{seas_adj}{The original series without calendar and seasonal effects}
#' \item{sfac7}{The day-of-the-week effect}
#' \item{sfac31}{The day-of-the-month effect}
#' \item{cfac}{The calendar effects}
#' \item{outlier}{The outlier effects}
#' }
#' @examples x=sim_daily(5, multiplicative=TRUE, outlier=list(k=5, type=c("AO", "LS"), effect_size=50))
#' ts.plot(x[,1])
#' @author Daniel Ollech
#' @details Standard deviation of the seasonal factor is in percent if a multiplicative time series model is assumed. Otherwise it is in unitless.
#' Using a non-seasonal ARIMA model for the initialization of the seasonal factor does not impact the seasonality of the time series. It can just make it easier for human eyes to grasp the seasonal nature of the series. The definition of the ar and ma parameter needs to be inline with the chosen model.
#' If only change_sd is specified, the change parameters for the single seasonal factors are set individually as change_sd/365*(length of seasonal cycle)
#' The parameters that can be set for calendar and outlier are those defined in sim_outlier and sim_calendar.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export

sim_daily <- function(N, sd=2.5, change_sd=0.05, week_sd=NA, month_sd=NA, year_sd=NA, week_change_sd=NA, month_change_sd=NA, year_change_sd=NA, innovations_sd=1, sa_sd=NA,model=list(order=c(3,1,1), ma=0.5, ar=c(0.2, -0.4, 0.1)), beta_1=0.9, beta_tau=0, start=c(2020,1), multiplicative=TRUE, extra_smooth=FALSE, calendar=list(which="Easter", from=-2, to=2), outlier=NULL, timewarping=TRUE, as_index=FALSE) {
  if (is.na(week_sd)) {week_sd <- sd}
  if (is.na(month_sd)) {month_sd <- sd}
  if (is.na(year_sd)) {year_sd <- sd}
  
  if (is.na(week_change_sd)) {week_change_sd <- change_sd/365*7}
  if (is.na(month_change_sd)) {month_change_sd <- change_sd/365*31}
  if (is.na(year_change_sd)) {year_change_sd <- change_sd}
  
  # Basic series
  series <- stats::ts(stats::arima.sim(n=365*N-1, model=model, sd=innovations_sd), start=start, frequency=365)

  if (!is.na(sa_sd)) {
    series <- series * sa_sd / sd(series, na.rm=FALSE)
  }
  
  if (multiplicative & min(series, na.rm=TRUE) < 0) series <- series - min(series, na.rm=TRUE) + 100
  
  if(as_index & N >= 2) series <- series / mean(series[366:730], na.rm=TRUE) * 100

  # Add day-of-the-year-effect
  sfac365 <- sim_sfac(length(series), freq=365, sd=year_sd, change_sd=year_change_sd, beta_1=beta_1, beta_tau=beta_tau, ar=0.99, ma=0.99, start=start, burnin=3, multiplicative=multiplicative, extra_smooth=extra_smooth)
  
  series <- dsa::ts2xts(series)
  times <- seq.Date(from=stats::start(series), to=stats::end(series), by="days") # Adding 29.2
  blank <- xts::xts(rep(NA,length(times)), times)
  series <- zoo::na.approx(xts::merge.xts(series,blank)$series)
  seas_adj <- series
  
  sfa <- dsa::ts2xts(sfac365)
  sfac365 <- zoo::na.approx(xts::merge.xts(sfa, blank)$sfa)
  
  if(year_sd > 0){
  if (multiplicative) {
    sfac365 <- (sfac365-1) / stats::sd(sfac365) * year_sd/100 
    sfac365 <- 1+sfac365 
    series <- series * as.numeric(sfac365) } else {
      sfac365 <- sfac365 / stats::sd(sfac365) * year_sd
      series <- series + as.numeric(sfac365)   
    }}

  # Add day-of-the-week-effect
  sfac7 <- sim_sfac(length(series), freq=7, sd=week_sd, change_sd=week_change_sd, beta_1=beta_1, beta_tau=beta_tau, start=start, multiplicative=multiplicative, extra_smooth=extra_smooth)
  if (multiplicative) {
    series <- series * as.numeric(sfac7) } else {
      series <- series + as.numeric(sfac7)   
    }


  # Add day-of-the-month-effect
  series31 <- dsa:::.fill31(series,"spline") 
  sfac31 <- sim_sfac(length(series31), freq=31, sd=month_sd, change_sd=month_change_sd, beta_1=beta_1, beta_tau=beta_tau, ar=0.9, ma=0.99, start=start, burnin=1, multiplicative=multiplicative, extra_smooth=extra_smooth) 
  if (timewarping) {
    sfac31 <- .stretch_re(sfac31)
  } else {
  sfac31 <- dsa:::.drop31(stats::ts(sfac31, start=start, frequency=372), 1, 365)} # if number of years is dividable by 4 then it should probably be 366 instead of 365
  if(month_sd > 0){
  if (multiplicative) {
    sfac31 <- (sfac31-1) / stats::sd(sfac31) * month_sd/100 
    sfac31 <- 1+sfac31 
    series <- series * as.numeric(sfac31) } else {
      sfac31 <- sfac31 / stats::sd(sfac31) * month_sd
      series <- series + as.numeric(sfac31)   
    }}
  
  # Add calendar effect
  if (is.null(calendar)) {cfac <- series * 0 + as.numeric(multiplicative)} else {
  calendar[["n"]]=length(series)
  calendar[["multiplicative"]] <- multiplicative
  calendar[["start"]] <- as.Date(paste(start, collapse="-"), "%Y-%j")
  calendar[["freq"]] <- 365.25
  cfac <- suppressMessages(do.call(sim_calendar, calendar))}
  
  if (multiplicative) {series <- series*cfac} else {series <- series + cfac}
  
  # Add outlier
  if (is.null(outlier)) {outlier_effect <- series * 0 + as.numeric(multiplicative)} else {
    outlier[["n"]]=length(series)
    outlier[["multiplicative"]] <- multiplicative
    outlier[["start"]] <- as.Date(paste(start, collapse="-"), "%Y-%j")
    outlier_effect <- apply(do.call(sim_outlier, outlier), 1, ifelse(multiplicative,prod, sum)) # RowProduct
  }
  
  if (multiplicative) {series <- series*outlier_effect; seas_adj <- seas_adj * outlier_effect} else {series <- series + outlier_effect; seas_adj <- seas_adj + outlier_effect}
  
  
  
  # Post adjustment
  out <- xts::merge.xts(series, seas_adj, xts::xts(sfac7, zoo::index(series)), xts::xts(sfac31, zoo::index(series)), sfac365, cfac, outlier_effect)
  colnames(out) <- c("original", "seas_adj", "sfac7", "sfac31", "sfac365", "cfac", "outlier")
  
  return(out)
  
}


#' Use time warping to reduce the number of observations in a month
#' 
#' Reduce the number of observations in a month using time warping / stretching. Only relevant if a daily time series is simulated
#' @param  seas_component Seasonal component for day-of-the-month
#' @author Daniel Ollech
#' @details Usually time warping would be used to stretch the number of observations of a time series in a given interval to more observations. Here it is used to reduce the number of observations (31) to the number of days in a given month while maintaining the underlying trajectory of the data. This is done by first creating a very long time series for each month, interpolating missing values by spline interpolation and then reducing the number of observations to the number suitable for a given month.
#' @return Returns a \code{xts} time series containing the day-of-the-month effect.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export


.stretch_re <- function(seas_component) {
  
  start <- stats::start(seas_component)
  
  # Change time series to matrix, where each column is a given month
  Mat <- matrix(c(seas_component, rep(NA, 31-length(seas_component)%%31)), nrow=31)
  colnames(Mat) <- as.character(seq.Date(as.Date(paste0(start[1], "-", start[2], "-01")), by="months", length.out=ncol(Mat)))
  
  .transformX <- function(x, to=30) { # Reduce the observations in a given month from 31 to X
    z <- matrix(NA, nrow=to, ncol=length(x))
    z[1,1] <- x[1]
    z[to,] <- x
    zz <- zoo::na.spline(as.vector(z)) # using Spline interpolation
    
    y <- c(matrix(zz, nrow=31)[31,], rep(NA, 31-to))
    return(y)
  }
  
  .is.leapyear <- function(Year) { # Identify leap years
    Year <- as.numeric(Year)
    if(Year %% 4 == 0) {
      if(Year %% 100 == 0) {
        if(Year %% 400 == 0) {return(TRUE)} else {
          return(FALSE)
        }} else {return(TRUE)}
    } else {return(FALSE)}}
  
  .not.last <- function(X) { # Identify last month (that shall not be changed)
    out <- rep(1, ncol(X))
    out[ncol(X)] <- 0
    return(as.logical(out))
  }
  
  # Months that need to be changed
  month30 <- c("04", "06", "09", "11")
  month28 <- c("02")
  
  # Shortens months with 30, with 28 and 29 days by time warping /stretching. Last month is not changed
  
  ## 30-day-months
  Mat[,format(as.Date(colnames(Mat)), "%m") %in% month30 & .not.last(Mat)] <- apply(Mat[,format(as.Date(colnames(Mat)), "%m") %in% month30 & .not.last(Mat)], 2, .transformX, 30)
  
  ## 28-day-months
  Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & !sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat) )] <- apply(Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & !sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))], 2, .transformX, 28)
  
  ## 29-day-months (more complicated, because apply does not work with vectors, i.e. problems if only one leap year is in the data)
  if (sum((format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat) )) > 1) { 
    Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))] <- apply(Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))], 2, .transformX, 29)} else {
      if (sum((format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))) == 1) { 
        Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))] <- .transformX(Mat[,(format(as.Date(colnames(Mat)), "%m") %in% month28 & sapply(format(as.Date(colnames(Mat)), "%Y"), .is.leapyear) & .not.last(Mat))], 29)
      }}
  
  # Output
  raw <- stats::na.omit(as.vector(Mat))
  
  out <- xts::xts(raw, seq.Date(as.Date(paste0(start[1], "-", start[2], "-01")), by="days", length.out=length(raw)))
  
  return(out)
}





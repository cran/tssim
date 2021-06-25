#' Simulate calendar effects
#' 
#' Simulate a time series containing specified calendar effects 
#' @param n Time series length
#' @param which Holidays to be used, functions from timeDate package used
#' @param from days before the Holiday to include
#' @param to days after the Holiday to include
#' @param freq Frequency of the time series
#' @param effect_size Mean size of calendar effect
#' @param start Start Date of output time series
#' @param multiplicative Boolean. Is multiplicative time series model assumed?
#' @param time_dynamic Should the calendar effect change over time
#' @param center Should calendar variable be center, i.e. mean=0 
#' @return The function returns a time series of class \code{xts}
#' @examples plot(sim_calendar(60, from=0, to=4, freq=12))
#' @author Daniel Ollech
#' @details 
#' If multiplicative is true, the effect size is measured in percentage. If is not true, the effect size is unit less and thus adopts the unit of the time series the calendars are added to.
#' The time_dynamic parameter controls the change of the calendar effect. The effect of the previous year is multiplied by the time_dynamic factor.
#' @references Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
#' @export


sim_calendar <- function(n, which=c("Easter", "Ascension"), from=0, to=0, freq=12, effect_size=3, start="2020-01-01", multiplicative=TRUE, time_dynamic=1, center=TRUE) {
  # Create empty calendar variable
 Cal <- xts::xts(rep(0,ceiling(n/freq*366)), seq.Date(from=as.Date(start), by="days", length.out=ceiling(n/freq*366)))
 # Find dates 
  s <- as.numeric(format(utils::head(zoo::index(Cal),1)  , "%Y"))
  e <- as.numeric(format(utils::tail(zoo::index(Cal),1)  , "%Y"))
  dts <- sapply(which, function(x) timeDate::as.Date.timeDate(eval(parse(text=paste0("timeDate::", x, "(", s, ":",e, ")")))))
  # Move dates by from and to
  dates <- as.Date(sapply(from:to, function(x) dts + x), origin="1970-01-01")
  dates <- sort(dates)
  # Create the calendar effect
  Cal[dates] <- (rep(stats::rnorm(length(dates)/length(s:e), effect_size * (stats::rbinom(length(s:e), 1, 0.5)*2-1), effect_size/4), length(s:e))*rep(time_dynamic^(1:length(s:e)), each=length(which)))[1:length(Cal[dates])]
  
  # Potentially aggregate series
  if (freq==12) calendar <- xts::apply.monthly(Cal, sum)
  if (freq==4) calendar <- xts::apply.quarterly(Cal, sum)
  if (freq >= 52 && freq <= 53){calendar <- xts::apply.weekly(Cal, sum)
                                #zoo::index(calendar) <- format(zoo::index(calendar), "%Y-%V")
                                message("It is assumed that this is supposed to be a weekly series")}
  if (freq >= 365 && freq <= 366){calendar <- Cal
                                message("It is assumed that this is supposed to be a daily series")}
  if (!any(freq==12,freq==4,freq >= 52 & freq <= 53, freq >= 365 & freq <= 366)){ calendar <- Cal
                                message("This function can only handle daily, weekly, monthly and quarterly data. The output is now a daily time series. You have to convert it yourself.")}
  
  calendar <- calendar[seq_len(n)]
  
  if (center) {calendar <- .centeruser.xts(calendar)}
  
  if (multiplicative) calendar <- 1+calendar/100
  
  return(calendar)
}



.centeruser.xts <- function(x) {
  .substr_from<- function(x, k) {return(substr(x, k, nchar(x)))}
  yearless_dates <-  .substr_from(as.character(zoo::index(x)),6) # 10x faster than gsub solution
  
  Means <- stats::aggregate(x=x, by=yearless_dates, FUN=mean, na.rm=TRUE)
  
  .subs <- function(y, Means) {
    as.numeric(y) - as.numeric(Means[grep(.substr_from(as.character(zoo::index(y)),6), zoo::index(Means))])}
  
  out <- xts::xts(sapply(seq_len(length(x)), function(y) .subs(x[y], Means=Means)), zoo::index(x))
  
  return(out)
}

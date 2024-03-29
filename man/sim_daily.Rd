% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_daily.R
\name{sim_daily}
\alias{sim_daily}
\title{Simulate a daily seasonal series}
\usage{
sim_daily(
  N,
  sd = 2.5,
  change_sd = 0.05,
  week_sd = NA,
  month_sd = NA,
  year_sd = NA,
  week_change_sd = NA,
  month_change_sd = NA,
  year_change_sd = NA,
  innovations_sd = 1,
  sa_sd = NA,
  model = list(order = c(3, 1, 1), ma = 0.5, ar = c(0.2, -0.4, 0.1)),
  beta_1 = 0.9,
  beta_tau = 0,
  start = c(2020, 1),
  multiplicative = TRUE,
  extra_smooth = FALSE,
  calendar = list(which = "Easter", from = -2, to = 2),
  outlier = NULL,
  timewarping = TRUE,
  as_index = FALSE
)
}
\arguments{
\item{N}{length in years}

\item{sd}{Standard deviation for all seasonal factors}

\item{change_sd}{Standard deviation of simulated change for all seasonal factors}

\item{week_sd}{Standard deviation of the seasonal factor for day-of-the-week}

\item{month_sd}{Standard deviation of the seasonal factor for day-of-the-month}

\item{year_sd}{Standard deviation of the seasonal factor for day-of-the-year}

\item{week_change_sd}{Standard deviation of simulated change to seasonal factor for day-of-the-week}

\item{month_change_sd}{Standard deviation of simulated change to seasonal factor for month-of-the-week}

\item{year_change_sd}{Standard deviation of simulated change to seasonal factor for year-of-the-week}

\item{innovations_sd}{Standard deviation of the innovations used in the non-seasonal regarima model}

\item{sa_sd}{Standard deviation of the non-seasonal time series}

\item{model}{Model for non-seasonal time series. A list.}

\item{beta_1}{Persistance wrt to previous period of the seasonal change}

\item{beta_tau}{Persistance wrt to one year/cycle before of the seasonal change}

\item{start}{Start date of output time series}

\item{multiplicative}{Boolean. Should multiplicative seasonal factors be simulated}

\item{extra_smooth}{Boolean. Should the seasonal factors be smooth on a period-by-period basis}

\item{calendar}{Parameters for calendar effect, a list, see sim_calendar}

\item{outlier}{Parameters for outlier effect, a list, see sim_outlier}

\item{timewarping}{Should timewarping be used to obtain the day-of-the-month factors}

\item{as_index}{Shall series be made to look like an index (i.e. shall values be relative to reference year = second year)}
}
\value{
Multiple simulated daily time series of class xts including:
\describe{
\item{original}{The original series}
\item{seas_adj}{The original series without calendar and seasonal effects}
\item{sfac7}{The day-of-the-week effect}
\item{sfac31}{The day-of-the-month effect}
\item{cfac}{The calendar effects}
\item{outlier}{The outlier effects}
}
}
\description{
Simulate a daily seasonal series as described in Ollech (2021).
}
\details{
Standard deviation of the seasonal factor is in percent if a multiplicative time series model is assumed. Otherwise it is in unitless.
Using a non-seasonal ARIMA model for the initialization of the seasonal factor does not impact the seasonality of the time series. It can just make it easier for human eyes to grasp the seasonal nature of the series. The definition of the ar and ma parameter needs to be inline with the chosen model.
If only change_sd is specified, the change parameters for the single seasonal factors are set individually as change_sd/365*(length of seasonal cycle)
The parameters that can be set for calendar and outlier are those defined in sim_outlier and sim_calendar.
}
\examples{
x=sim_daily(5, multiplicative=TRUE, outlier=list(k=5, type=c("AO", "LS"), effect_size=50))
ts.plot(x[,1])
}
\references{
Ollech, D. (2021). Seasonal adjustment of daily time series. Journal of Time Series Econometrics. \doi{10.1515/jtse-2020-0028}
}
\author{
Daniel Ollech
}

## percent_outbound.R
##


percout <- function(time_series, windows) {
  p_out <- sum(time_series < windows[,1] | time_series > windows[,2]) / length(time_series)
  return(p_out)
}


# A function to generate a periodic pulse
# A = amplitude
# delta = shift in location of pulse
# P = periodicity
# Pu = pulse width
# m0 = intercept
periodic_pulse_predict <- function(
    A, delta, P, Pu, time_steps, m0
){
  pulse <- matrix(0, ncol = 500, nrow = length(time_steps))
  
  for(n in 1:500){
    pulse[,n] <- (2*A/(n*pi))*sin((pi*n*Pu)/P)*
      cos((2 *pi*n * (time_steps - delta))/(P))
  }
  
  ans <- m0 + rowSums(pulse)
  return(ans)
}

time <- 1:200
pulse_vals <- periodic_pulse_predict(
  A = 60,
  delta = 25,
  P = 50,
  Pu = 28,
  time_steps = time,
  m0 = 34
)

pulse_vals <- floor(pulse_vals)
plot(pulse_vals ~ time, type = "o")

library(bbplot)


# get a sequence of repeating deployments
pulse_vals <- rep(c(0,60, 0), times = c(12,28, 12))
pulse_vals <- rep(pulse_vals, 4)

# want nice square edges
time <- rep(NA, length(pulse_vals))
time[1] <- 1
to_sub <- 0
for(i in 2:length(time)){
  if(pulse_vals[i] == 60 & pulse_vals[i-1] == 0){
    to_sub <- to_sub - 1
  }
  if(pulse_vals[i] == 0 & pulse_vals[i-1] == 60){
    to_sub <- to_sub - 1
  }
  time[i] <- i + to_sub
}


one_time <- pulse_vals
my_60 <- which(one_time == 60)
one_time[my_60[1]:my_60[length(my_60)]] <- 60

sum(one_time == 60) - (28*4)
one_time[13:(13+39)] <- 0
one_time[(196-39):196] <- 0
# add more leading and trailing zeroes to better fit 4 distinct sampling
#  periods.


plot(one_time ~ time, type = "l")

windows(4,4)

svg("./plots/temp_pulse.svg", width = 4, height = 4)
par(mar = c(2,2,0.5,0.5), lend = 2)
bbplot::blank(xlim = c(0,200), ylim = c(0,60), bty = "l")
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(
  text = c(0,100,200),at = c(0,100,200),
  side = 1, line = 0.8, cex = 1.4)
bbplot::axis_text(
  text = c(0,30,60), at = c(0,30,60),
  side = 2, line = 0.6, cex = 1.4, las = 1)
lines(
  x = time,
  y = pulse_vals,
  lwd = 3
)

dev.off()
svg("./plots/temp_long_sample.svg", width = 4, height = 4)
par(mar = c(2,2,0.5,0.5), lend = 2)
bbplot::blank(xlim = c(0,200), ylim = c(0,60), bty = "l")
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(
  text = c(0,100,200),at = c(0,100,200),
  side = 1, line = 0.8, cex = 1.4)
bbplot::axis_text(
  text = c(0,30,60), at = c(0,30,60),
  side = 2, line = 0.6, cex = 1.4, las = 1)
lines(
  x = time,
  y = one_time,
  lwd = 3
)

dev.off()







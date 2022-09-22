rm(list = ls())
# Demonstrate the superposition of sin/consine waves
t = 1:100
x1 = 2*cos(2*pi*t*6/100) + 3*sin(2*pi*t*6/100)
x2 = 4*cos(2*pi*t*10/100) + 5*sin(2*pi*t*10/100)
x3 = 6*cos(2*pi*t*40/100) + 7*sin(2*pi*t*40/100)

x = x1 + x2 + x3


plot.ts(x1, ylim = c(-16,16), main = 'frequency = 6/100, amp^2 = 13')
plot.ts(x2, ylim = c(-16,16), main = 'frequency = 10/100, amp^2 = 41')
plot.ts(x3, ylim = c(-16,16), main = 'frequency = 40/100, amp^2 = 85')
plot.ts(x, ylim = c(-16,16), main = "Simulated stationary time series")

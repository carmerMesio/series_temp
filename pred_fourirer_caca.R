install.packages("TSA")
library(TSA)

# compute the Fourier Transform
p = periodogram(superm$`General [kWh]`)
dd = data.frame(freq=p$freq, spec=p$spec)
order = dd[order(-dd$spec),]
top2 = head(order, 2)

# display the 2 highest "power" frequencies
top2
time = 1/top2$f
time

require(readr)
require(lubridate)
require(mgcv)
require(sfsmisc)
require(xgboost)
require(doRNG)
require(doSNOW)
source("multi_spec.R")
source("fitpred.R")

mitvo <- read_csv("data/Metro_Interstate_Traffic_Volume.csv")
mitv <- mitvo[,c("temp", "rain_1h", "clouds_all", "date_time", "traffic_volume")]
# remove some outliers
wr <- which(mitv$rain_1h > 100)
wt <- which(mitv$temp == 0)
mitv <- mitv[-c(wr, wt), ]
# remove duplicate observations
mitv <- aggregate(mitv, list(mitv$date_time), mean)
# day of the year
mitv$doy = yday(as.Date(mitv$date_time))
# day of the week
mitv$dow = wday(as.Date(mitv$date_time))
mitv$hour <- as.numeric(format(mitv$date_time, "%H"))
mitv$year <- as.numeric(format(mitv$date_time, "%Y"))
# encode target variable
mitv$y <- mitv$traffic_volume
mitv <- mitv[,c("y", "temp", "rain_1h", "clouds_all", "doy", "dow", "hour", "year")]

# find time difference between observation i and i-1
mitv.diff <- (mitv[-1, ] -mitv[-nrow(mitv),])[, c("doy", "hour", "year")]

good <- list(c(0, 1, 0), c(1, -23, 0), c(-365, -23, 1), c(-364, -23, 1))
ind <- numeric()
for (go in good)
  ind <- c(ind, which(apply(t(mitv.diff) == go, 2, sum) == 3))
ind <- sort(ind)

mitv$y.pre <-c(NA, mitv$y[-nrow(mitv)])

# only keep observations where on hour before information is available
mitv <- mitv[ind + 1,]

if(!exists("analysis"))
  analysis <- list()

cols <- c("y", "hour", "dow", "temp", "rain_1h", "clouds_all")
name <- paste(cols[-1], collapse = "+")
analysis[[name]] <- multi.spec(mitv[,cols], B = 25, fitting = fitxg, predicting = predxg, return.predictor = TRUE,
                 return.residual = TRUE, return.indices = TRUE,
                 parallel = TRUE, sockets = 5)

# ms <- multi.spec(data.frame(mitv), B = 25, fitting = fitxg, predicting = predxg, return.predictor = TRUE)

save(analysis, file = paste("data/analysis ", format(Sys.time(), "%d-%b-%Y %H.%M"), ".RData", sep = ""))


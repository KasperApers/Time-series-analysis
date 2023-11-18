rm(list = ls())   # Clear workspace

setwd("C:\\Users\\apers\\OneDrive\\Documenten\\R\\time series")
output= "C:\\Users\\apers\\OneDrive\\Documenten\\R\\time series"

if(!require(urca)){install.packages("urca")}
if(!require(tools)){install.packages("tools")}
if(!require(stargazer)){install.packages("stargazer")}
if(!require(dynlm)){install.packages("forecast")}
if(!require(Hmisc)){install.packages("Hmisc")}
if(!require(dynlm)){install.packages("dynlm")}
if(!require(lmtest)){install.packages("lmtest")}
if(!require(vars)){install.packages("vars")}
if(!require(VAR.etp)){install.packages("VAR.etp")}
if(!require(forecast)){install.packages("forecast")}
if(!require(TSstudio)){install.packages("TSstudio")}

library(urca)            
library(stargazer)
library(tools)
library(Hmisc)
library(dynlm)
library(lmtest)  
library(vars)
library(VAR.etp)
library(forecast)
library(TSstudio)

source("TimeSeriesFunctions_v3.R") # includes a number of additional functions not available online (download from Ufora and save in your input directory)

data = read.table("data.csv", header=TRUE, sep=",")
View(data)

#Q1
# Create a time series object for the US_IP
ts_lnYt00 <- ts(data$US_IP, start=c(1974, 1), end=c(2017, 12), frequency=12)
diff_ts_lnYt00 <- diff(ts_lnYt00) # Calculate the differences
par(mfrow=c(1, 1))
# Plot ln Yt
plot.ts(ts_lnYt00,main="Log of U.S. Industrial Production Over Time", ylab="ln Yt", xlab="Time", col="blue", lwd=2)
grid()
# Plot ∆ ln Yt
plot.ts(diff_ts_lnYt00,main="Change in U.S. Industrial Production Over Time", ylab="Delta ln Yt", xlab="Time")
abline(h = 0, col = "blue", lty = 2)

#Q2 ACF PACF for ΔlnYt
ts_lnYt <- ts(data$US_IP, start=c(1974, 1), end=c(2015, 12), frequency=12)
diff_ts_lnYt <- diff(ts_lnYt)

#ACF level: MA(6)
diff_ts_lnYt_acf <- acf(diff_ts_lnYt, main="ACF of ΔlnYt for 1974.1-2015.12",ylab = NULL, xlab = "Lag",lwd=2, col="blue")
plot(diff_ts_lnYt_acf, main="ACF of ΔlnYt for 1974.1-2015.12",
     ylab=NULL, bty="l", xlab="lag", lwd=7.5, mar=c(2,2,0,0), col="darkblue", col.axis="black")

#PACF level:AR(3)
diff_ts_lnYt_pacf <- pacf(diff_ts_lnYt, main="PACF of ΔlnYt for 1974.1-2015.12",ylab = NULL, xlab = "Lag",lwd=2, col="red")
plot(diff_ts_lnYt_pacf, main="PACF of ΔlnYt for 1974.1-2015.12",
     ylab=NULL, bty="l",xlab="lag", lwd=7.5, mar=c(2,2,0,0), col="darkblue", col.axis="black")



#Q3 Construct and estimate an ARMA model for ∆lnYt


#AR
AIC_matrix_AR <- aic_table_css(diff_ts_lnYt, 10, 0)
AIC_matrix_AR
min(AIC_matrix_AR) #AR(5)

BIC_matrix_AR <- bic_table_css(diff_ts_lnYt, 10, 0)
BIC_matrix_AR
min(BIC_matrix_AR) #AR(3)

AR3 = arima(diff_ts_lnYt, order = c(3,0,0), method = "CSS-ML")
AR3
AR3BIC = round(AIC(AR3,k = log(length(diff_ts_lnYt))),digits=4)
AR3BIC

AR4 = arima(diff_ts_lnYt, order = c(4,0,0), method = "CSS-ML")
AR4
AR4BIC = round(AIC(AR4,k = log(length(diff_ts_lnYt))),digits=4)
AR4BIC

ARMA41 = arima(diff_ts_lnYt, order = c(4,0,1), method = "CSS-ML")
ARMA41

#Q4 Residuals.
resARMA41 <- ARMA41$residuals
resARMA41_ACF <- acf(resARMA41, lag.max = 50)
resARMA41_PACF <- pacf(resARMA41, lag.max = 50)
plot(resARMA41_ACF, main="ACF of residuals of ARMA(4,1).",
     ylab=NULL, bty="l", xlab="lag", lwd=7.5, mar=c(2,2,0,0), col="darkblue", col.axis="black")
plot(resARMA41_PACF, main="PACF of residuals of ARMA(4,1).",
     ylab=NULL, bty="l", xlab="lag", lwd=7.5, mar=c(2,2,0,0), col="darkblue", col.axis="black")
QstatresresARMA41 = LjungBox(resARMA41,5,0)
stargazer(QstatresresARMA41, type = "text", summary = FALSE) #According to the p values no autocorr in the error terms. (white noise)


#Q5 Stability.
## AR4
res_AR4 <- residuals(AR4)-mean(residuals(AR4))
sigma_AR4 <- sqrt(AR4$sigma2)
cumsum_AR4 <- cumsum(res_AR4) / (sqrt(length(diff_ts_lnYt))*sigma_AR4)
cumsum_AR4.ts = ts(cumsum_AR4,start = c(1974, 1), end = c(2015, 12), frequency = 12)

par(mar=c(3,3,3,3))
plot(cumsum_AR4.ts, ylab = NULL, xlab = NULL,lwd = 3, main= "Cusum AR(4).", col="darkblue", bty = 'l', col.axis = "black", type = "l",
     ylim = c(min(cumsum_AR4.ts, -1.5), max(cumsum_AR4.ts,1.5))) #doesn't step outside the bounds: model should be stable
abline(h = c(-1.36,1.36), col = "darkgrey", lty = 2,lwd = 2)
abline(h = c(0), col = "black", lwd = 2)

#Q6 construct out-of-sample forecasts
noholdout1 <- window(diff_ts_lnYt00 ,start=c(1974,1),end=c(2015,12))
holdout1 <- window(diff_ts_lnYt00,start=c(2016,1),end=c(2016,12))
noholdout2 <- window(diff_ts_lnYt00,start=c(1974,1),end=c(2016,12))
holdout2 <- window(diff_ts_lnYt00,start=c(2017,1),end=c(2017,12))

#forecast for 2016
AR4_noholdout1 = arima(noholdout1, order=c(4,0,0), method = "CSS-ML", include.mean = TRUE)
Fcast_AR4a <- forecast(AR4_noholdout1,h=12)
summary(Fcast_AR4a)
plot_forecast(Fcast_AR4a, title = "Forecast for 2016.", Xtitle = "Year")
accuracy(Fcast_AR4a,holdout1)


#forecast for 2017
AR4_noholdout2 = arima(noholdout2, order = c(4,0,0), method = "CSS-ML", include.mean = TRUE)
Fcast_AR4b <- forecast(AR4_noholdout2,h=12)
summary(Fcast_AR4b)
plot_forecast(Fcast_AR4b, title = "Forecast for 2017.", Xtitle = "Year")
accuracy(Fcast_AR4b, holdout2)


#Q7 diff_lnYt stationary?
## ADF test
### With trend and drift.
ADF_trend_drift = ur.df(diff_ts_lnYt00, lags = 3, type = "trend")
summary(ADF_trend_drift)

  ## tau  : gamma = 0 (Reject)
  ## phi3 : gamma = phi = 0 (Reject)
  ## phi2 : alpha = gamma = phi = 0 (Reject)

### Test with only drift.
ADF_drift = ur.df(diff_ts_lnYt00, lags = 3, type = "drift")
summary(ADF_drift)

  ## tau  : gamma = 0 (Reject)
  ## phi1 : alpha = gamma  = 0 (Reject)

### Test without trend or drift.
ADF_test = ur.df(diff_ts_lnYt00, lags = 3, type = "none")
summary(ADF_test)

  ## tau  : gamma = 0 (Reject)

#Q8 lnYt stationary?
## ADF test
### With trend and drift.
ADF_trend_drift_2 = ur.df(ts_lnYt00, lags = 3, type = "trend")
summary(ADF_trend_drift_2)

## tau  : gamma = 0 (Cant Reject)
## phi3 : gamma = phi = 0 (Cant Reject)
## phi2 : alpha = gamma = phi = 0 (Cant Reject)

### Test with only drift.
ADF_drift_2 = ur.df(ts_lnYt00, lags = 3, type = "drift")
summary(ADF_drift_2)

## tau  : gamma = 0 (Cant Reject)
## phi1 : alpha = gamma  = 0 (Cant Reject)

### Test without trend or drift.
ADF_test_2 = ur.df(ts_lnYt00, lags = 3, type = "none")
summary(ADF_test_2)

## tau  : gamma = 0 (Cant Reject)




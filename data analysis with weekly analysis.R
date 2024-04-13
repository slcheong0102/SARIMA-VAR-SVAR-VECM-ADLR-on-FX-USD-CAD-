# ----- Setup -----
rm(list=ls());
library(readxl);
library(logr);
library(tidyverse);
library(dplyr);
library(ggplot2);
library(lmtest);
library(forecast);
library(tseries);
library(seasonal);
library(vars);
library(urca);
library(tsDyn);
library(TSstudio);
library(lmtest);
library(mFilter);
library(sensitivity);
library(dynamac);
library(strucchange);

# Import forex data
setwd("C:/Users/soulu/OneDrive - University of Toronto/University of Toronto/UTM 2021-2023/8_2024 Winter/ECO475H5S/Term Papaer/Data Analysis");
forex <- read_excel("FOREX.xlsx");
colnames(forex)[colnames(forex) == "week"] <- "date";
forex$pct_change_us_cpi <- replace_na(forex$pct_change_us_cpi, 0);
forex$pct_change_ca_cpi <- replace_na(forex$pct_change_ca_cpi, 0);
forex$pct_change_us_trade <- replace_na(forex$pct_change_us_trade, 0);
forex$pct_change_ca_trade <- replace_na(forex$pct_change_ca_trade, 0);
forex$us_unemployment_rate <- replace_na(forex$us_unemployment_rate, 0);
forex$ca_unemployment_rate <- replace_na(forex$ca_unemployment_rate, 0);

# Define time series by week of year
date <- format(forex$date, "%Y-%m-%d");
wofy <- format(forex$date, "%Y-W%V");
year <- as.numeric(substr(forex$date,1,4));
week <- as.numeric(substr(wofy,7,8));
month <- as.numeric(substr(forex$date,6,7));

# Define variables
forex$rate_differential <- forex$us_rate - forex$ca_rate;
forex$year_month <- format(forex$date,"%Y-%m");
forex$pct_change_cpi_differential <- forex$pct_change_us_cpi - forex$pct_change_ca_cpi;
forex$pct_change_trade_differential <- forex$pct_change_us_trade - forex$pct_change_ca_trade;
forex$unemployment_differential <- forex$us_unemployment_rate - forex$ca_unemployment_rate;

# Create forex_month dataframe
forex_month <- forex %>% group_by(year_month) %>% summarise(price.month=mean(price),
                                                            rate.differential.month=mean(rate_differential),                                             
                                                            cpi.differential.month=mean(pct_change_cpi_differential),
                                                            trade.differential.month=mean(pct_change_trade_differential),
                                                            unemployment.differential.month=mean(unemployment_differential)
);
forex_month$unemployment.differential.month[176:187] <- NA; #no obs for first and last month
forex_month$unemployment.differential.month[1:55] <- NA;
forex_month <- na.omit(forex_month);
year2 <- as.numeric(substr(forex_month$year_month,1,4));
month2 <- as.numeric(substr(forex_month$year_month,6,7));
time=seq(1,length(forex$date),by=1);
month=seq(1,length(forex_month$year_month),by=1);

# Create forex_month_VAR for running the VAR/VECM model
forex_month_VAR <- forex_month;
forex_month_VAR$year_month <- NULL;

## ----- Basic Summary Statistics -----
summary(forex);
summary(forex_month);
cor(forex_month_VAR);


## ----- 1.1 ARIMA of Average Forex Closing Price (forex$price) -----

############################ Weekly Analysis
# Plot forex$price
price.ts=ts(forex$price, start=c(year[1],week[1]), frequency=52);
ts_plot(price.ts, title="Time Series Plot of Weekly Average Fx Closing Price", Ytitle="Fx Price", Xtitle="Time");
plot(price.ts, main="Time Series Plot of Weekly Average Fx Closing Price", ylab="Fx Price", xlab="Time");
adf.test(price.ts);
acf(price.ts,lag=100);
pacf(price.ts);

## Take first difference of forex$price
price.diff <- diff(unlist(forex$price));
price.diff.ts=ts(price.diff, start=c(year[2],week[2]), frequency=52);
plot(price.diff.ts, main="Time Series Plot of First-Differenced Fx Closing Price", ylab="First-Diff Price")
abline(h=mean(price.diff.ts),lty=2,col='red');
ts_plot(price.diff.ts, title="Time Series Plot of First-Differenced Fx Closing Price", Ytitle="First-Diff Fx Price", Xtitle="Time");
adf.test(price.diff);
acf(price.diff);
pacf(price.diff);

## Fit ARIMA model and forecast
price.diff.fit <- auto.arima(price.diff);
price.diff.fit;
checkresiduals(price.diff.fit);
check_res(price.diff.fit);
autoplot(forecast(price.diff.fit));
forecast(price.diff.fit);
# The first-difference process is MA(2)

## AR decomposition on price.ts using Yule Walker
price.ar <- ar.yw(price.ts);
price.ar;
plot.ts(price.ar$aic, main="Plot of Akaike's Criterion", xlab="Order", ylab="Akaike's Criterion (AIC)");
price.ar$aic;
plot.ts(price.ar$resid, main="Plot of Residuals from AR Model", xlab="Time", ylab="Residuals from AR Model")
abline(h=mean(price.ar$resid[3:810]),lty=2, col='red');
acf(price.ar$resid[3:810]);
price.ar$ar;
price.ar$var.pred;
price.ar$asy.var.coef;
checkresiduals(price.ar);
autoplot(forecast(price.ar));
forecast(price.ar);
# Forex prices can also be fitted with an AR(2) model

######################## Monthly Analysis
price.month.ts=ts(forex_month$price.month, start=c(year2[1],month2[1]), frequency=12);
plot(price.month.ts, main="Time Series Plot of Monthly Average Fx Closing Price", ylab="Fx Price", xlab="Time");
ts_plot(price.month.ts, title="Time Series Plot of Monthly Average Fx Closing Price", Ytitle="Fx Price", Xtitle="Time");
adf.test(price.month.ts);
acf(price.month.ts);
pacf(price.month.ts);

## Fit ARIMA Model and Forecast
fit1.1 <- auto.arima(price.month.ts);
summary(fit1.1);
autoplot(fit1.1);
checkresiduals(fit1.1);
autoplot(forecast(fit1.1));
forecast(fit1.1);

## Take the first difference
price.diff.month.ts <- diff(price.month.ts);
adf.test(price.diff.month.ts);
plot(price.diff.month.ts, main="Time Series Plot of First-Differenced Monthly Fx Closing Price", ylab="Fx Price", xlab="Time");
abline(h=mean(price.diff.month.ts),lty=2,col='red');
acf(price.diff.month.ts);
pacf(price.diff.month.ts);
fit1.1 <- auto.arima(price.diff.month.ts);
summary(fit1.1);

## Spectral Analysis/STL and ETS Model (Advanced Method)
price.month.z <- spec.pgram(price.month.ts, fast=FALSE, taper=0.0);
stl1.1 <- stl(price.month.ts, s.window="periodic");
summary(stl1.1);
autoplot(stl1.1);
stl1.1.seasonal <- stl1.1$time.series[1:185];
stl1.1.trend <- stl1.1$time.series[186:370];
stl1.1.remainder <- stl1.1$time.series[371:555];
plot.ts(stl1.1.seasonal);
plot.ts(stl1.1.trend);
plot.ts(stl1.1.remainder)
abline(h=mean(stl1.1.remainder),lty=2,col='red');
adf.test(stl1.1.remainder);
fit1.1 <- stl1.1.seasonal + stl1.1.trend;
autoplot(cbind(price.month.ts,fit1.1));
autoplot(stlf(price.month.ts));
stlf(price.month.ts);
# Find expression for seasonal component using fourier decomposition/spectral analysis
fit1.1 <- lm(stl1.1.seasonal~sin(2*pi*month/12)+cos(2*pi*month/12)
             +sin(4*pi*month/12)+cos(4*pi*month/12)
             +sin(6*pi*month/12)+cos(6*pi*month/12));
summary(fit1.1);
stl1.1.fit <- stl1.1.trend+fit1.1$fitted.values;
autoplot(cbind(price.month.ts,stl1.1.fit));

## Verification using In-Sample Forecasting and Diagnostics (Construct training and testing set, i.e. train-test split)
adf.test(price.month.ts); # Non-stationary --> needs to take difference (verified that first difference is stationary)
adf.test(diff(price.month.ts)); # Stationary --> price.month is first-difference stationary, i.e. price.month.ts is I(1)
split <- ts_split(price.month.ts, sample.out=12);
price.train <- split$train;
price.test <- split$test;

# Use an ARIMA Diagnostic Plot on the training set
arima_diag(price.train);
# Based on the diagnostic plot, we guess that our data follow SARIMA Model

# Build a SARIMA(1,1,1)(1,0,0)[12] model with training data set
sarima111.100 <- arima(price.train, order=c(1,1,1), seasonal=c(1,0,0));
autoplot(sarima111.100);
check_res(sarima111.100);

# Build a SARIMA(2,1,2)(1,0,0)[12] model with training data set
sarima212.100 <- arima(price.train, order=c(2,1,2), seasonal=c(1,0,0));
autoplot(sarima212.100);
check_res(sarima212.100);

# Build a SARIMA(2,1,2)(1,0,1)[12] model with training data set (suggested by auto.arima(price.train))
auto.arima(price.train); # Use r code to see which model the computer thinks is best fit
sarima212.101 <- arima(price.train, order=c(2,1,2), seasonal=c(1,0,1));
autoplot(sarima212.101);
check_res(sarima212.101);

# Model selection based on forecast values and diagnostics
sarima111.100.fit <- forecast(sarima111.100, h=12); # Forecast 12 periods ahead (since we split the price.test to be the last 12 observations)
test_forecast(actual=price.month.ts, forecast.obj=sarima111.100.fit, test=price.test);
accuracy(sarima111.100.fit, price.test);

sarima212.100.fit <- forecast(sarima212.100, h=12);
test_forecast(actual=price.month.ts, forecast.obj=sarima212.100.fit, test=price.test);
accuracy(sarima212.100.fit, price.test);

sarima212.101.fit <- forecast(sarima212.101, h=12);
test_forecast(actual=price.month.ts, forecast.obj=sarima212.101.fit, test=price.test);
accuracy(sarima212.101.fit, price.test);

# Define models with arguments (Construct train-validation-test split)
methods <- list(model1=list(method="arima",
                            method_arg=list(order=c(1,1,1), seasonal=list(order=c(1,0,0))),
                            notes="SARIMA(1,1,1)(1,0,0)"),
                model2=list(method="arima",
                            method_arg=list(order=c(2,1,2), seasonal=list(order=c(1,0,0))),
                            notes="SARIMA(2,1,2)(1,0,0)"),
                model3=list(method="arima",
                            method_arg=list(order=c(2,1,2), seasonal=list(order=c(1,0,1))),
                            notes="SARIMA(2,1,2)(1,0,1)"));

# Train the models using backtesting
obj <- train_model(input=price.month.ts, methods=methods, 
                   train_method=list(partitions=3, sample.out=12, space=3),
                   horizon=12, error="RMSE");

# Plot the model
plot_model(model.obj=obj);
# Seems like SARIMA(1,1,0)(0,0,1)[12] is the best model considering prediction power and model complexity

# Generate the optimal fit based on our training data SARIMA(1,1,0)(0,0,1)[12]
fit1.1 <- auto.arima(price.month.ts, seasonal=TRUE);
summary(fit1.1);
autoplot(fit1.1);
check_res(fit1.1);

# Generate final forecast based on our whole dataset (for period t+1)
forecast <- forecast(price.month.ts, model=fit1.1, h=12);
plot_forecast(forecast);
summary(forecast);


## ----- 1.2 ARIMA of Treasury Bill Rate Differential (forex$rate_differential) -----

########################### Weekly Analysis 
# Plot rate.differential.ts
rate.differential.ts=ts(forex$rate_differential, start=c(year[1],week[1]), frequency=52);
ts_plot(rate.differential.ts);
mod1.2 <- lm(forex$rate_differential~time);
intercept1.2 <- mod1.2$coefficients[1];
slope1.2 <- mod1.2$coefficients[2];
y1.2 <- intercept1.2+slope1.2*length(time);
plot(rate.differential.ts, main="Time Series Plot of TBill Rate Differential", ylab="Rate Differential", xlab="Time")
lines(c(time(rate.differential.ts)[1],time(rate.differential.ts)[810]), c(intercept1.2,y1.2),lty=2,col='red');
# TBill rate differential "seems" to be non-stationary and follows a positive linear trend

# Augmented Dickey Fuller test of rate.differential.ts
adf.test(rate.differential.ts);
# Since the p-value is 0.5777 > 0.05, we fail to reject H0: non-stationarity, i.e. rate.differential.ts is a unit root process.

# ACF and PACF of rate.differential.ts
acf(rate.differential.ts,100);
pacf(rate.differential.ts);

# STL and ETS method of Model Decomposition
rate.differential.z <- spec.pgram(rate.differential.ts, fast=FALSE, taper=0.0);
stl1.2 <- stl(rate.differential.ts, s.window="periodic");
summary(stl1.2);
autoplot(stl1.2);
stl1.2.seasonal <- stl1.2$time.series[1:810];
stl1.2.trend <- stl1.2$time.series[811:1620];
stl1.2.remainder <- stl1.2$time.series[1621:2430];
plot.ts(stl1.2.seasonal);
plot.ts(stl1.2.trend);
plot.ts(stl1.2.remainder);
abline(h=mean(stl1.2.remainder),lty=2,col='red');
adf.test(stl1.2.remainder);
fit1.2 <- stl1.2.seasonal + stl1.2.trend;
autoplot(cbind(rate.differential.ts,fit1.2));
autoplot(stlf(rate.differential.ts));
stlf(price.month.ts);

# Take first difference of rate.differential.ts
rate.differential.diff.ts <- diff(rate.differential.ts);

# Plot the first difference
plot(rate.differential.diff.ts, main="Time Series Plot of First-Differenced TBill Rate Differential", ylab="First-Diff TBill Rate")
abline(h=mean(rate.differential.diff.ts),lty=2,col='red');
ts_plot(rate.differential.diff.ts, title="Time Series Plot of First-Differenced TBill Rate Differential", Ytitle="First-Diff TBill Rate", Xtitle="Time");

# Augmented Dickey Fuller test of first difference
adf.test(rate.differential.diff.ts);
# The p-value is less than 0.01, we can reject H0: non-stationarity, i.e. the first-difference is a sufficient model.

# ACF and PACF of first difference
acf(rate.differential.diff.ts);
pacf(rate.differential.diff.ts);

# Find expression for seasonal component using fourier decomposition/spectral analysis
fit1.2 <- lm(stl1.2.seasonal~sin(2*pi*time/52)+cos(2*pi*time/52)
             +sin(4*pi*time/52)+cos(4*pi*time/52)
             +sin(6*pi*time/52)+cos(6*pi*time/52)
             +sin(8*pi*time/52)+cos(8*pi*time/52)
             +sin(10*pi*time/52)+cos(10*pi*time/52));

summary(fit1.2);
stl1.2.fit <- stl1.2.trend+fit1.2$fitted.values;
autoplot(cbind(rate.differential.ts,stl1.2.fit));

# Find an expression for both rate.differential.ts and rate.differential.diff.ts
auto.arima(rate.differential.ts, method="CSS");
auto.arima(rate.differential.diff.ts, method="CSS"); # Note that "ML" method here has an extremely long run time

## AR decomposition on rate.differential.ts using Yule Walker
plot(rate.differential.ts);
plot(rate.differential.diff.ts);
rate.differential.ar <- ar.yw(rate.differential.ts);
rate.differential.ar;
names(rate.differential.ar);
plot.ts(rate.differential.ar$aic, main="Plot of Akaike's Criterion", xlab="Order", ylab="Akaike's Criterion (AIC)");
rate.differential.ar$aic;
plot.ts(rate.differential.ar$resid, main="Plot of Residuals from AR Model", xlab="Time", ylab="Residuals from AR Model")
abline(h=mean(rate.differential.ar$resid[8:810]),lty=2, col='red');
acf(rate.differential.ar$resid[8:810]);
rate.differential.ar$ar;
rate.differential.ar$var.pred;
rate.differential.ar$asy.var.coef;
checkresiduals(rate.differential.ar);
autoplot(forecast(rate.differential.ar));


############################ Monthly Analysis 

## Fit an ARIMA model
rate.differential.month.ts=ts(forex_month$rate.differential.month, start=c(year2[1],month2[1]), frequency=12);
adf.test(rate.differential.month.ts);
adf.test(diff(rate.differential.month.ts)); # rate.differential.month.ts is I(1)
acf(rate.differential.month.ts,lag.max = 50);
pacf(rate.differential.month.ts,lag.max = 50);
auto.arima(rate.differential.month.ts);

## Take the first difference
rate.differential.diff.month.ts <- diff(rate.differential.month.ts);
adf.test(rate.differential.diff.month.ts);
acf(rate.differential.diff.month.ts);
pacf(rate.differential.diff.month.ts);
fit1.2 <- auto.arima(rate.differential.diff.month.ts);
summary(fit1.2);
checkresiduals(fit1.2);
autoplot(fit1.2);
autoplot(forecast(fit1.2));

## Spectral Analysis/STL and ETS Method (Advanced Method)
stl1.2 <- stl(rate.differential.month.ts,s.window="periodic");
summary(stl1.2);
autoplot(stl1.2);
autoplot(forecast(stl1.2));
forecast(stl1.2);
# Perform fourier analysis on the seasonal trend
stl1.2.seasonal <- stl1.2$time.series[1:185];
stl1.2.trend <- stl1.2$time.series[186:370];
stl1.2.remainder <- stl1.2$time.series[371:555];
fit1.2 <- lm(stl1.2.seasonal~sin(2*pi*month/12)+cos(2*pi*month/12)
             +sin(4*pi*month/12)+cos(4*pi*month/12)
             +sin(6*pi*month/12)+cos(6*pi*month/12));
summary(fit1.2);
stl1.2.fit <- stl1.2.trend+fit1.2$fitted.values;
autoplot(cbind(rate.differential.month.ts,stl1.2.fit));

## In-Sample Forecasting and Diagnostics (Construct train-test split)
adf.test(rate.differential.month.ts);
split <- ts_split(rate.differential.month.ts, sample.out=12);
rate.train <- split$train;
rate.test <- split$test;

# Use an ARIMA Diagnostic Plot on the training set
arima_diag(rate.train);
# Based on the diagnostic plot, we guess that our data follow SARIMA Model

# Build a SARIMA(1,1,1)(1,0,1)[12] model with training data set
sarima111.101 <- arima(rate.train, order=c(1,1,1), seasonal=c(1,0,1));
autoplot(sarima111.101);
check_res(sarima111.101);

# Build a SARIMA(2,1,1)(1,0,1)[12] model with training data set
sarima211.101 <- arima(rate.train, order=c(2,1,1), seasonal=c(1,0,1));
autoplot(sarima211.101);
check_res(sarima211.101);

# Build a SARIMA(2,1,0)(2,0,0)[12] model with training data set (suggested by auto.arima(rate.train))
auto.arima(rate.train);
sarima210.200 <- arima(rate.train, order=c(2,1,0), seasonal=c(2,0,0));
autoplot(sarima210.200);
check_res(sarima210.200);

# Build a SARIMA(2,1,1)(2,0,1)[12] model with training data set
auto.arima(rate.train);
sarima211.201 <- arima(rate.train, order=c(2,1,1), seasonal=c(2,0,1));
autoplot(sarima211.201);
check_res(sarima211.201);

# Model selection based on forecast values and diagnostics
sarima111.101.fit <- forecast(sarima111.101, h=12);
test_forecast(actual=rate.differential.month.ts, forecast.obj=sarima111.101.fit, test=rate.test);
accuracy(sarima111.101.fit, rate.test);

sarima211.101.fit <- forecast(sarima211.101, h=12);
test_forecast(actual=rate.differential.month.ts, forecast.obj=sarima211.101.fit, test=rate.test);
accuracy(sarima211.101.fit, rate.test);

sarima210.200.fit <- forecast(sarima210.200, h=12);
test_forecast(actual=rate.differential.month.ts, forecast.obj=sarima210.200.fit, test=rate.test);
accuracy(sarima210.200.fit, rate.test);

sarima211.201.fit <- forecast(sarima211.201, h=12);
test_forecast(actual=rate.differential.month.ts, forecast.obj=sarima211.201.fit, test=rate.test);
accuracy(sarima211.201.fit, rate.test);

# Define models with arguments (Construct train-validation-test split)
methods <- list(model1=list(method="arima",
                            method_arg=list(order=c(1,1,1), seasonal=list(order=c(1,0,1))),
                            notes="SARIMA(1,1,1)(1,0,1)"),
                model2=list(method="arima",
                            method_arg=list(order=c(2,1,1), seasonal=list(order=c(1,0,1))),
                            notes="SARIMA(2,1,1)(1,0,1)"),
                model3=list(method="arima",
                            method_arg=list(order=c(2,1,0), seasonal=list(order=c(2,0,0))),
                            notes="SARIMA(2,1,0)(2,0,0)"),
                model4=list(method="arima",
                            method_arg=list(order=c(2,1,0), seasonal=list(order=c(2,0,1))),
                            notes="SARIMA(2,1,0)(2,0,1)")
);


# Train the models using backtesting
obj <- train_model(input=rate.differential.month.ts, methods=methods, 
                   train_method=list(partitions=3, sample.out=12, space=3),
                   horizon=12, error="RMSE");

# Plot the model
plot_model(model.obj=obj);
# Seems like SARIMA(2,1,1)(1,0,1)[12] is the best model considering prediction power and model complexity (using mape and rmse)

# Generate the optimal fit based on our training data SARIMA(2,1,1)(1,0,1)[12]
fit1.2 <- arima(rate.differential.month.ts, order=c(2,1,1), seasonal=c(1,0,1));
summary(fit1.2);
autoplot(fit1.2);
check_res(fit1.2);

# Generate final forecast based on our whole dataset (for period t+1)
forecast <- forecast(rate.differential.month.ts, model=fit1.2, h=12);
plot_forecast(forecast);
summary(forecast);


## ----- 1.3 ARIMA of CPI Percentage Change Differential (forex_month$cpi.differential.month) -----
cpi.differential.month.ts=ts(forex_month$cpi.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(cpi.differential.month.ts, main="Time Series Plot of Monthly CPI Pct Change Differential", ylab="CPI Pct Change Differential", xlab="Time");
ts_plot(cpi.differential.month.ts, title="Time Series Plot of Monthly CPI Pct Change Differential", Ytitle="CPI Pct Change Differential", Xtitle="Time");
# We can see there are four significant structural breaks (spikes/plunges) before 2010 and after 2020 (identify the periods)
# (2009.167, -71.42857) is (01/02/2009, -71.42857)
# (2009.417, 980769) is (01/05/2009, 98.0769)
# (2020.417, 56.25) is (01/05/2020, 56.25)
# (2020.667, -78.46154) is (01/08/2020, -78.45154)

# Identify the structural break periods
plot(efp(forex_month$cpi.differential.month~month, type="OLS-CUSUM"));
plot(Fstats(forex_month$cpi.differential.month~month));
cpi.breakpoints <- breakpoints(forex_month$cpi.differential.month~month, h=30);
cpi.breakpoints;
coef(cpi.breakpoints);
# It seems like R suggests that there are no structural breaks in cpi.differential.month

adf.test(cpi.differential.month.ts);
acf(cpi.differential.month.ts);
pacf(cpi.differential.month.ts);
# cpi.differential.month.ts is stationary

# Fit ARIMA model
auto.arima(cpi.differential.month.ts);
# cpi.differential.month.ts follows ARIMA(1,0,0)(1,0,0)[12]


# STL and ETS method of Model Decomposition
cpi.differential.month.z <- spec.pgram(cpi.differential.month.ts, fast=FALSE, taper=0.0);
stl1.3 <- stl(cpi.differential.month.ts, s.window="periodic");
summary(stl1.3);
autoplot(stl1.3);
forecast(stl1.3);

# Perform fourier analysis on the seasonal trend
stl1.3.seasonal <- stl1.3$time.series[1:185];
stl1.3.trend <- stl1.3$time.series[186:370];
stl1.3.remainder <- stl1.3$time.series[371:555];
fit1.3 <- lm(stl1.3.seasonal~sin(2*pi*month/12)+cos(2*pi*month/12)
             +sin(4*pi*month/12)+cos(4*pi*month/12)
             +sin(6*pi*month/12)+cos(6*pi*month/12));
summary(fit1.3);
stl1.3.fit <- stl1.3.trend+fit1.3$fitted.values;
autoplot(cbind(cpi.differential.month.ts,stl1.3.fit));

## In-Sample Forecasting and Diagnostics (Construct train-test split)
adf.test(cpi.differential.month.ts);
split <- ts_split(cpi.differential.month.ts, sample.out=12);
cpi.train <- split$train;
cpi.test <- split$test;

# Use an ARIMA Diagnostic Plot on the training set
arima_diag(cpi.train);
# Based on the diagnostic plot, we guess that our data follow SARIMA Model

# Build a SARIMA(2,0,1)(1,0,0)[12] model with training data set
sarima201.100 <- arima(cpi.train, order=c(2,0,1), seasonal=c(1,0,0));
autoplot(sarima201.100);
check_res(sarima201.100);

# Build a SARIMA(2,0,1)(1,0,1)[12] model with training data set
sarima201.101 <- arima(cpi.train, order=c(2,0,1), seasonal=c(1,0,1));
autoplot(sarima201.101);
check_res(sarima201.101);

# Build a SARIMA(1,0,0)(1,0,0)[12] model with training data set (suggested by auto.arima(cpi.train))
auto.arima(cpi.train);
sarima100.100 <- arima(cpi.train, order=c(1,0,0), seasonal=c(1,0,0));
autoplot(sarima100.100);
check_res(sarima100.100);

# Model selection based on forecast values and diagnostics
sarima201.100.fit <- forecast(sarima201.100 , h=12);
test_forecast(actual=cpi.differential.month.ts, forecast.obj=sarima201.100.fit, test=cpi.test);
accuracy(sarima201.100.fit, cpi.test);

sarima201.101.fit <- forecast(sarima201.101 , h=12);
test_forecast(actual=cpi.differential.month.ts, forecast.obj=sarima201.101.fit, test=cpi.test);
accuracy(sarima201.101.fit, cpi.test);

sarima100.100.fit <- forecast(sarima100.100 , h=12);
test_forecast(actual=cpi.differential.month.ts, forecast.obj=sarima100.100.fit, test=cpi.test);
accuracy(sarima100.100.fit, cpi.test);

# Define models with arguments (Construct train-validation-test split)
methods <- list(model1=list(method="arima",
                            method_arg=list(order=c(2,0,1), seasonal=list(order=c(1,0,0))),
                            notes="SARIMA(2,0,1)(1,0,0)"),
                model2=list(method="arima",
                            method_arg=list(order=c(2,0,1), seasonal=list(order=c(1,0,1))),
                            notes="SARIMA(2,0,1)(1,0,1)"),
                model3=list(method="arima",
                            method_arg=list(order=c(1,0,0), seasonal=list(order=c(1,0,0))),
                            notes="SARIMA(1,0,0)(1,0,0)")
);


# Train the models using backtesting
obj <- train_model(input=cpi.differential.month.ts, methods=methods, 
                   train_method=list(partitions=3, sample.out=12, space=3),
                   horizon=12, error="RMSE");

# Plot the model
plot_model(model.obj=obj);
# Seems like SARIMA(2,1,1)(1,0,1)[12] is the best model considering prediction power and model complexity (using mape and rmse)

# Generate the optimal fit based on our training data SARIMA(2,1,1)(1,0,1)[12]
fit1.3 <- auto.arima(cpi.differential.month.ts, seasonal=TRUE);
summary(fit1.3);
autoplot(fit1.3);
check_res(fit1.3);

# Generate final forecast based on our whole dataset (for period t+1)
forecast <- forecast(cpi.differential.month.ts, model=fit1.3, h=12);
plot_forecast(forecast);
summary(forecast);


## ----- 1.4 ARIMA of Trade Balance Percentage Change Differential (forex_month$trade.differential.month) -----
trade.differential.month.ts=ts(forex_month$trade.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(trade.differential.month.ts, main="Time Series Plot of Monthly Trade Balance Pct Change Differential");
ts_plot(trade.differential.month.ts, title="Time Series Plot of Monthly Trade Balance Pct Change Differential", Ytitle="Trade Pct Change", Xtitle="Time");
# We have two unusual spikes/plunges:
# (2009.167, -4717.896) is (01/02/2009, -4717.896)
# (2011.167, -4675.246) is (01/02/2011, 4675.246)
adf.test(trade.differential.month.ts);
# Since p-value is less than 0.01, we can reject H0, i.e. trade.differential.month.ts is stationary
auto.arima(trade.differential.month.ts);
# "Seems" like trade.differential.month.ts follows a ARIMA(0,0,0)(2,0,1)[12] model

# Identify the structural break periods
plot(efp(forex_month$trade.differential.month~month, type="OLS-CUSUM"));
plot(Fstats(forex_month$trade.differential.month~month));
trade.breakpoints <- breakpoints(forex_month$trade.differential.month~month, h=30);
trade.breakpoints;
coef(trade.breakpoints);
# There are 33 observations with structural breaks

# Replace structural break periods with NA
forex_month1 <- forex_month;
forex_month1$trade.differential.month[9] <- NA;
forex_month1$trade.differential.month[33] <- NA;
forex_month1 <- na.omit(forex_month1);
year3 <- year2;
month3 <- month2;
year3[9] <- NA;
year3[33] <- NA;
month3[9] <- NA;
month3[33] <- NA;
year3 <- na.omit(year3);
month3 <- na.omit(month3);

# Create new time series for trade without the structural break point (we determined that this structural break is an outlier that should be removed)
trade.differential.month.ts2 <- ts(forex_month1$trade.differential.month, start=c(year3[1],month3[1]), frequency=12);

# Perform the structural break test again
x <- seq(1,183,by=1);
plot(efp(forex_month1$trade.differential.month~x, type="OLS-CUSUM"));
plot(Fstats(forex_month1$trade.differential.month~x));
trade.breakpoints <- breakpoints(forex_month1$trade.differential.month~x, h=30);
trade.breakpoints;
coef(trade.breakpoints);
# This time we passed the structural break test

# Stationarity test and ARIMA decomposition
adf.test(trade.differential.month.ts2);
auto.arima(trade.differential.month.ts2);
# Once we removed the structural break, now trade does not follow SARIMA model, and looks more like white noise.
acf(trade.differential.month.ts2);
pacf(trade.differential.month.ts2);
checkresiduals(trade.differential.month.ts2);
autoplot(forecast(trade.differential.month.ts2));


## ----- 1.5 ARIMA of Unemployment Rate Differential (forex_month$unemployment.differential.month) ------
unemployment.differential.month.ts=ts(forex_month$unemployment.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(unemployment.differential.month.ts, main="Time Series Plot of Monthly Unemployment Rate Differential");
ts_plot(unemployment.differential.month.ts, title="Time Series Plot of Monthly Unemployment Rate Differential", Ytitle="Unemployment Rate Differential", Xtitle="Time");

autoplot(stl(unemployment.differential.month.ts, s.window="periodic"));
adf.test(unemployment.differential.month.ts)

a <- unemployment.differential.month.ts;


a <- a[59:185];
a <- a[55:176];
plot.ts(a)


adf.test(a)

















## ----- 2.1 Cointegration Tests (Engle Granger Test) -----

## Cointegration tests for weekly time series
# Cointegration test between price.ts and rate.differential.ts
autoplot(cbind(price.ts,rate.differential.ts));
adf.test(price.ts);
adf.test(rate.differential.ts);
mod1.5 <- lm(price.ts~rate.differential.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);
plot.ts(res1.5, start=c(year[1],week[1]), frequency=52)
abline(h=mean(res1.5),lty=2,col='red');
auto.arima(res1.5);

mod1.5 <- lm(rate.differential.ts~price.ts);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);
# 1. Both Forex price and TBill rate differential are unit root processes
# 2. Forex price and TBill rate differential are cointegrated

## Cointegration tests for monthly time series
# Cointegration test between price.month.ts and rate.differential.month.ts
autoplot(cbind(price.month.ts,rate.differential.month.ts));
adf.test(price.month.ts);
adf.test(rate.differential.month.ts);
mod1.5 <- lm(price.month.ts~rate.differential.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);
plot.ts(res1.5, start=c(year2[1],month2[1]), frequency=12)
abline(h=mean(res1.5),lty=2,col='red');
auto.arima(res1.5);
res1.5.diff <- diff(res1.5); 
adf.test(res1.5.diff);
plot.ts(res1.5.diff, start=c(year2[2],month2[2]), frequency=12)
abline(h=mean(res1.5.diff), lty=2, col='red');
acf(res1.5.diff);
pacf(res1.5.diff);
auto.arima(res1.5.diff);
autoplot(forecast(res1.5.diff));

mod1.5 <- lm(rate.differential.month.ts~price.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);
# 1. Both price.month.ts and rate.differential.month.ts are unit root processes
# 2. rate.differential.month.ts does not cointegrate with price.month.ts 
# 3. price.month.ts cointegrates with rate.differential.month.ts (first cointegrated relationship)

# Cointegration test between price.month.ts and cpi.differential.month.ts
autoplot(cbind(price.month.ts,cpi.differential.month.ts));
adf.test(price.month.ts);
adf.test(cpi.differential.month.ts);
# 1. price.month.ts is a unit root process but cpi.month.differential.ts is stationary
# 2. p-value is 0.3863 > 0.05, fail to reject H0: non-stationarity, i.e. no cointegration
# 3. In monthly time series, Forex price and CPI differential are not cointegrated

# Cointegration test between price.diff.month.ts and rate.differential.diff.month.ts
price.diff.month.ts <- diff(price.month.ts);
rate.differential.diff.month.ts <- diff(rate.differential.month.ts);
adf.test(price.diff.month.ts);
adf.test(rate.differential.diff.month.ts);
# 1. The first difference of both Forex price and TBill are stationary
# 2. The residuals of the first differences are stationary
# 3. The residuals follow an non-invertible MA(1) process

# Cointegration test between price.month.ts and trade.differential.month.ts
autoplot(cbind(price.month.ts, trade.differential.month.ts));
adf.test(price.month.ts);
adf.test(trade.differential.month.ts);
# 1. Trade balance percentage change differential is stationary

# Section 1 Summary
# 1. Forex Price (weekly/monthly) is a unit root process
# 2. TBill Rate Differential (weekly/monthly) is a unit root process
# 3. CPI percentage change differential (monthly) is stationary
# 4. Trade balance percentage change differential (monthly) is stationary
# 5. Forex Price and TBill Rate Differential (weekly) are cointegrated
# 6. The first difference of the residuals from regressing Forex Price against TBill Rate Differential (monthly) is stationary


## ----- 2.1 Granger Causality Tests -----

# Granger Causality Test between Forex Price and TBill Rate Differential (weekly)
grangertest(price.ts~rate.differential.ts,order=10);
grangertest(rate.differential.ts~price.ts,order=10);
# TBill rate differential granger causes Forex Price

# Granger Causality Test between Forex Price and TBill Rate Differential (monthly)
grangertest(price.month.ts~rate.differential.month.ts,order=10);
grangertest(rate.differential.month.ts~price.month.ts,order=10);
# rate.differential.month.ts granger causes price.month.ts
# price.month.ts granger causes rate.differential.month.ts

# Granger Causality Test between Forex Price and CPI Percentage Change Differential (monthly)
grangertest(price.month.ts~cpi.differential.month.ts,order=10);
grangertest(cpi.differential.month.ts~price.month.ts,order=10);
# Foex price granger causes CPI differential

# Granger Causality Test between Forex Price and Trade Balance Percentage Change Differential (monthly)
grangertest(price.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~price.month.ts,order=10);

# Granger Causality Test between TBill and CPI Percentage Change Differential (monthly)
grangertest(rate.differential.month.ts~cpi.differential.month.ts,order=10);
grangertest(cpi.differential.month.ts~rate.differential.month.ts,order=10);

# Granger Causality Test between TBill and Trade Balance Percentage Change Differential (monthly)
grangertest(rate.differential.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~rate.differential.month.ts,order=10);
# TBill rate differential granger causes trade balance percentage change differential (borderline)

# Granger Causality Test between CPI and Trade Balance Percentage Change Differential (monthly)
grangertest(cpi.differential.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~cpi.differential.month.ts,order=10);
# CPI granger causes trade differential


## ----- 2.2 VAR Model -----
# Overall plot of all variables
autoplot(cbind(price.diff.month.ts, rate.differential.month.ts, cpi.differential.month.ts, trade.differential.month.ts));

# Model Selection
VARselect(forex_month_VAR, lag.max=10, type="const")[["criteria"]];
VARselect(forex_month_VAR, lag.max=10, type="const")[["selection"]];
# Using the BIC, we fit a VAR(1) model (denoted by SC(n))

# VAR(1) Model
VAR1 <- VAR(forex_month_VAR, p=1, season=NULL, exog=NULL, type="const");
summary(VAR1);
AIC(VAR1);
BIC(VAR1);

# Serial Correlation Test
serial.test(VAR1,lags.pt=12,type="PT.asymptotic");
# There exists serial correlation in this VAR(1) model

# Heteroskedasticity Test
arch.test(VAR1, lags.multi=12, multivariate.only=TRUE);
# The model does not suffer from heteroskedasticity

# Normal Distribution Test of Residuals
normality.test(VAR1, multivariate.only=TRUE);
# From the Jarque-Bera Test, the residuals are not Normally distributed

# Structural Break/Stability Test of Residuals
stability <- stability(VAR1, type="OLS-CUSUM");
plot(stability);
# There is no structural breaks in the residuals

# VAR(3) Model
VAR3 <- VAR(forex_month_VAR,p=3, season=NULL, exog=NULL, type="const");
summary(VAR3);
AIC(VAR3);
BIC(VAR3);

# Serial Correlation Test
serial.test(VAR3,lags.pt=12,type="PT.asymptotic");
# There exists serial correlation in this VAR(1) model

# Heteroskedasticity Test
arch.test(VAR3, lags.multi=12, multivariate.only=TRUE);
# The model does not suffer from heteroskedasticity

# Normal Distribution Test of Residuals
normality.test(VAR3, multivariate.only=TRUE);
# From the Jarque-Bera Test, the residuals are not Normally distributed

# Structural Break/Stability Test of Residuals
stability <- stability(VAR3, type="OLS-CUSUM");
plot(stability);
# There is no structural breaks in the residuals


## ----- 2.3 Structural VAR Model -----

# Plot the time series objects
ts_plot(price.month.ts);
ts_plot(rate.differential.month.ts);
ts_plot(cpi.differential.month.ts);
ts_plot(trade.differential.month.ts);

# SVAR Restrictions
amat <- diag(4);
amat[1,3] <- NA;
amat[1,4] <- NA;
amat[2,3] <- NA;
amat[2,4] <- NA;
amat[3,2] <- NA;
amat[3,4] <- NA;
amat[4,1] <- NA;

# Construct the Model
VARselect(forex_month_VAR, lag.max=10, type="const")[["criteria"]];
VARselect(forex_month_VAR, lag.max=10, type="const")[["selection"]];
SVAR1 <- SVAR(VAR1, Amat=amat, Bmat=NULL, hessian=TRUE, estmethod=c("scoring","direct"));
SVAR1;

# Impulse Response Function
SVAR1.price.price <- irf(SVAR1, impulse="price.month", response="price.month");
SVAR1.price.price;
plot(SVAR1.price.price);

SVAR1.price.rate <- irf(SVAR1, impulse="rate.differential.month", response="price.month");
SVAR1.price.rate;
plot(SVAR1.price.rate);

SVAR1.rate.rate <- irf(SVAR1, impulse="rate.differential.month", response="rate.differential.month");
SVAR1.rate.rate;
plot(SVAR1.rate.rate);

SVAR1.rate.price <- irf(SVAR1, impulse="price.month", response="rate.differential.month");
SVAR1.rate.price;
plot(SVAR1.rate.price);

SVAR1.cpi.cpi <- irf(SVAR1, impulse="cpi.differential.month", response="cpi.differential.month");
SVAR1.cpi.cpi;
plot(SVAR1.cpi.cpi);

SVAR1.trade.trade <- irf(SVAR1, impulse="trade.differential.month", response="trade.differential.month");
SVAR1.trade.trade;
plot(SVAR1.trade.trade);

SVAR1.trade.rate <- irf(SVAR1, impulse="rate.differential.month", response="trade.differential.month");
SVAR1.trade.rate;
plot(SVAR1.trade.rate);

SVAR1.trade.cpi <- irf(SVAR1, impulse="cpi.differential.month", response="trade.differential.month");
SVAR1.trade.cpi;
plot(SVAR1.trade.cpi);


## ----- 2.4 VECM/Cointegrated VAR Model Model -----
# Overall plot of all variables
autoplot(cbind(price.diff.month.ts, rate.differential.month.ts, cpi.differential.month.ts, trade.differential.month.ts));

# Model Selection
VARselect(forex_month_VAR, lag.max=10, type="const")[["criteria"]];
VARselect(forex_month_VAR, lag.max=10, type="const")[["selection"]];

# Johansen Test
ca.jo.trace <- ca.jo(forex_month_VAR,ecdet="const",type=c("trace"),K=2);
summary(ca.jo.trace);
# Since we are able to reject rank <= 1 but fail to reject rank <= 2, 
#so we can conclude there are at most 2 cointegrated relationships 
#(which is also verified by the Engle-Granger Cointegration Tests)

ca.jo.eigen <- ca.jo(forex_month_VAR,ecdet="const",type=c("eigen"),K=2);
summary(ca.jo.eigen);
# The eigenvalue test yields the same results as the trace test

# Fit a VECM Model
## VECM with 1 lag and 2 cointegrated relationships
VECM1 <- VECM(forex_month_VAR,lag=2,r=2,estim='ML');
summary(VECM1);

# Diagnostic Tests
VECM1.VAR <- vec2var(ca.jo.trace, r=2);

# Serial Correlation Test
serial.test(VECM1.VAR, lags.pt=5, type="PT.asymptotic");
# We fail to reject H0, i.e. there is no serial correlation

# Heteroskedasticity Test
arch.test(VECM1.VAR, lags.multi=10, multivariate.only=TRUE);
# There is no varch effects in this model

# Normality Test on Residuals
normality.test(VECM1.VAR, multivariate.only=TRUE);
# Residuals are not Nomally distributed

# Impulse Functions
irf.price.rate <- irf(VECM1.VAR, impulse="rate.differential.month", response="price.month", n.ahead=12, boot=TRUE);
plot(irf.price.rate, ylab="Fx Price (Monthly)", main="Monthly TBill rate differential shock/impulse to Fx price");

irf.rate.price <- irf(VECM1.VAR, impulse="price.month", response="rate.differential.month", n.ahead=12, boot=TRUE);
plot(irf.rate.price, ylab="TBill Rate Differential (Monthly)", main="Monthly Fx price shock/impulse to TBill rate differential");

irf.cpi.price <- irf(VECM1.VAR, impulse="price.month", response="cpi.differential.month", n.ahead=12, boot=TRUE);
plot(irf.cpi.price, ylab="CPI Pct Change Differential (Monthly)", main="Monthly Fx price shock/impulse to CPI pct change differential");

irf.trade.rate <- irf(VECM1.VAR, impulse="rate.differential.month", response="trade.differential.month", n.ahead=12, boot=TRUE);
plot(irf.trade.rate, ylab="Trade Balance Pct Change Differential (Monthly)", main="TBill rate differential shock/impulse to Trade Balance Differential");

irf.trade.cpi <- irf(VECM1.VAR, impulse="cpi.differential.month", response="trade.differential.month", n.ahead=12, boot=TRUE);
plot(irf.trade.cpi, ylab="Trade Balance Pct Change Differential (Monthly)", main="CPI differential shock/impulse to Trade Balance Differential");

# Variance Decomposition
fevd1 <- fevd(VECM1.VAR, n.ahead=12);
plot(fevd1);

## VECM with 9 lags and 2 cointegrated relationships
# Johansen Test
ca.jo.trace <- ca.jo(forex_month_VAR,ecdet="const",type=c("trace"),K=9);
summary(ca.jo.trace);
# Since we are able to reject rank <= 1 but fail to reject rank <= 2, 
#so we can conclude there are at most 2 cointegrated relationships 
#(which is also verified by the Engle-Granger Cointegration Tests)

ca.jo.eigen <- ca.jo(forex_month_VAR,ecdet="const",type=c("eigen"),K=9);
summary(ca.jo.eigen);
# The eigenvalue test yields the same results as the trace test

VECM2 <- VECM(forex_month_VAR,lag=9,r=2,estim='ML');
summary(VECM2);

# Diagnostic Tests
VECM2.VAR <- vec2var(ca.jo.trace, r=2);

# Serial Correlation Test
serial.test(VECM2.VAR, lags.pt=5, type="PT.asymptotic");
# We fail to reject H0, i.e. there is no serial correlation

# Heteroskedasticity Test
arch.test(VECM2.VAR, lags.multi=10, multivariate.only=TRUE);
# There is no varch effects in this model

# Normality Test on Residuals
normality.test(VECM2.VAR, multivariate.only=TRUE);
# Residuals are not Nomally distributed

# Impulse Functions
irf.price.rate <- irf(VECM2.VAR, impulse="rate.differential.month", response="price.month", n.ahead=12, boot=TRUE);
plot(irf.price.rate, ylab="Fx Price (Monthly)", main="Monthly TBill rate differential shock/impulse to Fx price");

irf.rate.price <- irf(VECM2.VAR, impulse="price.month", response="rate.differential.month", n.ahead=12, boot=TRUE);
plot(irf.rate.price, ylab="TBill Rate Differential (Monthly)", main="Monthly Fx price shock/impulse to TBill rate differential");

irf.cpi.price <- irf(VECM2.VAR, impulse="price.month", response="cpi.differential.month", n.ahead=12, boot=TRUE);
plot(irf.cpi.price, ylab="CPI Pct Change Differential (Monthly)", main="Monthly Fx price shock/impulse to CPI pct change differential");

irf.trade.rate <- irf(VECM2.VAR, impulse="rate.differential.month", response="trade.differential.month", n.ahead=12, boot=TRUE);
plot(irf.trade.rate, ylab="Trade Balance Pct Change Differential (Monthly)", main="TBill rate differential shock/impulse to Trade Balance Differential");

irf.trade.cpi <- irf(VECM2.VAR, impulse="cpi.differential.month", response="trade.differential.month", n.ahead=12, boot=TRUE);
plot(irf.trade.cpi, ylab="Trade Balance Pct Change Differential (Monthly)", main="CPI differential shock/impulse to Trade Balance Differential");

# Variance Decomposition
fevd2 <- fevd(VECM2.VAR, n.ahead=12);
plot(fevd2);


AIC(VECM1.VAR)
AIC(VECM2.VAR)
BIC(VECM1.VAR)
BIC(VECM2.VAR)


## ----- 2.5 ARDL Model -----
# Check for stationarity
pp.test(price.month.ts); #non-stationary
pp.test(rate.differential.month.ts); #non-stationary
pp.test(cpi.differential.month.ts); #stationary
pp.test(trade.differential.month.ts); #stationary

# Estimate the ARDL Model
ARDL1 <- dynardl(price.month~rate.differential.month+cpi.differential.month+trade.differential.month,
                 data=forex_month_VAR,
                 lags=list("price.month"=1, "rate.differential.month"=1),
                 diffs="rate.differential.month",
                 ec=TRUE, simulate=FALSE);
summary(ARDL1);
dynardl.auto.correlated(ARDL1);


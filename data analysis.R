# ----- Setup -----
rm(list=ls());
library(readxl);
library(logr);
library(stargazer);
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
setwd("C:/Users/soulu/OneDrive - University of Toronto/University of Toronto/UTM 2021-2023/8_2024 Winter/ECO475H5S/Term Papaer/Data Analysis/Submission");
forex <- read_excel("FX.xlsx");
colnames(forex)[colnames(forex) == "week"] <- "date";
forex$pct_change_us_cpi <- replace_na(forex$pct_change_us_cpi, 0);
forex$pct_change_ca_cpi <- replace_na(forex$pct_change_ca_cpi, 0);
forex$pct_change_us_trade <- replace_na(forex$pct_change_us_trade, 0);
forex$pct_change_ca_trade <- replace_na(forex$pct_change_ca_trade, 0);
forex$us_employ_rate <- replace_na(forex$us_employ_rate, 0);
forex$ca_employ_rate <- replace_na(forex$ca_employ_rate, 0);

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
forex$unemployment_differential <- forex$us_employ_rate - forex$ca_employ_rate;

# Create forex_month dataframe
forex_month <- forex %>% group_by(year_month) %>% summarise(price.month=mean(price),
               rate.differential.month=mean(rate_differential),                                             
               cpi.differential.month=mean(pct_change_cpi_differential),
               trade.differential.month=mean(pct_change_trade_differential),
               unemployment.differential.month=mean(unemployment_differential)
);
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
adf.test(price.diff.month.ts); # price.month.ts is I(1)
plot(price.diff.month.ts, main="Time Series Plot of First-Differenced Monthly Fx Closing Price", ylab="Fx Price", xlab="Time");
abline(h=mean(price.diff.month.ts),lty=2,col='red');
acf(price.diff.month.ts);
pacf(price.diff.month.ts);
fit1.1 <- auto.arima(price.diff.month.ts);
summary(fit1.1);

## AR decomposition on price.ts using Yule Walker
price.month.ar <- ar.yw(price.month.ts);
price.month.ar;
plot.ts(price.month.ar$aic, main="Plot of Akaike's Criterion", xlab="Order", ylab="Akaike's Criterion (AIC)");
price.month.ar$aic;
plot.ts(price.month.ar$resid, main="Plot of Residuals from AR Model", xlab="Time", ylab="Residuals from AR Model")
abline(h=mean(price.month.ar$resid[2:132]),lty=2, col='red');
acf(price.month.ar$resid[2:132]);
price.month.ar$ar;
price.month.ar$var.pred;
price.month.ar$asy.var.coef;
autoplot(price.month.ar);
checkresiduals(price.month.ar);
autoplot(forecast(price.month.ar));
forecast(price.month.ar);
# Forex prices can also be fitted with an AR(2) model

## Verification using In-Sample Forecasting and Diagnostics (Construct training and testing set, i.e. train-test split)
adf.test(price.month.ts); # Non-stationary --> needs to take difference (verified that first difference is stationary)
adf.test(diff(price.month.ts)); # Stationary --> price.month is first-difference stationary, i.e. price.month.ts is I(1)
split <- ts_split(price.month.ts, sample.out=12);
price.train <- split$train;
price.test <- split$test;

# Use an ARIMA Diagnostic Plot on the training set
arima_diag(price.train);
# Based on the diagnostic plot, we guess that our data follow SARIMA Model

# Build a SARIMA(1,1,1)(0,0,1) model with training data set
sarima111.001 <- arima(price.train, order=c(1,1,1), seasonal=c(0,0,1));
autoplot(sarima111.001);
check_res(sarima111.001);

# Build a ARIMA(1,1,1) model with training data set
arima111 <- arima(price.train, order=c(1,1,1));
autoplot(arima111);
check_res(arima111);

# Build a ARIMA(2,1,2) model with training data set (suggested by auto.arima(price.train))
auto.arima(price.train); # Use r code to see which model the computer thinks is best fit
arima212 <- arima(price.train, order=c(2,1,2));
autoplot(arima212);
check_res(arima212);

# Build a ARIMA(2,1,2) model (with drift) with training data set 
arima212.drift <- arima(price.train, order=c(2,1,2), include.mean=TRUE);
autoplot(arima212.drift);
check_res(arima212.drift);

# Model selection based on forecast values and diagnostics
sarima111.001.fit <- forecast(sarima111.001, h=12); # Forecast 12 periods ahead (since we split the price.test to be the last 12 observations)
test_forecast(actual=price.month.ts, forecast.obj=sarima111.001.fit, test=price.test);
accuracy(sarima111.001.fit, price.test);

arima111.fit <- forecast(arima111, h=12);
test_forecast(actual=price.month.ts, forecast.obj=arima111.fit, test=price.test);
accuracy(arima111.fit, price.test);

arima212.fit <- forecast(arima212, h=12);
test_forecast(actual=price.month.ts, forecast.obj=arima212.fit, test=price.test);
accuracy(arima212.fit, price.test);

arima212.drift.fit <- forecast(arima212.drift, h=12);
test_forecast(actual=price.month.ts, forecast.obj=arima212.drift.fit, test=price.test);
accuracy(arima212.drift.fit, price.test);

# Define models with arguments (Construct train-validation-test split)
methods <- list(model1=list(method="arima",
                            method_arg=list(order=c(1,1,1), seasonal=list(order=c(0,0,1))),
                            notes="SARIMA(1,1,1)(0,0,1)"),
                model2=list(method="arima",
                            method_arg=list(order=c(1,1,1)),
                            notes="ARIMA(1,1,1)"),
                model3=list(method="arima",
                            method_arg=list(order=c(2,1,2)),
                            notes="ARIMA(2,1,2)"),
                model4=list(method="arima",
                            method_arg=list(order=c(2,1,2), include.mean=TRUE),
                            notes="ARIMA(2,1,2)(drift)")
                );

# Train the models using backtesting
obj <- train_model(input=price.month.ts, methods=methods, 
                   train_method=list(partitions=3, sample.out=12, space=3),
                   horizon=12, error="RMSE");

# Plot the model
plot_model(model.obj=obj);
# Seems like SARIMA(1,1,1)(0,0,1)[12] is the best model considering prediction power and model complexity

# Generate the optimal fit based on our training data ARIMA(1,1,1)
fit1.1 <- arima(price.month.ts, order=c(1,1,1), seasonal=c(0,0,1));
summary(fit1.1);
autoplot(fit1.1);
check_res(fit1.1);

# Generate final forecast based on our whole dataset (for period t+1)
forecast <- forecast(price.month.ts, model=fit1.1, h=12);
plot_forecast(forecast);
summary(forecast);


## ----- 1.2 ARIMA of Treasury Bill Rate Differential (forex$rate_differential) -----
## Fit an ARIMA model
rate.differential.month.ts=ts(forex_month$rate.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(rate.differential.month.ts, main="Time Series Plot of Monthly TBill Rate Differential", ylab="TBill Rate Differential", xlab="Time");
ts_plot(price.month.ts, title="Time Series Plot of Monthly TBill Rate Differential", Ytitle="TBill Rate Differential", Xtitle="Time");
adf.test(rate.differential.month.ts);
adf.test(diff(rate.differential.month.ts)); # rate.differential.month.ts is I(1)
acf(rate.differential.month.ts,lag.max = 50);
pacf(rate.differential.month.ts,lag.max = 50);
fit1.2 <- auto.arima(rate.differential.month.ts);
summary(fit1.2);
autoplot(fit1.2);
check_res(fit1.2);
checkresiduals(fit1.2);
autoplot(forecast(fit1.2));

## Take the first difference
rate.differential.diff.month.ts <- diff(rate.differential.month.ts);
adf.test(rate.differential.diff.month.ts); # rate.differential.month is I(1)
acf(rate.differential.diff.month.ts);
pacf(rate.differential.diff.month.ts);
fit1.2 <- auto.arima(rate.differential.diff.month.ts);
summary(fit1.2);
checkresiduals(fit1.2);
autoplot(fit1.2);
autoplot(forecast(fit1.2));
forecast(fit1.2);


## ----- 1.3 ARIMA of CPI Percentage Change Differential (forex_month$cpi.differential.month) -----
cpi.differential.month.ts=ts(forex_month$cpi.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(cpi.differential.month.ts, main="Time Series Plot of Monthly CPI Pct Change Differential", ylab="CPI Pct Change Differential", xlab="Time");
ts_plot(cpi.differential.month.ts, title="Time Series Plot of Monthly CPI Pct Change Differential", Ytitle="CPI Pct Change Differential", Xtitle="Time");
# We suspect structural breaks here around 2020 May and 2020 Aug

adf.test(cpi.differential.month.ts);
acf(cpi.differential.month.ts);
pacf(cpi.differential.month.ts);
# cpi.differential.month.ts is stationary

# Fit ARIMA model
fit1.3 <- auto.arima(cpi.differential.month.ts);
summary(fit1.3);
autoplot(fit1.3);
check_res(fit1.3);
checkresiduals(fit1.3);
autoplot(forecast(fit1.3));
forecast(fit1.3);


## ----- 1.4 ARIMA of Trade Balance Percentage Change Differential (forex_month$trade.differential.month) -----
trade.differential.month.ts=ts(forex_month$trade.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(trade.differential.month.ts, main="Time Series Plot of Monthly Trade Balance Pct Change Differential");
ts_plot(trade.differential.month.ts, title="Time Series Plot of Monthly Trade Balance Pct Change Differential", Ytitle="Trade Pct Change", Xtitle="Time");
# Here we suspect structural breaks again

adf.test(trade.differential.month.ts);
# Since p-value is less than 0.01, we can reject H0, i.e. trade.differential.month.ts is stationary
auto.arima(trade.differential.month.ts);

# Stationarity test and ARIMA decomposition
adf.test(trade.differential.month.ts);
fit1.4 <- auto.arima(trade.differential.month.ts);
acf(trade.differential.month.ts);
pacf(trade.differential.month.ts);
checkresiduals(fit1.4);
check_res(fit1.4);
autoplot(fit1.4);
autoplot(forecast(fit1.4));
forecast(fit1.4);


## ----- 1.5 ARIMA of Unemployment Rate Differential (forex_month$unemployment.differential.month) ------
unemployment.differential.month.ts=ts(forex_month$unemployment.differential.month, start=c(year2[1],month2[1]), frequency=12);
plot(unemployment.differential.month.ts, main="Time Series Plot of Monthly Unemployment Rate Differential");
ts_plot(unemployment.differential.month.ts, title="Time Series Plot of Monthly Unemployment Rate Differential", Ytitle="Unemployment Rate Differential", Xtitle="Time");

# Unit Root Test
adf.test(unemployment.differential.month.ts);
adf.test(diff(unemployment.differential.month.ts)); # unemployment.differential.month.ts is I(1)
acf(unemployment.differential.month.ts);
pacf(unemployment.differential.month.ts);

# Fit ARIMA Model
fit1.5 <- auto.arima(unemployment.differential.month.ts);
summary(fit1.5);
autoplot(fit1.5);
checkresiduals(fit1.5);
check_res(fit1.5);
autoplot(forecast(fit1.5));
forecast(fit1.5);


## ----- 2.1 Cointegration Tests (Engle Granger Test) -----

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

# Cointegration test between price.month.ts and unemployment.differential.month.ts
autoplot(cbind(price.month.ts, unemployment.differential.month.ts));
adf.test(price.month.ts);
adf.test(unemployment.differential.month.ts);
mod1.5 <- lm(price.month.ts~unemployment.differential.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);

mod1.5 <- lm(unemployment.differential.month.ts~price.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);

# Cointegration test between trade.differential.month.ts and unemployment.differential.month.ts
autoplot(cbind(rate.differential.month.ts, unemployment.differential.month.ts));
adf.test(rate.differential.month.ts);
adf.test(unemployment.differential.month.ts);
mod1.5 <- lm(rate.differential.month.ts~unemployment.differential.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);

mod1.5 <- lm(unemployment.differential.month.ts~rate.differential.month.ts);
summary(mod1.5);
res1.5 <- mod1.5$residuals;
adf.test(res1.5);


## ----- 2.1 Granger Causality Tests -----
# Granger Causality Test between Forex Price and TBill Rate Differential 
grangertest(price.month.ts~rate.differential.month.ts,order=10);
grangertest(rate.differential.month.ts~price.month.ts,order=10);
grangertest(price.diff.month.ts~rate.differential.diff.month.ts,order=10);
grangertest(rate.differential.diff.month.ts~price.diff.month.ts,order=10);

# Granger Causality Test between Forex Price and CPI Percentage Change Differential 
grangertest(price.month.ts~cpi.differential.month.ts,order=10);
grangertest(cpi.differential.month.ts~price.month.ts,order=10);
grangertest(price.diff.month.ts~cpi.differential.month.ts[2:132],order=10);
grangertest(cpi.differential.month.ts[2:132]~price.diff.month.ts,order=10);
# Forex price granger causes CPI differential
# First-difference forex price granger causes CPI differential at 10% significance level

# Granger Causality Test between Forex Price and Trade Balance Percentage Change Differential 
grangertest(price.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~price.month.ts,order=10);
grangertest(price.diff.month.ts~trade.differential.month.ts[2:132],order=10);
grangertest(trade.differential.month.ts[2:132]~price.diff.month.ts,order=10);

# Granger Causality Test between TBill and CPI Percentage Change Differential 
grangertest(rate.differential.month.ts~cpi.differential.month.ts,order=10);
grangertest(cpi.differential.month.ts~rate.differential.month.ts,order=10);
grangertest(rate.differential.diff.month.ts~cpi.differential.month.ts[2:132],order=10);
grangertest(cpi.differential.month.ts[2:132]~rate.differential.diff.month.ts,order=10);

# Granger Causality Test between TBill and Trade Balance Percentage Change Differential 
grangertest(rate.differential.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~rate.differential.month.ts,order=10);
grangertest(rate.differential.diff.month.ts~trade.differential.month.ts[2:132],order=10);
grangertest(trade.differential.month.ts[2:132]~rate.differential.diff.month.ts,order=10);
# Trade balance pct change differential granger causes TBill rate
# Trade balance pct change differential granger causes first difference of TBill rate

# Granger Causality Test between CPI and Trade Balance Percentage Change Differential 
grangertest(cpi.differential.month.ts~trade.differential.month.ts,order=10);
grangertest(trade.differential.month.ts~cpi.differential.month.ts,order=10);
# CPI granger causes trade differential

# Granger Causality Test between Forex Price and Unemployment Rate Differential
grangertest(price.month.ts~unemployment.differential.month.ts,order=10);
grangertest(unemployment.differential.month.ts~price.month.ts,order=10);
unemployment.differential.diff.month.ts <- diff(unemployment.differential.month.ts);
grangertest(price.diff.month.ts~unemployment.differential.diff.month.ts,order=10);
grangertest(unemployment.differential.diff.month.ts~price.diff.month.ts,order=10);

# Granger Causality Test between TBill and Unemployment Rate Differential
grangertest(rate.differential.month.ts~unemployment.differential.month.ts,order=10);
grangertest(unemployment.differential.month.ts~rate.differential.month.ts,order=10);
grangertest(rate.differential.diff.month.ts~unemployment.differential.diff.month.ts,order=10);
grangertest(unemployment.differential.diff.month.ts~rate.differential.diff.month.ts,order=10);
# TBill granger causes umemployment rate differential at 10% significance level

# Granger Causality Test between CPI and Unemployment Rate Differential
grangertest(cpi.differential.month.ts~unemployment.differential.month.ts,order=10);
grangertest(unemployment.differential.month.ts~cpi.differential.month.ts,order=10);
grangertest(unemployment.differential.diff.month.ts~cpi.differential.month.ts[2:132],order=10);
grangertest(cpi.differential.month.ts[2:132]~unemployment.differential.diff.month.ts,order=10);
# CPI granger causes unemployment rate differential
# CPI granger causes first difference of unemployment rate differential

# Granger Causality Test between Trade Balance and Unemployment Rate Differential
grangertest(trade.differential.month.ts~unemployment.differential.month.ts,order=10);
grangertest(unemployment.differential.month.ts~trade.differential.month.ts,order=10);
grangertest(unemployment.differential.diff.month.ts~trade.differential.month.ts[2:132],order=10);
grangertest(trade.differential.month.ts[2:132]~unemployment.differential.diff.month.ts,order=10);


## ----- 2.2 VAR Model -----
# Overall plot of all variables
autoplot(cbind(price.diff.month.ts, rate.differential.diff.month.ts, cpi.differential.month.ts, 
               trade.differential.month.ts, unemployment.differential.diff.month.ts), 
               main="Time Series Plot of All 5 Variables Overlaid", 
               ylab="All Variables", xlab="Time");

# Construct the matrix for the VAR model
Y <- cbind(price.diff.month.ts, rate.differential.diff.month.ts, cpi.differential.month.ts[2:132], 
           trade.differential.month.ts[2:132], unemployment.differential.diff.month.ts);

# Model Selection
VARselect(Y, lag.max=6, type="const")[["criteria"]];
VARselect(Y, lag.max=6, type="const")[["selection"]];
# Using the BIC, we fit a VAR(1) model (denoted by SC(n))

# VAR(1) Model
VAR1 <- VAR(Y, p=1, season=NULL, exog=NULL, type="const");
summary(VAR1);
stargazer(VAR1[["varresult"]], type="latex");
AIC(VAR1);
BIC(VAR1);

# Serial Correlation Test
serial.test(VAR1,lags.pt=1,type="PT.asymptotic");
# There exists serial correlation in this VAR(1) model

# Heteroskedasticity Test
arch.test(VAR1, lags.multi=1, multivariate.only=TRUE);
# The model does not suffer from heteroskedasticity

# Normal Distribution Test of Residuals
normality.test(VAR1, multivariate.only=TRUE);
# From the Jarque-Bera Test, the residuals are not Normally distributed

# Structural Break/Stability Test of Residuals
stability <- stability(VAR1, type="OLS-CUSUM");
plot(stability);
# There is no structural breaks in the residuals

# VAR(5) Model
VAR5 <- VAR(Y,p=5, season=NULL, exog=NULL, type="const");
summary(VAR5);
stargazer(VAR5[["varresult"]], type="latex");
AIC(VAR5);
BIC(VAR5);

# Serial Correlation Test
serial.test(VAR5,lags.pt=5,type="PT.asymptotic");
serial.test(VAR5,lags.pt=17,type="PT.asymptotic");
# There exists serial correlation in this VAR(1) model

# Heteroskedasticity Test
arch.test(VAR5, lags.multi=5, multivariate.only=TRUE);
# The model does not suffer from heteroskedasticity

# Normal Distribution Test of Residuals
normality.test(VAR5, multivariate.only=TRUE);
# From the Jarque-Bera Test, the residuals are not Normally distributed

# Structural Break/Stability Test of Residuals
stability <- stability(VAR5, type="OLS-CUSUM");
plot(stability);
# There is no structural breaks in the residuals



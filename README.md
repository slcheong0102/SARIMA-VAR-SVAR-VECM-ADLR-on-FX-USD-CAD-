This project analyses macroeconomic forces and dynamics of FX Price with Inflation Rate Differential (captured by CPI Percentage Change Differential), 
Interest Rate Differential (captured by Treasury Bill Rate Percentage Change Differential), Trade Balance Percentage Change Differenital, and Unemployment Rate Differential. 
We performed time series model fitting on each variable using the Augmented Dickey Fuller's Test, train-validate-test split method (machine learning), spectral/fourier analysis,
etc. and performed model selection based on AIC/BIC information criterion minimization with exact likelihood estimates of the parameters.
We also performed the Engle Granger Cointegration Test on the unit root processes (random walks), which are FX Price, Interest Rate Differential, and Unemployment Rate Differential; 
as well as the Johansen Test to see the order of cointegration in constructing a Vector Error Correction Model (VECM). In addition, we performed the Granger Causality Test, 
and Impulse Response Funciton on pairwise variables to investigate causality between variables such that we can potentially construct a Structured Vector Autoregressive Model (SVAR) or an 
Autoregressive Distributed Lag Model (ARDL). Note that each model is put under a series of robustness tests, such as serial correlation test, autoregressive conditional heteroskedasticity test, 
Normality test, stability test, etc. to verify the stability (stationarity/invertibility) and goodness of fit of the model. A 12-period (months) forecast is also performed for each variable and 
plots are produced to visualize the seasonality and trend of each stochastic processes.

An additonal place we can look into in future projects is to do per-second time interval analysis for FX Price or other financial instrutments/derivatives in the spot, futures, or forwards market. In 
this case, we can apply Brownian Motion, stochastic differential equations, Ito's Integral to analyze short-term arbitraging opportunity, and hence explore algorithmic trading possibilities.

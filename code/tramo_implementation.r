library(tseries)
library(seasonal)
library(forecast)
source('code/helper_functions.r')
# we start with loading the data 
data = c(105.4, 106.4, 106.5, 107.5, 109.1, 112.4, 115.2, 117.3, 119.0, 121.1, 123.9, 118.7,
    115.6, 114.8, 115.7, 117.4, 118.8, 120.5, 125.4, 127.5, 126.4, 125.8, 125.3, 123.4,
    122.7, 122.7, 122.8, 123.4, 124.5, 127.1, 128.1, 130.3, 131.3, 132.4, 132.9, 131.3, 
    131.1, 129.2, 129.2, 131.3, 133.8, 137.0, 138.8, 138.0, 136.5, 136.8, 135.6, 133.1,
    131.9, 131.8, 131.8, 132.1, 132.4, 134.1, 138.3, 140.1, 138.2, 139.4, 141.5, 139.7,
    138.1, 136.1, 135.5, 135.8, 136.5, 138.0, 140.1, 140.5, 138.9, 138.2, 137.8, 136.0,
    135.0, 135.1, 135.9, 137.3, 139.0, 141.1, 143.4, 144.7, 146.0, 149.1, 151.6, 155.3, 
    153.4, 149.7, 147.8, 153.4, 151.8, 153.4, 156.7, 157.8, 161.6, 165.5, 166.0, 160.6, 
    156.4, 155.5, 155.0, 156.4, 159.4, 161.3, 162.9, 162.7, 162.7, 166.9, 169.1, 167.1, 
    164.9, 164.6, 166.9, 169.4, 172.1, 173.8, 173.8, 175.1, 176.7, 178.6, 177.0, 174.1, 
    174.8, 174.4, 174.9, 175.9, 177.2, 181.7, 193.8, 192.5, 188.4, 190.4, 192.4, 190.7, 
    189.3, 189.5, 189.8, 191.2, 192.6, 198.7, 204.3, 203.4)
diwali_ind = c(0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0405, -0.0405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.6595,  0.6595,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.3595,  0.3595  ,0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.6595,  0.6595  ,0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.0595,  0.0595  ,0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -0.6595,  0.6595  ,0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.3405, -0.3405,  0.0000)

# create  timeseries object

y = ts(data, start = c(2013, 1), frequency = 12)
diwali = ts(diwali_ind[1:length(y)], start = c(2013, 1), frequency =  12)
print(y)
# we fit the default model on y and log(y) an compare the aic test to select. transformation
# This corresponds to the following in the spec file:
#           transform{function  = auto}
no_transformation_test_fit = Arima(y , order =c(0,1,1), seasonal = c(0,1,1),include.mean = F, method = 'ML', lambda = NULL)
log_transformation_test_fit = Arima(log(y), order =c(0,1,1), seasonal = c(0,1,1),include.mean = F,  method = 'ML',  lambda = NULL)

plot_model_output(no_transformation_test_fit, series_name = "No Transformation")
plot_model_output(log_transformation_test_fit,series_name = "Log Transformed Series", adj_inf_criteria_for_log = T)

# note that the values of AICC values are matching those we get from the X13-ArimaSeats 
# if AICC_nolog - AICC_log < -2 (default of delta AICC) then prefer no log transform . See page 218 of manual
# 529.65 - 514.66 = 14.99 : not less than -2 so Prefer Log transform

Z = log(y)


## Default Model Estimation : initial outlier identification, tests for trading day and easter

# including the interfears with the log level tests so i am not concidering them right now

# test for the presence of constant term by a t-test on the residialsof default model with no intercept


const = const_term(y, 1, 1)
const = rep(1, 140)

default_model_fit1 = Arima(Z ,order =c(0,1,1), seasonal = c(0,1,1),include.mean = F,  method = 'ML',  lambda = NULL)

plot_model_output( default_model_fit1, adj_inf_criteria_for_log = T)

t.test(residuals(default_model_fit1))
# the p value of 0.5394 indicates that the constant term is not significant
# the Null Hypothesis of absence of constant term is rejected if |t| < 1.96: See page 72


# Now we fit the default model and try to find outliers
# see page 40 for robust estimate formula for sig
# Deafault critical value is 3.88 for the test See table 7.22


Arima(default_model_fit1$x,xreg = LS(Z, 45),  order = as.numeric(paste0(default_model_fit1$call$order)[2:4]),
      seasonal = as.numeric(paste0(default_model_fit1$call$seasonal)[2:4]), 
      include.mean = FALSE, method = 'ML', 
      fixed = c(coef(default_model_fit1), NA))

for(i in 3:(length(Z)-1)){
    ress_t = detect_outlier(default_model_fit1, i, "LS")
    if(ress_t>3.88){
        print(ress_t)
    }
}

for(i in 2:(length(Z)-1)){
    ress_t = detect_outlier(default_model_fit1, i, "AO")
    if(ress_t>3.88){
        print(ress_t)
    }
}
### values not matching 
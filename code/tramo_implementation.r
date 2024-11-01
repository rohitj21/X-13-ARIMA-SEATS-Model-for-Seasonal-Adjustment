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

# for diwali regressor we have taken the impact period to be 10 days before diwali (inclusive)
td1nolpyear = td1nolpyear_regressor(start = c(2013,1), n = length(data))
# It matches exactly with the regressor being used in the program 
# here is the image of the regressor of the program : <insert image>


# create  timeseries object

y = ts(data, start = c(2013, 1), frequency = 12)
diwali = ts(diwali_ind[1:length(y)], start = c(2013, 1), frequency =  12)
td1nolpyear = ts(td1nolpyear, start = c(2013, 1), frequency = 12)
print(td1nolpyear)
print(y)
# we fit the default model on y and log(y) an compare the aic test to select. transformation
# This corresponds to the following in the spec file:
#           transform{function  = auto}

# regression matrix 
Xreg = cbind(td1nolpyear, diwali)
no_transformation_test_fit = Arima(y , xreg = Xreg,  order =c(0,1,1), seasonal = c(0,1,1),include.mean = F, method = 'ML', lambda = NULL)
log_transformation_test_fit = Arima(log(y), xreg= Xreg, order =c(0,1,1), seasonal = c(0,1,1),include.mean = F,  method = 'ML',  lambda = NULL)

plot_model_output(no_transformation_test_fit, series_name = "No Transformation")
plot_model_output(log_transformation_test_fit,series_name = "Log Transformed Series", adj_inf_criteria_for_log = T)

# note that the values of AICC values are matching those we get from the X13-ArimaSeats 
# if AICC_nolog - AICC_log < -2 (default of delta AICC) then prefer no log transform . See page 218 of manual
# 529.65 - 514.66 = 14.99 : not less than -2 so Prefer Log transform

Z = log(y)


## Default Model Estimation : initial outlier identification, tests for trading day


# test for the presence of constant term by a t-test on the residialsof default model with no intercept


const = const_term(y, 1, 1) # funtion to create constant term


default_model_fit1 = Arima(Z ,xreg = Xreg ,order =c(0,1,1), seasonal = c(0,1,1),include.mean = F,  method = 'ML',  lambda = NULL)

plot_model_output( default_model_fit1, adj_inf_criteria_for_log = T)

t.test(residuals(default_model_fit1))
# the p value of 0.5391 indicates that the constant term is not significant
# the Null Hypothesis of absence of constant term is rejected if |t| < 1.96: See page 72

# We will not use the constant term in the model

# AIC test for td1nolpyear and diwali
plot_model_output(Arima(Z,xreg = td1nolpyear,  order = c(0,1,1),seasonal = c(0,1,1), include.mean = F,  method = 'ML',  lambda = NULL) , T)
plot_model_output(Arima(Z,  order = c(0,1,1),seasonal = c(0,1,1), include.mean = F,  method = 'ML',  lambda = NULL) , T)
# AICC_with - AICC_without = -754.38 + 754.18  = 0.2 
# 0.2 + aiccdiff(defaul = 0) = 0.2 > 0 
# trading day regressor will not be included in the model

plot_model_output(Arima(Z,xreg = diwali,  order = c(0,1,1),seasonal = c(0,1,1), include.mean = F,  method = 'ML',  lambda = NULL) , T)
plot_model_output(Arima(Z,  order = c(0,1,1),seasonal = c(0,1,1), include.mean = F,  method = 'ML',  lambda = NULL) , T)
# AICC_with - AICC_without = -754.57 + 754.18  = 0.39
# 0.39 + aiccdiff(defaul = 0) = 0.39 > 0 
# diwali regressor will not be included in the model


# Now we fit the default model and try to find outliers
# see page 40 of the manual for robust estimate formula for sig
# Deafault critical value is 3.88 for the test See table 7.22

Xreg = NULL

default_model = Arima(Z,  order = c(0,1,1),seasonal = c(0,1,1), xreg = Xreg,
                     include.mean = F,  method = 'ML',  lambda = NULL) 

default_model

# note that the coeff are negative of what we get from X13 arima seats 
#this is because of different parameterizations of the arima model used by Arima()
# ref: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/arima

print("forward pass 1")
curr_outlier = forward_pass(default_model,xreg = Xreg, types = c("AO", "LS"),
                            tcritical = 3.88)
k = 2
while(! is.null(curr_outlier)){
  Xreg = cbind(Xreg, curr_outlier)
  print(paste("forward pass", k))
  k = k +1
  curr_outlier = forward_pass(default_model,xreg = Xreg, types = c("AO", "LS"),
                              tcritical = 3.88)
}


k = 1

while(T){
  print(paste("backward_pass", k))
  k = k+1
  ind = backward_pass(default_model, Xreg, tcritical = 3.88)
  if(is.null(ind)){
    print("outlier detection is done")
    break
  }
  else{
    Xreg = Xreg[,-ind]
  }

}
# outliers identified from the default model
print(colnames(Xreg))
# these match with those from the program

# Our next goul is to understand and implement the IGLS algorithm so that 
#we can match the results of paramenter estimates with the program

library(ggplot2)
library(gridExtra)
library(forecast)
td1nolpyear_regressor <- function(start = c(2023, 1), n) {
    # Initialize an empty vector to store the regressor values
    regressor <- numeric(n)
    
    # Loop through each month from the start date for n months
    for (i in 0:(n - 1)) {
        # Calculate the year and month for the current iteration
        year <- start[1] + (start[2] + i - 1) %/% 12
        month <- (start[2] + i - 1) %% 12 + 1
        
        # Start and end dates for the current month
        start_date <- as.Date(paste(year, sprintf("%02d", month), "01", sep = "-"))
        end_date <- seq(start_date, by = "month", length.out = 2)[2] - 1
        
        # Get all dates in the current month
        dates <- seq(start_date, end_date, by = "day")
        
        # Count weekdays and weekends (Saturday and Sunday)
        weekdays_count <- sum(!weekdays(dates) %in% c("Saturday", "Sunday"))
        weekends_count <- sum(weekdays(dates) %in% c("Saturday", "Sunday"))
        
        # Calculate the regressor for the month
        regressor[i + 1] <- weekdays_count - (5 / 2) * weekends_count
    }
    
    return(regressor)
}


selection_criteria <- function(model, adjust_for_log_transform = F) {
    
    # Extract log-likelihood from the model
    loglik <- as.numeric(logLik(model))  # Extract log-likelihood
    
    # Extract original series and frequency from model
    original_series <- model$x
    freq <- frequency(model$x)
    
    # Extract non-seasonal and seasonal orders from model$call
    arima_order <- eval(model$call$order)      # Non-seasonal order (p, d, q)
    seasonal_order <- eval(model$call$seasonal) # Seasonal order (P, D, Q)
    
    number_of_regressors = dim(eval(model$call$xreg))[2]
    if(is.null(number_of_regressors)){
        number_of_regressors = 1
    }
    
    # Non-seasonal and seasonal differencing
    order_diff <- arima_order[2]         # Non-seasonal differencing (d)
    seasonal_diff <- seasonal_order[2]   # Seasonal differencing (D)
    
    # Calculate the effective number of observations
    effective_n <- length(original_series) - order_diff - freq * seasonal_diff
    
    if(adjust_for_log_transform){
        loglik = loglik - sum(original_series[(length(original_series) - effective_n + 1):(length(original_series))])
    }
    # Calculate the number of parameters (including the variance term)
    # The number of parameters includes AR (p), MA (q), Seasonal AR (P), Seasonal MA (Q), and variance.
    n_p <- number_of_regressors +  arima_order[1] + arima_order[3] +
        seasonal_order[1] + seasonal_order[3] + 1  # Add 1 for variance
    
    # Calculate AIC
    aic <- -2 * loglik + 2 * n_p
    
    # Calculate BIC
    bic <- -2 * loglik + n_p * log(effective_n)
    
    # Calculate AICc (corrected AIC)
    # Calculate AICc (corrected AIC)
    aicc <- aic + (2 * n_p * (n_p + 1)) / (effective_n - n_p - 1)
    
    # Calculate Hannan-Quinn Criterion
    hq <- -2 * loglik + 2 * n_p * log(log(effective_n))
    
    # Return the results as a list
    return(list(AIC = aic, BIC = bic, AICc = aicc, HQ = hq, npar = n_p))
}

    
plot_model_output <- function(model, series_name = "Series", adj_inf_criteria_for_log = F) {
    
    # Extract the original series, fitted values, and residuals from the model
    original_series <- model$x
    fitted_values <- fitted(model)
    residuals <- residuals(model)
    
    # Perform Ljung-Box test for residuals independence
    ljung_box_test <- Box.test(residuals, lag = 20, type = "Ljung-Box")
    lb_pvalue <- round(ljung_box_test$p.value, 4)
    
    # Extract model coefficients and standard errors
    coef_vals <- round(coef(model), 4)
    coef_errors <- round(sqrt(diag(vcov(model))), 4)
    
    # Create the model equation as a string
    model_eq <- paste0("ARIMA",substr(paste(model$call)[3], 2, 10), 
                       " x ", substr(paste(model$call)[4], 2, 10), " Method = ", paste(model$call$method))
    
    # Use the selection_criteria function to get AIC, BIC, AICc, and HQ
    criteria <- selection_criteria(model,adjust_for_log_transform = adj_inf_criteria_for_log)
    
    
    n_p <- criteria$npar
    
    # Create a dataframe to store the original series, fitted values, and residuals
    df <- data.frame(
        Time = as.numeric(time(original_series)),
        Original = as.numeric(original_series),
        Fitted = as.numeric(fitted_values),
        Residuals = as.numeric(residuals)
    )
    
    # Plot the original series and fitted values
    p1 <- ggplot(df, aes(x = Time)) +
        geom_line(aes(y = Original), color = "blue", linewidth = 1, show.legend = TRUE) +
        geom_line(aes(y = Fitted), color = "red", linetype = "dashed", linewidth = 1, show.legend = TRUE) +
        ggtitle(paste(series_name, ": Original vs Fitted")) +
        xlab("Time") + ylab("Values") +
        theme_minimal() +
        labs(subtitle = paste("Model Summary: AIC =", round(criteria$AIC, 2), 
                              ", npar= ", n_p,
                              ", BIC =", round(criteria$BIC, 2), 
                              ", AICc =", round(criteria$AICc, 2),
                              ", HQ =", round(criteria$HQ, 2),
                              ", Log-Likelihood =", round(logLik(model), 2),
                              "\nLjung-Box Test p-value:", lb_pvalue,
                              "\nModel Equation:", model_eq,
                              "\nCoefficients (Estimate ± Std. Error):\n", 
                              paste(names(coef_vals), "=", coef_vals, "±", coef_errors, collapse = ", ")))
    
    # Plot the residuals
    p2 <- ggplot(df, aes(x = Time, y = Residuals)) +
        geom_line(color = "purple", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        ggtitle(paste(series_name, ": Residuals")) +
        xlab("Time") + ylab("Residuals") +
        theme_minimal()
    
    # Print the model summary and Ljung-Box test result on the console
    print(summary(model))
    print(ljung_box_test)
    
    # Combine the plots
    grid.arrange(p1, p2, nrow = 2)
}

const_term <-  function(y, d, D){
    return(diffinv(diffinv(rep(1, length(y)), lag = 1, differences =d), lag = 12, differences = D)[-(1:(d+D*frequency(y)))])
}
AO <- function(y, t) {
    time_point = start(y)
    time_point[2] = time_point[2] + (t-1)
    time_point[1]  = time_point[1] + ((time_point[2]-1) %/% 12)  
    time_point[2] = ifelse( time_point[2] %% 12 == 0 , 12, time_point[2]%%12)
    new_var_name = paste("AO",time_point[1] , ".", time_point[2], sep  = "")
    if(time_point[2] <10){
        new_var_name = paste("AO",time_point[1] , ".0", time_point[2], sep  = "")
    }

    x = as.integer(seq_along(y) == t)
    x = matrix(x, ncol = 1)
    colnames(x)[1] <- new_var_name
    return(x)    
    
}
LS <- function(y, t) {
    time_point = start(y)
    time_point[2] = time_point[2] + (t-1)
    time_point[1]  = time_point[1] + ((time_point[2]-1) %/% 12)  
    time_point[2] = ifelse( time_point[2] %% 12 == 0 , 12, time_point[2]%%12)
    new_var_name = paste("LS", time_point[1] , ".", time_point[2], sep  = "")
    if(time_point[2] <10){
        new_var_name = paste("LS",time_point[1] , ".0", time_point[2], sep  = "")
    }
    x = as.integer(seq_along(y) >= t)-1

    x = matrix(x, ncol = 1)
    colnames(x)[1] <- new_var_name
    return(x)
    
}


outlier_t_value <- function(model, t, xreg = NULL, type = c("AO, LS")){
    y = model$x
    X = xreg
    nreg = ifelse(is.null(xreg) ,  0, dim(xreg)[2])
    ncoeff = sum(c(eval(model$call$order)[-2] , eval(model$call$seasonal)[-2]))
    model_order = eval(model$call$order)
    model_seasonal = eval(model$call$seasonal)
    new_variable = NULL
    if(type == 'AO'){
        new_variable = AO(y, t)
    }
    else if(type == "LS"){
        new_variable = LS(y, t)
    }
    X = cbind(X, new_variable)
    if(abs(det(t(X) %*% X)) < 1e-10){
        return(list(tvalue = 0, time_point =colnames(X)[dim(X)[2]] ))
    }
    fit = Arima(y, order = model_order,
                seasonal = model_seasonal,
                xreg = X,
                include.mean = FALSE, method = 'ML', 
                fixed = c(coef(model)[1:ncoeff], rep(NA, nreg + 1)))
    ind = length(fit$coef)
    ind2 = dim(fit$var.coef)[1]
    return(list(tvalue = abs((fit$coef)[ind]) / sqrt((fit$var.coef)[ind2, ind2]), time_point = colnames(X)[dim(X)[2]]) )
}

forward_pass <- function(model, xreg=NULL, types = c("AO", "LS"), tcritical=3.88){
    y = model$x
    n = length(y)
    
    max_tvalue = 0
    ind_max_tvalue = NULL
    type_max_tvalue = NULL
    for(outliertype in types){
        for(i in 1:n){
            foo = outlier_t_value(model, i,xreg, outliertype )
            tvalue = foo$tvalue
            time_point = foo$time_point
            if(tvalue > tcritical){
                print(paste(time_point, tvalue))
                if(tvalue>max_tvalue){
                    max_tvalue = tvalue
                    ind_max_tvalue = i
                    type_max_tvalue = outliertype
                }
            }
        }
    }
    if(is.null(type_max_tvalue)){
        return(NULL)
    }
    if(type_max_tvalue == "AO"){
        return(AO(y, ind_max_tvalue))
    }
    if(type_max_tvalue == "LS"){
        return(LS(y, ind_max_tvalue))
    }
    return(NULL)
}

backward_pass <-  function(model, xreg=NULL, tcritical=3.88){
    y = model$x
    X = xreg
    nreg = ifelse(is.null(xreg) ,  0, dim(xreg)[2])
    ncoeff = sum(c(eval(model$call$order)[-2] , eval(model$call$seasonal)[-2]))
    model_order = eval(model$call$order)
    model_seasonal = eval(model$call$seasonal)

    fit = Arima(y, order = model_order,
                seasonal = model_seasonal,
                xreg = X,
                include.mean = FALSE, method = 'ML', 
                fixed = c(coef(model)[1:ncoeff], rep(NA, nreg)))
    betas = abs(fit$coef[-c(1:ncoeff)])
    se_betas = sqrt(diag(fit$var.coef))
    tvalues  = betas/se_betas
    
    min_tvalue = 1e8
    min_ind = NULL
    for(i in 1:length(tvalues)){
        if(grepl("^(AO|LS)(1[6-9][0-9]|20[0-9][0-9])\\.(0[1-9]|1[0-2])$",
                 colnames(X)[i])){
            if(tvalues[i] < tcritical ){
                print(paste(colnames(X)[i], tvalues[i]))
                if( tvalues[i] < min_tvalue ){
                    min_tvalue = tvalues[i]
                    min_ind = i                    
                }
                
            }
        }
    }
    
    return(min_ind)
}



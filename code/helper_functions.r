library(ggplot2)
library(gridExtra)
library(forecast)
selection_criteria <- function(model, adjust_for_log_transform = F) {
    
    # Extract log-likelihood from the model
    loglik <- as.numeric(logLik(model))  # Extract log-likelihood
    
    # Extract original series and frequency from model
    original_series <- model$x
    freq <- frequency(model$x)
    
    # Extract non-seasonal and seasonal orders from model$call
    arima_order <- eval(model$call$order)      # Non-seasonal order (p, d, q)
    seasonal_order <- eval(model$call$seasonal) # Seasonal order (P, D, Q)
    
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
    n_p <- arima_order[1] + arima_order[3] + seasonal_order[1] + seasonal_order[3] + 1  # Add 1 for variance
    
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
    return(list(AIC = aic, BIC = bic, AICc = aicc, HQ = hq))
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
AO <- function(y, t_oa) {
    as.integer(seq_along(y) == t_oa)
}
LS <- function(y, t_ls) {
    as.integer(seq_along(y) >= t_ls)-1
}

get_time <- function(t){
    return(c(floor(t), (t-floor(t))*12))
}
detect_outlier <- function(default_model, t_outlier, type = c("AO", "LS")) {
    y <- default_model$x
    outlier_time <- time(y)[t_outlier]  # Find the time corresponding to the outlier index
    
    # Get the year and month using the get_time function
    time_parts <- get_time(outlier_time)
    year <- time_parts[1]
    month <- sprintf("%02g", time_parts[2])  # Format month to be 2 digits
    
    # Create the name for the new variable based on outlier type and time
    var_name <- paste0(type, year, ".", month)
    
    # Generate the regression variable based on outlier type
    reg_var <- if (type == "AO") AO(y, t_outlier) else LS(y, t_outlier)
    reg_var <- cbind(reg_var) # Assign the dynamically created name to the variable
    colnames(reg_var) <- var_name
    # Fit the RegARIMA model with the new regression variable
    fit <- Arima(y, order = as.numeric(paste0(default_model$call$order)[2:4]), 
                 seasonal =as.numeric(paste0(default_model$call$order)[2:4]), 
                 include.mean = FALSE, xreg = reg_var,
                 method = 'ML', 
                 fixed = c(coef(default_model_fit1), NA))
    
    # Extract the coefficient and standard error for the outlier variable
    coef_outlier <- coef(fit)[var_name]
    se_outlier <- sqrt(diag(fit$var.coef))[var_name]
    
    # Calculate the t-value
    t_value <- coef_outlier / se_outlier
    names(t_value) <- var_name
    return(t_value)
}



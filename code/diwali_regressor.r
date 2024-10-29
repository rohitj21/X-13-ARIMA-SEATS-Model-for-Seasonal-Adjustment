ddates <- read.csv('data/diwali.csv')
# Assuming ddates is a dataframe with columns: year, month, and day for Diwali dates
diwali_w <- function(year, w) {
    
    # Find Diwali month and day for the given year
    m <- ddates$month[which(ddates$year == year)]
    d <- ddates$day[which(ddates$year == year)]
    
    # Create a Date object for Diwali
    diwali_date <- as.Date(paste(year, m, d, sep = "-"))
    
    # Calculate the date 'w' days before Diwali
    start_date <- diwali_date - w
    
    # Initialize counts for each month
    sep_count <- 0
    oct_count <- 0
    nov_count <- 0
    
    # Iterate over each day in the range
    for (i in 1:w) {
        current_date <- diwali_date - i + 1
        current_month <- format(current_date, "%m")
        
        # Count days in each month
        if (current_month == "09") {
            sep_count <- sep_count + 1
        } else if (current_month == "10") {
            oct_count <- oct_count + 1
        } else if (current_month == "11") {
            nov_count <- nov_count + 1
        }
    }
    
    total_days <- w
    return(c("sep" = sep_count / total_days, "oct" = oct_count / total_days, "nov" = nov_count / total_days))
}
P_diwali <- matrix(0,ncol  =3, nrow = 200)
for(i in 1:200){
    P_diwali[i,] <- diwali_w(1899+i, 20)
}
colMeans(P_diwali)

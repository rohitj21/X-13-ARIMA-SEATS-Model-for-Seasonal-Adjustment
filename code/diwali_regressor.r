
# Assuming ddates is a dataframe with columns: year, month, and day for Diwali dates
ddates <- read.csv('data/diwali.csv')
diwali_w <- function(year, w) {

    m <- ddates$month[which(ddates$year == year)]
    d <- ddates$day[which(ddates$year == year)]
    

    diwali_date <- as.Date(paste(year, m, d, sep = "-"))
    sep_count <- 0
    oct_count <- 0
    nov_count <- 0

    for (i in 1:w) {
        current_date <- diwali_date - i+1
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
    
    return(c("sep" = sep_count / w, "oct" = oct_count / w, "nov" = nov_count / w))
}
Dp <- list()
for(w in 1:25){
    Dpw <- matrix(0, 200, 3)
    for(i in 1:200){
        Dpw[i,] <- diwali_w(1899 + i, w)
    }
    print(colMeans(Dpw))
    Dp[[w]] = Dpw
}
diwali_reg_values <- list()
for(w in 1:25){
    diwali_reg_values[[w]] <- t(t(Dp[[w]]) - colMeans(Dp[[w]]))
}
diwali_reg <- list()
for(w in 1:25){
    diwali_reg[[w]] <- matrix(0, 200*12, 3)
    colnames(diwali_reg[[w]]) <- c("year", "month", "reg")
    i = 1
    for(y in 1900:2099){
        for(m in 1:12){
            diwali_reg[[w]][i, 1] <- y
            diwali_reg[[w]][i, 2] <- m
            if(m >=9 & m<=11){
                diwali_reg[[w]][i, 3] <- diwali_reg_values[[w]][y-1899, m - 8]    
            }
            i = i+1
        }
    }
    
    diwali_reg[[w]] = data.frame(diwali_reg[[w]])
}
ll <- diwali_reg[[10]]

ll[which(ll$year >= 2013 & ll$year <= 2025), 3]

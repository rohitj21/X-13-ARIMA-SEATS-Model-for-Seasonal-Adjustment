library(seasonal)

# Seasonal adjustment
m <- seas(AirPassengers)
summary(m)
?seas
tramo_output <- series(m, "ref")

print(tramo_output)
# Install and load DiagrammeR if not already installed
install.packages("DiagrammeR")
library(DiagrammeR)


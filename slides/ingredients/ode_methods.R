library(deSolve)

# Define model function
SIR_model <- function(t, state, parms)
{
    # Get variables
    S <- state["S"]
    I <- state["I"]
    R <- state["R"]
    N <- S + I + R
    # Get parameters
    beta <- parms["beta"]
    gamma <- parms["gamma"]
    # Define differential equations
    dS <- -(beta * I / N) * S
    dI <- (beta * I / N) * S - gamma * I
    dR <- gamma * I
    res <- list(c(dS, dI, dR))
    return (res)
}

# Define parameter values
parms <- c(beta = 0.4, gamma = 0.2)

# Define time to solve equations
times <- seq(from = 0, to = 100, by = 0.1)

# Define initial conditions
N <- 100
I_0 <- 0.1
R_0 <- 0
S_0 <- N - I_0 - R_0
y <- c(S = S_0, I = I_0, R = R_0)

DT = 1
nh = 10

# Solve equations with RK1 (Euler)
y <- c(S = S_0, I = I_0, R = R_0)
output <- c(t = 0, y)
for (t in seq(0, 100-DT, DT))
{
	k <- SIR_model(t, y, parms)[[1]]
	for (dt in seq(DT/nh, DT, DT/nh)) {
		yp <- y + dt * k
		output <- rbind(output, c(t = t + dt, yp))
	}
	y <- yp
}
output1 <- output

# Solve equations with RK2
y <- c(S = S_0, I = I_0, R = R_0)
output <- c(t = 0, y)
for (t in seq(0, 100-DT, DT))
{
	k1 <- SIR_model(t, y, parms)[[1]]
	for (dt in seq(DT/nh, DT, DT/nh)) {
		k2 <- SIR_model(t + dt, y + dt * k1, parms)[[1]]
		yp <- y + dt * (k1 + k2) / 2
		output <- rbind(output, c(t = t + dt, yp))
	}
	y <- yp
}
output2 <- output

# Solve equations with RK3
y <- c(S = S_0, I = I_0, R = R_0)
output <- c(t = 0, y)
for (t in seq(0, 100-DT, DT))
{
	k1 <- SIR_model(t, y, parms)[[1]]
	for (dt in seq(DT/nh, DT, DT/nh)) {
		k2 <- SIR_model(t + 0.5 * dt, y + dt * 0.5 * k1, parms)[[1]]
		k3 <- SIR_model(t + 1.0 * dt, y + dt * -1. * k1 + dt * 2 * k2, parms)[[1]]
		yp <- y + dt * (k1/6 + 2*k2/3 + k3/6)
		output <- rbind(output, c(t = t + dt, yp))
	}
	y <- yp
}
output3 <- output


# Solve equations with deSolve
y <- c(S = S_0, I = I_0, R = R_0)
output4 <- ode(y = y, times = times, func = SIR_model, parms = parms)

pt <- seq(1, nrow(output1), by = nh)
cex <- 1

plot(output4[, 1], output4[, 3], type = "l", col = "darkgrey", lwd = 8, 
	xlab = "Time (days)", ylab = "Infectious people", ylim = c(-8, 25))
legend("topright", bt = "n", col = "darkgrey", lwd = 8, legend = "Exact solution")

plot(output4[, 1], output4[, 3], type = "l", col = "lightgrey", lwd = 8, 
	xlab = "Time (days)", ylab = "Infectious people", ylim = c(-8, 25))
lines(output1[, 1], output1[, 3], type = "l", col = "red", lwd = 2)
points(output1[pt, 1], output1[pt, 3], col = "red", pch = 16, cex = cex)
legend("topright", bt = "n", col = c("lightgrey", "red"), lwd = c(8, 2), 
	legend = c("Exact solution", "Linear approx."))

plot(output4[, 1], output4[, 3], type = "l", col = "lightgrey", lwd = 8, 
	xlab = "Time (days)", ylab = "Infectious people", ylim = c(-8, 25))
lines(output1[, 1], output1[, 3], type = "l", col = "pink", lwd = 2)
points(output1[pt, 1], output1[pt, 3], col = "pink", pch = 16, cex = cex)
lines(output2[, 1], output2[, 3], type = "l", col = "blue", lwd = 2)
points(output2[pt, 1], output2[pt, 3], col = "blue", pch = 16, cex = cex)
legend("topright", bt = "n", col = c("lightgrey", "pink", "blue"), lwd = c(8, 2, 2), 
	legend = c("Exact solution", "Linear approx.", "Quadr. approx."))

plot(output4[, 1], output4[, 3], type = "l", col = "lightgrey", lwd = 8, 
	xlab = "Time (days)", ylab = "Infectious people", ylim = c(-8, 25))
lines(output1[, 1], output1[, 3], type = "l", col = "pink", lwd = 2)
points(output1[pt, 1], output1[pt, 3], col = "pink", pch = 16, cex = cex)
lines(output2[, 1], output2[, 3], type = "l", col = "lightblue", lwd = 2)
points(output2[pt, 1], output2[pt, 3], col = "lightblue", pch = 16, cex = cex)
lines(output3[, 1], output3[, 3], type = "l", col = "purple", lwd = 2)
points(output3[pt, 1], output3[pt, 3], col = "purple", pch = 16, cex = cex)
legend("topright", bt = "n", col = c("lightgrey", "pink", "lightblue", "purple"), lwd = c(8, 2, 2, 2), 
	legend = c("Exact solution", "Linear approx.", "Quadr. approx.", "Cubic approx."))


library(tweak)
library(data.table)
library(ggplot2)

curve(0.4 + 0.3 * cos(x * 2 * pi / 365), 0, 365, n = 1000, ylim = c(0, 1), xaxt = "n",
	xlab = "Month", ylab = expression(beta))
axis(1, at = seq(0, 330, by = 30), labels = month.abb)

curve(ifelse(x > 50 & x < 80, 0.1, 0.8), from = 0, to = 100, n = 1001, ylim = c(0, 1),
	xlab = "Time (days)", ylab = expression(beta))
highlight(x = c(50, 80), col = "red", label = "Lockdown")

mob = fread("~/Documents/uk_covid_data_sensitive/google-mobility-old/Global_Mobility_Report-2021-01-03.csv")
mob = mob[country_region_code == "GB" & sub_region_1 == ""]
plot(mob$date, 1 + mob$transit_stations_percent_change_from_baseline / 100, type = "l", xlim = c(as.Date("2020-01-01"), as.Date("2020-12-31")),
	xlab = NA, ylab = "Transit stations mobility", ylim = c(0, 1.1))




library(deSolve) # For solving systems of ODEs

# Define model function
SIIR_model <- function(times, state, parms)
{
    # Get variables
    S <- state["S"]
    I1 <- state["I1"]
    I2 <- state["I2"]
    R <- state["R"]
    N <- S + I1 + I2 + R
    # Get parameters
    beta <- parms["beta"]
    gamma <- parms["gamma"]
    alpha <- parms["alpha"]
    # Define differential equations
    dS <- -(beta * (I1 + alpha * I2) / N) * S
    dI1 <- (beta * I1 / N) * S - gamma * I1
    dI2 <- (beta * alpha * I2 / N) * S - gamma * I2
    dR <- gamma * (I1 + I2)
    res <- list(c(dS, dI1, dI2, dR))
    return (res)
}

# Define parameter values
parms <- c(beta = 0.25, gamma = 0.15, alpha = 1.5)

# Define time to solve equations
times <- seq(from = 0, to = 100, by = 1)

# Define initial conditions
N <- 500
I_0 <- 1
R_0 <- 0
S_0 <- N - I_0 - R_0
y <- c(S = S_0, I1 = I_0, I2 = 0, R = R_0)

# Solve equations
output_raw <- ode(y = y, times = seq(0, 100, 0.1), func = SIIR_model, parms = parms,
	events = list(data = data.frame(var = "I2", time = 20, value = 1, method = "replace")))

# Convert matrix to data frame for easier manipulation
output <- as.data.frame(output_raw)

# Plot model output
plot(output$time, output$I1, lwd = 2, col = "red", type = "l", xlab = "Time", ylab = "Number", ylim = c(0, 100))
lines(output$time, output$I2, lwd = 2, col = "pink", type = "l")
legend("topright", legend = c("Strain 1", "Strain 2"),
       col = c("red", "pink"), lwd = 2, bty = "n")


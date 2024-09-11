# Load in required packages, deSolve and socialmixr
# If either package is not installed, install using the install.packages() function
library(deSolve)
library(socialmixr)

polymod <- contact_matrix(polymod, countries = "United Kingdom", 
	age.limits = seq(0, 70, by = 10))
mat <- polymod$matrix
n_groups <- nrow(mat)

# vignette("socialmixr")

# Define model function 
SIR_conmat_model <- function(times, state, parms) {
    # Get variables
    S <- state[1:n_groups]
    I <- state[(n_groups + 1):(2 * n_groups)]
    R <- state[(2 * n_groups + 1):(3 * n_groups)]
    N <- S + I + R
    # Extract parameters
    p <- parms[["p"]]
    cm <- parms[["cm"]]
    gamma <- parms[["gamma"]]
    # Build force of infection
    lambda = rep(0, n_groups)
    for (i in 1:n_groups) {
    	lambda[i] <- sum(p * cm[i, ] * I/N)
    }
    lambda <- (p * cm %*% (I/N))[,1]
    dS <- -lambda * S
    dI <- lambda * S - gamma * I
    dR <- gamma * I
    res <- list(c(dS, dI, dR))
    return (res)
}

# Define parameters  
parms <- list(p = 0.05, cm = mat, gamma = 0.2)

# Define time to run model
times <- seq(from = 0, to = 50, by = 1)

# Define initial conditions
N0 <- seq(1000, 600, length.out = n_groups)
I0 <- rep(1, n_groups)
S0 <- N0 - I0
R0 <- rep(0, n_groups)
names(S0) <- paste0("S", 1:n_groups)
names(I0) <- paste0("I", 1:n_groups)
names(R0) <- paste0("R", 1:n_groups)
state <- c(S0, I0, R0)

# Solve equations
output_raw <- ode(y = state, 
                  times = times, 
                  func = SIR_conmat_model, 
                  parms = parms)

# Convert to data frame for easy extraction of columns
output <- as.data.frame(output_raw)

# Plot output
par(mfrow = c(1, 1))
matplot(output$time, output[, 2:ncol(output)], type = "l", xlab = "Time (years)", ylab = "Number of people")
# Cool, but not easy to see... let's split it out

# Make a 2x2 grid
par(mfrow = c(2, 2))
# Plot S, I, R on separate plots
matplot(output$time, output[, (0 * n_groups + 2):(1 * n_groups + 1)], type = "l", 
	xlab = "Time (years)", ylab = "Susceptible", col = hcl.colors(n_groups), lty = 1)
matplot(output$time, output[, (1 * n_groups + 2):(2 * n_groups + 1)], type = "l", 
	xlab = "Time (years)", ylab = "Infectious", col = hcl.colors(n_groups), lty = 1)
matplot(output$time, output[, (2 * n_groups + 2):(3 * n_groups + 1)], type = "l", 
	xlab = "Time (years)", ylab = "Recovered", col = hcl.colors(n_groups), lty = 1)
# Put a legend in the bottom right of the grid
plot.new()
legend("center", legend = rownames(mat), col = hcl.colors(n_groups), lty = 1)

# Put plot grid back to 1x1
par(mfrow = c(1, 1))

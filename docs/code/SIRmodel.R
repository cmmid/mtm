# Solve an ODE equation 
# Arguments are the parameter values
solveODE <- function(parameters, plot.all.results = TRUE, num.row = 1, num.col = 1){

  # Define time to solve equations
  times <- seq(from = 0, to = 50, by = 1)
  
  # Define initial conditions
  N <- 100
  I_0 <- 1
  S_0 <- N - I_0
  R_0 <- 0
  state <- c( S = S_0, I = I_0, R = R_0)
  
  # Solve equations
  output_raw <- ode(y = state, times = times, func = SIR_model, parms = parameters,
                    method = "rk4")
  # Convert to data frame for easy extraction of columns
  output <- as.data.frame(output_raw)
  
  if(plot.all.results == TRUE){
    par(new=TRUE)
    par( mfrow = c(num.row, num.col))
    plot( output$time, output$S, type = "l", col = 4, lwd = 2, ylim = c(0, N),
          xlab = "Time", ylab = "Number of People", main = "")
    lines( output$time, output$I, lwd = 2, col = 2, type = "l")
    legend("topright", legend = c("Susceptible", "Infected"),
           lty = rep(1, 2), col = c(4, 2), lwd = 2, bty = "n")
  } 
  
  return(max(output[,"I"]))
}

### Define the SIR model
SIR_model <- function(times, state, parms){
  ## Define variables
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S + I + R
  
  # Extract parameters
  beta <- parms["beta"]
  gamma <- parms["gamma"]
  
  # Define differential equations
  dS <- - (beta * S * I) / N
  dI <- (beta * S * I) / N - gamma * I
  dR <- gamma * I
  res <- list(c(dS, dI, dR))
  return(res)
}


# Solve an ODE equation 
# Arguments are the parameter values
solveODE_2 <- function(parameters){
  
  # Define time to solve equations
  times <- seq(from = 0, to = 50, by = 1)
  
  # Define initial conditions
  N <- 100
  I_0 <- 1
  S_0 <- N - I_0
  R_0 <- 0
  state <- c( S = S_0, I = I_0, R = R_0)
  
  # Solve equations
  output_raw <- ode(y = state, times = times, func = SIR_model_2, parms = parameters,
                    method = "rk4")
  # Convert to data frame for easy extraction of columns
  output <- as.data.frame(output_raw)
  return(max(output[,"I"]))
}

### Define the SIR model
SIR_model_2 <- function(times, state, parms){
  ## Define variables
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S + I + R
  
  # Extract parameters
  R0 = parms["R0"]
  gamma <- parms["gamma"]
  beta = gamma * R0
  
  # Define differential equations
  dS <- - (beta * S * I) / N
  dI <- (beta * S * I) / N - gamma * I
  dR <- gamma * I
  res <- list(c(dS, dI, dR))
  return(res)
}


# Function to plot the prevalence as gamma changes with uncertainty around R0

MCplot <- function(data){
  ggplot2::ggplot(data = data) + 
    geom_violin(aes(x = factor(1/gamma), y = max.prev)) + 
    xlab("Infectious duration (days)") + 
    ylab("Maximum prevalence") + 
    title("Monte Carlo Sampling over R0")
  
}

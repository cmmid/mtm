library(ggplot2) ## for plotting
library(data.table) ## for manipulation of results

## Function SIR_gillespie.
## This takes three arguments:
## - init_state: the initial state
##   (a named vector containing the number in S, I and R)
## - parms: the parameters
##   (a named vector containing the rates beta and gamma)
## - tf: the end time
SIR_gillespie <- function(init_state, parms, tf) {

  time <- 0 ## initialise time to 0

  ## assign parameters to easy-access variables
  beta <- parms["beta"]
  gamma <- parms["gamma"]

  ## assign states to easy-access variables
  S <- init_state["S"]
  I <- init_state["I"]
  R <- init_state["R"]
  N <- S + I + R

  ## create results data frame
  results_df <- data.frame(time=0, t(init_state))

  ## loop until end time is reached
  while (time < tf) {
    ## update current rates
    rates <- c()
    rates["infection"] <- beta * S * I / N
    rates["recovery"] <- gamma * I

    if (sum(rates) > 0) { ## check if any event can happen
      ## time of next event
      time <- time + rexp(n=1, rate=sum(rates))
      ## check if next event is supposed to happen before end time
      if (time <= tf) {
        ## generate cumulative sum of rates, to determine the type of the next
        ## event
        cumulative_rates <- cumsum(rates)
        ## determine type of next event
        type <- runif(n=1, min=0, max=sum(rates))
        if (type < cumulative_rates["infection"]) {
          ## infection
          S <- S - 1
          I <- I + 1
        } else if (type < cumulative_rates["recovery"]){
          ## recovery
          I <- I - 1
          R <- R + 1
        }
      } else { ## next event happens after end time
        time <- tf
      }
    } else { ## no event can happen - go straight to end time
      time <- tf
    }
    ## add new row to results data frame
    results_df <- rbind(results_df, c(time=time, S=S, I=I, R=R))
  }
  ## return results data frame
  return(results_df)
}

init.values <- c(S=249, I=1, R=0) ## initial state
parms <- c(beta=1, gamma=0.5) ## parameter vector
tmax <- 20 ## end time

## run Gillespie simulation
r <- SIR_gillespie(init_state=init.values, parms=parms, tf=tmax)

## Plot the result
ggplot(r) +
	geom_line(aes(time, I, colour = "I"))

## Re-run this from the line `r <- SIR_gillespie(...)` a few times to 
## convince yourself the output is different every time.

## Run multiple simulation runs and plot a few of them
nsim <- 100 ## number of trial simulations

## We store the simulations in a data frame, traj, which
## contains the results from multiple simulation runs and an additional column
## that represents the simulation index
trajlist <- list()

for (i in 1:nsim) {
	trajlist[[i]] <- SIR_gillespie(init.values, parms, tmax)
}

## convert to long data.table
traj <- rbindlist(trajlist, idcol = "run")
traj <- melt(traj, id.vars = c("run", "time"))

## Next, plot the multiple simulation runs
ggplot(traj, aes(x = time, y = value, group = run, color = variable)) +
    geom_line() +
    facet_wrap(~variable)

## Plot the distribution of overall outbreak sizes.

## create data frame of outbreak sizes
outbreak_sizes <- traj %>%
  filter(compartment=="R") %>%
  group_by(i) %>%
  filter(time==max(time)) %>%
  select(i, size=value)

## plot it as a histogram
ggplot(outbreak_sizes, aes(x=size)) +
  geom_histogram()

## determine number of large outbreaks (defined as larger than 50)
outbreak_sizes %>%
  filter(size > 50) %>%
  nrow

## Extract the values of the trajectory at integer time points
timeTraj <- mlr %>%
  group_by(i, compartment) %>%
  summarise(traj=list(data.frame(
              time=seq(0, tmax, by=0.1),
              value=approx(x=time, y=value, xout=seq(0, tmax, by=0.1),
                           method="constant")$y))) %>%
  unnest(traj)

## Calculate summary trajectory with mean & sd of infectious people over time
sumTraj <- timeTraj %>%
  filter(compartment=="I") %>%
  group_by(time) %>%
  summarise(mean=mean(value),
            sd=sd(value))

## plot
ggplot(sumTraj, aes(x=time, y=mean, ymin=pmax(0, mean-sd), ymax=mean+sd)) +
  geom_line() +
  geom_ribbon(alpha=0.3)

## Only consider trajectories that have not gone extinct:
sumTrajAll <- sumTraj %>%
  mutate(trajectories="all")

sumTrajGr0 <- timeTraj %>%
  filter(compartment=="I", value > 0) %>%
  group_by(time) %>%
  summarise(mean=mean(value),
            sd=sd(value)) %>%
  mutate(trajectories="greater_than_zero")

iTraj <- bind_rows(sumTrajAll, sumTrajGr0)

## plot
ggplot(iTraj, aes(x=time, y=mean, ymin=pmax(0, mean-sd), ymax=mean+sd,
                  colour=trajectories, fill=trajectories)) +
  geom_line() +
  geom_ribbon(alpha=0.3) +
  scale_color_brewer(palette="Set1")

## We now compare the two averages to the deterministic trajectory

## Model function (from ODE practical)
library(deSolve)
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
  res <- list(c(dS, dI, dR ))
  return(res)
}

ode_output_raw <-
  ode(y = init.values, times = seq(0, tmax), func = SIR_model, parms = parms,
      method = "rk4")
## Convert to data frame for easy extraction of columns
ode_output <- as.data.frame(ode_output_raw)

## Combine into one big data frame
allTraj <- ode_output %>%
  gather(compartment, value, 2:ncol(.)) %>% ## convert to long format
  filter(compartment=="I") %>%
  rename(mean=value) %>% ## in deterministic, mean=value
  mutate(trajectories="deterministic", ## label trajectories
         sd=0) %>% ## in deterministic, sd=0
  bind_rows(iTraj)

## plot
ggplot(allTraj, aes(x=time, y=mean, colour=trajectories)) +
  geom_line() +
  scale_color_brewer(palette="Set1")
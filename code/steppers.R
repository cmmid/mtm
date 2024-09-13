library(deSolve)
library(tweak)

methods = list(
	"Euler" = "euler",
	"Runge-Kutta 2" = rkMethod("rk2"),
	"Runge-Kutta 4" = "rk4",
	"LSODA" = "lsoda"
)

spring = function(t, y, p)
{
	return (list(c(y[2], -p[1] * y[1])))
}

tweak({
	t = seq(0, 10, by = `∆t`)
	sol = ode(func = spring, y = c(1, 0), times = t, 
		parms = k, method = methods[[method]])
	curve(cos(sqrt(k) * x), from = 0, to = 10, n = 1001, col = "lightblue", lwd = 8,
		ylim = range(c(1, -1, sol[, 2])))
	lines(sol[, 1:2], col = "red")
	legend("topright", legend = c("Exact", method), col = c("lightblue", "red"), lwd = c(8, 1))
	},
	k = c(1, 0, 20),
	`∆t` = c(0.0001, 1),
	method = names(methods)
)



grav = function(t, y, p)
{
	return (list(c(y[2], p[1])))
}

tweak({
	t = seq(0, 10, by = `∆t`)
	sol = ode(func = grav, y = c(0, 0), times = t, 
		parms = a, method = methods[[method]])
	curve(0.5 * a * x^2, from = 0, to = 10, n = 1001, col = "lightblue", lwd = 8, 
		ylim = range(c(0, 0.5 * a * 100, sol[, 2])))
	lines(sol[, 1:2], col = "red")
	legend("topright", legend = c("Exact", method), col = c("lightblue", "red"), lwd = c(8, 1))
	},
	a = c(-9.81, -20, 0),
	`∆t` = c(0.0001, 1),
	method = names(methods)
)

###############################################################################
# This file simulates the host-parasite dynamics over a range of transmission 
# rates and rate of parasite-induced stopping, capturing migratory stalling.
# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>
# Date: April 7, 2020
###############################################################################

library(animation)
library(fBasics) # for Heaviside function
library(doParallel)

######################################################################################
# Functions
######################################################################################

# This code uses a simplification of the model with:
# 1) No natural host mortality (mu = 0)
# 2) YES host stopping/starting (omega = 0, theta > 0)
# 3) No parasite-induced host mortality

#-------------------------------------------------------------------------------
# Function of *time* derivatives
#-------------------------------------------------------------------------------
partial_t.stall <- function(V, params){ 
	
	# Calculate time derivatives
	d.N <- params['theta']*V[5,]*V[4,]
	
	d.m <- params['lambda']*V[7,] - params['sigma']*V[2,] + (V[4,]/V[1,])*params['theta']*V[5,]*(V[6,] + V[5,] - V[2,])
	
	d.A <- (1 - V[3,])*(params['lambda']*V[7,]/V[2,] + params['sigma']) + (V[4,]*V[5,])/(V[1,]*V[2,])*params['theta']*(V[6,]*(3*V[5,] + 2*V[6,] - 1 - V[3,] - 2*V[2,]) + (V[5,] - V[2,])^2 - V[3,]*V[5,])
	
	d.N_hat<- - params['theta']*V[5,]*V[4,]
	
	d.m_hat<- params['lambda']*V[7,] - V[5,]*(params['sigma'] + params['theta']*V[6,])
	
	d.A_hat <- (1 - V[6,])*(params['lambda']*V[7,]/V[5,] + params['sigma'] + V[6,]*params['theta'])
	
	d.L <- params['kappa']*(V[1,]*V[2,] + V[4,]*V[5,]) - params['mu.L']*V[7,] - params['lambda']*V[7,]*(V[1,] + V[4,])
	
	return(rbind(d.N, d.m, d.A, d.N_hat, d.m_hat, d.A_hat, d.L))
} #end function

#-------------------------------------------------------------------------------
# Function to simulate model on space-time grid
#-------------------------------------------------------------------------------
sim <- function(params, grid, init, partial_t){
	# grid is a named list of: xmax, dx, nx, x, tmax, dt, nt, t
	
	# 1) Set up matrices to store solutions:
	V <- matrix(rep(NA, dim(init)[1] * (grid$nx + 3) * (grid$nt + 1)))
	dim(V) <- c(dim(init)[1], (grid$nx + 3), (grid$nt + 1))
	V[, , 1] <- init
		
	# Advection speed for each variable
	u <- as.numeric(c(rep(0, 3), rep(params['c.'], 3), 0))
	if(is.whole(u * grid$dt / grid$dx) == FALSE) stop("Step size not integer!")
	
	# 2) Loop through timesteps	
	for (n in 1:grid$nt) { # For each timestep, n
		
		# Calculate boundary conditions
		Vn <- V[, , n]
		Vn[, 1] <- Vn[, 3]
		Vn[, grid$nx + 3] <- Vn[, grid$nx + 1]
		
		# Spatial advection (upstream differencing)
		Vnp1 <- Vn
		Vnp1[, c(2:(grid$nx + 2))] <- Vn[, c(2:(grid$nx + 2))] - (u * grid$dt / grid$dx) * (Vn[, c(2:(grid$nx + 2))] - Vn[, c(1:(grid$nx + 1))])
		
		# Temporal dynamics (4th order Runge Kutta)
		k1 <- partial_t(Vnp1, params)
		k2 <- partial_t(Vnp1 + grid$dt / 2 * k1, params)
		k3 <- partial_t(Vnp1 + grid$dt / 2 * k2, params)
		k4 <- partial_t(Vnp1 + grid$dt * k3, params)
		
		V[, , n + 1] <- Vnp1 + grid$dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)#, 0) # Added in max to avoid negative values (Mar 14, 2019)
	} #end timestep n
	
	dimnames(V) <- list(c("N", "m", "A", "N_hat", "m_hat", "A_hat", "L"), grid$x, c(0,grid$t))
	return(V)
	}

#------------------------------------------------------------------------------
# Base parameters
#------------------------------------------------------------------------------
# I commented out the ones that don't matter ( = 0) for stall simulations.

params <- c(
	#	beta = 0, 			# Host birth rate
	#	mu = 0,	 				# Natural host death rate (other base_params.R = 0)
	# alpha = 0, 	 	  # Per-parasite rate of parasite-induced host death
	#	rho = 0,				# Within-host parasite reproduction
	sigma = 0.01,		  # Parasite death rate
	kappa = 0.03, 	  # Production of parasite larvae (1)
	lambda = 0.0005,  # Successful transmission probability (0.01)
	mu.L = 0.015, 	  # Larval parasite death rate
	#	omega = 0, 			# Starting rate for stationary hosts
	#	gamma = 0, 			# Stopping rate for migrating, parasite-free hosts
	theta = 0, 		    # Per-parasite increase in stopping rate (0.1)
	c. = 25		        # Migration speed (km per year) # 20000 km/yr*1yr/365 days = 55 km/day
)

k.start <- 0.8 # Starting overdispersion parameter (constant)
m.start <- 5 # Starting mean parasite burden (constant)

H0.tot <- 10^4 # Total host population size (integrated over space)

mean.start <- 0
sd.start <- 30

#------------------------------------------------------------------------------
# Grid and initial conditions
#------------------------------------------------------------------------------
# For culling, migration speed is constant at 25 km/day so we can fix initial
# conditions and grid

speed <- as.numeric(params['c.'])

# What is the best grid for each speed?
tmax <- 365
xmax <- speed * tmax + sd.start * 8
xmin <- - sd.start * 8

# First try
nx0 <- 4000
dx0 <- (xmax - xmin)/nx0
dt0 <- dx0/speed

dt <- round(dt0, 2)
dx <- dt * speed
nx <- floor((xmax - xmin)/dx)

# Check: should move a whole number each timestep
speed * dt / dx 

grid <- list(
	# space (x) in kms
	xmin = xmin,
	xmax = xmax,
	dx = dx,
	x = seq(xmin - dx, xmax + dx, dx),
	nx = floor((xmax - xmin)/dx),
	# time (t) in days
	tmax = tmax,
	dt = dt,
	nt = floor(tmax / dt),
	t = seq(dt, tmax, dt)
)

# Gaussian starting distribution
H0.Gaus <- H0.tot * (1 / sqrt(2*pi*(sd.start^2))) * exp(-(grid$x - mean.start)^2 / (2*sd.start^2))

init.stall <- matrix(c(
	rep(0.1 , grid$nx + 3),                        # Stationary host density (N)
	rep(m.start, grid$nx + 3),                    # Stationary parasite burden (m)
	rep((m.start+k.start)/k.start, grid$nx + 3),  # Stationary variance-to-mean ratio (A)
	H0.Gaus,                                      # Mig. host density (N_hat)
	rep(m.start, grid$nx + 3),                    # Mig. Parasite burden (m_hat)
	rep((m.start+k.start)/k.start, grid$nx + 3),  # Mig. variance-to-mean ratio (A_hat)
	H0.Gaus*10),                                  # Stationary larvae
	nrow=7, ncol=grid$nx + 3, byrow=TRUE)

# Figure of initial conditions
xR <- c(1 : (100+min(which(H0.Gaus[100:length(H0.Gaus)]<10^-10))))
quartz(width = 6.3, height=5, pointsize = 10)
par(mfrow=c(2,2), mar=c(3,4,2,1), oma=c(2,0,0,0))

plot(grid$x[xR], init.stall[7, xR], "l", lwd=2, las=1, ylab="Larval density", xlab="", col=grey(0.5), lty=c(1,2))
mtext(side=3, adj=0, line=0.5, "a)")
text(-150, par('usr')[4] - 0.15*(par('usr')[4]-par('usr')[3]), expression(paste(italic(L(x, t[0])))), cex=1.5, col=grey(0.5))

plot(grid$x[xR], init.stall[4, xR], "l", lwd=2, las=1, ylab="Host density", xlab="", col=grey(0.7))
mtext(side=3, adj=0, line=0.5, "b)")
lines(grid$x[xR], init.stall[1, xR], lty=2, lwd=1.5, col='#b30000')
text(-150, par('usr')[4] - 0.15*(par('usr')[4]-par('usr')[3]), expression(paste(italic(hat(H)(x, t[0])))), cex=1.5, col=grey(0.5))
text(-150, par('usr')[4] - 0.30*(par('usr')[4]-par('usr')[3]), expression(paste(italic(H(x, t[0])))), cex=1.5, col='#b30000')

plot(grid$x[xR], init.stall[5, xR], "l", lwd=2, las=1, ylab="Parasite burden", xlab="", col=grey(0.7), lty=c(1,2))
mtext(side=3, adj=0, line=0.5, "c)")
lines(grid$x[xR], init.stall[2, xR], lty=2, lwd=1.5, col='#b30000')
text(-150, par('usr')[4] - 0.15*(par('usr')[4]-par('usr')[3]), expression(paste(italic(hat(P)(x, t[0])))), cex=1.5, col=grey(0.5))
text(-150, par('usr')[4] - 0.30*(par('usr')[4]-par('usr')[3]), expression(paste(italic(P(x, t[0])))), cex=1.5, col='#b30000')

plot(grid$x[xR], init.stall[6, xR], "l", lwd=2, las=1, ylab="Variance-to-mean ratio", xlab="", col=grey(0.7), lty=c(1,2))
mtext(side=3, adj=0, line=0.5, "d)")
lines(grid$x[xR], init.stall[3, xR], lty=2, lwd=1.5, col='#b30000')
# legend("topright", lwd=c(2,1.5), col=c(grey(0.7), '#b30000'), c("Migrating", "Stationary (stalling case only)"), bty="n")
text(-150, par('usr')[4] - 0.15*(par('usr')[4]-par('usr')[3]), expression(paste(italic(hat(A)(x, t[0])))), cex=1.5, col=grey(0.5))
text(-150, par('usr')[4] - 0.30*(par('usr')[4]-par('usr')[3]), expression(paste(italic(A(x, t[0])))), cex=1.5, col='#b30000')

mtext(side=1, outer=TRUE, "Distance along migration (km)")

###############################################################################################
# Migratory stalling as a function of lambda and theta
###############################################################################################

#------------------------------------------------
# Changing uptake rate lambda and parasite-induced stopping theta

lambda.all <- c(0.0001, seq(0.001, 0.05, 0.002))
theta.all <- c(0:40)*10^-4 

# Create vector of all parameter combinations
sensPar.stall <- cbind(theta = rep(theta.all, each=length(lambda.all)), lambda = rep(lambda.all, length(theta.all)))
nPar <- dim(sensPar.stall)[1]

#------------------------------------------------
registerDoParallel(cores=10)

ptime <- system.time({ 
	
	out.stall<-foreach(i=1:nPar) %dopar%{
		params.i <- params
		params.i['lambda'] <- sensPar.stall[i, 'lambda']
		params.i['theta'] <- sensPar.stall[i, 'theta']
		
		V <- sim(params = params.i, grid = grid, init.stall, partial_t = partial_t.stall)
		
		# Calculate mean parasite burden
		m.avg <- matrix(NA, nrow = dim(V)[3], ncol = 2, dimnames = list(NULL, c("stat", "mig")))
		for(j in 1:(dim(V)[3])){ # for each timestep
			m.avg[j, 'stat'] <- (t(V['N',,j]) %*% V['m',,j]) * grid$dx / (sum(V['N', , j]) * grid$dx)
			m.avg[j, 'mig'] <- (t(V['N_hat',,j]) %*% V['m_hat',,j]) * grid$dx / (sum(V['N_hat', , j]) * grid$dx)
		}
		
		N.tot <- cbind(stat = apply(V['N', ,], 2, sum) * grid$dx, mig = apply(V['N_hat', ,], 2, sum) * grid$dx) 
		
		days <- as.numeric(dimnames(V)[[3]])
		ind.peak <- which(m.avg[, 'mig'] == max(m.avg[, 'mig']))
		peak.m <- days[ind.peak]
		start.m <- days[which.min(abs(m.avg[2:length(days), 'mig'] - m.start)) + 1]
		if(start.m < peak.m) start.m <- max(days)
		
		rm(V)
		list(peak.m, start.m, m.avg, N.tot, days)
		
	}
})[3]

ptime/60 # 7 mins

# saveRDS(object = out.stall, file="Workspaces/outStall.rds")
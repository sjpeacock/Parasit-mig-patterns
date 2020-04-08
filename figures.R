# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>

# Code to produce figures for migratory escape, culling, and stalling

colPalette <- read.csv("colourPalette1.csv")

colMet <- list( # List of colours for the four metrics
	c(as.character(colPalette[, 2]), "#FFFFFF"),
	c(as.character(colPalette[, 3]), "#FFFFFF"),
	c(as.character(colPalette[, 4]), "#FFFFFF"),
	rev(c(as.character(colPalette[, 5]), "#FFFFFF"))
) 

cex.lab <- 0.75

###############################################################################
# Fig 2: Example of migratory escape dynamics
###############################################################################

load("Workspaces/Fig2_ExampleSim.RData")

quartz(width = 3.42, height = 3, pointsize=8)
layout(matrix(c(1:3, rep(4, 3)), nrow=2, byrow=TRUE))
par(mar=c(3,0,1,0), oma=c(0,3,0,3), bg="white", cex = 1, cex.main = 1, cex.axis = 0.75, mgp = c(2.5, 0.4, 0))

for(d in c(1, 80, 360)){
	i <- findInterval(d, days)
	x <- which(V.eg[1, , i] > 10^-10)
	
	# Host density
	
	plot(dist[x], V.eg[1, x, i], "l", yaxt="n", ylab="", col=grey(0.6), xlab="", xaxt="n")
	if(d==1){
		axis(side = 1, at = c(-100, 0, 100), tck = -0.05)
	}
	if(d==80) axis(side=1, at = c(3800, 4000, 4200))
	if(d==360){
		axis(side = 1, at= c(17900, 18100), tck = -0.05)
		axis(side = 4, las=1, tck = -0.05)
		mtext(side = 4, line=1.5, "Host density", col=grey(0.6))
	}
	if(d == 80) mtext(side = 1, "Distance migrated (km)", line=1.5)
	# Parasite burden
	par(new=TRUE)
	plot(dist[x], V.eg[2, x, i], "l", col=1, yaxt="n", xaxt="n", ylim = c(0, 12), lwd=1.5, xlab="")

		if(d == 1){
		axis(side=2, las=1, tck = -0.05)
		mtext(side = 2, line=1.5, "Parasite burden", col=1)
	}
	
	H.scaled <- V.eg[1, x, i]*12/max(V.eg[1, x, i])
	y <- pmin(V.eg[2, x, i], H.scaled)
	abline( h = m.start, col=1, lty=3)
	text(dist[x[length(x)]], 11.5, paste("day", d), pos=2, cex=0.75)
	
}

#-------------
plot(days, m.avg, "n", las = 1, xlab="", ylab="", ylim=c(2, 9), tck = -0.05)
mtext(side = 2, line=1.5, "Mean parasite burden")
mtext(side = 1, line=1.5, "Duration of migration (days)")
abline(h = m.start, lty=3)
u <- par('usr')

abline(v = c(peak.m, start.m), col=c(colMet[[1]][1], colMet[[2]][1]))
text(c(peak.m, start.m), (u[4]) - 0.05*(u[4] - u[3]), pos=3, xpd=NA, c(expression(paste(italic(t), "(", bar(italic(P))[max], ")", sep="")), expression(paste(italic(t), "(", bar(italic(P))[0], ")", sep=""))), col=c(colMet[[1]][1], colMet[[2]][1]), cex=1)

abline(h=m.avg[ind.peak], col=colMet[[3]][1])
text(u[2], m.avg[ind.peak],  expression(bar(italic(P))[max]), col=colMet[[3]][1], pos=4, xpd=NA, cex=1)
lines(days, m.avg, lwd=1.5)

text(u[1] - 0.1*(u[2]-u[1]), (u[4]) + 0.15*(u[4] - u[3]), "B", xpd=NA)
text(u[1] - 0.1*(u[2]-u[1]), (u[4]) + 1.6*(u[4] - u[3]), "A", xpd=NA)


###############################################################################
# Figure S1: Initial conditions
###############################################################################
quartz(width = 4.33, height = 3, pointsize=8)
par(mfrow=c(2,2), oma = c(2,0,0,0), mar = c(1,4,1,1), bg="white", cex = 1, cex.main = 1, cex.axis = 0.75, mgp = c(2.5, 0.8, 0))

x.ind <- which(V.eg[1,,1] > 10^-10)

for(i in 1:4){
	plot(dist[x.ind], V.eg[c(4,1,2,3)[i], x.ind, 1], "l", las=1, xlab='', ylab=c("Larval density", "Host density", "Parasite burden", "Variance-to-mean ratio")[i])
	mtext(side = 3, adj=0, line = -1.2, paste("  ", LETTERS[i]))
	if(i == 2) abline(h = 0, lty = 3, col = grey(0.7))
	if(i == 3) abline(h = 5, lty = 3, col = grey(0.7))
	if(i == 4) abline(h = V.eg[3,1,1], lty = 3, col = grey(0.7))
	
}
mtext(side = 1, outer=TRUE, "Distance along migration (km)", line=1)
legend("topright", lty = c(1,3), col=c(1, grey(0.7)), c("Moving", "Stationary"), bty="n")

###############################################################################
# Figure 3: Migratory escape and recovery
###############################################################################

#------------------------------------------------------------------------------
# Base parameters
#------------------------------------------------------------------------------
# I commented out the ones that don't matter ( = 0) for escape simulations.

params <- c(
	#	beta=0, 			# Host birth rate
	#	mu=0,	 				# Natural host death rate (other base_params.R = 0)
	#	alpha=0, 	 		# Per-parasite rate of parasite-induced host death
	#	rho=0,				# Within-host parasite reproduction
	sigma=0.01,		# Parasite death rate
	kappa=0.03, 	# Production of parasite larvae (1)
	lambda=0.0005,# Successful transmission probability (0.01)
	mu.L=0.015, 	# Larval parasite death rate
	#	omega=0, 			# Starting rate for stationary hosts
	#	gamma=0, 			# Stopping rate for migrating, parasite-free hosts
	#	theta=0, 		  # Per-parasite increase in stopping rate (0.1)
	c.=25		      # Migration speed (km per year) # 20000 km/yr*1yr/365 days = 55 km/day
)

# params['sigma'] <- 0.05
k.start <- 0.8 # Starting overdispersion parameter (constant)
m.start <- 5 # Starting mean parasite burden (constant)

H0.tot <- 10^4 # Total host population size (integrated over space)

mean.start <- 0
sd.start <- 30

#------------------------------------------------------------------------------
# Initial conditions
#------------------------------------------------------------------------------

# Even grid for final simulations
speed.all <- c(1, seq(5, 100, 5))
lambda.all <- c(0.0001, seq(0.001, 0.022, 0.001))
rat <- data.frame(speed = rep(speed.all, each = length(lambda.all)), lambda = rep(lambda.all, length(speed.all)))

nSim <- dim(rat)[1]

eg <- c(
	which(round(rat[,'speed']) == 20 & round(rat[,'lambda'], 3) == 0.005), 
	which(round(rat[,'speed']) == 80 & round(rat[,'lambda'], 3) == 0.005))

params.eg<-c(sigma = 0.01, kappa = 0.03, lambda = 0.004, mu.L = 0.015, c. = 50)

# What is the best grid for each speed?
tmax <- 365
xmax <- speed.all * tmax + sd.start * 8
xmin <- - sd.start * 8

# First try
nx0 <- 4000
dx0 <- (xmax - xmin)/nx0
dt0 <- dx0/speed.all

dt <- round(dt0, 2)
dx <- dt * speed.all
nx <- floor((xmax - xmin)/dx)

speed.all * dt / dx

speed.ind <- match(rat$speed, unique(rat$speed))

# Over four different values of sigma (within-host parasite mortality)
sigma.all <- c(0.001, 0.005, 0.01, 0.02, 0.05)

#------------------------------------------------------------------------------
# Load escape data
#------------------------------------------------------------------------------
out.sigma <- readRDS("Workspaces/outSigma_20190325.rds")

#------------------------------------------------------------------------------
# Summarize output
#------------------------------------------------------------------------------

peak.m <- matrix(NA, nrow=nSim, ncol=length(sigma.all)) # timing of peak m.avg within 365 days
start.m <- matrix(NA, nrow=nSim, ncol=length(sigma.all)) # timing of m.avg = m.start (if within 365 days)
max.m <-matrix(NA, nrow=nSim, ncol=length(sigma.all)) # maximum m.avg over 365 days
for(s in 1:length(sigma.all)){
	for(i in 1:nSim){
		peak.m[i, s] <- out.sigma[[s]][[i]][[1]]
		start.m[i, s] <- out.sigma[[s]][[i]][[2]]
		max.m[i, s] <- max(out.sigma[[s]][[i]][[3]])
	}}

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

quartz(width = 4.5, height = 3.8, pointsize = 8)
# pdf(file = "Fig3.pdf", width = 4.5, height = 3.8, pointsize = 8)
layout(matrix(c(1:((length(sigma.all) + 1)*3)), nrow=3, byrow=TRUE))
par(oma=c(3,4,2,0), bg="white", cex = 1, cex.main = 1, cex.axis = 0.75, mgp = c(2.5, 0.4, 0), tck = -0.05)

#------------------------------------------------------------------------------
# Days until peak mean parasite burden
#------------------------------------------------------------------------------
par(mar=c(0,0,2,0))

for(s in 1:length(sigma.all)){
	plot(rat[,'speed'], rat[,'lambda'], "n", yaxs="i", xaxs="i", yaxt="n", xaxt="n", xlab="")
	axis(side = 1, at = c(25, 50, 75), labels = FALSE)
	
	if(s==1){
		axis(side = 2, las = 1)
		mtext(side=3, expression(paste("Days until peak mean parasite burden - ", italic(t), "(", italic(bar(P))[max], ")")), line=0, adj=0, col=colMet[[1]][1])
		}
	
	image(matrix(peak.m[, s], nrow=length(speed.all), ncol=length(lambda.all), byrow=TRUE), x=speed.all, y=lambda.all, col=colorRampPalette(colMet[[1]])(n=12), ylim=c(0, 0.022), breaks=seq(0, 360, 30), add=TRUE)
	
	mtext(side = 3, line = -0.8, adj = 0, paste("", LETTERS[s]), cex = cex.lab, col = c(1,1,1,1,"white")[s])

	text(10, 0.03, substitute(paste(mu[P] == S), list(S = sigma.all[s])), pos=4, xpd=NA)
	if(s==3) mtext(side=3, line=2.7, expression(paste("Parasite natural death (", italic(mu[P]), " ,", d^-1,")")))
	
	if(sigma.all[s] == params.eg['sigma']) points(params.eg['c.'], params.eg['lambda'], pch=8, lwd=0.8, col="white")
} # end s

par(mar=c(0,1,2,2))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 12), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[1]])(n=12)[i], border=NA)
}
axis(side=4, at=seq(0, 12, 2), labels=c(seq(0, 300, 60), ">360"), las=1)

#------------------------------------------------------------------------------
# Days until m.start
#------------------------------------------------------------------------------
par(mar=c(0,0,2,0))

for(s in 1:length(sigma.all)){
	plot(rat[,'speed'], rat[,'lambda'], "n", yaxs="i", xaxs="i", yaxt="n", xaxt="n", xlab="")
	if(s==1){
		axis(side = 2, las=1)
		mtext(side=3, expression(paste("Days until mean parasite burden declines to initial - ", italic(t), "(", italic(bar(P))[0], ")")), line=0, adj=0, col=colMet[[2]][1])
	}
	
	axis(side = 1, at = c(25, 50, 75), labels = FALSE)
	
	image(matrix(start.m[, s], nrow=length(speed.all), ncol=length(lambda.all), byrow=TRUE), x=speed.all, y=lambda.all, col=colorRampPalette(colMet[[2]])(n=12), ylim=c(0, 0.022), breaks=seq(0, 360, 30), add=TRUE)
	
	mtext(side = 3, line = -0.8, adj = 0, paste("", LETTERS[s+5]), cex = cex.lab, col = c(1,1,1,1,"white")[s])
	
	if(sigma.all[s] == params.eg['sigma']) points(params.eg['c.'], params.eg['lambda'], pch=8, lwd=0.8)
	
} # end s

par(mar=c(0,1,2,2))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 12), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[2]])(n=12)[i], border=NA)
}
axis(side=4, at=seq(0, 12, 2), labels=c(seq(0, 300, 60), ">360"), las=1)


#------------------------------------------------------------------------------
# Maximum parasite burden
#------------------------------------------------------------------------------

par(mar=c(0,0,2,0))
for(s in 1:length(sigma.all)){
	plot(rat[,'speed'], rat[,'lambda'], "n", yaxs="i", xaxs="i", yaxt="n", xaxt="n", xlab="")
	if(s==1){
		axis(side = 2, las=1)
		mtext(side=3, expression(paste("Peak mean parasite burden - ", italic(bar(P))[max])), line=0, adj=0, col=colMet[[3]][1])
	}
	
	axis(side = 1, at = c(25, 50, 75), labels = FALSE)
	axis(side = 1, at = c(20, 50, 80), labels = c(25, 50, 75), tck = 0)
	
	
	u <- par('usr')
	# polygon(x = c(u[1], u[2], u[2], u[1]), y = c(0, 0, 0.03, 0.03), col=tail(colMet[[3]], 1), border=NA)
	
	image(matrix(max.m[, s], nrow=length(speed.all), ncol=length(lambda.all), byrow=TRUE), x=speed.all, y=lambda.all, col=colorRampPalette(colMet[[3]])(n=10), yaxt="n", bty="o", ylim=c(0, 0.022), breaks=seq(0, 100, 10), add=TRUE)
	
	mtext(side = 3, line = -0.8, adj = 0, paste("", LETTERS[s+10]), cex = cex.lab, col = 1)
	
	if(sigma.all[s] == params.eg['sigma']) points(params.eg['c.'], params.eg['lambda'], pch=8, lwd=0.8, col=1)
} # end s

par(mar=c(0,1,2,2))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 10), xlim=c(0,1))
# polygon(x=c(0,1,1,0), y = c(12, 12, 13, 13), col=tail(colMet[[3]], 1), border=NA)
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[3]])(n=10)[i], border=NA)
}
axis(side=4, at=seq(0, 10, 2), labels=c(seq(0, 80, 20), ">100"), las=1)

#------------------------------------------------------------------------------
# Axis lables
#------------------------------------------------------------------------------

mtext(side=1, outer=TRUE, expression(paste("Migration speed (", italic(c), ", km ", d^-1,")")), line=2)
mtext(side=2, outer=TRUE, expression(paste("Transmission rate (",beta, ", ", d^-1, ")")), line=2.5)

###############################################################################
# Figure 4: Combined migratory culling and migratory stalling
###############################################################################

#------------------------------------------------------------------------------
# Load culling data
#------------------------------------------------------------------------------
lambda.cull <- c(0.0001, seq(0.001, 0.05, 0.001))
alpha.all<-seq(0, 0.003, 0.0001) 
sensPar.culling<-cbind(alpha = rep(alpha.all, each=length(lambda.cull)), lambda = rep(lambda.cull, length(alpha.all)))
nPar.cull <- dim(sensPar.culling)[1]

out.cull <- readRDS("Workspaces/outCull_1581_20190421.rds")



#------------------------------------------------------------------------------
# Load stalling data
#------------------------------------------------------------------------------
lambda.stall <- c(0.0001, seq(0.001, 0.05, 0.002))
theta.all <- c(0:40)*10^-4 
sensPar.stall <- cbind(theta = rep(theta.all, each=length(lambda.stall)), lambda = rep(lambda.stall, length(theta.all)))
nPar.stall <- dim(sensPar.stall)[1]

theta.eg <- c(0, 0.001, 0.001)
lambda.eg <- c(0.031, 0.031, 0.0001)

out.stall <- readRDS("Workspaces/outStall_20190412.rds")

#------------------------------------------------------------------------------
# Create list to store output from both culling and stalling
#------------------------------------------------------------------------------

peak.m <- list(cull = numeric(nPar.cull), stall = numeric(nPar.stall)) # timing of peak m.avg within 365 days
start.m <- list(cull = numeric(nPar.cull), stall = numeric(nPar.stall)) # timing of m.avg = m.start (if within 365 days)
max.m <- list(cull = numeric(nPar.cull), stall = numeric(nPar.stall))# maximum m.avg over 365 days
ppnMig.start.m <- list(cull = numeric(nPar.cull), stall = numeric(nPar.stall)) # maximum m.avg over 365 days

# Fill in culling
for(i in 1:nPar.cull){
	peak.m[[1]][i] <- out.cull[[i]][[1]]
	start.m[[1]][i] <- out.cull[[i]][[2]]
	max.m[[1]][i] <- max(out.cull[[i]][[3]])
	j <- which(out.cull[[i]][[5]] == out.cull[[i]][[2]])
	if(length(j) > 0) ppnMig.start.m[[1]][i] <- out.cull[[i]][[4]][j] else ppnMig.start.m[[1]][i] <- NA
}


# Fill in stalling
for(i in 1:nPar.stall){
	peak.m[[2]][i] <- out.stall[[i]][[1]]
	start.m[[2]][i] <- out.stall[[i]][[2]]
	max.m[[2]][i] <- max(out.stall[[i]][[3]][, 'mig'])
	j <- which(out.stall[[i]][[5]] == out.stall[[i]][[2]])
	ppnMig.start.m[[2]][i] <- out.stall[[i]][[4]][j, 'mig']/sum(out.stall[[i]][[4]][j, ])
}

ppnMig.start.m[[2]][which(start.m[[2]] == 365)] <- NA
#------------------------------------------------------------------------------
# Figure setup
#------------------------------------------------------------------------------
quartz(width = 3.42, height = 5, pointsize = 8)
# pdf(file = "Fig4.pdf", width = 3.42, height = 5, pointsize = 8)
layout(matrix(c(rep(c(1,2),each=4), rep(3, 2), rep(c(4,5), each=4), rep(6, 2), rep(c(7,8), each=4), rep(9, 2), rep(c(10,11), each=4), rep(12, 2)), nrow=4, byrow=TRUE))

par(oma=c(3,3.2,2,0), bg="white", cex = 1, cex.main = 1, cex.axis = 0.75, mgp = c(2.5, 0.4, 0), tck = -0.05)

#------------------------------------------------------------------------------
# peak.m
#------------------------------------------------------------------------------
par(mar=c(1,0,1,1))

# A - Culling
plot(sensPar.culling[,'alpha'], sensPar.culling[,'lambda'], "n", xlab="", ylab="", xaxs="i", yaxs="i", las=1, xaxt="n", ylim=c(0.0001, 0.0490)) 
axis(side=1, at = seq(0, 0.003, 0.001), labels=FALSE)
image(matrix(peak.m[[1]], nrow=length(alpha.all), ncol=length(lambda.cull), byrow=TRUE), x=alpha.all, y=lambda.cull, col=colorRampPalette(colMet[[1]])(n=12), ylim=c(0, 0.022), breaks=seq(0, 360, 30), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[1]), cex = cex.lab, col = "white")
mtext(side=3, line=2, "Culling")
mtext(side=3, expression(paste("Days until peak mean parasite burden - ", italic(t), "(", italic(bar(P))[max], ")")), line=0, adj=0, col=colMet[[1]][1])

# B = Stalling
plot(sensPar.stall[,'theta'], sensPar.stall[,'lambda'], "n",xlab="", ylab="", xaxs="i", yaxs="i", las=1, yaxt="n", xaxt="n")
axis(side=1, at=c(0, 0.001, 0.002, 0.003, 0.004), labels=FALSE)
axis(side=2, at=c(0.01, 0.02, 0.03, 0.04), labels=FALSE)
image(matrix(peak.m[[2]], nrow=length(theta.all), ncol=length(lambda.stall), byrow=TRUE), x=theta.all, y=lambda.stall, col=colorRampPalette(colMet[[1]])(n=12), ylim=c(0, 0.022), breaks=seq(0, 360, 30), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[2]), cex = cex.lab, col = "white")
mtext(side=3, line=2, "Stalling")
points(theta.eg, lambda.eg, pch=c(1,2,6), lwd=0.8, col=1, xpd=NA)

# Scale bar
par(mar=c(1,0,1,2.5))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 12), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[1]])(n=12)[i], border=NA)
}
axis(side=4, at=seq(0, 12, 2), labels=c(seq(0, 300, 60), ">360"), las=1)

#------------------------------------------------------------------------------
# start.m
#------------------------------------------------------------------------------
par(mar=c(1,0,1,1))

# Culling
plot(sensPar.culling[,'alpha'], sensPar.culling[,'lambda'], "n", xlab="", ylab="", xaxs="i", yaxs="i", las=1, xaxt="n", ylim=c(0.0001, 0.0490)) 
axis(side=1, at = seq(0, 0.003, 0.001), labels=FALSE)
image(matrix(start.m[[1]], nrow=length(alpha.all), ncol=length(lambda.cull), byrow=TRUE), x=alpha.all, y=lambda.cull, col=colorRampPalette(colMet[[2]])(n=12), breaks=seq(0, 360, 30), add=TRUE)
mtext(side=3, expression(paste("Days until mean parasite burden declines to initial - ", italic(t), "(", italic(bar(P))[0], ")")), line=0, adj=0, col=colMet[[2]][1])
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[3]), cex = cex.lab, col = 1)

# Stalling
plot(sensPar.stall[,'theta'], sensPar.stall[,'lambda'], "n",xlab="", ylab="", xaxs="i", yaxs="i", las=1, yaxt="n", xaxt="n")
axis(side=1, at=c(0, 0.001, 0.002, 0.003, 0.004), labels=FALSE)
axis(side=2, at=c(0.01, 0.02, 0.03, 0.04), labels=FALSE)
image(matrix(start.m[[2]], nrow=length(theta.all), ncol=length(lambda.stall), byrow=TRUE), x=theta.all, y=lambda.stall, col=colorRampPalette(colMet[[2]])(n=12), breaks=seq(0, 360, 30), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[4]), cex = cex.lab, col = 1)
points(theta.eg, lambda.eg, pch=c(1,2,6), lwd=0.8, col=1, xpd=NA)

# Scale bar
par(mar=c(1,0,1,2.5))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 12), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[2]])(n=12)[i], border=NA)
}
axis(side=4, at=seq(0, 12, 2), labels=c(seq(0, 300, 60), ">360"), las=1)

#------------------------------------------------------------------------------
# max.m
#------------------------------------------------------------------------------
par(mar=c(1,0,1,1))

# Culling
plot(sensPar.culling[,'alpha'], sensPar.culling[,'lambda'], "n", xlab="", ylab="", xaxs="i", yaxs="i", las=1, xaxt="n", ylim=c(0.0001, 0.0490)) 
axis(side=1, at = seq(0, 0.003, 0.001), labels=FALSE)
u <- par('usr')
# polygon(x = c(u[1], u[2], u[2], u[1]), y = c(u[3], u[3], u[4], u[4]), col=tail(colMet[[3]], 1), border=NA)
image(matrix(max.m[[1]], nrow=length(alpha.all), ncol=length(lambda.cull), byrow=TRUE), x=alpha.all, y=lambda.cull, col=colorRampPalette(colMet[[3]])(n=10), yaxt="n", bty="o", breaks=seq(0, 100, 10), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[5]), cex = cex.lab, col = 1)
mtext(side=3, expression(paste("Peak mean parasite burden - ", italic(bar(P))[max])), line=0, adj=0, col=colMet[[3]][1])#tail(colMet[[3]],1))

# Stalling
plot(sensPar.stall[,'theta'], sensPar.stall[,'lambda'], "n",xlab="", ylab="", xaxs="i", yaxs="i", las=1, yaxt="n", xaxt="n")
axis(side=1, at=c(0, 0.001, 0.002, 0.003, 0.004), labels=FALSE)
axis(side=2, at=c(0.01, 0.02, 0.03, 0.04), labels=FALSE)
u <- par('usr')
# polygon(x = c(u[1], u[2], u[2], u[1]), y = c(u[3], u[3], u[4], u[4]), col=tail(colMet[[3]], 1), border=NA)
image(matrix(max.m[[2]], nrow=length(theta.all), ncol=length(lambda.stall), byrow=TRUE), x=theta.all, y=lambda.stall, col=colorRampPalette(colMet[[3]])(n=10), yaxt="n", bty="o", breaks=seq(0, 100, 10), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[6]), cex = cex.lab, col = 1)
points(theta.eg, lambda.eg, pch=c(1,2,6), lwd=0.8, col=1, xpd=NA)

# Scale bar
par(mar=c(1,0,1,2.5))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 10), xlim=c(0,1))
polygon(x=c(0,1,1,0), y = c(12, 12, 13, 13), col=tail(colMet[[3]], 1), border=NA)
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[3]])(n=10)[i], border=NA)
}
axis(side=4, at=seq(0, 10, 2), labels=c(seq(0, 80, 20), ">100"), las=1)

#------------------------------------------------------------------------------
# ppnMig.start.m
#------------------------------------------------------------------------------
par(mar=c(1,0,1,1))

# Culling
plot(sensPar.culling[,'alpha'], sensPar.culling[,'lambda'], "n", xlab="", ylab="", xaxs="i", yaxs="i", las=1, xaxt="n", ylim=c(0.0001, 0.0490)) 

polygon(x = c(-1, 1, 1, -1), y = c(-1, -1, 1, 1), col = 1)

axis(side=1, at = seq(0, 0.003, 0.001))
image(matrix(ppnMig.start.m[[1]], nrow=length(alpha.all), ncol=length(lambda.cull), byrow=TRUE), x=alpha.all, y=lambda.cull, col=colorRampPalette(colMet[[4]])(n=10), yaxt="n", bty="o", ylim=c(0, 1), breaks=seq(0, 1, 0.1), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[7]), cex = cex.lab, col = 1)
mtext(side=3, expression(paste("Fraction of hosts migrating at ", italic(t), "(", italic(bar(P))[0], ")")), line=0, adj=0, col=tail(colMet[[4]],1))

mtext(side=1, line=1.8, "Parasite-induced")
mtext(side=1, line=3, expression(paste("mortality (", alpha, ", ", d^-1,")")))

# Stalling
plot(sensPar.stall[,'theta'], sensPar.stall[,'lambda'], "n",xlab="", ylab="", xaxs="i", yaxs="i", las=1, yaxt="n", xaxt="n")
polygon(x = c(-1, 1, 1, -1), y = c(-1, -1, 1, 1), col = 1)
axis(side=1, at=c(0, 0.001, 0.002, 0.003, 0.004), labels=FALSE)
axis(side=1, at=c(0.001, 0.003))
axis(side=2, at=c(0.01, 0.02, 0.03, 0.04), labels=FALSE)
image(matrix(ppnMig.start.m[[2]], nrow=length(theta.all), ncol=length(lambda.stall), byrow=TRUE), x=theta.all, y=lambda.stall, col=colorRampPalette(colMet[[4]])(n=10), yaxt="n", bty="o", ylim=c(0, 1), breaks=seq(0, 1, 0.1), add=TRUE)
mtext(side = 3, line = -0.8, adj = 0, paste(" ", LETTERS[8]), cex = cex.lab, col = "white")
points(theta.eg, lambda.eg, pch=c(1,2,6), lwd=0.8, col=1, xpd=NA)

mtext(side=1, line=1.8, "Parasite-induced")
mtext(side=1, line=3, expression(paste("stopping (", theta, ", ", d^-1,")")))

# Scale bar
par(mar=c(1,0,1,2.5))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 10), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(colMet[[4]])(n=10)[i], border=NA)
}
axis(side=4, at=seq(0, 10, 2), labels=seq(0, 1, 0.2), las=1)


mtext(side=2, outer=TRUE, expression(paste("Transmission rate (", beta, ", ", d^-1, ")")), line=2)


###############################################################################
# Figure S2: What's going on with N.tot in Fig. 4d??
###############################################################################
cols <- c("orange",2)
eg1 <- c(
	which(round(sensPar.culling[,'alpha'], 4)==0.0005 & round(sensPar.culling[,'lambda'], 4) == 0.0400),
	which(round(sensPar.culling[,'alpha'], 4)==0.0005 & round(sensPar.culling[,'lambda'], 4) == 0.0150),
	which(round(sensPar.culling[,'alpha'], 4)==0.0005 & round(sensPar.culling[,'lambda'], 4) == 0.0020)
)

eg2 <- c(
	which(round(sensPar.culling[,'alpha'], 4)==0.0015 & round(sensPar.culling[,'lambda'], 4) == 0.0400),
	which(round(sensPar.culling[,'alpha'], 4)==0.0015 & round(sensPar.culling[,'lambda'], 4) == 0.0150),
	which(round(sensPar.culling[,'alpha'], 4)==0.0015 & round(sensPar.culling[,'lambda'], 4) == 0.0020)
)

quartz(width = 7, height = 7, pointsize=10)
# pdf(file = "FigS2.pdf", width = 7, height = 7, pointsize = 10)

#--------------------------------
layout(matrix(c(rep(c(1,1,2), 3), rep(c(3,5,7), 2), rep(c(4,6,8), 2)), nrow=7, byrow=TRUE))
par(oma=c(4,2,0,0))

par(mar=c(6,3,2,1))
plot(sensPar.culling[,'alpha'], sensPar.culling[,'lambda'], "n", xlab="", ylab="", xaxs="i", yaxs="i", las=1, xaxt="n", xlim=c(0, 0.002)) 
axis(side=1, at = seq(0, 0.003, 0.001))
polygon(x = c(-1, 1, 1, -1), y = c(-1, -1, 1, 1), col = 1)

image(matrix(ppnMig.start.m[[1]], nrow=length(alpha.all), ncol=length(lambda.cull), byrow=TRUE), x=alpha.all, y=lambda.cull, col=colorRampPalette(rev(colMet[[4]]))(n=10), yaxt="n", bty="o", ylim=c(0, 1), breaks=seq(0, 1, 0.1), add=TRUE)

mtext(side=1, line=2.5, expression(paste("Parasite-induced mortality (", alpha, ", ", d^-1,")")), cex=cex.lab)
mtext(side=2, expression(paste("Transmission rate (", beta, ", ", d^-1, ")")), cex=cex.lab, line=3)


points(sensPar.culling[,'alpha'][eg1], sensPar.culling[,'lambda'][eg1], pch = c(24,21,25), col=1, lwd = 1.2, bg = "white")
points(sensPar.culling[,'alpha'][eg2], sensPar.culling[,'lambda'][eg2], pch = c(24,21,25), col=grey(0.7), lwd = 1.2, bg="white")
mtext(side=3, adj=0, line=0, expression(paste("A) Fraction of hosts migrating at ", italic(t)==italic(t), "(", italic(bar(P))[0], ")")), col=col4[1], cex=cex.lab)

par(mar=c(6,1,2,15))
plot(1, 1, "n", bty="o", xaxt="n", yaxt="n", ylab="", xlab="", yaxs="i", xaxs="i", ylim=c(0, 10), xlim=c(0,1))
for(i in 1:12){
	polygon(x=c(0,1,1,0), y = c(i-1, i-1, i, i), col=colorRampPalette(rev(colMet[[4]]))(n=10)[i], border=NA)
}
axis(side=4, at=seq(0, 10, 2), labels=seq(0, 1, 0.2), las=1)

text(3, 8, "Low parasite-\ninduced mortality", xpd=NA, font=2, col=1, cex=1.2)
text(3, 6.5, expression(paste(alpha, " = 0.0005")), font = 2, xpd = NA, cex=1.2)
text(3, 4, "High parasite-\ninduced mortality", xpd=NA, cex=1.2, font=2, col=grey(0.7))
text(3, 2.5, expression(alpha == 0.0015), font = 2, xpd=NA, cex=1.2, col=grey(0.7))


xind <- seq(1, length(out.cull[[eg1[1]]][[5]]), 30)
par(mar=c(1,3,3,1))

for(i in 1:3){
	plot(out.cull[[eg1[1]]][[5]], out.cull[[eg1[1]]][[3]], "n", las=1, xlab="", ylab="m.avg")
	u <- par('usr')
	abline(h = 5)
	mtext(side=3, line=1.2, cex=cex.lab, c(
		expression(paste("High transmission (", beta == 0.040, ")", sep="")), 
		expression(paste("Moderate transmission (", beta == 0.015, ")", sep="")), 
		expression(paste("Low transmission (", beta == 0.002, ")", sep="")))[i])
	points(mean(u[1:2]), u[4] + 0.3*(u[4] - u[3]), pch = c(24,21,25)[i], col=1, lwd = 1.2, bg = "white", cex=2, xpd=NA)
	lines(out.cull[[eg1[i]]][[5]][xind], out.cull[[eg1[i]]][[3]][xind], col=1)
	lines(out.cull[[eg2[i]]][[5]][xind], out.cull[[eg2[i]]][[3]][xind], col=grey(0.7))
	
	if(i==1) mtext(side=2, line=3, expression(paste("Average parasite burden (", bar(italic(P)), ")")), cex=cex.lab)
	mtext(side=3, adj=0, line=0, cex=cex.lab, LETTERS[i+1])
	
	if(i==3) legend("topright", lty = 1, col = c(1, grey(0.7)), legend = c(expression(paste(alpha, " = 0.0005")), expression(alpha == 0.0015)), bty="n")
	plot(out.cull[[eg1[1]]][[5]], out.cull[[eg1[1]]][[4]], "n", las=1, xlab="", ylab="N.frac", ylim=c(0,1))
	
	lines(out.cull[[eg1[i]]][[5]][xind], out.cull[[eg1[i]]][[4]][xind],  col=1)
	lines(out.cull[[eg2[i]]][[5]][xind], out.cull[[eg2[i]]][[4]][xind], col=grey(0.7))
	
	segments(x0=out.cull[[eg1[i]]][[2]], x1=out.cull[[eg1[i]]][[2]], y0=1.613743, y1=out.cull[[eg1[i]]][[4]][findInterval(out.cull[[eg1[i]]][[2]], out.cull[[eg1[i]]][[5]])], xpd=NA, col=1)
	segments(x0=out.cull[[eg2[i]]][[2]], x1=out.cull[[eg2[i]]][[2]], y0=1.613743, y1=out.cull[[eg2[i]]][[4]][findInterval(out.cull[[eg2[i]]][[2]], out.cull[[eg2[i]]][[5]])], xpd=NA, col=grey(0.7))
	segments(x0=out.cull[[eg1[i]]][[2]], x1=-5, y0=out.cull[[eg1[i]]][[4]][findInterval(out.cull[[eg1[i]]][[2]], out.cull[[eg1[i]]][[5]])], col=1)
	segments(x0=out.cull[[eg2[i]]][[2]], x1=-5, y0=out.cull[[eg2[i]]][[4]][findInterval(out.cull[[eg2[i]]][[2]], out.cull[[eg2[i]]][[5]])], col=grey(0.7))
	
	if(i==1) mtext(side=2, line=3, "Fraction of hosts migrating", cex=cex.lab)
	mtext(side=1, line=3, "Time (days)", cex=cex.lab)
	mtext(side=3, adj=0, line=0, cex=cex.lab, LETTERS[i+4])
	
}

###############################################################################
# Figure 5: Stalling example
###############################################################################

theta.eg <- c(0, 0.001, 0.001)
lambda.eg <- c(0.031, 0.031, 0.0001)
eg <- numeric(length(theta.eg))
for(i in 1:length(theta.eg)){
	eg[i] <- which(sensPar.stall[,'theta'] == theta.eg[i] & round(sensPar.stall[,'lambda'], 4) == round(lambda.eg[i], 4))}


days <- out.stall[[1]][[5]]
xind <- seq(1, length(days), 30)


#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

quartz(width = 4.5, height = 3.4, pointsize = 8)
# pdf(file = "Fig5.pdf", width = 4.5, height = 3.4, pointsize = 8)
par(mfrow = c(2, 3), mar=c(1,1,0,0), oma=c(1.5, 2, 5, 0.5), bg="white", cex = 1, cex.main = 1, cex.axis = 0.75, mgp = c(2.5, 0.4, 0), tck = -0.03)

for(i in 1:length(eg)){
	
	plot(days[xind], out.stall[[eg[i]]][[3]][xind, 'mig'], "l", xlab="", ylab = "", ylim=c(0, 100), las=1, yaxt="n", xaxt="n")
	axis(side=1, labels=FALSE)
	if(i==1) axis(side=2, las=1) else axis(side=2, labels=FALSE)
	lines(days[xind], out.stall[[eg[i]]][[3]][xind, 'stat'], col=grey(0.6), lty=2)
	mtext(side=3, adj=0, line = -1, paste("  ", LETTERS[i], sep=""))
	if(i==1) mtext(side=2, "Mean parasite burden", line=2)
	
	# text on top
	u <- par('usr')
	tx <- c(0.14, 0.5, 0.4)
	points(mean(days[xind]), u[4] + tx[1]*(u[4] - u[3]), pch = c(21, 24, 25)[i], cex=2, xpd=NA)
	# Stopping
	if(i == 1) text(mean(days[xind]), u[4] + tx[2]*(u[4] - u[3]), expression(paste("No stopping (", theta == 0, ")", sep="")), xpd = NA)
	if(i==2){
		segments(x0=u[1], x1 = 800, y0 =u[4] + tx[2]*(u[4] - u[3]), y1= u[4] + tx[2]*(u[4] - u[3]), xpd = NA, col=colMet[[3]][1], lwd=2)
		segments(x0 = 204, x1 = 536, y0 = u[4] + tx[2]*(u[4] - u[3]), y1 = u[4] + tx[2]*(u[4] - u[3]), col="white", lwd=3, xpd=NA)
		text(mean(c(u[1], 2*u[2])), u[4] + tx[2]*(u[4] - u[3]), expression(paste("Stopping (", theta == 0.001, ")", sep="")), xpd = NA, col=colMet[[3]][1])
	}
	# Transmission
	if(i == 3){
		text(mean(days[xind]), u[4] + tx[3]*(u[4] - u[3]), expression(paste("Low transmission")), xpd = NA)
		text(mean(days[xind]), u[4] + (tx[3] - 0.1)*(u[4] - u[3]), expression(paste("(", beta == 0.0001, ")", sep="")), xpd = NA)
	}
	if(i == 2){
		segments(x0=-444, x1 = 372, y0 =u[4] + (tx[3]-0.08)*(u[4] - u[3]), y1= u[4] + (tx[3]-0.08)*(u[4] - u[3]), xpd = NA, col=colMet[[3]][1], lwd=2)
		segments(x0 = -278, x1 = 202, y0 = u[4] + (tx[3]-0.08)*(u[4] - u[3]), y1 = u[4] + (tx[3]-0.08)*(u[4] - u[3]), col="white", lwd=3, xpd=NA)
		text(mean(c(-444, 372)), u[4] + (tx[3]-0.08)*(u[4] - u[3]), expression(paste("High transmisison (", beta == 0.031, ")", sep="")), xpd = NA, col=colMet[[3]][1])
	}
}
legend("topright", lty=c(1,2), c("migrating", "stationary"), col=c(1, grey(0.6)), bty="n", lwd=1)

for(i in 1:length(eg)){
	plot(days[xind], out.stall[[eg[i]]][[4]][xind, 'mig']/sum(out.stall[[eg[i]]][[4]][1,]), "l", xlab="", ylab = "Mean parasite burden", ylim=c(0, 1), las=1,  yaxt="n", bg="white")
	if(i==1) axis(side=2, las=1) else axis(side=2, labels=FALSE)
	lines(days[xind], out.stall[[eg[i]]][[4]][xind, 'stat']/sum(out.stall[[eg[i]]][[4]][1,]), col=grey(0.6), lty=2)
	if(i==1) mtext(side=2, "Fraction of hosts", line = 2)
	mtext(side=3, adj=0, line = -1, paste("  ", LETTERS[3 + i], sep=""))
}

mtext(side=1, outer=TRUE, "Days migrating", line = 0.5)
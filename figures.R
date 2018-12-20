# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>

# Code to produce figures for migratory escape, culling, and stalling

# Color palettes used
col.pal<-c(red="#ED2124", grey(0.6), grey(0.4), 1)
col.pal2<-c("#f4a582", grey(0.9), grey(0.7), grey(0.5))

###############################################################################
# Migratory escape as a function of lambda/c
###############################################################################

rat.esc <- cbind(rep(c.all, each=length(lambda.all)), rep(lambda.all, length(c.all)))
SpreadCat_escape <- esc.cat

# Which parameter combinations to use as an example?
eg <- c(
	spread=which(rat.esc[,1]==11&rat.esc[,2]==0.0036), 
	spread.here=which(rat.esc[,1]==31&rat.esc[,2]==0.0020), 
	escape.here=which(rat.esc[,1]==45&rat.esc[,2]==0.0008), 
	escape.always=which(rat.esc[,1]==57&rat.esc[,2]==0.0000))

# quartz(width=3.4252, height=5, pointsize=10)
pdf(file="Figures/Escape.pdf", width=3.4252, height=5, pointsize=10)
par(mfrow=c(2,1), mar=c(4,5,2,2), mgp=c(2.5, 1, 0))
plot(time[eg[1],], m.avg.all[eg[1],]/5, "l", ylim=range(m.avg.all[eg,]/5), xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5, yaxt="n")
axis(side=2, at=c(1:3))
points(seq(15, 60, 15), m.avg.all[eg[1],seq(1, 500, length.out=5)[2:5]]/5, pch=1, col=col.pal[1])
mtext(side=2, "Change in parasite burden", line=2.5)
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(time[eg[j],], m.avg.all[eg[j],]/5, col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m.avg.all[eg[j], seq(1, 500, length.out=5)[2:5]]/5, pch=j, col=col.pal[j])
}
mtext(side=3, adj=0, line=0.5, "a)")

image(matrix(esc.cat, nrow=length(c.all), ncol=length(lambda.all), byrow=TRUE), x=c.all, y=lambda.all, xlab=expression(paste("Migration speed (", italic(c), ", km ", d^-1,")")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=col.pal2, yaxt="n", bty="l")
axis(side=2, at=c(0, 0.002, 0.004))
# lines(cec[[1]][,'x'], cec[[1]][,'y'])

points(rat.esc[eg,1], rat.esc[eg,2], pch=1:4, col=col.pal, lwd=1.5)

abline(v=base.params['c.'], lty=3)
abline(h=base.params['lambda'], lty=3)

text(15, 0.0028, "Spread", col=col.pal[1])
text(60, 0.002, "Escape", col=col.pal[4])
mtext(side=3, adj=0, line=0.5, "b)")
dev.off()

###############################################################################
# Migratory culling 
###############################################################################
eg<-c(
	spread=which(rat[,1]==0.0002&rat[,2]==0.0034), 
	spread.here=which(rat[,1]==0.0010&rat[,2]==0.0034), 
	escape.here=which(rat[,1]==0.0036&rat[,2]==0.0034), 
	escape.always=which(rat[,1]==0.0036&rat[,2]==0.0002))

m.frac_culling<-matrix(nrow=500, ncol=length(eg))
N.frac_culling<-matrix(nrow=500, ncol= length(eg))
for(i in 1:length(eg)){
	m.frac_culling[,i]<-output.both[[1]][[eg[i]]][[1]]/5
	N.frac_culling[,i]<-output.both[[1]][[eg[i]]][[2]]
}

pdf(file="Figures/Culling.pdf", width=3.4252, height=7, pointsize=10)
par(mfcol=c(3,1), mar=c(4,5,2,2), mgp=c(2.5, 1, 0))

#----------------------------
# PARASITE
#----------------------------
plot(V4plot[[2]]*dx/base.params['c.'], m.frac_culling[,1], "l", ylim=range(m.frac_culling), xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5)
points(seq(15, 60, 15), m.frac_culling[seq(1, 500, length.out=5)[2:5],1], pch=1, col=col.pal[1])
mtext(side=2, "Change in parasite burden", line=2.5, cex=par('cex'))
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(V4plot[[2]]*dx/base.params['c.'], m.frac_culling[,j], col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), m.frac_culling[seq(1, 500, length.out=5)[2:5],j], pch=j, col=col.pal[j])
}
mtext(side=3, adj=0, line=0.5, "a)", cex=par('cex'))
# mtext(side=3, line=1.5, "Change in parasite burden")

##############
plot(V4plot[[2]]*dx/base.params['c.'], N.frac_culling[,1], "l", xlab="Time in transient phase (days)", bty="l", ylab="", col=col.pal[1], lwd=1.5, yaxt="n", ylim=range(N.frac_culling))
axis(side=2, at=c(0.4, 0.6, 0.8, 1.0))
points(seq(15, 60, 15), N.frac_culling[seq(1, 500, length.out=5)[2:5],1], pch=1, col=col.pal[1])
mtext(side=2, "Change in host population", line=2.5, cex=par('cex'))
abline(h=1, lty=2)

for(j in 2:length(eg)){
	lines(V4plot[[2]]*dx/base.params['c.'], N.frac_culling[,j], col=col.pal[j], lwd=1.5)
	points(seq(15, 60, 15), N.frac_culling[seq(1, 500, length.out=5)[2:5],j], pch=j, col=col.pal[j])
}
mtext(side=3, adj=0, line=0.5, "b)", cex=par('cex'))
# mtext(side=3, line=1.5, "Change in parasite burden")

################
image(matrix(esc.cat[[1]], nrow=length(alpha.all), ncol=length(lambda.all), byrow=TRUE)[which(alpha.all<=0.004),], x=alpha.all[which(alpha.all<=0.004)], y=lambda.all, xlab=expression(paste("Parasite-induced mortality (", alpha, ")")), ylab=expression(paste("Transmission rate (", lambda, ")")), col=col.pal2, bty="n", yaxt="n", xaxt="n")

axis(side=2, at=c(0, 0.002, 0.004))
axis(side=1, at=seq(0, 0.004, 0.001))

# text(rat[,1], rat[,2], esc.cat, cex=0.8, col=c(col.pal['red'], "#1E90FF", '#104E8B', col.pal['green'])[esc.cat])

points(rat[eg,1], rat[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(v=base.params['alpha'], lty=3)
abline(h=base.params['lambda'], lty=3)

text(0.0018, 0.0025, "Culling", col=col.pal['blue'])
mtext(side=3, adj=0, line=0.5, "c)", cex=par('cex'))
dev.off()

###############################################################################
# Migratory culling and stalling: parameter space only
###############################################################################
SpreadCat_culling <- esc.cat[[1]]
HostPop_culling <- N_hat.frac[,500]

barx<-c(0.009, 0.01, 0.01, 0.009)

pdf(file="Figures/Culling-Stalling.pdf", width=7.008, height=5, pointsize=10)
par(mfrow=c(2,2), mar=c(4,4,2,5), mgp=c(2.5, 1, 0), oma=c(0,2,2,0))

#------------------------------------------------------------------------------
# Culling
#------------------------------------------------------------------------------

# Change in parasite burden

image(matrix(SpreadCat_culling, nrow=length(alpha.all), ncol=length(lambda.all), byrow=TRUE), x=alpha.all, y=lambda.all, xlab=expression(paste("Parasite-induced mortality (", alpha, ")")), ylab=expression(paste("Transmission rate (", lambda, ")")), col=col.pal2, bty="n", yaxt="n", xaxt="n")
axis(side=2, at=c(0, 0.002, 0.004))
axis(side=1, at=seq(0, 0.006, 0.002))

abline(v=base.params['alpha'], lty=3)
abline(h=base.params['lambda'], lty=3)

# text(0.003, 0.0025, "Culling")
mtext(side=3, adj=0, line=0.5, "a)", cex=par('cex'))

mtext(side=3, line=2, "Change in parasite burden")
mtext(side=2, line=4.5, "Migratory culling")

# Change in host population

image(matrix(HostPop_culling, nrow=length(alpha.all), ncol=length(lambda.all), byrow=TRUE), x=alpha.all, y=lambda.all, xlab=expression(paste("Parasite-induced mortality (", alpha, ")")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=colorRampPalette(c(1, "white"))(n=21), yaxt="n", bty="l", zlim=c(0, 1), xaxt="n")

axis(side=2, at=c(0, 0.002, 0.004))
axis(side=1, at=seq(0, 0.006, 0.002))
abline(h=base.params['alpha'], lty=3)
abline(v=base.params['lambda'], lty=3)

mtext(side=3, adj=0, line=0.5, "b)", cex=par('cex'))
mtext(side=3, line=2, "Change in moving hosts")


#------------------------------------------------------------------------------
# Stalling
#------------------------------------------------------------------------------

# Change in parasite burden of moving hosts

image(matrix(SpreadCat_stalling, nrow=length(theta.all), ncol=length(lambda.all), byrow=TRUE), x=theta.all, y=lambda.all, xlab=expression(paste("Parasite-induced stopping (", theta, ", parasite",{}^-1,d^-1, ")", sep="")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=col.pal2, yaxt="n", bty="l")
axis(side=2, at=c(0, 0.002, 0.004))

# text(0.004, 0.0025, "Stalling")
mtext(side=3, adj=0, line=0.5, "c)", cex=par('cex'))
mtext(side=2, line=4.5, "Migratory stalling")

# points(rat.stalling[eg,1], rat.stalling[eg,2], pch=1:4, col=col.pal, lwd=1.5)
abline(h=base.params['lambda'], lty=3)
abline(v=base.params['theta'], lty=3)

polygon(x=barx, y = c(0.0112, 0.0112, 0.0112*0.75, 0.0112*0.75), border=NA, col=col.pal2[1], xpd=NA)
text(mean(barx), 0.0098, "Spread", col=col.pal[1], xpd=NA, srt=90)
polygon(x=barx, y = c(0.0112*0.5, 0.0112*0.5, 0.0112*0.75, 0.0112*0.75), border=NA, col=col.pal2[2], xpd=NA)
text(mean(barx), 0.007, "Decl. above", col=col.pal[2], xpd=NA, srt=90)
polygon(x=barx, y = c(0.0112*0.5, 0.0112*0.5, 0.0112*0.25, 0.0112*0.25), border=NA, col=col.pal2[3], xpd=NA)
text(mean(barx), 0.0042, "Decl. below", col=col.pal[3], xpd=NA, srt=90)
polygon(x=barx, y = c(0, 0, 0.0112*0.25, 0.0112*0.25), border=NA, col=col.pal2[4], xpd=NA)
text(mean(barx), 0.0014, "Always decl.", col=col.pal[4], xpd=NA, srt=90)
polygon(x=barx, y=c(0,0,0.0112, 0.0112), col=NA, xpd=NA)


# Change in moving host population size

image(matrix(HostPop_stalling, nrow=length(theta.all), ncol=length(lambda.all), byrow=TRUE), x=theta.all, y=lambda.all, xlab=expression(paste("Parasite-induced stopping (", theta, ", parasite",{}^-1,d^-1, ")", sep="")), ylab=expression(paste("Transmission rate (",lambda, ", ", d^-1, ")")), col=colorRampPalette(c(1, "white"))(n=21), yaxt="n", bty="l", zlim=c(0, 1))

axis(side=2, at=c(0, 0.002, 0.004))
mtext(side=3, adj=0, line=0.5, "d)", cex=par('cex'))

for(i in 1:21) polygon(x=barx, y=c(0.0112*i*(1/21), 0.0112*i*(1/21), 0.0112*(i-1)*(1/21), 0.0112*(i-1)*(1/21)), col=colorRampPalette(c(1, "white"))(n=21)[i], xpd=NA, border=NA)
polygon(x=barx, y=c(0,0,0.0112, 0.0112), col=NA, xpd=NA)
text(mean(barx), 0.0112*(c(1, 6, 11, 16, 21)-0.5)*(1/21), c(0, 0.25, 0.5, 0.75, 1), xpd=NA, col=c(rep("white", 2), rep(1, 3)))


dev.off()

#
# Main function, contains:
# (1) one simulation example;
# (2) call CNV dection functions;
#     (i) One uses constant prior for position parameter k and lambdas;
#     (ii)One uses constant prior for position parameter k and jeffreys prior on lambdas;
# (3) plot CNV detection results.
#
#

### simulate poisson read counts data with two change points.
sim.x <- c(rpois(30, 30), rpois(40, 50), rpois(30, 25))


### call CNV detection functions.
### (i) constant prior on parameter k and lambdas.
###
source("/Users/Tieming/Dropbox/github/CNV-Detection/mvcm-uniform.R")
mvcm.uniform.res <- mvcm.uniform(x=sim.x, window.size=50, threshold=0.3)
mvcm.uniform.res

### plot
par(mar=c(5,5,1,1), mfrow=c(2,1))
plot(1:length(sim.x), sim.x, pch=19, col="blue", xlab="Position", ylab="Read Count", cex.axis=1.5, cex.lab=1.5)
abline(v=mvcm.uniform.res$est.pnt, col="darkgray", lwd=3, cex.axis=1.5, cex.lab=1.5)
for(mean.index in 1:length(mvcm.uniform.res$est.lambdas)){
  lines(x=c(mvcm.uniform.res$region.start[mean.index], 
            mvcm.uniform.res$region.end[mean.index]),
        y=rep(mvcm.uniform.res$est.lambdas[mean.index], 2), col="darkgray", lty="dashed", lwd=3)
}

barplot(mvcm.uniform.res$all.pp, names.arg=1:length(sim.x), cex.axis=1.5, cex.lab=1.5, cex.names=1.5, xlab="Position", ylab="Posterior Probability")


### (ii) constant prior on parameter k and jeffreys prior on lambdas.
source("/Users/Tieming/Dropbox/github/CNV-Detection/mvcm-jeffreys.R")
mvcm.jeffreys.res <- mvcm.jeffreys(x=sim.x, window.size=50, threshold=0.5)
mvcm.jeffreys.res


### plot
par(mar=c(5,5,1,1), mfrow=c(2,1))
plot(1:length(sim.x), sim.x, pch=19, col="blue", xlab="Position", ylab="Read Count", cex.axis=1.5, cex.lab=1.5)
abline(v=mvcm.jeffreys.res$est.pnt, col="darkgray", lwd=3, cex.axis=1.5, cex.lab=1.5)
for(mean.index in 1:length(mvcm.jeffreys.res$est.lambdas)){
  lines(x=c(mvcm.jeffreys.res$region.start[mean.index], 
            mvcm.jeffreys.res$region.end[mean.index]),
        y=rep(mvcm.jeffreys.res$est.lambdas[mean.index], 2), col="darkgray", lty="dashed", lwd=3)
}

barplot(mvcm.jeffreys.res$all.pp, names.arg=1:length(sim.x), cex.axis=1.5, cex.lab=1.5, cex.names=1.5, xlab="Position", ylab="Posterior Probability")







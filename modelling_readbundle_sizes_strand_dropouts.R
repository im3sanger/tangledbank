# Inigo Martincorena - 03/11/2022
# Playing with RB sizes and strand dropout estimation for nanoseq

# 1. Fit zero-truncated Poisson and zero-truncated NB to observed RB sizes
# 2. Estimate strand dropout

rbfile = "Example_readbundles.RBs"

rb = as.matrix(read.table(rbfile, header=0, sep="\t", stringsAsFactors = F)[,2:3])

# Simulation
if (0) {
  n = 1e6; mu = 1; size = 1; p = 0.3
  m = rnbinom(n=n, mu=mu, size=size)
  m[runif(n)<p] = 0
  m = array(m, dim=c(n/2,2)) 
  rb = m[rowSums(m)>0,] # Removing 0,0 read bundles
}
## 1. Modelling the distribution of single-strand bundle size as zero-truncated Poisson or truncated Negative Binomial

# Zero-truncated Poisson density function
lik_ztp = function(par, x) { 
  -sum(dpois(x=x, lambda=par, log=T)-ppois(q=0, lambda=par, lower.tail = F, log.p = T))
}

# Zero-truncated Negative Binomial density function
lik_ztnb = function(par, x) { 
  -sum(dnbinom(x=x, mu=par[1], size=par[2], log=T) - pnbinom(q=0, mu=par[1], size=par[2], lower.tail = F, log.p = T))
}

# Fitting ZTP and ZTNB to the observed strand bundle distribution
v = rb[rb>0]
mle_pois = optim(par=mean(v), fn=lik_ztp, method="Brent", lower=c(1e-4), upper=c(100), x=v) # Truncated Poisson
mle_nb = optim(par=c(mean(v),2), fn=lik_ztnb, method="L-BFGS-B", lower=c(1e-2,0.1), upper=c(100,100), x=v) # Trunc NB

# Plot 1: fitting truncated Poisson and NB distributions to single-strand RBs
dev.new(width=8, height=8)
par(mfrow=c(2,2))

max_size = 10
fobs = (table(v)/length(v))[1:max_size]; fobs[is.na(fobs)] = 0 # Observed size distribution
fexp_pois = dpois(x=1:10, lambda=mle_pois$par[1]) / ppois(q=0, lambda=mle_nb$par[1], lower.tail = F)
fexp_nb = dnbinom(x=1:10, mu=mle_nb$par[1], size=mle_nb$par[2]) / pnbinom(q=0, mu=mle_nb$par[1], size=mle_nb$par[2], lower.tail = F)

ymax = max(c(fobs,fexp_nb,fexp_pois))*1.2
bc = barplot(fobs, ylim=c(0,ymax), las=1, border=NA, ylab="Frequency", xlab="Single-strand bundle size")
lines(bc, fexp_pois, pch=19, col="cadetblue4"); points(bc, fexp_pois, pch=19, col="cadetblue4", cex=0.5)
lines(bc, fexp_nb, pch=19, col="darkorchid4"); points(bc, fexp_nb, pch=19, col="darkorchid4", cex=0.5)

## 2. Estimating the strand dropout rate (currently minimising the difference in P(0), consider loglik fitting)

pvec = d = seq(0,1,by=0.01) # Values to test
for (j in 1:length(pvec)) {
  p1 = dnbinom(x=0, mu=mle_nb$par[1], size=mle_nb$par[2])
  p2 = p1*(1-pvec[j])+pvec[j]
  d[j] = p2*(1-p2)/(1-p2^2)
}
sdo = pvec[which.min(abs(d-mean(rb==0)))]
text(x=max(bc), y=ymax/1.1, adj = c(1,1), labels=sprintf("NB: m=%0.3g,d=%0.3g\nPois: m=%0.3g\nSDO=%0.3g",mle_nb$par[1],mle_nb$par[2],mle_pois$par[1],sdo))
title(gsub(".RBs","",rbfile))

## 3. Efficiency plots

# Simulating RBs with desired parameters
rbsim = function(mu, size, sdo, n=1e6) {
  m = rnbinom(n=2*n, mu=mu, size=size) # Negative binomial distribution
  m[runif(2*n)<sdo] = 0 # Strand dropout
  m = array(m, dim=c(n,2)) # RB (assuming independence of both strands)
  m = m[rowSums(m)>0,] # Removing (0,0) RBs (never seen)
  return(m)
}

muvec = seq(0.1,10, by=0.1) # mu of the NB
eff = duprate = rvec = xvec = dxvec = array(NA, length(muvec))
for (j in 1:length(muvec)) {
  rbj = rbsim(mu=muvec[j], size=mle_nb$par[2], sdo=sdo)
  rvec[j] = mean(2*rbj)
  duprate[j] = (rvec[j]-1)/rvec[j]
  xvec[j] = sum(rbj)
  dxvec[j] = sum(rowSums(rbj>1)==2)
  eff[j] = dxvec[j] / xvec[j]
}

# Observed metrics
truedx = sum(rowSums(rb>1)==2)
truex = sum(rb)
truer = mean(2*rb)
truedr = (truer-1)/truer
trueeff = truedx/truex
  
# Plot 2: Duplicate rate vs efficiency
plot(duprate, eff, type="l", xlim=c(0,1), xlab="Duplicate rate", ylab="Efficiency (duplex bp / total bp)", las=1)
opt_duprate = duprate[which.max(eff)] # Optimal duplicate rate appears to be 0.805
ind = which(eff > max(eff, na.rm=T)*0.8)
abline(v=opt_duprate, col="red")
abline(v=duprate[c(min(ind),max(ind))], col="red", lty=3)
semiopt_duprate = duprate[c(min(ind),max(ind))]
abline(v=truedr, col="cadetblue"); points(truedr, trueeff, pch=19, col="cadetblue", cex=0.5)
text(x=0, y=max(eff)/1.1, adj = c(0,1), labels=sprintf("obs:%0.3g\nopt:%0.3g\n80%%:%0.3g-%0.3g",truedr,opt_duprate,semiopt_duprate[1],semiopt_duprate[2]))

# Plot 3: Sequencing ratio vs efficiency
plot(rvec, eff, type="l", xlab="Seq ratio (x per molecule)", ylab="Efficiency (duplex bp / total bp)", las=1)
opt_ratio = rvec[which.max(eff)] # Optimal duplicate rate appears to be 5.1
ind = which(eff > max(eff, na.rm=T)*0.8)
abline(v=opt_ratio, col="red")
abline(v=rvec[c(min(ind),max(ind))], col="red", lty=3)
semiopt_ratio = rvec[c(min(ind),max(ind))]
abline(v=truer, col="cadetblue"); points(truer, trueeff, pch=19, col="cadetblue", cex=0.5)
text(x=max(rvec), y=max(eff)/1.1, adj = c(1,1), labels=sprintf("obs:%0.3g\nopt:%0.3g\n80%%:%0.3g-%0.3g",truer,opt_ratio,semiopt_ratio[1],semiopt_ratio[2]))

# Plot 4: Increase in dx per increase in x
truedx_sim = dxvec[which.min(abs(truer-rvec))]
plot(rvec/truer, dxvec/truedx_sim, type="l", las=1, ylab="dx relative to current", xlab="x relative to current")
abline(a=0, b=1, col="grey50")
abline(h=1, col="grey50")
abline(v=1, col="grey50")

dev.copy(pdf, file=sprintf("RBplot_%s.pdf",gsub(".RBs","",rbfile)), width=8, height=8, useDingbats=F); dev.off(); dev.off()

## OTHER USEFUL CODE
if (0) { 

  # fa8
  rbs = data.frame(x=rb[,1], y=rb[,2], size=rowSums(rb))
  rbs = rbs[rbs$size>=4, ]
  total_missed = 0
  for (size in c(4:300)) {
    exp_orphan = (0.5 ** size) * 2
    total_this_size = nrow(rbs[which(rbs$size == size),])
    if (total_this_size > 0) {
      with_both_strands = nrow(rbs[which(rbs$size == size & rbs$x > 0 & rbs$y > 0),])
      obs_orphan = 1 - with_both_strands / total_this_size
      missed = (obs_orphan - exp_orphan) * total_this_size
      total_missed = total_missed + missed
    }
  }
  total_missed_fraction = total_missed / nrow(rbs[which(rbs$size >= 4),])
  total_missed_fraction
  
} # end of if (0)

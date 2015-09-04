#
#
# MVCM algorithm to find CNV change points using
# a constant prior for lambda parameters in the 
# Poisson distributions.
#
# Tieming, 05/19/2015
#
#

mvcm.uniform <- function(x, window.size, threshold){
# Function parameters:
# x: normalized read counts for each bin, recommended bin size: 10kbp or 100kbp.
# window.size: number of bins per window, recommended window size: 40 to 100.
# threshold: posterior probability threshold for claiming a change point,
#            recommended threshold: 0.5 to 0.7 depending on bin and window sizes.
#
#
  chrom.size <- length(x)
  search.start <- 1
  search.end   <- window.size
  est.change.pnt <- NULL
  est.change.pnt.prob <- NULL
  all.pp <- 0
  transform.x <- sqrt(x+3/8)
  
  while(search.start + 1 < chrom.size){
    	tmp.x <- transform.x[search.start:search.end]
  	tmp.window.size <- search.end - search.start + 1
  	k <- 1:(tmp.window.size-1)
  	
  	# y1 is the vector of sum((x_i-sqrt(1/8))^2).
  	y1 <- cumsum((tmp.x-sqrt(1/8))^2)
  	# y2 is the vector of sum((x_[-i]-sqrt(1/8))^2).
  	y2 <- rev(cumsum((rev(tmp.x)-sqrt(1/8))^2))
  	x1.bar <- cumsum(tmp.x)/(1:tmp.window.size)
  	x2.bar <- rev(cumsum(rev(tmp.x))/(1:tmp.window.size))
  	x1.2   <- cumsum(tmp.x^2)
  	x2.2   <- rev(cumsum(rev(tmp.x)^2))
  	
  	delta <- mean(diff(tmp.x))^2
  	prod1 <- (1/2*k)*exp(-2*y1[k]+delta)-
  	         x1.bar[k]*sqrt(pi)/sqrt(2*k)*
  	         exp(-2*x1.2[k]+2*k*(x1.bar[k]^2)+delta)*
  	         (1-pnorm(q=sqrt(1/8), mean=x1.bar[k], sd=sqrt(1/4*k)))
  	prod2 <- (1/2*(tmp.window.size-k))*exp(-2*y2[k+1]+delta)-
  	         x2.bar[k+1]*sqrt(pi)/sqrt(2*(tmp.window.size-k))*
  	         exp(-2*x2.2[k+1]+2*(tmp.window.size-k)*(x2.bar[k+1]^2)+delta)*
  	         (1-pnorm(q=sqrt(1/8), mean=x2.bar[k+1], 
  	         sd=sqrt(1/4*(tmp.window.size-k))))
  	tmp.pi.k <- prod1*prod2
  	if(sum(tmp.pi.k)==0){
      pi.k <- rep(0, tmp.window.size-2)
  	}else{
  	  pi.k <- tmp.pi.k/sum(tmp.pi.k)
  	}
  	all.pp <- c(all.pp, pi.k)

  	if(max(pi.k) >= threshold){
  		est.change.pnt <- c(est.change.pnt, 
  		                    search.start + which.max(pi.k))
  		est.change.pnt.prob <- c(est.change.pnt.prob, max(pi.k))
  	}
  	search.start <- search.end
  	search.end   <- search.start + window.size -1
  	if(search.end > chrom.size){
  		search.end = chrom.size
  	}
  }
  if(search.start + 1 == chrom.size){all.pp <- c(all.pp,0)}
  
  ## compute estimated mean read count in each region given change points.
  num.region   <- length(est.change.pnt) + 1
  lambdas      <- NULL
  region.start <- c(1, est.change.pnt)
  region.end   <- c(est.change.pnt-1, length(x))
  for(region.index in 1:num.region){
  	lambdas <- c(lambdas, median(x[region.start[region.index]:region.end[region.index]]))
  }

  # est.pnt: estimated change point positions;
  # est.pnt.prob: estimated posterior probability 
  #               of a change point at the corresponding position;           
  # If no change point in any bin has posterior probability higher
  # than the threshold, then output will be a NULL vector.  
  # all.pp: posterior probabilities of all positions.       
  return(list(est.pnt=est.change.pnt,
              est.pnt.prob=est.change.pnt.prob,
              all.pp=all.pp,
              est.lambdas=lambdas,
              region.start=region.start,
              region.end=region.end))
}










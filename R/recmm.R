#functions actually used for fitting
#em.run - plain em
#em.cens - using the censored likelihood. 
#em.rreml - with robust = FALSE and TRUE
#em.runbi - for uncensored but not nec normal
#Sbi.ureda 
#LaxA
source("wt_med.r") #Henrik Bengtsson code

em.run <- function(y, i0.mat, w0, m0, sd0, niter=30) {
 	em2 <- em.run1(y, i0.mat, w0, m0, sd0)
 	for (j in 1:niter) {
   		tmp2 <- em.run1(y, em2$id, em2$wt, em2$mean, em2$sd)
   		em2 <- tmp2
 	}
 	list(id=em2$id, wt=em2$wt, mean=em2$mean, sd=em2$sd)
}

em.run1 <- function(y, i0.mat, w0, m0, sd0) {
# 1 run of EM with current proportions, means, SDs
 	k0 <- ncol(i0.mat)
 	nn <- nrow(i0.mat)
 	i1.mat <- i0.mat
 	i1.den <- rep(0, nn)
 	m1 <- m0
    sd1 <- sd0
 	for (j in 1:k0) {
    	i1.mat[, j] <- w0[j]*dnorm(y, m0[j], sd0[j])  #P{obs come from pop j}
    	i1.den <- i1.den + i1.mat[, j]
 	}
 	for (j in 1:k0) {
   		i1.mat[, j] <- i1.mat[, j]/i1.den
		i1.mat[, j][is.nan(i1.mat[, j])] <- 0
    	m1[j] <- sum(i1.mat[, j]*y)/sum(i1.mat[, j])
   		sd1[j] <- sqrt(sum(i1.mat[, j]*(y-m1[j])^2)/sum(i1.mat[, j]))
 	}
 	w1 <- apply(i1.mat, 2, sum)/nn
 	list(id=i1.mat, wt=w1, mean=m1, sd=sd1)
}

em.cens <- function(y, i0.mat, w0, m0, sd0, niter=30, verbose=FALSE) {
 	em2 <- em.run1(y, i0.mat, w0, m0, sd0)
	if(verbose) {
		print('starting means and sds')
		print(m0)
		print(sd0)
	}
 	for (j in 1:niter) {
   		tmp2 <- em.run1c(y, em2$id, em2$wt, em2$mean, em2$sd, verbose=verbose)
   		em2 <- tmp2
 	}
    list(id=em2$id, wt=em2$wt, mean=em2$mean, sd=em2$sd)
}

em.run1c <- function(y, i0.mat, w0, m0, sd0, verbose=verbose) {
# 1 run of EM with current proportions, means, SDs
 	k0 <- ncol(i0.mat)
 	nn <- nrow(i0.mat)
 	i1.mat <- i0.mat
 	i1.den <- rep(0, nn)
 	m1 <- m0
    sd1 <- sd0
 	for (j in 1:k0) {
    	i1.mat[, j] <- w0[j]*dnorm(y, m0[j], sd0[j])  #P{obs come from pop j}
    	i1.den <- i1.den + i1.mat[, j]
 	}
 	for (j in 1:k0) {
   		i1.mat[, j] <- i1.mat[, j]/i1.den
		i1.mat[, j][is.nan(i1.mat[, j])] <- 0

		#set up the maximization
		wtleft <- i1.mat[1, j]
		wtright <- i1.mat[length(y), j]
		xl <- min(y)
		xu <- max(y)
		iwtc <- i1.mat[y>xl&y<xu, j]
		l <- sum(y==xl)
		u <- sum(y==xu)
		unc <- y[y>xl&y<xu]

		mstep <- nlminb(c(m1[j], sd1[j]), negloglikemc, wtleft=wtleft, wtright=wtright, l=l, u=u, xl=xl, xu=xu, iwtc=iwtc, unc=unc)

    	m1[j] <- mstep$par[1]
   		sd1[j] <- mstep$par[2]
		if (verbose & mstep$convergence ==0) {
			print(mstep)
			mu <- m1[j]
			sig  <- sd1[j]
			print(wtleft)
			print(wtright)
			print(l)
			print(u)
			print(xl)
			print(xu)
			print(wtleft*l*log(pnorm(xl, mean=mu, sd=sig)))
			print(wtright*u*log((1-pnorm(xu, mean=mu, sd=sig))))
			midstep <- dnorm(unc, mean=mu, sd=sig)
			print(sum(iwtc[midstep!=0]*log(midstep[midstep!=0])))
			print(mstep$objective)
		}

 	}
 	w1 <- apply(i1.mat, 2, sum)/nn
	if (verbose){
			print("w1/w.true   "); print(rbind( w1, p.true))
			print("m1/m.true   "); print(rbind( m1, m.true))
			print("sd1/s.true  "); print(rbind(sd1, s.true))
	}
 	list(id=i1.mat, wt=w1, mean=m1, sd=sd1)
}

negloglikemc <- function(theta, wtleft, wtright, l, u, xl, xu, iwtc, unc) {
	mu <- theta[1]
	sig <- theta[2]
	
	left <- ifelse(pnorm(xl, mean=mu, sd=sig) >0, wtleft*l*log(pnorm(xl, mean=mu, sd=sig)), 0)
	right <- ifelse(pnorm(xu, mean=mu, sd=sig) ==1, 0, wtright*u*log((1-pnorm(xu, mean=mu, sd=sig))))
	midstep <- dnorm(unc, mean=mu, sd=sig)
	mid <- sum(iwtc[midstep!=0]*log(midstep[midstep!=0]))

	ll <- left + right + mid
	nll <- -1*ll
	nll
}

em.rreml <- function(y, i0.mat, w0, m0, sd0, niter=30, verbose=FALSE, robust=FALSE) {
 	em2 <- em.run1(y, i0.mat, w0, m0, sd0)
	if(verbose) {
		print('starting means and sds')
		print(m0)
		print(sd0)
	}
 	for (j in 1:niter) {
   		tmp2 <- em.run1r(y, em2$id, em2$wt, em2$mean, em2$sd, verbose=verbose, robust=robust)
   		em2 <- tmp2
 	}
	if(verbose) {
		print('ending means and sds')
		print(em2$mean)
		print(em2$sd)
	}
    
    list(id=em2$id, wt=em2$wt, mean=em2$mean, sd=em2$sd)
}

em.run1r <- function(y, i0.mat, w0, m0, sd0, robust=FALSE, verbose=FALSE) {
# 1 run of EM with current proportions, means, SDs
 	k0 <- ncol(i0.mat)
 	nn <- nrow(i0.mat)
 	i1.mat <- i0.mat
 	i1.den <- rep(0, nn)
 	m1 <- m0
    sd1 <- sd0
 	for (j in 1:k0) {
    	i1.mat[, j] <- w0[j]*dnorm(y, m0[j], sd0[j])  #P{obs come from pop j}
    	i1.den <- i1.den + i1.mat[, j]
		if(verbose & sum(i1.mat[, j])==0) {
			print("all wts went to 0: m, sd, sum(w)")
			print(m0[j])
			print(sd0[j])
			print(sum(w0[j]))
		}
 	}
 	for (j in 1:k0) {
        if (verbose & sum(is.nan(i1.mat[, j]))>100) {
            print(sum(is.nan(i1.mat[, j])))
            print("^ number of NaN about to be turned to 0")    
        }
   		i1.mat[, j] <- i1.mat[, j]/i1.den
		i1.mat[, j][is.nan(i1.mat[, j])] <- 0
		if (robust) {
            yl <- min(y)
            yr <- max(y)
            uncens <- y[y>yl&y<yr]
            wts <- i1.mat[y>yl&y<yr, j]
			tmp <- mybiwt.wt2(uncens, 6, wts, sbi.real=TRUE, verbose=verbose)
			tmp <- reml.mix(y, i1.mat[, j], k.m=tmp$t, k.s=tmp$s, verbose=verbose)
			m1[j] <- tmp$m
			sd1[j] <- tmp$s
		} else {
    		m1[j] <- sum(i1.mat[, j]*y)/sum(i1.mat[, j])
   			sd1[j] <- sqrt(sum(i1.mat[, j]*(y-m1[j])^2)/sum(i1.mat[, j]))
			tmp <- reml.mix(y, i1.mat[, j], k.m=m1[j], k.s=sd1[j], verbose=verbose)
			m1[j] <- tmp$m
			sd1[j] <- tmp$s
		}
 	}
 	w1 <- apply(i1.mat, 2, sum)/nn
 	list(id=i1.mat, wt=w1, mean=m1, sd=sd1)
}

mybiwt.wt2 <- function(y, tune.const = 4, owt, verbose=FALSE, sbi.real=FALSE) {
	x <- y[!is.na(y)]
	n <- length(x)
	c0 <- 3*tune.const/2

	if(missing(owt))
	{
		iwt <- rep(1, n)
	} else 
	{
		iwt <- owt
	}
	# the following calculation is the MLE estimates as stated in HTF p 238
	t0 <- sum(iwt*x)/sum(iwt)
	if(verbose) {
		print('first few x and iwt entries, t0, sum(iwt)')
		print(x[1:20])
		print(iwt[1:20])
        print(t0)
        print(sum(iwt))
	}
    ad <- abs(x-t0)
    s0 <- weighted.median(ad, iwt)
    #if it STILL is 0, we have a ton of censoring going on, so try on just the uncensored obs
    if(s0==0) {
        xr <- max(x)
        xl <- min(x)
        ad <- abs(x[x!=xr&x!=xl] -t0)
        s0 <- weighted.median(ad, iwt[x!=xr&x!=xl])
    }

	u <- (x-t0)/(tune.const*s0)
	u2 <- u*u
	wt1 <- ifelse(abs(u) <= 1, 1, 0)
	w <- wt1*(1-u2)*(1-u2)

	t1 <- sum(x*w)/sum(w)
	sbi <- sqrt(sum(w*(x-t0)^2)/sum(w))
	eps <- 1e-4
	dd <- abs((t0-t1)/t1)
    if (is.nan(dd) | is.na(dd)) {
        print('dd was NaN or NA. it was divided by t1, which was divided by sum(w)')
        print(dd)
        print(t1)
        print(sum(w))
    }
	if(verbose) {
		print('change from MLE in one biweight iter, dd:')
		print(dd)
		print(t0)
		print(t1)
	}
	while (dd > eps) {
   		t0 <- t1
   		u <- (x-t0)/(tune.const*sbi)
   		w <- ifelse(abs(u) < 1, (1-u*u)^2, 0) * iwt
   		t1 <- sum(x*w)/sum(w)
   		dd <- abs((t0-t1)/t1)
	}
	if(sbi.real) {
        ad <- abs(x-t1)
        s0 <- weighted.median(ad, iwt)

        u <- (x-t1)/(tune.const*s0)
        w <- ifelse(abs(u) < 1, (1-u*u)^2, 0) * iwt
        p <- w*u
        sipp <- sum(iwt*psiprime(u))
        q <- sum(iwt*p^2)/(sipp*(-1+sipp))
        sqrtq <- sqrt(q)
        sqrtiwt <- sqrt(sum(iwt))
        sbi <- sqrtiwt*tune.const*s0*sqrtq
	}
    if(verbose) {
        print('sbi')
        print(sbi)
    }

	list(t=t1, s=sbi, w=w)
}

psi <- function(u)
{
    y <- ifelse(abs(u) < 1, u*(1-u*u)^2, 0)
}

psiprime <- function(u)
{
	# following Lax... 
   	y <- ifelse(abs(u) < 1, (1-u^2)*(1-5*u^2), 0) 
}

reml.mix <- function(x, iwt, censored, bc=FALSE, k.m=0, k.s=0, verbose=FALSE) {
    n <- sum(iwt)
	xl <- min(x[iwt>0])
	xr <- max(x[iwt>0])
	nl <- sum(iwt[x==xl])
	nr <- sum(iwt[x==xr])

	n0 <- n - nl - nr

	#when restricting the likelihood, we make the assumption that F(u_r) = 1-nr/n & F(u_l) = nl/n
	ul <- qnorm(min(nl/n, 1))
	ur <- qnorm(1 - min(nr/n, 1))
	Zr <- dnorm(ur)/(pnorm(ur)-pnorm(ul))
	Zl <- dnorm(ul)/(pnorm(ur)-pnorm(ul))
	#premultiply to save from Inf u in uncensored case
	uZr <- ifelse(Zr > 0, ur*Zr, 0) 
	uZl <- ifelse(Zl > 0, ul*Zl, 0) 

	y <- (x-xl)*iwt
	d <- xr-xl
	ybar <- sum(y)/n0
	ybarsq <- sum(y^2)/n0
    if (verbose) {
            print (ul)
            print (n1)
            print (n)
    }
    if (Zl == 0) {
        C <- dnorm(ur)*(xr-xl)*n/n0
    } else {
        C <- ybar*ul+dnorm(ur)*(xr-xl)*n/n0
    }
	sstar <- (1/2)*C + (1/2)*sqrt(C^2 + 4*ybarsq)
	mrml <- k.m - sstar*(Zl-Zr)
	srml <- sqrt(k.s^2 -sstar^2*(uZl-uZr-(Zl-Zr)^2))

	if(bc) {
		#see Schneider 86 page 108
		n1 <- 2*n0 - n
		pm <- n0/(n1+1)
		ps <- n0/(n+1)
		Bm <- (-1)*exp(2.692 - 5.439*pm)
		Bs <- (-1)*(0.312 + 0.859*ps)^(-2)
		mc <- mrml - srml*Bm*sign(nr-nl)/(n+1)
		sc <- srml - srml*Bs/(n+1)
		mrml <- mc
		srml <- sc
	}
	
	if (verbose & (is.nan(srml) | is.na(srml))) {
		print('srml is NaN or NA')
		print(k.s^2 -sstar^2*(ul*Zl-ur*Zr-(Zl-Zr)^2))
		print(k.s)
		print(sstar)
        print(nl)
		print(ul)
		print(Zl)
		print(ur)
		print(Zr)
        print(nr)
	}

	list(m=mrml, s=srml)
}

em.runbi <- function(y, i0.mat, w0, m0, sd0, niter=30, verbose=FALSE) {
 	em2 <- em.run1(y, i0.mat, w0, m0, sd0)
	if(verbose) {
		print('starting means and sds')
		print(m0)
		print(sd0)
	}
 	for (j in 1:niter) {
   		tmp2 <- em.run1b(y, em2$id, em2$wt, em2$mean, em2$sd, verbose=verbose)
   		em2 <- tmp2
 	}
    list(id=em2$id, wt=em2$wt, mean=em2$mean, sd=em2$sd)
}

em.run1b <- function(y, i0.mat, w0, m0, sd0, verbose=FALSE) {
# 1 run of EM with current proportions, means, SDs
 	k0 <- ncol(i0.mat)
 	nn <- nrow(i0.mat)
 	i1.mat <- i0.mat
 	i1.den <- rep(0, nn)
 	m1 <- m0
    sd1 <- sd0
 	for (j in 1:k0) {
    	i1.mat[, j] <- w0[j]*dnorm(y, m0[j], sd0[j])  #P{obs come from pop j}
    	i1.den <- i1.den + i1.mat[, j]
 	}
 	for (j in 1:k0) {
   		i1.mat[, j] <- i1.mat[, j]/i1.den
		i1.mat[, j][is.nan(i1.mat[, j])] <- 0
    	tmp   <- mybiwt.wt2(y, 4, i1.mat[, j], verbose=verbose)
    	m1[j] <- tmp$t
   		sd1[j] <- tmp$s
 	}
 	w1 <- apply(i1.mat, 2, sum)/nn
 	list(id=i1.mat, wt=w1, mean=m1, sd=sd1)
}

Sbi.ureda <- function(x, bi, c0, wt, s0=0, verbose=FALSE)
{
    iwt <- wt[!is.nan(wt)]
	if(s0==0)
	{
		ad <- abs(x-bi)
		s0 <- weighted.median(ad, wt)
	}
    #if it STILL is 0, we have a ton of censoring going on, so try on just the uncensored obs
    if(s0==0) {
        xr <- max(x)
        xl <- min(x)
        ad <- abs(x[x!=xr&x!=xl] -bi)
        s0 <- weighted.median(ad, wt[x!=xr&x!=xl])
    }

    u=(x-bi)/(c0*s0)
	
	# UREDA p417 claims KKs thesis says to use sqrt((denom)*(-1+denom)) of normal A est
	# I think thats remove the max from the technometrics formulation
	q <- sum(iwt*psi(u)^2)/(sum(iwt*psiprime(u))*(-1+sum(iwt*psiprime(u))))
	if(verbose) {
		print('s0')
		print(s0)
		print('q numerator')
		print(sum(iwt*psi(u)^2))
		print('q denom part, denom=part*(-1+part)')
		print(sum(iwt*psiprime(u)))
		print('q')
		print(q)
		print('sqrt(sum(iwt))')
		print(sqrt(sum(iwt)))
		print(sqrt(q))
	}

	sqrtq <- sqrt(q)
	sqrtiwt <- sqrt(sum(iwt))
    sbi <- sqrtiwt*c0*s0*sqrtq
    return(sbi)
}

LaxA <- function(x, bi, c0, wt, s0=0) {
	if(s0==0)
	{
		ad <- abs(x-bi)
		s0 <- weighted.median(ad, wt)
	}
    u <- (x-bi)/(c0*s0)

	numerator <- sqrt(sum(wt))*sqrt(sum(wt*(x-bi)^2*(1-u^2)^4))
	denominator <- abs(sum(wt*(1-u^2)*(1-5*u^2)))
	A <- numerator/denominator
	return(A)
}


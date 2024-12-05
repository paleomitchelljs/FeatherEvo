# This code is from Liam Revell's phytool's package.
# A modified version of the phytools::spreadlabels() function from, eg, phenogram()
ff <- function(zz, yy, cost, mm, mo = 1, ms = 1) {           
    ZZ <- cbind(zz - mm/2, zz + mm/2)
    ZZ <- ZZ[order(zz),]
    oo <- 0
    for (i in 2:nrow(ZZ)) oo <- if (ZZ[i - 1, 2] > ZZ[i,1]) {
       oo <- oo + ZZ[i - 1, 2] - ZZ[i, 1]
    }
    else { 
        oo <- oo
    }
    pp <- sum((zz - yy)^2)       
    oo <- if (oo < (1e-06 * diff(par()$usr[3:4]))) {
        0
    }
    else    {
        oo 
    }         
    pp <- if (pp < (1e-06 * diff(par()$usr[3:4])))  {
        0
    }
    else    {
        pp 
    }         
    oo/mo * cost[1] + pp/ms * cost[2]                
}    

# x = values vector
spreadlabels <- function(x, fsize = 0.5, cost=c(1, 0)) {
        rng <- range(x)
        yy <- x
        zz <- setNames((rank(yy, ties.method = "random") - 1)/(length(yy) - 1) * diff(range(yy)) + range(yy)[1], names(yy))      
        mm <- max(fsize * strheight(names(x)))
        mo <- ff(yy, zz, cost = c(1, 0), mm=mm) 
        ms <- ff(yy, zz, cost = c(0, 1), mm=mm) 
        if (mo == 0 && ms == 0){        
            return(yy)    
        }         
        else {
            mo <- ifelse(mo == 0, mo+0.001, mo)
            ms <- ifelse(ms == 0, ms+0.001, ms)
            rr <- optim(zz, ff, yy = yy, mo = mo, ms = ms, mm=mm, cost = cost, method = "L-BFGS-B", lower = rep(rng[1], length(yy)), upper = rep(rng[2], length(yy)))
            return(rr$par)
        }
}

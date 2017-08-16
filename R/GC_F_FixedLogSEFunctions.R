## author: Leonhard Wachutka

## TODO: Leo, the F invokes a warning with biocCheck as it thinks you mean FALSE. maybe you can rename it?
.guess_initial2 <- function(counts, ti, gc, sF)
{
    sF <- sF/sF[1]
    sF <- sF[-1]
    gen.num <- nrow(counts)
    numFinite <- sum(is.finite(ti)) #TODO this is not used
    params.initial <- list(gl=runif(gen.num, 0.01, 0.1),
                            gs=runif(gen.num, 0.01, 0.1),
                            gm=runif(gen.num, 10, 500),
                            sF=sF,
                            gc=gc)
    return (params.initial)
}


.parametrize2 <- function(p.initial)
{
    p <- p.initial
    p$gl <- log(p$gl)
    p$gs <- log(p$gs)
    p$gm <- log(p$gm)
    p$sF <- log(p$sF)
    p$gc <- log(1/p$gc-1)
    return(p)
}


.merge_and_repar2 <- function(p.free, p.fixed, p.free.skeleton)
{
    p <- relist(p.free, skeleton=p.free.skeleton)
    p <- append(p, p.fixed)
    p$gl <- exp(p$gl)
    p$gs <- exp(p$gs)
    p$gm <- exp(p$gm)
    p$sF <- exp(p$sF)
    p$sF <- append(p$sF, p$F0, 0)#the first F is set to F0
    p$gc <- 1/(1+exp(p$gc)) #Fermifunction
    return(p)
}


.SElikelihood2 <- function(param.free, param.fixed, param.free.skeleton, k)
{
    #Resynthesize and reparametrisize our parameters to a single list
    p <- .merge_and_repar2(param.free, param.fixed, param.free.skeleton)
    out <- 0 #the
    #Exon-Intron ist SE-Reads * effective-length
    for(i in 1:length(p$t))
    {
        if(p$t[i] != +Inf)#Steady state has diffrent count formulas due to missing purification
        {
            #Exon-Intron junctions
            EI <- p$sF[i]*p$gm/p$gs*(1-exp(-p$t[i]*p$gs) + p$gc) * p$l * p$uc
        }else{
            #steady state without Purification
            EI <- p$sF[i] * p$gm / p$gs * p$l
        }
        
        out <- out + sum(dnbinom(k[, i], mu=EI, size=p$ga, log=TRUE))
    }
    return (log(-out))
}


.SElikelihood2.grad <- function(param.free, param.fixed, param.free.skeleton, k)
{
    p <- .merge_and_repar2(param.free, param.fixed, param.free.skeleton)
    d_gm <- rep(0, p$N)
    d_gl <- rep(0, p$N)
    d_gs <- rep(0, p$N)
    d_F <- rep(0, length(p$t)-1)
    d_gc <- 0
    #Per Gen values
    for(i in 1:length(p$t))
    {
        
        if(p$t[i] != +Inf)#Steady state has diffrent count formulas due to missing purification
        {
            #Exon-Intron junctions
            EI <- p$sF[i] * p$gm/p$gs * (1-exp(-p$t[i]*p$gs) + p$gc) * p$l * p$uc
            
            d_gm <- d_gm + (k[, i]/EI-1) * p$ga/(EI+p$ga) * EI#First EI
            
            d_gs <- d_gs + (k[, i]/EI-1) * p$ga/(EI+p$ga) * (-p$sF[i]*p$gm/p$gs* (1-exp(-p$t[i]*p$gs) + p$gc)+p$sF[i]*p$gm*exp(-p$t[i]*p$gs)*p$t[i])#First EI
            
        }else{
            #steady state without Purification
            EI <- p$sF[i] * p$gm/p$gs * p$l
            
            d_gm <- d_gm + (k[, i]/EI-1) * p$ga/(EI+p$ga) * EI#First EI
            
            d_gs <- d_gs + (k[, i]/EI-1) * p$ga/(EI+p$ga) * (-EI)#First EI
            
        }
    }
    
    
    out <- 0 #the
    #Exon-Intron ist SE-Reads * effective-length
    for(i in 1:length(p$t))
    {
        if(p$t[i] != +Inf)#Steady state has diffrent count formulas due to missing purification
        {
            #Exon-Intron junctions
            EI <- p$sF[i] * p$gm/p$gs * (1-exp(-p$t[i]*p$gs) + p$gc) * p$l * p$uc
        }else{
            #steady state without Purification
            EI <- p$sF[i] * p$gm / p$gs * p$l
        }
        
        out <- out + sum(dnbinom(k[, i], mu=EI, size=p$ga, log=TRUE))
    }
    
    #dont forget the "minus" because of minimization
    return(c(-d_gs, -d_gm)/-out)
}
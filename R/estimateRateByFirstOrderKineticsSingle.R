## author: Carina Demel

.estimateRateByFirstOrderKineticsSingle <- function(featureCounts, 
                                                    replicate=NULL, 
                                                    BPPARAM=NULL)
{
    ## sample information
    conditions <- featureCounts@colData$condition
    uniqueConditions <- unique(conditions)
    conditionsLabeling <- featureCounts@colData$LT
    labelingTime <- featureCounts@colData$labelingTime
    replicates <- as.character(featureCounts@colData$replicate)
    sequencingDepths <- featureCounts@colData$sequencing.depth
    crossContamination <- featureCounts@colData$cross.contamination
    
    ## gene information
    rows <- rowRanges(featureCounts)
    genes <- names(rows)
    if(is.null(genes)){
        genes <- 1:length(rows)
    }
    lengths <- width(rows)
    names(lengths) <- genes
    counts <- assay(featureCounts)
    rownames(counts) <- genes
    dispersion_L <- rowRanges(featureCounts)$dispersion_L
    dispersion_T <- rowRanges(featureCounts)$dispersion_T
    dispersion <- cbind(dispersion_L, dispersion_T)
    rownames(dispersion) <- genes
    
    stopifnot(!is.null(conditions) &
                  !is.null(conditionsLabeling) & 
                  !is.null(labelingTime) &
                  !is.null(replicates) &
                  !is.null(sequencingDepths) &
                  !is.null(crossContamination) &
                  !is.null(dispersion) & 
                  !is.null(rows))
    
    if(is.null(replicate)){
        ur <- unique(replicates)
        replicate <- c(as.character(ur), paste(ur, collapse=':'))
    }
    
    N <- 1 ## number of cells, my model does not take into account other numbers
    numberOfIterations <- elementMetadata(featureCounts)$numberOfIterations
    if(is.null(numberOfIterations)){
        numberOfIterations <- 3 ##TODO (leo calls it rep) integrate into metadata?
    }
    
    
    singleConditions <- expand.grid(condition=unique(conditions),
                                    rep=replicate,
                                    gene=genes,
                                    iteration=1:numberOfIterations)
    
    ## give worker or use MultiCoreParam
    if(is.null(BPPARAM)){
        BPPARAM <- BiocParallel::MulticoreParam()
    }
    BiocParallel::register(BPPARAM, default=TRUE)
    
    res <- BiocParallel::bplapply(1:nrow(singleConditions), .estimate,
                                  singleConditions=singleConditions, 
                                  counts=counts, dispersion=dispersion, 
                                  conditions=conditions, 
                                  conditionsLabeling=conditionsLabeling, 
                                  replicates=replicates, 
                                  labelingTime=labelingTime, lengths=lengths,
                                  crossContamination=crossContamination, 
                                  sequencingDepths=sequencingDepths)
    
    res <- data.table::rbindlist(res)
    #Take the median of all fits
    res <- res[,.(synthesis=median(synthesis, na.rm=TRUE), degradation=median(degradation, na.rm=TRUE)), 
                by=.(gene=singleConditions$gene, rep=singleConditions$rep, cond=singleConditions$condition)]
    
    labTimeCond = sapply(res$cond, function(c) { subset(featureCounts@colData, condition == c & LT == "L" & replicate == "1")$labelingTime})
    res$labeledAmount = res$synthesis/res$degradation * (1-exp(-res$degradation * labTimeCond))
    res$unlabeledAmount <- res$synthesis/res$degradation - res$labeledAmount
    res$half.life <- log(2)/res$degradation
    
    # ## empty results object code by LEO in estimate...Kinetics.R
    # ## createRResultCubeRates(featureCounts, replicate)
    rRates <- createRResultCubeRatesExtended(featureCounts, replicate)
    assay(rRates[,rRates$rate == 'synthesis']) <- matrix(res$synthesis,
                                                         nrow=length(genes), 
                                                         ncol=length(unique(conditions))*length(replicate), 
                                                         byrow = TRUE)
    assay(rRates[,rRates$rate == 'degradation']) <- matrix(res$degradation,
                                                           nrow=length(genes),
                                                           ncol=length(unique(conditions))*length(replicate),
                                                           byrow = TRUE)
    assay(rRates[,rRates$rate == 'half.life']) <- matrix(res$half.life,
                                                         nrow=length(genes),
                                                         ncol=length(unique(conditions))*length(replicate),
                                                         byrow = TRUE)
    assay(rRates[,rRates$rate == 'labeled.amount']) <- matrix(res$labeledAmount,
                                                              nrow=length(genes),
                                                              ncol=length(unique(conditions))*length(replicate),
                                                              byrow = TRUE)
    assay(rRates[,rRates$rate == 'unlabeled.amount']) <- matrix(res$unlabeledAmount,
                                                                nrow=length(genes),
                                                                ncol=length(unique(conditions))*length(replicate),
                                                                byrow = TRUE)
    return(rRates)
}





## the underlying function to fit labeled and unlabeled RNA amounts to observed
## read count data
.fitAlphaBeta_log_reparam <- function(
    counts,
    labeledSamples,
    totalSamples,
    L,
    N,
    crossContamination,
    sequencingDepths,
    disp,
    alphaInitial,
    betaInitial,
    maxit = 10000,
    trace = FALSE,
    logloglik = TRUE
){
    ## due to reparametrization adjust initial guess values
    alphaInitial <- log(alphaInitial)
    betaInitial <- log(betaInitial)
    
    ## minimize negative log likelihood == maximize likelihood
    fit <- optim(par=c(alphaInitial, betaInitial), fn=.cost, gr=.grad, 
                 counts=counts, L=L, N=N, crossContamination=crossContamination, 
                 sequencingDepths=sequencingDepths, disp=disp,
                 control=list(maxit=maxit, trace=trace), method="BFGS",
                 logloglik=logloglik, labeledSamples=labeledSamples,
                 totalSamples=totalSamples)
    return(fit)
}


## estimation of synthesis and degradation rates for single condition, replicate, and gene
.estimate <- function(index, singleConditions, counts, dispersion, conditions, 
                      conditionsLabeling, replicates, labelingTime, lengths,
                      crossContamination, sequencingDepths){
    message(index)
    c <- singleConditions$condition[index]
    rep <- unlist(strsplit(as.character(singleConditions$rep[index]), ':'))
    gene <- singleConditions$gene[index]
    labeledSamples <- which(conditions == c & conditionsLabeling == "L" & 
                                replicates %in% rep)
    totalSamples <- which(conditions == c & conditionsLabeling == "T" &
                              replicates %in% rep)
    samples <- c(labeledSamples, totalSamples)
    nl <- length(labeledSamples)
    nt <- length(totalSamples)
    
    if (any(counts[gene, samples] > 0) & all(!is.na(dispersion[gene, ]))) {
        muInitialGene <- runif(1, 0.01, 0.1)
        lambdaInitialGene <- runif(1, 0.01, 0.1)
        
        alphaInitialGene <- muInitialGene/lambdaInitialGene * 
            (1-exp(-lambdaInitialGene * labelingTime[labeledSamples[1]]))
        betaInitialGene <- muInitialGene/lambdaInitialGene - alphaInitialGene
        
        tryCatch({
            fab = .fitAlphaBeta_log_reparam(counts=counts[gene, ], 
                                            labeledSamples=labeledSamples, 
                                            totalSamples=totalSamples, 
                                            L=lengths[gene], N=N, 
                                            crossContamination=crossContamination, 
                                            sequencingDepths=sequencingDepths,
                                            disp=dispersion[gene, ], 
                                            alphaInitial=alphaInitialGene,
                                            betaInitial=betaInitialGene)
            
            fittingResult <- exp(fab$par)
            alpha <- fittingResult[1]
            beta <- fittingResult[2]
            lambda <- -1 / labelingTime[labeledSamples[1]] * log(beta / (alpha + beta))
            mu <- (alpha + beta) * lambda
            hl <- log(2)/lambda
            #TODO maybe don't return extended information
            return(data.table('synthesis'=mu, 'degradation'=lambda,
                              'half.life'=hl, 
                              'labeled.amount'=alpha, 'unlabeled.amount'=beta))
        }, error=function(e) {
            return(data.table('synthesis'=NA, 'degradation'=NA,  'half.life'=NA,
                              'labeled.amount'=NA, 'unlabeled.amount'=NA))
        })
        
    }else{
        return(data.table('synthesis'=NA, 'degradation'=NA, 'half.life'=NA,
                          'labeled.amount'=NA, 'unlabeled.amount'=NA))
    }
}


## underlying function to calculated expected counts based on labeled and unlabeld
## RNA amounts
.getExpectedCounts <- function(L, labeledAmount, unlabeledAmount, N, crossCont, seqDepths)
{
    # if(is.na(labeledAmount) & is.na(unlabeledAmount)){
    if(all(is.na(labeledAmount)) & all(is.na(unlabeledAmount))){
        return(rep(NA, length(crossCont))) #used if multiple samples are calculated
    }else{
        # exp.counts <- t(N) * L * t(seqDepths * (t(labeledAmount) + crossCont * t(unlabeledAmount)))
        # #TODO I removed N because of ERror with "  non-conformable arrays" and I am not using it anywhere else (only when calculating for rep 1 and 2 (2 seqDepth and crossCont values))
        exp.counts <- L * t(seqDepths * (t(labeledAmount) + crossCont * t(unlabeledAmount)))
        # L * t(seqDepths * t(as.vector(labeledAmount) + t(outer(crossCont, unlabeledAmount)[,,1])))
        # L * t(seqDepths * t(as.vector(labeledAmount) + t(matrix(outer(crossCont, unlabeledAmount),nrow=length(unlabeledAmount),byrow = T))))
        
        return(exp.counts)
    }
}


.vgetExpectedCounts <- Vectorize(.getExpectedCounts, vectorize.args=c("crossCont", "seqDepths"))


.cost <- function(
    theta,
    counts,
    L,
    N,
    crossContamination,
    sequencingDepths,
    disp,
    logloglik,
    labeledSamples,
    totalSamples
){
    l <- theta[1]
    u <- theta[2]
    labeledDispersion <- disp[1]
    totalDispersion <- disp[2]
    
    ## reparametrization
    labeledAmount <- exp(l) ## labeled amount alpha
    unlabeledAmount <- exp(u) ## unlabeled amount beta
    
    expectedCountsLabeled <- .getExpectedCounts(L=L, labeledAmount=labeledAmount,
                                                unlabeledAmount=unlabeledAmount, 
                                                N=N,
                                                crossCont=crossContamination[labeledSamples],
                                                seqDepths=sequencingDepths[labeledSamples])
    expectedCountsTotal <- .getExpectedCounts(L=L, labeledAmount=labeledAmount,
                                              unlabeledAmount=unlabeledAmount,
                                              N=N,
                                              crossCont=crossContamination[totalSamples],
                                              seqDepths=sequencingDepths[totalSamples])
    
    neg_LL <- -sum(stats::dnbinom(x=counts[labeledSamples], size=labeledDispersion,
                                  mu=expectedCountsLabeled, log=TRUE)) -
        sum(stats::dnbinom(x=counts[totalSamples], size=totalDispersion,
                           mu=expectedCountsTotal, log=TRUE)) 
    
    if (logloglik) {
        res <- log(neg_LL)
    } else {
        res <- neg_LL
    }
    return(res)
}

## gradient function
.grad <- function(
    theta,
    counts,
    L,
    N,
    crossContamination,
    sequencingDepths,
    disp,
    logloglik,
    labeledSamples,
    totalSamples
){
    l <- theta[1]
    u <- theta[2]
    samples <- c(labeledSamples, totalSamples)
    
    ## reparametrization
    labeledAmount <- exp(l) ## labeled amount alpha
    unlabeledAmount <- exp(u) ## unlabeled amount beta
    
    expectedCounts <- N * L * sequencingDepths[samples] *
        (labeledAmount + crossContamination[samples] * unlabeledAmount)
    
    if (length(samples) > 2) { ## replicates
        disp <- c(rep(disp[1], length(labeledSamples)), 
                  rep(disp[2], length(totalSamples)))
    }
    ## vectorial solution
    d_a <- - sum((counts[samples] - expectedCounts) * N * L * sequencingDepths[samples] * 
                     labeledAmount/(expectedCounts * (1 + expectedCounts / disp)))
    d_b <- - sum((counts[samples] - expectedCounts) * N * L * sequencingDepths[samples] *
                     (crossContamination[samples] * unlabeledAmount) / 
                     (expectedCounts * (1 + expectedCounts / disp)) )
    
    if (logloglik) {
        neg_LL <- -sum(stats::dnbinom(x=counts[samples], size=disp, 
                                      mu=expectedCounts, log=TRUE))
        res <- c(d_a, d_b) / neg_LL
    } else {
        res <- c(d_a, d_b)
    }
    return(res)
}

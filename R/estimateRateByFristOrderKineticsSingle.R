## author: Carina Demel

.estimateRateByFirstOrderKineticsSingle <- function(featureCounts, replicate=NULL, BPPARAM=NULL)
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
    minAlphaBeta <- 1e-4
    #TODO What about initial values
    alphaInitial = betaInitial = NULL
    
    singleConditions <- expand.grid(condition=unique(conditions),
                                    rep=replicate,
                                    gene=genes)

    estimate <- function(index){
        c <- singleConditions$condition[index]
        rep <- unlist(strsplit(as.character(singleConditions$rep[index]), ':')) #unlist(singleConditions$rep[index])
        
        gene <- singleConditions$gene[index]
        labeledSamples <- which(conditions == c & conditionsLabeling == "L" & replicates %in% rep)
        totalSamples <- which(conditions == c & conditionsLabeling == "T" & replicates %in% rep)
        samples <- c(labeledSamples, totalSamples)
        nl <- length(labeledSamples)
        nt <- length(totalSamples)
        
        #TODO initialization: random, run 10 times, take average result?
        if (any(counts[gene, samples] > 0) & all(!is.na(dispersion[gene, ]))) {	
            if (is.null(alphaInitial)) {
                alphaInitialGene <- 1 / (N * lengths[gene]) * 
                    1 / (nl * sum(crossContamination[totalSamples]) -
                            nt * sum(crossContamination[labeledSamples])) * 
                    (sum(crossContamination[totalSamples]) * 
                         sum(counts[gene, labeledSamples] / sequencingDepths[labeledSamples])	-
                         sum(crossContamination[labeledSamples]) * 
                         sum(counts[gene, totalSamples] / (sequencingDepths[totalSamples])))
            } else if (length(alphaInitial) == 1) {
                alphaInitialGene <- alphaInitial
            } else {
                alphaInitialGene <- alphaInitial[gene, t]
            }
            
            if (is.null(betaInitial)) {
                betaInitialGene <- 1 / (N * lengths[gene]) * 
                    1 / (sum(crossContamination[totalSamples]) / nt - 
                             sum(crossContamination[labeledSamples]) / nl) * 
                    ( sum(counts[gene, totalSamples] / 
                              (nt * sequencingDepths[totalSamples])) -
                          sum(counts[gene, labeledSamples] /
                                  (nl * sequencingDepths[labeledSamples])))
            } else if (length(betaInitial) == 1) {
                betaInitialGene <- betaInitial
            } else {
                betaInitialGene <- betaInitial[gene, t]
            }
            ## set alpha and beta to a minimum value, so that the estimation on log
            ## scale can be performed. (should be non-negative alpha and beta values)
            alphaInitialGene <- max(alphaInitialGene, minAlphaBeta)
            betaInitialGene <- max(betaInitialGene, minAlphaBeta)				
            
            tryCatch({
                fab = .fitAlphaBeta_log_reparam(counts[gene, ], labeledSamples, 
                                               totalSamples, lengths[gene], N, 
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
                ## exp counts are a vector of the lenght of individual labeled/total samples under consideration
                ## therefore do.call("rbind", res) does not work if looped over reps and combinations of reps
                # expectedCountsLabeled <- .getExpectedCounts(lengths[gene], alpha, beta, N,
                #                                   crossContamination[labeledSamples],
                #                                   sequencingDepths[labeledSamples])
                # expectedCountsTotal <- .getExpectedCounts(lengths[gene], alpha, beta, N,
                #                                   crossContamination[totalSamples],
                #                                   sequencingDepths[totalSamples])

                # loss = fab$value # not return when we initialize 10 times
                return(c('synthesis'=mu, 'degradation'=lambda, 'half.life'=hl, 'labeled.amount'=alpha, 'unlabeled.amount'=beta))#, expectedCountsLabeled, expectedCountsTotal))
            }, error=function(e) {
                return(c('synthesis'=NA, 'degradation'=NA,  'half.life'=NA, 'labeled.amount'=NA, 'unlabeled.amount'=NA))#, rep(NA, nl), rep(NA,nt)))
            })	
            
        }else{
            return(c('synthesis'=NA, 'degradation'=NA, 'half.life'=NA, 'labeled.amount'=NA, 'unlabeled.amount'=NA))#, rep(NA, nl), rep(NA,nt)))
        }
        
    }
    
    ## give worker or use MultiCoreParam
    if(is.null(BPPARAM)){
        BPPARAM <- MulticoreParam(workers = ncores)
    }
    BiocParallel::register(BPPARAM, default = TRUE)
    ## BiocParallel::registered()
    
    res <- BiocParallel::bplapply(1:nrow(singleConditions), estimate)
    res <- do.call("rbind", res)

    # resTable <- cbind(singleConditions, as.data.frame(res))
    # resTableReshaped <- data.table::melt(resTable, id.vars=c("condition","rep","gene"))
    # resTableReshaped$merge_var <- paste(resTableReshaped$condition, resTableReshaped$variable, resTableReshaped$rep, sep="_")
    # resTableReshaped <- reshape(resTableReshaped, idvar=c("gene"), timevar=c("merge_var"),
    #                              direction="wide", drop=c("condition", "rep", "variable"))
    # colnames(resTableReshaped) <- gsub("value.", "", colnames(resTableReshaped))
    # rownames(resTableReshaped) <- resTableReshaped$gene
    #    
    # # match.L <- match(conditions[labeledSamples], uniqueConditions)
    # # match.T <- match(conditions[totalSamples], uniqueConditions)
    #     
    # # expectedCountsLabeled <- getExpectedCounts(lengths[gene.indices], alpha[, match.L], 
    # #                                   beta[, match.L], N, crossContamination[labeledSamples], 
    # #                                   sequencingDepths[labeledSamples]) 	
    # # expectedCountsTotal <- getExpectedCounts(lengths[gene.indices], alpha[, match.T], 
    # #                                   beta[, match.T], N, crossContamination[totalSamples], sequencingDepths[totalSamples])
    #     
    # ## empty results object code by LEO in estimate...Kinetics.R
    # ## createRResultCubeRates(featureCounts, replicate)
    rRates <- createRResultCubeRatesExtended(featureCounts, replicate)
    # assay(rRates) <- resTableReshaped[ ,colnames(assay(rates))]

    assay(rRates[,rRates$rate == 'synthesis']) = matrix(res[,'synthesis'], nrow=length(genes), 
                                                        ncol=length(unique(conditions))*length(replicate), byrow = TRUE)
    assay(rRates[,rRates$rate == 'degradation']) = matrix(res[,'degradation'], nrow=length(genes), 
                                                        ncol=length(unique(conditions))*length(replicate), byrow = TRUE)
    assay(rRates[,rRates$rate == 'half.life']) = matrix(res[,'half.life'], nrow=length(genes), 
                                                        ncol=length(unique(conditions))*length(replicate), byrow = TRUE)
    assay(rRates[,rRates$rate == 'labeled.amount']) = matrix(res[,'labeled.amount'], nrow=length(genes), 
                                                        ncol=length(unique(conditions))*length(replicate), byrow = TRUE)
    assay(rRates[,rRates$rate == 'unlabeled.amount']) = matrix(res[,'unlabeled.amount'], nrow=length(genes), 
                                                        ncol=length(unique(conditions))*length(replicate), byrow = TRUE)
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
    
    cost <- function(
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
        
        expectedCountsLabeled <- .getExpectedCounts(L, labeledAmount, unlabeledAmount, N,
                                          crossContamination[labeledSamples], sequencingDepths[labeledSamples])
        expectedCountsTotal <- .getExpectedCounts(L, labeledAmount, unlabeledAmount, N,
                                          crossContamination[totalSamples], sequencingDepths[totalSamples])
        
        neg_LL <- -sum(dnbinom(x=counts[labeledSamples], size=labeledDispersion,
                               mu=expectedCountsLabeled, log=TRUE)) -
            sum(dnbinom(x=counts[totalSamples], size=totalDispersion,
                        mu=expectedCountsTotal, log=TRUE)) 
        
        if (logloglik) {
            res <- log(neg_LL)
        } else {
            res <- neg_LL
        }
        return(res)
    }
    
    grad <- function(
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
            neg_LL <- -sum(dnbinom(x=counts[samples], size=disp, 
                                   mu=expectedCounts, log=TRUE))
            res <- c(d_a, d_b) / neg_LL
        } else {
            res <- c(d_a, d_b)
        }
        return(res)
    }
    
    ## minimize negative log likelihood == maximize likelihood
    fit <- optim(par=c(alphaInitial, betaInitial), fn=cost, gr=grad, 
                 counts=counts, L=L, N=N, crossContamination=crossContamination, 
                 sequencingDepths=sequencingDepths, disp=disp,
                 control=list(maxit=maxit, trace=trace), method="BFGS",
                 logloglik=logloglik, labeledSamples=labeledSamples,
                 totalSamples=totalSamples)
    return(fit)
}


## underlying function to calculated expected counts based on labeled and unlabeld
## RNA amounts
.getExpectedCounts <- function(L, labeledAmount, unlabeledAmount, N, crossCont, seqDepths)
{
    exp.counts <- t(N) * L * t(seqDepths * (t(labeledAmount) + crossCont * t(unlabeledAmount)))
    return(exp.counts)
}




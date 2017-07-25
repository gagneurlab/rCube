context("Fitting Functions")


test_that("costfunction gives correct result", {
    expect_equal(.cost(theta=c(10,10),
                       counts=c(1000,500),
                       L=1000,
                       N=1,
                       crossContamination=c(0.01, 1),
                       sequencingDepths=c(0.98, 0.52),
                       disp=c(0.2, 0.002),
                       logloglik=FALSE,
                       labeledSamples=1,
                       totalSamples=2),
                 23.2137064364781)
})


test_that("gradient is correct", {
    require(numDeriv)
    expect_equal(.grad(theta=c(10,10),
                       counts=c(1000,500),
                       L=1000,
                       N=1,
                       crossContamination=c(0.01, 1),
                       sequencingDepths=c(0.98, 0.52),
                       disp=c(0.2, 0.002),
                       logloglik=FALSE,
                       labeledSamples=1,
                       totalSamples=2),
                 numDeriv::grad(function(u) .cost(u, counts=c(1000,500),
                                                  L=1000,
                                                  N=1,
                                                  crossContamination=c(0.01, 1),
                                                  sequencingDepths=c(0.98, 0.52),
                                                  disp=c(0.2, 0.002),
                                                  logloglik=FALSE,
                                                  labeledSamples=1,
                                                  totalSamples=2), c(10,10)))
})


test_that("optim function", {
    expect_equal(.fitAlphaBeta_log_reparam(counts=c(9898, 10400),
                                           labeledSamples=1,
                                           totalSamples=2,
                                           L=1000,
                                           N=1,
                                           crossContamination=c(0.01, 1),
                                           sequencingDepths=c(0.98, 0.52),
                                           disp=c(0.2,0.002),
                                           alphaInitial=10,
                                           betaInitial=10,
                                           maxit=10000,
                                           trace=FALSE,
                                           logloglik=TRUE),
                 list("par"=c(2.30258509299405, 2.30258509299405),"value"=3.28554321041966,"counts"=c("function"=1,"gradient"=1),"convergence"=0,"message"=NULL))
})

context("Expected Counts")

test_that(".getExpectedCounts gives correct numbers", {
    ### Single-gene, single condition test
    expect_equal(as.vector(.getExpectedCounts(L=1000, labeledAmount=10, 
                                              unlabeledAmount=10,
                                              N=1,
                                              crossCont=0.01,
                                              seqDepths=0.98)),
                 9898)
    ### single gene, multi condition
    # expect_equal(as.vector(.getExpectedCounts(L=1000,
    #                                           labeledAmount=10,
    #                                           unlabeledAmount=10,
    #                                           N=1,
    #                                           crossCont=c(0.01, 1),
    #                                           seqDepths=c(0.98,0.52))),
    #              c(9898,10400))
### Multi-gene, single condition test
    expect_equal(as.vector(.getExpectedCounts(L=c(1000, 2000),
                                              labeledAmount=c(10, 10),
                                              unlabeledAmount=c(10, 10),
                                              N=1, crossCont=0.01, 
                                              seqDepths=0.98)),
                 c(9898, 19796))
    
})


test_that(".vgetExpectedCounts gives correct numbers", {
    expect_equal(.vgetExpectedCounts(L=1000,
                                     labeledAmount=10,
                                     unlabeledAmount=10,
                                     N=1,
                                     crossCont=c(0.01, 1),
                                     seqDepths=c(0.98, 0.52)),
                 c(9898, 10400))
    expect_equal(.vgetExpectedCounts(L=c(1000, 2000),
                                     labeledAmount=c(10, 10),
                                     unlabeledAmount=c(10, 10),
                                     N=1,
                                     crossCont=c(0.01, 1),
                                     seqDepths=c(0.98 ,0.52)),
                 matrix(c(9898, 10400, 19796, 20800), nrow=2, byrow=TRUE))
})

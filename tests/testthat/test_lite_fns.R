library('M3D')
context('lite functions')

test_that('M3D_Wrapper_lite works properly', {
    data(rrbsDemo)
    data(CpGsDemo)
    data(MMDlistDemo)
    CpGs1 <- CpGsDemo[1:10]
    CpGs2 <- CpGsDemo[110:120]
    overlaps1 <- findOverlaps(CpGs1,rrbsDemo)
    overlaps2 <- findOverlaps(CpGs2,rrbsDemo)
    M3D_stat_lite1 <- M3D_Wrapper_lite(rrbsDemo, overlaps1)
    M3D_stat_lite2 <- M3D_Wrapper_lite(rrbsDemo, overlaps2)
    stat <- MMDlistDemo[[1]]-MMDlistDemo[[2]]
    res <- cbind(rowMeans(stat[,2:5],na.rm=TRUE),stat[,1],stat[,6])
    colnames(res) <- c('MeanBetween','H1-hESC1  vs  H1-hESC2','K562-1  vs  K562-2')
    expect_equal(res[1:10,],M3D_stat_lite1)
    expect_equal(res[110:120,],M3D_stat_lite2)
})

test_that('pvals_lite works properly', {
    data('rrbsDemo')
    data('CpGsDemo')
    data('MMDlistDemo')
    data('PDemo')
    group1 <- 'H1-hESC'
    group2 <- 'K562'
    M3D <- MMDlistDemo$Full-MMDlistDemo$Coverage
    M3D <- cbind(rowMeans(M3D[,2:5],na.rm=TRUE),M3D[,1],M3D[,6])
    colnames(M3D) <- c('MeanBetween','H1-hESC1  vs  H1-hESC2','K562-1  vs  K562-2')
    expect_equal(pvals_lite(rrbsDemo, CpGsDemo, M3D, group1, group2,
                           smaller=FALSE,comparison='allReps',method='empirical',
                           closePara=0.005),PDemo, tolerance=1.0e-4)
})
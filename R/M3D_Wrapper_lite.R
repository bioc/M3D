#' Computes the components of the M3D test-statistic over all regions for all
#'  sample-pairs.
#'
#' Returns the M3D test-statistic, without the compenents, for all regions
#' with the cross group sample pairs averaged to save memory (in column one),
#'  as a matrix.
#'
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param overlaps The overlaps between the list of testing regions and the
#'  methylation data. This is obtained using the
#' function findOverlaps(CpGs,rrbs) for a GRanges object CpGs detailing the
#'  testing regions.
#' @param group1 The name of the first group for the comparison. This is stored
#'  in colData(rrbs). Default finds first unique group in colData(rrbs).
#' @param group2 The name of the second group for the comparison. This is stored
#'  in colData(rrbs). Default finds second unique group in colData(rrbs).
#' @param verbose Logical vector. If true, the function prints a progress bar.
#' @return This returns the two components of the M3D test-statistic for each
#'  region over all sample pairs as a matrix.
#' Subtracting them gives the M3D test-statistic. This is processed with the
#'  function pvals.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @references Gretton, A., Borgwardt, K. M., Rasch, M., Scholkopf, B., Smola,
#'  A. J. (2006). A kernel method for the two-sample-problem. In Advances in
#'   neural information processing systems (pp. 513-520).
#' @examples
#' data(rrbsDemo)
#' data(CpGsDemo)
#' CpGsDemo <- CpGsDemo[1:5]
#' overlaps <- GenomicRanges::findOverlaps(CpGsDemo,rrbsDemo)
#' M3D_list <- M3D_Wrapper_lite(rrbsDemo,overlaps)
#' head(M3D_list)
#' @export

M3D_Wrapper_lite <- function(rrbs, overlaps, group1=NaN, group2=NaN, verbose = TRUE){
    if (is.na(group1)){
        group1 <- as.character(unique(colData(rrbs)$group)[1])
    }
    if (is.na(group2)){
        group2 <- as.character(unique(colData(rrbs)$group)[2])
    }
    nSamples = sum(colData(rrbs)$group==group1,colData(rrbs)$group==group2)
    if (nSamples==2){
        samplesIdx <- c(1,2)
        numPairs <- 1
    } else {
        a <- unlist(lapply(1:(nSamples-1), function(i) rep(i,(nSamples-i))))
        b <- unlist(lapply(1:(nSamples-1), function(i) (i+1):nSamples))
        samplesIdx <- cbind(a,b)
        numPairs <- length(samplesIdx[,1])
    }

    # store objects here to avoid slow lookups, potential memory hazard
    m_reads <- methReads(rrbs)
    t_reads <- totalReads(rrbs)
    locs_all <- start(ranges(rowRanges(rrbs)))
    q_hits <- queryHits(overlaps)
    sub_hits <- subjectHits((overlaps))

    ### make colnames - vectorize
    all_sample_names <- rownames(colData(rrbs))[colData(rrbs)$group==group1]
    all_sample_names <- c(all_sample_names,
                          rownames(colData(rrbs))[colData(rrbs)$group==group2])
    ColumnNames <- unlist(lapply(1:numPairs, function(pairInd){
        if (numPairs==1){
            pair <- c(1,2)
        } else {
            pair <- samplesIdx[pairInd,]
        }
        sample1 <- all_sample_names[pair[1]]
        sample2 <- all_sample_names[pair[2]]
        return(paste(sample1, ' vs ', sample2))
    }))

    l <- length(ColumnNames)
    temp <- matrix(rep(0,2*l),nrow=2,ncol=l) # to work with findcomps (hack)
    colnames(temp) <- ColumnNames
    ### code from pvals
    samples1 <- rownames(colData(rrbs))[colData(rrbs)[,]==group1]
    samples2 <- rownames(colData(rrbs))[colData(rrbs)[,]==group2]
    within1 <- determineGroupComps(samples1,type='within')
    within2 <- determineGroupComps(samples2,type='within')
    within <- c(within1,within2)
    between <- determineGroupComps(samples1,samples2,type='between')

    idsWithin <- c(findComps(temp,within1),findComps(temp,within2))
    idsBetween <- findComps(temp,between)
    ##### end code from pvals

    islands <- unique(queryHits(overlaps))
    M3D_stat <- matrix(NA,nrow=length(islands),ncol=1+length(idsWithin))
    col_sample_names <- c(samples1, samples2)

    # loop over islands, then over samples
    if(verbose){
        pb <- txtProgressBar(min=1,max=length(islands),style=3)
    }
    for (i in 1:length(islands)) {
        island <- islands[i]
        methIndices <- sub_hits[q_hits==island]
        meth_isl <- m_reads[methIndices,]
        total_isl <- t_reads[methIndices,]
        unmeth_isl <- total_isl - meth_isl

        # compute the location matrix
        ## define easy labels
        locs <- locs_all[methIndices]
        locs <- locs - min(locs)
        # make location matrix
        G <- locs%*%t(locs)
        L <- ncol(G)
        nor <- rep(G[seq(1,L^2,L+1)],L)
        # locMx is the squared distance between sites
        locMx <- -2*G + matrix(nor, nrow=L, byrow=TRUE) + matrix(nor, nrow=L)
        locInds <- which(locMx!=0)

        temp <- rep(NA,numPairs)
        for (pairInd in 1:numPairs){
            if (numPairs==1){
                pair <- c(1,2)
            } else {
                pair <- samplesIdx[pairInd,]
            }
            sample1 <- col_sample_names[pair[1]]
            sample2 <- col_sample_names[pair[2]]

            # this is the meth data
            methData <- meth_isl[,c(sample1,sample2)]
            # this is the total data
            unmethData <- unmeth_isl[,c(sample1,sample2)]
            testData <- matrix(0, nrow = dim(methData)[1], ncol = 4)
            testData[,1:2] <- methData
            testData[,3:4] <- unmethData
            # compute the MMDs
            res <- M3D_Single(testData,locMx, locInds, method='MinusCovMMD')
            temp[pairInd] <- res[[1]] - res[[2]] # MMD - MMDCoverage
            # store as vector for the island
        }
        M3D_stat[i,1] <- mean(temp[idsBetween],na.rm=TRUE)
        M3D_stat[i,2:(1+length(idsWithin))] <- temp[idsWithin]
        if(verbose){
            setTxtProgressBar(pb,i)
        }
    }
    if(verbose){
        close(pb)
    }
    colnames(M3D_stat) <- c('MeanBetween',ColumnNames[idsWithin])
    return(M3D_stat)
}

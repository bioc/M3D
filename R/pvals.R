#' Computes p-values
#' 
#' Returns p-values for each region reflecting the probability of observing the
#'  mean test-statistic 
#' of the between group comparisons among the inter-replicate comparisons.
#' 
#' @param rrbs An rrbs object containing methylation and coverage data as
#'  created using the BiSeq pacakge
#' @param CpGs A GRanges object with each row being a testing region
#' @param MMD A matrix containing the M3D test-statistic, the difference the
#'  full and methylation blind metrics,
#' for each region in the CpGs object. Each column is a comparison between two
#'  samples, which are described in the column names.
#' @param group1 The name of the first group for the comparison. This is stored
#'  in colData(rrbs)
#' @param group2 The name of the second group for the comparison. This is stored
#'  in colData(rrbs)
#' @param smaller Determines whether the p-value is computed whether the
#'  test-statistic is greater or lesser than inter-replicate 
#' values. For our purposes, it should be set to FALSE.
#' @param comparison Details which groups we are using to define our empirical
#'  testing distribution. The default is to 
#' use all of them, however, should the user find one group contains unusually
#'  high variability, then that group can be selected. 
#' Values are either 'allReps', 'Group1' or 'Group2'. 
#' @return Returns a list P, with 2 entries. 'FDRmean' is the Benjamini-Hochberg
#'  adjusted p-values. The unadjusted p-values are stored in 'Pmean'.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @export
#' @import BiSeq
#' @import GenomicRanges
#' @import IRanges
#' @examples
#' \donttest{data(rrbsDemo)
#' data(CpGsDemo)
#' data(MMDlistDemo)
#' M3Dstat <- MMDlistDemo$Full-MMDlistDemo$Coverage
#' group1 <- unique(colData(rrbsDemo)$group)[1]
#' group2 <-unique(colData(rrbsDemo)$group)[2]
#' PDemo <- pvals(rrbsDemo, CpGsDemo, M3Dstat,
#'             group1, group2, smaller=FALSE,comparison='allReps')}


pvals <- function(rrbs, CpGs, MMD, group1, group2, 
                  smaller=FALSE,comparison='allReps'){
  # get the names of the columns to compare from MMD object
  samples1 <- rownames(colData(rrbs))[colData(rrbs)[,]==group1]
  samples2 <- rownames(colData(rrbs))[colData(rrbs)[,]==group2]
  within1 <- determineGroupComps(samples1,type='within')
  within2 <- determineGroupComps(samples2,type='within')
  within <- c(within1,within2)
  between <- determineGroupComps(samples1,samples2,type='between')
  
  ids1 <- findComps(MMD,within1)
  ids2 <- findComps(MMD,within2)
  ids3 <- findComps(MMD,between)
  
  # get the MMD in order of between groups, then within groups
  Within <- cbind(MMD[,within1],MMD[,within2])
  colnames(Within) <- c(within1,within2)
  Between <- MMD[,between]
  AvsB <- cbind(Between,Within)
  
  # find the total counts over each island, per sample
  nCpGs <- length(CpGs)
  ovlaps <- findOverlaps(CpGs,rowData(rrbs))
  islandList <- unique(queryHits(ovlaps))
  nIslands <- length(islandList)
  
  Pmean <- lapply(1:nIslands, function(j) {
    if (smaller==TRUE){
      if (comparison=='Group1'){
        length(which(MMD[,within1] <= mean(
          Between[j,],na.rm=TRUE))) / length(MMD[,within1])
      } else if (comparison=='Group2')  {
        length(which(MMD[,within2] <= mean(
          Between[j,],na.rm=TRUE)))/length(MMD[,within2])      
      } else {
        length(which(Within <= mean(
          Between[j,],na.rm=TRUE)))/length(Within)
      }
      
    } else {
      
      if (comparison=='Group1'){
        length(which(MMD[,within1] >= mean(
          Between[j,],na.rm=TRUE)))/length(MMD[,within1])
      } else if (comparison=='Group2')  {
        length(which(MMD[,within2] >= mean(
          Between[j,],na.rm=TRUE)))/length(MMD[,within2])     
      } else {
        length(which(Within >= mean(
          Between[j,],na.rm=TRUE)))/length(Within)
      }
    }
  })
  Pmean <- unlist(Pmean)
  FDRmean <- p.adjust(as.vector(Pmean),method='BH')
  P <- list(Pmean,FDRmean)
  names(P) <- c('Pmean','FDRmean')
  return(P)
}

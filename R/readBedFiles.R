#' Reads in bed file to make RRBS data
#'
#' Reads in RRBS data in bed file format from various styles and outputs
#'  an rrbs data structure. Adapted from readBismark in the BiSeq package.
#'
#' @param files A character pointing the the rrbs files downloads from the
#'  ENCODE database.
#' @param colData Samples' names plus additional sample information as
#'  character, data.frame or DataFrame.
#' @param bed_type Character string representing the style of the bed file.
#' Options are
#' @param eData Experiment data to describe the work. This is used to create
#' the BSraw object as in the BiSeq package.
#' @return Returns a BSraw object storing methylation and coverage data -
#'  the underlying structure for this package.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @export
#' @examples
#' \donttest{
#' # download the files and change the working directory
#' # to that location
#' files <- c('wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz',
#' 'wgEncodeHaibMethylRrbsH1hescHaibSitesRep2.bed.gz',
#' 'wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed.gz',
#' 'wgEncodeHaibMethylRrbsK562HaibSitesRep2.bed.gz')
#' group <- factor(c('H1-hESC','H1-hESC','K562','K562'))
#' samples <- c('H1-hESC1','H1-hESC2','K562-1','K562-2')
#' colData <- DataFrame(group,row.names= samples)
#' rrbs <- readBedFiles(files,colData)}

readBedFiles <- function(files, colData, bed_type = 'Encode', eData=NaN) {
    if (nrow(colData) != length(files)) {
        stop("Row number of colData must equal length of files.")
    }

    if (!test_colData_structure(colData)){
        stop("Files must be specified group by group. Please re-make colData to
             show all of one group, then all of the next, etc.")
    }

    message('Reading in files\n')
    methData <- readBedRaw(files, colData, bed_type='Encode')

    cat("Building BSraw object.\n")

    fData <- methData[[1]]

    if(length(methData) > 1){
        for(i in seq(along=methData)[-1]){
            fData <- unique(c(fData, methData[[i]]))
        }
    }
    elementMetadata(fData) <- NULL
    names(fData) <- as.character(1:length(fData))

    tReads <- matrix(integer(length = length(fData) *
                                 length(methData)), nrow=length(fData))
    mReads <- matrix(integer(length = length(fData) *
                                 length(methData)), nrow=length(fData))

    for(i in seq(along=methData)){
        mtch <- findOverlaps(fData, methData[[i]])
        mtch.m <-  as.matrix(mtch)
        ind1 <- mtch.m[, 1]
        ind2 <- mtch.m[, 2]
        tReads[ind1, i] <- elementMetadata(methData[[i]])$reads[ind2]
        mReads[ind1, i] <- elementMetadata(methData[[i]])$methylated[ind2]
    }

    colnames(tReads) <- rownames(colData)
    colnames(mReads) <- rownames(colData)
    rownames(tReads) <- names(fData)
    rownames(mReads) <- names(fData)
    if(is.na(eData)){
        eData <- SimpleList(exp='Demo')
    }

    rrbs = BSraw(
        exptData=eData,
        colData = colData,
        rowRanges = fData,
        totalReads = tReads,
        methReads = mReads)

    return(rrbs)
}

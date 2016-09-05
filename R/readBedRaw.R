#' Reads in RRBS data in raw format to be processed later
#'
#' Reads in RRBS data in bed file format from a variety of formats and outputs
#' a list that is turned into an rrbs data structure in the wrapper function.
#' Adapted from readBismark in the BiSeq package.
#'
#' @param files A character pointing the the rrbs files downloads from the
#'  ENCODE database.
#' @param colData Samples' names plus additional sample information as
#'  character, data.frame or DataFrame.
#' @param bed_type Character string representing the style of the bed file.
#' Options are
#' @return Returns a list with all the data that will be turned into a BSraw
#' object -  the underlying structure for this package.
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
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
#' rrbs <- readBedRaw(files,colData)}

readBedRaw <- function(files, colData, bed_type='Encode') {

    methData = list()

    for (i in 1:length(files)) {
        cat(paste("Processing sample ", rownames(colData)[i], " ... \n", sep=""))
        if (bed_type == 'Encode'){
            bismark <- scan(files[i], skip=1, sep="\t",
                            what=list("character", integer(), NULL, NULL,
                                      integer(), 'character',
                                      NULL,NULL,NULL,NULL,integer()))
            message('Begin sorting data')
            bismarkEntries <- c(1,2,5,6,11)
            Order <- with(bismark,order(bismark[[6]],bismark[[1]],bismark[[2]]))
            for (j in 1:length(bismarkEntries)){
                entry <- bismarkEntries[j]
                bismark[[entry]] <- bismark[[entry]][Order]
            }
            rm(Order)
            message('Data sort complete')
            methData[[i]] = GRanges(
                seqnames=bismark[[1]],
                strand=bismark[[6]],
                ranges=IRanges(start=bismark[[2]], width=1),
                methylated=as.integer(round(0.01*bismark[[5]]*bismark[[11]])),
                reads=bismark[[5]])
            rm(bismark)
        } else if (bed_type == 'bismarkCytosine') {
            bismark <- scan(files[i], skip=1, sep="\t",
                            what=list("character", integer(), "character", integer(),
                                      integer(), NULL, NULL))
            message('Begin sorting data')
            bismarkEntries <- c(1,2,3,4,5)
            Order <- with(bismark,order(bismark[[3]],bismark[[1]],bismark[[2]]))
            for (j in 1:length(bismarkEntries)){
                entry <- bismarkEntries[j]
                bismark[[entry]] <- bismark[[entry]][Order]
            }
            rm(Order)
            message('Data sort complete')
            # co-ordinates are 1-based, shift to 0-based
            methData[[i]] = GRanges(
                seqnames=bismark[[1]],
                strand=bismark[[3]],
                ranges=IRanges(start=bismark[[2]] - 1, width=1),
                methylated = bismark[[4]],
                reads = bismark[[4]] + bismark[[5]])
            rm(bismark)
        } else if (bed_type == 'bismarkCov') {
            bismark <- scan(files[i], skip=1, sep="\t",
                            what=list("character", integer(), NULL, NULL,
                                      integer(), integer()))
            message('Begin sorting data')
            bismarkEntries <- c(1,2,5,6)
            Order <- with(bismark,order(bismark[[1]],bismark[[2]]))
            for (j in 1:length(bismarkEntries)){
                entry <- bismarkEntries[j]
                bismark[[entry]] <- bismark[[entry]][Order]
            }
            rm(Order)
            message('Data sort complete')
            # co-ordinates are 1-based, shift to 0-based
            methData[[i]] = GRanges(
                seqnames=bismark[[1]],
                strand = rep('*', length(bismark[[1]])),
                ranges=IRanges(start=bismark[[2]] - 1, width=1),
                methylated = bismark[[5]],
                reads = bismark[[5]] + bismark[[6]])
            rm(bismark)
        } else if (bed_type == 'EPP') {
            bismark <- scan(files[i], skip=1, sep="\t",
                            what=list("character", integer(), NULL, "character",
                                      NULL, "character"))
            message('Begin sorting data')
            bismarkEntries <- c(1,2,4,6)
            Order <- with(bismark,order(bismark[[6]], bismark[[1]],bismark[[2]]))
            for (j in 1:length(bismarkEntries)){
                entry <- bismarkEntries[j]
                bismark[[entry]] <- bismark[[entry]][Order]
            }
            rm(Order)
            message('Data sort complete')
            # co-ordinates are 0-based
            methreads <- vector(length = length(bismark[[4]]))
            totalreads <- vector(length = length(bismark[[4]]))
            for (j in 1:length(bismark[[4]])){
                x <- bismark[[4]][j]
                temp <- strsplit(x, split='/')[[1]]
                methreads[j] <- as.integer(temp[1])
                totalreads[j] <- as.integer(temp[2])
            }
            methData[[i]] = GRanges(
                seqnames=bismark[[1]],
                strand = bismark[[6]],
                ranges=IRanges(start=bismark[[2]], width=1),
                methylated = methreads,
                reads = totalreads)
            rm(bismark, methreads, totalreads, temp, x)

        } else if (bed_type == 'bisSNP') {
            bismark <- scan(files[i], skip=1, sep="\t",
                            what=list("character", integer(), NULL, numeric(),
                                      integer(), 'character',
                                      NULL,NULL,NULL,NULL,NULL))
            message('Begin sorting data')
            bismarkEntries <- c(1,2,4,5,6)
            Order <- with(bismark,order(bismark[[6]],bismark[[1]],bismark[[2]]))
            for (j in 1:length(bismarkEntries)){
                entry <- bismarkEntries[j]
                bismark[[entry]] <- bismark[[entry]][Order]
            }
            rm(Order)
            message('Data sort complete')
            methData[[i]] = GRanges(
                seqnames=bismark[[1]],
                strand=bismark[[6]],
                ranges=IRanges(start=bismark[[2]], width=1),
                methylated=as.integer(round(0.01*bismark[[4]]*bismark[[5]])),
                reads=bismark[[5]])
            rm(bismark)
        }
    }
    return(methData)
}

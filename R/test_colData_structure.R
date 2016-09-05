#' Tests that the files are specified in a workable order for M3D to run
#'
#' M3D requires that the files be specified group by group (all of condition 1,
#' then all of condition 2 etc.). This function tests this, and is used before
#' loading so that the error is caught before time is wasted.
#'
#' @param colData Samples' names plus additional sample information as
#'  character, data.frame or DataFrame.
#' @return Returns TRUE if the colData is structured correctly for downstream
#' analysis, FALSE otherwise
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}


test_colData_structure <- function(colData){
    grps <- as.character(colData$group)
    grps_uniq <- unique(grps)
    test_per_group <- sapply(grps_uniq, function(grp){
        inds <- which(grps==grp)
        all(diff(inds)==1)
    })
    test <- all(test_per_group)
    return(test)
}

saveExpList <- function(scelist, postfix, dir.prefix="RDS/scelist")
{
    n <- paste0(dir.prefix, "_", postfix)
    nms <- lapply(seq_along(scelist), function(i){
        nm <- paste0(n, "_", i)
        saveHDF5SummarizedExperiment(scelist[[i]],
                            dir=nm)
        nm})
    message(n, " wrote on disk (for ", length(scelist), " elements)")
    return(nms)
}

loadExpList <- function(postfix, dir.prefix="RDS/scelist")
{
    scelist <- lapply(seq_along(scelist), function(i){
        loadHDF5SummarizedExperiment(dir=paste0(dir.prefix, "_", postfix, "_", i))})
    return(scelist)
}

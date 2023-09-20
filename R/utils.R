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
    message("Reading ", dir.prefix, "_", postfix, " files")
    files <- list.files(recursive=TRUE)
    dirs <- unique(dirname(files[grep(
            paste0(dir.prefix, "_", postfix, "_[1-9]"), files)]))
    scelist <- lapply(seq_along(dirs), function(i){
        loadHDF5SummarizedExperiment(dir=dir[i])})
    return(scelist)
}

sce <- read10xMultiome("multiome_10xweb_Alzeheimer/outs/", sample.name="Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_", type="HDF5")
sce$replicate <- unlist(lapply(strsplit(sce$Barcode,"-"), function(x){ return(x[[2]])}))

ss <- read.csv("multiome_10xweb_Alzeheimer/outs/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_summary.csv")
ncs <- ss[,grep("Number.of.cells",colnames(ss))]

tr <- table(sce$replicate)
ts <- unlist(lapply(tr, function(t){names(ncs)[ncs==t]}))
cond <- gsub("[.]","",gsub("Number.of.cells..", "",ts))

sce$condition <- NA
for(i in seq_along(cond))
{

    sce$condition[sce$replicate==names(cond)[i]] <- cond[i]
}

## checking numerosity of cells
# ncs
# table(sce$condition)
## ok
allen <- readRDS("allen_23.rds")
sce1 <- assignLabels(sce, allen, "subclass_label")

sce$genotype <-  unlist(lapply(strsplit(sce$condition,"_"), function(x) x[[1]]))
sce$timep <-  unlist(lapply(strsplit(sce$condition,"_"), function(x) x[[2]]))
sce$mouseid <- sce$replicate
sce$replicate <-  unlist(lapply(strsplit(sce$condition,"_"), function(x) x[[3]]))
table(sce$genotype, sce$replicate, sce$timep)
table(sce$genotype, sce$timep)

sce57 <- sce[,sce$timep=="5p7"]
sce57 <- swapAltExp(sce57, name="ATAC")
scelist57 <- splitSCE(sce57, by="genotype")
scelistl57 <- lapply(scelist57, splitSCE, by="replicate")
# scelistlp <- lapply(scelistl, lapply, splitSCE, by="timep")

# scelistl57 <- lapply(scelistl57, swapAltExp, name="ATAC")
scelistl57a <- lapply(seq_along(scelistl57), function(j)
{
    print(j)
    scelistl <- scelistl57[[j]]
    print(scelistl)
    lapply(seq_along(scelistl), function(i)
    {
        sce <- aggregateAcrossCells(scelistl[[i]], use.assay.type="counts",
            ids=DataFrame(label=scelistl[[i]]$SingleR))
        sce$replicate <- names(scelistl)[i]
        sce$genotype <- names(scelist57)[j]
        sce$cnames <- paste0(sce$SingleR, "_", sce$genotype, "_", sce$replicate)
        sce
    })
})


sce57a <- cbind(scelistl57a[[1]][[1]], scelistl57a[[1]][[2]], scelistl57a[[2]][[1]], scelistl57a[[2]][[2]])
##################
##################
##################


scetp <- sce[,sce$timep%in%c("5p7", "17p9")]
scetp <- swapAltExp(scetp, name="ATAC")
scelisttp <- splitSCE(scetp, by="genotype")
scelistltp <- lapply(scelisttp, splitSCE, by="replicate")
## Removing rep2 and rep6 because at the same time point of WT
lapply(scelisttp, function(sce) {table(sce$replicate, sce$timep)})
scelistltp[[2]] <- scelistltp[[2]][-(which(names(scelistltp[[2]]) %in% c("rep2", "rep6")) )]
###

scelistltpa <- lapply(seq_along(scelistltp), function(j)
{
    print(j)
    scelistl <- scelistltp[[j]]
    print(scelistl)
    lapply(seq_along(scelistl), function(i)
    {
        sce <- aggregateAcrossCells(scelistl[[i]], use.assay.type="counts",
                                    ids=DataFrame(label=scelistl[[i]]$SingleR))
        sce$replicate <- names(scelistl)[i]
        sce$genotype <- names(scelisttp)[j]
        sce$cnames <- paste0(sce$SingleR, "_", sce$genotype, "_", sce$replicate)
        sce
    })
})


scetpa <- cbind(scelistltpa[[1]][[1]], scelistltpa[[1]][[2]], scelistltpa[[2]][[1]], scelistltpa[[2]][[2]])

dir.create('G:/project/pregnant_mouse/10x/10x_v3')
setwd('G:/project/pregnant_mouse/10x/10x_v3')
dir.create('rawdata')

library(Seurat)
###################################################################3
#                          load
###################################################################
prefix.dir <- 'H:/serve_backup/project/pregnant_mouse/10x_v3/'
src.raw.G0 <- Read10X(data.dir = paste(prefix.dir,"G0_20200416/outs/raw_feature_bc_matrix/",sep = ""))
src.raw.G14.5 <- Read10X(data.dir = paste(prefix.dir,"G14.5_20200416/outs/raw_feature_bc_matrix/",sep = ""))


colnames(src.raw.G0) <- paste0("G0_20200416",colnames(src.raw.G0))
colnames(src.raw.G14.5) <- paste0("G14.5_20200416",colnames(src.raw.G14.5))


save.image('rawdata/read.10x.G0G14.5G14.5_2nd.RData')



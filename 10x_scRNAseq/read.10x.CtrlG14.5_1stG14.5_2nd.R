dir.create('G:/project/pregnant_mouse/10x/10x_v3')
setwd('G:/project/pregnant_mouse/10x/10x_v3')
dir.create('rawdata')

library(Seurat)
###################################################################3
#                          load
###################################################################
prefix.dir <- 'H:/serve_backup/project/pregnant_mouse/10x_v3/'
src.raw.Ctrl <- Read10X(data.dir = paste(prefix.dir,"Ctrl_20200416/outs/raw_feature_bc_matrix/",sep = ""))
src.raw.G14.5_1st <- Read10X(data.dir = paste(prefix.dir,"G14.5_1st_20200416/outs/raw_feature_bc_matrix/",sep = ""))
src.raw.G14.5_2nd <- Read10X(data.dir = "H:/FileHistory/shuan/DESKTOP-S56EULO/Data/G/project/pregnant_mouse/10x/rawdata/20200308/G14.5_2nd/")

colnames(src.raw.Ctrl) <- paste0("Ctrl_20200416",colnames(src.raw.Ctrl))
colnames(src.raw.G14.5_1st) <- paste0("G14.5_1st_20200416",colnames(src.raw.G14.5_1st))
colnames(src.raw.G14.5_2nd) <- paste0("G14.5_2nd_20200308",colnames(src.raw.G14.5_2nd))

save.image('rawdata/read.10x.CtrlG14.5_1stG14.5_2nd.RData')
# src.raw.Ctrl=CreateSeuratObject(src.raw.Ctrl)
# src.raw.G14.5_1st=CreateSeuratObject(src.raw.G14.5_1st)
# src.raw.G14.5_2nd=CreateSeuratObject(src.raw.G14.5_2nd)

# src.raw.Ctrl$Time="Ctrl"
# src.raw.G14.5_1st$Time="G14.5_1st"
# src.raw.G14.5_2nd$Time="G14.5_2nd"
# src.raw.Ctrl$FACS="Sox9GFP"
# src.raw.G14.5_1st$FACS="Sox9GFP"
# src.raw.G14.5_2nd$FACS="Sox9GFP"

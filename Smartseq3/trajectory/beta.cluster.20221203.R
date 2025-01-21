dir.create('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/')
dir.create('cluster')
setwd('cluster')
load('seu.ref.raw.beta.Glut2LH.RData')
library(future)
plan("multisession", workers = 10)
plan()


source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)

seu.ref.raw <- readRDS('../20220913/cluster/seu.ref.raw.beta.celltype.rds')
pre.ambigous.sym <- readRDS('../20220913/cluster/pre.ambigous.sym.rds')
exclude.cell <- c('m220608_hs_M_zy_20220606_1_lib15_03_044','m220608_hs_M_zy_20220606_1_lib15_03_004',#G0
                  'm220608_hs_M_zy_20220606_1_lib15_03_005',#G0Acss2
                  'm220523_hs_H_zy_20220520_1_lib13_04_082',
                  'm220523_hs_H_zy_20220520_1_lib13_04_093',
                  'm220523_hs_H_zy_20220520_1_lib08_02_132',
                  'm220523_hs_H_zy_20220520_1_lib08_02_094',#G12.5
                  'm220913_hs_L03_M_zy_20220909_04_079',
                  'm220913_hs_L03_M_zy_20220909_04_061',#G6.5
                  "m220420_hs_M_ZY_20220418_3_117","m220420_hs_M_ZY_20220418_3_120",'m220620_hs_G14.5_C57BL6J_4_islet_67' #G14.5 Acss2
                  )
G10.5.cell <- c('m220523_hs_H_zy_20220520_1_lib08_02_031',
                'm220523_hs_H_zy_20220520_1_lib13_04_066',
                'm220523_hs_H_zy_20220520_1_lib13_04_170',
                'm220523_hs_H_zy_20220520_1_lib13_04_090',
                'm220523_hs_H_zy_20220520_1_lib13_04_170')
seu.ref.raw <- subset(seu.ref.raw,cells = colnames(seu.ref.raw)[!seu.ref.raw$SampleName %in% exclude.cell])
seu.ref.raw$Type <- as.character(seu.ref.raw$Type)
seu.ref.raw@meta.data[G10.5.cell,'Type'] <- 'G10.5'
seu.ref.raw$Type <- factor(seu.ref.raw$Type,levels = names(ref.time.colors))
seu.ref.raw <- subset(seu.ref.raw,cells = colnames(seu.ref.raw)[!seu.ref.raw$Type_rep %in% c('Virgin_rep1','Virgin_rep2')])
seu.ref.raw@meta.data[seu.ref.raw$Type_rep=='Virgin_rep3','Rep'] <- 'rep1'
seu.ref.raw@meta.data[seu.ref.raw$Type_rep=='Virgin_rep4','Rep'] <- 'rep2'

seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
seu.ref.raw$Type <- as.character(seu.ref.raw$Type)
seu.ref.raw@meta.data[c('m220523_hs_H_zy_20220520_1_lib13_04_072',
                        'm220523_hs_H_zy_20220520_1_lib13_04_144',
                        'm220523_hs_H_zy_20220520_1_lib13_04_119',
                        'm220523_hs_H_zy_20220520_1_lib08_02_125'
),'Type'] <- 'G14.5'

seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
############
seu.ref.raw <- FindVariableFeatures(seu.ref.raw,selection.method = "vst", nfeatures = 2000)
VariableFeatures(seu.ref.raw) <- setdiff(VariableFeatures(seu.ref.raw),pre.ambigous.sym)

c1Color <- MyName2Col(seu.ref.raw@meta.data[,"Type"],
                      ref.time.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.ref.raw@meta.data[,"Rep"],
                      rep.colors[-1])
c2Color <- as.matrix(c2Color)
# c3Color <- MyName2Col(seu.ref.raw@meta.data[,"SeqDate"],
#                       colors.group)
# c3Color <- as.matrix(c3Color)

cColor <- cbind(c1Color,c2Color)

png('excludeambigous.Vstvar.2000.geneover6000.raw.png',2000,3000)
var2000.row.tree.exclude <-
  MyHeatmap(as.matrix(exp(seu.ref.raw@assays$RNA@data))[VariableFeatures(seu.ref.raw),colnames(seu.ref.raw)],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = cColor,
            ColSideColorsSize = 2,
            return.tree = "row",
            graph = T)
dev.off()

# var2000.row.tree.den <- as.dendrogram(var2000.row.tree)
# pre.ambigous.sym.list <- list()
# pre.ambigous.sym.list[['Alb']] <- labels(var2000.row.tree.den[[2]][[1]])
# pre.ambigous.sym.list[['jun']] <- labels(var2000.row.tree.den[[2]][[2]][[2]][[2]][[1]][[1]])
# pre.ambigous.sym.list[['exo']] <- labels(var2000.row.tree.den[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]])
# pre.ambigous.sym <- unlist(pre.ambigous.sym.list)
# 
# grep('Pnlip',labels(var2000.row.tree.den[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]))
# 
# var2000.row.tree.den <- as.dendrogram(var2000.row.tree.exclude)
pre.cc.sym <- c(labels(as.dendrogram(var2000.row.tree.exclude)[[1]]))
#length(pre.cc.sym)#241

saveRDS(pre.cc.sym,'pre.cc.sym.rds')

VariableFeatures(seu.ref.raw) <- setdiff(VariableFeatures(seu.ref.raw),c(pre.ambigous.sym,pre.cc.sym))
###############
all.var.co <- rownames(MyCo(as.matrix(seu.ref.raw@assays$RNA@data),
                            var.gene = setdiff(VariableFeatures(seu.ref.raw),c(pre.cc.sym,pre.ambigous.sym)),
                            exp.prop.whole.max = 1,
                            exp.prop.whole.min = 0.01,
                            # vector.group = samples.inf.qc$GroupNameFig1,
                            # exp.prop.group.min = 0.1,
                            # exp.prop.group.max = 0.5,
                            cor.method = "rho",
                            cor.cutoff = 0.15,
                            partner.cutoff = 10,
                            refine.cor = T))
length(all.var.co)#233
saveRDS(all.var.co,'all.var.co0.15.10.1.0.01.pc5.rds')


grep('Mt2',all.var.co)
grep('Ucn3',all.var.co)

png('vst2000.rho.rmcc.cor0.15.10.1.0.01.rho.png',2000,3000)
var.co0.15.row.tree <-
  MyHeatmap(as.matrix(exp(seu.ref.raw@assays$RNA@data))[all.var.co,],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = cColor,
            ColSideColorsSize = 2,
            return.tree = "row",
            graph = T)
dev.off()
# 
# var.co0.2.row.tree.den <- as.dendrogram(var.co0.15.row.tree)
# G14.5.heter.sym <- labels(var.co0.2.row.tree.den[[1]][[2]][[1]])
##########
seu.ref.raw <- ScaleData(seu.ref.raw,features = VariableFeatures(seu.ref.raw))
seu.ref.raw <- RunPCA(seu.ref.raw, features = all.var.co)
#seu.ref.raw <- RunPCA(seu.ref.raw,features = VariableFeatures(seu.ref.raw))

pdf('rmcc.cor0.15.1.DimHeatmap.pdf',
    10,40)
DimHeatmap(seu.ref.raw, dims = 1:30, cells = 200, balanced = TRUE)
dev.off()

PCAPlot(seu.ref.raw,dims = c(1,5))


pc.use <- 1:5;#cor0.1
seu.ref.raw <- FindNeighbors(seu.ref.raw, dims = pc.use)
#seu.ref.raw <- FindClusters(seu.ref.raw, resolution = 1.5)
seu.ref.raw <- FindClusters(seu.ref.raw, resolution = 0.2)
seu.ref.raw <- FindClusters(seu.ref.raw, resolution = 0.15)
#seu.ref.raw <- FindClusters(seu.ref.raw, resolution = 3.5)
#seu.ref.raw <- RunUMAP(seu.ref.raw,dims = pc.use)
# 
# pdf('umap.rmcc.vst.cor0.15.1.0.01.10pc13.res0.2.pdf',
#     6,5)
# seu.ref.raw <- SetIdent(seu.ref.raw,value = seu.ref.raw$RNA_snn_res.0.2)
# DimPlot(seu.ref.raw,
#         reduction = "umap",
#         cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#         label.size = 6,
#         sizes.highlight = 4,
#         pt.size = 1,
#         label = T)  
# 
# seu.ref.raw <- SetIdent(seu.ref.raw,value = seu.ref.raw$Type)
# DimPlot(seu.ref.raw,
#         reduction = "umap",
#         cols = c(ref.time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#         label.size = 6,
#         sizes.highlight = 4,
#         pt.size = 1,
#         label = F)
# seu.ref.raw <- SetIdent(seu.ref.raw,value = seu.ref.raw$beta)
# DimPlot(seu.ref.raw,
#         reduction = "umap",
#         cols = c(ref.time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#         label.size = 6,
#         sizes.highlight = 4,
#         pt.size = 1,
#         label = F)
# dev.off()
# for(tmp.perplexity in seq(50,150,10)){
#   seu.ref.raw <- RunTSNE(seu.ref.raw,
#                          dims = pc.use,
#                          #perplexity= round((30+ncol(seu.ref.raw)/100)),
#                          perplexity= tmp.perplexity,
#                          check_duplicates = F)
#   pdf(paste('perplexity',tmp.perplexity,'tsne.rmcc.vst.cor0.15.1.0.01.10.pc8.res0.1.res0.2.pdf',sep = ''),
#       6,5)
#   seu.ref.raw <-
#     SetIdent(seu.ref.raw, value = seu.ref.raw$RNA_snn_res.0.1)
#   print(DimPlot(seu.ref.raw,
#           reduction = "tsne",
#           cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#           label.size = 6,
#           sizes.highlight = 4,
#           pt.size = 2,
#           label = T))
#   seu.ref.raw <-
#     SetIdent(seu.ref.raw, value = seu.ref.raw$RNA_snn_res.0.2)
#   print(DimPlot(seu.ref.raw,
#           reduction = "tsne",
#           cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#           label.size = 6,
#           sizes.highlight = 4,
#           pt.size = 2,
#           label = T))
#   seu.ref.raw <-
#     SetIdent(seu.ref.raw, value = seu.ref.raw$beta)
#   print(DimPlot(seu.ref.raw,
#           reduction = "tsne",
#           cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#           label.size = 6,
#           sizes.highlight = 4,
#           pt.size = 2,
#           label = T))
#   seu.ref.raw <-
#     SetIdent(seu.ref.raw, value = seu.ref.raw$Type)
#   print(DimPlot(seu.ref.raw,
#           reduction = "tsne",
#           cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#           label.size = 6,
#           sizes.highlight = 4,
#           pt.size = 2,
#           label = F))
#   dev.off()
#   
# }
seu.ref.raw <- RunTSNE(seu.ref.raw,
                       dims = pc.use,
                       perplexity= round((30+ncol(seu.ref.raw)/100)),
                    #   perplexity= 100,
                       check_duplicates = F)
 
pdf('tsne.rmcc.vst.cor0.15.1.0.01.10.pc5.res0.1.res0.2.pdf',
    6,5)
seu.ref.raw <-
  SetIdent(seu.ref.raw, value = seu.ref.raw$RNA_snn_res.0.15)
DimPlot(seu.ref.raw,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = T)
seu.ref.raw <-
  SetIdent(seu.ref.raw, value = seu.ref.raw$RNA_snn_res.0.2)
DimPlot(seu.ref.raw,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = T)
seu.ref.raw <-
  SetIdent(seu.ref.raw, value = seu.ref.raw$beta)
DimPlot(seu.ref.raw,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = T)
seu.ref.raw <-
  SetIdent(seu.ref.raw, value = seu.ref.raw$Type)
DimPlot(seu.ref.raw,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = F)
dev.off()

Mybarplot(seu.ref.raw@meta.data,c1 = 'Type',c2 = 'RNA_snn_res.0.2',xlim = 20,cols = time.colors)
########marker########
marker.sym <- c(
  "Neurod1",#endocrine
  "Chga",
  "Chgb",
  "Mki67",
  "Ins1",
  'Ins2',
  "Ppy",
  "Ghrl",
  "Nkx6-1",
  "Ucn3",
  "Mafa",
  "Slc2a2",
  "Arx",
  "Gcg",
  "Mafb",
  "Sst",
  "Hhex",
  'Prlr',
  'Oxtr',
  'Ovol2',
  'Id4',
  'Tph1','Tph2','Gbp8','Sftpd',
  'Il1r1',
  'Mafb',
  'Hsbp1',
  'Tspan8',
  'Cd81',
  'Gpx3'
)
length(marker.sym)#
marker.sym <- marker.sym[marker.sym %in% rownames(seu.ref.raw)]
length(marker.sym)#
marker.sym <- unique(marker.sym)
length(marker.sym)#57
total.count <- length(marker.sym)
run.count <- 1
Glut2LH.sym <- c('Chgb','Ovol2','Hsbp1','Ttr','Ucn3','Slc2a2','Cd81','Tspan8')
seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')

pdf("FigS3a.tsne.log.marker.rmcc.vst.cor0.15.1.0.01.10pc5.pdf",
    6,
    7)
Myseuratmarker(seu.ref.raw,marker.sym = Glut2LH.sym,reduction = 'tsne',pt.size= 1)
dev.off()

##########cell Type#######
# seu.ref.raw@meta.data$subgroup <- '/'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==3,'subgroup'] <- 'beta1'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==1,'subgroup'] <- 'beta2'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==2,'subgroup'] <- 'beta3'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==0,'subgroup'] <- 'beta4'

seu.ref.raw@meta.data$beta <- '/'
seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==3,'beta'] <- 'Glut2L'
seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2 %in% c(2,1,0),'beta'] <- 'Glut2H'

# seu.ref.raw@meta.data$subgroup2 <- '/'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==3,'subgroup2'] <- 'beta3'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==1,'subgroup2'] <- 'beta1'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2 %in% c(0,2),'subgroup2'] <- 'beta2'
# 
# 
# seu.ref.raw@meta.data$subgroup3 <- '/'
# seu.ref.raw@meta.data[seu.ref.L$SampleName[seu.ref.L$subgroup=='group1'],'subgroup3'] <- 'beta3'
# seu.ref.raw@meta.data[seu.ref.L$SampleName[seu.ref.L$subgroup=='group2'],'subgroup3'] <- 'beta4'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2==1,'subgroup3'] <- 'beta1'
# seu.ref.raw@meta.data[seu.ref.raw$RNA_snn_res.0.2%in% c(2,0),'subgroup3'] <- 'beta2'


# group.col <- time.colors[c(8,5,6,9)]
# names(group.col) <- c('beta1','beta2','beta3','beta4')
# 
# pdf('FS1.all.subgroup.barplot.pdf',8,8)
# Mybarplot(seu.ref.raw@meta.data,c1 = 'Type',c2 = 'subgroup2',xlim = 20,cols = group.col)
# dev.off()
# pdf('FS1.all.subgroup3.barplot.pdf',8,8)
# Mybarplot(seu.ref.raw@meta.data,c1 = 'Type',c2 = 'subgroup3',xlim = 20,cols = group.col)
# dev.off()  
##############
celltype.col <- c('#1f78b4','#ff7f00')
names(celltype.col) <- c('Glut2H','Glut2L')

celltype.col2 <- c("#E08698","#6a3d9a","#5F9EA0", "#A65628")
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2')

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/FS3e.all.Glut2LH.ratio.barplot.pdf',8,8)
Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'beta',xlim = 20,cols = c('gray70','red3'))
dev.off()

# group.col <- colors.group[1:4]
# names(group.col) <- c('beta1','beta2','beta3','beta4')
seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
seu.ref.L <- readRDS('Glut2L/seu.ref.L.endo.celltype.rds')

seu.ref.L$subgroup <- factor(seu.ref.L$subgroup,levels = c('group2','group1'))
pdf('Glut2L/FS3.all.Glut2L.ratio.barplot.pdf',12,8)
Mybarplot(seu.ref.L@meta.data,c1 = 'Type',c2 = 'subgroup',xlim = 20,cols = celltype.col2[2:1])
dev.off()


seu.ref.raw$beta2 <- as.character(seu.ref.raw$beta)
seu.ref.raw@meta.data[seu.ref.L$SampleName[seu.ref.L$subgroup=='group1'],'beta2'] <- 'Glut2L_1'
seu.ref.raw@meta.data[seu.ref.L$SampleName[seu.ref.L$subgroup=='group2'],'beta2'] <- 'Glut2L_2'

table(seu.ref.raw$beta2)
time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")
celltype.col2 <- time.col[c(2,9,17,12)]
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2')


time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")
MyPlotColor(time.col,length(time.col))
ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)

MyPlotColor(ref.time.colors2,length(ref.time.colors2))


seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
p.umap <- MySeuratv3TSNE10x2Gg(seu.ref.raw,
                               seu.ref.raw@meta.data
)
pdf("G:/lab/Article/heshuang/BYLW/sm3/ref/F1a.left.tsne.Glut2LH.cellType.pdf",
    14.5,
    8)
print(p.umap +
        scale_color_manual(values = celltype.col2) +
        geom_point(aes(x = x.pos,
                       y = y.pos,
                       col =beta2#,
                       #shape = State
        ),
        size = 2
        )+ theme(aspect.ratio=1))

dev.off()

pdf("G:/lab/Article/heshuang/BYLW/sm3/ref/F1a.right.tsne.Type.pdf",
    14.5,
    8)
print(p.umap +
        scale_color_manual(values = ref.time.colors2) +
        geom_point(aes(x = x.pos,
                       y = y.pos,
                       col =Type#,
                       #shape = State
        ),
        size = 2
        )+theme(aspect.ratio=1))

dev.off()
#########rep#########
seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
sample.inf.keep.cp <- seu.ref.raw@meta.data
age.list <- names(table(sample.inf.keep.cp$Type))


pdf("FigS2d.pc5.cor0.15.1.0.01.res0.1.tsne.rep.batch.age.pdf",
    14.5,
    8)
for (i in age.list){
  print(i)
  sample.inf.keep.cp$Rep1 <- "*"
  sample.inf.keep.cp[sample.inf.keep.cp$Type %in% i,"Rep1"] <- as.character(sample.inf.keep.cp[sample.inf.keep.cp$Type == i,"Rep"])
  
  sample.inf.keep.cp$Rep1 <- factor(sample.inf.keep.cp$Rep1,
                                    levels = c("*",
                                               "rep1",
                                               "rep2",
                                               "rep3",
                                               "rep4"#,
                                               # "Fucci-rep1"
                                    ))
  
  
  p.tmp <- MySeuratv3TSNE10x2Gg(seu.ref.raw,
                                sample.inf.keep.cp)
  p.tmp_ <- p.tmp
  p.tmp_$data$Rep1_ <- p.tmp_$data$Rep1
  p.tmp_$data$Rep1_ <-  factor(p.tmp_$data$Rep1_,
                               levels = c("*",
                                          #"Fucci-rep1",
                                          
                                          "rep4",
                                          "rep3",
                                          "rep2",
                                          "rep1"))
  
  p.tmp_$data <- p.tmp_$data[order(p.tmp_$data$Rep1_),]
  
  
  print(p.tmp_ +
          scale_color_manual(values = rep.colors)+
          geom_point(aes(x = x.pos,
                         y = y.pos,
                         col = Rep1#,
                         #shape = State
          ),
          size = 5
          )+
          labs(title = i)+theme(aspect.ratio=1))
}
dev.off()
########
source('G:/pcatest/MyFunction.R')
beta.si.tab <- MyReadDelim('seu.ref.beta.meta.tab')
pdf('SF2bc.QC.over6000.genecount.rawcount.pdf')

hist(seu.ref.raw$Genecount_all,breaks = 120,
     lwd = 4,xlim = c(6000,12000),ylim = c(0,120)
)


hist(seu.ref.raw$Raw_count,breaks = 120,xlim = c(50000,2000000),ylim = c(0,140),
     lwd=4#,
     #axes=FALSE
)

dev.off()
########
MyWriteTable(table(seu.ref.raw$Type_rep,seu.ref.raw$beta),row.names = T,
             'celltype.Glut2LH.rep.count.tab')
MyWriteTable(table(seu.ref.raw$Type,seu.ref.raw$beta),row.names = T,
             'celltype.Glut2LH.count.tab')
MyWriteTable(table(seu.ref.raw$Type,seu.ref.raw$Rep),row.names = T,
             'type.rep.count.tab')

MyWriteTable(table(seu.ref.raw$Type),row.names = T,
             'celltype.count.tab')
MyWriteTable(table(seu.ref.raw$beta2),row.names = T,
             'celltype.beta2.count.tab')

MyWriteTable(seu.ref.raw@meta.data,'seu.ref.beta.meta.tab')
########
pre.alb.sym.tab <- genes.inf.input[as.character(unlist(pre.ambigous.sym.list)),]
pre.alb.sym.tab$Type <- 'Alb'
pre.alb.sym.tab[pre.ambigous.sym.list$exo,'Type'] <- 'exo'
pre.alb.sym.tab[pre.ambigous.sym.list$jun,'Type'] <- 'jun'

MyWriteTable(pre.alb.sym.tab,
             'pre.ambigous.sym.tab')

saveRDS(pre.ambigous.sym,'pre.ambigous.sym.rds')

##########
seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
saveRDS(seu.ref.raw,'seu.ref.raw.beta.celltype.rds')
##########
seu.ref.raw <- readRDS('seu.ref.raw.beta.celltype.rds')
rm(seu.ref.raw,seu.ref.L)
save.image('seu.ref.raw.beta.Glut2LH.RData')




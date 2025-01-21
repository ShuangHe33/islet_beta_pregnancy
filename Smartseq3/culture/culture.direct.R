dir.create('direct')
dir.create('direct/20221120/')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/project/20220913/direct/20221120/')
load('direct.20230704.RData')
load('direct.20221120.RData')
load('direct.20230606.RData')

preg.c1.gene <- sub('_','-',preg.gene.si.tab$SymbolDedu[preg.gene.si.tab$cluster=='cc'])

rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
seu.vitro.nor.list$Nif <- subset(seu.vitro.nor.list$Nif,cells = setdiff(colnames(seu.vitro.nor.list$Nif),exclude.cell))
seu.vitro.nor.list$Nif@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalTp0.1mExcludeGene(as.matrix(seu.vitro.nor.list$Nif@assays$RNA@counts),
                                                                                                     genes.inf.input[rownames(seu.vitro.nor.list$Nif@assays$RNA@counts),'GeneLength'],
                                                                                                     c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/project/20220913/seu.vitro.nor.list.rds')
seu.vitro.nor.list[['sm3.Nif.A485.only']] <- subset(seu.vitro.nor.list$sm3.Nif.A485,cells = colnames(seu.vitro.nor.list$sm3.Nif.A485)[
  seu.vitro.nor.list$sm3.Nif.A485$Type %in% c('Ctrl_DMSO','Ctrl_A485_Nif','Preg_DMSO','Preg_A485_Nif')
])

seu.vitro.nor.list[['sm3.NifA485.all']] <- merge(seu.vitro.nor.list$Nif,seu.vitro.nor.list[c('sm3.A485','sm3.Nif.A485')])
table(seu.vitro.nor.list[['sm3.NifA485.all']]$Type)
seu.vitro.nor.list[['sm3.NifA485.all']]$Type <- factor(seu.vitro.nor.list[['sm3.NifA485.all']]$Type,levels = c('Ctrl_DMSO','Ctrl_A485','Ctrl_20um_Nif','Ctrl_A485_Nif','Preg_DMSO','Preg_A485','Preg_20um_Nif','Preg_A485_Nif'))


setwd('../../')

row.tree.list <- list()
seu <- 'sm3.A485'
seu <- 'sm3.Nif.A485.only'
seu <- 'sm3.NifA485.all'
seu <- 'sm3.A485.add'
seu <- 'sm3.Nif.add'
seu <- 'sm3.A485.addpregA485'
seu <- 'sm3.Nif.add'
seu.vitro.nor.list <- readRDS('../../seu.vitro.nor.list.rds')
seu.vitro.nor.list[['sm3.Nif.batch1']] <- subset(seu.vitro.nor.list$sm3.Nif.add,cells=colnames(seu.vitro.nor.list$sm3.Nif.add)[
  seu.vitro.nor.list$sm3.Nif.add$SeqDate=='20220613'
])

seu.vitro.nor.list[['sm3.Nif.batch2']] <- subset(seu.vitro.nor.list$sm3.Nif.add,cells=colnames(seu.vitro.nor.list$sm3.Nif.add)[
  seu.vitro.nor.list$sm3.Nif.add$SeqDate=='20230606'
  ])


seu <- 'sm3.Nif.batch1'
seu <- 'sm3.Nif.batch2'
for (seu in names(seu.vitro.nor.list)[c(7)]) {
  seu.vitro.nor.list[[seu]] <- 
  FindVariableFeatures(seu.vitro.nor.list[[seu]], selection.method = "vst", nfeatures = 2000)
  VariableFeatures(seu.vitro.nor.list[[seu]]) <- setdiff(VariableFeatures(seu.vitro.nor.list[[seu]]),c(pre.ambigous.sym,
                                                                                                        beta.cc.sym
                                                                                                       ))
  c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Rep"],
                        rep.colors[-1])
  c2Color <- as.matrix(c2Color)

  cColor <- cbind(c1Color,c2Color)
  
  png(paste('direct/20221120/',seu,"vst.2000.heatmap.png",sep = ""),2000,2000)
  row.tree.list[[seu]] <-
  MyHeatmap(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[VariableFeatures(seu.vitro.nor.list[[seu]]),colnames(seu.vitro.nor.list[[seu]])],
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
}
labels(as.dendrogram(row.tree.list$Nif)[[2]][[1]][[1]][[1]])
sm3.A485.cc.gene <- labels(as.dendrogram(row.tree.list$sm3.A485)[[2]][[2]][[1]])
sm3.A485.cc.gene2 <- labels(as.dendrogram(row.tree.list$sm3.A485)[[2]][[2]][[2]])

# VariableFeatures(seu.vitro.nor.list$sm3.Nif.A485.only) <- setdiff(VariableFeatures(seu.vitro.nor.list[[seu]]),
#                                                                   labels(as.dendrogram(row.tree.list$sm3.Nif.A485.only)[[2]][[1]])
#                                                                   )

dir.create('go_kegg')

MyGOwritetable(genes.inf.input[sm3.A485.cc.gene,1],pvalue = 1,'go_kegg/go.sm3.A485.cc.gene.221.tab')
MyGOwritetable(genes.inf.input[sm3.A485.cc.gene2,1],pvalue = 1,'go_kegg/go.sm3.A485.cc.gene.222.tab')

var.co.list <- list()
for(seu in names(seu.vitro.nor.list)[c(1,6)]){
  print(seu)
  var.co.list[[seu]] <- 
    rownames(MyCo(as.matrix(seu.vitro.nor.list[[seu]]@assays$RNA@data),
                  var.gene = VariableFeatures(seu.vitro.nor.list[[seu]]),
                                    
                  exp.prop.whole.max = 1,
                  exp.prop.whole.min = 0.03,
                  # vector.group = samples.inf.qc$GroupNameFig1,
                  # exp.prop.group.min = 0.1,
                  # exp.prop.group.max = 0.5,
                  cor.method = "rho",
                  cor.cutoff = 0.2,
                  partner.cutoff = 10,
                  refine.cor = T))
  c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Rep"],
                        rep.colors[-1])
  cColor <- cbind(c1Color,c2Color)
  
  png(paste('direct/20221120/',seu,".cor0.2.10.1.0.03.heatmap.png",sep = ""),2000,2000)
  row.tree.list[[paste(seu,'cor0.2',sep = '')]] <- 
         MyHeatmap(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data)[var.co.list[[seu]],]),
                   type = "log.row.relat",
                   hc.c.data.type = "log.row.relat",
                   hc.r.data.type = "log.row.relat",
                   c.cov.method = "s",
                   r.cov.method = "s",
                   ColSideColorsSize = 4,
                   ColSideColors = cColor,
                   c.hc.method = "ward.D2",
                   r.hc.method = "ward.D2",
                   return.tree = "row",
                   graph = T)
  dev.off()
}

# var.co.list$Nif <- setdiff(var.co.list$Nif,c(labels(as.dendrogram(row.tree.list$Nifcor0.2)[[2]][[2]][[1]])
#                                              ))

var.co.list$sm2.A485 <- setdiff(var.co.list$sm2.A485,c(
                                             labels(as.dendrogram(row.tree.list$sm2.A485cor0.2)[[2]][[2]][[2]][[2]])
))


var.co.list$sm3.A485 <- setdiff(var.co.list$sm3.A485,c(
  labels(as.dendrogram(row.tree.list$sm3.A485cor0.2)[[1]])
))


var.co.list$sm3.Nif.A485.only <- setdiff(var.co.list$sm3.Nif.A485.only,c(
  labels(as.dendrogram(row.tree.list$sm3.Nif.A485.onlycor0.2)[[2]][[1]])
))

sm3.A485.gene <- labels(as.dendrogram(row.tree.list$sm3.A485cor0.2)[[1]])
sm3.A485.gene <- labels(as.dendrogram(row.tree.list$sm3.A485cor0.2)[[1]])

labels(as.dendrogram(row.tree.list$sm3.NifA485.allcor0.2)[[2]][[1]][[1]])


MyGOwritetable(genes.inf.input[sm3.A485.gene,1],pvalue = 1,'direct/20221120/go_kegg/go.sm3.A485.gene.cor.1.tab')

seu.vitro.nor.list <- readRDS('../../seu.vitro.nor.list.rds')


var.co.list$sm3.A485.add <- setdiff(var.co.list$sm3.A485.add,c(
  labels(as.dendrogram(row.tree.list$sm3.A485.addcor0.2)[[1]][[1]])
))

var.co.list$sm3.A485.addpregA485 <- setdiff(var.co.list$sm3.A485.addpregA485,c(
  labels(as.dendrogram(row.tree.list$sm3.A485.addpregA485cor0.2)[[1]])
))


var.co.list$sm3.Nif.add <- setdiff(var.co.list$sm3.Nif.add,c(
  labels(as.dendrogram(row.tree.list$sm3.Nif.addcor0.2)[[1]][[1]])
))

#########cc########
cc.col.tree.list <- list()
for (seu in names(seu.vitro.nor.list)[c(1,3)]) {
  
  c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[,"Rep"],
                        rep.colors[-1])
  c2Color <- as.matrix(c2Color)
  
  cColor <- cbind(c1Color,c2Color)
  
  png(paste('direct/20221120/',seu,".col.cc.heatmap.png",sep = ""),2000,2000)
  cc.col.tree.list[[seu]] <-
    MyHeatmap(as.matrix(MyGeneExp(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data)[rownames(seu.vitro.nor.list[[seu]]@assays$RNA@data) %in% preg.c1.gene,colnames(seu.vitro.nor.list[[seu]])],1,1)),
              type = "log.row.relat",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = cColor,
              ColSideColorsSize = 2,
              return.tree = "col",
              graph = T)
  dev.off()
}
cc.cell.list <- list()
cc.cell.list[['sm3.A485']] <- c(labels(as.dendrogram(cc.col.tree.list$sm3.A485)[[1]]),labels(as.dendrogram(cc.col.tree.list$sm3.A485)[[2]][[2]][[1]]))
cc.cell.list[['Nif']] <- c(labels(as.dendrogram(cc.col.tree.list$Nif)[[1]]),labels(as.dendrogram(cc.col.tree.list$Nif)[[2]][[1]]))

for (seu in names(seu.vitro.nor.list)[c(1,3)]) {
  seu.vitro.nor.list[[seu]]$proliferation <- 'quiescent'
  seu.vitro.nor.list[[seu]]@meta.data[cc.cell.list[[seu]],'proliferation'] <- 'proliferative'
}

########fmnn########
library(batchelor)
fMNN.res <- fastMNN(as.matrix(seu.vitro.nor.list$sm3.Nif.add@assays$RNA@data[,seu.vitro.nor.list$sm3.Nif.add$SeqDate%in%"20220613"]),
                    as.matrix(seu.vitro.nor.list$sm3.Nif.add@assays$RNA@data[,seu.vitro.nor.list$sm3.Nif.add$SeqDate%in%"20230606"]),
                    d=30,
                    subset.row = unique(c(var.co.list$sm3.Nif.batch1
                                          #,var.co.list$sm3.Nif.batch2
                                          )))

seu.vitro.nor.list$sm3.Nif.add@reductions$pca@cell.embeddings <- fMNN.res@assays@data$reconstructed@seed@components[seu.vitro.nor.list$sm3.Nif.add$SampleName,]
colnames(seu.vitro.nor.list$sm3.Nif.add@reductions$pca@cell.embeddings) <- paste("PC_",1:30,sep="")

seu <- 'sm3.Nif.add'
###########
seu <- 'Nif'
marker.sym <- c('Ucn3','Cd81','Tspan8','Gpx3','Prlr','Pdx1','Nkx6.1','Ttr','Chgb','Tph1','Ivd')
for(seu in names(seu.vitro.nor.list)[c(1,6)]){
  print(seu)
# pdf(paste('direct/20221120/',seu,'.pref.c1.cor0.2.1.0.03.pc6.res0.5.pdf',sep = ''))
seu.vitro.nor.list[[seu]] <- ScaleData(seu.vitro.nor.list[[seu]],features = c(var.co.list[[seu]]#,preg.c1.gene
                                                                              ))
 # seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = c(var.co.list[[seu]],preg.c1.gene))
  seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = var.co.list[[seu]])
  seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = setdiff(var.co.list[[seu]],sm3.A485.cc.gene2))
  seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = setdiff(var.co.list[[seu]],sm3.A485.gene))
  
  # seu.vitro.nor.list[[seu]] <- ScaleData(seu.vitro.nor.list[[seu]],features = VariableFeatures(seu.vitro.nor.list[[seu]]))
  # seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features =setdiff(VariableFeatures(seu.vitro.nor.list[[seu]]),
  #   labels(as.dendrogram(row.tree.list$sm2.A485cor0.2)[[2]][[2]][[2]][[2]])
  # ))
  # 
  # # 
  # seu.vitro.nor.list[[seu]] <- ScaleData(seu.vitro.nor.list[[seu]],features = DEG.mrge.list[[seu]])
  # seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = DEG.mrge.list[[seu]])
  # # 
  pc.use <- 1:6
  seu.vitro.nor.list[[seu]] <- FindNeighbors(seu.vitro.nor.list[[seu]], dims = pc.use)
  seu.vitro.nor.list[[seu]] <- FindClusters(seu.vitro.nor.list[[seu]], resolution = 0.5)

  seu.vitro.nor.list[[seu]] <- RunUMAP(seu.vitro.nor.list[[seu]],dims = pc.use)
#   #seu.vitro.nor.list[[seu]] <- RunTSNE(seu.vitro.nor.list[[seu]],dims = pc.use)
#   seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$RNA_snn_res.0.5)
#   # 
#   # print(DimPlot(seu.vitro.nor.list[[seu]],
#   #               reduction = "umap",
#   #               cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
#   #               label.size = 6,
#   #               sizes.highlight = 4,
#   #               pt.size = 2,
#   #               label = F))
seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$Type)
#   print(DimPlot(seu.vitro.nor.list[[seu]],
#                 reduction = "umap",
#                 cols = c(time.colors,brewer.pal(8,"Set1")),
#                 label.size = 6,
#                 sizes.highlight = 4,
#                 pt.size = 2,
#                 label = F))
  print(DimPlot(seu.vitro.nor.list[[seu]],
                reduction = "umap",
                cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 3,
                label = F))
# 
# # print(DimPlot(seu.vitro.nor.list[[seu]],
# #                 reduction = "tsne",
# #                 cols = c(ko.time.col[as.character(unique(seu.vitro.nor.list[[seu]]$Type))],'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
# #                 label.size = 6,
# #                 sizes.highlight = 4,
# #                 pt.size = 2,
# #                 label = F))
#   
#   
  seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$SeqDate)
  print(DimPlot(seu.vitro.nor.list[[seu]],
                reduction = "pca",
                cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 3,
                label = F))
   
#   Myseuratmarker(seu.vitro.nor.list[[seu]],marker.sym = c('Mki67','Cdt1','Gmnn','Ucn3','Slc2a2','Cd81'),reduction = 'pca')
#   dev.off()
}
pdf('Ctrl.A485.tmp.marker.pdf',6,7)
Myseuratmarker(subset(seu.vitro.nor.list[[seu]],cells = colnames(seu.vitro.nor.list[[seu]])[seu.vitro.nor.list[[seu]]$Type %in% c('Ctrl_DMSO','Ctrl_A485')]),c('Ins1','Prlr','Gcg','Irx1','Mafb','Gpx3'))
dev.off()
UMAPPlot(subset(seu.vitro.nor.list[[seu]],cells = colnames(seu.vitro.nor.list[[seu]])[seu.vitro.nor.list[[seu]]$Type %in% c('Ctrl_DMSO','Ctrl_A485')])
         ,cols=time.colors,pt.size=3)
###########
seu.vitro.nor.list$Nif$Type <- factor(seu.vitro.nor.list$Nif$Type,levels = names(table(seu.vitro.nor.list$Nif$Type))[table(seu.vitro.nor.list$Nif$Type)!=0])
seu.vitro.nor.list$sm3.A485$Type <- factor(seu.vitro.nor.list$sm3.A485$Type,levels = c('Ctrl_DMSO','Ctrl_A485','Preg_DMSO','Preg_A485'))
seu.vitro.nor.list$sm3.Nif.A485.only$Type <- factor(seu.vitro.nor.list$sm3.Nif.A485.only$Type,levels = c('Ctrl_DMSO','Ctrl_A485_Nif','Preg_DMSO','Preg_A485_Nif'))

seu.vitro.nor.list <- readRDS('seu.vitro.nor.list.rds')

seu <- "sm3.A485.add"
seu <- 'sm3.A485.addpregA485'

seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,1] <- -seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,1]
seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,2] <- -seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,2]
seu <- 'sm3.Nif.add'
seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,3] <- -seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,3]
seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,2] <- -seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,2]

for(seu in names(seu.vitro.nor.list)[1:6]){
  pca.df.list[[seu]] <- cbind(Embeddings(seu.vitro.nor.list[[seu]],'pca'),seu.vitro.nor.list[[seu]]@meta.data[,c('SampleName','Type')])
  #pca.df.list[[seu]] <- cbind(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,],seu.ref.KO.list2[[seu]]@meta.data[,c('SampleName','Type')])
  pca.df.list[[seu]] <- pca.df.list[[seu]][pca.df.list[[seu]]$Type %in% names(table(seu.vitro.nor.list[[seu]]$Type)),]
  
  ko.coord.list[[seu]] <- Mycoord(pca.df.list[[seu]],dims=c(1,2),#seu <- "sm3.A485.add" 1,3# sm3.Nif.add 2,3
                                  class = 'Type',
                                  names(table(seu.vitro.nor.list[[seu]]$Type)))
}

#######QC#####
pdf('G:/lab/Article/heshuang/BYLW/sm3/culture/genecount.boxplot.all.pdf',4,5)
boxplot(c(seu.vitro.nor.list$Nif$Genecount_all,seu.vitro.nor.list$sm3.A485$Genecount_all),lwd=2,main='gene number')
box(lwd=3)
boxplot(c(seu.vitro.nor.list$Nif$Raw_count,seu.vitro.nor.list$sm3.A485$Raw_count),lwd=2,main='read count')
box(lwd=3)
dev.off()
###########

width=15;height=12
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.circle0.7.0.6.shape.dirct.rmcc.pc12.pdf',sep = ''),width,height)
  seu.vitro.nor.list[[seu]]$Type2 <- seu.vitro.nor.list[[seu]]$Type
  seu.vitro.nor.list[[seu]]$Type2 <- factor(seu.vitro.nor.list[[seu]]$Type2,levels = names(table(seu.vitro.nor.list[[seu]]$Type))[1:3])
  p.pca <- MySeuratDR2Gg2(seu.vitro.nor.list[[seu]],seu.vitro.nor.list[[seu]]@meta.data,reduction.use = 'pca',
                          reduction.key = 'PC',estimate.variation.explain.percentage = T,
                          x.dim = 1,y.dim = 2)
  p.pca$data_ <- p.pca$data
  p.pca$data_$Type_ <- p.pca$data_$Type
  p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(seu.vitro.nor.list[[seu]]$Type)))
  p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]
  
  plot(p.pca+
         scale_color_manual(values =c('skyblue1','sienna1','deeppink1','purple1')[c(2,1,4,3)]) +
         scale_shape_manual(values = c(19,17,8,15,1,6))+
         geom_point(aes(x = -x.pos,
                        y = y.pos,
                        col = Type,
                        shape = Rep
         ),
         size =4)
          +theme(aspect.ratio = 1)+stat_ellipse(aes(x = -x.pos,
                                                     y = y.pos,col=Type2,size=I(2.5)),level=0.6,type='norm')+
          geom_point(aes(x=-ko.coord.list[[seu]][1,1],y=ko.coord.list[[seu]][1,2],size=6),colour='black')+
          geom_point(aes(x=-ko.coord.list[[seu]][2,1],y=ko.coord.list[[seu]][2,2],size=6),colour='black')+
          geom_point(aes(x=-ko.coord.list[[seu]][3,1],y=ko.coord.list[[seu]][3,2],size=6),colour='black')
        )
  
  dev.off()
}

# library(ggord)
# library(ggplot2)
# for(seu in names(seu.vitro.nor.list)){
#   pca.res <- FactoMineR::PCA(t(as.matrix(seu.vitro.nor.list[[seu]]@assays$RNA@data[setdiff(var.co.list[[seu]],sm3.A485.gene),])))
#   pca.res$ind$coord[,1] <- seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,'PC_1']
#   pca.res$ind$coord[,2] <- -seu.vitro.nor.list[[seu]]@reductions$pca@cell.embeddings[,'PC_2']
#   pdf(paste(seu,'circle.0.5.pdf'),15,12)
#   print(ggord(pca.res,seu.vitro.nor.list[[seu]]$Type,col = c(time.colors[1:4]),arrow=NULL,vec_ext =0,txt=NULL,ellipse=T,ellipse_pro=0.5,
#               poly = F,alpha=0.5,size=I(8))+ylim(c(-10,13))+xlim(c(-18,13))+scale_shape_manual(values = c(19,17,17,17))+
#           theme_bw() + 
#           theme(panel.grid.major = element_blank(), 
#                 panel.grid.minor = element_blank())+theme(aspect.ratio = 1)+aes(size=I(2.5))+
#           theme(panel.border = element_rect(size = 4,
#                                             colour = "black"))+
#           theme(axis.text = element_text(size = 20, colour = "black")) +
#           theme(axis.title.x = element_text(size = 40, colour = "black")) +
#           theme(axis.title.y = element_text(size = 40, colour = "black")) +
#           theme(legend.text = element_text(size = 45, colour = "black")) +
#           theme(legend.title = element_text(size = 45, colour = "black")) +
#           theme(panel.border = element_rect(size = 4,
#                                             colour = "black")))
#   dev.off()
# }
##########KOscore#######
for(seu in names(seu.ref.KO.list2)[1:3]){
  x <- ko.coord.list[[seu]][3,]-ko.coord.list[[seu]][1,]
  x_norm <- x/sqrt(sum(x^2))
  
  seu.vitro.nor.list[[seu]]$cg.koscore <- as.matrix(Embeddings(seu.vitro.nor.list[[seu]],'pca')[,1:2]) %*% as.matrix(x_norm)
  #seu.ref.KO.list2[[seu]]$pseudotime <- as.matrix(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,1:2]) %*% as.matrix(x_norm)
}

seu.vitro.nor.list[[seu]]$pseudotime <- seu.vitro.nor.list[[seu]]$cg.koscore

direct.select.si.list <- list()
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
 # seu.vitro.nor.list[[seu]]$pseudotime <- -Embeddings(seu.vitro.nor.list[[seu]],'pca')[,2]
  width=10;height=6
  pdf(paste('direct/20221120/',seu,'.231002.direct.pseudobox.2.pdf',sep = ''),width,height)
  print(MyPseudotimebox(seu.vitro.nor.list[[seu]]@meta.data,time.colors = time.colors,size.point = 2.5,box.lwd = 1.5)+theme(aspect.ratio = 1.5))
  dev.off()
  #direct.select.si.list[[seu]] <- seu.vitro.nor.list[[seu]]@meta.data
}
saveRDS(seu.vitro.nor.list,'seu.vitro.nor.list.rds')


direct.projct.pval.list <- list()
for(seu in names(direct.select.si.list)){
  age.list <- names(table(direct.select.si.list[[seu]]$Type))
  direct.projct.pval.list[[seu]] <- matrix("-",
                                           ncol = length(age.list),
                                           nrow = length(age.list)
  )
  colnames(direct.projct.pval.list[[seu]]) <- age.list
  rownames(direct.projct.pval.list[[seu]]) <- age.list
}

for(seu in names(direct.select.si.list)){
  age.list <- names(table(direct.select.si.list[[seu]]$Type))
  rm <- age.list [1]
  for(a1 in age.list){
    rm <- unique(c(rm,a1))
    age.list2 <- age.list[!age.list %in% rm]
    for(a2 in age.list2){
      tmp.p <- wilcox.test(direct.select.si.list[[seu]][direct.select.si.list[[seu]]$Type==a1,'pseudotime'],
                           direct.select.si.list[[seu]][direct.select.si.list[[seu]]$Type==a2,'pseudotime'],
                           paired = F
      )
      direct.projct.pval.list[[seu]][a1,a2] <- tmp.p$p.value
    }
  }
}

##########
for(seu in names(direct.select.si.list)){
  MyWriteTable(direct.projct.pval.list[[seu]],row.names = T,paste('direct/20221120/',seu,'.231002.direct.cor.preg.c1.pval.tab',sep = ''))
  #MyWriteTable(direct.select.si.list[[seu]],row.names = T,paste(seu,'.si.tab',sep = ''))
  MyWriteTable(table(direct.select.si.list[[seu]]$Type),row.names = T,paste(seu,'.type.count.tab',sep = ''))
}

##########
# set.seed(5)
# graph.list <- list()
# seu <- names(seu.vitro.nor.list)[2]
# graph.list[[seu]] = graph.adjacency(as.matrix(seu.vitro.nor.list[[seu]]@graphs$RNA_snn),mode = "undirected",weighted = T)
# layout.3d = layout_with_fr(graph.list[[seu]],dim = 3)
# 
# pan.data.par3d <- par3d()
# open3d(zoom = pan.data.par3d$zoom, userMatrix = pan.data.par3d$userMatrix, windowRect = pan.data.par3d$windowRect)
# 
# size.point <- 5
# plot3d(layout.3d,
#        # aspect = 'iso',
#        xlab='FDL1',ylab='FDL2',zlab='FDL3',
#        size = size.point,
#        col = time.colors[seu.vitro.nor.list[[seu]]@meta.data$Type])
# 
# ##########3d#########
library(igraph)
library(plotly)
library(rgl)
p.3d.pca <- plot_ly(x = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,1],
                    y = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,2],
                    z = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,3],
                    type = "scatter3d",
                    mode = "markers",
                    color = seu.vitro.nor.list[[seu]]@meta.data$Type,
                    colors = time.colors,
                    text = seu.vitro.nor.list[[seu]]@meta.data$SampleName,
                    size =30,
                    sizes = c(30,200))

htmlwidgets::saveWidget(as.widget(p.3d.pca), "direct/20221120/sm2.A485.time.3d.pca.html")

plot3d(x = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,1],
       y = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,2],
       z = Embeddings(seu.vitro.nor.list[[seu]],'pca')[,3],
       xlab='PC1',ylab='PC2',zlab='PC3',
       size = 8,
       col=MyName2Col(seu.vitro.nor.list[[seu]]$Type,
                      time.colors))
##########
rm(seu.vitro.nor.list)
save.image('direct/20221120/direct.20221120.RData')
save.image('direct/20221120/direct.20230606.RData')
save.image('direct/20221120/direct.20230704.RData')
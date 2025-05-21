dir.create('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20221203')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20221203/')
load('P7NLHFD.RData')
load('HFD.RData')
library(future)
plan("multisession", workers = 10)
plan()

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(Vennerable)

rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)

ko.time.col <- c('gray80',time.colors[c(1,2)],gene.tree.colors[3],time.colors[c(20,
                                        5,20,8)],gene.tree.colors[6],time.colors[13]
                 ,"#8ECFC9",'#FA7F6F',ref.time.colors['P7NL'],ref.time.colors['P7L']
                 )

time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")

ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)

ko.time.col <- c('gray80',ref.time.colors2['Virgin'],time.col[c(2)],time.col[3],
                 time.col[c(21,5,23)],ref.time.colors2['G14.5'],time.col[6],time.col[13],"#8ECFC9",'#FA7F6F',ref.time.colors2['P7NL'],
                 ref.time.colors2['P7L']
)

names(ko.time.col) <- c('/','Virgin','Virgin_HFD','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD',
                        'G0Acss2HFDWT','G0Acss2HFDHMKO',
                        'G14.5','G14.5CDWT','G14.5_HFD','G14.5HFDWT','G14.5_Acss2HMKO_HFD','P7NL','P7NL_HFD'
                        
)

MyPlotColor(ko.time.col,length(ko.time.col))

pre.ambigous.sym <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/ref/CD_HFD/20220913/pre.ambigous.sym.rds')
TAM.sym <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/KOdirect/220801/TAM.sym.rds')
exo.sym <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/ref/CD_HFD/220801/beta.exo.sym.rds')


seu.beta <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/cluster/seu.beta.ges.post.rds')
seu.HFD.beta <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20220913/G0G14.5HMKO/addG18.5/final/seu.WTHFD.rds')
seu.beta$KOtype <- '/'
seu.beta$KOtype2 <- seu.beta$Type

seu.ko.list <- list()
seu.ko.list[['ref_HMKOHFD']] <- merge(subset(seu.beta,cells = colnames(seu.beta)[seu.beta$Type %in% c('Virgin','G14.5')]),
                                      subset(seu.HFD.beta,cells = colnames(seu.HFD.beta)[!seu.HFD.beta$Type %in% c('Virgin','G14.5')])
                                      )
table(seu.ko.list[['ref_HMKOHFD']]$Type)  
seu.ko.list[['ref_HMKOHFD']]$Type <- factor(seu.ko.list[['ref_HMKOHFD']]$Type,levels = c(c('Virgin','Virgin_HFD','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD','G14.5','G14.5CDWT','G14.5_HFD','G14.5HFDWT','G14.5_Acss2HMKO_HFD')))
seu.ko.list[['ref_HMKOHFD']]$KOtype2 <- factor(seu.ko.list[['ref_HMKOHFD']]$KOtype2,levels = c('/',c('Virgin','Virgin_HFD','G14.5','G14.5_HFD')))
seu.ko.list[['ref_HMKOHFD']]$KOtype <- factor(seu.ko.list[['ref_HMKOHFD']]$KOtype,levels = c('/','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD','G14.5CDWT','G14.5HFDWT','G14.5_Acss2HMKO_HFD'))
i <- 'ref_HMKOHFD'

outlier.sn <- colnames(seu.ko.list[[i]])[(exp(seu.ko.list[[i]]@assays$RNA@data)-1)['Acss2',]>=20]
seu.ko.list[[i]] <- subset(seu.ko.list[[i]],cells = setdiff(colnames(seu.ko.list[[i]]),outlier.sn))

gene.input <- rownames(MyGeneExp(seu.ko.list$ref_HMKOHFD@assays$RNA@data,log(2),10))
gene.input <- setdiff(gene.input,c(TAM.sym,exo.sym,pre.ambigous.sym))
length(gene.input)#129009

pre.ambigous.sym <- unique(c(pre.ambigous.sym,TAM.sym,exo.sym))


i <- 'B6'
seu.ko.list[[i]] <- subset(seu.ko.list[['ref_HMKOHFD']],cells = colnames(seu.ko.list[['ref_HMKOHFD']])[seu.ko.list[['ref_HMKOHFD']]$KOtype=='/'])
seu.ko.list[[i]]$Type <- factor(seu.ko.list[[i]]$Type,levels = c('Virgin','Virgin_HFD','G14.5','G14.5_HFD'))

seu.ko.list <- readRDS('seu.ko.list.rds')
table(seu.ko.list[[i]]$Type)

seu.ko.list <- readRDS('seu.ko.list.rds')
i <- 'ref_HMKOHFD'

seu.ko.list <- readRDS('seu.ko.list.rds')


seu.P7NL <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/20220726/cluster/seu.ko.HFD.list.GLH.rds')
seu.P7NL <- subset(seu.P7NL$HFD,cells = colnames(seu.P7NL$HFD)[seu.P7NL$HFD$Type=='P7NL_HFD'&seu.P7NL$HFD$knn.beta=='Glut2H'])
table(seu.P7NL$Type,seu.P7NL$Rep)
seu.beta.all <- readRDS('../../../beta/20221203/Glut2H/cluster/seu.beta.ges.post.rds')
seu.beta.all <- subset(seu.beta.all,cells = colnames(seu.beta.all)[seu.beta.all$Type=='P7NL'])

pdf('P7NLHFD.genenumber.pdf')
hist(seu.P7NL$Genecount_all,breaks = 30,ylim = c(0,40))
dev.off()

seu.ko.list[['P7NLHFD']] <- merge(seu.P7NL,c(seu.beta.all,seu.ko.list$B6))
table(seu.ko.list[['P7NLHFD']]$Type)

seu.ko.list[['P7NLHFD']]$Type <- factor(seu.ko.list[['P7NLHFD']]$Type,levels = c('Virgin','Virgin_HFD','G14.5','G14.5_HFD','P7NL','P7NL_HFD'))
i <- 'P7NLHFD'
rm(seu.P7NL)
gc()
save.image('P7NLHFD.RData')
seu.ko.list <- readRDS('seu.ko.list.rds')
names(seu.ko.list)

seu.ko.list[['ref_HMKOHFD_P7NLHFD']] <- merge(seu.ko.list$ref_HMKOHFD,subset(seu.ko.list$P7NLHFD,cells = colnames(seu.ko.list$P7NLHFD)[seu.ko.list$P7NLHFD$Type %in% c('P7NL','P7NL_HFD')]))
table(seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type)

seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type <- as.character(seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type)
seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type <- factor(seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type,levels = c(c('Virgin','Virgin_HFD','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD','G14.5','G14.5CDWT','G14.5_HFD','G14.5HFDWT','G14.5_Acss2HMKO_HFD','P7NL','P7NL_HFD')))
seu.ko.list[['ref_HMKOHFD_P7NLHFD']]@meta.data[seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type=='P7NL','KOtype2'] <- 'P7NL'
seu.ko.list[['ref_HMKOHFD_P7NLHFD']]@meta.data[seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type=='P7NL_HFD','KOtype2'] <- 'P7NL_HFD'

seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$KOtype2 <- factor(seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$KOtype2,levels = c('/',c('Virgin','Virgin_HFD','G14.5','G14.5_HFD','P7NL','P7NL_HFD')))

seu.ko.list[['ref_HMKOHFD_P7NLHFD']]@meta.data[seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$Type %in% c('P7NL','P7NL_HFD'),'KOtype'] <- '/'

seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$KOtype <- factor(seu.ko.list[['ref_HMKOHFD_P7NLHFD']]$KOtype,levels = c('/','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD','G14.5CDWT','G14.5HFDWT','G14.5_Acss2HMKO_HFD'))

i <- 'ref_HMKOHFD_P7NLHFD'
#########
seu.ko.list[[i]] <- FindVariableFeatures(seu.ko.list[[i]], selection.method = "vst", nfeatures = 2000)
VariableFeatures(seu.ko.list[[i]]) <- setdiff(VariableFeatures(seu.ko.list[[i]]),c(pre.ambigous.sym))
c1Color <- MyName2Col(seu.ko.list[[i]]@meta.data[,"Type"],
                      time.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.ko.list[[i]]@meta.data[,"Rep"],
                      colors.group)

c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)

png(paste(i,'.Vstvar.2000.geneover6000.raw.png',sep = ''),2000,3000)
assign(paste(i,'var2000.row.tree',sep = '.'),
       MyHeatmap(as.matrix(exp(seu.ko.list[[i]]@assays$RNA@data))[VariableFeatures(seu.ko.list[[i]]),colnames(seu.ko.list[[i]])],
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
                 graph = T))
dev.off()

##########
VariableFeatures(seu.ko.list[[i]]) <- setdiff(VariableFeatures(seu.ko.list[[i]]),labels(as.dendrogram(get(paste(i,'var2000.row.tree',sep = '.')))[[1]]))
seu.ko.list[[i]] <- ScaleData(seu.ko.list[[i]],features = VariableFeatures(seu.ko.list[[i]]))
seu.ko.list[[i]] <- RunPCA(seu.ko.list[[i]],features = VariableFeatures(seu.ko.list[[i]]))
#seu.ko.list[[i]] <- RunPCA(seu.ko.list[[i]],features = rownames(get(paste(i,'beta.var.co15',sep = ''))))

pc.use <- 1:5;#pc5 ref KO P7NLHFD:8
seu.ko.list[[i]] <- FindNeighbors(seu.ko.list[[i]], dims = pc.use)
seu.ko.list[[i]] <- FindClusters(seu.ko.list[[i]], resolution = 0.3)
seu.ko.list[[i]] <- RunUMAP(seu.ko.list[[i]],dims = pc.use)
seu.ko.list[[i]] <- RunTSNE(seu.ko.list[[i]],dims = pc.use)

seu.ko.list[[i]] <- SetIdent(seu.ko.list[[i]],value = seu.ko.list[[i]]$Type)
DimPlot(seu.ko.list[[i]],
        reduction = "umap",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 3,
        label = F)
DimPlot(seu.ko.list[[i]],
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 3,
        label = F)
seu.ko.list[[i]] <- SetIdent(seu.ko.list[[i]],value = seu.ko.list[[i]]$Type)
DimPlot(seu.ko.list[[i]],
        reduction = "pca",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = F)

pdf('B6.pc13.pdf',15,12)
p.pca <- MySeuratDR2Gg2(seu.ko.list[[i]],seu.ko.list[[i]]@meta.data,reduction.use = 'pca',
                        reduction.key = 'PC',x.dim = 1,y.dim = 2,
                        estimate.variation.explain.percentage = T)



plot(p.pca+
       scale_color_manual(values = c(ko.time.col)) +
       geom_point(aes(x = -x.pos,
                      y = y.pos,
                      col = Type#,
                      # shape = State
       ),
       size = 5
       )+theme(aspect.ratio = 1))
# 
dev.off()

pdf('G:/lab/Article/heshuang/BYLW/sm3/HFD/B6.only.umap12.pdf',15,12)
p.pca <- MySeuratDR2Gg2(seu.ko.list[[i]],seu.ko.list[[i]]@meta.data,reduction.use = 'umap',
                        reduction.key = 'UMAP',x.dim = 1,y.dim = 2,
                        estimate.variation.explain.percentage = F)

p.pca$data_ <- p.pca$data
p.pca$data_$Type_ <- p.pca$data_$KOtype2
p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(p.pca$data_$KOtype2)))
p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ko.time.col[names(table(seu.ko.list[[i]]$KOtype2))]) +
       geom_point(aes(x = x.pos,
                      y = y.pos,
                      col = KOtype2#,
                      # shape = State
       ),
       size = 6
       )+theme(aspect.ratio = 1))
dev.off()

pdf('P7NLHFD.pc5.only.pca.pdf',15,12)
p.pca <- MySeuratDR2Gg2(seu.ko.list[[i]],seu.ko.list[[i]]@meta.data,reduction.use = 'pca',
                        reduction.key = 'PC',x.dim = 1,y.dim = 2,
                        estimate.variation.explain.percentage = T)

# p.pca$data_ <- p.pca$data
# p.pca$data_$Type_ <- p.pca$data_$KOtype2
# p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(p.pca$data_$KOtype2)))
# p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ko.time.col[names(table(seu.ko.list[[i]]$Type))]) +
       geom_point(aes(x = x.pos,
                      y = y.pos,
                      col = Type#,
                      # shape = State
       ),
       size = 6
       )+theme(aspect.ratio = 1))
p.pca <- MySeuratDR2Gg2(seu.ko.list[[i]],seu.ko.list[[i]]@meta.data,reduction.use = 'pca',
                        reduction.key = 'PC',x.dim = 1,y.dim = 3,
                        estimate.variation.explain.percentage = T)

plot(p.pca+
       scale_color_manual(values = ko.time.col[names(table(seu.ko.list[[i]]$Type))]) +
       geom_point(aes(x = x.pos,
                      y = y.pos,
                      col = Type#,
                      # shape = State
       ),
       size = 6
       )+theme(aspect.ratio = 1))
dev.off()

seu.ko.list <- readRDS('seu.ko.list.rds')
pdf('G:/lab/Article/heshuang/BYLW/sm3/HFD/WT.umap12.pdf',18,12)
p.pca <- MySeuratDR2Gg2(seu.ko.list[[i]],seu.ko.list[[i]]@meta.data,reduction.use = 'umap',
                        reduction.key = 'UMAP',x.dim = 1,y.dim = 2,
                        estimate.variation.explain.percentage = F)

p.pca$data_ <- p.pca$data
p.pca$data_$Type_ <- p.pca$data_$KOtype
p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(p.pca$data_$KOtype)))
p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ko.time.col[names(table(seu.ko.list[[i]]$KOtype))]) +
       geom_point(aes(x = x.pos,
                      y = y.pos,
                      col = KOtype#,
                      # shape = State
       ),
       size = 5
       )+theme(aspect.ratio = 1))


dev.off()





seu <- 'P7NLHFD'
seu.ko.list[[seu]]$pseudotime <- -Embeddings(seu.ko.list[[seu]],'pca')[,1]
select.si.list2 <- list()
for(seu in names(seu.ko.list)[2]){
  # seu.ko.list[[seu]]$pseudotime <- -seu.ko.list[[seu]]$pseudotime
  select.si.list2[[seu]] <- seu.ko.list[[seu]]@meta.data
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/HFD/',seu,'.pc1.pseudobox.pdf',sep = ''),10,6)
  print(MyPseudotimebox(select.si.list2[[seu]],time.colors = ko.time.col[names(table(select.si.list2[[seu]]$Type))],size.point = 2,box.lwd = 1.5)+theme(aspect.ratio = 1.5))
  dev.off()
}
#########pval#######
projct.pval.list2 <- list()
for(seu in names(select.si.list2)){
  
  age.list <- names(table(select.si.list2[[seu]]$Type))
  projct.pval.list2[[seu]] <- matrix("-",
                                     ncol = length(age.list),
                                     nrow = length(age.list)
  )
  colnames(projct.pval.list2[[seu]]) <- age.list
  rownames(projct.pval.list2[[seu]]) <- age.list
}

for(seu in names(select.si.list2)){
  age.list <- names(table(select.si.list2[[seu]]$Type))
  rm <- age.list [1]
  for(a1 in age.list){
    rm <- unique(c(rm,a1))
    age.list2 <- age.list[!age.list %in% rm]
    for(a2 in age.list2){
      tmp.p <- wilcox.test(select.si.list2[[seu]][select.si.list2[[seu]]$Type==a1,'pseudotime'],
                           select.si.list2[[seu]][select.si.list2[[seu]]$Type==a2,'pseudotime'],
                           paired = F
      )
      projct.pval.list2[[seu]][a1,a2] <- tmp.p$p.value
    }
  }
}

for(seu in names(select.si.list2)){
  MyWriteTable(projct.pval.list2[[seu]],row.names = T,paste(seu,'.B6.pc1.pval.tab',sep = ''))
  #MyWriteTable(select.si.list[[seu]],row.names = T,paste(seu,'si.tab',sep = ''))
  #MyWriteTable(table(select.si.list[[seu]]$Type),row.names = T,paste(seu,'type.count.tab',sep = ''))
}
#############
seu.ko.list <- readRDS('seu.ko.list.rds')
library(rgl)
library(plotly)
#par3d()
P7NLHFD.res.3dPCA <- readRDS('P7NLHFD.res.3dPCA.rds')
open3d(zoom = P7NLHFD.res.3dPCA$zoom, userMatrix = P7NLHFD.res.3dPCA$userMatrix, windowRect = P7NLHFD.res.3dPCA$windowRect)
plot3d(x=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,1],
       y=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,2],
       z=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,3],
       size=12,
       #sizes = c(300,200),
       xlab = "", ylab = "", zlab = "",
       box=F,lwd = 1,
       axes = F,
       col=MyName2Col(seu.ko.list[[i]]$Type,
                      c(ko.time.col[names(table(seu.ko.list[[i]]$Type))])))
box3d(lwd=2)


rgl.snapshot("G:/lab/Article/heshuang/BYLW/sm3/HFD/P7NLHFD.3d.PCA.box.png")

P7NLHFD.res <- par3d()
saveRDS(P7NLHFD.res,'P7NLHFD.res.3dPCA.rds')
saveRDS(seu.ko.list,'seu.ko.list.rds')

plot3d(x=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,1],
       y=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,2],
       z=seu.ko.list[[i]]@reductions$pca@cell.embeddings[,3],
       size=12,
       #sizes = c(300,200),
       xlab = "", ylab = "", zlab = "",
       box=F,lwd = 1,
       axes = F,
       col=MyName2Col(seu.ko.list[[i]]$KOtype,
                      c(ko.time.col[names(table(seu.ko.list[[i]]$KOtype))])))

box3d(lwd=2)
title3d(xlab = paste("PC1 (",'3.9',"%)",sep = ""),
        ylab = paste("PC2 (",'1.4',"%)",sep = ""),
        zlab = paste("PC3 (",'0.9',"%)",sep = ""))

title3d(xlab = paste("PC1 (",'4.1',"%)",sep = ""),
        ylab = paste("PC2 (",'2',"%)",sep = ""),
        zlab = paste("PC3 (",'1',"%)",sep = ""))

###########
seu.ko.list[[seu]]$PC_1 <- seu.ko.list[[seu]]@reductions$pca@cell.embeddings[,1]
seu.ko.list[[seu]]$PC_2 <- seu.ko.list[[seu]]@reductions$pca@cell.embeddings[,2]
seu.ko.list[[seu]]$PC_3 <- seu.ko.list[[seu]]@reductions$pca@cell.embeddings[,3]
pca.df <- as.data.frame(seu.ko.list[[seu]]@meta.data)

Type3dList <- split(pca.df, pca.df$KOtype)
size.list <- list(5,12,12,12,12,12)
ko.time.col2 <- ko.time.col[names(table(pca.df$KOtype))]
open3d()
with(pca.df, plot3d(PC_1, PC_2, PC_3, size=0,aspect = 1,pch=19,xlab = '',ylab = '',zlab = '',#axes = F,
                    col=ko.time.col2[as.factor(KOtype)]))


for(i in seq_along(Type3dList)) {
  print(i)
  with(Type3dList[[i]], points3d(PC_1, PC_2, PC_3,size=size.list[[i]],col=ko.time.col[as.factor(KOtype)]))
}


rgl.snapshot("G:/lab/Article/heshuang/BYLW/sm3/HFD/B6.G14.5WTKO.pca.3d.png")

##########vioplot######
seu.ko.list <- readRDS('seu.ko.list.rds')
i <- 'B6'
seu.tmp <- seu.ko.list[[i]]
seu.tmp <- subset(seu.ko.list[[i]],cells = colnames(seu.ko.list[[i]])[seu.ko.list[[i]]$KOtype =='/'])
seu.tmp$Type <- factor(seu.tmp$Type,levels = c('Virgin','Virgin_HFD','G14.5','G14.5_HFD'))  
p.list <- list()
p.list[['G14.5CD_HFD']] <- wilcox.test(seu.tmp@assays$RNA@data['Acss2',colnames(seu.tmp)[seu.tmp$Type=='G14.5']],
                                       seu.tmp@assays$RNA@data['Acss2',colnames(seu.tmp)[seu.tmp$Type=='G14.5_HFD']],
                                       paired = F
)$p.value
p.list[['G0CD_HFD']] <- wilcox.test(seu.tmp@assays$RNA@data['Acss2',colnames(seu.tmp)[seu.tmp$Type=='Virgin']],
                                    seu.tmp@assays$RNA@data['Acss2',colnames(seu.tmp)[seu.tmp$Type=='Virgin_HFD']],
                                    paired = F
)$p.value
pdf('G:/lab/Article/heshuang/BYLW/sm3/HFD/Fig6e.B6vio.Acss2.tp0.1m.pdf',6,8)
Myseuvioplot(seu.tmp,c('Acss2'),type = 'Type',type.colors = ko.time.col[names(table(seu.tmp$Type))],log = F,ext=3,med.lwd=8)
MyText(paste(
             'G0CD_HFD: pval',p.list[['G0CD_HFD']],'\n',
             'G14.5CD_HFD: pval',p.list[['G14.5CD_HFD']],sep=''
             ),text.cex = 0.5)
dev.off()

rm(seu.tmp)
#########P7NLHFD heatmap######
G14.5HFD.gene.tab <- MyReadDelim('DEG/gene_inf/padj0.01.fc1.5.G14.5_G14.5_HFD.si.tab')
pseu.df <- do.call(cbind,by(t(expm1(seu.ko.list$P7NLHFD@assays$RNA@data[G14.5HFD.gene.tab$gene[G14.5HFD.gene.tab$cluster=='G14.5_HFD'],])),seu.ko.list$P7NLHFD$Type,colMeans))
boxplot(pseu.df,outline=F)
VlnPlot(seu.ko.list$P7NLHFD,'nFeature_RNA')
MyReadDelim('../')


c1Color <- MyName2Col(seu.ko.list[[i]]@meta.data[,"Type"],
                      time.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.ko.list[[i]]@meta.data[,"Rep"],
                      colors.group)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)

png(paste(i,'.Vstvar.2000.geneover6000.raw.png',sep = ''),2000,3000)
       MyHeatmap(as.matrix(pseu.df),
                 type = "log.row.relat",
                 hc.c.data.type = "log.row.relat",
                 hc.r.data.type = "log.row.relat",
                 c.cov.method = "s",
                 r.cov.method = "s",
                 c.hc.method = "ward.D2",
                 r.hc.method = "ward.D2",
                # ColSideColors = cColor,
                 ColSideColorsSize = 2,
                 Colv = 'none',
                Rowv = 'do',
                dendrogram='row',
                # return.tree = "row",
                 graph = T)
dev.off()
rm(seu.ko.list)
save.image('P7NLHFD.RData')
#######
rm(seu.ko.list,seu.tmp)
saveRDS(seu.ko.list,'seu.ko.list.rds')
   rm(seu.ko.list)
save.image('G14.5Acss2HMKOHFD.RData')

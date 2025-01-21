setwd('G:/project/pregnant_mouse/10x/10x_v3')
dir.create('Ctl_G14.5_1st/')
setwd('Ctl_G14.5_1st/')
load('pre.cluster.RData')
load('beta.10x.DEG.RData')
source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
# library(dplyr)
# library(rgl)
# library(scran)
colors.group <- c("#E08698",
                  "#999999",
                  "#A9AF62",
                  "#41BAB0",
                  "#8D9AC6",
                  "#E579E5",
                  # "black",
                  "#ff7f00",
                  "#6a3d9a",
                  "#e31a1c",
                  "#b15928",
                  "#33a02c",
                  "#a6cee3",
                  "#67001f",
                  "#969696"
                  
)
shape.group <- c(16,8,2)

colors.exp <- colorRampPalette(c("#0000FF","white","#FF0000"),
                               space="Lab")(50)
time.colors <- c("#4DAF4A",
                 "#752a78",
                 #"black",
                 "#F01414",
                 "#5F9EA0",#P21-2-B-con
                 "#FF7F00",
                 "#7570B3", 
                 "#452cff",#G5.5
                 "#A65628", 
                 "#F781BF",
                 "#999999",
                 #"#ffe478",#P2-B
                 "#cab2d6",
                 "#ba984d",
                 "#1f78b4",
                 "#ff4a21",#P7-14-21,P21-B
                 "#752a78",
                 "#a6cee3",
                 "#53751c",
                 "#C71585",
                 "#d15126",
                 #"#6bff9c",#G14.5-2nd
                 "#5ce6ba",
                 "#ff8e9d"
)
MyPlotColor(time.colors,length(time.colors))

gene.tree.colors <- brewer.pal(8,"Set2")


preg.colors <- c("#4DAF4A",
                 "#F781BF"
)

names(preg.colors) <- c("Ctrl",
                        "G14.5_1st")
MyPlotColor(preg.colors,length(preg.colors))


colors.exp <- colorRampPalette(c("#0000FF","white","#FF0000"),
                               space="Lab")(50)
##########load data#########
Ctrl.seu <- readRDS('../QC/src.duct.Ctrl.rds')
Ctrl.seu.meta.tab <- MyReadDelim('../QC/ctrl.celltype.new.meta.data')
#Ctrl.seu$SampleName <- paste0("Ctrl_20200416",Ctrl.seu$SampleName)
#colnames(Ctrl.seu) <- Ctrl.seu$SampleName
Ctrl.seu <- subset(Ctrl.seu,cells = Ctrl.seu.meta.tab$SampleName[!Ctrl.seu.meta.tab$cluster %in% c( 'Doublet/multi-hormone cell',
                                                                                                   'Low quality cell')])

G14.5_1st.seu <- readRDS('../QC/pre.seu.rds')
G14.5_1st.seu.meta.tab <- MyReadDelim('../QC/G14.5_1st.celltype.new.meta.data')
G14.5_1st.seu <- subset(G14.5_1st.seu,cells = G14.5_1st.seu.meta.tab$SampleName[G14.5_1st.seu.meta.tab$cluster != 'Doublet/multi-hormone cell'])

Ctrl.seu <- RenameCells(Ctrl.seu,new.names = paste0('Ctrl_20200416',colnames(Ctrl.seu)))
G14.5_1st.seu <- RenameCells(G14.5_1st.seu,new.names = paste0('G14.5_1st_20200416',colnames(G14.5_1st.seu)))

pre.seu <- merge(x=Ctrl.seu,
                 y=G14.5_1st.seu
                 )
rm(Ctrl.seu)
rm(G14.5_1st.seu)


pre.seu <- readRDS('Ctl_G14.5_1st/pre.seu.rds')

mean(pre.seu$nCount_RNA)
mean(pre.seu$nFeature_RNA)
table(pre.seu$merge.endo.cluster)
pre.seu <- SetIdent(pre.seu,value = pre.seu$merge.endo.cluster)
VlnPlot(pre.seu,'Slc2a5',split.by = 'merge.endo.cluster')

######gene count#####
pre.seu <- readRDS('pre.seu.rds')
dir.create('G:/lab/Article/heshuang')
dir.create('G:/lab/Article/heshuang/BYLW')
dir.create('G:/lab/Article/heshuang/10x')
pdf('G:/lab/Article/heshuang/BYLW/10x/genecount.2.pdf',4,5)
boxplot(pre.seu$nFeature_RNA,lwd=2,main='gene number')
box(lwd=3)
boxplot(pre.seu$nCount_RNA,lwd=2,main='read count')
box(lwd=3)
dev.off()
########cluster######
dir.create('all')
pre.seu <- 
  FindVariableFeatures(pre.seu, selection.method = "vst", nfeatures = 2000)

all.var.co <- rownames(MyCo(as.matrix(pre.seu@assays$RNA@data),
                            var.gene = VariableFeatures(pre.seu),
                            exp.prop.whole.max = 0.9,
                            exp.prop.whole.min = 0.005,
                            # vector.group = samples.inf.qc$GroupNameFig1,
                            # exp.prop.group.min = 0.1,
                            # exp.prop.group.max = 0.5,
                            cor.method = "rho",
                            cor.cutoff = 0.15,
                            partner.cutoff = 5,
                            refine.cor = T))
length(all.var.co)#1264
# c1Color <- MyName2Col(pre.seu@meta.data[,"RNA_snn_res.0.5"],
#                       brewer.pal(9,"Set1"))
# c1Color <- as.matrix(c1Color)

png('all/varcor0.15.png',2000,3000)
var.row.tree <-
  MyHeatmap(as.matrix(exp(pre.seu@assays$RNA@data)[all.var.co,]),
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
           # ColSideColors = c1Color,
            return.tree = "row",
            graph = T)
dev.off()

var.row.tree <- as.dendrogram(var.row.tree)
all.cc.sym <- labels(var.row.tree[[2]][[2]][[2]][[1]][[2]][[1]])


pre.seu <- 
  ScaleData(pre.seu)
pre.seu <- RunPCA(pre.seu, features = setdiff(all.var.co,
                                              all.cc.sym))

MyWriteTable(setdiff(all.var.co,
                     all.cc.sym),
             'pre.input.var.tab')
######eound2 exclude cc gene#########
pdf('all/beta.DimHeatmap.rmdoubelt.rmcc.pdf',
    10,40)
DimHeatmap(pre.seu, dims = 1:30, cells = 200, balanced = TRUE)
dev.off()
pre.seu <- readRDS('pre.seu.rds')

pc.use <- 1:12
pre.seu <- FindNeighbors(pre.seu, dims = pc.use)
pre.seu <- FindClusters(pre.seu, resolution = 2)
pre.seu <- FindClusters(pre.seu, resolution = 1)
pre.seu <- FindClusters(pre.seu, resolution = 0.5)
pre.seu <- FindClusters(pre.seu, resolution = 0.2)
pre.seu <- FindClusters(pre.seu, resolution = 0.1)


time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")

ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)


pre.seu <- RunUMAP(pre.seu,dims = pc.use)

pre.seu <- RunTSNE(pre.seu,
                   dims = pc.use,
                   perplexity= round((30+ncol(pre.seu)/100)),
                   check_duplicates = F)
pre.seu <-
  SetIdent(pre.seu, value = pre.seu$RNA_snn_res.2)


pdf('all/rmcc.beta.pure.tsne.pc12.res2.pdf',
    10,8)
DimPlot(pre.seu,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 1,
        label = T)
dev.off()

p.tsne <- MySeuratv3TSNE10x2Gg(pre.seu,
                               pre.seu@meta.data
                               #pp.meta.data
)

p.tsne <- MySeuratv3UMAP10x2Gg(pre.seu,
                               pre.seu@meta.data
                               #pp.meta.data
)


p.age <- p.tsne+
  #scale_y_reverse() +
  scale_color_manual(values = time.colors) +
  #scale_shape_manual(values = shape.group) +
  guides(colour = guide_legend(title = "",
                               order = 1)) +
  theme(axis.text = element_text(size = 20, colour = "black")) +
  theme(axis.title.x = element_text(size = 40, colour = "black")) +
  theme(axis.title.y = element_text(size = 40, colour = "black")) +
  theme(legend.text = element_text(size = 45, colour = "black")) +
  theme(legend.title = element_text(size = 45, colour = "black")) +
  theme(panel.border = element_rect(size = 4,
                                    colour = "black")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(title = "",
                               keywidth = 3,
                               keyheight = 3,
                               override.aes = list(size = 12),
                               order = 1,
                               ncol = 1)) +
  guides(shape = guide_legend(title = "",
                              keywidth = 3,
                              keyheight = 3,
                              override.aes = list(size = 12),
                              order = 2,
                              ncol = 1))

pdf("G:/lab/Article/heshuang/BYLW/10x/umap.cellType.pc12old.pdf",
    18,
    10)
p.plot <- p.age +
  scale_color_manual(values = time.colors) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =merge.endo.cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot+theme(aspect.ratio = 1))

dev.off()

pdf("G:/lab/Article/heshuang/BYLW/10x/umap.Time.pc12.pdf",
    15,
    10)
p.plot <- p.age +
  scale_color_manual(values = c("#b6d4a8","#b89fc1")) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =Time#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot+theme(aspect.ratio = 1))
dev.off()

########marker########
marker.sym <- c("Spi1",#°×Ï¸°û
                "Procr",
                "Fcgr1",
                "Neurod1",#endocrine
                "Chga",
                "Chgb",
                "Prox1",
                "Mki67",
                "Col1a1",#mesenchymal
                "Col5a1",#mesenchymal
                "Col3a1",#mesenchymal
                "Mgp",
                "Esm1",
                "Vim",
                "Cdh11",
                #"S100a4",
                "Itga5",
                "Sdc1",
                "Pecam1",#Ñª¹Ü
                "Gypa",#ºìÏ¸°û
                "Sox17",#¸ÎÍâµ¨¹Ü
                "Sox9",#duct
                "Spp1",#duct
                "Rbpjl",#acinar
                # "Amy2b",#acinar
                "Pnlip",#exo
                "Ctrb1",#exo
                "Prss1",
                "Neurog3",
                "Pax4",
                "Ins1",
                "Gcg",
                "Sst",
                "Ppy",
                "Ghrl",
                "Nkx6-1",
                "Ucn3",
                "Mafb",
                "Mafa",
                "Slc2a2",
                "Meg3",
                "Cacna1a",
                "Xist",
                "Kcnq1ot1",
                "Arx",
                "Hhex",
                "Kdr",
                "Ascl1",
                "Ptprc"
                
                
)
length(marker.sym)#
marker.sym <- marker.sym[marker.sym %in% rownames(pre.seu)]
length(marker.sym)#
marker.sym <- unique(marker.sym)
length(marker.sym)#45
total.count <- length(marker.sym)
run.count <- 1

marker.sym <- c('Neurod1','Pecam1','Col3a1','Fcgr1','Spp1')

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
pre.seu$merge.endo.cluster <- factor(pre.seu$merge.endo.cluster,levels = )
pre.seu <- SetIdent(pre.seu,value = pre.seu$merge.endo.cluster)

pdf('G:/lab/Article/heshuang/BYLW/10x/marker.vioplot.pdf',8,8)
StackedVlnPlot(pre.seu,features = marker.sym,cols=time.colors,pt.size = 1)+theme(aspect.ratio = 1)
dev.off()


pdf(paste('all/pc12.marker',".color.exp.pdf",sep = ""),6,7)

for (gene in c('percent.mt','nFeature_RNA',marker.sym)) {
  cat(paste(run.count, "/",total.count,"\n"))
  print(FeaturePlot(pre.seu,gene,reduction = 'tsne')+scale_color_gradientn(colours = colors.exp)+
          theme(plot.title = element_text(size = rel(3.5))) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.title = element_blank()) +
          theme(axis.text  = element_blank()) +
          theme(axis.ticks = element_blank()) +
          theme(legend.position="bottom")+
          theme(panel.border = element_rect(size = 4,
                                            colour = "black"))+
          guides(colour = guide_colorbar(title = "ln(TP10K+1)",
                                         title.position = "top",
                                         barwidth = 27,
                                         title.hjust = 0.5,
                                         title.theme = element_text(angle = 0,
                                                                    size = 20),
                                         label.theme = element_text(angle = 0,
                                                                    size = 20),
                                         ticks = T)) )
  run.count = run.count + 1
}
dev.off()
#######cell type######

pre.seu <- readRDS('pre.seu.rds')
pre.seu@meta.data$merge.cluster=as.character(pre.seu@meta.data$RNA_snn_res.2)
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(3,14,9,13,22,15,12,10,7,21,6,8,11,6,16,4),]$merge.cluster="beta"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(0,1),]$merge.cluster="alpha"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(18),]$merge.cluster="pp"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(2),]$merge.cluster="delta"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(17,5,25),]$merge.cluster="Endothelial cell"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(19),]$merge.cluster="Macrophage cell"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(23),]$merge.cluster="Myeloid cell"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(20,24,26),]$merge.cluster="Mesenchymal cell"
pre.seu@meta.data[pre.seu@meta.data$merge.cluster%in%c(27),]$merge.cluster="Duct cell"

pre.seu@meta.data$merge.cluster <- factor(pre.seu@meta.data$merge.cluster,
                                          levels = c('beta',
                                                     'alpha',
                                                     'delta',
                                                     'pp',
                                                     'Endothelial cell',
                                                     'Macrophage cell',
                                                     "Myeloid cell",
                                                     'Mesenchymal cell',
                                                     'Duct cell'))

table(pre.seu@meta.data$merge.cluster)

pre.seu@meta.data$merge.endo.cluster <- as.character(pre.seu$merge.cluster)
pre.seu@meta.data[pre.seu$merge.endo.cluster %in% c('beta',
                                                    'alpha',
                                                    'delta',
                                                    'pp'),'merge.endo.cluster'] <- 'Endocrine cell'


pre.seu@meta.data$merge.endo.cluster <- factor(pre.seu@meta.data$merge.endo.cluster,
                                          levels = c('Endocrine cell',
                                                     'Endothelial cell',
                                                     'Mesenchymal cell',
                                                     'Immune cell',
                                                     'Duct cell'))
table(pre.seu@meta.data$merge.endo.cluster)
pre.seu <- readRDS('pre.seu.rds')
p.tsne <- MySeuratv3TSNE10x2Gg(pre.seu,
                               pre.seu@meta.data
                               #pp.meta.data
)
p.age <- p.tsne+
  #scale_y_reverse() +
  scale_color_manual(values = time.colors) +
  #scale_shape_manual(values = shape.group) +
  guides(colour = guide_legend(title = "",
                               order = 1)) +
  theme(axis.text = element_text(size = 20, colour = "black")) +
  theme(axis.title.x = element_text(size = 40, colour = "black")) +
  theme(axis.title.y = element_text(size = 40, colour = "black")) +
  theme(legend.text = element_text(size = 45, colour = "black")) +
  theme(legend.title = element_text(size = 45, colour = "black")) +
  theme(panel.border = element_rect(size = 4,
                                    colour = "black")) +
  theme(legend.key = element_blank()) +
  guides(colour = guide_legend(title = "",
                               keywidth = 3,
                               keyheight = 3,
                               override.aes = list(size = 12),
                               order = 1,
                               ncol = 1)) +
  guides(shape = guide_legend(title = "",
                              keywidth = 3,
                              keyheight = 3,
                              override.aes = list(size = 12),
                              order = 2,
                              ncol = 1))

time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")

pdf("G:/lab/Article/heshuang/BYLW/10x/tsne.time.pc12.res2.pdf",
    18,
    10)
# p.plot <- p.age +
#   scale_color_manual(values = colors.group) +
#   geom_point(aes(x = x.pos,
#                  y = y.pos,
#                  col =merge.cluster#,
#                  #shape = State
#   ),
#   size = 1
#   )+theme(aspect.ratio = 1)
# 
# print(p.plot)
p.plot <- p.age +
  scale_color_manual(values = c('cyan3','red3')) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =Time#,
                 #shape = State
  ),
  size = 1
  )+theme(aspect.ratio = 1)

print(p.plot)
dev.off()

celltype.col <- colors.group[c(1:4,8)]
names(celltype.col) <- c('Endocrine cell',
                         'Endothelial cell',
                         'Mesenchymal cell',
                         "Immune cell",
                         'Duct cell')
pdf("FigS1.tsne.merge.endo.cellType.pc12.res2.pdf",
    14.5,
    8)
print(p.age +
  scale_color_manual(values = celltype.col) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =merge.endo.cluster#,
                 #shape = State
  ),
  size = 1
  )+theme(aspect.ratio = 1))

dev.off()
#########marker######
pre.seu <- readRDS('pre.seu.rds')
marker.meanexp.tab <-  t(data.frame('Neurod1' = 1:5,
                                    # 'Chga' = 1:5,
                                    'Pecam1' = 1:5,
                                    'Col3a1' = 1:5,
                                    'Fcgr1' = 1:5,
                                    'Spp1' = 1:5))

colnames(marker.meanexp.tab) <- names(celltype.col)
for(i in rownames(marker.meanexp.tab)){
  for (j in colnames(marker.meanexp.tab)){
    marker.meanexp.tab[i,j] <- mean(as.matrix(pre.seu@assays$RNA@data[i,pre.seu$merge.endo.cluster==j]))
  }
}
marker.meanexp.relative.tab <- t(apply(marker.meanexp.tab, 1, function(x){x/max(x)}))
saveRDS(marker.meanexp.relative.tab,'marker.meanexp.relative.tab.rds')

colors.marker <- colorRampPalette(c("#BEBEBE","white","#E72703"
                                    ),
                                  space="Lab")(50)
pdf('all/Figs1d.marker.celltype.pdf',10,10)
par(oma= c(3,0,0,7))
MyHeatmap(as.matrix(marker.meanexp.relative.tab),
          type = "raw",
          #hc.c.data.type = "log.row.relat",
          #hc.r.data.type = "log.row.relat",
          #c.cov.method = "p",
          ## r.cov.method = "p",
          c.hc.method = "ward.D2",
          #color.palette = colors.heat,
          color.palette = colors.marker,
          #r.hc.method = "ward.D2",
          #ColSideColors = cColor,
          #ColSideColorsSize = 2,
          Colv = 'none',
          Rowv = 'none',
          labRow = rownames(marker.meanexp.relative.tab),
          cexRow = 2#,
          # return.tree = "row"
)
dev.off()
########dot plot#######
marker.meanexp.relative.tab <- readRDS('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/marker.meanexp.relative.tab.rds')

pseudoexp.select.list <- list()

pseudoexp.select.list[['all.exp.ratio']] <- marker.meanexp.relative.tab

for(gene in rownames(pseudoexp.select.list[['all.exp.ratio']])){
  for (group in colnames(pseudoexp.select.list[['all.exp.ratio']])) {
    pseudoexp.select.list[['all.exp.ratio']][gene,group] <- sum(pre.seu@assays$RNA@data[gene,pre.seu$merge.endo.cluster==group]>0)/sum(pre.seu$merge.endo.cluster==group)
  }
}

pseudoexp.select.list[['all.exp.ratio.group']] <- reshape2::melt(pseudoexp.select.list[['all.exp.ratio']][5:1,])
colnames(pseudoexp.select.list[['all.exp.ratio.group']]) <- c('gene','group','exp.ratio')

marker.meanexp.relative.tab2 <- reshape2::melt(marker.meanexp.relative.tab[5:1,])
colnames(marker.meanexp.relative.tab2) <- c('gene2','group2','relat.exp')

pseudoexp.gene.exp.ratio <- cbind(pseudoexp.select.list[['all.exp.ratio.group']],marker.meanexp.relative.tab2)
pseudoexp.gene.exp.ratio$gene <- as.character(pseudoexp.gene.exp.ratio$gene)
pseudoexp.gene.exp.ratio$gene <- factor(pseudoexp.gene.exp.ratio$gene,levels = rev(rownames(marker.meanexp.relative.tab)))

pseudoexp.gene.exp.ratio <- pseudoexp.gene.exp.ratio[order(pseudoexp.gene.exp.ratio$gene,decreasing = F),]

color.dot <- colorRampPalette(c("gray50","#CBABA3","#DDB6AC","#FF0000"),
                              space="Lab")(50)

p <- ggplot(pseudoexp.gene.exp.ratio, 
            aes(group,
                gene))+
  theme_minimal()+
  xlab(NULL)+
  ylab(NULL)
p<- p+scale_color_gradientn(colours = color.dot)+
  geom_point(aes(size=exp.ratio,
                 col=relat.exp))+
  scale_size_continuous(range = c(1, 15))

p <- p + guides(color=guide_colorbar(title ="relative ln(TP0.1M)"), 
                size=guide_legend("ratio"))+
  theme(panel.border = element_rect(size = 2,
                                    colour = "black",
                                    fill=NA))

p.sub.merge <- p+ 
  theme(axis.text = element_text(size = 20, colour = "black")) +
  theme(axis.title.x = element_text(size = 20, colour = "black")) +
  theme(axis.title.y = element_text(size = 20, colour = "black")) +
  theme(legend.text = element_text(size = 20, colour = "black")) +
  theme(legend.title = element_text(size = 20, colour = "black"))

pdf("marekrgene.dotplot.ratio.pdf",
    8,5)
plot(p.sub.merge)
dev.off()


###########
table(pre.seu@meta.data[,c('Time','merge.cluster')])



saveRDS(pre.seu,
        'pre.seu.rds')

rm(pre.seu)
save.image('pre.cluster.RData')
pre.seu <- readRDS('pre.seu.rds')
load('pre.cluster.RData')


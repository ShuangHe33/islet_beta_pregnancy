setwd('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/')
source('G:/pcatest/MyFunction.R')
#######g1g23#########
load('beta/beta.cluster.RData')
load('beta.10x.DEG.RData')

beta.seu <- readRDS('beta/src.beta.rds')
source('G:/pcatest/MyFunction.R')
endoType.col <- c(colors.group[1],"#5F9EA0", "#A65628",colors.group[c(5,6,3)])
names(endoType.col) <- c('beta cell1',
                         'beta cell2',
                         'beta cell3',
                         'alpha cell',
                         'delta cell',
                         'PP cell')
library(Seurat)
pc12.heatmap.col <- colorRampPalette(c( "midnightblue","white","darkorange2","red","darkred"), 
                                     space="Lab")(20)

keep.gene <- rownames(beta.seu)[rowSums(as.matrix(beta.seu@assays$RNA@data))>0]
keep.gene <- setdiff(keep.gene,c(beta.cc.sym,
                                 beta.ex.sym))
length(keep.gene)#17995
beta.seu <- SetIdent(beta.seu,value = beta.seu$betagroup)

beta.seu.g23.diffmarker <- FindMarkers(beta.seu,
                                       # only.pos = T,
                                       features = keep.gene,
                                       ident.1 = 'beta cell2',
                                       ident.2 = c('beta cell3'),
                                      # test.use = 'roc'#,
                                       #return.thresh = 0.7
                                       logfc.threshold = log(1.2)
)
dim(beta.seu.g23.diffmarker)
beta.seu.g23.diffmarker <- MyReadDelim('beta/DEG/beta.seu.g23.diffmarker.tab')
rownames(beta.seu.g23.diffmarker) <- beta.seu.g23.diffmarker$gene

g2g3.preg.high.sym <- rownames(beta.seu.g23.diffmarker)[beta.seu.g23.diffmarker$avg_logFC<0]
g2g3.preg.low.sym <- rownames(beta.seu.g23.diffmarker)[beta.seu.g23.diffmarker$avg_logFC>0]

beta.seu.g23.diffmarker <- na.omit(beta.seu.g23.diffmarker)

beta.seu.g23.diffmarker$cluster <- NA
beta.seu.g23.diffmarker$cluster <- 'beta2'
beta.seu.g23.diffmarker[beta.seu.g23.diffmarker$avg_logFC<0,'cluster'] <- 'beta3'
beta.seu.g23.diffmarker$gene <- rownames(beta.seu.g23.diffmarker)
beta.seu.g23.diffmarker.filter <- beta.seu.g23.diffmarker[beta.seu.g23.diffmarker$p_val_adj<=0.05&
                                                           # abs(beta.seu.g23.diffmarker$avg_logFC)>=log(1.6) & 
                                                            (beta.seu.g23.diffmarker$pct.1>=0.3&beta.seu.g23.diffmarker$pct.2>=0.3),]
table(beta.seu.g23.diffmarker.filter$cluster)


MyWriteTable(cbind(beta.seu.g23.diffmarker,
                   genes.inf.input[rownames(beta.seu.g23.diffmarker),]),'beta/DEG/beta.seu.g23.diffmarker.tab')


c1Color <- MyName2Col(beta.seu@meta.data[,"Time"],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[,"betagroup"],
                      endoType.col)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,
                c2Color)
pdf("wilcox.heatmap.log0.25.2.pdf",10,12)
#par(oma = c(0,0,0,4))
beta.g2g3.DEGs.row.tree <- MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[rownames(beta.seu.g23.diffmarker),],
                                     type = "log.row.zscore",
                                     hc.c.data.type = "log.row.relat",
                                     hc.r.data.type = "log.row.relat",
                                     c.cov.method = "p",
                                     r.cov.method = "p",
                                     c.hc.method = "ward.D2",
                                     # color.palette = colors.heat,
                                     r.hc.method = "ward.D2",
                                     ColSideColors = cColor,
                                     ColSideColorsSize = 2,
                                     # labRow = gene.var.co,
                                     #cexRow = 0.3,
                                     return.tree = "row"
)

r1Color <- MyName2Col(cutree(beta.g2g3.DEGs.row.tree,
                             2),
                      gene.tree.colors,
                      is.row = T)
r1Color <- as.matrix(r1Color)
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[rownames(beta.seu.g23.diffmarker),],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          # color.palette = colors.heat,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 2,
          RowSideColors = r1Color
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()
beta.g2g3.DEGs.row.tree.den <- as.dendrogram(beta.g2g3.DEGs.row.tree)

c1Color <- MyName2Col(beta.seu@meta.data[colnames(beta.seu)[order(beta.seu$Time)],"Time"],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[colnames(beta.seu)[order(beta.seu$Time)],"betagroup"],
                      endoType.col)
c2Color <- as.matrix(c2Color)
cColor <- cbind(c1Color,
                c2Color)


pdf('g2g3.heatmap.orderbyTime.3.pdf',10,12)
r1Color <- c(rep(endoType.col[3],length(g2g3.preg.high.sym)),
           #  rep(gene.tree.colors[2],length(beta.high.sym)),
           #  rep(gene.tree.colors[3],length(preg.low.sym)),
             rep(endoType.col[2],length(g2g3.preg.low.sym)))
r1Color <- as.matrix(t(r1Color))
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[c(g2g3.preg.high.sym,
                                                     g2g3.preg.low.sym),colnames(beta.seu)[order(beta.seu$Time)]],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          Rowv = 'none',
          Colv = 'none',
          dendrogram = 'none',
          #color.palette = colors.exp,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 1.5,
          RowSideColors = r1Color,
          RowSideColorsSize = 1.5
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()

pdf('g2g3.heatmap.orderbyhccluster.pdf',10,12)
r1Color <- c(rep(gene.tree.colors[1],length(g2g3.preg.low.sym)),
             #  rep(gene.tree.colors[2],length(beta.high.sym)),
             #  rep(gene.tree.colors[3],length(preg.low.sym)),
             rep(gene.tree.colors[2],length(g2g3.preg.high.sym)))
r1Color <- as.matrix(t(r1Color))
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[c(g2g3.preg.low.sym,
                                                     g2g3.preg.high.sym),colnames(beta.seu)[order(beta.seu$Time)]],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          Rowv = 'none',
          Colv = 'none',
          dendrogram = 'none',
          color.palette = colors.exp,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 1.5,
          RowSideColors = r1Color,
          RowSideColorsSize = 1.5
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()
######ctrl g1/g2g3#####
beta.seu <- readRDS('beta/src.beta.rds')
beta.seu$time.group <- paste(beta.seu$Time,beta.seu$betagroup,sep = '_')
beta.seu <- SetIdent(beta.seu,value = beta.seu$time.group)

beta.seu.g1g23.ctrl.diffmarker <- FindMarkers(beta.seu,
                                              features = keep.gene,
                                              # only.pos = T,
                                              ident.1 = 'Ctrl_beta cell1',
                                              ident.2 = c('Ctrl_beta cell2','Ctrl_beta cell3')#,
                                              #  test.use = 'roc'#,
                                              # return.thresh = 0.7
                                              #logfc.threshold = log(1.5),
)

dim(beta.seu.g1g23.ctrl.diffmarker)


ctrl.g1.g23.g1.low.sym <- rownames(beta.seu.g1g23.ctrl.diffmarker)[beta.seu.g1g23.ctrl.diffmarker$avg_logFC<0]
ctrl.g1.g23.g1.high.sym <- rownames(beta.seu.g1g23.ctrl.diffmarker)[beta.seu.g1g23.ctrl.diffmarker$avg_logFC>0]

beta.seu.g1g23.ctrl.diffmarker$cluster <- NA
beta.seu.g1g23.ctrl.diffmarker[ctrl.g1.g23.g1.low.sym,'cluster'] <- 'beta23'
beta.seu.g1g23.ctrl.diffmarker[ctrl.g1.g23.g1.high.sym,'cluster'] <- 'beta1'
beta.seu.g1g23.ctrl.diffmarker$gene <- rownames(beta.seu.g1g23.ctrl.diffmarker)
rownames(genes.inf.input) <- genes.inf.input$SymbolDedu

rownames(genes.inf.input) <- sub('-','_',rownames(genes.inf.input))

MyWriteTable(cbind(beta.seu.g1g23.ctrl.diffmarker,
                   genes.inf.input[rownames(beta.seu.g1g23.ctrl.diffmarker),]),'beta/DEG/beta.seu.g1g23.ctrl.diffmarker.tab')

c1Color <- MyName2Col(beta.seu@meta.data[colnames(beta.seu)[beta.seu$Time=='Ctrl'],"Time"],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[colnames(beta.seu)[beta.seu$Time=='Ctrl'],"betagroup"],
                      endoType.col)
c2Color <- as.matrix(c2Color)
cColor <- cbind(c1Color,
                c2Color)



pdf('ctrl.g1.g2g3.heatmap.pdf',10,12)
r1Color <- c(rep(gene.tree.colors[1],length(ctrl.g1.g23.g1.high.sym)),
             #  rep(gene.tree.colors[2],length(beta.high.sym)),
             #  rep(gene.tree.colors[3],length(preg.low.sym)),
             rep(gene.tree.colors[2],length(ctrl.g1.g23.g1.low.sym)))
r1Color <- as.matrix(t(r1Color))
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[c(ctrl.g1.g23.g1.high.sym,
                                                     ctrl.g1.g23.g1.low.sym),colnames(beta.seu)[beta.seu$time.group %in% c('Ctrl_beta cell1',
                                                                                                                           'Ctrl_beta cell2',
                                                                                                                           'Ctrl_beta cell3')]],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          Rowv = 'none',
          Colv = 'none',
          dendrogram = 'none',
          #color.palette = colors.exp,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 2,
          RowSideColors = r1Color,
          RowSideColorsSize = 1.5
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()

ctrl.meta.data <- beta.seu@meta.data[beta.seu$Time=='Ctrl',]

c1Color <- MyName2Col(ctrl.meta.data[rownames(ctrl.meta.data)[order(ctrl.meta.data$betagroup)],"Time"],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(ctrl.meta.data[rownames(ctrl.meta.data)[order(ctrl.meta.data$betagroup)],"betagroup"],
                      endoType.col)
c2Color <- as.matrix(c2Color)
cColor <- cbind(c1Color,
                c2Color)

c2Color <- c(rep(gene.tree.colors[1],length(colnames(beta.seu)[beta.seu$time.group %in% c('Ctrl_beta cell1')])),
             rep(gene.tree.colors[2],length(colnames(beta.seu)[beta.seu$time.group %in% c('Ctrl_beta cell2','Ctrl_beta cell3')])))
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,
                c2Color)

pdf('reorder.ctrl.g1.g2g3.heatmap.ctrl.pdf',10,12)
r1Color <- c(rep(gene.tree.colors[1],length(ctrl.g1.g23.g1.high.sym)),
             #  rep(gene.tree.colors[2],length(beta.high.sym)),
             #  rep(gene.tree.colors[3],length(preg.low.sym)),
             rep(gene.tree.colors[2],length(ctrl.g1.g23.g1.low.sym)))
r1Color <- as.matrix(t(r1Color))
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[c(ctrl.g1.g23.g1.high.sym,
                                                     ctrl.g1.g23.g1.low.sym),rownames(ctrl.meta.data)[order(ctrl.meta.data$betagroup)]],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          Rowv = 'none',
          Colv = 'none',
          dendrogram = 'none',
          #color.palette = colors.exp,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 2,
          RowSideColors = r1Color,
          RowSideColorsSize = 1.5
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()
#########G0 G14.5#####
beta.seu <- readRDS('../')
beta.seu <- SetIdent(beta.seu,value = beta.seu$time.group)
G0G14.5.Glut2H.DEG.tab <- Myseufindmarker(seu.ob = beta.seu,gene.include = keep.gene,ident.1 = c('Ctrl_beta cell2',
                                                                                                  'Ctrl_beta cell3'),ident.2 = c('G14.5_1st_beta cell2','G14.5_1st_beta cell3'),c1 = 'Virgin_Glut2H',c2 = 'G14.5_Glut2H')
MyWriteTable(cbind(G0G14.5.Glut2H.DEG.tab,
                   genes.inf.input[rownames(G0G14.5.Glut2H.DEG.tab),]),'beta/DEG/G0G14.5.Glut2H.DEG.tab')
G0G14.5.DEG.10x.filter <-  G0G14.5.Glut2H.DEG.tab[ G0G14.5.Glut2H.DEG.tab$p_val_adj<=0.05&
                                                     abs( G0G14.5.Glut2H.DEG.tab$avg_logFC) >= log(1.2) & 
                                                     ( G0G14.5.Glut2H.DEG.tab$pct.1>=0.3 |  G0G14.5.Glut2H.DEG.tab$pct.2>=0.3),]

G0G14.5.DEG.10x.filter <- G0G14.5.DEG.10x.filter[order(G0G14.5.DEG.10x.filter$avg_logFC,decreasing = T),]
G0G14.5.DEG.10x.G14.5.tab <- G0G14.5.DEG.10x.filter[G0G14.5.DEG.10x.filter$cluster=='G14.5_Glut2H',]
G0G14.5.DEG.10x.G14.5.tab <- G0G14.5.DEG.10x.G14.5.tab[order(G0G14.5.DEG.10x.G14.5.tab$avg_logFC,decreasing = F),]

G0G14.5.DEG.10x.G0.tab <- G0G14.5.DEG.10x.filter[G0G14.5.DEG.10x.filter$cluster=='Virgin_Glut2H',]
G0G14.5.DEG.10x.filter <- rbind(G0G14.5.DEG.10x.G0.tab,G0G14.5.DEG.10x.G14.5.tab)
MyWriteTable(cbind(genes.inf.input[G0G14.5.DEG.10x.filter$gene,c('EnsemblGeneID','Symbol')],G0G14.5.DEG.10x.filter),'beta/DEG/Figs1g.G0G14.5.DEG.10x.filter.fc.tab')
G0G14.5.Glut2H.DEG.tab$cluster <- factor(G0G14.5.Glut2H.DEG.tab$cluster,levels = c('Virgin_Glut2H','G14.5_Glut2H'))
G0G14.5.DEG.10x.filter$cluster <- factor(G0G14.5.DEG.10x.filter$cluster,levels = c('Virgin_Glut2H','G14.5_Glut2H'))


MyDEGfilterplot(seu,G0G14.5.DEG.10x.filter,G0G14.5.Glut2H.DEG.tab,prefix = paste('beta/DEG/Figs1g.','padj0.05','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = preg.colors,text.cex=1,pval.gap = 50,x.gap = 0.2,x.ext = -1,ext = -150,point.size = 2, width.input=7,height.input=8)
MyDEGfilterplot(seu,G0G14.5.DEG.10x.filter,G0G14.5.Glut2H.DEG.tab,prefix = paste('beta/DEG/text.','padj0.05','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = preg.colors,text.cex=1,y.cex = 0.01,pval.gap = 50,x.gap = 0.2,x.ext = -1,ext = -150,point.size = 2, width.input=7,height.input=8)



beta.seu <- SetIdent(beta.seu,value = beta.seu$beta)
Glut2H.1.2.DEG.tab <- Myseufindmarker(seu.ob = beta.seu,gene.include = keep.gene,ident.1 = c('Glut2H_1'),ident.2 = c('Glut2H_2'),c1 = 'Glut2H_1',c2 = 'Glut2H_2')
Glut2H.1.2.DEG.tab$cluster <- as.character(Glut2H.1.2.DEG.tab$cluster)
Glut2H.1.2.DEG.tab[Glut2H.1.2.DEG.tab$avg_logFC>0,'cluster'] <- 'Glut2H_1'
Glut2H.1.2.DEG.tab[Glut2H.1.2.DEG.tab$avg_logFC<0,'cluster'] <- 'Glut2H_2'


Glut2H.1.2.DEG.tab.filter <-  Glut2H.1.2.DEG.tab[Glut2H.1.2.DEG.tab$p_val_adj<=0.01&
                                                     abs( Glut2H.1.2.DEG.tab$avg_logFC) >= log(1.2) & 
                                                     ( Glut2H.1.2.DEG.tab$pct.1>=0.3 |  Glut2H.1.2.DEG.tab$pct.2>=0.3),]

table(Glut2H.1.2.DEG.tab.filter$cluster)

Glut2H.1.2.DEG.tab$cluster <- factor(Glut2H.1.2.DEG.tab$cluster,levels = c('Virgin_Glut2H','G14.5_Glut2H'))
Glut2H.1.2.DEG.tab$cluster <- factor(Glut2H.1.2.DEG.tab$cluster,levels = c('Virgin_Glut2H','G14.5_Glut2H'))

Glut2H.1.2.DEG.tab.filter.cp <- Glut2H.1.2.DEG.tab.filter
Glut2H.1.2.DEG.tab.cp <- Glut2H.1.2.DEG.tab
Glut2H.1.2.DEG.tab.cp[Glut2H.1.2.DEG.tab.cp$p_val==0,'p_val'] <- 10e-300
Glut2H.1.2.DEG.tab.filter.cp[Glut2H.1.2.DEG.tab.filter.cp$p_val==0,'p_val'] <- 10e-300


MyDEGfilterplot(seu,Glut2H.1.2.DEG.tab.filter.cp,Glut2H.1.2.DEG.tab.cp,prefix = paste('cp.Figs1g.GLUT2H_12.','padj0.01','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = celltype.col2[3:4],text.cex=0.5,pval.gap = 50,x.gap = 0.5,x.ext = 0,ext = 50,point.size = 2)
MyDEGfilterplot(seu,Glut2H.1.2.DEG.tab.filter.cp,Glut2H.1.2.DEG.tab.cp,prefix = paste('cp.Figs1g.text.GLUT2H_12.','padj0.05','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = celltype.col2[3:4],text.cex=0.5,y.cex = 0.01,pval.gap = 50,x.gap = 0.5,x.ext = 0,ext = 50,point.size = 2,width.input=7,height.input=6)

MyWriteTable(cbind(genes.inf.input[Glut2H.1.2.DEG.tab.filter$gene,c('EnsemblGeneID','Symbol')],Glut2H.1.2.DEG.tab.filter),'beta/DEG/Figs1g.Gltu2H1_2.DEG.10x.filter.fc.tab')
rm(beta.seu)


Glut2L.1.2.DEG.tab <- Myseufindmarker(seu.ob = beta.seu,gene.include = keep.gene,ident.1 = c('Glut2L_1'),ident.2 = c('Glut2L_2'),c1 = 'Glut2L_1',c2 = 'Glut2L_2')

Glut2L.1.2.DEG.tab.filter <-  Glut2L.1.2.DEG.tab[Glut2L.1.2.DEG.tab$p_val_adj<=0.01&
                                                   abs( Glut2L.1.2.DEG.tab$avg_logFC) >= log(1.2) & 
                                                   ( Glut2L.1.2.DEG.tab$pct.1>=0.3 |  Glut2L.1.2.DEG.tab$pct.2>=0.3),]
table(Glut2L.1.2.DEG.tab.filter$cluster)

Glut2L.1.2.DEG.tab$cluster <- factor(Glut2L.1.2.DEG.tab$cluster,levels = c('Glut2L_1','Glut2L_2'))
Glut2L.1.2.DEG.tab$cluster <- factor(Glut2L.1.2.DEG.tab$cluster,levels = c('Glut2L_1','Glut2L_2'))

MyDEGfilterplot(seu,Glut2L.1.2.DEG.tab.filter,Glut2L.1.2.DEG.tab,prefix = paste('G:/lab/Article/heshuang/BYLW/10x/Figs1g.GLUT2L_12.2.','padj0.01','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = celltype.col2[c(1,2)],text.cex=1,pval.gap = 5,x.gap = 0.2,x.ext = -1.2,ext = 0,point.size = 2)
MyDEGfilterplot(seu,Glut2L.1.2.DEG.tab.filter,Glut2L.1.2.DEG.tab,prefix = paste('G:/lab/Article/heshuang/BYLW/10x/Figs1g.text.GLUT2_12.2.','padj0.01','fc1.2',sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = celltype.col2[c(1,2)],text.cex=1,y.cex = 0.01,pval.gap = 5,x.gap = 0.2,x.ext = -1.2,ext = 0,point.size = 2, width.input=7,height.input=6)

# MyDEGfilterplot(beta.seu,DEG.raw.filter = seu,Glut2L.1.2.DEG.tab.filter,DEG.raw = Glut2L.1.2.DEG.tab,'beta/DEG/Glut2L1_2.',plotheatmap = T,display.cells = colnames(beta.seu)[beta.seu$beta %in% c('Glut2L_1','Glut2L_2')],plotfireplot = F,
#                 group = c('Time','beta')
# )

MyWriteTable(cbind(genes.inf.input[Glut2L.1.2.DEG.tab.filter$gene,c('EnsemblGeneID','Symbol')],Glut2L.1.2.DEG.tab.filter),'beta/DEG/Figs1g.Gltu2L1_2.DEG.10x.filter.fc.tab')
rm(beta.seu)

save.image('beta.10x.DEG.RData')
##########
seu.beta.L <- readRDS('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/endo/Glut2L/seu.Glut2L.rds')
beta.seu$beta <- '/'
beta.seu@meta.data[colnames(seu.beta.L)[seu.beta.L$beta=='beta cell2_1'],'beta'] <- 'Glut2L_1'
beta.seu@meta.data[colnames(seu.beta.L)[seu.beta.L$beta=='beta cell2_2'],'beta'] <- 'Glut2L_2'
beta.seu@meta.data[beta.seu$betagroup=='beta cell3','beta'] <- 'Glut2H_2'
beta.seu@meta.data[beta.seu$betagroup=='beta cell2','beta'] <- 'Glut2H_1'

table(endo.seu$beta)
beta.seu$beta <- factor(beta.seu$beta,levels = c('Glut2H_1','Glut2H_2','Glut2L_1','Glut2L_2'))
table(beta.seu$beta,beta.seu$Time)
celltype.col2 <- c("#5F9EA0", "#A65628","#6a3d9a",colors.group[1])
names(celltype.col2) <- c('Glut2H_1','Glut2H_2','Glut2L_1','Glut2L_2')
MyPlotColor(celltype.col2,4)

c1Color <- MyName2Col(beta.seu@meta.data[order(beta.seu$Time,decreasing = F),'Time'],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[order(beta.seu$Time,decreasing = F),'beta'],
                      celltype.col2)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)


Virgin.gene.tab <- G0G14.5.DEG.10x.filter[G0G14.5.DEG.10x.filter$cluster=='Virgin_Glut2H',]
G14.5.gene.tab <- G0G14.5.DEG.10x.filter[G0G14.5.DEG.10x.filter$cluster=='G14.5_Glut2H',]
Virgin.gene <- Virgin.gene.tab$gene[order(Virgin.gene.tab$p_val,decreasing = F)]
G14.5.gene <- G14.5.gene.tab$gene[order(G14.5.gene.tab$p_val,decreasing = F)]

colname.order <- c(labels(as.dendrogram(DEG.col.tree)[[1]]),labels(as.dendrogram(DEG.col.tree)[[2]][[1]]),
                  labels(as.dendrogram(DEG.col.tree)[[2]][[2]]))

c1Color <- MyName2Col(beta.seu@meta.data[colname.order,'Time'],
                      preg.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[colname.order,'beta'],
                      celltype.col2)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)


Virgin.gene.order<- MyordergenewithPseudotime(as.matrix(exp(beta.seu@assays$RNA@scale.data))[Virgin.gene,colname.order],
                                          graph = T,
                                          Virgin.gene)
G14.5.gene.order<- MyordergenewithPseudotime(as.matrix(exp(beta.seu@assays$RNA@scale.data))[G14.5.gene,colname.order],
                                              graph = T,
                                             G14.5.gene)

library(Seurat)
#beta.seu <- ScaleData(beta.seu,features = G0G14.5.DEG.10x.filter$gene)
pdf('beta/DEG/FigS1g.10x.G0G14.5.GlutH.DEG.padj.0.05.fc1.2.heatmap.scale.5.pdf',10,12)
#DEG.col.tree <- 
  MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@scale.data))[c(rev(Virgin.gene.order)[sample(1:length(Virgin.gene.order))],rev(G14.5.gene.order)),colname.order],
          type = 'log.row.relat',
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          color.palette = pc12.heatmap.col,
          c.cov.method = "s",
          r.cov.method = "s",
          Colv = 'do',
          Rowv='none',
          dendrogram='col',
          c.hc.method = "ward.D2",
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 1.5,
          #RowSideColorsSize = 1.5,
          RowSideColors = MyName2Col(sort(G0G14.5.DEG.10x.filter$cluster,decreasing = T),rev(preg.colors),is.row = T),
         # return.tree = "col",
          graph = T
)
dev.off()

#######g1g23#########
beta.seu.g1g23.diffmarker <- FindMarkers(beta.seu,
                                         features = keep.gene,
                                         # only.pos = T,
                                         ident.1 = 'beta cell1',
                                         ident.2 = c('beta cell2','beta cell3')#,
                                       #  test.use = 'roc'#,
                                         # return.thresh = 0.7
                                         #logfc.threshold = log(1.5),
)

dim(beta.seu.g1g23.diffmarker)

pdf("wilcox.g1vsg2g3.heatmap.pdf",10,12)
#par(oma = c(0,0,0,4))
beta.g1.g2g3.DEGs.filter.row.tree <- MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[rownames(beta.seu.g1g23.diffmarker),],
                                               type = "log.row.relat",
                                               hc.c.data.type = "log.row.relat",
                                               hc.r.data.type = "log.row.relat",
                                               c.cov.method = "p",
                                               r.cov.method = "p",
                                               c.hc.method = "ward.D2",
                                               # color.palette = colors.heat,
                                               r.hc.method = "ward.D2",
                                               ColSideColors = c2Color,
                                               ColSideColorsSize = 2,
                                               # labRow = gene.var.co,
                                               #cexRow = 0.3,
                                               return.tree = "row"
)

# r1Color <- MyName2Col(cutree(beta.g1.g2g3.DEGs.filter.row.tree,
#                              2),
#                       gene.tree.colors,
#                       is.row = T)
# r1Color <- as.matrix(r1Color)
# 
# MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[rownames(beta.seu.g1g23.diffmarker.filter),beta.seu@meta.data$SampleName],
#           type = "log.row.zscore",
#           hc.c.data.type = "log.row.relat",
#           hc.r.data.type = "log.row.relat",
#           c.cov.method = "p",
#           r.cov.method = "p",
#           c.hc.method = "ward.D2",
#           # color.palette = colors.heat,
#           r.hc.method = "ward.D2",
#           ColSideColors = cColor,
#           ColSideColorsSize = 2,
#           RowSideColors = r1Color
#           #labRow = gene.var.co,
#           #cexRow = 0.3,
#           # return.tree = "row"
# )
dev.off()
##########g1 g2 g3#####
g1g2g3.gene <- unique(c(rownames(beta.seu.g1g23.diffmarker),
                        rownames(beta.seu.g23.diffmarker)))
length(g1g2g3.gene)#124

pdf("g1.g2g3.g2g3.wilcox.heatmap.2.pdf",10,12)
#par(oma = c(0,0,0,4))
beta.g1g2g3.DEGs.filter.row.tree <- MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[g1g2g3.gene,],
                                              type = "log.row.relat",
                                              hc.c.data.type = "log.row.relat",
                                              hc.r.data.type = "log.row.relat",
                                              c.cov.method = "p",
                                              r.cov.method = "p",
                                              c.hc.method = "ward.D2",
                                              # color.palette = colors.heat,
                                              r.hc.method = "ward.D2",
                                              ColSideColors = cColor,
                                              ColSideColorsSize = 2,
                                              # labRow = gene.var.co,
                                              #cexRow = 0.3,
                                              return.tree = "row"
)
r1Color <- MyName2Col(cutree(beta.g1g2g3.DEGs.filter.row.tree,
                             4),
                      gene.tree.colors,
                      is.row = T)
r1Color <- as.matrix(r1Color)

MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@data))[g1g2g3.gene,],
          type = "log.row.zscore",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          # color.palette = colors.heat,
          r.hc.method = "ward.D2",
          ColSideColors = cColor,
          ColSideColorsSize = 2,
          RowSideColors = r1Color
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()

beta.g1g2g3.DEGs.filter.row.tree.den <- as.dendrogram(beta.g1g2g3.DEGs.filter.row.tree)
labels(beta.g1g2g3.DEGs.filter.row.tree.den[[2]][[1]])
labels(beta.g1g2g3.DEGs.filter.row.tree.den[[2]][[2]])
#########reorder heatmap#########
preg.high.sym <- labels(beta.g1g2g3.DEGs.filter.row.tree.den[[1]])
preg.low.sym <- labels(beta.g1g2g3.DEGs.filter.row.tree.den[[2]][[2]][[2]])
beta.low.sym <- labels(beta.g1g2g3.DEGs.filter.row.tree.den[[2]][[1]])
beta.high.sym <- labels(beta.g1g2g3.DEGs.filter.row.tree.den[[2]][[2]][[1]])

beta.reorder.sym <- c(beta.low.sym,beta.high.sym,preg.low.sym,preg.high.sym)
# c1Color <- MyName2Col(beta.seu@meta.data[,"Type3"],
#                       time.colors)
# c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(beta.seu@meta.data[colnames(beta.seu)[order(beta.seu$betagroup)],"betagroup"],
                      endoType.col)
c2Color <- as.matrix(c2Color)
# cColor <- cbind(c1Color,
#                 c2Color)


pdf('../g1.g2g3.g2g3.wilcox.log1.5.pct.2.0.2.padj5.heatmap.pdf')
r1Color <- c(rep(gene.tree.colors[1],length(beta.low.sym)),
             rep(gene.tree.colors[2],length(beta.high.sym)),
             rep(gene.tree.colors[3],length(preg.low.sym)),
             rep(gene.tree.colors[4],length(preg.high.sym)))
r1Color <- as.matrix(t(r1Color))
MyHeatmap(as.matrix(exp(beta.seu@assays$RNA@scale.data))[beta.reorder.sym,colnames(beta.seu)[order(beta.seu$betagroup)]],
          type = "log.row.relat",
          hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          c.cov.method = "p",
          r.cov.method = "p",
          c.hc.method = "ward.D2",
          Rowv = 'none',
          Colv = 'none',
          dendrogram = 'none',
          color.palette = pc12.heatmap.col,
          r.hc.method = "ward.D2",
          ColSideColors = c2Color,
          ColSideColorsSize = 1.5,
          RowSideColors = r1Color,
          RowSideColorsSize = 1.5
          #labRow = gene.var.co,
          #cexRow = 0.3,
          # return.tree = "row"
)
dev.off()
#########
rm(beta.seu)
saveRDS(beta.seu,'beta/src.beta.rds')
rm(beta.seu)
save.image('beta.10x.DEG.RData')

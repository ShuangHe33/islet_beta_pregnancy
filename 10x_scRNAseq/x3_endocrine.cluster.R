setwd('G:/project/pregnant_mouse/10x/10x_v3')
dir.create('Ctl_G14.5_1st/')
setwd('Ctl_G14.5_1st/')
endo.seu <- readRDS('endo/endo.seu.rds')


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
dir.create('endo')
#############
pre.seu <- readRDS('pre.seu.rds')
endo.seu <- subset(pre.seu,cells = colnames(pre.seu)[pre.seu$merge.cluster %in% c('beta','alpha','delta','pp')])
rm(pre.seu)

endo.seu <- 
  FindVariableFeatures(endo.seu, selection.method = "vst", nfeatures = 2000)

endo.var.co <- rownames(MyCo(as.matrix(endo.seu@assays$RNA@data),
                            var.gene = VariableFeatures(endo.seu),
                            exp.prop.whole.max = 0.9,
                            exp.prop.whole.min = 0.005,
                            # vector.group = samples.inf.qc$GroupNameFig1,
                            # exp.prop.group.min = 0.1,
                            # exp.prop.group.max = 0.5,
                            cor.method = "rho",
                            cor.cutoff = 0.15,
                            partner.cutoff = 5,
                            refine.cor = T))
length(endo.var.co)#1264
# c1Color <- MyName2Col(endo.seu@meta.data[,"RNA_snn_res.0.5"],
#                       brewer.pal(9,"Set1"))
# c1Color <- as.matrix(c1Color)

png('endo/varcor0.15.png',2000,3000)
var.row.tree <-
  MyHeatmap(as.matrix(exp(endo.seu@assays$RNA@data)[endo.var.co,]),
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
grep('Mki67',labels(var.row.tree[[1]][[2]][[2]][[1]]))

endo.cc.sym <- labels(var.row.tree[[1]][[2]][[2]])
MyGOwritetable(genes.inf.input[endo.cc.sym,1],
               pvalue = 1,
               'endo/go.endo.cc.2.tsv')


endo.seu <- 
  ScaleData(endo.seu)
endo.seu <- RunPCA(endo.seu, features = setdiff(endo.var.co,
                                              endo.cc.sym))
######eound2 exclude cc gene#########
pdf('endo/DimHeatmap.rmcc.pdf',
    10,40)
DimHeatmap(endo.seu, dims = 1:30, cells = 200, balanced = TRUE)
dev.off()

pc.use <- 1:15
endo.seu <- FindNeighbors(endo.seu, dims = pc.use)

endo.seu <- FindClusters(endo.seu, resolution = 1)
endo.seu <- FindClusters(endo.seu, resolution = 0.5)
endo.seu <- FindClusters(endo.seu, resolution = 0.3)
endo.seu <- FindClusters(endo.seu, resolution = 0.1)

endo.seu <- RunTSNE(endo.seu,
                   dims = pc.use,
                   perplexity= round((30+ncol(endo.seu)/100)),
                   check_duplicates = F)
endo.seu <-
  SetIdent(endo.seu, value = endo.seu$RNA_snn_res.0.4)


pdf('endo/rmcc.tsne.pc15.res0.4.pdf',
    10,8)
DimPlot(endo.seu,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 1,
        label = T)
dev.off()

p.tsne <- MySeuratv3TSNE10x2Gg(endo.seu,
                               endo.seu@meta.data
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

pdf("endo/tsne.cellType.pc15.pdf",
    14,
    10)
p.plot <- p.age +
  scale_color_manual(values = colors.group) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =merge.cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()

pdf("endo/tsne.Time.pc15.pdf",
    15,
    10)
p.plot <- p.age +
  scale_color_manual(values = preg.colors) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =Time#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()

endo.seu <- RunUMAP(endo.seu,
                    dims = pc.use#,
                    #perplexity= round((30+ncol(endo.seu)/100)),
                    #check_duplicates = F
                    )

p.tsne <- MySeuratv3UMAP10x2Gg(endo.seu,
                               endo.seu@meta.data
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

pdf("endo/umap.cellType.pc15.pdf",
    14,
    10)
p.plot <- p.age +
  scale_color_manual(values = colors.group) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =merge.cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()
pdf("endo/umap.Time.pc15.pdf",
    15,
    10)
p.plot <- p.age +
  scale_color_manual(values = preg.colors) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =Time#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()
########marker########
marker.sym <- c(
                # "Spi1",#°×Ï¸°û
                # "Procr",
                # "Fcgr1",
                # "Neurod1",#endocrine
                # "Chga",
                # "Chgb",
                # "Prox1",
                # "Mki67",
                # "Col1a1",#mesenchymal
                # "Col5a1",#mesenchymal
                # "Col3a1",#mesenchymal
                # "Mgp",
                # "Esm1",
                # "Vim",
                # "Cdh11",
                # #"S100a4",
                # "Itga5",
                # "Sdc1",
                # "Pecam1",#Ñª¹Ü
                # "Gypa",#ºìÏ¸°û
                # "Sox17",#¸ÎÍâµ¨¹Ü
                # "Sox9",#duct
                # "Spp1",#duct
                # "Rbpjl",#acinar
                # # "Amy2b",#acinar
                # "Pnlip",#exo
                # "Ctrb1",#exo
                # "Prss1",
                # "Neurog3",
                # "Pax4",
                "Ins1",
                "Gcg",
                "Sst",
                "Ppy",
                "Ghrl",
                "Nkx6-1",
                "Ucn3",
                "Slc2a2",
                'G6pc2',
                "Mafb",
                "Mafa",
                'Chgb',
                'Fmo1',
                'Ttr',
                # "Meg3",
                # "Cacna1a",
                # "Xist",
                # "Kcnq1ot1",
                "Arx",
                "Hhex"#,
                # "Kdr",
                # "Ascl1",
                # "Ptprc"
                
                
)
length(marker.sym)#
marker.sym <- marker.sym[marker.sym %in% rownames(endo.seu)]
length(marker.sym)#
marker.sym <- unique(marker.sym)
length(marker.sym)#45
total.count <- length(marker.sym)
run.count <- 1


marker.sym <- c('Ucn3','Slc2a2','Chgb','Il1r1','Ttr','Ivd','Tph1')

marker.sym <- c('Ins1','Gcg','Sst','Ppy','Cd81','Tspan8','Ucn3','Slc2a2','Il1r1','Ttr','Rgs2','Ovol2','Chgb','Ivd','Tph1','Tph2')


pdf(paste('G:/lab/Article/heshuang/BYLW/10x/pc15.marker',"beta.umap.pdf",sep = ""),6,7)
Myseuratmarker(beta.seu,
               marker.sym = marker.sym,reduction = 'umap',pt.size = 1,col = colors.exp)
dev.off()


pdf(paste('endo/pc15.marker',".tsne.2.pdf",sep = ""),6,7)

for (gene in c('percent.mt','nFeature_RNA',marker.sym)) {
  cat(paste(run.count, "/",total.count,"\n"))
  print(FeaturePlot(endo.seu,gene,reduction = 'tsne')+scale_color_gradientn(colours = colors.exp)+
          theme(plot.title = element_text(size = rel(3.5))) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.title = element_blank()) +
          theme(axis.text  = element_blank()) +
          theme(axis.ticks = element_blank()) +
          theme(legend.position="bottom")+
          theme(panel.border = element_rect(size = 4,
                                            colour = "black"))+
          guides(colour = guide_colorbar(title = "scaled(expression)",
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
######celltype######
#######cell type#######
endo.seu$endoType <- as.character(endo.seu$RNA_snn_res.0.4)
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(7),'endoType'] <- 'PP cell'
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(1),'endoType'] <- 'alpha cell'
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(5),'endoType'] <- 'delta cell'
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(6),'endoType'] <- 'beta cell1'
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(2,3),'endoType'] <- 'beta cell3'
endo.seu@meta.data[endo.seu$RNA_snn_res.0.4 %in% c(0,4),'endoType'] <- 'beta cell2'

table(endo.seu$endoType)


seu.Glut2L <- readRDS('endo/Glut2L/seu.Glut2L.rds')
endo.seu$beta <- as.character(endo.seu$endoType)
endo.seu@meta.data[colnames(seu.Glut2L)[seu.Glut2L$beta=='beta cell2_1'],'beta'] <- 'Glut2L_1'
endo.seu@meta.data[colnames(seu.Glut2L)[seu.Glut2L$beta=='beta cell2_2'],'beta'] <- 'Glut2L_2'
endo.seu@meta.data[endo.seu$endoType=='beta cell2','beta'] <- 'Glut2H_1'
endo.seu@meta.data[endo.seu$endoType=='beta cell3','beta'] <- 'Glut2H_2'

table(endo.seu$beta)
endo.seu <- readRDS('endo/endo.seu.rds')
endo.seu$beta <- factor(endo.seu$beta,levels = c('Glut2H_1','Glut2H_2','Glut2L_1','Glut2L_2',
                                                 'alpha cell',
                                                 'delta cell',
                                                 'PP cell'))
celltype.col2 <- c("#6a3d9a","#E08698","#5F9EA0", "#A65628",colors.group[c(5,6,3)])
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2','alpha cell',
                          'delta cell',
                          'PP cell')
MyPlotColor(celltype.col2,10)


time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")


endo.seu$endoType <- factor(endo.seu$endoType,
                            levels = c('beta cell1',
                                       'beta cell2',
                                       'beta cell3',
                                       'alpha cell',
                                       'delta cell',
                                       'PP cell'))

endoType.col <- c(colors.group[9],"#5F9EA0", "#A65628",colors.group[c(5,6,3)])




names(endoType.col) <- c('beta cell1',
                         'beta cell2',
                         'beta cell3',
                         'alpha cell',
                         'delta cell',
                         'PP cell')
MyPlotColor(endoType.col,6)
endo.seu <- readRDS('endo/endo.seu.rds')



beta.seu <- subset(endo.seu,cells = colnames(endo.seu)[!endo.seu$endoType %in% c('alpha cell','delta cell','PP cell')])
beta.seu <- subset(beta.seu,cells = colnames(beta.seu)[!names(beta.seu$SampleName) %in% 'G14.5_1st_20200416TAACGACCATCGGATT???1'])

p.tsne <- MySeuratv3UMAP10x2Gg(beta.seu,
                               beta.seu@meta.data)
p.age <- p.tsne+
  #scale_y_reverse() +
 # scale_color_manual(values = endoType.col) +
  scale_shape_manual(values = shape.group) +
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
                               order = 1)) +
  guides(shape = guide_legend(title = "",
                              keywidth = 3,
                              keyheight = 3,
                              override.aes = list(size = 12),
                              order = 2))

pdf("G:/lab/Article/heshuang/BYLW/10x/Figs1.umap.beta.2.pdf",
    14.5,
    8)
p.plot <- p.age +
  scale_color_manual(values = endoType.col2) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col = endoType#,
                 # shape = State
  ),
  size = 2.5
  )+theme(aspect.ratio = 1)

print(p.plot)
dev.off()

setwd('G:/lab/Article/heshuang/BYLW/10x/')
pdf("G:/lab/Article/heshuang/BYLW/10x/Figs1.umap.endoType.Time.pdf",
    14.5,
    8)
p.plot <- p.age +
  scale_color_manual(values = c("#b6d4a8","#b89fc1")) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col = Time#,
                 # shape = State
  ),
  size = 2.5
  )+theme(aspect.ratio = 1)

print(p.plot)
dev.off()

celltype.col2 <- time.col[c(2,9,17,12)]
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2')
MyPlotColor(celltype.col2,10)


pdf("G:/lab/Article/heshuang/BYLW/10x/Figs1.umap.betaType.t.1.pdf",
    14.5,
    8)
p.plot <- p.age +
  scale_color_manual(values = celltype.col2) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col = beta#,
                 # shape = State
  ),
  size = 2.5
  )+theme(aspect.ratio = 1)


# p.plot <- p.age +
#   scale_color_manual(values = celltype.col2) +
#   geom_text(aes(x = x.pos,
#                  y = y.pos,
#                  col = beta#,
#                  # shape = State
#   ),
#   size = 1
#   )+theme(aspect.ratio = 1)


print(p.plot)
dev.off()

##############
save(endo.var.co,
     var.row.tree,
     #endo.seu,
     file = 'endo/endo.cluster.RData')

saveRDS(endo.seu,
        'endo/endo.seu.rds')
load('endo/endo.cluster.RData')

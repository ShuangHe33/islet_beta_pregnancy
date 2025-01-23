setwd('G:/project/pregnant_mouse/10x/10x_v3')
dir.create('QC')
load('QC/pre.QC.RData')

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(rgl)
library(scran)

genes.inf.input <- MyReadDelim("G:/lab/gene inf/mm10/mm10.gene.inf.merge.v1.1.tab")
rownames(genes.inf.input) <- genes.inf.input$SymbolDedu

time.colors <- c("#4DAF4A",
                 "#FF7F00",
                 "#7570B3", 
                 "#452cff",
                 "#A65628", 
                 "#F781BF",
                 "#999999",
                 "#ffe478",
                 "#cab2d6",
                 "#ba984d",
                 "#1f78b4",
                 "#ff4a21",
                 "#752a78",
                 "#a6cee3",
                 "#53751c",
                 "#C71585",
                 "#F01414",
                 "#ff8e9d",
                 "#234682",
                 "#ff723f",
                 "#d96ca8",
                 "#0f8575",
                 "#5277cf",
                 "#de562c",
                 "#5a7368",
                 "#c77fc1",
                 "#3fa0a3",
                 "#c7a96c",
                 "#805e2b",
                 "#F01414",
                 "#d15126",
                 "#6343c4"                             
)
marker.sym <- c("Spi1",
                "Procr",
                "Fcgr1",
                "Neurod1",
                "Chga",
                "Chgb",
                "Prox1",
                "Mki67",
                "Col1a1",
                "Col5a1",
                "Col3a1",
                "Mgp",
                "Esm1",
                "Vim",
                "Cdh11",
                "Itga5",
                "Sdc1",
                "Pecam1",
                "Gypa",
                "Sox17",
                "Sox9",
                "Spp1",
                "Rbpjl",
                "Amy2b",
                "Pnlip",
                "Ctrb1",
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

###################################################################
#                          load
###################################################################
load('rawdata/read.10x.G0G14.5G14.5_2nd.RData')
src.G0=CreateSeuratObject(src.raw.G0)
src.G14.5=CreateSeuratObject(src.raw.G14.5)

src.G0$Time="G0"
src.G14.5$Time="G14.5"


rm(src.raw.G0)
rm(src.raw.G14.5)

########################################################################
#                            QC
########################################################################
src.G0[["percent.mt"]] <- PercentageFeatureSet(src.G0, pattern = "^mt-")
src.G14.5[["percent.mt"]] <- PercentageFeatureSet(src.G14.5, pattern = "^mt-")

pdf('QC/features.pdf')
for(seu in c("src.G0",
             "src.G14.5"
             )){
  hist(get(seu)$nFeature_RNA,
       1000,
       xlim = c(0,10000),
       ylim = c(0,100),
       xlab = 'gene number',
       main = seu)
  
  hist(get(seu)$nCount_RNA,
       5000,
       xlim = c(0,30000),
       ylim = c(0,100),
       xlab = 'UMI count',
       main = seu)
  hist(get(seu)$percent.mt,
       500,
       xlim = c(0,100),
       ylim = c(0,2000),
       xlab = 'percent.mt',
       main = seu)
}
print(seu)
dev.off()


for (seurat in c("src.G0",
                 "src.G14.5"
                 
                 )) {
  assign(seurat,
         NormalizeData(get(seurat),
                       normalization.method = "LogNormalize",
                       scale.factor = 10000))
  # mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  # mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  # assign(seurat,
  #        AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
  #assign(seurat,
  #seurat[["percent.mt"]] <- PercentageFeatureSet(get(seurat), pattern = "^mt-"))
  assign(seurat,
         subset(get(seurat), subset = nFeature_RNA > 1500  & percent.mt < 10))
  assign(seurat,
         FindVariableFeatures(get(seurat), selection.method = "vst", nfeatures = 2000))
  VariableFeaturePlot(get(seurat))
  assign(seurat,
         ScaleData(get(seurat),features = VariableFeatures(get(seurat))))
  assign(seurat,
         MyDoubletFinder(get(seurat),round(ncol(get(seurat)) * 0.1,0)))
  # assign(paste(seurat,".row.tree",sep = ""),
  #        MyHeatmap(as.matrix(exp(get(seurat)@assays$RNA@data)[VariableFeatures(get(seurat)),]),
  #                  type = "log.row.relat",
  #                  hc.c.data.type = "log.row.relat",
  #                  hc.r.data.type = "log.row.relat",
  #                  c.cov.method = "s",
  #                  r.cov.method = "s",
  #                  c.hc.method = "ward.D2",
  #                  r.hc.method = "ward.D2",
  #                  return.tree = "row",
  #                  graph = F))
  # assign(paste(seurat,".row.tree",sep = ""),
  #        as.dendrogram(get(paste(seurat,".row.tree",sep = ""))))
  print(seurat)
}


for (seurat in c("src.G0",
                 "src.G14.5"
                 )){
  png(paste('QC/',seurat,'.pre.var.heatmap.2000.png'),
      2000,3000)
  assign(paste(seurat,".row.tree",sep = ""),
         MyHeatmap(as.matrix(exp(get(seurat)@assays$RNA@data)[VariableFeatures(get(seurat)),]),
                   type = "log.row.relat",
                   hc.c.data.type = "log.row.relat",
                   hc.r.data.type = "log.row.relat",
                   c.cov.method = "s",
                   r.cov.method = "s",
                   c.hc.method = "ward.D2",
                   r.hc.method = "ward.D2",
                   return.tree = "row",
                   graph = T)
  )
  assign(paste(seurat,".row.tree",sep = ""),
         as.dendrogram(get(paste(seurat,".row.tree",sep = ""))))
  print(seurat)
  dev.off()
}


VariableFeatures(src.G0) <- setdiff(VariableFeatures(src.G0),
                                            c(labels(src.G0.row.tree[[2]][[2]][[2]][[2]][[2]][[1]])))
VariableFeatures(src.G14.5) <- setdiff(VariableFeatures(src.G14.5),
                                            c(labels(src.G14.5.row.tree[[2]][[2]][[1]][[2]])))

for (seurat in c("src.G0",
                 "src.G14.5"
                )) {
  assign(seurat,
         RunPCA(get(seurat), features = VariableFeatures(get(seurat))))
  
  assign(seurat,
         FindNeighbors(get(seurat), dims = 1:30))
  assign(seurat,
         FindClusters(get(seurat), resolution = 3))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.5))
  assign(seurat,
         RunTSNE(get(seurat),
                 dims = 1:30,
                 perplexity= round((30+ncol(get(seurat))/100)),
                 check_duplicates = F))
  pdf(paste('QC/',seurat,".cluster.pdf",sep = ""),8,7)
  assign(seurat,
         SetIdent(get(seurat), value = get(seurat)$RNA_snn_res.1.5))
  print(TSNEPlot(get(seurat),label = T)+
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
          theme_bw() +
          scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"),
                                        brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"))) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1))
  assign(seurat,
         SetIdent(get(seurat), value = get(seurat)$RNA_snn_res.3))
  print(TSNEPlot(get(seurat),label = T)+
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
          theme_bw() +
          scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"),
                                        brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"))) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1))
  assign(seurat,
         SetIdent(get(seurat), value = get(seurat)$pANNPredictions))
  print(TSNEPlot(get(seurat))+
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
          theme_bw() +
          scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1))
  dev.off()
  
}

for (seurat in c("src.G0",
                  "src.G14.5"
                
)) {

pdf(paste('QC/',seurat,".marker.pdf",sep = ""),6,7)
for (gene in c("nFeature_RNA",marker.sym)) {
  print(FeaturePlot(get(seurat),gene,cols = c("#C7E8CC","#FF0000"))+
          theme(plot.title = element_text(size = rel(3.5))) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.title = element_blank()) +
          theme(axis.text  = element_blank()) +
          theme(axis.ticks = element_blank()) +
          theme(legend.position="bottom")+
          theme(panel.border = element_rect(size = 4,
                                            colour = "black"))+
          guides(colour = guide_colorbar(title = "ln(TP10K + 1)",
                                         title.position = "top",
                                         barwidth = 27,
                                         title.hjust = 0.5,
                                         title.theme = element_text(angle = 0,
                                                                    size = 20),
                                         label.theme = element_text(angle = 0,
                                                                    size = 20),
                                         ticks = T)) )
}
  dev.off()
}

#######################################################################
src.G0@meta.data$SampleName <- rownames(src.G0@meta.data)
src.G0@meta.data$cluster=as.character(src.G0@meta.data$RNA_snn_res.1.5)
src.G0@meta.data[src.G0@meta.data$cluster%in%c(1,14,5,6,2,3,12,10,8),]$cluster="beta"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(0),]$cluster="alpha"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(9),]$cluster="pp"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(7),]$cluster="delta"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(4),]$cluster="Endothelial cell"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(15,19,20,22),]$cluster="Immune cell"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(16,17),]$cluster="Mesenchymal cell"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(18),]$cluster="Duct cell"
src.G0@meta.data[src.G0@meta.data$cluster%in%c(21,11,23,13)|
                           src.G0@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet/multi-hormone cell"


src.G0@meta.data$cluster <- factor(src.G0@meta.data$cluster,
                                          levels = c('beta',
                                                     'alpha',
                                                     'delta',
                                                     'pp',
                                                     'Endothelial cell',
                                                     'Immune cell',
                                                     'Mesenchymal cell',
                                                     'Duct cell',
                                                     'Doublet/multi-hormone cell',
                                                     'Low quality cell'))


pp.meta.data <- MyReadDelim('G0_endo/pp/pp.gcgpp.meta.tab')
rownames(pp.meta.data) <- pp.meta.data$SampleName

sum(pp.meta.data[pp.meta.data$group == 'Gcg;Ppy','SampleName'] %in% src.G0$SampleName[src.G0$RNA_snn_res.3 == 21])

src.G0@meta.data[pp.meta.data[pp.meta.data$group == 'Gcg;Ppy','SampleName'],'cluster'] <- "Doublet/multi-hormone cell"

beta.low.quality.tab <- MyReadDelim('G0_endo/beta.low.quality.cell.tab')
src.G0@meta.data$cluster <- as.character(src.G0@meta.data$cluster)

src.G0@meta.data[beta.low.quality.tab$x,'cluster'] <- "Low quality cell"



p.tsne <- MySeuratv3TSNE10x2Gg(src.G0,
                               src.G0@meta.data
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

pdf("QC/G0.tsne.cellType.pc30.rmGcgPpy.rmbetalowquality.pdf",
    19,
    10)
p.plot <- p.age +
  scale_color_manual(values = time.colors) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()

src.G14.5 <- SetIdent(src.G14.5,value = src.G14.5$RNA_snn_res.1.5)
src.G14.5.res1.5.marker.res <- FindMarkers(src.G14.5,
                                                    ident.1 = 19,
                                                    ident.2 = c(1,2,5,6,8,9,3,11,14)
)

res0.15.19.marker <- src.G14.5.res1.5.marker.res[src.G14.5.res1.5.marker.res$avg_logFC > 1 &
                                                            src.G14.5.res1.5.marker.res$p_val_adj < 10^-5 ,]

rownames(res0.15.19.marker)

seurat <- 'src.G14.5'
pdf(paste('QC/',seurat,"res0.15.19.marker.pdf",sep = ""),6,7)
for (gene in rownames(res0.15.19.marker)) {
  print(FeaturePlot(get(seurat),gene,cols = c("#C7E8CC","#FF0000"))+
          theme(plot.title = element_text(size = rel(3.5))) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.title = element_blank()) +
          theme(axis.text  = element_blank()) +
          theme(axis.ticks = element_blank()) +
          theme(legend.position="bottom")+
          theme(panel.border = element_rect(size = 4,
                                            colour = "black"))+
          guides(colour = guide_colorbar(title = "ln(TP10K + 1)",
                                         title.position = "top",
                                         barwidth = 27,
                                         title.hjust = 0.5,
                                         title.theme = element_text(angle = 0,
                                                                    size = 20),
                                         label.theme = element_text(angle = 0,
                                                                    size = 20),
                                         ticks = T)) )
}
dev.off()

src.G14.5@meta.data$SampleName <-  rownames(src.G14.5@meta.data)

#src.G14.5@meta.data$SampleName <-  paste0("G14.5_20200416",rownames(src.G14.5@meta.data))
src.G14.5@meta.data$cluster=as.character(src.G14.5@meta.data$RNA_snn_res.1.5)
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(1,2,5,6,8,9,3,11,14),]$cluster="beta"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(0),]$cluster="alpha"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(12),]$cluster="pp"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(7),]$cluster="delta"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(4),]$cluster="Endothelial cell"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(15,20),]$cluster="Immune cell"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(13,17),]$cluster="Mesenchymal cell"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(21),]$cluster="Duct cell"
src.G14.5@meta.data[src.G14.5@meta.data$cluster%in%c(16,10,18,22,19,23)|
                               src.G14.5@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet/multi-hormone cell"


src.G14.5@meta.data$cluster <- factor(src.G14.5@meta.data$cluster,
                                               levels = c('beta',
                                                          'alpha',
                                                          'delta',
                                                          'pp',
                                                          'Endothelial cell',
                                                          'Immune cell',
                                                          'Mesenchymal cell',
                                                          'Duct cell',
                                                          'Doublet/multi-hormone cell'))

table(src.G14.5@meta.data$cluster)

src.G14.5@meta.data$cluster <- as.character(src.G14.5@meta.data$cluster)
G14.5.pp.gcg.cell <- MyReadDelim('G14.5_endo/rmcc/pp.Gcg.cell.tab')
src.G14.5@meta.data[G14.5.pp.gcg.cell$x,'cluster'] <- 'Doublet/multi-hormone cell'


p.tsne <- MySeuratv3TSNE10x2Gg(src.G14.5,
                               src.G14.5@meta.data
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

pdf("QC/rmppGcgG14.5.tsne.cellType.pc30.pdf",
    19,
    10)
p.plot <- p.age +
  scale_color_manual(values = time.colors) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()
##write out########
MyWriteTable(src.G0@meta.data,
             'QC/G0.celltype.new.meta.data')
saveRDS(src.G0,
        'QC/src.G0.rds')

MyWriteTable(src.G14.5@meta.data,
             'QC/G14.5.celltype.new.meta.data')
saveRDS(src.G14.5,
        'QC/src.G14.5.rds')

save.image('QC/pre.QC.RData')

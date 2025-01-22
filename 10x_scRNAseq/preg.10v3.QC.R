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
MyNorm.new=function(raw.data){
  size.factor <- colSums(raw.data)/mean(colSums(raw.data))
  norm.data <-  t(t(raw.data) / size.factor)
  return(norm.data)
}
genes.inf.input <- MyReadDelim("G:/lab/gene inf/mm10/mm10.gene.inf.merge.v1.1.tab")
rownames(genes.inf.input) <- genes.inf.input$SymbolDedu

time.colors <- c("#4DAF4A",
                 # "black",
                 # "#F01414",
                 # "#5F9EA0",#P21-2-B-con
                 "#FF7F00",
                 "#7570B3", 
                 "#452cff",#G5.5
                 "#A65628", 
                 "#F781BF",
                 "#999999",
                 "#ffe478",#P2-B
                 "#cab2d6",
                 "#ba984d",
                 "#1f78b4",
                 "#ff4a21",#P7-14-21,P21-B
                 "#752a78",
                 "#a6cee3",
                 "#53751c",
                 "#C71585",
                 "#F01414",
                 "#ff8e9d",
                 "#234682",#G14.5-HFD,
                 "#ff723f",#control culture for 4 days
                 "#d96ca8",#preg culture for 4 days
                 "#0f8575",#G18.5-KO
                 "#5277cf",#Nif-unpreg-Ctrl
                 "#de562c",#Nif-preg-Ctrl
                 "#5a7368",#Nif-preg
                 "#c77fc1",#C646-unpreg-Ctrl,
                 "#3fa0a3",#C646-preg-Ctrl
                 "#c7a96c",#C646-perg
                 "#805e2b",#human serum ctrl
                 "#F01414",
                 "#d15126",
                 #"#6bff9c",#G14.5-2nd
                 #"#5ce6ba"#,
                 # "#39a866",#human serum preg
                 "#6343c4"#human serum GDM
                 
                 
)
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

###################################################################
#                          load
###################################################################
load('rawdata/read.10x.CtrlG14.5_1stG14.5_2nd.RData')
src.duct.Ctrl=CreateSeuratObject(src.raw.Ctrl)
src.duct.G14.5_1st=CreateSeuratObject(src.raw.G14.5_1st)
src.duct.G14.5_2nd=CreateSeuratObject(src.raw.G14.5_2nd)
src.duct.Ctrl$Time="Ctrl"
src.duct.G14.5_1st$Time="G14.5_1st"
src.duct.G14.5_2nd$Time="G14.5_2nd"

rm(src.raw.Ctrl)
rm(src.raw.G14.5_1st)
rm(src.raw.G14.5_2nd)

########################################################################
#                            QC
########################################################################
src.duct.Ctrl[["percent.mt"]] <- PercentageFeatureSet(src.duct.Ctrl, pattern = "^mt-")
src.duct.G14.5_1st[["percent.mt"]] <- PercentageFeatureSet(src.duct.G14.5_1st, pattern = "^mt-")
src.duct.G14.5_2nd[["percent.mt"]] <- PercentageFeatureSet(src.duct.G14.5_2nd, pattern = "^mt-")


pdf('QC/features.pdf')
for(seu in c("src.duct.Ctrl",
             "src.duct.G14.5_1st",
             "src.duct.G14.5_2nd"
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


for (seurat in c("src.duct.Ctrl",
                 "src.duct.G14.5_1st",
                 "src.duct.G14.5_2nd"
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


for (seurat in c("src.duct.Ctrl",
                 "src.duct.G14.5_1st",
                 "src.duct.G14.5_2nd")){
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


VariableFeatures(src.duct.Ctrl) <- setdiff(VariableFeatures(src.duct.Ctrl),
                                            c(labels(src.duct.Ctrl.row.tree[[2]][[2]][[2]][[2]][[2]][[1]])))
VariableFeatures(src.duct.G14.5_1st) <- setdiff(VariableFeatures(src.duct.G14.5_1st),
                                            c(labels(src.duct.G14.5_1st.row.tree[[2]][[2]][[1]][[2]])))
VariableFeatures(src.duct.G14.5_2nd) <- setdiff(VariableFeatures(src.duct.G14.5_2nd),
                                            c(labels(src.duct.G14.5_2nd.row.tree[[2]][[2]][[2]][[2]][[1]])))

for (seurat in c(#"src.duct.Ctrl",
                 "src.duct.G14.5_1st",
                 "src.duct.G14.5_2nd"
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

for (seurat in c(#"src.duct.Ctrl"#,
                  "src.duct.G14.5_1st"#,
                 # "src.duct.G14.5_2nd"
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
src.duct.Ctrl@meta.data$SampleName <- rownames(src.duct.Ctrl@meta.data)
src.duct.Ctrl@meta.data$cluster=as.character(src.duct.Ctrl@meta.data$RNA_snn_res.1.5)
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(1,14,5,6,2,3,12,10,8),]$cluster="beta"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(0),]$cluster="alpha"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(9),]$cluster="pp"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(7),]$cluster="delta"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(4),]$cluster="Endothelial cell"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(15,19,20,22),]$cluster="Immune cell"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(16,17),]$cluster="Mesenchymal cell"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(18),]$cluster="Duct cell"
src.duct.Ctrl@meta.data[src.duct.Ctrl@meta.data$cluster%in%c(21,11,23,13)|
                           src.duct.Ctrl@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet/multi-hormone cell"


src.duct.Ctrl@meta.data$cluster <- factor(src.duct.Ctrl@meta.data$cluster,
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


pp.meta.data <- MyReadDelim('Ctrl_endo/pp/pp.gcgpp.meta.tab')
rownames(pp.meta.data) <- pp.meta.data$SampleName

sum(pp.meta.data[pp.meta.data$group == 'Gcg;Ppy','SampleName'] %in% src.duct.Ctrl$SampleName[src.duct.Ctrl$RNA_snn_res.3 == 21])

src.duct.Ctrl@meta.data[pp.meta.data[pp.meta.data$group == 'Gcg;Ppy','SampleName'],'cluster'] <- "Doublet/multi-hormone cell"

beta.low.quality.tab <- MyReadDelim('Ctrl_endo/beta.low.quality.cell.tab')
src.duct.Ctrl@meta.data$cluster <- as.character(src.duct.Ctrl@meta.data$cluster)

src.duct.Ctrl@meta.data[beta.low.quality.tab$x,'cluster'] <- "Low quality cell"



p.tsne <- MySeuratv3TSNE10x2Gg(src.duct.Ctrl,
                               src.duct.Ctrl@meta.data
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

pdf("QC/ctrl.tsne.cellType.pc30.rmGcgPpy.rmbetalowquality.pdf",
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

src.duct.G14.5_1st <- SetIdent(src.duct.G14.5_1st,value = src.duct.G14.5_1st$RNA_snn_res.1.5)
src.duct.G14.5_1st.res1.5.marker.res <- FindMarkers(src.duct.G14.5_1st,
                                                    ident.1 = 19,
                                                    ident.2 = c(1,2,5,6,8,9,3,11,14)
)

res0.15.19.marker <- src.duct.G14.5_1st.res1.5.marker.res[src.duct.G14.5_1st.res1.5.marker.res$avg_logFC > 1 &
                                                            src.duct.G14.5_1st.res1.5.marker.res$p_val_adj < 10^-5 ,]

rownames(res0.15.19.marker)

seurat <- 'src.duct.G14.5_1st'
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

src.duct.G14.5_1st@meta.data$SampleName <-  rownames(src.duct.G14.5_1st@meta.data)

#src.duct.G14.5_1st@meta.data$SampleName <-  paste0("G14.5_1st_20200416",rownames(src.duct.G14.5_1st@meta.data))
src.duct.G14.5_1st@meta.data$cluster=as.character(src.duct.G14.5_1st@meta.data$RNA_snn_res.1.5)
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(1,2,5,6,8,9,3,11,14),]$cluster="beta"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(0),]$cluster="alpha"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(12),]$cluster="pp"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(7),]$cluster="delta"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(4),]$cluster="Endothelial cell"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(15,20),]$cluster="Immune cell"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(13,17),]$cluster="Mesenchymal cell"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(21),]$cluster="Duct cell"
src.duct.G14.5_1st@meta.data[src.duct.G14.5_1st@meta.data$cluster%in%c(16,10,18,22,19,23)|
                               src.duct.G14.5_1st@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet/multi-hormone cell"


src.duct.G14.5_1st@meta.data$cluster <- factor(src.duct.G14.5_1st@meta.data$cluster,
                                               levels = c('beta',
                                                          'alpha',
                                                          'delta',
                                                          'pp',
                                                          'Endothelial cell',
                                                          'Immune cell',
                                                          'Mesenchymal cell',
                                                          'Duct cell',
                                                          'Doublet/multi-hormone cell'))

table(src.duct.G14.5_1st@meta.data$cluster)

src.duct.G14.5_1st@meta.data$cluster <- as.character(src.duct.G14.5_1st@meta.data$cluster)
G14.5_1st.pp.gcg.cell <- MyReadDelim('G14.5_1st_endo/rmcc/pp.Gcg.cell.tab')
src.duct.G14.5_1st@meta.data[G14.5_1st.pp.gcg.cell$x,'cluster'] <- 'Doublet/multi-hormone cell'


p.tsne <- MySeuratv3TSNE10x2Gg(src.duct.G14.5_1st,
                               src.duct.G14.5_1st@meta.data
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

pdf("QC/rmppGcgG14.5_1st.tsne.cellType.pc30.pdf",
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

src.duct.G14.5_2nd <- SetIdent(src.duct.G14.5_2nd,value = src.duct.G14.5_2nd$RNA_snn_res.1.5)
src.duct.G14.5_2nd.res1.5.marker.res <- FindMarkers(src.duct.G14.5_2nd,
                                                    ident.1 = 19,
                                                    ident.2 = c(1,2,5,6,8,9,3,11,14)
)

res0.15.19.marker <- src.duct.G14.5_2nd.res1.5.marker.res[src.duct.G14.5_2nd.res1.5.marker.res$avg_logFC > 1 &
                                                            src.duct.G14.5_2nd.res1.5.marker.res$p_val_adj < 10^-5 ,]

rownames(res0.15.19.marker)

seurat <- 'src.duct.G14.5_2nd'
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


src.duct.G14.5_2nd@meta.data$SampleName <-  paste0("G14.5_2nd_20200308",rownames(src.duct.G14.5_2nd@meta.data))
src.duct.G14.5_2nd@meta.data$cluster=as.character(src.duct.G14.5_2nd@meta.data$RNA_snn_res.3)
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(19,7,5,8,10,3,9,0,18,16,30,6,14,1,21,15,22),]$cluster="beta"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(17,12,13,25),]$cluster="alpha"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(11),]$cluster="pp"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(4),]$cluster="delta"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(2),]$cluster="Endothelial cell"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(24),]$cluster="Immune cell"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(29,23),]$cluster="Mesenchymal cell"
#src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(21),]$cluster="Duct cell"
src.duct.G14.5_2nd@meta.data[src.duct.G14.5_2nd@meta.data$cluster%in%c(20,27,28,32,31,26)|
                               src.duct.G14.5_2nd@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet/multi-hormone cell"


G14.5_2nd.doublet.tab <- MyReadDelim('G14.5_2nd_endo/G14.5_2nd.Immunedoublet.tab')

src.duct.G14.5_2nd@meta.data[G14.5_2nd.doublet.tab$x,'cluster'] <- 'Doublet/multi-hormone cell'



src.duct.G14.5_2nd@meta.data$cluster <- factor(src.duct.G14.5_2nd@meta.data$cluster,
                                               levels = c('beta',
                                                          'alpha',
                                                          'delta',
                                                          'pp',
                                                          'Endothelial cell',
                                                          'Immune cell',
                                                          'Mesenchymal cell',
                                                          # 'Duct cell',
                                                          'Doublet/multi-hormone cell'))

p.tsne <- MySeuratv3TSNE10x2Gg(src.duct.G14.5_2nd,
                               src.duct.G14.5_2nd@meta.data
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

pdf("QC/G14.5_2nd.tsne.cellType.pc30.new.pdf",
    19,
    10)
p.plot <- p.age +
  scale_color_manual(values = time.colors[c(1:7,9)]) +
  geom_point(aes(x = x.pos,
                 y = y.pos,
                 col =cluster#,
                 #shape = State
  ),
  size = 1
  )

print(p.plot)
dev.off()

MyWriteTable(src.duct.Ctrl@meta.data,
             'QC/ctrl.celltype.new.meta.data')
saveRDS(src.duct.Ctrl,
        'QC/src.duct.Ctrl.rds')

MyWriteTable(src.duct.G14.5_1st@meta.data,
             'QC/G14.5_1st.celltype.new.meta.data')
saveRDS(src.duct.G14.5_1st,
        'QC/src.duct.G14.5_1st.rds')

MyWriteTable(src.duct.G14.5_2nd@meta.data,
             'QC/G14.5_2nd.celltype.new.meta.data')
saveRDS(src.duct.G14.5_2nd,
        'QC/src.duct.G14.5_2nd.rds')

save.image('QC/pre.QC.RData')

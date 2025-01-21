setwd('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/beta/')
dir.create('proliferation')
setwd('proliferation/')
load('10x.beta.10x.cc.gene.RData')

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)

seu.beta.si.tab <- MyReadDelim('seu.beta.hc.knn.cc.si.tab')
table(seu.beta.si.tab$betagroup,seu.beta.si.tab$hc.knn.proliferative)

seu.beta.si.tab$betagroup.time <- paste(seu.beta.si.tab$betagroup,seu.beta.si.tab$Time,sep = '_')
table(seu.beta.si.tab$betagroup.time)
table(seu.beta.si.tab$betagroup.time,seu.beta.si.tab$hc.knn.proliferative)


setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/Glut2H/cluster/excludebatch/proliferation/')

library(future)
plan("multisession", workers = 10)
plan()

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)

preg.colors <- c("#4DAF4A",
                 "#F781BF"
)

names(preg.colors) <- c("Ctrl",
                        "G14.5_1st")
#########
seu.beta <- readRDS('../src.beta.rds')

beta.cc.sym2 <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/Glut2H/cluster/excludebatch/pre.cc.sym.rds')
beta.cc.sym <- readRDS('../beta.pre.cc.sym.rds')
beta.pre.cc.sym <- rownames(MyGeneExp(as.matrix(exp(seu.beta@assays$RNA@data)[beta.cc.sym,]),
                                      2,exp = 'yes',5))
length(beta.pre.cc.sym)

sum(beta.cc.sym2 %in% beta.pre.cc.sym)

beta.cc.merge <- c(beta.cc.sym2,beta.pre.cc.sym)
length(beta.cc.merge)#505

beta.pre.cc.sym2 <- rownames(MyGeneExp(as.matrix(exp(seu.beta@assays$RNA@data)[beta.cc.merge[beta.cc.merge %in% rownames(seu.beta)],]),
                                      2,exp = 'yes',5))
length(beta.pre.cc.sym2)#376

sum(beta.pre.cc.sym %in% beta.pre.cc.sym2)


c1Color <- MyName2Col(seu.beta@meta.data[,"Time"],
                      preg.colors)
c1Color <- as.matrix(c1Color)

png('vst.heatmap.merge.png',2000,3000)
cc.col.tree <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[beta.pre.cc.sym,],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = c1Color,
            ColSideColorsSize = 2,
            return.tree = "col",
            graph = T)
dev.off()
png('vst.heatmap.merge.png',2000,3000)
cc.col.tree2 <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[beta.pre.cc.sym2,],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = c1Color,
            ColSideColorsSize = 2,
            return.tree = "col",
            graph = T)
dev.off()

cc.col.tree.den <- as.dendrogram(cc.col.tree)
beta.cc.sn <- c(labels(cc.col.tree.den[[2]]),labels(cc.col.tree.den[[1]][[1]]),
                labels(cc.col.tree.den[[1]][[1]]),
                labels(cc.col.tree.den[[1]][[2]][[1]])
                )

cc.col.tree.den2 <- as.dendrogram(cc.col.tree2)
beta.cc.sn2 <- c(labels(cc.col.tree.den2[[1]]),labels(cc.col.tree.den2[[2]][[2]]))


seu.beta$hc.proliferative <- 'quiescent'
seu.beta@meta.data[beta.cc.sn,'hc.proliferative'] <- 'proliferative'
table(seu.beta$hc.proliferative,seu.beta$Time)

seu.beta$hc.proliferative2 <- 'quiescent'
seu.beta@meta.data[beta.cc.sn2,'hc.proliferative2'] <- 'proliferative'

seu.beta@meta.data['G14.5_1st_20200416GGTTAACTCCCAAGTA???1','hc.proliferative2'] <- 'proliferative'
table(seu.beta$hc.proliferative2,seu.beta$Time)

cc.type.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Time',c2 = 'hc.proliferative',cols = time.colors)
cc.group.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Time',c2 = 'knn.proliferative',cols = preg.colors)


pdf('proliferative/cc.bar.pdf',8,6)
barplot(cc.type.ratio[1,])
barplot(cc.group.ratio[1,])
dev.off()

table(seu.beta$Type_rep,seu.beta$hc.proliferative)
hc.knn.cc.sn <- unique(c(colnames(seu.beta)[seu.beta$knn.proliferative=='proliferative'],
                         colnames(seu.beta)[seu.beta$hc.proliferative=='proliferative']))
length(hc.knn.cc.sn)#227
hc.knn.cc.sn <- na.omit(hc.knn.cc.sn)

seu.beta$hc.knn.proliferative <- 'quiescent'
seu.beta@meta.data[hc.knn.cc.sn,'hc.knn.proliferative'] <- 'proliferative'
cc.type.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Time',c2 = 'hc.knn.proliferative',cols = time.colors)
barplot(cc.type.ratio[1,])

table(seu.beta$knn.proliferative,seu.beta$hc.proliferative)

###########
seu.beta <- RunPCA(seu.beta,features = beta.pre.cc.sym)
Myseuratmarker(seu.beta,marker.sym = c('Mki67','Cdt1','Gmnn'),reduction = 'pca')
setwd('G:/lab/Article/heshuang/BYLW/10x/')
pdf('knn.hc.pc12.proliferative.pdf',13,12)
p.pca <- MySeuratDR2Gg2(seu.beta,seu.beta@meta.data,
                        reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
p.pca$data_ <- p.pca$data
p.pca$data_$hc.proliferative_ <- p.pca$data_$Time
p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(seu.beta$Time))))
p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]


plot(p.pca+
       scale_color_manual(values = c('cyan3','red3')) +
       scale_shape_manual(values = c(17,19)) +
       geom_point(aes(x =-x.pos,
                      y = y.pos,
                      col = Time,
                      shape = hc.knn.proliferative
       ),
       size = 6
       )+
       theme(aspect.ratio=1))

# p.pca$data_ <- p.pca$data
# p.pca$data_$hc.proliferative_ <- p.pca$data_$knn.hc.proliferative
# p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(seu.beta$knn.hc.proliferative))))
# p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]

# plot(p.pca+
#        scale_color_manual(values = c("#ba984d",'gray90')) +
#        geom_point(aes(x =-x.pos,
#                       y = y.pos,
#                       col = hc.proliferative#,
#                       # shape = State
#        ),
#        size =6
#        )+
#        theme(aspect.ratio=1))

# plot(p.pca+
#        
#        scale_color_manual(values = c("#ba984d",'gray90')) +
#        geom_point(aes(x =-x.pos,
#                       y = y.pos,
#                       col = hc.knn.proliferative#,
#                       # shape = State
#        ),
#        size =6
#        )+
#        theme(aspect.ratio=1))
# plot(p.pca+
#        scale_color_manual(values = c('orange','gray90')) +
#        geom_text(aes(x =-x.pos,
#                       y = y.pos,
#                       col = hc.proliferative#,
#                       # shape = State
#        ),
#        size = 0.5
#        ))
# plot(p.pca+
#        scale_color_manual(values = c('orange','gray90','black')) +
#        geom_point(aes(x =x.pos,
#                       y = y.pos,
#                       col = cc.old#,
#                       # shape = State
#        ),
#        size = 4
#        ))
dev.off()
##########
p.pca <- MySeuratDR2Gg2(seu.beta,seu.beta@meta.data,
                        reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
p.pca$data_ <- p.pca$data
p.pca$data_$hc.proliferative_ <- p.pca$data_$Time
p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(seu.beta$Time))))
p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]

p.loop <- p.pca +
  scale_color_gradientn(colours = colors.exp) +
  #scale_shape_manual(values = 1:5)
  theme(plot.title = element_text(size = rel(3.5))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_blank()) +
  theme(axis.text  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(legend.position="bottom")+
  theme(panel.border = element_rect(size = 4,
                                    colour = "black"))+
  guides(colour = guide_colorbar(title = "ln(TP0.1M + 1)",
                                 title.position = "top",
                                 barwidth = 20,
                                 title.hjust = 0.5,
                                 title.theme = element_text(angle = 0,
                                                            size = 20),
                                 label.theme = element_text(angle = 0,
                                                            size = 20),
                                 ticks = T)) 
cellcycle.marker.set=c("Cdk2","Ccne1","Mcm2","E2f7",      #G1/S
                       "Cdk1","Ccnb1","Kif4",             #G2/M
                       "Cdkn1c","Cdkn2c","Cdkn1b",
                       'Cdt1','Gmnn','Top2a','Mki67')   #p57(Cdkn1c)p18(Cdkn2c)p16(Cdkn2a)p27(Cdkn1b) inhibitor

pdf("hc.knn.cc.marker.cellcycle.pdf",
    6,
    7)
for (plot.ens in cellcycle.marker.set){
  print(plot.ens)
  print(p.plot <- p.loop +
          geom_point(aes(x = -x.pos,
                         y = y.pos,
                         
                         color = seu.beta@assays$RNA@data[plot.ens,rownames(p.pca$data)]
                         #shape = Age
          ),size = 5) +
          labs(title = plot.ens)+
          theme(aspect.ratio=1))
  
}

dev.off()

MyWriteTable(seu.beta@meta.data,'seu.beta.hc.knn.cc.si.tab')
#MyWriteTable(table(seu.beta$Type,seu.beta$hc.proliferative),row.names = T,'seu.beta.type.cc.count.tab')
MyWriteTable(table(seu.beta$Time,seu.beta$hc.knn.proliferative),row.names = T,'10x.hc.knn.seu.beta.group.cc.count.tab')

##########
pc.use <- 1:5
seu.beta <- FindNeighbors(seu.beta, dims = pc.use)
seu.beta <- FindClusters(seu.beta, resolution = 1.5)
seu.beta <- RunUMAP(seu.beta,dims = pc.use)

seu.beta <- RunTSNE(seu.beta,
                    dims = pc.use,
                    perplexity= round((30+ncol(seu.beta)/100)),
                    check_duplicates = F)

pdf('cc.tsne.marker.res1.5.pdf')
seu.beta <- SetIdent(seu.beta,value = seu.beta$hc.proliferative)
DimPlot(seu.beta,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = F)

seu.beta <- SetIdent(seu.beta,value = seu.beta$RNA_snn_res.1.5)
DimPlot(seu.beta,
        reduction = "tsne",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 2,
        label = T)

dev.off()


pdf('knn.tsne.pdf')
Myseuratmarker(seu.beta,marker.sym = c('Mki67','Gmnn','Cdt1','Top2a','Mcm2'),reduction = 'tsne')
dev.off()

p.umap <- MySeuratv3TSNE10x2Gg(seu.beta,
                               seu.beta@meta.data
)

pdf("FS1.umap.subgroup.pdf",
    14.5,
    8)
print(p.umap +
        scale_color_manual(values = preg.colors) +
        geom_point(aes(x = x.pos,
                       y = y.pos,
                       col =Time#,
                       #shape = State
        ),
        size = 3
        )+ theme(aspect.ratio=1))

dev.off()
# 
seu.beta$knn.proliferative <- 'quiescent'
seu.beta@meta.data[seu.beta$RNA_snn_res.1.5 %in% c(17,18),'knn.proliferative'] <- 'proliferative'
# 
# 
cc.type.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Time',c2 = 'knn.proliferative',cols = time.colors)

pdf('proliferative/cc.bar.pdf',8,6)
barplot(cc.type.ratio[1,])
dev.off()
# 
# 
# MyWriteTable(table(seu.beta@meta.data$Type,seu.beta@meta.data$hc.proliferative3),row.names = T,'proliferative/type.cc.count.tab')
# MyWriteTable(table(seu.beta@meta.data$Type_rep,seu.beta@meta.data$hc.proliferative3),row.names = T,'proliferative/type_rep.cc.count.tab')
# 
# 
# MyWriteTable(seu.beta@meta.data,'proliferative/seu.beta.cc.meta.tab')
# MyWriteTable(table(seu.beta$Type_rep,seu.beta$hc.proliferative),row.names = T,'proliferative/seu.beta.cc.meta.tab')

###########pc2 genes#########
PCA.supp <- FactoMineR::PCA(t(as.matrix(seu.beta@assays$RNA@data[beta.cc.sym,])),
                            graph = T
)
library(FactoMineR)
beta.PC12.dim.res <- dimdesc(PCA.supp,
                             axes = c(1,2),
                             proba = 0.01
)
beta.pc12.dim2 <- as.data.frame(beta.PC12.dim.res$Dim.2$quanti)
cc.pc2.genes <- rownames(beta.pc12.dim2)[-log10(beta.pc12.dim2$p.value)>=5]

length(cc.pc2.genes)#205

beta.pc12.dim1 <- as.data.frame(beta.PC12.dim.res$Dim.1$quanti)
cc.pc1.genes <- rownames(beta.pc12.dim1)[-log10(beta.pc12.dim1$p.value)>=5]
length(cc.pc1.genes)#313

beta.pc12.cc.gene <- unique(c(cc.pc2.genes,cc.pc1.genes))
length(beta.pc12.cc.gene)#255

# c1Color <- MyName2Col(seu.beta@meta.data[,"Type"],
#                       time.colors)
# c1Color <- as.matrix(c1Color)
# c2Color <- MyName2Col(seu.beta@meta.data[,"hc.proliferative"],
#                       c('orange','gray90'))
# c2Color <- as.matrix(c2Color)
# 
# cColor <- cbind(c1Color,c2Color)
# 
# png('proliferative/pc12.filter.cc.heatmap.png',2000,3000)
# cc.col.pc12.tree <-
#   MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[beta.pc12.cc.gene,],
#             type = "log.row.relat",
#             hc.c.data.type = "log.row.relat",
#             hc.r.data.type = "log.row.relat",
#             c.cov.method = "s",
#             r.cov.method = "s",
#             c.hc.method = "ward.D2",
#             r.hc.method = "ward.D2",
#             ColSideColors = cColor,
#             ColSideColorsSize = 2,
#             return.tree = "col",
#             graph = T)
# dev.off()
# 
# cc.col.pc12.tree.den <- as.dendrogram(cc.col.pc12.tree)
# beta.pc12.cc.sn <- c(labels(cc.col.pc12.tree.den[[2]]),labels(cc.col.pc12.tree.den[[1]][[1]]),
#                      labels(cc.col.pc12.tree.den[[1]][[2]][[1]][[1]]))
# 
# seu.beta$pc12.hc.proliferative <- 'quiescent'
# seu.beta@meta.data[beta.pc12.cc.sn,'pc12.hc.proliferative'] <- 'proliferative'
# table(seu.beta$Type,seu.beta$pc12.hc.proliferative)
# table(seu.beta$Type_rep,seu.beta$pc12.hc.proliferative)
# 
# 
# cc.type.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'pc12.hc.proliferative',cols = time.colors)
# cc.typerep.ratio <- Mybarplot(seu.beta@meta.data,c1 = 'Type_rep',c2 = 'pc12.hc.proliferative',cols = time.colors)
# pdf('proliferative/cc.bar.pdf',8,6)
# barplot(cc.type.ratio[1,])
# barplot(cc.typerep.ratio[1,])
# dev.off()
########
rm(seu.beta)
save.image('10x.beta.10x.cc.gene.RData')

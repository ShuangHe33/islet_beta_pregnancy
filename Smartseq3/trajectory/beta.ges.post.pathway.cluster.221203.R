
dir.create('cluster')
setwd('cluster')


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
time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")
celltype.col2 <- time.col[c(2,9,17,12)]
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2')

ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)

length(pre.ambigous.sym)#198

pre.mt.gene <- Virgin.LH.DEG.filter$gene
pre.ambigous.sym <- c(pre.ambigous.sym,pre.mt.gene)
length(pre.ambigous.sym)#227

seu.beta$beta2 <- factor(seu.beta$beta2,levels = c('Glut2H_2','Glut2H_1'))
seu.beta$beta2 <- factor(seu.beta$beta2,levels = c('Glut2L_2','Glut2L_1'))
celltype.col2 <- c("#A65628","#5F9EA0", "#A65628")
names(celltype.col2) <- c('Glut2L','Glut2H_1','Glut2H_2')

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/Fs2j.Glut2H.subgroup.barplot.pdf',8,8)
Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'beta2',xlim = 20,cols = celltype.col2[c('Glut2H_2','Glut2H_1')])
dev.off()

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/Fs2j.Glut2L.subgroup.barplot.pdf',8,8)
Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'beta2',xlim = 20,cols = celltype.col2[c('Glut2L_2','Glut2L_1')])
dev.off()

pre.ambigous.sym <- c(pre.ambigous.sym,Virgin.LH.DEG.filter$gene)
pre.ambigous.sym <- setdiff(pre.ambigous.sym,'Chgb')
pre.ambigous.sym <- readRDS('pre.ambigous.sym.mt.rds')
dir.create('rmG0batch1')
setwd('rmG0batch1/')
seu.beta <- readRDS('seu.beta.ges.post.rds')
seu.beta$Type <- as.character(seu.beta$Type)

seu.beta$Type <- factor(seu.beta$Type,levels = names(ref.time.colors))
seu.beta <- readRDS('seu.beta.ges.post.rds')

seu.beta <- merge(seu.beta.new,seu.beta)
seu.beta$Type <- factor(seu.beta$Type,levels = c(names(ref.time.colors)[1:10],'P0',names(ref.time.colors)[11:13],'P1NL',names(ref.time.colors)[14:16]))
table(seu.beta$Type)
#######QC#####
pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/QC.genecount.box.pdf',4,5)
boxplot(seu.beta$Genecount_all,lwd=2,main='gene number')
box(lwd=3)
boxplot(seu.beta$Raw_count,lwd=2,main='read count')
box(lwd=3)
dev.off()
############
seu.beta <- FindVariableFeatures(seu.beta,selection.method = "vst", nfeatures = 2000)
VariableFeatures(seu.beta) <- setdiff(VariableFeatures(seu.beta),pre.ambigous.sym)#

c1Color <- MyName2Col(seu.beta@meta.data[,"Type"],
                      ref.time.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.beta@meta.data[,"Rep"],
                      rep.colors[-1])
c2Color <- as.matrix(c2Color)
c3Color <- MyName2Col(seu.beta@meta.data[,"beta2"],
                      celltype.col2)
c3Color <- as.matrix(c3Color)

cColor <- cbind(c1Color,c2Color,c3Color)

png('Vstvar.2000.geneover6000.raw.png',2000,3000)
var2000.row.tree <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[VariableFeatures(seu.beta),colnames(seu.beta)],
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

var2000.row.tree.den <- as.dendrogram(var2000.row.tree)
grep('Mt2',VariableFeatures(seu.beta))

pre.cc.sym <- c(labels(as.dendrogram(var2000.row.tree)[[1]]))
length(pre.cc.sym)#257

VariableFeatures(seu.beta) <- setdiff(VariableFeatures(seu.beta),c(pre.ambigous.sym,pre.cc.sym))
saveRDS(pre.cc.sym,'pre.cc.sym.rds')
###############
all.var.co2.cc <- rownames(MyCo(as.matrix(seu.beta@assays$RNA@data),
                                var.gene = c(VariableFeatures(seu.beta),pre.cc.sym),
                                exp.prop.whole.max = 1,
                                exp.prop.whole.min = 0.01,
                                # vector.group = samples.inf.qc$GroupNameFig1,
                                # exp.prop.group.min = 0.1,
                                # exp.prop.group.max = 0.5,
                                cor.method = "rho",
                                cor.cutoff = 0.2,
                                partner.cutoff = 10,
                                refine.cor = T))
length(all.var.co2.cc)#cor0.2 254

saveRDS(all.var.co2.cc,'all.var.co2.cc.final.rds')

all.var.co2 <- rownames(MyCo(as.matrix(seu.beta@assays$RNA@data),
                                 var.gene = VariableFeatures(seu.beta),
                                 exp.prop.whole.max = 1,
                                 exp.prop.whole.min = 0.01,
                                 # vector.group = samples.inf.qc$GroupNameFig1,
                                 # exp.prop.group.min = 0.1,
                                 # exp.prop.group.max = 0.5,
                                 cor.method = "rho",
                                 cor.cutoff = 0.2,
                                 partner.cutoff = 10,
                                 refine.cor = T))
length(all.var.co2)#cor0.2 82


saveRDS(all.var.co2,'all.var.co2.10.0.01.rmcc.pathway.rds')
saveRDS(all.var.co2.cc,'all.var.co2.10.0.01.cc.pathway.rds')

png('vst2000.rho.cor0.2.1.0.01.10.rho.png',2000,3000)
var.co0.2.row.tree <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[all.var.co2,],
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

beta.heter.sym <- labels(as.dendrogram(var.co0.2.row.tree)[[1]][[2]])

#var.co0.2.row.tree.den <- as.dendrogram(var.co0.2.row.tree)
png('cc.vst2000.rho.cor0.2.1.0.01.10.rho.png',2000,3000)
cc.var.co0.2.row.tree <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[all.var.co2.cc,],
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

grep('Mt1',all.var.co2.cc)
##########
seu.beta <- readRDS('seu.beta.ges.post.rds')

seu.beta <- ScaleData(seu.beta,features = c(VariableFeatures(seu.beta),pre.cc.sym))
#seu.beta <- RunPCA(seu.beta, features = vst.gene.input)
#seu.beta <- RunPCA(seu.beta, features = all.var.co)
seu.beta <- RunPCA(seu.beta, features = all.var.co2)

seu.beta <- RunPCA(seu.beta, features = all.var.co2.cc)

# pdf('rmcc.cor0.2.1.DimHeatmap.pdf',
#     10,40)
# DimHeatmap(seu.beta, dims = 1:30, cells = 200, balanced = TRUE)
# dev.off()
# 
# #pc.use <- 1:13;#vst
pc.use <- 1:15;#cor0.2
seu.beta <- FindNeighbors(seu.beta, dims = pc.use)
#seu.beta <- FindClusters(seu.beta, resolution = 1.5)
seu.beta <- FindClusters(seu.beta, resolution = 0.2)
#seu.beta <- FindClusters(seu.beta, resolution = 3.5)
seu.beta <- RunUMAP(seu.beta,dims = pc.use)

#saveRDS(all.var.co,'all.var.co0.15.1.0.01.10.pc10.beforefilterhetergene.rds')
Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'RNA_snn_res.0.3',xlim = 20,cols = time.colors)
pdf('umap.rmcc.vst.cor0.15.1.0.01.10pc15.res0.2.pdf',
    6,5)
seu.beta <- SetIdent(seu.beta,value = seu.beta$RNA_snn_res.0.2)
DimPlot(seu.beta,
        reduction = "umap",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 1,
        label = T)
seu.beta <- SetIdent(seu.beta,value = seu.beta$subgroup)
DimPlot(seu.beta,
        reduction = "umap",
        cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 1,
        label = F)

seu.beta <- SetIdent(seu.beta,value = seu.beta$Type)
DimPlot(seu.beta,
        reduction = "umap",
        cols = c(ref.time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
        label.size = 6,
        sizes.highlight = 4,
        pt.size = 1,
        label = F)
dev.off()

seu.beta$RNA_snn_res.0.2 <- factor(seu.beta$RNA_snn_res.0.2,levels = c('1','2','0'))

seu.beta$sub.time <- paste(seu.beta$Type,seu.beta$RNA_snn_res.0.2,sep = '_')

seu.beta$subgroup <- 'Glut2H_1'
seu.beta@meta.data[seu.beta$RNA_snn_res.0.2==2,'subgroup'] <- 'Glut2H_2_1'
seu.beta@meta.data[seu.beta$RNA_snn_res.0.2==0,'subgroup'] <- 'Glut2H_2_2'

MyWriteTable(table(seu.beta$Type,seu.beta$subgroup),'subgroup.count.tab')

seu.meta.tab <- MyReadDelim('../../cluster/proliferation/seu.beta.meta.tab')
rownames(seu.meta.tab) <- seu.meta.tab$SampleName
seu.beta@meta.data$proliferation <- seu.meta.tab[seu.beta$SampleName,]$hc.proliferative

MyWriteTable(seu.beta@meta.data,'seu.beta.meta.data.tab')

MyWriteTable(table(seu.beta$Type,seu.beta$subgroup,seu.beta$proliferation),'subgroup.cc.count.tab')


seu.beta$sub.time <- factor(seu.beta$sub.time,levels = paste(rep(names(ref.time.colors),each=3),c('1','2','0'),sep = '_'))

pdf('vio.Acss2.sub.split.pdf',18,8)
Myseuvioplot(seu.beta,sym.list = c('Acss2','Oxtr'),type = 'sub.time',type.colors = rep(ref.time.colors,each=3))
dev.off()


pdf('vio.Acss2.sub.pdf',6,8)
Myseuvioplot(seu.beta,sym.list = c('Acss2','Oxtr'),type = 'RNA_snn_res.0.2')
dev.off()

tmp.p <- wilcox.test(seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$RNA_snn_res.0.2=='2']],
                     seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$RNA_snn_res.0.2=='0']]
                     )
tmp.p$p.value


tmp.p <- wilcox.test(seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G12.5_2']],
                     seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G12.5_0']]
)
tmp.p$p.value#0.1

tmp.p <- wilcox.test(seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G14.5_2']],
                     seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G14.5_0']]
)
tmp.p$p.value#0.2

tmp.p <- wilcox.test(seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G16.5_2']],
                     seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G16.5_0']]
)
tmp.p$p.value#0.2

tmp.p <- wilcox.test(seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G18.5_2']],
                     seu.beta@assays$RNA@data['Acss2',colnames(seu.beta)[seu.beta$sub.time=='G18.5_0']]
)
tmp.p$p.value#0.44
###########marker#######
seu.beta <- SetIdent(seu.beta,value = seu.beta$RNA_snn_res.0.2)
glut2H.sub.DEG <- Myseufindmarker(seu.beta,gene.include = gene.include,ident.1 = '2',ident.2 = '0',c1 = 'Glut2H2_1',c2 = 'Glut2H2_2')
Myseuratmarker(seu.beta,marker.sym = c('Oxtr','Ovol2','Acss2','Gbp8','Ikzf4'),reduction = 'umap',pt.size= 1)
glut2H.sub.DEG.filter <- 


pdf("hetermaketr.umap.log.marker.rmcc.vst.cor0.15.1.0.01.10.pc10.pdf",
    6,
    7)
Myseuratmarker(seu.beta,marker.sym = 'Mt1',reduction = 'umap',pt.size= 1)
dev.off()
Myseuratmarker(seu.beta,marker.sym = c('Oxtr','Ovol2','Acss2','Gbp8','Ikzf4'),reduction = 'umap',pt.size= 1)

#########pseudotime######
seu.beta <- readRDS('seu.beta.ges.post.rds')
seu.beta <- RunPCA(seu.beta, features = all.var.co2)
seu.beta <- RunPCA(seu.beta, features = all.var.co2.cc)

seu.beta$pc1 <- Embeddings(seu.beta,'pca')[,1]
seu.beta$pseudotime <- seu.beta$pc1

pdf(paste('pc1.pseudobox.pdf',sep = ''),10,6)
print(MyPseudotimebox(seu.beta@meta.data,time.colors = ref.time.colors,size.point = 2,box.lwd = 1.5))
dev.off()


pdf('.pc12.pdf',15,12)
p.pca <- MySeuratDR2Gg2(seu.beta,seu.beta@meta.data,reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
p.pca$data_ <- p.pca$data
p.pca$data_$Type_ <- p.pca$data_$Type
p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = rev(names(table(seu.beta$Type))))
p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ref.time.colors2) +
       geom_point(aes(x = x.pos,
                      y = y.pos,
                      col = Type#,
                      # shape = State
       ),
       size = 4
       )+theme(aspect.ratio=1))


dev.off()
pdf('Figs2k.ref.rmcccc.cor0.2.1.pc12.2.pdf',15,12)
p.pca <- MySeuratDR2Gg2(seu.beta,seu.beta@meta.data,reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
p.pca$data_ <- p.pca$data
p.pca$data_$Type_ <- p.pca$data_$Type
p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = rev(names(table(seu.beta$Type))))
p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ref.time.colors) +
       geom_point(aes(x =x.pos,
                      y = y.pos,
                      col = Type
                      # shape = State
       ),
       size = 3
       )+theme(aspect.ratio=1))



dev.off()

si.select <- seu.beta@meta.data
pregScore=c()
for (time in names(table(seu.beta$Type))) {
  regress=lm(-Embeddings(seu.beta,'pca')[si.select$Type==time,2]~
               Embeddings(seu.beta,'pca')[si.select$Type==time,1])$coefficients[2]
  temp=-c(1,regress)%*%rbind(Embeddings(seu.beta,'pca')[si.select$Type==time,1],
                             Embeddings(seu.beta,'pca')[si.select$Type==time,2])
  # temp=temp-mean(temp)
  names(temp)=si.select[si.select$Type==time,]$SampleName
  pregScore=c(pregScore,temp)
}

seu.beta@meta.data[names(pregScore),'pregScore'] <- pregScore
seu.beta$pseudotime <- seu.beta$pregScore

#############cc 124 24 preg########
si.select <- seu.beta@meta.data
pregScore=c()
for (time in names(table(seu.beta$Type))) {
  regress=lm(-Embeddings(seu.beta,'pca')[si.select$Type==time,4]~
               Embeddings(seu.beta,'pca')[si.select$Type==time,2])$coefficients[2]
  temp=-c(1,regress)%*%rbind(Embeddings(seu.beta,'pca')[si.select$Type==time,2],
                             Embeddings(seu.beta,'pca')[si.select$Type==time,4])
  # temp=temp-mean(temp)
  names(temp)=si.select[si.select$Type==time,]$SampleName
  pregScore=c(pregScore,temp)
}

seu.beta@meta.data[names(pregScore),'pregScore'] <- pregScore
seu.beta$pseudotime <- -seu.beta$pregScore

########
seu.beta <- readRDS('seu.beta.ges.post.rds')
pdf('PCA.cor0.1.1.10,0.01.preg.merge.loess.final.pdf',12,8)
#MyPseudotimebox2(seu.beta@meta.data,time.colors = group.col)
MyPseudotimebox(seu.beta@meta.data,time.colors = ref.time.colors)
MyPseudotimebox(seu.beta@meta.data,time.colors = ref.time.colors[c(1,1,1,1,rep(2:6,each=2),7,7,7,8,8,9,9,9,rep(10:11,each=2),12,12,12,rep(13:16,each=2))],rep = T)
dev.off()

seu.beta.si.tab <- seu.beta@meta.data
seu.beta.si.tab$Type <- factor(seu.beta.si.tab$Type,levels = c('Virgin','G3.5','G5.5','G6.5','G8.5','G10.5','G12.5','G14.5','G16.5','G18.5','Gap2',
                                                               'P2L','P7L','P14L','P2NL','P7NL','P14NL'))

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/Fig1c.centergravity.PCA.preg.merge.loess.3.pdf',15,8)
#MyPseudotimebox2(seu.beta@meta.data,time.colors = group.col)
print(MyPseudotimebox(seu.beta.si.tab,time.colors = ref.time.colors2,box.lwd = 1.8,size.point = 3))
dev.off()

saveRDS(seu.beta,'seu.beta.ges.post.rds')
##########cell Type#######
group.col <- time.colors[c(8,5,7)]
names(group.col) <- c('beta1','beta2','beta3')
seu.beta <- readRDS('seu.beta.ges.post.rds')
seu.beta@meta.data$heter.group.merge <- factor(seu.beta@meta.data$heter.group.merge,levels = c('beta3','beta2','beta1'))

pdf('Glut2H.subgroup.barplot.pdf',8,8)
Mybarplot(seu.beta@meta.data,c1 = 'Type',c2 = 'heter.group.merge',xlim = 20,cols = group.col[c(3,2,1)])
dev.off()

##############
p.umap <- MySeuratv3UMAP10x2Gg(seu.beta,
                               seu.beta@meta.data
)

pdf("F1B.umap.subgroup.pdf",
    14.5,
    8)
print(p.umap +
        scale_color_manual(values = group.col) +
        geom_point(aes(x = x.pos,
                       y = y.pos,
                       col =heter.group.merge#,
                       #shape = State
        ),
        size = 2
        )+ theme(aspect.ratio=1))

dev.off()

pdf("F1B.umap.Type.final.pdf",
    14.5,
    8)
print(p.umap +
        scale_color_manual(values = ref.time.colors) +
        geom_point(aes(x = x.pos,
                       y = y.pos,
                       col =Type#,
                       #shape = State
        ),
        size = 2
        )+theme(aspect.ratio=1))

dev.off()

#########rep#########
sample.inf.keep.cp <- seu.beta@meta.data
age.list <- names(table(sample.inf.keep.cp$Type))

pdf("G:/lab/Article/heshuang/BYLW/sm3/ref/rep.batch.age.pdf",
    13,
    11)
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
  
  
  p.tmp <- MySeuratv3TSNE10x2Gg(seu.beta,
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
          scale_color_manual(values = c("gray90",'skyblue1','sienna1','deeppink1','purple1'))+
          geom_point(aes(x = x.pos,
                         y = y.pos,
                         col = Rep1#,
                         #shape = State
          ),
          size = 5
          )+
          labs(title = i))
}
dev.off()
########heatmap######
seu.beta <- readRDS('seu.beta.ges.post.rds')
seu.beta$Type2 <- '/'
seu.beta@meta.data[seu.beta$Type %in% c('P2L','P7L','P14L','P2NL','P7NL','P14NL'),'Type2'] <- as.character(seu.beta@meta.data[seu.beta$Type %in% c('P2L','P7L','P14L','P2NL','P7NL','P14NL'),'Type'])
seu.beta$Type2 <- factor(seu.beta$Type2,levels = c('/','P2L','P7L','P14L','P2NL','P7NL','P14NL'))


seu.beta$Type3 <- 'ges'
seu.beta@meta.data[seu.beta$Type %in% c('P2L','P7L','P14L'),'Type3'] <- 'L'
seu.beta@meta.data[seu.beta$Type %in% c('P2NL','P7NL','P14NL'),'Type3'] <- 'NL'
seu.beta$Type3 <- factor(seu.beta$Type3,levels = c('ges','L','NL'))


preg.gene <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/ges/LMgene/gene_inf/preg.gene.si.1e20.tab')
table(preg.gene$cluster)
preg.1.sym <- sub('_','-',preg.gene$SymbolDedu[preg.gene$cluster=='preg.up'])
preg.2.sym <- sub('_','-',preg.gene$SymbolDedu[preg.gene$cluster=='preg.down'])
beta.cc.sym <- sub('_','-',preg.gene$SymbolDedu[preg.gene$cluster=='cc'])

c1Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"Type2"],
                      c('gray90',ref.time.colors[10:15]))
c1Color <- as.matrix(c1Color)
# c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pregScore)],"Type2"],
#                       c('gray90',ref.time.colors[c('P2L','P7L','P14L','P2NL','P7NL','P14NL')]))
# c2Color <- as.matrix(c2Color)

c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"Type3"],
                      c('gray90',time.colors[1:2]))
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)

preg.1.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[preg.1.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                          graph = T,
                                          preg.1.sym)
preg.2.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[preg.2.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                          graph = T,
                                          preg.2.sym)
beta.cc.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[beta.cc.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                           graph = T,
                                           beta.cc.sym)

seu.beta <- ScaleData(seu.beta,features = c(beta.cc.order,preg.1.order,preg.2.order))

png("FigS2g.postpartum.preg.gene.scaledata.final.png",
    2000,2500)
# row.tree.count <- 2
r1Color <- c(rep(gene.tree.colors[1],length(beta.cc.order)),
             rep(gene.tree.colors[2],length(preg.1.order)),
             rep(gene.tree.colors[3],length(preg.2.order))
)
r1Color <- as.matrix(t(r1Color))

MyHeatmap(as.matrix(as.matrix(exp(seu.beta@assays$RNA@scale.data))[c(beta.cc.order,
                                                                     preg.1.order,
                                                                     preg.2.order),colnames(seu.beta)[order(seu.beta$pseudotime)]]),
          type = "log.row.relat",
          #hc.c.data.type = "log.row.relat",
          hc.r.data.type = "log.row.relat",
          color.palette = pc12.heatmap.col,
          #c.cov.method = "s",
          # r.cov.method = "s",
          Colv = "none",
          Rowv = "none",
          dendrogram = "none",
          #c.hc.method = "ward.D2",
          #r.hc.method = "ward.D2",
          RowSideColors = r1Color,
          ColSideColors = cColor,
          ColSideColorsSize = 2#,
          #return.tree = "none"
)
dev.off()
########
MyWriteTable(table(seu.beta$Type_rep,seu.beta$heter.group.merge),row.names = T,
             'celltype.Glut2H.rep.count.tab')

MyWriteTable(table(seu.beta$Type),row.names = T,
             'celltype.Glut2H.ges.post.type.count.tab')
MyWriteTable(seu.beta@meta.data,'seu.ref.beta.meta.tab')
##########
library(rgl)
library(plotly)
#par3d()
plot3d(x=seu.beta@reductions$pca@cell.embeddings[,1],
       y=seu.beta@reductions$pca@cell.embeddings[,2],
       z=seu.beta@reductions$pca@cell.embeddings[,4],
       size=12,
       #sizes = c(300,200),
       xlab = "PC1", ylab = "PC2", zlab = "PC4",
       box=T,lwd = 1,
       axes = T,
       col=MyName2Col(seu.beta$Type,ref.time.colors))

########
saveRDS(seu.beta,'seu.beta.ges.post.rds')
rm(seu.beta,seu.beta.V)
save.image('seu.beta.ges.post.RData')


library(future)
plan("multisession", workers = 10)
plan()

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
#############cc genes#######
seu.ref.L <- readRDS('Glut2L/seu.ref.L.endo.celltype.rds')
gene.input <- rownames(MyGeneExp(as.matrix(seu.ref.L@assays$RNA@data),log(2),10))
gene.input <- setdiff(gene.input,c(pre.ambigous.sym,Glut2L.ambigous.sym))
length(gene.input)#11281
seu.ref.L <- ScaleData(seu.ref.L,features = gene.input)

pdf("Glut2L/LMgene/pca.scaledata.gene.sn.supp.cc.pdf")
PCA.supp <- FactoMineR::PCA(t(as.matrix(seu.ref.L@assays$RNA@data[gene.input,])),
                            graph = T,
                            quanti.sup = which(!gene.input %in% c(Glut2L.cc.pca.input.sym))
)

dev.off()

plot(PCA.supp$ind$coord)
library(FactoMineR)
beta.PC12.dim.res <- dimdesc(PCA.supp,
                             axes = c(1,2),
                             proba = 0.01
)

saveRDS(beta.PC12.dim.res,'Glut2L/LMgene/beta.PC12.dim.res.rds')

beta.pc12.dim1 <- as.data.frame(beta.PC12.dim.res$Dim.1$quanti) 
pos.pc12.dim1.ens <- na.omit(beta.pc12.dim1[-log10(beta.pc12.dim1$p.value) >= 4, ])
dim(pos.pc12.dim1.ens)#164

Glut2L.keep.sym <- rownames(MyGeneExp(as.matrix(seu.ref.L@assays$RNA@data)[rownames(seu.ref.L) %in%  rownames(pos.pc12.dim1.ens),],
                                      log(2),
                                      1
))
length(Glut2L.keep.sym)#164

Glut2L.keep.sym2 <- rownames(MyGeneExp(as.matrix(seu.ref.L@assays$RNA@data)[Glut2L.keep.sym,],
                                       log(2),
                                       exp = 'no',
                                       10
))
length(Glut2L.keep.sym2)#142


png('Glut2L/LMgene/exp1.noexp10.pc1.cc.logp4.png',2000,3000)
cc.row.treeL <-
  MyHeatmap(as.matrix(exp(seu.ref.L@assays$RNA@data))[Glut2L.keep.sym2,],
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

cc.row.treeL.den <- as.dendrogram(cc.row.treeL)
length(labels(cc.row.treeL.den[[2]]))#

Glut2L.keep.cc.sym <- Glut2L.keep.sym2#
saveRDS(Glut2L.keep.cc.sym,'Glut2L/LMgene/Glut2L.keep.cc.sym.rds')

#########preg##########
cl<-makeCluster(10)
gene.input <- setdiff(gene.input,c(pre.cc.sym,Glut2L.cc.sym))
preg.gene <- MyLM.parallel(seu.ref.L@assays$RNA@data[gene.input,colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
                           variable = sort(seu.ref.L$pseudotime))

# preg.gene.1 <- MyLM.parallel(seu.ref.L@assays$RNA@data[gene.input,colnames(seu.ref.L)[order(seu.ref.L$ges.pregScore1)]],
#                              variable = sort(seu.ref.L$ges.pregScore1))
saveRDS(preg.gene,'Glut2L/LMgene/preg.gene.rds')

preg.gene['Acss2']
preg.gene['Stat3']

# pdf('LMgene/Acss2.Stat3.pseudotime.log.pdf',7,7)
# Mygene2pseudotime(genes.inf.input = genes.inf.input,gene.list = c('Acss2','Stat3'),Time = seu.ref.L$Type,cols = ref.time.colors,tpm.data = as.matrix(seu.ref.L@assays$RNA@data),Pseudotime = seu.ref.L$pseudotime,log = T)
# dev.off()

#######ges post####
preg.gene.L <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/cluster/Glut2L/LMgene/gene_inf/preg.gene.si.1e5.addcc.tab')
preg.gene.H <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/LMgene/gene_inf/cg.preg.gene.si.1e28.addcc.tab')
library(Vennerable)
Glut2HL.gene.venn <- Venn(list('Glut2H' = sub('_','-',preg.gene.H$SymbolDedu[preg.gene.H$cluster == 'preg.up']),
                               'Glut2L' = sub('_','-',preg.gene.L$SymbolDedu[preg.gene.L$cluster=='preg.up'])))


Glut2HL.low.gene.venn <- Venn(list('Glut2H' = sub('_','-',preg.gene.H$SymbolDedu[preg.gene.H$cluster == 'preg.down']),
                                   'Glut2L' = sub('_','-',preg.gene.L$SymbolDedu[preg.gene.L$cluster=='preg.down'])))

Glut2HL.gene.venn.p <- Myvenn(Glut2HL.gene.venn,col.list = c('black','black'),lwd = 8)
Glut2HL.gene.venn.low.p <- Myvenn(Glut2HL.low.gene.venn,col.list = c('black','black'),lwd = 8)

pdf('Glut2L/LMgene/FigS3c.cg.all.Glut2LH.preg.down.gene.overlap.1e5.pdf'    )
plot(Glut2HL.low.gene.venn)
plot(Glut2HL.low.gene.venn,gp=Glut2HL.gene.venn.low.p)
dev.off()

pdf('Glut2L/LMgene/FigS3c.cg.all.Glut2LH.preg.up.gene.overlap.1e5.pdf')
plot(Glut2HL.gene.venn)
plot(Glut2HL.gene.venn,gp=Glut2HL.gene.venn.p)
dev.off()
#########
preg.gene.filter2 <- preg.gene[-log10(preg.gene)>5]
length(preg.gene.filter2)#316#471

group.col <- time.colors[c(8,5,7)]
names(group.col) <- c('beta1','beta2','beta3')
subgroup.col <- c("#E08698","#6a3d9a")
MyPlotColor(subgroup.col,3)

c1Color <- MyName2Col(seu.ref.L@meta.data[colnames(seu.ref.L)[order(seu.ref.L$pseudotime)],"Type"],
                      ref.time.colors)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.ref.L@meta.data[colnames(seu.ref.L)[order(seu.ref.L$pseudotime)],"subgroup"],
                      subgroup.col)
c2Color <- as.matrix(c2Color)
cColor <- cbind(c1Color,c2Color)
seu.ref.L <- readRDS('Glut2L/seu.ref.L.endo.celltype.rds')

# L.heatmap <- MyGeneExp(as.matrix(exp(seu.ref.L@assays$RNA@data))[names(preg.gene.filter2),colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
#                        2,exp = 'no',3
#                        )
# dim(L.heatmap)

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/cor2.1.8.pregTop1e4.heatmap.pdf',20,30)
preg.row.tree2 <-
  MyHeatmap(as.matrix(exp(seu.ref.L@assays$RNA@data))[names(preg.gene.filter2),colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
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
            Colv = 'none',
            dendrogram = 'row',
            graph = T)
dev.off()

preg.row.tree.den2 <- as.dendrogram(preg.row.tree2)

preg.2.sym <- labels(as.dendrogram(preg.row.tree2)[[1]])
preg.1.sym <- c(labels(as.dendrogram(preg.row.tree2)[[2]]))

length(preg.2.sym)#55
length(preg.1.sym)#416


preg.1.order <- MyordergenewithPseudotime(as.matrix(exp(seu.ref.L@assays$RNA@data))[preg.1.sym,colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
                                          graph = T,
                                          preg.1.sym)
preg.2.order <- MyordergenewithPseudotime(as.matrix(exp(seu.ref.L@assays$RNA@data))[preg.2.sym,colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
                                          graph = T,
                                          preg.2.sym)
beta.cc.order <- MyordergenewithPseudotime(as.matrix(exp(seu.ref.L@assays$RNA@data))[Glut2L.keep.cc.sym,colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]],
                                           graph = T,
                                           Glut2L.keep.cc.sym)
# length(beta.cc.sym2)#368
#########identify cc genes######
#beta.cc.sym <- pre.cc.sym#285
seu.ref.L <- ScaleData(seu.ref.L,features = c(beta.cc.order,preg.1.order,preg.2.order))

c1Color <- MyName2Col(seu.ref.L@meta.data[colnames(seu.ref.L)[order(seu.ref.L$pseudotime)],"Type"],
                      ref.time.colors2)
c1Color <- as.matrix(c1Color)
c2Color <- MyName2Col(seu.ref.L@meta.data[colnames(seu.ref.L)[order(seu.ref.L$pseudotime)],"subgroup"],
                      celltype.col2)
c2Color <- as.matrix(c2Color)
cColor <- cbind(c1Color,c2Color)

png("G:/lab/Article/heshuang/BYLW/sm3/ref/FigS3c.addcc.cor0.2.1.all.preg.reorder.10e4.scaledata.final.png",
    2000,2500)
# row.tree.count <- 2
r1Color <- c(#rep(gene.tree.colors[1],length(beta.cc.order)),
             rep(gene.tree.colors[2],length(preg.1.order)),
             rep(gene.tree.colors[3],length(preg.2.order))
)
r1Color <- as.matrix(t(r1Color))

MyHeatmap(as.matrix(as.matrix(exp(seu.ref.L@assays$RNA@scale.data))[c(#beta.cc.order,
                                                                     preg.1.order,
                                                                     preg.2.order),colnames(seu.ref.L)[order(seu.ref.L$pseudotime)]]),
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

#########
dir.create('Glut2L/LMgene/gene_inf')

preg.gene.si.tab <- genes.inf.input[c(#beta.cc.order,
                                      preg.1.sym,preg.2.sym),]

#preg.gene.si.tab$cluster <- 'cc'
preg.gene.si.tab[preg.1.sym,'cluster'] <- 'preg.up'
preg.gene.si.tab[preg.2.sym,'cluster'] <- 'preg.down'
preg.gene.si.tab$pval <- '/'
preg.gene.si.tab[c(preg.1.sym,preg.2.sym),]$pval <- preg.gene[c(preg.1.sym,preg.2.sym)]

preg.gene.cc.si.tab <- genes.inf.input[beta.cc.order,]
preg.gene.cc.si.tab$cluster <- 'cc'
preg.gene.cc.si.tab$pval <- beta.pc12.dim1[beta.cc.order,'p.value']


preg.gene.si.tab <- rbind(preg.gene.si.tab,preg.gene.cc.si.tab)
MyWriteTable(preg.gene.si.tab,'Glut2L/LMgene/gene_inf/preg.gene.si.1e4.addcc.tab')
###########
rm(seu.ref.L)
save.image('seu.ref.raw.beta.Glut2LH.RData')


source('G:/pcatest/MyFunction.R')
source('G:/r script/file_function.R')

##########cc gene######

celltype.col2 <- c("#E08698","#5F9EA0", "#A65628")
names(celltype.col2) <- c('Glut2L','Glut2H_1','Glut2H_2')

seu.beta$beta2 <- factor(seu.beta$beta2,levels = c('Glut2L','Glut2H_1','Glut2H_2'))


gene.input <- rownames(MyGeneExp(seu.beta@assays$RNA@data[setdiff(rownames(seu.beta)[rowSums(seu.beta@assays$RNA@data)>0],
                                                                  c(pre.ambigous.sym,pre.cc.sym)),],
                                 log(2),10))
length(gene.input)#13408

saveRDS(gene.input,'gene.input.exp2.10.rds')

pdf("LMgene/pca.gene.sn.supp.cc.pdf")
PCA.supp <- FactoMineR::PCA(t(as.matrix(seu.beta@assays$RNA@data[c(gene.input,pre.cc.sym),])),
                            graph = T,
                            quanti.sup = which(!c(gene.input,pre.cc.sym) %in% c(all.var.co2.cc))
)
dev.off()

plot(PCA.supp$ind$coord)

library(FactoMineR)
beta.PC12.dim.res <- dimdesc(PCA.supp,
                             axes = c(1,2),
                             proba = 0.01
)

saveRDS(beta.PC12.dim.res,'LMgene/beta.PC12.dim.res.rds')

beta.pc12.dim1 <- as.data.frame(beta.PC12.dim.res$Dim.1$quanti) 
pos.pc12.dim1.ens <- na.omit(beta.pc12.dim1[-log10(beta.pc12.dim1$p.value) >= 30, ])
dim(pos.pc12.dim1.ens)#421

sum(pre.cc.sym %in% rownames(pos.pc12.dim1.ens))


Glut2L.keep.sym <- rownames(MyGeneExp(as.matrix(seu.beta@assays$RNA@data)[rownames(seu.beta) %in%  rownames(pos.pc12.dim1.ens),],
                                      log(2),
                                      5
))
length(Glut2L.keep.sym)#419

Glut2L.keep.sym2 <- rownames(MyGeneExp(as.matrix(seu.beta@assays$RNA@data)[Glut2L.keep.sym,],
                                       log(2),
                                       exp = 'no',
                                       10
))
length(Glut2L.keep.sym2)#408

c1Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"Type"],
                      ref.time.colors)
c1Color <- as.matrix(c1Color)
# c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"heter.group.merge"],
#                       group.col)
# c2Color <- as.matrix(c2Color)
c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"beta2"],
                      celltype.col2)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)

png('LMgene/exp1.noexp10.pc1.cc.logp30.png',2000,3000)
cc.row.tree <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[Glut2L.keep.sym2,],
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

Glut2H.keep.cc.sym <- Glut2L.keep.sym2#411
saveRDS(Glut2L.keep.cc.sym,'LMgene/Glut2H.keep.cc.sym.rds')
###########
length(gene.input)
seu.beta <- SetIdent(seu.beta,value = seu.beta$Type)
DEG.G0_G14.5 <- Myseufindmarker(seu.beta,gene.include = gene.input,ident.1 = 'Virgin',ident.2 = 'G14.5',c1 = 'Virgin',c2 = 'G14.5')
MyWriteTable(DEG.G0_G14.5,'DEG.G0_G14.5.tab')


cl<-makeCluster(10)
preg.gene <- MyLM.parallel(seu.beta@assays$RNA@data[gene.input,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                           variable = sort(seu.beta$pseudotime))

preg.gene['Acss2']
preg.gene['Stat3']
preg.gene['Acly']

pdf('LMgene/Acss2.Stat3.pseudotime.log.pdf',7,7)
Mygene2pseudotime(genes.inf.input = genes.inf.input,gene.list = c('Acly','Acss2','Stat3'),Time = seu.beta$Type,cols = ref.time.colors,tpm.data = as.matrix(seu.beta@assays$RNA@data),Pseudotime = seu.beta$pseudotime,log = T)
dev.off()

###########preg##########
preg.gene.filter2 <- preg.gene[-log10(preg.gene)>=28]
length(preg.gene.filter2)#746#cg 769

group.col <- time.colors[c(8,5,7)]
names(group.col) <- c('beta1','beta2','beta3')


c1Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"Type"],
                      ref.time.colors2)
c1Color <- as.matrix(c1Color)
# c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"heter.group.merge"],
#                       group.col)
# c2Color <- as.matrix(c2Color)
c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"beta2"],
                      celltype.col2)
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)

png('cg.cor2.1.pregTop1e28.heatmap.png',2000,3000)
preg.row.tree27 <-
  MyHeatmap(as.matrix(exp(seu.beta@assays$RNA@data))[names(preg.gene.filter2),colnames(seu.beta)[order(seu.beta$pseudotime)]],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = cColor,
            ColSideColorsSize = 1.5,
            return.tree = "row",
            Colv = 'none',
            dendrogram = 'row',
            graph = T)
dev.off()

preg.1.sym <- labels(as.dendrogram(preg.row.tree27)[[2]])
preg.2.sym <- labels(as.dendrogram(preg.row.tree27)[[1]])
length(preg.2.sym)#180#188
length(preg.1.sym)#566#581

library(Vennerable)
cc.venn <- Venn(list('cc' = Glut2L.keep.sym2,
          'preg.1' = preg.1.sym,
          'preg.2' = preg.2.sym))
plot(cc.venn)
# 
# beta.cc.sym <- setdiff(Glut2L.keep.sym2,cc.venn@IntersectionSets$`110`)
# length(beta.cc.sym)#408




preg.1.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[preg.1.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                          graph = T,
                                          preg.1.sym)

preg.2.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[preg.2.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                          graph = T,
                                          preg.2.sym)
Glut2H.keep.cc.sym <- setdiff(Glut2H.keep.cc.sym,c(preg.1.sym,preg.2.sym))
beta.cc.order <- MyordergenewithPseudotime(as.matrix(exp(seu.beta@assays$RNA@data))[Glut2H.keep.cc.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
                                          graph = T,
                                          Glut2H.keep.cc.sym)
##########
seu.beta <- readRDS('cluster/seu.beta.ges.post.rds')
library(Seurat)
seu.beta <- ScaleData(seu.beta,features = c(#beta.cc.order,
                                            preg.1.order,preg.2.order))
########trend#########
MyPlotTrend(preg.2.sym,
            colnames(seu.beta)[order(seu.beta$pseudotime)],
            as.matrix(exp(seu.beta@assays$RNA@data)-1)[preg.2.sym,colnames(seu.beta)[order(seu.beta$Type)]],
            seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"pseudotime"])

MyPlotTrend(preg.1.sym,
            colnames(seu.beta)[order(seu.beta$pseudotime)],
            as.matrix(exp(seu.beta@assays$RNA@data)-1)[preg.1.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]],
            seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"pseudotime"])

MyGetTrend(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"pseudotime"],
           as.matrix(exp(seu.beta@assays$RNA@data)-1)[preg.1.sym,colnames(seu.beta)[order(seu.beta$pseudotime)]]
           )

debugonce(MyPlotTrend)

library(ggplot2)
ggplot(seu.beta@assays$RNA@data, aes(pseudotime, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # Line for each gene
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.3, color = "blue") + 
  # Trend line
  geom_smooth(size = 2, se = FALSE, color = "orange") +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  theme_classic()

#########heatmap######
load('../cg.LMgene.RData')
c1Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"Type"],
                      ref.time.colors2)
c1Color <- as.matrix(c1Color)
# c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"heter.group.merge"],
#                       group.col)
# c2Color <- as.matrix(c2Color)
c2Color <- MyName2Col(seu.beta@meta.data[colnames(seu.beta)[order(seu.beta$pseudotime)],"beta2"],
                      celltype.col2[c(3,4)])
c2Color <- as.matrix(c2Color)

cColor <- cbind(c1Color,c2Color)
png("G:/lab/Article/heshuang/BYLW/sm3/ref/Fig1d.cg.cor0.2.1.all.preg.reorder.10e28.scaledata.final.2.png",
    2000,2500)
# row.tree.count <- 2
r1Color <- c(#rep(gene.tree.colors[1],length(beta.cc.order)),
             rep(gene.tree.colors[2],length(preg.1.order)),
             rep(gene.tree.colors[3],length(preg.2.order))
)
r1Color <- as.matrix(t(r1Color))
MyText(paste('c1:',length(beta.cc.order),'\n','c2:',length(preg.1.order),'\n','c3:',length(preg.2.order),sep= ''))
MyHeatmap(as.matrix(as.matrix(exp(seu.beta@assays$RNA@scale.data))[c(#beta.cc.order,
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

#########
dir.create('LMgene/gene_inf')
dir.create('LMgene/go_kegg')
rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)

preg.gene.si.tab <- genes.inf.input[c(#beta.cc.order,
                                      preg.1.sym,preg.2.sym),]

preg.gene.si.tab$cluster <- '/'
preg.gene.si.tab[preg.1.sym,'cluster'] <- 'preg.up'
preg.gene.si.tab[preg.2.sym,'cluster'] <- 'preg.down'
preg.gene.si.tab$pval <- '/'

preg.gene.si.tab[c(preg.1.sym,preg.2.sym),]$pval <- preg.gene[c(preg.1.sym,preg.2.sym)]



preg.gene.cc.si.tab <- genes.inf.input[beta.cc.order,]
preg.gene.cc.si.tab$cluster <- 'cc'
preg.gene.cc.si.tab$pval <- beta.pc12.dim1[beta.cc.order,'p.value']

preg.gene.si.all.tab <- rbind(preg.gene.si.tab,preg.gene.cc.si.tab)
MyWriteTable(preg.gene.si.all.tab,'LMgene/gene_inf/cg.preg.gene.si.1e28.addcc.tab')

MyGOwritetable(preg.gene.si.all.tab$EnsemblGeneID[preg.gene.si.all.tab$cluster=='cc'],'LMgene/go_kegg/go.cc.tab')
#MyGOwritetable(preg.gene.si.tab$EnsemblGeneID[preg.gene.si.tab$cluster=='preg.down'],'LMgene/go.c3.tab')

####
rm(seu.beta)
load('../cg.LMgene.RData')
save.image('cg.LMgene.RData')

############GSEA##########
dir.create('LMgene/GSEA/')
setwd('LMgene/GSEA/')
seu.beta <- readRDS('../../cluster/seu.beta.ges.post.rds')
beta.sample.class <- sort(seu.beta$pseudotime,decreasing = F)
MyGSEA.prepare2(as.matrix(seu.beta@assays$RNA@data)[gene.input[!grepl('^-',gene.input)],colnames(seu.beta)[order(seu.beta$pseudotime,decreasing = F)]],
                sample.class = beta.sample.class,
                gct.name = "cg.beta.symbol.1.gct",
                cls.name = "cg.beta.cls",
                mode = "continuous",
                continuous_value_name = "preg"
)

#########
source('G:/pcatest/MyFunction.R')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/')
go_kegg.pos <- MyReadDelim('LMgene/GSEA/cor0.2.20221204.Gsea.1670128838895/gsea_report_for_preg_pos_1670128838895.tsv')
go_kegg.pos$GS.br..follow.link.to.MSigDB <- tolower(go_kegg.pos$GS.br..follow.link.to.MSigDB)
go_kegg.pos$NAME <- tolower(go_kegg.pos$NAME)
MyWriteTable(go_kegg.pos,'LMgene/GSEA/cor0.2.1.gokeggv7.5.1.preg.high.tab')


go_kegg.neg <- MyReadDelim('LMgene/GSEA/cor0.2.20221204.Gsea.1670128838895/gsea_report_for_preg_neg_1670128838895.tsv')
go_kegg.neg$GS.br..follow.link.to.MSigDB <- tolower(go_kegg.neg$GS.br..follow.link.to.MSigDB)
go_kegg.neg$NAME <- tolower(go_kegg.neg$NAME)
MyWriteTable(go_kegg.neg,'LMgene/GSEA/cor0.2.1.gokeggv7.5.1.preg.low.tab')



go_kegg.pos <- MyReadDelim('LMgene/GSEA/cg.cor0.2.20221214.Gsea.1671007908218/gsea_report_for_preg_pos_1671007908218.tsv')
go_kegg.pos$GS.br..follow.link.to.MSigDB <- tolower(go_kegg.pos$GS.br..follow.link.to.MSigDB)
go_kegg.pos$NAME <- tolower(go_kegg.pos$NAME)
MyWriteTable(go_kegg.pos,'LMgene/GSEA/cg.cor0.2.1.gokeggv7.5.1.preg.high.tab')


go_kegg.neg <- MyReadDelim('LMgene/GSEA/cg.cor0.2.20221214.Gsea.1671007908218/gsea_report_for_preg_neg_1671007908218.tsv')
go_kegg.neg$GS.br..follow.link.to.MSigDB <- tolower(go_kegg.neg$GS.br..follow.link.to.MSigDB)
go_kegg.neg$NAME <- tolower(go_kegg.neg$NAME)
MyWriteTable(go_kegg.neg,'LMgene/GSEA/cg.cor0.2.1.gokeggv7.5.1.preg.low.tab')

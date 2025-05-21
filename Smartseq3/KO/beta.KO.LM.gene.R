######pc2#######
pc.gene.list <- list()
seu.ref.KO.list2[[seu]]@reductions$pca@feature.loadings[,'PC_2']
print(seu.ref.KO.list2[[seu]][["pca"]], dims = 2, nfeatures = 10)
seu.ref.KO.list2[[seu]] <- JackStraw(seu.ref.KO.list2[[seu]], num.replicate = 100)
seu.ref.KO.list2[[seu]] <- ScoreJackStraw(seu.ref.KO.list2[[seu]], dims = 1:20)
pc.gene.list[[seu]] <- PCASigGenes(seu.ref.KO.list2[[seu]],pcs.use = 2,pval.cut = 0.001,use.full = FALSE,max.per.pc = NULL)

PCA.supp <- list()
PCA.supp[[seu]] <- FactoMineR::PCA(t(seu.ref.KO.list2[[seu]]@assays$RNA@data[VariableFeatures(seu.ref.KO.list2[[seu]]),]),graph = F)
library(FactoMineR)
beta.atac.dim.res <- list()
beta.atac.dim.res[[seu]] <- dimdesc(PCA.supp[[seu]],
                             axes = c(1,2),
                             proba = 1
)

beta.atac.pc12.dim2 <- na.omit(as.data.frame(beta.atac.dim.res[[seu]]$Dim.2$quanti))
beta.PC2.gene <-  beta.atac.pc12.dim2[-log10(beta.atac.pc12.dim2$p.value) >= -log10(0.01) , ]
dim(beta.PC2.gene)
pc.gene.list[[seu]] <- list()
pc.gene.list[[seu]][['pos']] <- rownames(beta.atac.pc12.dim2[-log10(beta.atac.pc12.dim2$p.value) >= -log10(0.001)& beta.atac.pc12.dim2$correlation>0, ])
pc.gene.list[[seu]][['neg']] <- rownames(beta.atac.pc12.dim2[-log10(beta.atac.pc12.dim2$p.value) >= -log10(0.001)& beta.atac.pc12.dim2$correlation<0, ])
sapply(pc.gene.list[[seu]], length)

#for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  for(tmp.cluster in names(pc.gene.list[[seu]])){
    MyGOwritetable(genes.inf.input[pc.gene.list[[seu]][[tmp.cluster]],1],universe.gene = c(genes.inf.input[unlist(pc.gene.list[[seu]]),1],preg.gene.si.tab$EnsemblGeneID[preg.gene.si.tab$cluster!='cc']),pvalue = 1,paste('go_kegg/go.bg.pc20.01',seu,tmp.cluster,'tab',sep = '.'))
    #MyGOwritetable(genes.inf.input[preg.gene.list[[seu]][[tmp.cluster]],1],pvalue = 1,paste('go_kegg/go.pval11.pc2.LM.Koscore',seu,tmp.cluster,'tab',sep = '.'))
    
  }
#}


###########
seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')

gene.input3 <- setdiff(gene.input2,c(exclude.P300.gene,exclude.stat3HMKO.gene,Novaseq.batch.gene))
preg.gene <- list()
cl<-makeCluster(10)
preg.gene.cluster <- list()
preg.gene.cluster[[seu]] <- list()
preg.gene.cluster[[seu]][['preg']] <- preg.gene[[seu]]
# preg.gene[[seu]] <- MyLM.parallel(seu.ref.KO.list2[[seu]]@assays$RNA@data[gene.input3,colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$pregScore)]],
#                            variable = sort(seu.ref.KO.list2[[seu]]$pregScore))
seu <- 'G0G14.5P300WTHMKO'
seu.ref.KO.list2[[seu]]$KOScore <- Embeddings(seu.ref.KO.list2[[seu]],'pca')[,2]
seu.ref.KO.list2[[seu]]$KOScore <- seu.ref.KO.list2[[seu]]$cg.koscore 
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  preg.gene[[seu]] <- MyLM.parallel(seu.ref.KO.list2[[seu]]@assays$RNA@data[gene.input3,colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore)]],
                                    variable = sort(seu.ref.KO.list2[[seu]]$KOScore))
}

preg.gene[[seu]]['Acss2']
preg.gene[[seu]]['Ep300']
preg.gene[[seu]]['Stat3']
preg.gene[[seu]]['Acly']
saveRDS( preg.gene,'20221106/20221121/p300/ preg.gene.rds')
# pdf('LMgene/Acss2.Stat3.KOScore.log.pdf',7,7)
# Mygene2KOScore(genes.inf.input = genes.inf.input,gene.list = c('Acly','Acss2','Stat3'),Time = seu.ref.KO.list2[[seu]]$Type,cols = ref.time.colors,tpm.data = as.matrix(seu.ref.KO.list2[[seu]]@assays$RNA@data),KOScore = seu.ref.KO.list2[[seu]]$KOScore,log = T)
# dev.off()
# saveRDS(preg.gene,'preg.gene.rds')
###########preg##########
preg.gene.filter <- list()
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  preg.gene.filter[[seu]] <- preg.gene[[seu]][-log10(preg.gene[[seu]])>=11]
}
sapply(preg.gene.filter, length)
dir.create('20221106/20221121/p300')
preg.row.tree <- list()
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore)],"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
png(paste('20221106/20221121/p300/',seu,'cg.koscore.pregTop.heatmap.png',sep = ''),2000,3000)
preg.row.tree[[seu]] <-
  MyHeatmap(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[names(preg.gene.filter[[seu]]),colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore)]],
            type = "log.row.relat",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = c1Color,
            ColSideColorsSize = 1.5,
            return.tree = "row",
            Colv = 'none',
            dendrogram = 'row',
            graph = T)
dev.off()
}

preg.gene.list <- list()
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  preg.gene.list[[seu]] <- list()
  preg.gene.list[[seu]][['WT']] <- labels(as.dendrogram(preg.row.tree[[seu]])[[2]])
  preg.gene.list[[seu]][['HMKO']] <- labels(as.dendrogram(preg.row.tree[[seu]])[[1]])
}
preg.gene.si.tab
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  for(tmp.cluster in names(preg.gene.list[[seu]])){
    #MyGOwritetable(genes.inf.input[preg.gene.list[[seu]][[tmp.cluster]],1],universe.gene = c(genes.inf.input[unlist(preg.gene.list[[seu]]),1],preg.gene.si.tab$EnsemblGeneID[preg.gene.si.tab$cluster!='cc']),pvalue = 1,paste('go_kegg/go.bg.pval11.pc2.LM.Koscore',seu,tmp.cluster,'tab',sep = '.'))
    MyGOwritetable(genes.inf.input[preg.gene.list[[seu]][[tmp.cluster]],1],pvalue = 1,paste('20221106/20221121/p300/go.pval11.cg.LM.Koscore',seu,tmp.cluster,'tab',sep = '.'))
    
  }
}

dir.create('gene_inf')
LM.gene.inf <- list()
for(seu in names(seu.ref.KO.list2)[c(2,3)]){
  LM.gene.inf[[seu]] <- cbind(genes.inf.input[names(preg.gene.filter[[seu]]),],preg.gene.filter[[seu]])
  LM.gene.inf[[seu]]$cluster <- 'WT'
  LM.gene.inf[[seu]][preg.gene.list[[seu]][['HMKO']],'cluster'] <- 'HMKO'
  MyWriteTable(LM.gene.inf[[seu]],paste('20221106/20221121/p300/',seu,'.cg.LMpval11.tab'))
}


for(seu in names(seu.ref.KO.list2)[c(2,3)]){
preg.1.order <- MyordergenewithPseudotime(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[preg.gene.list[[seu]][['WT']],colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore,decreasing = T)]],
                                          graph = T,
                                          preg.gene.list[[seu]][['WT']])

preg.2.order <- MyordergenewithPseudotime(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[preg.gene.list[[seu]][['HMKO']],colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore,decreasing = T)]],
                                          graph = T,
                                          preg.gene.list[[seu]][['HMKO']])
c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore)],"Type"],
                      colors.list[[seu]])
c1Color <- as.matrix(c1Color)

pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/KO/',seu,".cg.LM.zscore.pval11.pdf",sep =  ''),
    10,13)
r1Color <- c(#rep(gene.tree.colors[1],length(beta.cc.order)),
  rep(gene.tree.colors[2],length(preg.1.order)),
  rep(gene.tree.colors[3],length(preg.2.order))
)
r1Color <- as.matrix(t(r1Color))

MyHeatmap(as.matrix(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[c(#beta.cc.order,
  preg.1.order,
  preg.2.order),colnames(seu.ref.KO.list2[[seu]])[order(seu.ref.KO.list2[[seu]]$KOScore)]]),
  type = "log.row.zscore",
  #hc.c.data.type = "log.row.relat",
  hc.r.data.type = "log.row.relat",
  #  color.palette = pc12.heatmap.col,
  #c.cov.method = "s",
  # r.cov.method = "s",
  Colv = "none",
  Rowv = "none",
  dendrogram = "none",
  #c.hc.method = "ward.D2",
  #r.hc.method = "ward.D2",
  RowSideColors = r1Color,
  RowSideColorsSize = 1.5,
  ColSideColors = c1Color,
  ColSideColorsSize = 1.5#,
  #return.tree = "none"
)
dev.off()
}


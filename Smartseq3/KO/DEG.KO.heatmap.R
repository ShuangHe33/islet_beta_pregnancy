setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/KOdirect/20220913/20221106/20221121/')
load('KOdirect.RData')

seu <- 'G0G14.5P300WTHMKO'
seu.ref.KO.list2[[seu]]$Type_rep <- paste(seu.ref.KO.list2[[seu]]$Type,seu.ref.KO.list2[[seu]]$Rep,sep = '_')
seu.ref.KO.list2[[seu]]$Type_rep <- factor(seu.ref.KO.list2[[seu]]$Type_rep,levels = c('Virgin_P300KOWT_rep1','Virgin_P300KOWT_rep2',
                                                                                       'Virgin_P300HMKO_rep1','Virgin_P300HMKO_rep2',
                                                                                       'G14_5_P300WT_rep1','G14_5_P300WT_rep3',
                                                                                       'G14_5_P300HMKO_rep1','G14_5_P300HMKO_rep2'
                                                                                       ))
seu <- 'G0G18.5Stat3WTHMKO'
seu.ref.KO.list2[[seu]]$Type_rep <- paste(seu.ref.KO.list2[[seu]]$Type,seu.ref.KO.list2[[seu]]$Rep,sep = '_')
seu.ref.KO.list2[[seu]]$Type_rep <- factor(seu.ref.KO.list2[[seu]]$Type_rep,levels = c('G0_Stat3KOWT_rep1','G0_Stat3KOWT_rep2',
                                                                                       'Virgin_Stat3HMKO_rep1','Virgin_Stat3HMKO_rep2',
                                                                                       'G18.5_Stat3WT_rep1','G18.5_Stat3WT_rep3',
                                                                                       'G18.5_Stat3HMKO_rep1','G18.5_Stat3HMKO_rep2'
))

seu <- 'G0G18.5Acss2WTHMKO'
seu.ref.KO.list2[[seu]]$Type_rep <- paste(seu.ref.KO.list2[[seu]]$Type,seu.ref.KO.list2[[seu]]$Rep,sep = '_')
seu.ref.KO.list2[[seu]]$Type_rep <- factor(seu.ref.KO.list2[[seu]]$Type_rep,levels = c('G0_Acss2KOWT_rep1','G0_Acss2KOWT_rep2',
                                                                                       'G0_Acss2KOHMKO_rep1','G0_Acss2KOHMKO_rep2',
                                                                                       'G18.5_Acss2WT_rep1','G18.5_Acss2WT_rep2',
                                                                                       'G18.5_Acss2HMKO_rep1','G18.5_Acss2HMKO_rep2'
))

gene.input2 <- setdiff(gene.input,Novaseq.batch.gene)
seu <- "G0G18.5Acss2WTHMKO"

seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
KO.DEG.list <- readRDS('KO.DEG.list.20230704.rds')

seu <- "G0G18.5Acss2WTHMKO";type1 <- names(table(seu.ref.KO.list2[[seu]]$Type))[3];type2 <- names(table(seu.ref.KO.list2[[seu]]$Type))[4]
parame.list[[seu]] <- c(0.05,1.5,1.5)
parame.list[[seu]] <- c(0.05,1.5,1.5)
seu <- 'G0G18.5Stat3WTHMKO';type1 <- 'G0_Stat3KOWT';type2 <- 'Virgin_Stat3HMKO'
saveRDS(gene.input2,'gene.input2.rds')
KO.DEG.list <- KO.DEG.list[c(28:32,22,33)]
seu.ref.KO.list2 <- seu.ref.KO.list2[1:3]
KO.DEG.list.filter <- list()
fc.tmp <- 1.3
for(seu in names(seu.ref.KO.list2)){
  seu.ref.KO.list2[[seu]] <- SetIdent(seu.ref.KO.list2[[seu]],value = seu.ref.KO.list2[[seu]]$Type)
  parame.list[[seu]] <- c(0.05,fc.tmp,fc.tmp)
  type.list <- names(table(seu.ref.KO.list2[[seu]]$Type))
  for(type1 in type.list){
    print(type1)
    type.list <- type.list[-1]
    for(type2 in type.list){
      print(type2)
      # KO.DEG.list[[paste(type1,type2,sep = '_')]] <- Myseufindmarker(seu.ref.KO.list2[[seu]],gene.include = gene.input2,ident.1 = type1,
      #                                                                 ident.2 = type2,c1 = type1,c2 = type2)
       #KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster <- factor(KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster,levels = c(type1,type2))
      
       KO.DEG.list.filter[[paste(type1,type2,sep = '_')]] <-  KO.DEG.list[[paste(type1,type2,sep = '_')]][KO.DEG.list[[paste(type1,type2,sep = '_')]]$p_val_adj<=parame.list[[seu]][1]&
                                                                                                           abs( KO.DEG.list[[paste(type1,type2,sep = '_')]]$avg_logFC) >= log(parame.list[[seu]][3]) &
                                                                                                           ( KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.1>=0.3 |  KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.2>=0.3),]
       MyDEGfilterplot(seu.ref.KO.list2[[paste(type1,type2,sep = '_')]],KO.DEG.list.filter[[paste(type1,type2,sep = '_')]],KO.DEG.list[[paste(type1,type2,sep = '_')]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/KO/',type1,type2,'pval',parame.list[[seu]][1],'fc',parame.list[[seu]][2],sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = time.colors[7:8],text.cex=1,pval.gap = 5,x.gap = 0.5,x.ext = -0.5,ext =20,point.size = 2)
       MyDEGfilterplot(seu.ref.KO.list2[[paste(type1,type2,sep = '_')]],KO.DEG.list.filter[[paste(type1,type2,sep = '_')]],KO.DEG.list[[paste(type1,type2,sep = '_')]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/KO/text.',type1,type2,'pval',parame.list[[seu]][1],'fc',parame.list[[seu]][2],sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = colors.list[[seu]][1:2],text.cex=1,pval.gap = 5,y.cex = 0.01,x.gap = 0.5,x.ext = -0.5,ext = 20,point.size = 2)
      # DEG <- paste(type1,type2,sep = '_')
      # DEG.row.tree.list[[DEG]] <- MyDEGfilterplot(seu.ref.KO.list2[[seu]],KO.DEG.list.filter[[DEG]],KO.DEG.list[[DEG]],prefix = paste('DEG/padj0.05.fc1.3.tmp',DEG,sep = ''),plotheatmap = T,
      #                                             display.cells = colnames(seu.ref.KO.list2[[seu]])[seu.ref.KO.list2[[seu]]$Type %in% c(type1,type2)],plotfireplot = F,plotRep = T,heatmap.col1 = ko.time.col)
      # 
       
    }
  }
}

labels(as.dendrogram(DEG.row.tree.list[[DEG]])[[2]][[2]][[2]][[2]][[1]][[2]][[1]])
labels(as.dendrogram(DEG.row.tree.list[[DEG]])[[1]][[2]][[1]][[1]])
grep('Hist1h2al',labels(as.dendrogram(DEG.row.tree.list[[DEG]])[[1]][[2]][[1]][[2]][[1]]))
G18.5.WT.exclude.gene <- labels(as.dendrogram(DEG.row.tree.list[[DEG]])[[1]][[2]][[1]][[2]][[1]])
KO.DEG.list$G18.5_Acss2WT_G18.5_Acss2HMKO <- KO.DEG.list$G18.5_Acss2WT_G18.5_Acss2HMKO[!KO.DEG.list$G18.5_Acss2WT_G18.5_Acss2HMKO$gene %in% G18.5.WT.exclude.gene,]

rm(seu.ref.KO.list2)
save.image('KOdirect.RData')
saveRDS(KO.DEG.list,'KO.DEG.list.rds')
saveRDS(KO.DEG.list,'KO.DEG.list.20230704.rds')
KO.DEG.list <- readRDS('KO.DEG.list.rds')
MyWriteTable(cbind(genes.inf.input[KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO$gene,],KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO),
             'DEG/gene_inf/G18.5Acss2WTHMKO.padj0.05.fc1.3.si.tab')

MyWriteTable(cbind(genes.inf.input[KO.DEG.list.filter$G18.5_Stat3WT_G18.5_Stat3HMKO$gene,],KO.DEG.list.filter$G18.5_Stat3WT_G18.5_Stat3HMKO),
             'DEG/gene_inf/G18.5Stat3WTHMKO.padj0.05.fc1.3.si.tab')
MyWriteTable(cbind(genes.inf.input[KO.DEG.list.filter$G0_Stat3KOWT_Virgin_Stat3HMKO$gene,],KO.DEG.list.filter$G0_Stat3KOWT_Virgin_Stat3HMKO),
             'DEG/gene_inf/G0Stat3WTHMKO.padj0.05.fc1.3.si.tab')

MyWriteTable(cbind(genes.inf.input[KO.DEG.list.filter$G0_Stat3KOWT_Virgin_Stat3HMKO$gene,],KO.DEG.list.filter$G0_Stat3KOWT_Virgin_Stat3HMKO),
             'DEG/gene_inf/G0Stat3WTHMKO.padj0.05.fc1.5.si.tab')
for(tmp.cluster in names(table(KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO$cluster))){
  MyGOwritetable(genes.inf.input[KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO$gene[KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO$cluster==tmp.cluster],1],
                 paste('DEG/go_kegg/go.',tmp.cluster,'.tab',sep = ''))
}
############
genesKO.list <- list()
genesKO.list[[names(seu.ref.KO.list2)[1]]] <- 'Acss2'
genesKO.list[[names(seu.ref.KO.list2)[2]]] <- 'Stat3'
genesKO.list[[names(seu.ref.KO.list2)[3]]] <- 'P300'
DEG.mrge.list <- list()

for(seu in names(seu.ref.KO.list2)){
  DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- unique(do.call(rbind,KO.DEG.list.filter[grep(genesKO.list[[seu]],names(KO.DEG.list.filter))])$gene)
}
sapply(DEG.mrge.list, length)
rep.log.tp0.1m.list <- list()
merge.log.tp0.1m.list <- list()

fc.tmp <- 1.3;DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- KO.DEG.list.filter$G18.5_Acss2WT_G18.5_Acss2HMKO$gene

for(seu in names(seu.ref.KO.list2)){
  #rep.log.tp0.1m.list[[seu]] <- do.call(cbind,by(t(seu.ref.KO.list2[[seu]]@assays$RNA@data),seu.ref.KO.list2[[seu]]$Type_rep,colMeans))
  #merge.log.tp0.1m.list[[seu]] <- do.call(cbind,by(t(seu.ref.KO.list2[[seu]]@assays$RNA@data),seu.ref.KO.list2[[seu]]$Type,colMeans))
  
  png(paste('DEG/fc',fc.tmp,'.',seu,'.rep.G18.5WTKO.DEG.heatmap.raw.png',sep = ''),2000,3000)
  par(oma = c(8,0,0,0))
  DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]] <-
    MyHeatmap(as.matrix(rep.log.tp0.1m.list[[seu]])[DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]],],
              type = "row.zscore",
              hc.c.data.type = "row.relat",
              hc.r.data.type = "row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              #  ColSideColors = cColor,
              ColSideColorsSize = 2,
              Colv = 'none',
              dendrogram = 'row',
              return.tree = "row",
              graph = T)
  dev.off()
}

genes.list <- list()
#for(seu in names(seu.ref.KO.list2)){
seu <- 'G0G18.5Acss2WTHMKO';fc.tmp <- 1.5
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- list()
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster2']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[1]]),
                                                                       labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[2]][[1]])
                                                                       )
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster1']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[2]][[2]]))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]]))
  
  seu <- 'G0G18.5Stat3WTHMKO';fc.tmp <- 1.5
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- list()
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster2']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[1]]),
                                                                       labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[2]][[2]])
  )
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster1']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[2]][[1]]))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[2]]))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster4']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[1]]))
  
  
  seu <- 'G0G14.5P300WTHMKO';fc.tmp <- 1.5
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- list()
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster2']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[2]][[2]][[2]]))
  
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster1']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[2]][[1]]),
                                                                       labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[2]][[2]][[1]]
                                                                       ))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[1]][[2]]))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster4']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[1]][[1]]),labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]][[2]]))
  genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster5']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[1]]))
  #########
  
#}

sapply(genes.list[[seu]], length)
##########reorder######
select.merge.sample.list[[seu]] <- c('G0_Acss2KOWT','G18.5_Acss2WT','G18.5_Acss2HMKO')

colors.exp <- colorRampPalette(c("white","#FFFDFA","#FFF7E6","#FFE3A9","goldenrod1","#F08104","#EF7A01","darkorange2"),
                               space="Lab")(100)

for(seu in names(seu.ref.KO.list2)){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/KO/Fig5.DEGWTKO.color.merge.rowrelat.fc',fc.tmp,'.',seu,'.DEG.pdf',sep = ''),10,16)
  par(oma = c(8,0,0,0))
  ordr.tmp <- c()
  ordr.name <- c()
  for(tmp.cluster in sort(names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]]))){
    ordr.tmp <- c(ordr.tmp,genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]][sample(1:length(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]]))])
    ordr.name <- c(ordr.name,rep(tmp.cluster,length(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]])))
  }
  ordr.name <- factor(ordr.name,levels = sort(names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]])))
  
  par(oma  = c(5,0,0,0))
  MyHeatmap(as.matrix(merge.log.tp0.1m.list[[seu]][ordr.tmp,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
 # color.palette = colors.exp,
  c.cov.method = "s",
  r.cov.method = "s",
  c.hc.method = "ward.D2",
  r.hc.method = "ward.D2",
  ColSideColors = MyName2Col(1:4,add.alpha(c("#F781BF","#452cff","#A65628", "#5F9EA0"),0.9)),
  # ColSideColorsSize = 1.5,
  RowSideColors =  MyName2Col(sort(ordr.name),gene.tree.colors[c(1,2,2)],is.row = T),
  RowSideColorsSize = 1.5,
  Colv = 'none',
  Rowv = 'none',
  dendrogram='none',
  #  return.tree = "row",
  graph = T)
  dev.off()
}
######
dir.create('DEG/go_kegg')
dir.create('DEG/gene_inf')
for(seu in c("G0G18.5Acss2WTHMKO"
)){
  for(tmp.cluster in names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]])){
    MyGOwritetable(genes.inf.input[genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]],1],pvalue = 1,paste('DEG/go_kegg/DEG.go',seu,tmp.cluster,'tab',sep = '.'))
  }
}

MyGOwritetable(genes.inf.input[c(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster2']],
                                 genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']]
                                 ),1],pvalue = 1,paste('DEG/go_kegg/DEG.go',seu,'cluster2.3','tab',sep = '.'))


rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
WTHMKO.list <- list()
for(seu in c("G0G18.5Acss2WTHMKO"
)){
  
  WTHMKO.list[[seu]] <- genes.inf.input[DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]],]
  WTHMKO.list[[seu]]$hc.cluster <- '/'
  for(tmp.cluster in names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]])){
    WTHMKO.list[[seu]][genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]],'hc.cluster'] <- tmp.cluster
  }
  MyWriteTable(WTHMKO.list[[seu]],paste('DEG/gene_inf/G18.5WTKO.fc',fc.tmp,seu,'hc.cluster.tab',sep = '.'))
}
#######relative#####
merge.relative.tp0.1m.list <- list()
merge.relative.tp0.1m.list[[seu]] <- (exp(merge.log.tp0.1m.list[[seu]][,c('G0_Acss2KOWT','G18.5_Acss2WT','G18.5_Acss2HMKO')])-1)/rowMax(exp(merge.log.tp0.1m.list[[seu]][,c('G0_Acss2KOWT','G18.5_Acss2WT','G18.5_Acss2HMKO')])-1)
# pdf(paste('DEG/',seu,'relative.0.4.1.pdf'),5,5)
# count <- 0
# for(tmp.cluster in names(genes.list[[seu]])[5:9]){
#   count <- count+1
#   plot(colMeans(merge.relative.tp0.1m.list[[seu]][genes.list[[seu]][[tmp.cluster]],]),col = gene.tree.colors[count],main = tmp.cluster,
#        type='o',pch=20,cex=3,lwd=5,ylab = 'relative expression',ylim = c(0.4,1))
#   box(lwd=5)
# }
# dev.off()
# 
# rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
# 
# selected.gene.inf <- list()
# selected.gene.inf[[seu]] <- genes.inf.input[c(genes.list[[seu]][['preg.up']],genes.list[[seu]][['preg.down']],
#                                               genes.list[[seu]][['KO.up']]
#                                               
# ),]
# selected.gene.inf[[seu]]$cluster <- '/'
# count <- 0
# 
# for(tmp.cluster in c('preg.up','preg.down','KO.up')){
#   count <- count+1
#   selected.gene.inf[[seu]][genes.list[[seu]][[tmp.cluster]],'cluster'] <- paste('cluster',count,sep = '')
#   MyGOwritetable(genes.inf.input[genes.list[[seu]][[tmp.cluster]],1],pvalue = 1,paste('go_kegg/go.final.',seu,'.',paste('cluster',count,sep = ''),'.tab',sep = ''))
# }
# 
# MyWriteTable(selected.gene.inf[[seu]],'gene_inf/G18.5Acss2WTHMKO.DEG.heatmap.si.tab')
# 
# dir.create('gene_inf')
######selected genes######
go.cluster3.down <- MyReadDelim('../../overlap/final/go_kegg/go.Acss2HMKO.G18.5_Acss2HMKO.high.tab')
insulin.tem <- unlist(strsplit(go.cluster3.down[go.cluster3.down$ID=='GO:0030073','geneID'],'/'))
insulin.tem <- c(insulin.tem,'Mlxipl')
insulin.tem <- setdiff(insulin.tem,c('Gcg','Gipr','Rbp4','Cpt1a','Myt1','Glul','Baiap3','Ucp2'))
insulin.term.relative <- merge.relative.tp0.1m.list[[seu]][insulin.tem[insulin.tem %in% rownames(merge.relative.tp0.1m.list[[seu]])],]
insulin.term.relative[insulin.term.relative<0.5] <- 0.5
par(oma = c(5,0,0,5))
#MyHeatmap(insulin.term.relative,type = 'raw',Colv = 'none',Rowv = 'none',dendrogram='none',labRow=rownames(insulin.term.relative))
insulin.term.relative <- insulin.term.relative[order(insulin.term.relative[,1],decreasing = T),]
insulin.term.relative <- rbind(insulin.term.relative[1:7,][order(insulin.term.relative[1:7,3],decreasing = T),],insulin.term.relative[8:13,])
pdf('DEG/Figs8j.gene.heatmap.pdf',10,15)
par(oma = c(5,0,0,5))
pheatmap::pheatmap(insulin.term.relative,cluster_cols = F,cellwidth = 50,cellheight = 20,cluster_rows = F)
dev.off()
########pc2 genes#######
library(FactoMineR)
pca.res.list <- list()
dim.res.list <- list()
pc12.dim12.list <- list()
seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')

for(seu in names(seu.ref.KO.list2)){
  pca.res.list[[seu]] <- FactoMineR::PCA(t(seu.ref.KO.list2[[seu]]@assays$RNA@data[VariableFeatures(seu.ref.KO.list2[[seu]]),colnames(seu.ref.KO.list2[[seu]])]),graph = F)
  dim.res.list[[seu]] <- dimdesc(pca.res.list[[seu]],
                                 axes = c(1,2),
                                 proba = 0.05
  )
  pc12.dim12.list[[seu]][['pc1']] <- na.omit(as.data.frame(dim.res.list[[seu]]$Dim.1$quanti))
  pc12.dim12.list[[seu]][['pc2']] <- na.omit(as.data.frame(dim.res.list[[seu]]$Dim.2$quanti))
}

row.tree.list <- list()
pc12gene.list <- list()
pval=8
for (seu in names(seu.ref.KO.list2)) {
  print(seu)
  pc12gene.list[[seu]] <- list()
  for(pc.n in names(pc12.dim12.list[[seu]])){
    print(pc.n)
    order.sn <- rownames(pca.res.list[[seu]]$ind$coord)[order(pca.res.list[[seu]]$ind$coord[,grep(pc.n,names(pc12.dim12.list[[seu]]))],decreasing = F)]
    pc12gene.list[[seu]][[pc.n]] <- rownames(pc12.dim12.list[[seu]][[pc.n]])[-log10(pc12.dim12.list[[seu]][[pc.n]]$p.value)>=pval]
    
    c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[order.sn,"Type"],
                          ko.time.col[names(table(seu.ref.KO.list2[[seu]]@meta.data[order.sn,"Type"]))])
    c1Color <- as.matrix(c1Color)
    c2Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[order.sn,"Rep"],
                          rep.colors[-1])
    c2Color <- as.matrix(c2Color)
    
    cColor <- cbind(c1Color,c2Color)
    
    pdf(paste(seu,'.',pc.n,'.pval',pval,"order.col.heatmap.pdf",sep = ""),15,20)
    row.tree.list[[paste(seu,pc.n,sep = '')]] <-
      MyHeatmap(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[pc12gene.list[[seu]][[pc.n]],order.sn],
                type = "log.row.zscore",
                hc.c.data.type = "log.row.relat",
                hc.r.data.type = "log.row.relat",
                c.cov.method = "s",
                r.cov.method = "s",
                c.hc.method = "ward.D2",
                r.hc.method = "ward.D2",
                ColSideColors = cColor,
                ColSideColorsSize = 2,
                Colv = 'none',
                dendrogram='row',
                return.tree = "row",
                graph = T)
    dev.off()
  }
}

########
dir.create('go_kegg')
for (seu in names(seu.ref.KO.list2)) {
  print(seu)
  for(pc.n in names(pc12.dim12.list[[seu]])){
   MyGOwritetable(genes.inf.input[labels(as.dendrogram(row.tree.list[[paste(seu,pc.n,sep = '')]])[[1]]),1],pvalue = 1,paste('go_kegg/go.1.',seu,'.',pc.n,'.tab',sep = '')) 
   MyGOwritetable(genes.inf.input[labels(as.dendrogram(row.tree.list[[paste(seu,pc.n,sep = '')]])[[2]]),1],pvalue = 1,paste('go_kegg/go.2.',seu,'.',pc.n,'.tab',sep = '')) 
}
}

#############
rm(seu.ref.KO.list2)
save.image('KOdirect.RData')

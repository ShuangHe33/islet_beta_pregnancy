dir.create('DEG')
sn.DEG.list <- list()
sn.DEG.list[['Type']] <- c('Virgin','G14.5','G14.5_Acss2HMKO_HFD','G14.5_HFD',
                           'Virgin_HFD','G14.5HFDWT','G0Acss2HFDWT','G0Acss2HFDHMKO','G14.5CDWT',
                           'G14.5CDB6WT','G14.5HFDB6WT'
)

sn.DEG.list[['Type_rep']] <- list(c('Virgin_rep3','Virgin_rep4'),
                                  c('G14.5_rep1','G14.5_rep2','G14.5_rep3'),
                                  c('G14.5_Acss2HMKO_HFD_rep1','G14.5_Acss2HMKO_HFD_rep2'),
                                  c('G14.5_HFD_rep2','G14.5_HFD_rep3','G14.5_HFD_rep4'),
                                  c('Virgin_14.5HFD_rep2','Virgin_18.5HFD_rep1','Virgin_18.5HFD_rep2'),
                                  c('G14.5_Acss2OEWT_HFD_rep2','G14.5_Acss2KOWT_HFD_rep1'),
                                  c('G0_14d_Acss2WT_HFD_rep1','G0_Acss2KOWT_14d_HFD_rep2'),
                                  c('G0_14d_Acss2HMKO_HFD_rep1','G0_Acss2HMKO_14d_HFD_rep2','G0_Acss2HMKO_14d_HFD_rep3'),
                                  c('G14.5_Stat3WT_rep1','G14.5_Acss2OEWT_rep1'),
                                  c('G14.5_rep1','G14.5_rep2','G14.5_rep3','G14.5_Stat3WT_rep1','G14.5_Acss2OEWT_rep1'),
                                  c('G14.5_HFD_rep2','G14.5_HFD_rep3','G14.5_HFD_rep4','G14.5_Acss2OEWT_HFD_rep2','G14.5_Acss2KOWT_HFD_rep1')
)
names(sn.DEG.list[['Type_rep']]) <- sn.DEG.list[['Type']]




KO.DEG.list <- list()
KO.DEG.list.filter <- list()

seu <- 'ref_HMKOHFD'
seu.ko.list[[seu]] <- SetIdent(seu.ko.list[[seu]],value = seu.ko.list[[seu]]$Type_rep)
type.list <- c('G14.5CDB6WT','G14.5HFDB6WT','G14.5_Acss2HMKO_HFD')
type.list <- c('Virgin','Virgin_HFD')
type.list <- c('G14.5','G14.5_HFD')
type.list <- c('G0Acss2HFDWT','G0Acss2HFDHMKO')
type.list <- c('G14.5HFDWT','G14.5_Acss2HMKO_HFD')
type.list <- c('G14.5CDWT','G14.5HFDWT','G14.5_Acss2HMKO_HFD')

for(type1 in type.list){
  type.list <- type.list[-1]
  print(type1)
  for(type2 in type.list){
    print(type2)
    KO.DEG.list[[paste(type1,type2,sep = '_')]] <- Myseufindmarker(seu.ko.list[[seu]],gene.include = gene.input,ident.1 = sn.DEG.list[['Type_rep']][[type1]],ident.2 = sn.DEG.list[['Type_rep']][[type2]],c1 = type1,c2 = type2)
    KO.DEG.list.filter[[paste(type1,type2,sep = '_')]] <-  KO.DEG.list[[paste(type1,type2,sep = '_')]][KO.DEG.list[[paste(type1,type2,sep = '_')]]$p_val_adj<=0.01 &
                                                                                                         abs( KO.DEG.list[[paste(type1,type2,sep = '_')]]$avg_logFC) >= log(1.5) &
                                                                                                         ( KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.1>=0.3 |  KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.2>=0.3),]
    # row.tree.list[[paste(type1,type2,sep = '_')]] <- MyDEGfilterplot(seu.ko.list[[seu]],DEG.raw.filter = KO.DEG.list.filter[[paste(type1,type2,sep = '_')]],DEG.raw = KO.DEG.list[[paste(type1,type2,sep = '_')]],prefix = paste('',paste(type1,type2,sep = '_'),sep = ''),plotheatmap = T,display.cells = colnames(seu.ko.list[[seu]])[seu.ko.list[[seu]]$Type_rep %in% c(sn.DEG.list[['Type_rep']][[type2]],sn.DEG.list[['Type_rep']][[type1]])],plotfireplot = F,plotRep = T)
  }
}

type.list <- list()
type.list[['Virgin_Virgin_HFD']] <- list('Virgin','Virgin_HFD',0.01,1.5,10,40,2)
type.list[['Virgin_Virgin_HFD']] <- list('Virgin','Virgin_HFD',0.01,1.5,10,40,1)

type.list[['G14.5_G14.5_HFD']] <- list('G14.5','G14.5_HFD',0.01,1.5,10,10,1)

type.list[['G14.5HFDWT_G14.5_Acss2HMKO_HFD']] <- list('G14.5HFDWT','G14.5_Acss2HMKO_HFD',0.05,1.3,10,10,1)
type.list[['G0Acss2HFDWT_G0Acss2HFDHMKO']] <- list('G0Acss2HFDWT','G0Acss2HFDHMKO',0.05,1.3,10,40,2)




for(j in c(
  "Virgin_Virgin_HFD",
  "G14.5_G14.5_HFD",
  #'G14.5HFDWT_G14.5_Acss2HMKO_HFD',
  'G0Acss2HFDWT_G0Acss2HFDHMKO'
  )){
 # KO.DEG.list[[j]]$cluster <- factor(KO.DEG.list[[j]]$cluster,levels = unlist(type.list[[j]][c(1,2)]))
  # 
  # KO.DEG.list.filter[[j]] <-  KO.DEG.list[[j]][KO.DEG.list[[j]]$p_val_adj<=type.list[[j]][[3]]&
  #                                                abs( KO.DEG.list[[j]]$avg_logFC) >= log(type.list[[j]][[4]]) &
  #                                                ( KO.DEG.list[[j]]$pct.1>=0.3 |  KO.DEG.list[[j]]$pct.2>=0.3),]
  MyDEGfilterplot(seu.ko.list[[seu]],KO.DEG.list.filter[[j]],KO.DEG.list[[j]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/HFD/re.fc',type.list[[j]][[4]],'.padj',type.list[[j]][[3]],'.',j,sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = ko.time.col[unlist(type.list[[j]][c(1,2)])],text.cex=1,pval.gap = type.list[[j]][[5]],ext = type.list[[j]][[6]],point.size = type.list[[j]][[7]],x.gap = 0.5,x.ext = 0)
  MyDEGfilterplot(seu.ko.list[[seu]],KO.DEG.list.filter[[j]],KO.DEG.list[[j]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/HFD/re.text.fc',type.list[[j]][[4]],'.padj',type.list[[j]][[3]],'.',j,sep = ''),plotheatmap = F,plotfireplot = T,padj = F,fire.cols = ko.time.col[unlist(type.list[[j]][c(1,2)])],text.cex=1,pval.gap = type.list[[j]][[5]],y.cex = 0.01,ext = type.list[[j]][[6]],point.size = type.list[[j]][[7]],x.gap = 0.5,x.ext = 0)
}

saveRDS(KO.DEG.list,'DEG/KO.DEG.list.rds')
#######DEG heatmap######
dir.create('DEG/go_kegg')
dir.create('DEG/gene_inf')
rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
for(j in c(#"G14.5_B6G14.5_HFD",
 # "G14.5HFDWT_G14.5_Acss2HMKO_HFD"  ,
  "G14.5_G14.5_HFD"
)
){
  for(tmp.cluster in names(table(KO.DEG.list.filter[[j]]$cluster))){
    print(tmp.cluster)
    MyGOwritetable(genes.inf.input[KO.DEG.list.filter[[j]]$gene[KO.DEG.list.filter[[j]]$cluster==tmp.cluster],1],pvalue = 1,paste('DEG/go_kegg/go.padj0.01.fc1.5.',j,tmp.cluster,'tab',sep = '.'))
  }
}
for(j in c("Virgin_Virgin_HFD","G14.5_G14.5_HFD"
  #"G14.5HFDWT_G14.5_Acss2HMKO_HFD"  
)
){
  MyWriteTable(cbind(genes.inf.input[KO.DEG.list.filter[[j]]$gene,],KO.DEG.list.filter[[j]]),paste('DEG/gene_inf/padj0.01.fc1.5.',j,'.si.tab',sep = ''))
}
rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)

#######heatmap#######
DEG.mrge.list <- list()
fc.tmp <- 1.5
type.list <- c('G14.5CDB6WT','G14.5HFDB6WT','G14.5_Acss2HMKO_HFD')
for(type1 in type.list){
  type.list <- type.list[-1]
  print(type1)
  for(type2 in type.list){
    print(type2)
    KO.DEG.list.filter[[paste(type1,type2,sep = '_')]] <-  KO.DEG.list[[paste(type1,type2,sep = '_')]][KO.DEG.list[[paste(type1,type2,sep = '_')]]$p_val_adj<=0.01 &
                                                                                                         abs( KO.DEG.list[[paste(type1,type2,sep = '_')]]$avg_logFC) >= log(fc.tmp) &
                                                                                                         ( KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.1>=0.3 |  KO.DEG.list[[paste(type1,type2,sep = '_')]]$pct.2>=0.3),]
  }
}

for(seu in names(seu.ko.list)){
  DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- unique(do.call(rbind,KO.DEG.list.filter[grep('B6',names(KO.DEG.list.filter))])$gene)
}

sapply(DEG.mrge.list, length)#441
rep.log.tp0.1m.list <- list()
merge.log.tp0.1m.list <- list()
DEG.row.tree <- list()
for(seu in names(seu.ko.list)){
  #rep.log.tp0.1m.list[[seu]] <- do.call(cbind,by(t(seu.ko.list[[seu]]@assays$RNA@data),seu.ko.list[[seu]]$Type_rep,colMeans))
  #merge.log.tp0.1m.list[[seu]] <- do.call(cbind,by(t(seu.ko.list[[seu]]@assays$RNA@data),seu.ko.list[[seu]]$Type,colMeans))
  
  png(paste('DEG/fc',fc.tmp,'.',seu,'.rep.3conditions.DEG.heatmap.raw.png',sep = ''),2000,3000)
  par(oma = c(8,0,0,0))
  DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]] <-
    MyHeatmap(as.matrix(rep.log.tp0.1m.list[[seu]])[DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]],unlist(sn.DEG.list[['Type_rep']][c('G14.5CDB6WT','G14.5HFDB6WT','G14.5_Acss2HMKO_HFD')])],
              type = "row.relat",
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
#for(seu in names(seu.ko.list)){
fc.tmp <- 1.5
genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]] <- list()
genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster1']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[1]]))
genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster2']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[2]]))
genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']] <- c(labels(as.dendrogram(DEG.row.tree[[paste(seu,'_','fc',fc.tmp,sep = '')]])[[2]][[1]]))

sapply(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]], length)
##########reorder######
unlist(sn.DEG.list[['Type_rep']][c('G14.5CDB6WT','G14.5HFDB6WT','G14.5_Acss2HMKO_HFD')])
seu.ko.list[[seu]]$Type2 <- as.character(seu.ko.list[[seu]]$Type)
seu.ko.list[[seu]]@meta.data[seu.ko.list[[seu]]$Type %in% c('G14.5','G14.5CDWT'),'Type2'] <- 'G14.5CDB6WT'
seu.ko.list[[seu]]@meta.data[seu.ko.list[[seu]]$Type %in% c('G14.5_HFD','G14.5HFDWT'),'Type2'] <- 'G14.5HFDB6WT'

B6merge.log.tp0.1m.list <- list()
B6merge.log.tp0.1m.list[[seu]] <- do.call(cbind,by(t(seu.ko.list[[seu]]@assays$RNA@data),seu.ko.list[[seu]]$Type2,colMeans))


for(seu in names(seu.ko.list)){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/HFD/merge.color.merge.fc',fc.tmp,'.',seu,'.DEG.pdf',sep = ''),10,16)
  par(oma = c(8,0,0,0))
  # ordr.tmp <- c()
  # ordr.name <- c()
  # for(tmp.cluster in sort(names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]]))){
  #   ordr.tmp <- c(ordr.tmp,genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]][sample(1:length(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]]))])
  #   ordr.name <- c(ordr.name,rep(tmp.cluster,length(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]])))
  # }
  # ordr.name <- factor(ordr.name,levels = sort(names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]])))
  # 
  par(oma  = c(5,0,0,0))
  MyHeatmap(as.matrix(B6merge.log.tp0.1m.list[[seu]][ordr.tmp,
                                                   #c("G14.5","G14.5CDWT","G14.5_HFD","G14.5HFDWT","G14.5_Acss2HMKO_HFD")
                                                   c('G14.5CDB6WT','G14.5HFDB6WT','G14.5_Acss2HMKO_HFD')
                                                   ]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            #color.palette = colors.exp,
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = MyName2Col(1:3,add.alpha(ko.time.col[c("G14.5","G14.5_HFD","G14.5_Acss2HMKO_HFD")],0.9)),
            # ColSideColorsSize = 1.5,
            RowSideColors =  MyName2Col(sort(ordr.name),time.colors,is.row = T),
            RowSideColorsSize = 1.5,
            Colv = 'none',
            Rowv = 'none',
            dendrogram='none',
            #  return.tree = "row",
            graph = T)
  dev.off()
}
########
dir.create('DEG/go_kegg')
dir.create('DEG/gene_inf')

for(seu in c("G0G18.5Acss2WTHMKO"
)){
  for(tmp.cluster in names(genes.list[[seu]])[7:9]){
    MyGOwritetable(genes.inf.input[genes.list[[seu]][[tmp.cluster]],1],pvalue = 1,paste('DEG/go_kegg/go',seu,tmp.cluster,'tab',sep = '.'))
  }
}

WTHMKO.list <- list()
rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
for(seu in c("ref_HMKOHFD"
)){
  
  WTHMKO.list[[seu]] <- genes.inf.input[DEG.mrge.list[[paste(seu,'_','fc',fc.tmp,sep = '')]],]
  WTHMKO.list[[seu]]$hc.cluster <- '/'
  for(tmp.cluster in names(genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]])){
    WTHMKO.list[[seu]][genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][[tmp.cluster]],'hc.cluster'] <- tmp.cluster
  }
  MyWriteTable(WTHMKO.list[[seu]],paste('DEG/gene_inf/',seu,'hc.cluster.tab',sep = '.'))
}
#######
DEG.tab <- MyReadDelim('DEG/gene_inf/.ref_HMKOHFD.hc.cluster.tab')
MyGOwritetable(DEG.tab$EnsemblGeneID[DEG.tab$hc.cluster=='cluster1'],pvalue = 1,'DEG/go_kegg/go.cluster1.tab')
MyGOwritetable(DEG.tab$EnsemblGeneID[DEG.tab$hc.cluster!='cluster1'],pvalue = 1,'DEG/go_kegg/go.cluster2.3.tab')

######G18.5Acss2####
G18.5Acss2.DEG.tab <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/KOdirect/20220913/20221106/20221121/DEG/gene_inf/final.G0G18.5Acss2WTHMKO.hc.cluster.tab')
Acss2.up.ko <- Venn(list('G18.5Acss2KO' = sub('_','-',G18.5Acss2.DEG.tab$SymbolDedu[G18.5Acss2.DEG.tab$hc.cluster=='cluster2']),
          'G14.5Acss2HFDKO' = genes.list[[paste(seu,'_','fc',fc.tmp,sep = '')]][['cluster3']]))
plot(Acss2.up.ko)

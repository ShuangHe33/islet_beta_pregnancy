dir.create('G:/project/pregnant_mouse/chromatin/ATAC/20221106_addG0HFD_rep23')
setwd('G:/project/pregnant_mouse/chromatin/ATAC/20221106_addG0HFD_rep23/')
load('ATAC.RData')

#load('ATAC.G0G14.5G14.5HFD.RData')

source('G:/pcatest/MyFunction.R')
source('G:/r script/file_function.R')
library(ggplot2)

list.ATAC.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/ATAC/keep/list.ATAC.bdg.rds')
list.ATAC.bdg[["Ins1Ctrl_2_190829"]] <- NULL
list.ATAC.bdg[["Ins1G14.5_190710"]] <- NULL
list.ATAC.bdg[["Ins1G14.5_190807"]] <- NULL
list.ATAC.bdg$Ins1P7NB_190829 <- NULL
list.ATAC.bdg$Ins1P7NB_190808 <- NULL
list.ATAC.bdg$Ins1P7NB_190930 <- NULL


list.Acss2HFD.H3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/Acss2KO/HFD/new/list.bdg.Acss2HMKO.HFD.H3K27ac.rds')
list.Acss2HFD.H3K27ac.bdg$Ins1G14.5Acss2OEWT_HFD_N286_20220905 <- NULL

list.ATACHFD.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/ATAC/keep/HFD/list.ATAC.HFD.bdg.rds')
names(list.ATACHFD.bdg)
list.ATACHFD.bdg$Ins1G0HFD_220923_1 <- NULL
prefix.list <- list()
prefix.list[['ATAC']] <- c(names(list.ATAC.bdg)[1:2],names(list.ATACHFD.bdg)[c(1,4)],
                           names(list.ATAC.bdg)[3:4],names(list.ATACHFD.bdg)[2:3]
)

list.ATAC.bdg <- c(list.ATAC.bdg,list.ATACHFD.bdg)

list.ATAC.bed <-  readRDS('G:/project/pregnant_mouse/chromatin/ATAC/nbl.score.bed/list.consensus.peak.rds')
list.ATAC.bed$P7NB <- NULL
list.ATAC.bed$G0HFD <- NULL

sapply(list.ATAC.bed, length)

ATAC.nor.result <- MyCalUcscNormFactor2(list.ATAC.bed,list.ATAC.bdg,0.6)#370
saveRDS(ATAC.nor.result,'ATAC.nor.result.rds')
saveRDS(ATAC.nor.result$norm.factor,'ATAC.com.nor.factor.rds')

ATAC.prefix <-  prefix.list[['ATAC']]
rm(list.ATAC.bdg)
ATAC.merge.bed <- ATAC.nor.result$merge.bed

for(sn.tmp in ATAC.prefix){
  mcols(ATAC.merge.bed)[,sn.tmp] <- round(mcols(ATAC.merge.bed)[,sn.tmp],3)
  mcols(ATAC.merge.bed)[,sn.tmp][mcols(ATAC.merge.bed)[,sn.tmp]<0] <- 0
}
saveRDS(ATAC.merge.bed,'ATAC.merge.bed.rds')

list.G14.5HFD.H3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/HFD/list.H3K27ac.HFD.bdg.rds')
names(list.G14.5HFD.H3K27ac.bdg)[1] <- 'Ins1G14.5HFD_H3K27ac_220620'

list.Acss2HFD.H3K27ac.bdg <- c(list.G14.5HFD.H3K27ac.bdg[c(4:7)],list.Acss2HFD.H3K27ac.bdg[3:4])
prefix.list[['G14.5HFD.H3K27ac']] <- names(list.Acss2HFD.H3K27ac.bdg)

list.Acss2HFD.H3K27ac.bdg[['Ins1G0HFD_H3K27ac_2_220924']] <- list.G14.5HFD.H3K27ac.bdg$Ins1G0HFD_H3K27ac_2_220924
list.Acss2HFD.H3K27ac.bdg[['Ins1G0HFD_H3K27ac_3_220924']] <- list.G14.5HFD.H3K27ac.bdg$Ins1G0HFD_H3K27ac_3_220924

for(type in c(
  #'G14.5HFD.H3K27ac',
  'Acss2HFD.H3K27ac'
)){
  for(i in names(get(paste('list',type,'bdg',sep = '.')))[5:6]){
    print(i)
    mcols(ATAC.merge.bed)[,i] <- round(MyCalTd(ATAC.merge.bed,get(paste('list',type,'bdg',sep = '.'))[[i]],"td"),3)
    mcols(ATAC.merge.bed)[,i][mcols(ATAC.merge.bed)[,i]<0] <- 0
  }
}

boxplot(as.data.frame(ATAC.merge.bed)[,prefix.list$ATAC],outline=F)
#list.sf.factor <- readRDS('G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/G14.5Acss2HMKO/20221114/B6G0HFDG14.5CDHFD/list.sf.factor.rds')
list.sf.factor <- readRDS('G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/G14.5Acss2HMKO/20221114/WTG0HFDG14.5CDHFD/list.mean.sf.factor.20221203.rds')

###########
#list.sf.factor <- readRDS('G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/G14.5Acss2HMKO/20221114/CD_WT/list.sf.factor.rds')
prefix.list[['G14.5HFD.H3K27ac']] <- names(list.sf.factor$G14.5WTHFD.H3K27ac)
prefix.list[['G14.5HFD.H3K27ac']] <- names(list.Acss2HFD.H3K27ac.bdg)


list.nor.df <- list()
pdf('QC/H3K27ac.log.nor.pdf')
for(prefix in c(
  'G14.5HFD.H3K27ac'
  
)){
  print(prefix)
  list.nor.df[[prefix]] <- t(t(as.data.frame(ATAC.merge.bed)[,prefix.list[[prefix]]])/unlist(list.sf.factor$G14.5HFD.H3K27ac)[1:6])
  rownames(list.nor.df[[prefix]]) <- names(ATAC.merge.bed)
  for(sn.tmp in colnames(list.nor.df[[prefix]])){
    mcols(ATAC.merge.bed)[,paste(sn.tmp,'.nor',sep ='')] <- round(list.nor.df[[prefix]][,sn.tmp],3)
  }
  MyLogBoxPlot(list.nor.df[[prefix]],main=prefix,min.cutoff = 1,log = T)
}
dev.off()

list.sf.factor[['ATAC']] <- ATAC.nor.result$norm.factor

for(prefix in c(
 # 'G14.5HFD.H3K27ac',
  'ATAC'
  
)){
  print(prefix)
  list.nor.df[[prefix]] <- t(t(as.data.frame(ATAC.merge.bed)[,prefix.list[[prefix]]])/list.sf.factor$ATAC[prefix.list[[prefix]]])
  rownames(list.nor.df[[prefix]]) <- names(ATAC.merge.bed)
  for(sn.tmp in colnames(list.nor.df[[prefix]])){
    mcols(ATAC.merge.bed)[,paste(sn.tmp,'.nor',sep ='')] <- round(list.nor.df[[prefix]][,sn.tmp],3)
  }
  MyLogBoxPlot(list.nor.df[[prefix]],main=prefix,min.cutoff = 1,log = F)
}

MyLogBoxPlot(as.data.frame(ATAC.merge.bed)[,prefix.list$ATAC],main=prefix,min.cutoff = 1,log = F)

hist(log2(list.nor.df[[prefix]][,1]))

MyWriteBed(ATAC.merge.bed,output.col = paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'mean',sep = '.'),col.names = T,'merge.ATAC.bed')
###########
ATAC.si.tab <- data.frame('SampleName' = ATAC.prefix,
                          'Time' = rep(c('Ctrl','G0HFD','G14.5','G14.5HFD'),each=2))
rownames(ATAC.si.tab) <- ATAC.si.tab$SampleName
ATAC.si.tab$Time <- factor(ATAC.si.tab$Time,levels = c('Ctrl','G0HFD','G14.5','G14.5HFD'))

ATAC.mean.df <- round(do.call(cbind,by(t(as.data.frame(ATAC.merge.bed)[,ATAC.prefix]),ATAC.si.tab$Time,colMeans)),3)
head(ATAC.mean.df)
boxplot(ATAC.mean.df,outline = F)

for(i in colnames(ATAC.mean.df)){
  mcols(ATAC.merge.bed)[,paste(i,"mean",sep = ".")] <- round(ATAC.mean.df[,i],3)
  print(i)
}

ATAC.merge.bed$G0HFD.H3K27ac.nor.mean <- round((ATAC.merge.bed$Ins1G0HFD_H3K27ac_2_220924.nor+ATAC.merge.bed$Ins1G0HFD_H3K27ac_3_220924.nor)/2,3)
ATAC.merge.bed$CD.H3K27ac.nor.mean <- round((ATAC.merge.bed$Ins1G14.5CD_H3K27ac_221114_1.nor+ATAC.merge.bed$Ins1G14.5CD_H3K27ac_221114_2.nor)/2,3)
ATAC.merge.bed$HFD.H3K27ac.nor.mean <- round((ATAC.merge.bed$Ins1G14.5Acss2WT_HFD_30_H3K27ac_220802_2.nor+ATAC.merge.bed$Ins1G14.5Acss2OEWT_HFD_N275_20220905.nor)/2,3)    


ATAC.merge.bed$G0HFD.H3K27ac.mean <- round((ATAC.merge.bed$Ins1G0HFD_H3K27ac_2_220924+ATAC.merge.bed$Ins1G0HFD_H3K27ac_3_220924)/2,3)
ATAC.merge.bed$CD.H3K27ac.mean <- round((ATAC.merge.bed$Ins1G14.5CD_H3K27ac_221114_1+ATAC.merge.bed$Ins1G14.5CD_H3K27ac_221114_2)/2,3)
ATAC.merge.bed$HFD.H3K27ac.mean <- round((ATAC.merge.bed$Ins1G14.5Acss2WT_HFD_30_H3K27ac_220802_2+ATAC.merge.bed$Ins1G14.5Acss2OEWT_HFD_N275_20220905)/2,3)    



MyLogBoxPlot(as.data.frame(ATAC.merge.bed)[,ATAC.prefix],min.cutoff = 1)
######annotate######
KO.DEG.list <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20221203/DEG/KO.DEG.list.rds')
j <- 'G14.5_G14.5_HFD'

KO.DEG.list.filter <- list()
for(j in c(
           'G14.5_G14.5_HFD'#0.01 1.5
           )){
  KO.DEG.list.filter[[j]] <-  KO.DEG.list[[j]][KO.DEG.list[[j]]$p_val_adj<=0.01&
                                                 abs( KO.DEG.list[[j]]$avg_logFC) >= log(1.5) &
                                                 ( KO.DEG.list[[j]]$pct.1>=0.3 |  KO.DEG.list[[j]]$pct.2>=0.3),]
  
}


G14.5CDG14.5HFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$G14.5_G14.5_HFD$gene[KO.DEG.list.filter$G14.5_G14.5_HFD$cluster=='G14.5_HFD'],1]
mcols(ATAC.merge.bed)$G14.5HFD <- '-'
mcols(ATAC.merge.bed[ATAC.merge.bed$gene.id %in% G14.5CDG14.5HFD.up.gene])$G14.5HFD <- 'up'

G14.5HFD.up.bed <- ATAC.merge.bed[ATAC.merge.bed$G14.5HFD=='up']


# G14.5HFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$Virgin_B6G14.5_HFD$gene[KO.DEG.list.filter$Virgin_B6G14.5_HFD$cluster=='B6G14.5_HFD'],1]
# G14.5CD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$Virgin_G14.5$gene[KO.DEG.list.filter$Virgin_G14.5$cluster=='G14.5'],1]
# G0HFDG14.5HFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$Virgin_HFD_B6G14.5_HFD$gene[KO.DEG.list.filter$Virgin_HFD_B6G14.5_HFD$cluster=='B6G14.5_HFD'],1]

G14.5CDHFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$G14.5_G14.5_HFD$gene[KO.DEG.list.filter$G14.5_G14.5_HFD$cluster=='G14.5_HFD'],1]

length(G14.5HFD.up.gene)#374
length(G14.5CD.up.gene)#277
length(G0HFDG14.5HFD.up.gene)
length(G14.5CDHFD.up.gene)#412

gene.venn <- Venn(list('G0CD&G14.5HFD' = G14.5HFD.up.gene,
          'G0HFD&G14.5HFD' = G0HFDG14.5HFD.up.gene
          ))
plot(gene.venn)
G14.5.DEG.list <- list()
G14.5.DEG.list[['G0HFDG14.5HFD']] <- G0HFDG14.5HFD.up.gene
G14.5.DEG.list[['G0CDG14.5HFD']] <- G14.5HFD.up.gene
G14.5.DEG.list[['G0CDG14.5CD']] <- G14.5CD.up.gene
sapply(G14.5.DEG.list, length)


atac.o.gene.venn <- Venn(list('G0CDG14.5'= G14.5CD.up.gene[G14.5CD.up.gene %in% cluster1.bed$gene.id],
          'G0HFDG14.5HFD' = G0HFDG14.5HFD.up.gene[G0HFDG14.5HFD.up.gene %in% cluster1.bed$gene.id]))
plot(atac.o.gene.venn)


tmp.venn.list <- list()
tmp.bed.list <- list()
for(tmp.cluster in c('G0HFDG14.5HFD','G0CDG14.5HFD')){
  
  tmp.venn.list[[tmp.cluster]] <-  Venn(list('G14.5CD' = G14.5.DEG.list[['G0CDG14.5CD']],
            'G14.5HFD' = G14.5.DEG.list[[tmp.cluster]]
            ))
  tmp.bed.list[[tmp.cluster]] <- list()
  tmp.bed.list[[tmp.cluster]][['G14.5CDG14.5HFD']] <- ATAC.merge.bed[ATAC.merge.bed$gene.id %in% tmp.venn.list[[tmp.cluster]]@IntersectionSets$`11`]
  tmp.bed.list[[tmp.cluster]][['G14.5HFD']] <- ATAC.merge.bed[ATAC.merge.bed$gene.id %in% tmp.venn.list[[tmp.cluster]]@IntersectionSets$`01`]

  # tmp.bed.list[[tmp.cluster]][['cluster1.G14.5CDG14.5HFD']] <- cluster1.bed[cluster1.bed$gene.id %in% tmp.venn.list[[tmp.cluster]]@IntersectionSets$`11`]
  # tmp.bed.list[[tmp.cluster]][['cluster1.G14.5HFD']] <- cluster1.bed[cluster1.bed$gene.id %in% tmp.venn.list[[tmp.cluster]]@IntersectionSets$`01`]
}

pdf('Figure/ATAC.G0CDG14.5CDHFD.HFDCD.boxplot.pdf',3,5)
for(tmp.cluster in names(tmp.bed.list)[1]){
  for(tmp.con in names(tmp.bed.list[[tmp.cluster]])){
    print(tmp.con)
    boxplot(as.data.frame(tmp.bed.list[[tmp.cluster]][[tmp.con]])[,c('G14.5.mean','G14.5HFD.mean')],
            main = paste('ATAC',tmp.con,sep= '_'),outline=F,col=color.list$ATAC[['col']][c(50)],ylab='signal intensity',lwd=3,notch=T,
            ylim = c(0,26),
            xaxt="n")
    box(lwd=3)
    axis(1,at = 1:2,labels = c('CD','HFD'),cex.axis = 0.5)
    P300.p <- wilcox.test(as.data.frame(tmp.bed.list[[tmp.cluster]][[tmp.con]])[,c('G14.5.mean')],
                          as.data.frame(tmp.bed.list[[tmp.cluster]][[tmp.con]])[,c('G14.5HFD.mean')],
                          paired = T)
    
    MyText(paste("p-value: \n",
                 P300.p$p.value,"\n"),
           text.cex = 0.5)
  }
}
dev.off()


pdf('Figure/ATAC.413.G14.5G14.5HFD.box.pdf',3,5)
boxplot(as.data.frame(ATAC.merge500.bed)[,c('Ctrl.mean','G0HFD.mean','G14.5.mean','G14.5HFD.mean')],
        main =  paste('ATAC','G0CDG14.5CD',sep= '_'),outline=F,col='gray60',ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,25),
        xaxt="n")
box(lwd=3)
axis(1,at = 1:4,labels = c('G0CD','G0HFD','G14.5CD','G14.5HFD'),cex.axis = 1)

P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed)[,c('G14.5.mean')],
                      as.data.frame(G14.5HFD.up.bed)[,c('G14.5HFD.mean')],
                      paired = T)

MyText(paste("p-value: \n",
             P300.p$p.value,"\n"),
       text.cex = 0.5)
dev.off()

G14.5CD.up.bed <- ATAC.merge.bed[ ATAC.merge.bed$gene.id %in% G14.5CD.up.gene]
length(G14.5CD.up.bed)

pdf('Figure/H3K27ac.B630.413.G14.5G14.5HFD.box.pdf',3,5)
boxplot(as.data.frame(G14.5CD.up.bed)[,c('CD.H3K27ac.nor.mean','HFD.H3K27ac.nor.mean')],
        main =  paste('H3K27ac','G0G14.5CD',sep= '_'),outline=F,col='gray60',ylab='signal intensity',lwd=3,notch=T,
        # ylim = c(0,26),
        xaxt="n")
box(lwd=3)

axis(1,at = 1:2,labels = c('CD','HFD'),cex.axis = 0.5)

P300.p <- wilcox.test(as.data.frame(G14.5CD.up.bed)[,c('CD.H3K27ac.nor.mean')],
                      as.data.frame(G14.5CD.up.bed)[,c('HFD.H3K27ac.nor.mean')],
                      paired = T)

MyText(paste("p-value: \n",
             P300.p$p.value,"\n"),
       text.cex = 0.5)
dev.off()

cluster1.413.venn <- Venn(list('413GHFDup.genes' = G14.5CDG14.5HFD.up.gene,
          'cluster2' = unique(list.atac.element.bed$cluster3$gene.id)))

plot(cluster1.413.venn)


pdf('Figure/H3K27ac.B630.atac.2cluster.pdf',3,5)
for(tmp.cluster in names(list.atac.element.bed)[1:2]){
  for(his in c('H3K27ac')){
    boxplot(list.nor.mean.df[[his]][list.atac.element.bed[[tmp.cluster]]$rowname,],main= paste(tmp.cluster,his,sep = '.'),outline=F,col=color.list[[his]][['col']][c(50)],ylab='signal intensity',lwd=3,notch=T,
            ylim = c(0,20)
    )
    box(lwd=3)
    
    P300.p <- wilcox.test(list.nor.mean.df[[his]][list.atac.element.bed[[tmp.cluster]]$rowname,][,1],
                          list.nor.mean.df[[his]][list.atac.element.bed[[tmp.cluster]]$rowname,][,2],
                          paired = T)
    
    MyText(paste("p-value: \n",
                 P300.p$p.value,"\n"),
           text.cex = 0.5)
    
  }
}

dev.off()


pdf('Figure/ATAC.overall.G14.5G14.5HFD.box.pdf',3,5)
boxplot(as.data.frame(ATAC.merge500.bed)[,c('Ctrl.mean','G0HFD.mean','G14.5.mean','G14.5HFD.mean')],
        main =  paste('ATAC','G0CDG14.5CD',sep= '_'),outline=F,col='gray60',ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,25),
        xaxt="n")
box(lwd=3)
axis(1,at = 1:4,labels = c('G0CD','G0HFD','G14.5CD','G14.5HFD'),cex.axis = 1)
# 
# P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed)[,c('G14.5.mean')],
#                       as.data.frame(G14.5HFD.up.bed)[,c('G14.5HFD.mean')],
#                       paired = T)
# 
# MyText(paste("p-value: \n",
#              P300.p$p.value,"\n"),
#        text.cex = 0.5)
dev.off()

pdf('G:/lab/Article/heshuang/BYLW/chromatin/Figs10.ATAC.overall.G14.5G14.5HFD.box.pdf',3,5)
boxplot(as.data.frame(ATAC.merge500.bed)[,c('Ctrl.mean','G0HFD.mean','G14.5.mean','G14.5HFD.mean')],
        main =  paste('ATAC','G0CDG14.5CD',sep= '_'),outline=F,col=ko.time.col[c('Virgin','Virgin_HFD','G14.5','G14.5_HFD')],ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,25),
        xaxt="n")
box(lwd=3)
axis(1,at = 1:4,labels = c('G0CD','G0HFD','G14.5CD','G14.5HFD'),cex.axis = 1)
# 
# P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed)[,c('G14.5.mean')],
#                       as.data.frame(G14.5HFD.up.bed)[,c('G14.5HFD.mean')],
#                       paired = T)
# 
# MyText(paste("p-value: \n",
#              P300.p$p.value,"\n"),
#        text.cex = 0.5)
dev.off()


pdf('Figure/overall.H3K27ac.box.pdf',3,5)
boxplot(as.data.frame(ATAC.merge500.bed)[,c('G0HFD.H3K27ac.nor.mean','CD.H3K27ac.nor.mean','HFD.H3K27ac.nor.mean')],
        main =  'H3K27ac',outline=F,col='gray60',ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,22),
        xaxt="n")
box(lwd=3)

# axis(1,at = 1:2,labels = c('CD','HFD'),cex.axis = 0.5)
# 
# P300.p <- wilcox.test(as.data.frame(G14.5CD.up.bed)[,c('CD.H3K27ac.nor.mean')],
#                       as.data.frame(G14.5CD.up.bed)[,c('HFD.H3K27ac.nor.mean')],
#                       paired = T)
# 
# MyText(paste("p-value: \n",
#              P300.p$p.value,"\n"),
#        text.cex = 0.5)
dev.off()

################
rownames(genes.inf.input) <- genes.inf.input$EnsemblGeneID
gene.parts <-  MyReadTranscriptFeatures("G:/lab/genome/mm10/mm10_convert_ens_filter.refFlat.bed",
                                        up.flank = 2000,
                                        down.flank = 2000)
ATAC.merge.bed.anno <- annotateWithGeneParts( ATAC.merge.bed,
                                              gene.parts)
mcols(ATAC.merge.bed)$trans.id <-  ATAC.merge.bed.anno@dist.to.TSS$feature.name
mcols(ATAC.merge.bed)$gene.id <- mcols(tss.bed[mcols( ATAC.merge.bed)$trans.id])$gene.id
names(ATAC.merge.bed) <- paste("V",1:length( ATAC.merge.bed$trans.id),sep='')
mcols(ATAC.merge.bed)$Symbol <- genes.inf.input[mcols( ATAC.merge.bed)$gene.id,"Symbol"]
mcols(ATAC.merge.bed)$TF <- genes.inf.input[mcols( ATAC.merge.bed)$gene.id,"TF"]
mcols(ATAC.merge.bed)$promoter <-  ATAC.merge.bed.anno@members[,1]
mcols(ATAC.merge.bed)$exon <-  ATAC.merge.bed.anno@members[,2]
mcols(ATAC.merge.bed)$intron <-  ATAC.merge.bed.anno@members[,3]
mcols(ATAC.merge.bed)$distalregion <- ifelse(abs( ATAC.merge.bed.anno@dist.to.TSS$dist.to.feature) > 3000,1,0)
mcols(ATAC.merge.bed)$rowname <- names(ATAC.merge.bed)
table(mcols(ATAC.merge.bed)$promoter)
table(mcols(ATAC.merge.bed)$distalregion)

mcols(ATAC.merge.bed)$G14.5HFD <- '-'
mcols(ATAC.merge.bed[ATAC.merge.bed$gene.id %in% G14.5CDHFD.up.gene])$G14.5HFD <- 'up'
########
# G14.5HFD.gene <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20220913/G0G14.5HMKO/addG18.5/final/overlap/gene_inf/v2.c1_HFDgene.tab')
# G14.5HFD.bed <- ATAC.merge.bed[ATAC.merge.bed$gene.id %in% G14.5HFD.gene$EnsemblGeneID]
# length(G14.5HFD.bed)#437
# pdf('G14.5HFD.gene.pdf',3,5)
# boxplot(ATAC.mean.df[G14.5HFD.bed$rowname,],outline=F,notch=T,main = 'G14.5HFD upregulated genes',ylab = 'ATAC signal intensity',lwd=3,col=atac.color.image[50])
# box(lwd=3)
# wilcox.test.res1 <- wilcox.test(ATAC.mean.df[G14.5HFD.bed$rowname,'Ctrl'],
#                                 ATAC.mean.df[G14.5HFD.bed$rowname,'G14.5'],
#                                 paired = T)
# MyText(paste("G0CDG14.5CDp-value: \n",
#              wilcox.test.res1$p.value,"\n"),
#        text.cex = 0.5)
# 
# wilcox.test.res1 <- wilcox.test(ATAC.mean.df[G14.5HFD.bed$rowname,'G14.5HFD'],
#                                 ATAC.mean.df[G14.5HFD.bed$rowname,'G14.5'],
#                                 paired = T)
# MyText(paste("G14.5CDCDG14.5HFDp-value: \n",
#              wilcox.test.res1$p.value,"\n"),
#        text.cex = 0.5)
# dev.off()

KO.DEG.list.filter <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20220913/G0G14.5HMKO/addG18.5/final/20221105/KO.DEG.list.filter.rds')
G14.5HFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$G14.5_B6G14.5_HFD$gene[KO.DEG.list.filter$G14.5_B6G14.5_HFD$cluster=='B6G14.5_HFD'],1]
G14.5HFD.up.bed <- ATAC.merge.bed[ATAC.merge.bed$gene.id %in% G14.5HFD.up.gene]
length(G14.5HFD.up.bed)#423

MyLogBoxPlot(as.data.frame(ATAC.merge.bed)[,ATAC.prefix],min.cutoff = 1)


pdf('G14.5CD_G14.5HFD.up.gene.pdf',3,5)
boxplot(log2(ATAC.mean.df[G14.5HFD.up.bed$rowname,]+1),outline=F,notch=T,main = 'G14.5HFD.up upregulated genes',ylab = 'ATAC signal intensity',lwd=3,col=add.alpha(time.colors[c(1,11,8,13)],0.8))
box(lwd=3)
wilcox.test.res1 <- wilcox.test(ATAC.mean.df[G14.5HFD.up.bed$rowname,'Ctrl'],
                                ATAC.mean.df[G14.5HFD.up.bed$rowname,'G14.5'],
                                paired = T)
MyText(paste("G0CDG14.5CDp-value: \n",
             wilcox.test.res1$p.value,"\n"),
       text.cex = 0.5)

wilcox.test.res1 <- wilcox.test(ATAC.mean.df[G14.5HFD.up.bed$rowname,'G14.5HFD'],
                                ATAC.mean.df[G14.5HFD.up.bed$rowname,'G14.5'],
                                paired = T)
MyText(paste("G14.5CDCDG14.5HFD.upp-value: \n",
             wilcox.test.res1$p.value,"\n"),
       text.cex = 0.5)
dev.off()
#########
names(list.ATAC.bed) <- c(names(list.ATAC.bed)[1:3],'G0HFD')

pdf("ATAC.state.hist.raw.pdf",
    8,7)
for(i in names(list.ATAC.bed)){
  print(i)
  mcols(ATAC.merge.bed)[,paste(i,".state",sep = "")] <- 0
  mcols(ATAC.merge.bed[countOverlaps(ATAC.merge.bed,list.ATAC.bed[[i]])>0])[,paste(i,".state",sep = "")] <- 1
  p <- ggplot(data = as.data.frame(ATAC.merge.bed),
              aes(x = log2(get(paste(i,"mean",sep = "."))+1),
                  group = factor(get(paste(i,".state",sep = ""))),
                  fill = factor(get(paste(i,".state",sep = "")))))+
    # scale_x_continuous(breaks=seq(min(log2(get(i))), 5, 0.25), limits=c(min(log2(get(i))), 5))+
    geom_density(adjust=1.5,alpha=0.5)+
    labs(title = i)
  plot(p)
}
dev.off()
############
ATAC.cutoff <- round(2**2-1,3)
for(tmp.cluster in names(list.ATAC.bed)){
  mcols(ATAC.merge.bed)[,paste(tmp.cluster,'state2',sep = '.')] <- 0
  mcols(ATAC.merge.bed[mcols(ATAC.merge.bed)[,paste(tmp.cluster,'mean',sep = '.')]>ATAC.cutoff])[,paste(tmp.cluster,'state2',sep = '.')] <- 1
}


ATAC.merge.bed$G0G14.5CDHFD.state <- paste(ATAC.merge.bed$Ctrl.state2,ATAC.merge.bed$G0HFD.state2,
                                           ATAC.merge.bed$G14.5.state2,ATAC.merge.bed$G14.5HFD.state2,sep=  '')
table(ATAC.merge.bed$G0G14.5CDHFD.state)

state.peak.heatmap <- as.data.frame(ATAC.merge.bed[!ATAC.merge.bed$G0G14.5CDHFD.state  %in% c('0000','1111')])[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'state2',sep = '.')]   
MyHeatmap(state.peak.heatmap,Colv = 'none')

keep.state <- c('1010','0101',#HFD
                '0001','1110',#G14.5HFD
                '1101',#G14.5
                '0010',
                '0111',
                '0011',#G14.5
               
                '1011'#G0HFD
                )
keep.state <- c(
  '0111','0011',
   '1011','0001', '0101','1101','1110','0010','1010', '1000','1100'
)

ATAC.merge.keep.bed <- ATAC.merge.bed[ATAC.merge.bed$G0G14.5CDHFD.state  %in% keep.state]
ATAC.merge.keep.bed$G0G14.5CDHFD.state <- factor(ATAC.merge.keep.bed$G0G14.5CDHFD.state,levels = keep.state)

state.peak.heatmap <- as.data.frame(ATAC.merge.keep.bed)[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'state2',sep = '.')]   
state.peak.heatmap <- state.peak.heatmap[ATAC.merge.keep.bed$rowname[order(ATAC.merge.keep.bed$G0G14.5CDHFD.state)],]
MyHeatmap(state.peak.heatmap[,],Colv = 'none',
          Rowv = 'none',
          type = 'raw',labCol = c('G0CD','G0HFD','G14.5CD','G14.5HFD'))
signal.peak.heatmap <- as.data.frame(ATAC.merge.keep.bed)[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'mean',sep = '.')]   
signal.peak.heatmap <- signal.peak.heatmap[ATAC.merge.keep.bed$rowname[order(ATAC.merge.keep.bed$G0G14.5CDHFD.state)],]
MyHeatmap(as.matrix(signal.peak.heatmap[,]),Colv = 'none',
          Rowv = 'none',
          type = 'log.raw',
          labCol = c('G0CD','G0HFD','G14.5CD','G14.5HFD'))

signal.peak.rep.heatmap <- as.data.frame(ATAC.merge.keep.bed)[,ATAC.prefix]   
signal.peak.rep.heatmap <- signal.peak.rep.heatmap[ATAC.merge.keep.bed$rowname[order(ATAC.merge.keep.bed$G0G14.5CDHFD.state)],]
MyHeatmap(as.matrix(signal.peak.rep.heatmap[,]),Colv = 'none',
          Rowv = 'none',
          type = 'log.raw'#,
          #labCol = c('G0CD','G0HFD','G14.5CD','G14.5HFD')
)

###########limma#######
ATAC.df2 <- MyGeneExp(as.data.frame(ATAC.merge.bed)[,ATAC.prefix[3:6]],ATAC.cutoff,2)
ATAC.limma.res <- Mylimma(ATAC.df2,
                          colnames(ATAC.df2)[1:2],
                          colnames(ATAC.df2)[3:4]
)

ATAC.limma.res$cluster <- '-'
ATAC.limma.res[ATAC.limma.res$logFC>=0,'cluster'] <- 'G14.5HFD'
ATAC.limma.res[ATAC.limma.res$logFC<=0,'cluster'] <- 'G14.5'
table(ATAC.limma.res$cluster)#
fc.tmp <- 2
ATAC.limma.res.filter <- ATAC.limma.res[ATAC.limma.res$P.Value<=0.08&abs(ATAC.limma.res$logFC)>=log2(fc.tmp),]
table(ATAC.limma.res.filter$cluster)
ATAC.limma.res$state <-  '-'
ATAC.limma.res[ATAC.up.bed$rowname,'state'] <- 'up'
ATAC.limma.res[ATAC.down.bed$rowname,'state'] <- 'down'
fc.tmp <- 2
pdf('Figure/limma.com.nor.ATAC.cutoff.fc4.pval0.1.peak.pdf')
par(oma = c(0,1,0,0))
MyLimmaPlot3(ATAC.limma.res,P.Val.cutoff = 0.08,fc.max = 20,fc.cutoff = log2(fc.tmp),xlim = c(0,60),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             axes=FALSE,
             xlab = c('signal intensity'),
             count = T
             #ylim = c(-6,6)
)


MyLimmaPlot3(ATAC.limma.res,P.Val.cutoff = 0.08,fc.max = 20,fc.cutoff = log2(fc.tmp),xlim = c(0,40),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             axes=FALSE,
             xlab = c('signal intensity'),
             count = F
             #ylim = c(-6,6)
)

axis(2,at=seq(-20,20,by=4),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(0,60,by=5),cex.axis = 1.5,lwd = 3)
box(lwd = 3)
dev.off()
########preg.gene######
ATAC.up.bed <- ATAC.merge.bed[rownames(ATAC.limma.res.filter)[ATAC.limma.res.filter$cluster=='G14.5HFD']]
length(ATAC.up.bed)#1984
table(ATAC.up.bed$G14.5CDHFD.state)

ATAC.down.bed <- ATAC.merge.bed[rownames(ATAC.limma.res.filter)[ATAC.limma.res.filter$cluster=='G14.5']]
length(ATAC.down.bed)#545
table(ATAC.down.bed$G14.5CDHFD.state)

library(Vennerable)
preg.gene.tab <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20220913/Glut2H/LMgene/gene_inf/preg.gene.si.1e27.addcc.tab')
rownames(preg.gene.tab) <- preg.gene.tab$EnsemblGeneID

G14.5HFD.upcluster1.venn <- Venn(list('G14.5HFD.ATAC.up' = ATAC.up.bed$gene.id,
                                      'cluster1' = preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='cc']))
plot(G14.5HFD.upcluster1.venn)
G14.5HFD.upcluster2.venn <- Venn(list('G14.5HFD.ATAC.up' = ATAC.up.bed$gene.id,
                                      'cluster2' = preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up']))
plot(G14.5HFD.upcluster2.venn)
G14.5HFD.upcluster3.venn <- Venn(list('G14.5HFD.ATAC.up' = ATAC.up.bed$gene.id,
                                      'cluster3' = preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.down']))
plot(G14.5HFD.upcluster3.venn)


cluster2.bed <- ATAC.merge.bed[ATAC.merge.bed$gene.id %in% preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up']]
pdf('G14.5HFD.gene.pdf',3,5)
boxplot(ATAC.mean.df[cluster2.bed$rowname,2:3],outline=F,notch=T,main = 'cluster3',ylab = 'ATAC signal intensity',lwd=3,col=atac.color.image[50])
box(lwd=3)
wilcox.test.res1 <- wilcox.test(ATAC.mean.df[cluster2.bed$rowname,'Ctrl'],
                                ATAC.mean.df[cluster2.bed$rowname,'G14.5'],
                                paired = T)

MyText(paste("G0CDG14.5CDp-value: \n",
             wilcox.test.res1$p.value,"\n"),
       text.cex = 0.5)

wilcox.test.res1 <- wilcox.test(ATAC.mean.df[cluster2.bed$rowname,'G14.5HFD'],
                                ATAC.mean.df[cluster2.bed$rowname,'G14.5'],
                                paired = T)

MyText(paste("G14.5CDCDG14.5HFDp-value: \n",
             wilcox.test.res1$p.value,"\n"),
       text.cex = 0.5)
dev.off()

########
rm(list.ATAC.bdg,list.ATACHFD.bdg)
rm(seu.B6.CDHFD)
rm(sm.up.list)
rm(sm.up.mean.cluster.list,sm.up.mean.list,sm.up.sort.list)
rm(list.Acss2HFD.H3K27ac.bdg)
saveRDS(ATAC.merge.bed,'ATAC.merge.bed.rds')
save.image('ATAC.RData')


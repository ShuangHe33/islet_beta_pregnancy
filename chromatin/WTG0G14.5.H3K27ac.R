
source('G:/pcatest/MyFunction.R')
source("G:/r script/file_function.R")

library(genomation)
library(GenomicRanges)
library(RColorBrewer)
library(VennDiagram)
library(Vennerable)
library(ggplot2)

########load data#########
# list.G14.5Acss2HFD.consensuspeak <- readRDS("G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/G14.5Acss2HMKO/nbl.score.bed/list.G14.5Acss2HFD.H3K27ac.c5_8_c1.ConsensusPeaks.rds")
list.Acss2HFD.H3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/Acss2KO/HFD/new/list.bdg.Acss2HMKO.HFD.H3K27ac.rds')
list.Acss2HFD.H3K27ac.bdg$Ins1G14.5Acss2OEWT_HFD_N286_20220905 <- NULL

list.G14.5HFD.H3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/HFD/list.H3K27ac.HFD.bdg.rds')
names(list.G14.5HFD.H3K27ac.bdg)[1] <- 'Ins1G14.5HFD_H3K27ac_220620'
list.G14.5Acss2CDHFD.consensuspeak <- readRDS('G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/nbl.score.bed/list.G14.5HFD.H3K27ac.c5_8.ConsensusPeaks.rds')
list.G14.5Acss2HFD.consensuspeak <- readRDS("G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/G14.5Acss2HMKO/nbl.score.bed/list.G14.5Acss2HFD.H3K27ac.c5_8_c1.ConsensusPeaks.rds")


merge.bed <- reduce(c(list.G14.5Acss2CDHFD.consensuspeak$G0.c1,
                      list.G14.5Acss2HFD.consensuspeak$G14.5Acss2WTHFD.c1,
                      list.G14.5Acss2CDHFD.consensuspeak$G14.5CD.c1,
                      list.G14.5Acss2HFD.consensuspeak$G14.5Acss2HMKOHFD.c1
))
merge.bed <- merge.bed[!merge.bed@seqnames %in% c('chrM','chrY')]
length(merge.bed)#27704

list.G14.5HFD.H3K27ac.bdg <- c(list.G14.5HFD.H3K27ac.bdg[c(4:7)],list.Acss2HFD.H3K27ac.bdg[c(3,4,1,2)])


prefix.list <- list()
prefix.list[['G14.5HFD.H3K27ac']] <- names(list.G14.5HFD.H3K27ac.bdg)

for(type in c(
  'G14.5HFD.H3K27ac'
)){
  for(i in names(get(paste('list',type,'bdg',sep = '.')))){
    print(i)
    mcols(merge.bed)[,i] <- round(MyCalTd(merge.bed,get(paste('list',type,'bdg',sep = '.'))[[i]],"td"),3)
    mcols(merge.bed)[,i][mcols(merge.bed)[,i]<0] <- 0
  }
}

#########
dir.create('QC')
pdf('QC/hist.pdf')
for(prefix in names(prefix.list)){
  for(sn.tmp in prefix.list[[prefix]]){
    hist(log2(mcols(merge.bed)[,sn.tmp]),main = sn.tmp,breaks = 100)
  }
  MyLogBoxPlot(as.data.frame(merge.bed)[,prefix.list[[prefix]]],main = prefix,min.cutoff = 1)
}
dev.off()
####normalization at p#####
promoter.bed <- readRDS('../../20221011/promoter.bed.rds')
##########20221203####
promoter.bed <- readRDS('G:/project/pregnant_mouse/chromatin/promoter/20221209/promoter.bed.td.rds')
for(his in c(
  'G14.5HFD.H3K27ac'
  
)){
  for(i in names(get(paste('list',his,'bdg',sep = '.')))){
    print(i)
    mcols(promoter.bed)[,i] <- round(MyCalTd(promoter.bed,get(paste('list',his,'bdg',sep = '.'))[[i]],"td"),3)
    mcols(promoter.bed)[,i][mcols(promoter.bed)[,i]<0] <- 0
  }
}
saveRDS(promoter.bed,'promoter.bed.rds')

#########
dir.create('QC')
#promoter.bed <- readRDS('G:/project/pregnant_mouse/chromatin/H3K27ac/HFD/20220925_4/promoter.bed.rds')
list.close <- list()
for(his in c(
  "G14.5HFD.H3K27ac")){
  list.close[[his]] <- list()
  pdf(paste("QC/",his,".p.hist.according to atacpeak.pdf",sep = ''))
  for (tmp.time in prefix.list[[his]]){
    print(tmp.time)
    # con <- ifelse(grepl('CD',tmp.time),'G14.5','G14.5HFD')
    if(grepl('G14.5CD',tmp.time)){con <- 'G14.5'}
    if(grepl('G14.5Acss2',tmp.time)){con <- 'G14.5HFD'}
    if(grepl('G0HFD',tmp.time)){con <- 'G0HFD'}
    tmp.atac.0.bed <- promoter.bed[mcols(promoter.bed)[,paste(con,'.ATAC.peak.state',sep= '')]==0]
    tmp.atac.1.bed <- promoter.bed[mcols(promoter.bed)[,paste(con,'.ATAC.peak.state',sep= '')]==1]
    list.close[[his]][[tmp.time]] <- mcols(tmp.atac.0.bed)[,tmp.time]
    hist(log2(mcols(tmp.atac.0.bed)[,tmp.time]),
         breaks = seq(-100,100,0.1),
         xlim=c(-8,8),
         #ylim = c(0,350)
         main = paste(tmp.time,"atac.0",sep = "_")
    )
    
    abline(v = -1,col = "red")
    abline(v=-0.5,col = "red")
    abline(v=0.5,col = "red")
    abline(v=1,col = "red")
    hist(log2(mcols(tmp.atac.1.bed)[,tmp.time]),
         breaks = seq(-100,100,0.1),
         xlim=c(-8,8),
         #ylim = c(0,350)
         main = paste(tmp.time,"atac.1",sep = "_")
    )
    
    abline(v = -1,col = "red")
    abline(v=-0.5,col = "red")
    abline(v=0.5,col = "red")
    abline(v=1,col = "red")
    
  }
  dev.off()
}
boxplot(list.close$G14.5HFD.H3K27ac,outline=F,main=his)

###########select##########
list.sf.factor <- list()
# list.close <- readRDS('../list.close.rds')
# list.close2 <- readRDS('../B6G0HFDG14.5CDHFD/list.close.rds')
# list.close$G14.5CDHFD.H3K27ac[[names(list.close2$G14.5HFD.H3K27ac)[1]]] <- list.close2$G14.5HFD.H3K27ac$Ins1G0HFD_H3K27ac_2_220924
# list.close$G14.5CDHFD.H3K27ac[[names(list.close2$G14.5HFD.H3K27ac)[2]]] <- list.close2$G14.5HFD.H3K27ac$Ins1G0HFD_H3K27ac_3_220924

boxplot(list.close$G14.5CDHFD.H3K27ac[prefix.list$G14.5HFD.H3K27ac],outline=F)

for(his in c('G14.5HFD.H3K27ac'
)){
  list.sf.factor[[his]] <- list()
  for(sn.tmp in prefix.list[[his]]){
    print(sn.tmp)
    list.sf.factor[[his]][[sn.tmp]] <- mean(list.close$G14.5HFD.H3K27ac[[sn.tmp]])/mean(unlist(list.close$G14.5HFD.H3K27ac[prefix.list[[his]]]))
  }
}
saveRDS(list.sf.factor,'list.mean.sf.factor.20221203.rds')
unlist(list.sf.factor[[his]])
############
list.nor.df <- list()
pdf('QC/log.nor.pdf')
for(prefix in c(
  'G14.5HFD.H3K27ac'
  
)){
  print(prefix)
  list.nor.df[[prefix]] <- t(t(as.data.frame(merge.bed)[,prefix.list[[prefix]]])/unlist(list.sf.factor[[prefix]]))
  rownames(list.nor.df[[prefix]]) <- names(merge.bed)
  for(sn.tmp in colnames(list.nor.df[[prefix]])){
    mcols(merge.bed)[,paste(sn.tmp,'.nor',sep ='')] <- round(list.nor.df[[prefix]][,sn.tmp],3)
  }
  MyLogBoxPlot(list.nor.df[[prefix]],main=prefix,min.cutoff = 1,log = T)
}
dev.off()
#########save data#######
saveRDS(merge.bed,'merge.bed.rds')
##########annotation########
KO.DEG.list <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20220913/G0G14.5HMKO/addG18.5/final/20221105/KO.DEG.list.rds')
j <- 'G14.5_B6G14.5_HFD'
KO.DEG.list.filter <- list()
KO.DEG.list.filter[[j]] <-  KO.DEG.list[[j]][KO.DEG.list[[j]]$p_val_adj<=0.01&
                                               abs( KO.DEG.list[[j]]$avg_logFC) >= log(1.5) &
                                               ( KO.DEG.list[[j]]$pct.1>=0.3 |  KO.DEG.list[[j]]$pct.2>=0.3),]
G14.5HFD.up.gene <- genes.inf.input[sub('_','-',genes.inf.input$SymbolDedu) %in%  KO.DEG.list.filter$G14.5_B6G14.5_HFD$gene[KO.DEG.list.filter$G14.5_B6G14.5_HFD$cluster=='B6G14.5_HFD'],1]

gene.parts <-  MyReadTranscriptFeatures("G:/lab/genome/mm10/mm10_convert_ens_filter.refFlat.bed",
                                        up.flank = 2000,
                                        down.flank = 2000)
merge.bed <- merge.bed[!merge.bed@seqnames %in% c("chrM","chrY")]
merge.bed.anno <- annotateWithGeneParts(merge.bed,
                                        gene.parts)
mcols(merge.bed)$trans.id <- merge.bed.anno@dist.to.TSS$feature.name
mcols(merge.bed)$gene.id <- mcols(tss.bed[mcols(merge.bed)$trans.id])$gene.id
names(merge.bed) <- paste("V",1:length(merge.bed$trans.id),sep='')
mcols(merge.bed)$Symbol <- genes.inf.input[mcols(merge.bed)$gene.id,"Symbol"]
mcols(merge.bed)$TF <- genes.inf.input[mcols(merge.bed)$gene.id,"TF"]
mcols(merge.bed)$promoter <- merge.bed.anno@members[,1]
mcols(merge.bed)$exon <- merge.bed.anno@members[,2]
mcols(merge.bed)$intron <- merge.bed.anno@members[,3]
mcols(merge.bed)$distalregion <- ifelse(abs(merge.bed.anno@dist.to.TSS$dist.to.feature) > 3000,1,0)
mcols(merge.bed)$rowname <- names(merge.bed)
mcols(merge.bed)$G14.5HFD <- '-'
mcols(merge.bed[merge.bed$gene.id %in% G14.5HFD.up.gene])$G14.5HFD <- 'up'

mcols(merge.bed)$promoter2 <-  0
mcols(merge.bed)[countOverlaps(merge.bed,promoter.bed)>0,'promoter2'] <- 1
table(mcols(merge.bed)$promoter2)

# merge.bed$enhancer <- 0
# mcols(merge.bed[merge.bed$distalregion==1&(merge.bed$G14.5HFD.G14.5CDHFD.H3K4me1.nor.state==1|
#                                              merge.bed$G14.5CD.G14.5CDHFD.H3K4me1.nor.state==1)])[,'enhancer'] <- 1
# 
merge.bed <- merge.bed[merge.bed@ranges@width>=500]
G14.5HFD.bed <- merge.bed[merge.bed$G14.5HFD=='up']
boxplot(as.data.frame(G14.5HFD.bed)[,paste(prefix.list$G14.5HFD.H3K27ac,'.nor',sep = '')],outline=F)
boxplot(as.data.frame(G14.5HFD.bed)[,paste(c('G0HFD','G14.5CD','G14.5HFD'),his,'nor.mean',sep = '.')],outline=F,notch=T)
par(oma= c(8,0,0,0))
MyHeatmap(as.data.frame(G14.5HFD.bed)[,paste(prefix.list$G14.5HFD.H3K27ac[c(1:2,6:7,3:5)],'.nor',sep = '')],Colv = 'none')

########
merge.bed <- merge.bed[merge.bed@ranges@width>=500]
G14.5HFD.bed <- merge.bed[merge.bed$G14.5HFD=='up']

merge.bed$G0HFD.nor.mean <- round((merge.bed$Ins1G0HFD_H3K27ac_2_220924.nor+merge.bed$Ins1G0HFD_H3K27ac_3_220924.nor)/2,3)
merge.bed$G14.5CD.nor.mean <- round((merge.bed$Ins1G14.5CD_H3K27ac_221114_1.nor+merge.bed$Ins1G14.5CD_H3K27ac_221114_2.nor)/2,3)
merge.bed$G14.5HFD.nor.mean <- round((merge.bed$Ins1G14.5Acss2WT_HFD_30_H3K27ac_220802_2+merge.bed$Ins1G14.5Acss2OEWT_HFD_N275_20220905)/2,3)
merge.bed$G14.5HMKOHFD.nor.mean <- round((merge.bed$Ins1G14.5Acss2HMKO_HFD_7_H3K27ac_220802_2.nor+merge.bed$Ins1G14.5Acss2HMKO_HFD_99658_20220905.nor)/2,3)

boxplot(as.data.frame(merge.bed)[,paste(c('G0HFD','G14.5CD','G14.5HFD','G14.5HMKOHFD'),'nor.mean',sep = '.')],outline=F,notch=T)


##########cutoff######
split.bed.list <- list()
split.bed.list[["G14.5HFD.H3K27ac"]] <- list.G14.5Acss2CDHFD.consensuspeak[c("G0.c1","G14.5CD.c1",'c1')]
names(split.bed.list[["G14.5HFD.H3K27ac"]]) <- c('G0HFD','G14.5CD','G14.5HFD')

for(type in names(split.bed.list)){
  for(i in names(split.bed.list[[type]])){
    print(i)
    mcols(merge.bed)[,paste(i,type,"state",sep = ".")] <- 0
    mcols(merge.bed)[countOverlaps(merge.bed,split.bed.list[[type]][[i]]) > 0,paste(i,type,"state",sep = ".")] <- 1
  }
}

pdf("QC/state.hist.nor.pdf",
    8,7)
for(type in names(split.bed.list)){
  print(type)
  for(i in names(split.bed.list[[type]])){
    print(i)
    for(sn.tmp in grep(i,prefix.list[[type]],value = T)){
      print(sn.tmp)
      p <- ggplot(data = as.data.frame(merge.bed),
                  aes(x = log2(get(paste(sn.tmp,"nor",sep = "."))+1),
                      group = factor(get(paste(i,type,"state",sep = "."))),
                      fill = factor(get(paste(i,type,"state",sep = ".")))))+
        # scale_x_continuous(breaks=seq(min(log2(get(i))), 5, 0.25), limits=c(min(log2(get(i))), 5))+
        geom_density(adjust=1.5,alpha=0.5)+
        labs(title = paste(sn.tmp,i,type,"state",sep = '_'))
      plot(p)
    }
  }
}
dev.off()

################
cutoff.list <- list()
cutoff.list[["G14.5HFD.H3K27ac"]] <- round(2**2-1,3)

for(his in names(prefix.list)){
  print(his)
  for(time in c('G0HFD','G14.5CD'#,
                #'G14.5HFD'
  )){
    print(time)
    mcols(merge.bed)[,paste(time,his,'nor.mean',sep='.')] <- round((mcols(merge.bed)[,paste(prefix.list[[his]][grep(time,prefix.list[[his]],ignore.case = T)][1],'nor',sep = '.')]+mcols(merge.bed)[,paste(prefix.list[[his]][grep(time,prefix.list[[his]],ignore.case = T)][2],'nor',sep = '.')])/2,3)
  }
}


time <- 'G14.5HFD'
mcols(merge.bed)[,paste(time,his,'nor.mean',sep='.')] <- round((mcols(merge.bed)[,paste(prefix.list[[his]][grep(time,prefix.list[[his]],ignore.case = T)][2],'nor',sep = '.')]+mcols(merge.bed)[,paste(prefix.list[[his]][grep(time,prefix.list[[his]],ignore.case = T)][3],'nor',sep = '.')])/2,3)
boxplot(as.data.frame(merge.bed)[,paste(prefix.list$G14.5HFD.H3K27ac,'.nor',sep = '')])

boxplot(as.data.frame(merge.bed)[,paste(c('G0HFD','G14.5CD','G14.5HFD'),his,'nor.mean',sep = '.')])

for(his in c("G14.5HFD.H3K27ac"
)){
  for(time in c('G0HFD','G14.5CD',
                'G14.5HFD','G14.5HMKOHFD'
                
  )){
    print(time)
    mcols(merge.bed)[,paste(time,'nor.state',sep='.')] <- 0
    mcols(merge.bed[mcols(merge.bed)[,paste(time,'nor.mean',sep='.')]>=cutoff.list[[his]]])[,paste(time,'nor.state',sep='.')] <- 1
  }
}


####
merge.bed$G14.5CDHFD.H3K4me1.state <- paste(merge.bed$G14.5CD.G14.5CDHFD.H3K4me1.nor.state,merge.bed$G14.5HFD.G14.5CDHFD.H3K4me1.nor.state,sep = '')
table(merge.bed$G14.5CDHFD.H3K4me1.state)
table(merge.bed2$G14.5CDHFD.H3K4me1.state)

merge.bed$G14.5CDHFD.state <- paste(merge.bed$G14.5CD.G14.5CDHFD.nor.state,merge.bed$G14.5HFD.G14.5CDHFD.nor.state,sep = '')
table(merge.bed$G14.5CDHFD.state)
table(merge.bed2$G14.5CDHFD.state)



merge.bed$G14.5CDHFD.state <- paste(merge.bed$G14.5CD.G14.5CDHFD.nor.state,merge.bed$WT.Acss2HFD.nor.state,sep = '')
table(merge.bed$G14.5CDHFD.state)
table(merge.bed2$G14.5CDHFD.state)



merge.bed$G14.5CDHFD.peak.state <- paste(merge.bed$G14.5CD.c1.G14.5CDHFD.H3K27ac.state,merge.bed$c1.G14.5CDHFD.H3K27ac.state,sep = '')
table(merge.bed$G14.5CDHFD.peak.state)
table(merge.bed2$G14.5CDHFD.peak.state)

merge.bed$G14.5WTHMKO.state <- paste(merge.bed$WT.Acss2HFD.nor.state,merge.bed$HMKO.Acss2HFD.nor.state,sep = '')
merge.bed$G14.5WTHMKO.peak.state <- paste(merge.bed$G14.5Acss2WTHFD.c1.Acss2HFD.state,merge.bed$G14.5Acss2HMKOHFD.c1.Acss2HFD.state,sep = '')
table(merge.bed$G14.5WTHMKO.state)
table(merge.bed$G14.5WTHMKO.peak.state)

table(merge.bed2$G14.5WTHMKO.state)
table(merge.bed2$G14.5WTHMKO.peak.state)

hist(merge.bed[merge.bed$G14.5CDHFD.H3K4me1.state=='01']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5CDHFD.H3K4me1.state=='10']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5WTHMKO.peak.state=='01']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5WTHMKO.peak.state=='10']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5WTHMKO.state=='01']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5WTHMKO.state=='10']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5CDHFD.peak.state=='01']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5CDHFD.peak.state=='10']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5CDHFD.state=='01']@ranges@width)
abline(v=500)
hist(merge.bed[merge.bed$G14.5CDHFD.state=='10']@ranges@width)
abline(v=500)
hist(merge.bed@ranges@width,breaks=500)
abline(v=500)
##########QC of cutoff######
#######width over 500bp######
saveRDS(merge.bed,'merge.bed.rds')
merge.bed <- merge.bed[merge.bed@ranges@width>=500]

list.H3K27ac.element.bed <- list()
for(time in c(#'G14.5CDHFD.state',
  'G14.5CDHFD.H3K4me1.state'#,
  #'G14.5WTHMKO.state'
)){
  print(time)
  for(tmp.cluster in c('01','10')){
    print(tmp.cluster)
    list.H3K27ac.element.bed[[paste(time,tmp.cluster,sep = '.')]] <- merge.bed[mcols(merge.bed)[,time]==tmp.cluster]
  }
}
sapply(list.H3K27ac.element.bed, length)
########
pdf('QC/over500bp.nor.log.hist.pdf')
for(prefix in names(prefix.list)[c(1,3)]){
  for(sn.tmp in prefix.list[[prefix]]){
    hist(log2(mcols(merge.bed)[,sn.tmp]),main = sn.tmp,breaks = 100)
  }
  MyLogBoxPlot(as.data.frame(merge.bed)[,paste(prefix.list[[prefix]],'nor',sep = '.')],main = prefix,min.cutoff = 1.5)
}
dev.off()
########G14.5CD_HFD######
G14.5CD_HFD.limma.list <- list()
split.nor.df.filter.list <- list()
for(type in c(
  "G14.5HFD.H3K27ac"
)){
  print(type)
  split.nor.df.filter.list[[type]] <- MyGeneExp(as.data.frame(merge.bed)[,paste(prefix.list[[type]],'nor',sep = '.')],cutoff.list[[type]],1)
  colnames(split.nor.df.filter.list[[type]]) <- prefix.list[[type]]
  G14.5CD_HFD.limma.list[[type]] <- Mylimma(split.nor.df.filter.list[[type]],
                                            prefix.list[[type]][3:4],
                                            prefix.list[[type]][5:6]
  )
  G14.5CD_HFD.limma.list[[type]]$cluster <- '-'
  G14.5CD_HFD.limma.list[[type]][G14.5CD_HFD.limma.list[[type]]$logFC>0,'cluster'] <- 'G14.5HFD'
  G14.5CD_HFD.limma.list[[type]][G14.5CD_HFD.limma.list[[type]]$logFC<0,'cluster'] <- 'G14.5CD'
  table(G14.5CD_HFD.limma.list[[type]]$cluster)
  G14.5CD_HFD.limma.list[[type]]$Symbol <- mcols(merge.bed[rownames(G14.5CD_HFD.limma.list[[type]])])[,'Symbol']
}

G14.5CD_HFD.limma.filter.list <- list()
fc.tmp <- 4;type <- "G14.5HFD.H3K27ac"
G14.5CD_HFD.limma.filter.list[[type]] <- MyLimmaPlot2(G14.5CD_HFD.limma.list[[type]],P.Val.cutoff  = 0.1,fc.max = 12,fc.cutoff = log2(fc.tmp),xlim = c(0,20),main = paste('fold change:',fc.tmp,sep= ''),
                                                      cex.lab = 2,cex.axis = 2,cex.main = 2,
                                                      #axes=FALSE,
                                                      xlab = c('signal intensity')#,
                                                      #ylim = c(-6,6)
)

G14.5CD_HFD.limma.filter.list[[type]] <- G14.5CD_HFD.limma.list[[type]][G14.5CD_HFD.limma.list[[type]]$P.Value<=0.5&abs(G14.5CD_HFD.limma.list[[type]]$logFC)>=log2(fc.tmp),]
table(G14.5CD_HFD.limma.filter.list[[type]]$cluster)
G14.5HFD.up <- rownames(G14.5CD_HFD.limma.filter.list[[type]])[G14.5CD_HFD.limma.filter.list[[type]]$state=='up']
G14.5HFD.down <- rownames(G14.5CD_HFD.limma.filter.list[[type]])[G14.5CD_HFD.limma.filter.list[[type]]$cluster=='G14.5CD']

G14.5CD_HFD.limma.filter.list[[type]]$log2FC <- log2(rowMeans(split.nor.df.filter.list[[type]][rownames(G14.5CD_HFD.limma.filter.list[[type]]),prefix.list[[type]][5:6]])+0.01)-log2(rowMeans(split.nor.df.filter.list[[type]][rownames(G14.5CD_HFD.limma.filter.list[[type]]),prefix.list[[type]][3:4]])+0.01)
G14.5CD_HFD.limma.filter.list[[type]]$logFC <- G14.5CD_HFD.limma.filter.list[[type]]$log2FC
G14.5CD_HFD.limma.filter.list[[type]]$AveExpr <- log2(G14.5CD_HFD.limma.filter.list[[type]]$AveExpr+0.01)


G14.5CD_HFD.limma.filter.list[[type]]
fc.tmp <- 1
pdf('G:/lab/Article/heshuang/BYLW/chromatin/20230530.revised.Fig6d.CDHFDlimma.H3K27ac.cutoff.fc4.pva0.1.20221210.pdf')
par(oma = c(0,1,0,0))
MyLimmaPlot3(G14.5CD_HFD.limma.filter.list[[type]] ,P.Val.cutoff = 0.1,fc.max = 10,fc.cutoff = log2(1),xlim = c(-1,5),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             axes=T,
             xlab = c('log2(signal intensity)'),
             count = T
             #ylim = c(-6,6)
)

MyLimmaPlot3(G14.5CD_HFD.limma.filter.list[[type]] ,P.Val.cutoff = 0.1,fc.max = 10,fc.cutoff = log2(1),xlim = c(-1,5),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c("gray80","#79706A","#b89fc1"),
             axes=FALSE,
             xlab = c('log2(signal intensity)'),
             count = F
             #ylim = c(-6,6)
)
box(lwd=5)
axis(2,at=seq(-10,10,by=5),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(-1,5,by=1),cex.axis = 1.5,lwd = 3)

MyLimmaPlot3(G14.5CD_HFD.limma.filter.list[[type]] ,P.Val.cutoff = 0.1,fc.max = 10,fc.cutoff = log2(1),xlim = c(0,5),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c("gray80","#79706A","#b89fc1"),
             axes=FALSE,
             xlab = c('log2(signal intensity)'),
             count = F
             #ylim = c(-6,6)
)
box(lwd=5)
axis(2,at=seq(-10,10,by=5),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(0,5,by=1),cex.axis = 1.5,lwd = 3)


MyLimmaPlot3(G14.5CD_HFD.limma.filter.list[[type]] ,P.Val.cutoff = 0.1,fc.max = 6,fc.cutoff = log2(1),xlim = c(0,5),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c("gray80","#79706A","#b89fc1"),
             axes=FALSE,
             xlab = c('log2(signal intensity)'),
             count = F
             #ylim = c(-6,6)
)
box(lwd=5)
axis(2,at=seq(-6,6,by=2),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(0,5,by=1),cex.axis = 1.5,lwd = 3)

dev.off()
########
G14.5HFD.up.bed <- merge.bed[G14.5HFD.up]
dir.create('Figure')
pdf('G:/lab/Article/heshuang/BYLW/chromatin/Fig6h.CDHFD.up.box.pval0.1.fc4.pdf',3,5)

  boxplot(as.data.frame(G14.5HFD.up.bed)[,c('G14.5HFD.nor.mean','G14.5HMKOHFD.nor.mean')],
         outline=F,col=add.alpha(c("#8ECFC9","#FA7F6F"),1),ylab='signal intensity',lwd=3,notch=T,
          ylim = c(0,26),xaxt="n"
  )
  axis(1,at = 1:2,labels = c('WTHFD','HMKOHFD'),cex.axis = 0.5)
  box(lwd=3)
  P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed)[,'G14.5HFD.nor.mean'],
                        as.data.frame(G14.5HFD.up.bed)[,'G14.5HMKOHFD.nor.mean'],
                        paired = T)
  
  MyText(paste("G14.5WT KOp-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)
dev.off()

pdf('Figure/CDHFD.up.box.pval0.1.fc4.pdf',3,5)

boxplot(as.data.frame(merge.bed)[,c('G0HFD.nor.mean','G14.5CD.nor.mean','G14.5HFD.nor.mean')],
        outline=F,col=add.alpha(time.colors[c(13,14)],0.9),ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,20),xaxt="n"
)
axis(1,at = 1:3,labels = c('G0HFD','WTHFD','HMKOHFD'),cex.axis = 0.5)
box(lwd=3)
P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed)[,'G14.5HFD.nor.mean'],
                      as.data.frame(G14.5HFD.up.bed)[,'G14.5HMKOHFD.nor.mean'],
                      paired = T)

MyText(paste("G14.5WT KOp-value: \n",
             P300.p$p.value,"\n"),
       text.cex = 0.5)
dev.off()

###########
rm(list.Acss2HFD.H3K27ac.bdg)
saveRDS(list.close,'list.close.rds')
save.image('H3K27ac.RData')

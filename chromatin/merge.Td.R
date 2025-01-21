library(genomation)
library(GenomicRanges)
source('G:/pcatest/MyFunction.R')

list.p300.bdg <- readRDS("I:/serve_backup/preg_chromatin/bdg/P300/CST/list.P300.bdg.rds")
list.H3K27ac.bdg <- readRDS("I:/serve_backup/preg_chromatin/bdg/H3K27ac/list.H3K27ac.bdg.rds")
list.H3K27ac.bdg$G14.5_H3K27ac_190930 <- NULL
rm(list.p300.bdg)

list.ATAC.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/ATAC/keep/list.ATAC.bdg.rds')
list.ATAC.bdg <- list.ATAC.bdg[c("Ins1conG14.5_20190416", "Ins1Ctrl_190807_merge","Ins1G14.5_190930","Ins1G14.5_20190604","Ins1P7NB_190808","Ins1P7NB_190930")]


list.STAT3.bdg <- readRDS("I:/serve_backup/preg_chromatin/bdg/STAT3/CST/list.STAT3.bdg.rds")
list.STAT3.bdg <- list.STAT3.bdg[c(2:4,8:10)]
STAT3.com.nor.factor <- readRDS('G:/project/pregnant_mouse/chromatin/STAT3/20220518/nbl.score.bed/nor/STAT3.com.nor.factor.rds')


list.cultureH3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/culture/list.bdg.rds')
list.Stat3HMKOH3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/Stat3HMKO/CD/list.bdg.Stat3HMKO.H3K27ac.rds')
list.Acss2HMKOH3K27ac.bdg <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/Acss2KO/new/list.bdg.Acss2HMKO.H3K27ac.rds')

list.Stat3HMKOH3K27ac.bdg$Ins1G18.5_Acss2KOWT_N137_H3K27ac_221013 <- list.Acss2HMKOH3K27ac.bdg$Ins1G18.5_Acss2KOWT_N137_H3K27ac_221013
list.Stat3HMKOH3K27ac.bdg$Ins1G18.5_Acss2KOWT_N743_H3K27ac_221013 <- list.Acss2HMKOH3K27ac.bdg$Ins1G18.5_Acss2KOWT_N743_H3K27ac_221013


list.Acss2HMKOH3K27ac.bdg <- list.Acss2HMKOH3K27ac.bdg[c(1,2,3,5,6)]
Acss2.bdg.old.list <- readRDS('I:/serve_backup/preg_chromatin/bdg/H3K27ac/Acss2KO/old/list.bdg.Acss2HMKO.H3K27ac.rds')
list.Acss2HMKOH3K27ac.bdg$Ins1G18.5_Acss2HMKO_CD.36_H3K27ac_220718_1 <- Acss2.bdg.old.list$Ins1G18.5_Acss2HMKO_CD.36_H3K27ac_220718_1


Ctrl.p300.Acss2.bed <-  GRanges(seqnames = Rle('chr2'),ranges = IRanges(155516000:155564800,width = 50),strand = '*')#Acss2
Ctrl.p300.Oxtr.bed <-  GRanges(seqnames = Rle('chr6'),ranges = IRanges(112471500:112550000,width = 50),strand = '*')#oxtr
Ctrl.p300.Gbp8.bed <-  GRanges(seqnames = Rle('chr5'),ranges = IRanges(105000000:105080000,width = 50),strand = '*')#Gbp8
Ctrl.p300.Ovol2.bed <-  GRanges(seqnames = Rle('chr2'),ranges = IRanges(144300000:144340000,width = 50),strand = '*')#Gbp8

Ctrl.p300.Actb.bed <-  GRanges(seqnames = Rle('chr5'),ranges = IRanges(142900000:142910000,width = 50),strand = '*')#Acss2


type <- 'p300'
type <- 'H3K27ac'
type <- 'STAT3'
type <- 'ATAC'
type <- 'cultureH3K27ac'
type <- 'Stat3HMKOH3K27ac'
type <- 'Acss2HMKOH3K27ac'


gene <- 'Gbp8'
gene <- 'Ovol2'
gene <- 'Oxtr'
gene <- 'Acss2'
gene <- 'Actb'

G14.5.p300.bed2 <- Ctrl.p300.Oxtr.bed
G14.5.p300.bed2 <- Ctrl.p300.Gbp8.bed
G14.5.p300.bed2 <- Ctrl.p300.Actb.bed
G14.5.p300.bed2 <- Ctrl.p300.Ovol2.bed
G14.5.p300.bed2 <- Ctrl.p300.Acss2.bed


for(i in names(get(paste('list',type,'bdg',sep = '.')))[c(1:6)]){
  print(i)
  mcols(G14.5.p300.bed2)[,i] <- round(MyCalTd(G14.5.p300.bed2,get(paste('list',type,'bdg',sep = '.'))[[i]],"td"),3)
  mcols(G14.5.p300.bed2)[,i][mcols(G14.5.p300.bed2)[,i]<0] <- 0
}


G14.5.p300.bed2$G0.td <- round((G14.5.p300.bed2$Ins1Ctrl_P300_191107_1+G14.5.p300.bed2$Ins1Ctrl_P300_191107_2)/2,2)
G14.5.p300.bed2$G14.5.td <- round((G14.5.p300.bed2$Ins1G14.5_P300_190918+G14.5.p300.bed2$Ins1G14.5_P300_190930)/2,2)
G14.5.p300.bed2$P7NB.td <- round((G14.5.p300.bed2$Ins1P7NB_P300_191211_1+G14.5.p300.bed2$Ins1P7NB_P300_191211_2)/2,2)


G14.5.p300.bed2$G0.td <- round((G14.5.p300.bed2$Ctrl_H3K27ac_190528+G14.5.p300.bed2$Ctrl_H3K27ac_190621)/2,2)
G14.5.p300.bed2$G14.5.td <- round((G14.5.p300.bed2$G14.5_H3K27ac_190605+G14.5.p300.bed2$G14.5_H3K27ac_190621)/2,2)
G14.5.p300.bed2$P7NB.td <- round((G14.5.p300.bed2$P7NB_H3K27ac_190902+G14.5.p300.bed2$P7NB_H3K27ac_190930)/2,2)


G14.5.p300.bed2$G0.td <- round((G14.5.p300.bed2$Ins1conG14.5_20190416+G14.5.p300.bed2$Ins1Ctrl_190807_merge)/2,2)
G14.5.p300.bed2$G14.5.td <- round((G14.5.p300.bed2$Ins1G14.5_190930+G14.5.p300.bed2$Ins1G14.5_20190604)/2,2)
G14.5.p300.bed2$P7NB.td <- round((G14.5.p300.bed2$Ins1P7NB_190808+G14.5.p300.bed2$Ins1P7NB_190930)/2,2)


G14.5.p300.bed2$G0.td <- round((G14.5.p300.bed2$culture_Ctrl_200312+G14.5.p300.bed2$culture_Ctrl_200730)/2,2)
G14.5.p300.bed2$G14.5.td <- round((G14.5.p300.bed2$culture_Hpreg_200312+G14.5.p300.bed2$culture_Hpreg_200730)/2,2)

G14.5.p300.bed2$G0.td <- round((G14.5.p300.bed2$CTrl_STAT3_191107.2/STAT3.com.nor.factor['CTrl_STAT3_191107.2']+G14.5.p300.bed2$CTrl_STAT3_191211.1/STAT3.com.nor.factor['CTrl_STAT3_191211.1']+G14.5.p300.bed2$CTrl_STAT3_191211.2/STAT3.com.nor.factor['CTrl_STAT3_191211.2'])/3,2)
G14.5.p300.bed2$G14.5.td <- round((G14.5.p300.bed2$G14.5_STAT3_191016/STAT3.com.nor.factor['G14.5_STAT3_191016']+G14.5.p300.bed2$G14.5_STAT3_191211.1/STAT3.com.nor.factor['G14.5_STAT3_191211.1']+G14.5.p300.bed2$G14.5_STAT3_191211.2/STAT3.com.nor.factor['G14.5_STAT3_191211.2'])/3,2)

G14.5.p300.bed2$G18.5Stat3WT.td <- round((G14.5.p300.bed2$Ins1G18.5Stat3WT_M151_H3K27ac_221114+G14.5.p300.bed2$Ins1G18.5_Acss2KOWT_N137_H3K27ac_221013+G14.5.p300.bed2$Ins1G18.5_Acss2KOWT_N743_H3K27ac_221013)/3,2)
G14.5.p300.bed2$G18.5Stat3HMKO.td <- round((G14.5.p300.bed2$Ins1G18.5_Stat3HMKO_N809_H3K27ac_221013+G14.5.p300.bed2$Ins1G18.5Stat3HMKO_M145_H3K27ac_221114+G14.5.p300.bed2$Ins1G18.5Stat3HMKO_M146_H3K27ac_221114)/3,2)


G14.5.p300.bed2$G18.5Acss2WT.td <- round((G14.5.p300.bed2$Ins1G18.5_Acss2KOWT_CD_37_H3K27ac_220718_2+G14.5.p300.bed2$Ins1G18.5_Acss2KOWT_N735_H3K27ac_221013+G14.5.p300.bed2$Ins1G18.5_Acss2KOWT_N739_H3K27ac_221013)/3,2)
G14.5.p300.bed2$G18.5Acss2HMKO.td <- round((G14.5.p300.bed2$Ins1G18.5Acss2HMKO_CD_15_H3K27ac_220802+G14.5.p300.bed2$Ins1G18.5Acss2HMKO_CD_6_H3K27ac_220802+G14.5.p300.bed2$Ins1G18.5_Acss2HMKO_CD.36_H3K27ac_220718_1)/3,2)


for(i in c(#'G0','G14.5'#,'P7NB'
  'G18.5Acss2WT','G18.5Acss2HMKO'
           )){
  MyWriteBed(G14.5.p300.bed2,
           output.col = paste(i,".td",sep = ''),
           sample.name = paste(i,'.merge.',type,'.',gene,sep = ''),
           is.gz = T,
           paste("I:/serve_backup/preg_chromatin/bdg/H3K27ac/Acss2KO/new/",i,'.merge.',type,'.',gene,".bdg.gz",sep = ""))

}

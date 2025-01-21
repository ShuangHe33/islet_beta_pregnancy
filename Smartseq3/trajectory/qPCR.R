###########qPCR#########
source("G:/pcatest/MyFunction.R")
setwd("G:/weekly report2/pregnant/Figure/qPCR/")
load("qPCR.ovol2.RData")

qPCR.tab <- read.csv("Ovol2.sc.qPCR.csv")
qPCR.tab <- qPCR.tab[,c(1:3)]
qPCR.tab$Treatment <- factor(qPCR.tab$Treatment,
                             levels = c("P21-ctrl",
                                        "G14.5"))


Ovol2.p <- wilcox.test(qPCR.tab$Ovol2_a[qPCR.tab$Treatment=='G14.5'],qPCR.tab$Ovol2_a[qPCR.tab$Treatment=='P21-ctrl'],
            paired = F)
Ovol2.p$p.value#0.0003498919


Chgb.p <- wilcox.test(Ttr.qPCR.tab$ChgB[Ttr.qPCR.tab$Treatment=='G14.5'],Ttr.qPCR.tab$ChgB[Ttr.qPCR.tab$Treatment=='P7-con'],
                       paired = F)
Chgb.p$p.value#6.981662e-05

Ttr.p <- wilcox.test(Ttr.qPCR.tab$Ttr[Ttr.qPCR.tab$Treatment=='G14.5'],Ttr.qPCR.tab$Ttr[Ttr.qPCR.tab$Treatment=='P7-con'],
                      paired = F)
Ttr.p$p.value#5.293801e-06


gene.list <- c("Ovol2_a",
               "Ovol2_b")

Pdx1.tab <- read.csv('Pdx1Nkx6.1.sc.qPCR.csv',stringsAsFactors = F)
colnames(Pdx1.tab) <- c('Treatment','Pdx1','Nkx6.1','Mlxipl')
Pdx1.tab$Treatment <- factor(Pdx1.tab$Treatment,levels = c('G0','G18.5Acss2WT','G18.5Acss2HMKO'))

pdf("Pdx1.Nkx6.1.Mlxipl.q.pdf",
    12,15)
par(mfrow = c(3,4))
for(j in c('Pdx1','Nkx6.1','Mlxipl')){
  tmp.tpm <- as.numeric(Pdx1.tab[,j])
  names(tmp.tpm) <- Pdx1.tab$Treatment
  MyViolinBeeSwarmMed(Pdx1.tab$Treatment,
                      tmp.tpm,
                      color.violin = add.alpha(c("#b6d4a8","#752a78","#c7a96c"),alpha = 0.7),
                      ylab.plot  = "relative expression",
                      box.lwd = 3,
                      j
                      
                      
  )
}
dev.off()

setwd('G:/lab/Article/heshuang/BYLW/sm3/ref/')

pdf("Ovol2.q.pdf",
    15,15)
par(mfrow = c(3,4))
for(j in gene.list){
  tmp.tpm <-as.numeric(qPCR.tab[,j])
  names(tmp.tpm) <- qPCR.tab$Treatment
  
  MyViolinBeeSwarmMed(qPCR.tab$Treatment,
                      tmp.tpm,
                      color.violin = add.alpha(c("#b6d4a8","#b89fc1"),alpha = 0.7),
                      ylab.plot  = "relative expression",
                      box.lwd = 3,
                      j
                      
                      
  )
}
dev.off()

Ttr.qPCR.tab <- read.csv("G:/weekly report2/pregnant/Figure/qPCR/chgB.Ttr..sc.qPCR.csv")
Ttr.qPCR.tab <- Ttr.qPCR.tab[,c(1:3)]
Ttr.qPCR.tab$Treatment <- factor(Ttr.qPCR.tab$Treatment,
                                 levels = c("P7-con",
                                            "G14.5"))
gene.list <- c("ChgB",
               "Ttr")
pdf("Chgb.Ttr.marker.q.pdf",
    15,15)
par(mfrow = c(3,4))
for(j in gene.list){
  tmp.tpm <-as.numeric(Ttr.qPCR.tab[,j])
  names(tmp.tpm) <- Ttr.qPCR.tab$Treatment
  
  MyViolinBeeSwarmMed(Ttr.qPCR.tab$Treatment,
                      tmp.tpm,
                      color.violin = add.alpha(c("#b6d4a8","#b89fc1"),alpha = 0.7),
                      ylab.plot  = "relative expression",
                      box.lwd = 3,
                      j
                      
                      
  )
}
dev.off()
seu.ref.KO.list2 <- readRDS('G14.5G18.5KO/seu.ref.KO.list2.rds')
seu.ref.KO.list <- readRDS('seu.ref.KO.list2.rds')
pdf('G:/lab/Article/heshuang/BYLW/sm3/KO/QC.genecount.pdf',4,5)

boxplot(seu.ref.KO.list$G0G14.5P300WTHMKO$Genecount_all,lwd=2,main='p300 gene number')
box(lwd=3)
boxplot(seu.ref.KO.list$G0G14.5P300WTHMKO$Raw_count,lwd=2,main='p300 read count')
box(lwd=3)

boxplot(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Genecount_all,lwd=2,main='Acss2 gene number')
box(lwd=3)
boxplot(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Raw_count,lwd=2,main='Acss2 read count')
box(lwd=3)

boxplot(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Genecount_all,lwd=2,main='Stat3 gene number')
box(lwd=3)
boxplot(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Raw_count,lwd=2,main='Stat3 read count')
box(lwd=3)

dev.off()
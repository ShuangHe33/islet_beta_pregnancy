fc.tmp <- 1
G0_G14.limma.ATAC.filter <- MyLimmaPlot2(G0_G14.limma.list$ATAC,P.Val.cutoff = 0.2,fc.max = 10,fc.cutoff = log2(1),xlim = c(0,12),main = paste('fold change:',fc.tmp,sep= ''),
                                            cex.lab = 2,cex.axis = 2,cex.main = 2,
                                            #axes=FALSE,
                                            xlab = c('signal intensity')#,
                                            #ylim = c(-6,6)
)
G0_G14.limma.ATAC.filter$state <- as.character(G0_G14.limma.ATAC.filter$state)


ATAC.up <- rownames(G0_G14.limma.ATAC.filter)[G0_G14.limma.ATAC.filter$state=='up']
ATAC.down <- rownames(G0_G14.limma.ATAC.filter)[G0_G14.limma.ATAC.filter$state=='down']

#ATAC.up.final <- ATAC.up[merge.bed[ATAC.up]$==1]
ATAC.up.bed <- merge.bed[ATAC.up]
length(ATAC.up.bed)
min(ATAC.up.bed$G14.5.ATAC.mean)
ATAC.up.bed <- ATAC.up.bed[ATAC.up.bed$G14.5.ATAC.mean>=cutoff.list$ATAC]
table(ATAC.up.bed$G14.5.ATAC.state)
length(ATAC.up.bed)#4777

G0_G14.limma.ATAC.filter[c(ATAC.up[merge.bed[ATAC.up]$G14.5.ATAC.mean<cutoff.list$ATAC],
                           ATAC.down[merge.bed[ATAC.down]$Ctrl.ATAC.mean<cutoff.list$ATAC]),'state'] <- '-'


G0_G14.limma.ATAC.filter$log2FC <- log2(rowMeans(split.nor.df.filter.list[[type]][rownames(G0_G14.limma.ATAC.filter),prefix.list[[type]][3:4]])+0.01)-log2(rowMeans(split.nor.df.filter.list[[type]][rownames(G0_G14.limma.ATAC.filter),prefix.list[[type]][1:2]])+0.01)
G0_G14.limma.ATAC.filter$logFC <- G0_G14.limma.ATAC.filter$log2FC
G0_G14.limma.ATAC.filter$AveExpr <- log2(G0_G14.limma.ATAC.filter$AveExpr+0.01)


MyWriteBed(ATAC.up.bed,'ATAC.up.4777.bed')

pdf('G:/lab/Article/heshuang/BYLW/chromatin/20230530.revised.raw.F2c.limma.ATAC.cutoff.pval0.2.2.pdf')
par(oma = c(0,1,0,0))
MyLimmaPlot3(G0_G14.limma.ATAC.filter,P.Val.cutoff = 0.2,fc.max = 5,fc.cutoff = log2(fc.tmp),xlim = c(-2,8),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c('gray60','gold2','deeppink3'),
             axes=T,
             xlab = c('log2(signal intensity)'),
             count = T
             #ylim = c(-6,6)
)


MyLimmaPlot3(G0_G14.limma.ATAC.filter,P.Val.cutoff = 0.2,fc.max = 6,fc.cutoff = log2(fc.tmp),xlim = c(-2,8),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c('gray60','gold2','deeppink3'),
             axes=FALSE,
             xlab = c('log2(signal intensity)'),
            # axes=F,
             count = F
             #ylim = c(-6,6)
)

axis(2,at=seq(-6,6,by=2),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(-2,8,by=2),cex.axis = 1.5,lwd = 3)
box(lwd =3)
dev.off()


cluster2.venn <- Venn(list('cluster1.peak' = cluster2.bed$rowname,
                           'ATAC.up' = ATAC.up.bed$rowname))
plot(cluster2.venn)


cluster3.venn <- Venn(list('cluster2.peak' = cluster3.bed$rowname,
                           'ATAC.up' = ATAC.up.bed$rowname))
plot(cluster3.venn)


up.gene.venn <- Venn(list('cluster1.gene' = preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up'],
                          'ATAC.up.gene' = unique(ATAC.up.bed$gene.id)
))
plot(up.gene.venn)



up.gene.venn <- Venn(list('1510H3K27ac.up.gene' = unique(H3K27ac.up.bed$gene.id),
                          'ATAC.up.gene' = unique(ATAC.up.bed$gene.id)
))
plot(up.gene.venn)


up.peak.venn <- Venn(list('1510H3K27ac.up.peak' = H3K27ac.up.bed$rowname,
                          'ATAC.up.peak' = ATAC.up.bed$rowname
))
plot(up.peak.venn)

cluster2.venn <- Venn(list('cluster1.peak' = cluster2.bed$rowname,
                           '1510.H3K27ac.up' = H3K27ac.up.bed$rowname))
plot(cluster2.venn)


cluster3.venn <- Venn(list('cluster2.peak' = cluster3.bed$rowname,
                           '709.H3K27ac.down' = H3K27ac.down.bed$rowname))
plot(cluster3.venn)




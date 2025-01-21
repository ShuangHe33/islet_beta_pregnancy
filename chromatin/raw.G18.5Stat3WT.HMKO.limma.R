MyLimmaPlot3 <- 
  function(de.data,
           P.Val.cutoff = 0.05,
           fc.max1,
           fc.max2,
           fc.cutoff = 1,
           xlim,
           de.color=c("black","#F25362","#487CC2"),
           xlab = "mean of normalized counts",
           ylab = "log2 fold change",
           count = T,
           box.lwd=3,
           text.cex=1,
           # up.state,
           # down.state,
           ...){
    de.data <- as.data.frame(de.data)
    # de.data$state <- '-'
    # de.data[(!is.na(de.data$P.Value)) & 
    #           (de.data$P.Value <= P.Val.cutoff) & 
    #           (de.data$logFC <= -fc.cutoff),'state'] <- 'down'
    # de.data[(!is.na(de.data$P.Value)) & 
    #           (de.data$P.Value <= P.Val.cutoff) & 
    #           (de.data$logFC >= fc.cutoff),'state'] <- 'up'
    de.data$state <- factor(de.data$state,levels = c('-','up','down'))
    de.data <- de.data[order(de.data$state),]
    up.count <- sum(de.data$state=='up')
    down.count <- sum(de.data$state=='down')
    # shape.plot <- rep(16,length(de.data$logFC))
    # shape.plot[abs(de.data$logFC) >= fc.max] <- 17
    de.data$logFC[de.data$logFC >= fc.max2] <- fc.max2
    de.data$logFC[de.data$logFC <= -fc.max1] <- -fc.max1
    
    plot(de.data$AveExpr,
         de.data$logFC,
         xlab = xlab,
         ylab = ylab,
         ylim = c(-fc.max1, fc.max2),
         xlim = xlim,
         col = de.color[de.data$state],
         #log = "x",
         pch = 16,
         cex = 1,
         ...)
    box(lwd=box.lwd)
    if(count){
      text(1,
           fc.max2 / 2,
           as.character(up.count),
           #col = up.color,
           cex = text.cex)
      text(1,
           -(fc.max1 / 2),
           as.character(down.count),
           #col = down.color,
           cex = text.cex)
    }
  }


merge.WT.bed <- merge.bed3[merge.bed3$G18.5Stat3WTHMKO.state!='01']

split.nor.df.filter.list <- list()
WT_HMKO.limma.list <- list()
for(type in c(
  'H3K27ac')){
  print(type)
  split.nor.df.filter.list[[type]] <- MyGeneExp(as.data.frame(merge.WT.bed)[,prefix.list[[type]]],cutoff.list$WT,2)
  colnames(split.nor.df.filter.list[[type]]) <- prefix.list[[type]]
  WT_HMKO.limma.list[[type]] <- Mylimma(split.nor.df.filter.list[[type]],
                                       prefix.list[[type]][4:6],
                                       prefix.list[[type]][1:3]
  )
  WT_HMKO.limma.list[[type]]$cluster <- '-'
  WT_HMKO.limma.list[[type]][WT_HMKO.limma.list[[type]]$logFC>0,'cluster'] <- 'WT'
  WT_HMKO.limma.list[[type]][WT_HMKO.limma.list[[type]]$logFC<0,'cluster'] <- 'HMKO'
  table(WT_HMKO.limma.list[[type]]$cluster)
  WT_HMKO.limma.list[[type]]$Symbol <- mcols(merge.bed3[rownames(WT_HMKO.limma.list[[type]])])[,'Symbol']
}
##############
fc.tmp <- 1;type <- 'H3K27ac'
WT_HMKO.limma.H3K27ac.filter <- MyLimmaPlot2(WT_HMKO.limma.list[[type]],P.Val.cutoff  = 0.2,fc.max = 10,fc.cutoff = log2(fc.tmp),xlim = c(0,15),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             #axes=FALSE,
             xlab = c('signal intensity')#,
             #ylim = c(-6,6)
)
WT_HMKO.limma.H3K27ac.filter$state <- as.character(WT_HMKO.limma.H3K27ac.filter$state)


H3K27ac.up <- rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='up']
H3K27ac.up.final <- H3K27ac.up[merge.bed[H3K27ac.up]$WT.State==1]
H3K27ac.down <- rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='down']
H3K27ac.down[merge.bed[H3K27ac.down]$HMKO.State==1]
WT_HMKO.limma.H3K27ac.filter[c(H3K27ac.up[merge.bed[H3K27ac.up]$HMKO.State==0],
                               H3K27ac.down[merge.bed[H3K27ac.down]$WT.State==0]),'state'] <- '-'


WT_HMKO.limma.H3K27ac.filter$log2FC <- log2(rowMeans(split.nor.df.filter.list[[type]][rownames(WT_HMKO.limma.H3K27ac.filter),prefix.list[[type]][1:3]])+0.01)-log2(rowMeans(split.nor.df.filter.list[[type]][rownames(WT_HMKO.limma.H3K27ac.filter),prefix.list[[type]][4:6]])+0.01)
WT_HMKO.limma.H3K27ac.filter$logFC <- WT_HMKO.limma.H3K27ac.filter$log2FC
WT_HMKO.limma.H3K27ac.filter$AveExpr <- log2(WT_HMKO.limma.H3K27ac.filter$AveExpr+0.01)

fc.tmp <- 1
pdf('G:/lab/Article/heshuang/BYLW/chromatin/20230530.revised.Fig4h.raw.limma.H3K27ac.cutoff.pval0.2.2.pdf')
par(oma = c(0,1,0,0))
MyLimmaPlot3(WT_HMKO.limma.H3K27ac.filter,P.Val.cutoff = 0.2,fc.max1 = 4,fc.max2 = 10,fc.cutoff = log2(1),xlim = c(0,5),main = paste('fold change:',fc.tmp,sep= ''),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
            # axes=FALSE,
             xlab = c('log2(signal intensity)'),
             count = T
             #ylim = c(-6,6)
)
MyLimmaPlot3(WT_HMKO.limma.H3K27ac.filter,P.Val.cutoff = 0.2,fc.max1 = 4,fc.max2 = 10,fc.cutoff = log2(fc.tmp),xlim = c(0,5),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             axes=FALSE,de.color=c('black','cornflowerblue','violetred'),
             xlab = c('log2(signal intensity)'),
             count = F
             #ylim = c(-6,6)
)

axis(2,at=seq(-4,10,by=2),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(0,5,by=1),cex.axis = 1.5,lwd = 3)
box(lwd = 3)
dev.off()

MylimmafcDataplot(merge.bed3,rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='down'],
                  rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='up'],
                  'G18.5.WT.mean','G18.5.HMKO.mean','WT.HMKO.state',cutoff.list$WT,xlim=c(-5,5),ylim=c(-5,5)
                  )

########
preg.gene.tab <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/LMgene/gene_inf/cg.preg.gene.si.1e28.addcc.tab')
rownames(preg.gene.tab) <- preg.gene.tab$EnsemblGeneID

merge.bed3$G18.5.WT.mean <- round((merge.bed3$Ins1G18.5_Acss2KOWT_N137_H3K27ac_221013+merge.bed3$Ins1G18.5_Acss2KOWT_N743_H3K27ac_221013+merge.bed3$Ins1G18.5Stat3WT_M151_H3K27ac_221114)/3,3)
merge.bed3$G18.5.HMKO.mean <- round((merge.bed3$Ins1G18.5_Stat3HMKO_N809_H3K27ac_221013+merge.bed3$Ins1G18.5Stat3HMKO_M145_H3K27ac_221114+merge.bed3$Ins1G18.5Stat3HMKO_M146_H3K27ac_221114)/3,3)
merge.bed$G18.5.WT.mean <- round((merge.bed$Ins1G18.5_Acss2KOWT_N137_H3K27ac_221013+merge.bed$Ins1G18.5_Acss2KOWT_N743_H3K27ac_221013+merge.bed$Ins1G18.5Stat3WT_M151_H3K27ac_221114)/3,3)
merge.bed$G18.5.HMKO.mean <- round((merge.bed$Ins1G18.5_Stat3HMKO_N809_H3K27ac_221013+merge.bed$Ins1G18.5Stat3HMKO_M145_H3K27ac_221114+merge.bed$Ins1G18.5Stat3HMKO_M146_H3K27ac_221114)/3,3)

cluster2.bed <- merge.bed3[merge.bed3$gene.id %in% preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up']]
length(cluster2.bed)
cluster3.bed <- merge.bed3[merge.bed3$gene.id %in% preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.down']]

rm(list.Acss2.bdg)
list.cluster23.bed <- list()
list.cluster23.bed[['cluster2']] <- cluster2.bed
list.cluster23.bed[['cluster3']] <- cluster3.bed
for(tmp.cluster in names(list.cluster23.bed)){
  print(tmp.cluster)
  list.cluster23.bed[[paste(tmp.cluster,'.p',sep = '')]] <- list.cluster23.bed[[tmp.cluster]][list.cluster23.bed[[tmp.cluster]]$promoter2==1]
  list.cluster23.bed[[paste(tmp.cluster,'.e',sep = '')]] <- list.cluster23.bed[[tmp.cluster]][list.cluster23.bed[[tmp.cluster]]$distalregion==1]
}
sapply(list.cluster23.bed, length)
############
pdf('G:/lab/Article/heshuang/BYLW/chromatin/Figs7n.raw.cluster1.2.gene.signal.intensity.pdf',3,5)

MyLogBoxPlot(as.data.frame(cluster2.bed)[,c('G18.5.WT.mean','G18.5.HMKO.mean')],log = F,
             min.cutoff = 0,notch=T,col=c('cornflowerblue','violetred'),ylim = c(0,16),lwd=4,ylab='signal intensity', main = 'cluster1'
)
box(lwd=3)
WT.p <- wilcox.test(as.data.frame(cluster2.bed)[,c('G18.5.WT.mean')],
                    as.data.frame(cluster2.bed)[,c('G18.5.HMKO.mean')],
                    #  alternative = 'less',
                    pair=T)
MyText(paste("cluster1 WT:p-value: \n",
             WT.p$p.value,"\n"),
       text.cex = 0.6)

MyLogBoxPlot(as.data.frame(cluster3.bed)[,c('G18.5.WT.mean','G18.5.HMKO.mean')],log = F,
             min.cutoff = 0,notch=T,col=c('cornflowerblue','violetred'),ylim = c(0,16),lwd=4,ylab='signal intensity', main = 'cluster2'
)
box(lwd=3)
WT.p <- wilcox.test(as.data.frame(cluster3.bed)[,c('G18.5.WT.mean')],
                    as.data.frame(cluster3.bed)[,c('G18.5.HMKO.mean')],
                    pair=T#,alternative = 'greater'
)
MyText(paste("cluster2 WT:p-value: \n",
             WT.p$p.value,"\n"),
       text.cex = 1)

dev.off()
########
list.cluster23.bed[['cluster2.1kb']] <- list.cluster23.bed$cluster2[list.cluster23.bed$cluster2@ranges@width>=1000]
list.cluster23.bed[['cluster3.1kb']] <- list.cluster23.bed$cluster3[list.cluster23.bed$cluster3@ranges@width>=1000]


pdf('Figure/Figs7n.e.p.raw.cluster2.3.gene.signal.intensity.pdf',3,5)
for(tmp.cluster in names(list.cluster23.bed)){
  MyLogBoxPlot(as.data.frame(list.cluster23.bed[[tmp.cluster]])[,c('G18.5.WT.mean','G18.5.HMKO.mean')],log = F,
               min.cutoff = 0,notch=T,col=add.alpha(time.colors[c(7,8)],0.8),ylim = c(0,16),lwd=4,ylab='signal intensity', main = tmp.cluster
  )
  box(lwd=3)
  WT.p <- t.test(as.data.frame(list.cluster23.bed[[tmp.cluster]])[,c('G18.5.WT.mean')],
                      as.data.frame(list.cluster23.bed[[tmp.cluster]])[,c('G18.5.HMKO.mean')],
                      alternative = 'greater',
                      pair=T)
  MyText(paste("WT:p-value: \n",
               WT.p$p.value,"\n"),
         text.cex = 0.6)
}
dev.off()
######
H3K27ac.fc <- 2
color.fc <- c('gray60',time.colors[c(26,25)])
pdf('Figure/Fig4h.G18.5Stat3WTHFD.fc2.pdf',6,6.5)
merge.bed2 <- MyfcDataplot2(merge.bed3[merge.bed3$G18.5Stat3WTHMKO.state!='01'],
                            "G18.5.HMKO.mean",
                            "G18.5.WT.mean",
                            "H3K27ac.fc2.state",
                            "H3K27ac.fc2.up",
                            "H3K27ac.fc2.down",
                            H3K27ac.fc,
                            cutoff.list[1],
                            c(-6,5),
                            c(0,5)
)
dev.off()

########
save.image('G18.5Stat3HMKO.RData')

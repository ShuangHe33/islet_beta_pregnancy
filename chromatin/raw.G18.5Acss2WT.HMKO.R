
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


merge.WT.bed <- merge.bed3[merge.bed3$G18.5Acss2WTHMKO.state!='01']
prefix.list <- list()
prefix.list[[type]] <- G18.5Acss2.prefix
split.nor.df.filter.list <- list()
WT_HMKO.limma.list <- list()
for(type in c(
  'Acss2')){
  print(type)
  split.nor.df.filter.list[[type]] <- MyGeneExp(as.data.frame(merge.WT.bed)[,prefix.list[[type]]],round(cutoff.list$WT,1),2)
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
WT_HMKO.limma.list[[type]]$log2FC <- log2(rowMeans(split.nor.df.filter.list[[type]][rownames(WT_HMKO.limma.list[[type]]),prefix.list[[type]][1:3]])+0.01)-log2(rowMeans(split.nor.df.filter.list[[type]][rownames(WT_HMKO.limma.list[[type]]),prefix.list[[type]][4:6]])+0.01)
WT_HMKO.limma.list[[type]]$logFC <- WT_HMKO.limma.list[[type]]$log2FC
WT_HMKO.limma.list[[type]]$AveExpr <- log2(WT_HMKO.limma.list[[type]]$AveExpr+0.01)


fc.tmp <- 1;type <- 'Acss2'
WT_HMKO.limma.H3K27ac.filter <- MyLimmaPlot2(WT_HMKO.limma.list[[type]],P.Val.cutoff  = 0.2,fc.max = 10,fc.cutoff = log2(fc.tmp),xlim = c(0,20),main = paste('fold change:',fc.tmp,sep= ''),
                                             cex.lab = 2,cex.axis = 2,cex.main = 2,
                                             #axes=FALSE,
                                             xlab = c('signal intensity')#,
                                             #ylim = c(-6,6)
)
WT_HMKO.limma.H3K27ac.filter$state <- as.character(WT_HMKO.limma.H3K27ac.filter$state)

# H3K27ac.up <- rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='up']
# H3K27ac.up.final <- H3K27ac.up[merge.bed[H3K27ac.up]$HMKO.State==1]
# H3K27ac.down <- rownames(WT_HMKO.limma.H3K27ac.filter)[WT_HMKO.limma.H3K27ac.filter$state=='down']
# H3K27ac.down[merge.bed[H3K27ac.down]$WT.State==1]
WT_HMKO.limma.H3K27ac.filter[c(H3K27ac.up[merge.bed[H3K27ac.up]$HMKO.State==0],
                               H3K27ac.down[merge.bed[H3K27ac.down]$WT.State==0]),'state'] <- '-'

pdf('G:/lab/Article/heshuang/BYLW/chromatin/20230530.revised.F5h.raw.limma.H3K27ac.cutoff.pval0.2.2.pdf')
par(oma = c(0,1,0,0))
WT_HMKO.limma.H3K27ac.filter <- MyLimmaPlot2(WT_HMKO.limma.list[[type]],P.Val.cutoff  = 0.2,fc.max= 10,fc.cutoff = log2(fc.tmp),xlim = c(0,6),main = paste('fold change:',fc.tmp,sep= ''),
                                             cex.lab = 2,cex.axis = 2,cex.main = 2,de.color = c('gray60',"#A65628", "#5F9EA0"),
                                             #axes=FALSE,
                                             xlab = c('log2(signal intensity)')#,
                                             #ylim = c(-6,6)
)
MyLimmaPlot3(WT_HMKO.limma.H3K27ac.filter,P.Val.cutoff = 0.2,fc.max1 = 4,fc.max2 = 6,fc.cutoff = log2(fc.tmp),xlim = c(0,5),
             cex.lab = 2,cex.axis = 2,cex.main = 2,
             axes=FALSE,de.color=c('gray60',"#A65628", "#5F9EA0"),
             xlab = c('log2(signal intensity)'),
             count = F
             #ylim = c(-6,6)
)

axis(2,at=seq(-4,6,by=2),cex.axis=1.5,lwd = 3#,tck=-0.05
)
axis(1,at = seq(0,5,by=1),cex.axis = 1.5,lwd = 3)
box(lwd = 3)
# MyLimmaPlot3(WT_HMKO.limma.H3K27ac.filter,P.Val.cutoff = 0.2,fc.max1 = 4,fc.max2 = 6,fc.cutoff = log2(fc.tmp),xlim = c(0,6),
#              cex.lab = 2,cex.axis = 2,cex.main = 2,
#              axes=FALSE,de.color=c('gray60','royalblue1','magenta1'),
#              xlab = c('log2(signal intensity)'),
#              count = F
#              #ylim = c(-6,6)
# )
# 
# axis(2,at=seq(-4,6,by=2),cex.axis=1.5,lwd = 3#,tck=-0.05
# )
# axis(1,at = seq(0,6,by=1),cex.axis = 1.5,lwd = 3)
# box(lwd = 3)
dev.off()
########
KO.up.bed <- merge.bed3[merge.bed3$G18.5Acss2=='G18.5_Acss2HMKO']
KO.down.bed <- merge.bed3[merge.bed3$G18.5Acss2=='G18.5_Acss2WT']

list.KO.bed <- list()
list.KO.bed[['KO.down']] <- KO.down.bed
list.KO.bed[['KO.up']] <- KO.up.bed
sapply(list.KO.bed, length)

for(tmp.cluster in names(list.KO.bed)){
  print(tmp.cluster)
  list.KO.bed[[paste(tmp.cluster,'.p',sep = '')]] <- list.KO.bed[[tmp.cluster]][list.KO.bed[[tmp.cluster]]$promoter2==1]
  list.KO.bed[[paste(tmp.cluster,'.e',sep = '')]] <- list.KO.bed[[tmp.cluster]][list.KO.bed[[tmp.cluster]]$distalregion==1]
}
sapply(list.KO.bed, length)

list.KO.bed[[paste(tmp.cluster,'.1kb',sep = '')]] <- list.KO.bed[[tmp.cluster]][list.KO.bed[[tmp.cluster]]@ranges@width>1000]
list.KO.bed[[paste(tmp.cluster,'.1kb',sep = '')]] <- list.KO.bed[[tmp.cluster]][list.KO.bed[[tmp.cluster]]@ranges@width>1000]

dir.create('Figure')

pdf('G:/lab/Article/heshuang/BYLW/chromatin/Fig5l.raw.all.acss2HMKO.H3K27ac.signal.intensity.pdf',3,5)
MyLogBoxPlot(as.data.frame(merge.WT.bed)[,c('G18.5.WT.mean','G18.5.HMKO.mean')],log = F,
             min.cutoff = 0,notch=T,col=add.alpha("#A65628", "#5F9EA0",0.8),ylim = c(0,18),lwd=4,ylab='signal intensity', main = 'overall')

WT.p <- t.test(as.data.frame(merge.WT.bed)[,c('G18.5.WT.mean')],
               as.data.frame(merge.WT.bed)[,c('G18.5.HMKO.mean')],
               #alternative = 'greater',
               pair=T)
box(lwd=3)
MyText(paste("WT:p-value: \n",
             WT.p$p.value,"\n"),
       text.cex = 0.6)
dev.off()

pdf('G:/lab/Article/heshuang/BYLW/chromatin/Figs9m.raw.all.ko.up.down.signal.intensity.pdf',3,5)
for(tmp.cluster in names(list.KO.bed)){
  MyLogBoxPlot(as.data.frame(list.KO.bed[[tmp.cluster]])[,c('G18.5.WT.mean','G18.5.HMKO.mean')],log = F,
               min.cutoff = 0,notch=T,col=add.alpha(c("#A65628", "#5F9EA0"),0.8),ylim = c(0,18),lwd=4,ylab='signal intensity', main = tmp.cluster
  )
  box(lwd=3)
  WT.p <- t.test(as.data.frame(list.KO.bed[[tmp.cluster]])[,c('G18.5.WT.mean')],
                      as.data.frame(list.KO.bed[[tmp.cluster]])[,c('G18.5.HMKO.mean')],
                        #alternative = 'greater',
                      pair=T)
  MyText(paste("WT:p-value: \n",
               WT.p$p.value,"\n"),
         text.cex = 0.6)
}
dev.off()
###########reviewer gene######
gene.list <- c('Gck','Slc2a2','Hadh','Ucp2',
               'Pdx1','Rfx6','Neurod1','Glis3','Hnf4a',
               "Scamp5","Stxbp1",'Stxbp2',"Syt13","Cacna1a","Myo6","Rims3","Gnai2","Cplx2","Bcl2l1","Rph3al","Grik5","Braf","Baiap3",'Rab3b','Ncam1','Exoc6b'
               )
merge.function.bed <- merge.bed3[merge.bed3$Symbol %in% gene.list]
length(merge.function.bed)#56

gene.list <- list()
gene.list[['glucose_uptake']] <- c('Gck','Slc2a2','Hadh','Ucp2')
gene.list[['Vesicle_trafficking']] <- c(
                                        "Scamp5","Stxbp1",'Stxbp2',"Syt13","Cacna1a","Myo6","Rims3","Gnai2","Cplx2","Bcl2l1","Rph3al","Grik5","Braf","Baiap3",'Rab3b','Ncam1','Exoc6b'
)
gene.list[['TF']] <- c('Pdx1','Rfx6','Neurod1','Glis3','Hnf4a')

merge.function.bed.list <- list()

for(tmp.fun in names(gene.list)){
  merge.function.bed.list[[tmp.fun]] <- merge.bed3[merge.bed3$Symbol %in% gene.list[[tmp.fun]]]
}

mcols(merge.function.bed.list$glucose_uptake)[,c('Symbol','G18.5.WT.mean','G18.5.HMKO.mean')]
mcols(merge.function.bed.list$Vesicle_trafficking)[,c('Symbol','G18.5.WT.mean','G18.5.HMKO.mean')]
mcols(merge.function.bed.list$TF)[,c('Symbol','G18.5.WT.mean','G18.5.HMKO.mean')]


merge.function.df <- as.data.frame(merge.function.bed)[,c('seqnames','start','end','Symbol','G18.5.WT.mean','G18.5.HMKO.mean')]
MyWriteTable(merge.function.df,'merge.glucose.function.tab')

########
save.image('raw.G18.5Stat3HMKO.RData')

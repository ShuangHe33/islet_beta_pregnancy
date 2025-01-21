dir.create('Figure')

si.tab.list <- list()
list.nor.df2 <- list()
list.nor.mean.df <- list()

time.list <- list()
time.list[['H3K4me1']] <- c(rep(c('G0','G14.5',
                                  #'G14.5HFD',
                                  'P7NB'),each=2))
time.list[['ATAC']] <- time.list[['H3K4me1']]
time.list[['H3K27ac']] <- time.list[['H3K4me1']]
#time.list[['H3K27ac']] <- c(rep(c('G0','G14.5'),each=2),rep('G14.5HFD',3),rep('P7NB',2))
for(his in c('H3K4me1','H3K27ac','ATAC')){
  print(his)

  si.tab.list[[his]] <- data.frame('SampleName'=paste(prefix.list[[his]],'',sep = ''),
                                   'Time' = time.list[[his]]
                                   )
  si.tab.list[[his]][,'Time'] <- factor(si.tab.list[[his]][,'Time'],levels = c('G0','G14.5',
                                                                              # 'G14.5HFD',
                                                                               'P7NB'))
  rownames(si.tab.list[[his]]) <- si.tab.list[[his]]$SampleName
  
  list.nor.df2[[prefix]] <- round(as.data.frame(merge.bed)[,as.character(si.tab.list[[his]]$SampleName)],3)
  list.nor.mean.df[[his]] <- do.call(cbind,by(t(list.nor.df2[[prefix]]),si.tab.list[[his]][,'Time'],colMeans))
}

HFD.time.col <- time.colors[c(1,8,#12,
                              15)]

list.nor.mean.df[['ATAC']] <- do.call(cbind,by(t(as.data.frame(merge.bed)[,prefix.list$ATAC]),si.tab.list[[his]][,'Time'],colMeans))

merge.bed$Ctrl.H3K27ac.mean <- list.nor.mean.df$H3K27ac[merge.bed$rowname,1]
merge.bed$G14.5.H3K27ac.mean <- list.nor.mean.df$H3K27ac[merge.bed$rowname,2]


merge.bed$Ctrl.ATAC.mean <- list.nor.mean.df$ATAC[merge.bed$rowname,1]
merge.bed$G14.5.ATAC.mean <- list.nor.mean.df$ATAC[merge.bed$rowname,2]

saveRDS(list.nor.mean.df,'list.nor.mean.df.rds')
list.nor.mean.df
###########
preg.gene.tab <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/LMgene/gene_inf/cg.preg.gene.si.1e28.addcc.tab')
cluster2.bed <- merge.bed[merge.bed$gene.id %in% preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up']]
cluster3.bed <- merge.bed[merge.bed$gene.id %in% preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.down']]

list.H3K27ac.element.bed <- list()
list.H3K27ac.element.bed[['pc1.p.up']] <- H3K27ac.pc1.p.up.bed
list.H3K27ac.element.bed[['pc1.distal.up']] <- H3K27ac.pc1.distal.up.bed
list.H3K27ac.element.bed[['distal.up']] <- H3K27ac.up.bed[H3K27ac.up.bed$enhancer==1]
list.H3K27ac.element.bed[['p.up']] <- H3K27ac.up.bed[H3K27ac.up.bed$promoter2==1]
list.H3K27ac.element.bed[['distal.down']] <- H3K27ac.down.bed[H3K27ac.down.bed$enhancer==1]
list.H3K27ac.element.bed[['p.down']] <- H3K27ac.down.bed[H3K27ac.down.bed$promoter2==1]


#list.H3K27ac.element.bed[['up']] <- H3K27ac.up.bed
list.H3K27ac.element.bed[['cluster1.p']] <- cluster2.bed[cluster2.bed$promoter2==1]
list.H3K27ac.element.bed[['cluster1.e']] <- cluster2.bed[cluster2.bed$enhancer==1]
list.H3K27ac.element.bed[['cluster2.p']] <- cluster3.bed[cluster3.bed$promoter2==1]
list.H3K27ac.element.bed[['cluster2.e']] <- cluster3.bed[cluster3.bed$enhancer==1]


y.max.list <- list()
y.max.list[['H3K4me1']] <- rep(8,length(list.H3K27ac.element.bed))
names(y.max.list[['H3K4me1']]) <- names(list.H3K27ac.element.bed)
y.max.list[['H3K27ac']] <- rep(8,length(list.H3K27ac.element.bed))
names(y.max.list[['H3K27ac']]) <- names(list.H3K27ac.element.bed)
y.max.list[['ATAC']] <- c(23,8,23,8,25,8,15,28,15,28,15,28)
names(y.max.list[['ATAC']]) <- names(list.H3K27ac.element.bed)
#########!!!ATAC without normalization#######

time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")
MyPlotColor(time.col,length(time.col))
ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)



pdf('G:/lab/Article/heshuang/BYLW/chromatin/Fig2de.withcolor.cg.raw.cluster1.cluster2.color.3.pdf',3,5)
for(tmp.cluster in names(list.H3K27ac.element.bed)[1:4]){
  for(his in c('H3K4me1','H3K27ac',
               'ATAC')){
    # y.max = round(max(list.nor.mean.df[[his]][c(list.H3K27ac.element.bed[[tmp.cluster]]$rowname,H3K27ac.pc1.p.up.bed$rowname),]))
    boxplot(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,],main= paste(tmp.cluster,his,sep = '.'),outline=F,
            col=ref.time.colors2[c('Virgin','G14.5','P7NL')],ylab='signal intensity',lwd=3,notch=T,
            ylim = c(0,y.max.list[[his]][tmp.cluster])
    )
    
    box(lwd=3)
    #abline(h=log2(cutoff.list$H3K4me1+1))
    
    P300.p <- wilcox.test(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,1],
                          list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,2],
                          #alternative = 'less',
                          paired = T)
    MyText(paste("G0CD_G14.5CD p-value: \n",
                 P300.p$p.value,"\n"),
           text.cex = 0.3)
    
    P300.p <- wilcox.test(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,2],
                          list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,3],
                          #alternative = 'greater',
                          paired = T)
    
    MyText(paste("G14.5CD_P7NL p-value: \n",
                 P300.p$p.value,"\n"),
           text.cex = 0.3)
    
  }
}

dev.off()


pdf('Fig21510.ATAC.pdf',3,5)
his <- 'ATAC';tmp.cluster <- 'up'
boxplot(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,],main= paste(tmp.cluster,his,sep = '.'),outline=F,col='#CECCCC',ylab='signal intensity',lwd=3,notch=T,
        ylim = c(0,18)
)
box(lwd=3)
P300.p <- wilcox.test(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,1],
                      list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,2],
                      #alternative = 'less',
                      paired = T)
MyText(paste("G0CD_G14.5CD p-value: \n",
             P300.p$p.value,"\n"),
       text.cex = 0.3)

P300.p <- wilcox.test(list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,3],
                      list.nor.mean.df[[his]][list.H3K27ac.element.bed[[tmp.cluster]]$rowname,][,2],
                      #alternative = 'less',
                      paired = T)
MyText(paste("G14.5CD_P7NL p-value: \n",
             P300.p$p.value,"\n"),
       text.cex = 0.3)
dev.off()
save.image('H3K27ac.CD.RData')

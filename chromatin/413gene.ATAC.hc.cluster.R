ATAC.merge500.bed <- ATAC.merge.bed[ATAC.merge.bed@ranges@width>=500]
G14.5HFD.up.bed <- ATAC.merge500.bed[ATAC.merge500.bed$G14.5HFD=='up']#506
G14.5HFD.500.df <-  as.data.frame(G14.5HFD.up.bed)[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'mean',sep = '.')]  



png(paste('20221203.412',".hc.clusters.raw.png",sep = ''),2000,3000)
par(oma = c(10,0,0,0))
G14.5HFD.row.tree <- MyHeatmap(as.data.frame(G14.5HFD.up.bed)[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'mean',sep = '.')],
          type = "log.raw",
          Colv = "none",
         # dist='euclidean',
          Rowv = "do",
          #r.cov.method = 'p',
          dendrogram = 'row',
          ColSideColors =
            MyName2Col(histone.si.tab[colnames(tmp.heatmap),]$Time,time.colors,is.row = F),
          ColSideColorsSize = 1,
          return.tree = 'row'
         # RowSideColors = MyName2Col(sort(cluster.raw.list[['G14.5HFD.500.res1']]),time.colors,is.row = T)#,
          # margins = c(15,15)
)
dev.off()

library(ape)
library(ggtree)
tmp.tree <- as.phylo(G14.5HFD.row.tree)

pdf('Fig6c.row.hc.tree.pdf',5,8)
ggtree(tmp.tree,size=2,right = F,ladderize = F)
dev.off()

heatmap.data <- as.data.frame(G14.5HFD.up.bed)[,paste(c('Ctrl','G0HFD','G14.5','G14.5HFD'),'mean',sep = '.')]
heatmap.data <- as.matrix(heatmap.data)

MyCalSi <- 
function(data,
         k,
         graph = F){
  library(cluster)
 # data.log2.relat <- as.matrix(log2(data + 1) / rowMax(log2(data + 1)))
  data.cor <- cor(t(log2(data+1)),
                  method = "p")
  data.dist <- as.dist(abs(1 - data.cor))
  data.hc <- hclust(data.dist, 
                    method = "ward.D2") 
  si <- silhouette(cutree(data.hc, k = k), data.dist)
  if(graph){
    plot(si, nmax = 80, cex.names = 0.5)
  }
  return(mean(si[,3]))
}

data.si.10 <- c(MyCalSi(heatmap.data,2),
                MyCalSi(heatmap.data,3),
                MyCalSi(heatmap.data,4),
                MyCalSi(heatmap.data,5),
                MyCalSi(heatmap.data,6))

pdf('silhouette score.analysis.pdf',5,5)
plot(2:6,data.si.10,pch=20,cex=3,type='o',lwd=5,xlab='number of clusters',ylab = 'Silhouette score')
abline(v=3,lty=2,lwd=5)
box(lwd=3)
dev.off()

MyWriteTable(data.si.10,
             "data.si.10.txt",
             col.names = F)



G14.5HFD.bed.list <- list()
G14.5HFD.bed.list[['cluster1']] <- labels(as.dendrogram(G14.5HFD.row.tree)[[1]])
G14.5HFD.bed.list[['cluster2']] <- labels(as.dendrogram(G14.5HFD.row.tree)[[2]][[1]])
G14.5HFD.bed.list[['cluster3']] <- labels(as.dendrogram(G14.5HFD.row.tree)[[2]][[2]])

sapply(G14.5HFD.bed.list, length)


order.new.list <- list()
#list.atac.element.bed <- list()
for(tmp.cluster in c('cluster1','cluster2','cluster3')){
  print(tmp.cluster)#tmp.cluster3 order by  
  order.new.list[[tmp.cluster]] <- G14.5HFD.bed.list[[tmp.cluster]][order(merge.df[G14.5HFD.bed.list[[tmp.cluster]],"G14.5HFD.mean"],decreasing = T)]
  list.atac.element.bed[[tmp.cluster]] <- ATAC.merge.bed[order.new.list[[tmp.cluster]]]
  # boxplot(as.data.frame(list.atac.element.bed[[tmp.cluster]])[,c('CD.H3K27ac.nor.mean','HFD.H3K27ac.nor.mean')],ylim = c(0,20),outline=F)
}


col.1 <- rep(gene.tree.colors[1],length(order.new.list$cluster1))
col.2 <- rep(gene.tree.colors[2],length(order.new.list$cluster2))
col.3 <- rep(gene.tree.colors[3],length(order.new.list$cluster3))
rColor <- t(as.matrix(c(col.1,
                        col.2,
                        col.3
                        # col.4
)))


MyPlotColor(c("white","#FFFDFA","#FFF7E6","#FFE3A9","goldenrod1","#F08104","#EF7A01","darkorange2"),8)
colors.exp2 <- colorRampPalette(c("midnightblue","dodgerblue3","#6D94DA","#F3F5FC","#FFF7E6","#FFE3A9","#FFC53C","goldenrod1","#FCB31E","darkorange2"),
                               space="Lab")(100)


MyPlotColor(c("midnightblue","dodgerblue3","#6D94DA","#F3F5FC","#FFF7E6","#FFE3A9","#FFC53C","goldenrod1","#FCB31E","darkorange2"),
            11
            )

library(viridis)
cor.heat.color <- autumn(100)
MyPlotColor(cor.heat.color,100)

pdf(paste('Figure/Fig6b.rawhc.cluster.',"20221203.color.412hcclusters.raw.pdf",sep = ''),15,22)
par(oma = c(10,0,0,0))
MyHeatmap(merge.df[c(order.new.list$cluster1,order.new.list$cluster2,order.new.list$cluster3),],
          color.palette = colors.exp2,
          type = "log.raw",
          Colv = "none",
          Rowv = "none",
          dendrogram = 'none',
          ColSideColors =
            MyName2Col(histone.si.tab[colnames(tmp.heatmap),]$Time,time.colors[c(1,2,8,13)],is.row = F),
          ColSideColorsSize = 1,
          RowSideColors = rColor#,
          # margins = c(15,15)
)
dev.off()

pdf(paste('Figure/Fig6b.rawhc.cluster.',"20221203.color.412hcclusters.log.row.zscore.2.pdf",sep = ''),15,22)
par(oma = c(10,0,0,0))
MyHeatmap(as.matrix(merge.df[c(order.new.list$cluster1,order.new.list$cluster2,order.new.list$cluster3),]),
          color.palette = colors.exp2,
         # type = "log.row.zscore",
          Colv = "none",
          Rowv = "none",
          dendrogram = 'none',
          ColSideColors =
            MyName2Col(histone.si.tab[colnames(tmp.heatmap),]$Time,time.colors[c(1,2,8,13)],is.row = F),
          ColSideColorsSize = 1,
          RowSideColors = rColor#,
          # margins = c(15,15)
)
dev.off()

time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")

ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)

ko.time.col <- c('gray80',ref.time.colors2['Virgin'],time.col[c(2)],time.col[3],
                 time.col[c(21,5,23)],ref.time.colors2['G14.5'],time.col[6],time.col[13],"#8ECFC9",'#FA7F6F',ref.time.colors2['P7NL'],
                 ref.time.colors2['P7L']
)

names(ko.time.col) <- c('/','Virgin','Virgin_HFD','G0_Acss2HFDWT','G0_Acss2HMKO_14d_HFD',
                        'G0Acss2HFDWT','G0Acss2HFDHMKO',
                        'G14.5','G14.5CDWT','G14.5_HFD','G14.5HFDWT','G14.5_Acss2HMKO_HFD','P7NL','P7NL_HFD'
                        
)

pdf('G:/lab/Article/heshuang/BYLW/chromatin/raw.cluster.Fig6c.color.refc412all.gene.hccluster500bp.ATAC.pdf',3,5)
for(tmp.cluster in names(G14.5HFD.bed.list)[1:3]){
  boxplot(as.data.frame(G14.5HFD.up.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,c('Ctrl.mean','G0HFD.mean','G14.5.mean','G14.5HFD.mean')],
          main = tmp.cluster,outline=F,col=add.alpha(ko.time.col[c('Virgin','Virgin_HFD','G14.5','G14.5_HFD')],1),ylab='signal intensity',lwd=3,notch=T,
          ylim = c(0,32),
          xaxt="n")
  box(lwd=3)
  axis(1,at = 1:2,labels = c('CD','HFD'),cex.axis = 0.5)
  
  P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,c('G14.5.mean')],
                        as.data.frame(G14.5HFD.up.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,c('Ctrl.mean')],
                        paired = T)

  MyText(paste("G0CD G14.5CD p-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)

  P300.p <- wilcox.test(as.data.frame(G14.5HFD.up.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,c('G14.5.mean')],
                        as.data.frame(G14.5HFD.up.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,c('G14.5HFD.mean')],
                        paired = T)

  MyText(paste("G14.5CD HFD p-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)
  
  
}
dev.off()


G14.5HFD.bed.list[['413']] <- G14.5HFD.up.bed$rowname
pdf('G:/lab/Article/heshuang/BYLW/chromatin/20230207.raw.Figs10k.refc.412all.gene.hccluster500bp.WTG0HFDG14.5CDG14.5HFDH3K27ac.pdf',3,5)
for(tmp.cluster in names(G14.5HFD.bed.list)[1:3]){
  boxplot(as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('G0HFD','CD','HFD'),'H3K27ac.nor.mean',sep = '.')],
          main= paste(tmp.cluster,sep = '.'),outline=F,col=add.alpha(ko.time.col[c('Virgin_HFD','G14.5','G14.5_HFD')],1),ylab='signal intensity',lwd=3,notch=T,
          ylim = c(0,22),xaxt="n"
  )
  axis(1,at = 1:3,labels = c('G0HFD','G14.5CD','G14.5HFD'),cex.axis = 0.5)
  box(lwd=3)
  P300.p <- wilcox.test(as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('CD'),'H3K27ac.nor.mean',sep = '.')],
                        as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('G0HFD'),'H3K27ac.nor.mean',sep = '.')],
                        paired = T)
  
  MyText(paste("G0HFD G14.5CDp-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)
  
  P300.p <- wilcox.test(as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('CD'),'H3K27ac.nor.mean',sep = '.')],
                        as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('HFD'),'H3K27ac.nor.mean',sep = '.')],
                        paired = T)
  
  MyText(paste("G14.5CD HFDp-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)
  
  
  P300.p <- wilcox.test(as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('HFD'),'H3K27ac.nor.mean',sep = '.')],
                        as.data.frame(ATAC.merge.bed[G14.5HFD.bed.list[[tmp.cluster]]])[,paste(c('G0HFD'),'H3K27ac.nor.mean',sep = '.')],
                        paired = T)
  
  MyText(paste("G0HFD G14.5HFDp-value: \n",
               P300.p$p.value,"\n"),
         text.cex = 0.5)
}
dev.off()

##########gene########
Fig6d.cluster2.gene <- unique(ATAC.merge.bed[G14.5HFD.bed.list[['cluster2']]]$Symbol)
Fig6d.cluster3.gene <- unique(ATAC.merge.bed[G14.5HFD.bed.list[['cluster3']]]$Symbol)
length(Fig6d.cluster2.gene)#81
length(Fig6d.cluster3.gene)#90
Fig6d.cluster2.gi <- genes.inf.input[c(Fig6d.cluster2.gene), ]
Fig6d.cluster2.gi$cluster <- 'cluster2'

Fig6d.cluster3.gi <- genes.inf.input[c(Fig6d.cluster3.gene), ]
Fig6d.cluster3.gi$cluster <- 'cluster3'

Fig6d.gi <- rbind(Fig6d.cluster2.gi,Fig6d.cluster3.gi)
MyWriteTable(Fig6d.gi,'Figure/Fig6d.cluster2.cluster3.gi.tab')
MyGOwritetable(Fig6d.cluster3.gi$EnsemblGeneID,pvalue = 1,'Figure/go.Fig6d.cluster3.tab')
MyGOwritetable(Fig6d.cluster2.gi$EnsemblGeneID,pvalue = 1,'Figure/go.Fig6d.cluster2.tab')

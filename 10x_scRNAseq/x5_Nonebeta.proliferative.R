setwd('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/Nbeta/')
dir.create('proliferation')
setwd('proliferation/')
load('Nbeta.cc.RData')

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
preg.colors <- c("#4DAF4A",
                 "#F781BF"
)
names(preg.colors) <- c("Ctrl",
                        "G14.5_1st")

all.list <- readRDS('../all.list.rds')
beta.cc.sym <- readRDS('G:/project/pregnant_mouse/10x/10x_v3/Ctl_G14.5_1st/beta.pre.cc.sym.rds')
beta.pre.cc.sym.list <- list()
cc.col.tree.list <- list()

i <- 'pp';num <- 3
#num=1;
num <- 3;
for(i in names(all.list)){
  beta.pre.cc.sym.list[[i]] <- rownames(MyGeneExp(as.matrix(exp(all.list[[i]]@assays$RNA@data)[beta.cc.sym,]),
                                        2,exp = 'yes',num))
  c1Color <- MyName2Col(all.list[[i]]@meta.data[,"Time"],
                        preg.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(all.list[[i]]@meta.data[,"hc.proliferative"],
                        preg.colors)
  c2Color <- as.matrix(c2Color)
  cColor <- cbind(c1Color,c2Color)
  png(paste(i,'.expin3.beta.cc.png',sep=''),2000,3000)
  cc.col.tree.list[[i]] <-
    MyHeatmap(as.matrix(exp(all.list[[i]]@assays$RNA@data))[beta.pre.cc.sym.list[[i]],],
              type = "log.row.relat",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = cColor,
              ColSideColorsSize = 2,
              return.tree = "col",
              graph = T)
  dev.off()
}
hc.cc.list <- list()
hc.cc.list[['pp']] <- c()
alpha.cc.den <- as.dendrogram(cc.col.tree.list[['alpha']])
hc.cc.list[['alpha']] <- labels(alpha.cc.den[[1]][[1]])

delta.cc.den <- as.dendrogram(cc.col.tree.list[['delta']])
hc.cc.list[['delta']] <- c(labels(delta.cc.den[[1]]))

for(i in names(all.list)){
  all.list[[i]]$hc.proliferative <- 'quiescent'
  all.list[[i]]@meta.data[hc.cc.list[[i]],'hc.proliferative'] <- 'proliferative'
}
##########
for(i in c(
  'alpha',
  'pp',
  'delta'
)){
  all.list[[i]]$hc.proliferative <- factor(all.list[[i]]$hc.proliferative,levels = c('quiescent','proliferative'))
  
  all.list[[i]] <- ScaleData(all.list[[i]],features = beta.pre.cc.sym.list[[i]])
  all.list[[i]] <- RunPCA(all.list[[i]],features = beta.pre.cc.sym.list[[i]])
  # pdf(paste(i,'.pca.cc.marker.pdf'),6,7)
  # Myseuratmarker(all.list[[i]],marker.sym = c("Cdk2","Ccne1","Mcm2","E2f7",      #G1/S
  #                                             "Cdk1","Ccnb1","Kif4",             #G2/M
  #                                             "Cdkn1c","Cdkn2c","Cdkn1b",
  #                                             'Cdt1','Gmnn','Top2a','Mki67'),reduction = 'pca')
  # dev.off()
  
  pdf(paste(i,'.exp3.hc.pc12.proliferative.pdf',sep = ''),13,12)
  p.pca <- MySeuratDR2Gg2(all.list[[i]],all.list[[i]]@meta.data,
                          reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
  p.pca$data_ <- p.pca$data
  p.pca$data_$hc.proliferative_ <- p.pca$data_$Time
  p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(all.list[[i]]$Time))))
  p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]
  
  
  plot(p.pca+
         scale_color_manual(values =  c('cyan3','red3')) +
         scale_shape_manual(values = c(19,17)) +
         geom_point(aes(x =-x.pos,
                        y = y.pos,
                        col = Time,
                        shape = hc.proliferative
         ),
         size = 6
         )+
         theme(aspect.ratio=1))
  
  # p.pca$data_ <- p.pca$data
  # p.pca$data_$hc.proliferative_ <- p.pca$data_$knn.hc.proliferative
  # p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(all.list[[i]]$knn.hc.proliferative))))
  # p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]
  
  # plot(p.pca+
  #        scale_color_manual(values = c("#ba984d",'gray90')) +
  #        geom_point(aes(x =x.pos,
  #                       y = y.pos,
  #                       col = hc.proliferative#,
  #                       # shape = State
  #        ),
  #        size =6
  #        )+
  #        theme(aspect.ratio=1))
   dev.off()
}
########plot marker######
for(i in c(
 # 'alpha',
  'pp',
  'delta'
)){
  p.pca <- MySeuratDR2Gg2(all.list[[i]],all.list[[i]]@meta.data,
                          reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
  p.pca$data_ <- p.pca$data
  p.pca$data_$hc.proliferative_ <- p.pca$data_$Time
  p.pca$data_$hc.proliferative_ <- factor(p.pca$data_$hc.proliferative_,levels = rev(names(table(all.list[[i]]$Time))))
  p.pca$data <- p.pca$data_[order(p.pca$data_$hc.proliferative_),]
  
  p.loop <- p.pca +
    scale_color_gradientn(colours = colors.exp) +
    #scale_shape_manual(values = 1:5)
    theme(plot.title = element_text(size = rel(3.5))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_blank()) +
    theme(axis.text  = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(legend.position="bottom")+
    theme(panel.border = element_rect(size = 4,
                                      colour = "black"))+
    guides(colour = guide_colorbar(title = "ln(TP10K + 1)",
                                   title.position = "top",
                                   barwidth = 20,
                                   title.hjust = 0.5,
                                   title.theme = element_text(angle = 0,
                                                              size = 20),
                                   label.theme = element_text(angle = 0,
                                                              size = 20),
                                   ticks = T)) 
  cellcycle.marker.set=c("Cdk2","Ccne1","Mcm2","E2f7",      #G1/S
                         "Cdk1","Ccnb1","Kif4",             #G2/M
                         "Cdkn1c","Cdkn2c","Cdkn1b",
                         'Cdt1','Gmnn','Top2a','Mki67')   #p57(Cdkn1c)p18(Cdkn2c)p16(Cdkn2a)p27(Cdkn1b) inhibitor
  
  pdf(paste(i,".exp3.hc.knn.cc.marker.cellcycle.pdf",sep = ''),
      6,
      7)
  for (plot.ens in cellcycle.marker.set){
    print(plot.ens)
    print(p.plot <- p.loop +
            geom_point(aes(x = x.pos,
                           y = y.pos,
                           
                           color = all.list[[i]]@assays$RNA@data[plot.ens,rownames(p.pca$data)]
                           #shape = Age
            ),size = 5) +
            labs(title = plot.ens)+
            theme(aspect.ratio=1))
    
  }
  
  dev.off()
}
######
rm(all.list)
save.image('Nbeta.cc.RData')

setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/KOdirect/20220913/20221106/20221121/')
source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")
dir.create('vio')
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(Vennerable)
ko.time.col <- c('gray90',time.colors[c(1,8,10,4,c(25:26,25:26,15,30),25:28,25:28)])

names(ko.time.col) <- c('/','Virgin','G14.5','G18.5',
                        'Virgin_Stat3HMKO',
                        'G18.5_Stat3WT','G18.5_Stat3HMKO',
                        'G14.5_Acss2WT','G14.5_Acss2HMKO',
                        'G18.5_Acss2WT','G18.5_Acss2HMKO',
                        'Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif',
                        'DMSO-human-serum-Ctrl','A485-human-serum-Ctrl','DMSO-human-serum-Preg','A485-human-serum-Preg'
                        
)
seu.ref.KO.list <- readRDS('seu.ref.KO.list.rds')
seu.G18.5Stat3 <- subset(seu.ref.KO.list$Stat3HMKO,cells = colnames(seu.ref.KO.list$Stat3HMKO)[seu.ref.KO.list$Stat3HMKO$Type %in% c('G18.5_Stat3WT','G18.5_Stat3HMKO')])
seu.G18.5Acss2 <- subset(seu.ref.KO.list$Acss2HMKO,cells = colnames(seu.ref.KO.list$Acss2HMKO)[seu.ref.KO.list$Acss2HMKO$Type %in% c('G18.5_Acss2WT','G18.5_Acss2HMKO')])
seu.G18.5Acss2$Type <- factor(seu.G18.5Acss2$Type,levels = c('G18.5_Acss2WT','G18.5_Acss2HMKO')
                              )

seu.G18.5Acss2 <- subset(seu.ref.KO.list$Acss2HMKO,cells = colnames(seu.ref.KO.list$Acss2HMKO)[seu.ref.KO.list$Acss2HMKO$Type %in% c('Virgin','G18.5_Acss2WT','G18.5_Acss2HMKO')])
seu.G18.5Acss2$Type <- factor(seu.G18.5Acss2$Type,levels = c('Virgin','G18.5_Acss2WT','G18.5_Acss2HMKO')
)

rm(seu.ref.KO.list)

sym.list <- c('Rfx6','Nr1d1','Pdx1','Nkx6-1','Neurod1','Gck','Slc2a2','Hadh','Ucp2')

seu.G18.5Stat3 <- readRDS('vio/seu.G18.5Stat3.rds')
seu.G18.5Stat3$Type <- factor(seu.G18.5Stat3$Type,levels = c('G18.5_Stat3WT','G18.5_Stat3HMKO'))


pdf('vio/vioplot.Chgb.Stat5.tp0.1m.log.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G18.5Stat3,sym.list = c('Stat5a','Stat5b','Chgb','Ovol2','Ttr'),type = 'Type',type.colors = time.colors[c(10,8)],log = T)
dev.off()

seu.G18.5Acss2 <- readRDS('../../vio/seu.G18.5Acss2.rds')

pdf('G:/lab/Article/heshuang/BYLW/sm3/KO/Acss2KO.vioplot.tp0.1m.log.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G18.5Acss2,sym.list = sym.list,type = 'Type',type.colors = c("#A65628", "#5F9EA0"),log = T)
dev.off()


saveRDS(seu.G18.5Acss2,'vio/seu.G18.5Acss2.rds')
saveRDS(seu.G18.5Stat3,'vio/seu.G18.5Stat3.rds')
###########plot pseudotime#########
MyGetTrend2 <- 
  function(data.pseudotime,
           data.rpkm,
           min.rpkm = 1,
           trend.formula = "log2rpkm ~ sm.ns(pseudotime, df=2)",...){
    data <- data.frame(pseudotime = data.pseudotime,
                       rpkm = data.rpkm)
    data$log2rpkm <- data$rpkm
    expectation <- tryCatch({vg <- suppressWarnings(vgam(formula = as.formula(trend.formula), 
                                                         family = tobit(), 
                                                         data = data, 
                                                         maxit = 30, 
                                                         checkwz = FALSE))
    res <-  predict(vg, 
                    type = "response")
    res
    }, 
    error = function(e) {
      print("Error!")
      print(e)
      res <- rep(NA, nrow(data))
      res
    })
    data.out <- data.frame(pseudotime = data.pseudotime, 
                           expectation = expectation)
    return(data.out)
  }



Mygene2pseudotime <- 
  function(genes.inf.input = genes.inf.input ,
           gene.list,
           tpm.data,
           Pseudotime,
           Time,
           sample.list,
           log = T,
           point.size=3,
           cols=liver.col,
           continu=F,...
           
  ){
    library(VGAM)
    tmp.list <- gene.list
    
    if(log==T){
      time.data <- data.frame(Pseudotime = Pseudotime,
                              Time = Time,
                              log2(1+t(tpm.data[tmp.list,sample.list])))
    }else{
      time.data <- data.frame(Pseudotime = Pseudotime,
                              Time = Time,
                              t(tpm.data[tmp.list,sample.list]))
      
    }
    ylab <- ""
    xlab <- ""
    tmp.list <- sub('-','.',tmp.list)
    total.count = length(tmp.list)
    run.count = 1
    if(continu==T){
      p.loop.1 <- ggplot(data = time.data) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradientn(colours = cols) 
    }
    else{
      p.loop.1 <- ggplot(data = time.data) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_manual(values = cols) 
    }
    p.loop <-  p.loop.1+
      scale_shape_manual(values = c(16,1)) +
      guides(colour = "none") +
      guides(shape = "none") +
      theme(axis.text = element_text(size = 20, colour = "black")) +
      theme(axis.title.x = element_text(size = 25, colour = "black")) +
      theme(axis.title.y = element_text(size = 25, colour = "black")) +
      theme(plot.title = element_text(size = 30)) +
      theme(panel.border = element_rect(size = 2,
                                        colour = "black")) +
      theme(axis.line = element_line(colour = "black")) +
      ylab(ylab) +
      xlab(xlab) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    
    for(ens in tmp.list){
      cat(paste(run.count, "/",total.count,"\n"))
      p.plot <- p.loop +
        geom_point(aes(x = as.numeric(Pseudotime),
                       y = as.numeric(time.data[,ens]),
                       col = Time#,
                       # shape = CellCycle
        ),
        size = point.size) +
        labs(title = ens) +
        ylim(min(as.numeric(time.data[,ens])),NA) +
        geom_line(aes(x = pseudotime,
                      y = expectation + 0),
                  col = "black",
                  size = 2,
                  data = MyGetTrend2(time.data$Pseudotime,
                                     time.data[,ens],
                                     min.rpkm = 1))
      print(p.plot)
      
      run.count = run.count + 1
    }
  }

##########
seu.ref.KO.list <- readRDS('seu.ref.KO.list2.rds')
dir.create('gene_trend')

pdf('gene_trend/G0G18.5Stat3WTHMKOtype.color.logTP0.1M.pdf',6,6)
#par(mfrow = c(3,4))
Mygene2pseudotime(gene.list = c('Chgb','Ovol2','Ttr'),tpm.data = seu.ref.KO.list$G0G18.5Stat3WTHMKO@assays$RNA@data,
                  Pseudotime = seu.ref.KO.list$G0G18.5Stat3WTHMKO$pseudotime,
                  Time = seu.ref.KO.list$G0G18.5Stat3WTHMKO$Type,sample.list = colnames(seu.ref.KO.list$G0G18.5Stat3WTHMKO),log = F,cols = time.colors[c(5,6,10,8)],continu = F,point.size = 4)
dev.off()

pdf('gene_trend/G0G14.5P300WTHMKO.type.color.logTP0.1M.pdf',6,6)
#par(mfrow = c(3,4))
Mygene2pseudotime(gene.list = c('Chgb','Ovol2','Ttr'),tpm.data = seu.ref.KO.list$G0G14.5P300WTHMKO@assays$RNA@data,
                  Pseudotime = seu.ref.KO.list$G0G14.5P300WTHMKO$pseudotime,
                  Time = seu.ref.KO.list$G0G14.5P300WTHMKO$Type,sample.list = colnames(seu.ref.KO.list$G0G14.5P300WTHMKO),log = F,cols = time.colors[c(5,6,7,8)],continu = F,point.size = 4)
dev.off()

pdf('gene_trend/G0G18.5Acss2WTHMKO.type.color.logTP0.1M.pdf',6,6)
#par(mfrow = c(3,4))
Mygene2pseudotime(gene.list = c('Chgb','Ovol2','Ttr'),tpm.data = seu.ref.KO.list$G0G18.5Acss2WTHMKO@assays$RNA@data,
                  Pseudotime = seu.ref.KO.list$G0G18.5Acss2WTHMKO$pseudotime,
                  Time = seu.ref.KO.list$G0G18.5Acss2WTHMKO$Type,sample.list = colnames(seu.ref.KO.list$G0G18.5Acss2WTHMKO),log = F,cols = time.colors[c(1,2,15,12)],continu = F,point.size = 4)
dev.off()

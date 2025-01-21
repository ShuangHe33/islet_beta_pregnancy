preg.gene.tab <- MyReadDelim('../../Glut2H/LMgene/gene_inf/cg.preg.gene.si.1e28.addcc.tab')
MyGOwritetable(preg.gene.tab$EnsemblGeneID[preg.gene.tab$cluster=='preg.up'],'../../Glut2H/LMgene/go_kegg/go.preg.up.tab')


dir.create('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/cluster/vio')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/cluster/vio')
library(future)
plan("multisession", workers = 10)
plan()

source("G:/pcatest/MyFunction.R")
source('G:/r script/file_function.R')
source("G:/r script/MyFunction.Seurat3.R")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(Vennerable)

####
Myseuviopvalplot <- 
  function(seu.ob,sym.list,type = 'Type',type.colors = time.colors,log=F,gap=1,...){
    beta.preg.meta.tab <- seu.ob@meta.data
    tmp.median <- matrix(ncol = (length(names(table(seu.ob@meta.data[,type])))+2),nrow = length(sym.list))
    rownames(tmp.median) <- sym.list
    colnames(tmp.median) <- c('p-value.1','p-value.2',names(table(seu.ob@meta.data[,type])))
    
    if(log==F){
      for(j in sym.list){
        print(j)
        beta.preg.meta.tab$tmp <- as.matrix(seu.ob@assays$RNA@data)[j,beta.preg.meta.tab$SampleName]
        tmp.tpm <- exp(beta.preg.meta.tab$tmp)-1
        names(tmp.tpm) <- beta.preg.meta.tab[,type]
        MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                            tmp.tpm,
                            color.violin = add.alpha(c(#"#1f78b4",
                              type.colors),alpha = 0.7),
                            box.lwd = 3,
                            ylab.plot  = "TP0.1M",
                            gap=gap,
                            j,...
                            
                            
        )
        tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'],
                             beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
                             paired = F)
        tmp.median[j,'p-value.1'] <- tmp.p$p.value
        tmp.median[j,names(table(seu.ob@meta.data[,type]))[1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'])-1),2)
        tmp.median[j,names(table(seu.ob@meta.data[,type]))[2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'])-1),2)
        
        tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'],
                             beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
                             paired = F)
        tmp.median[j,'p-value.2'] <- tmp.p$p.value
        tmp.median[j,names(table(seu.ob@meta.data[,type]))[3]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'])-1),2)
        
      }
    }
    
    if(log==T){
      for(j in sym.list){
        print(j)
        beta.preg.meta.tab$tmp <- as.matrix(seu.ob@assays$RNA@data)[j,beta.preg.meta.tab$SampleName]
        tmp.tpm <- beta.preg.meta.tab$tmp
        names(tmp.tpm) <- beta.preg.meta.tab[,type]
        MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                            tmp.tpm,
                            color.violin = add.alpha(c(#"#1f78b4",
                              type.colors),alpha = 0.7),
                            box.lwd = 3,
                            ylab.plot  = "ln(TP0.1M+1)",
                            gap=gap,
                            j,...
                            
                            
        )
        tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'],
                             beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
                             paired = F)
        tmp.median[j,'p-value.1'] <- tmp.p$p.value
        tmp.median[j,names(table(seu.ob@meta.data[,type]))[1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'])-1),2)
        tmp.median[j,names(table(seu.ob@meta.data[,type]))[2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'])-1),2)
        
        # tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'],
        #                      beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
        #                      paired = F)
        # tmp.median[j,'p-value.2'] <- tmp.p$p.value
        # tmp.median[j,names(table(seu.ob@meta.data[,type]))[3]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'])-1),2)
        
      }
    }
    return(tmp.median)
  }

#####
seu.G0G14.5 <- readRDS('seu.G0G14.5.rds')
seu.beta <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/cluster/seu.beta.ges.post.rds')
seu.G0G14.5 <- subset(seu.beta,cells = colnames(seu.beta)[seu.beta$Type %in% c('Virgin','G14.5')])



seu.ref.raw <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/cluster/seu.ref.raw.beta.celltype.rds')
seu.ref.raw$Type <- as.character(seu.ref.raw$Type)

seu.G0G14.5 <- subset(seu.ref.raw,cells = colnames(seu.ref.raw)[seu.ref.raw$Type %in% c('Virgin','G14.5')])
table(seu.G0G14.5$Type,seu.G0G14.5$Rep)
seu.G0G14.5$Type <- factor(seu.G0G14.5$Type,levels = c('Virgin','G14.5'))


time.col <- c("#b6d4a8","#8c9864",
              "#adba7d","#9CBD13","#e9de9b","#c2bf3c","#e1d554","#e8bb78",
              "#e19368","#b08061","#c2bb9f","#B2A59C","#79706A",
              "#b89fc1","#d8c0ef","#7CBAB5","#e29edd","#919ac5","#b7c9e7","#c3aeca","#919ac5",
              "#ecc44e","#87afcc","#e09694","#88b494","#9ED3AD","#4EA25D","#76876F","#b09694","#c08677")
celltype.col2 <- time.col[c(2,9,17,12)]
names(celltype.col2) <- c('Glut2L_1','Glut2L_2','Glut2H_1','Glut2H_2')

ref.time.colors2 <- time.col[c(1,4,5,9:12,14,13,15,16:18,2,27,25)]
names(ref.time.colors2) <- names(ref.time.colors)

Myseuviopvalplot <- 
function(seu.ob,sym.list,type = 'Type',type.colors = time.colors,log=F,gap=1,alpha.col=0.7,...){
  beta.preg.meta.tab <- seu.ob@meta.data
  tmp.median <- matrix(ncol = (length(names(table(seu.ob@meta.data[,type])))+2),nrow = length(sym.list))
  rownames(tmp.median) <- sym.list
  colnames(tmp.median) <- c('p-value.1','p-value.2',names(table(seu.ob@meta.data[,type])))
  
  if(log==F){
    for(j in sym.list){
      print(j)
      beta.preg.meta.tab$tmp <- as.matrix(seu.ob@assays$RNA@data)[j,beta.preg.meta.tab$SampleName]
      tmp.tpm <- exp(beta.preg.meta.tab$tmp)-1
      names(tmp.tpm) <- beta.preg.meta.tab[,type]
      MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                          tmp.tpm,
                          color.violin = add.alpha(c(#"#1f78b4",
                            type.colors),alpha = alpha.col),
                          box.lwd = 3,
                          ylab.plot  = "TP0.1M",
                          gap=gap,
                          j,...
                          
                          
      )
      tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'],
                           beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
                           paired = F)
      tmp.median[j,'p-value.1'] <- tmp.p$p.value
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'])-1),2)
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'])-1),2)
      
      # tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'],
      #                      beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
      #                      paired = F)
      # tmp.median[j,'p-value.2'] <- tmp.p$p.value
      # tmp.median[j,names(table(seu.ob@meta.data[,type]))[3]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'])-1),2)
      
    }
  }
  
  if(log==T){
    for(j in sym.list){
      print(j)
      beta.preg.meta.tab$tmp <- as.matrix(seu.ob@assays$RNA@data)[j,beta.preg.meta.tab$SampleName]
      tmp.tpm <- beta.preg.meta.tab$tmp
      names(tmp.tpm) <- beta.preg.meta.tab[,type]
      MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                          tmp.tpm,
                          color.violin = add.alpha(c(#"#1f78b4",
                            type.colors),alpha = alpha.col),
                          box.lwd = 3,
                          ylab.plot  = "ln(TP0.1M+1)",
                          gap=gap,
                          j,...
                          
                          
      )
      tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'],
                           beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
                           paired = F)
      tmp.median[j,'p-value.1'] <- tmp.p$p.value
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[1],'tmp'])-1),2)
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'])-1),2)
      
      # tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'],
      #                      beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[2],'tmp'],
      #                      paired = F)
      # tmp.median[j,'p-value.2'] <- tmp.p$p.value
      # tmp.median[j,names(table(seu.ob@meta.data[,type]))[3]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[3],'tmp'])-1),2)
      
    }
  }
  return(tmp.median)
}

sym.list <- c('Chgb','Ovol2','Ttr')
pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/vioplot.gene.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Eif4ebp1','Akt1','Gsk3b','Hspa5','Hspa13','Eif2s1','Atf6','Ern1','Chgb','Ovol2','Ttr'),type = 'Type',type.colors = ref.time.colors2[ c('Virgin','G14.5')],log = T,alpha.col = 1,med.lwd=8)
dev.off()

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/GDMgene.vioplot.gene.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Slc30a8','Rbp4','Gck','Ttr'),type = 'Type',type.colors = ref.time.colors2[ c('Virgin','G14.5')],log = T,alpha.col = 1,med.lwd=8)
dev.off()

pdf('vioplot.Akt1.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Eif4ebp1','Akt1','Gsk3b'),type = 'Type',type.colors = ref.time.colors2[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()



pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/Acss2.vioplot.gene.logF.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Acss2','Acly','Pdhb','Pdha1','Pdhx'),type = 'Type',type.colors = ref.time.colors2[ c('Virgin','G14.5')],log = F,alpha.col = 1,med.lwd=8)
dev.off()


tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Eif4ebp1','Akt1','Gsk3b'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)

sym.list <- c('Ep300','Crebbp','Gcn5')
pdf('vioplot.Akt1.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Eif4ebp1','Akt1','Gsk3b'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()


pdf('Glut2H.G0G14.5.vioplot.Glut2.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Slc2a2'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()

seu.G0G14.5 <- SetIdent(seu.G0G14.5,value = seu.G0G14.5$Type)

tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Slc2a2'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)




GDM.symbol <- c("Cdkal1",
                "Tcf7l2",
                "Gck",
                "Kcnj11",
                "Kcnq1",
                "Igf2bp2",
                "Irs1",
                "Tspan8",
                "Gckr",
                "Mtnr1b",
                "Ide")

sym.list <- c('Cdkn2a','Cdkn2b',#
              'Slc30a8','Arap1',#'Centd2',
              'Rbp4','Slc2a4','Pck1','Pik3r1',
              'Synpr','Cdh18','Ctif','Ptgis'
              )

sym.list <-c('Lipf','Lipk','Lipn',#'Lipj',
             'Tph1','Tph2','Apc','Gsk3b','Brca1','Sycp2','Clock','Grin3b','Nr3c1','Prdm16','Rnls','Sall3','Slc12a8','Tmem259')


sym.list <- c('Slc30a8','Rbp4',#"Kcnj11",
              "Gck",#'Tmem259',
              'Tph1','Tph2','Synpr','Gsk3b'
              )

sym.list <- c(#'Ddit3',#'Eif2ak3',
              #'Eif2a',
              'Eif2s1',
              'Atf4','Atf6',#'Atf6b',
              'Ern1',#'Xbp1',
              'Hspa5',#'Bhlha15''
              'Hspa13'
              )

pdf('selected.Figs2e.GWAS.vioplot.GDMgene.log.2.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = sym.list,type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()
pdf('selected.Erstress.vio.log.2.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = sym.list,type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()

c('Gckr','Spc25','Adcy5','Pcsk1','Esr1','Ccnd2','Nedd1','Cmip','Map3k15')[! c('Gckr','Spc25','Adcy5','Pcsk1','Esr1','Ccnd2','Nedd1','Cmip','Map3k15') %in% rownames(seu.G0G14.5)]

pdf('GDM.NG2024.vio.log.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c('Gckr','Spc25','Adcy5','Pcsk1','Esr1','Ccnd2','Nedd1','Cmip','Map3k15'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = T,med.lwd=8)
dev.off()


MyViolinBeeSwarmMed <- 
function(x,
         y,
         color.violin = "white",
         title.plot = "",
         xlab.plot = "",
         ylab.plot = "",
         box.lwd = 4,
         cex.main = 3,
         corral.beeswarm = "wrap",
         pch.beeswarm = 20,
         point.col = "#666666",
         med.length = 0.2,
         med.lwd = 5,
         x.cex.axis = 1.5,
         y.cex.axis = 3,
         ylim = c(0,0),
         ext=0,
         ...){
  library(vioplot)
  library(beeswarm)
  group.all <- levels(x)
  # par(oma = c(0,3,0,0))
  plot(1:2,1:2,col = "white",
       xlab = xlab.plot,
       ylab = ylab.plot,
       xlim = c(0.5,length(group.all)+0.5),
       ylim = c(0,15),
       xaxt="n",
       yaxt="n",
       main = title.plot,
       cex.main = cex.main,
       
       ...)
  box(lwd = box.lwd)
  axis(1,lwd = 0,lwd.ticks = 0, at = 1:length(group.all),cex.axis = y.cex.axis,labels = group.all)
  axis(2,cex.axis = y.cex.axis)
  
  for (n.plot in 1:length(group.all)){
    if (length(unique(y[x == group.all[n.plot]])) > 1){
      vioplot(y[x == group.all[n.plot]],
              at = n.plot,
              drawRect = F,
              # rectCol = "white",
              # colMed = "black",
              # pchMed = 1,
              col = color.violin[n.plot],
              add = T,
              ...)
    }
  }
  beeswarm(y ~ x,
           col = point.col,
           corral = corral.beeswarm,
           pch = pch.beeswarm,
           add = TRUE,
           ...)
  for (n.plot in 1:length(group.all)){
    if (length(unique(y[x == group.all[n.plot]])) > 1){
      data.med <- median(y[x == group.all[n.plot]])
      # plot(c(n.plot-0.2,n.plot+0.2),
      #      c(data.med,data.med),
      #      type = "l",
      #      add = T)
      lines(c(n.plot-med.length/2,n.plot+med.length/2),
            c(data.med,data.med),
            lwd = med.lwd)
    }
  }
}

pdf('STAT.vioplot.2.pdf',12,15)
par(mfrow = c(3,4))
Myseuviopvalplot(seu.G0G14.5,sym.list = c('Stat1','Stat2','Stat3','Stat4','Stat5a','Stat5b','Stat6'),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = F,med.lwd=8)
dev.off()


pdf('mTOR.ref.vioplot.2.pdf',12,15)
par(mfrow = c(3,4))
Myseuviopvalplot(seu.G0G14.5,sym.list = c('Foxo1','Gsk3a','Irs1','Irs2','Gsk3b','Akt1','Akt2','Akt3',
                                          grep('^Pik3',genes.inf.input$Symbol,value = T)
                                          
                                          ),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = F,med.lwd=8)
dev.off()
pdf('mTOR.ref.vioplot.down.pdf',12,15)
par(mfrow = c(3,4))
Myseuviopvalplot(seu.G0G14.5,sym.list = c("Cab39",
                                          "Pik3cb",
                                          "Cab39l",
                                          "Pik3r3",
                                          "Pik3r2",
                                          "Ulk3",
                                          "Hif1a",
                                          "Braf",
                                          "Eif4ebp1",
                                          "Tsc2",
                                          "Ulk1",
                                          "Ddit4",
                                          "Pdpk1",
                                          "Akt2",
                                          "Rictor",
                                          "Prkaa1",
                                          "Tsc1",
                                          "Vegfb",
                                          "Pgf",
                                          "Vegfa"),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = F,med.lwd=8)


dev.off()

grep('^Pik3',genes.inf.input$Symbol,value = T)
grep('Gsk3',genes.inf.input$Symbol,value = T)






MyWriteTable(tp0.1M.pval,row.names = T,'Chgb.tp0.1.m.pval.tab')

sym.list <- c('Acss2','Acly')
pdf('Fig5a.vioplot.Acss2.Acly.tp0.1m.pdf',12,15)
par(mfrow = c(3,4))
tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = sym.list,type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = F,ext=1,med.lwd=8)
dev.off()
MyWriteTable(tp0.1M.pval,row.names = T,'Acss2Acly.tp0.1.m.pval.tab')
saveRDS(seu.G0G14.5,'seu.G0G14.5.rds')

# 
# ensyme.list <- list()
# ensyme.list[['KAT.list']] <- c('Ep300','Crebbp',grep('Kat',genes.inf.input$Symbol,value = T))
# ensyme.list[['KAT.list']] <- ensyme.list[['KAT.list']][1:13]
# ensyme.list[['HDAC.list']] <- grep('Hdac',genes.inf.input$Symbol,value = T)
# 
# pdf('vioplot.kat.tp0.1m.pdf',15,15)
# par(mfrow = c(3,4))
# tp0.1M.pval <- Myseuviopvalplot(seu.G0G14.5,sym.list = c(ensyme.list$KAT.list,ensyme.list$HDAC.list),type = 'Type',type.colors = ref.time.colors[ c('Virgin','G14.5')],log = F)
# dev.off()
# MyWriteTable(tp0.1M.pval,row.names = T,'Kat.tp0.1.m.pval.tab')
#######post#########
seu.ref.raw <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/cluster/seu.ref.raw.beta.celltype.rds')
seu.ref.raw <- subset(seu.ref.raw,cells = colnames(seu.ref.raw)[seu.ref.raw$Type %in% c('Virgin','G14.5','P7NL')])
seu.ref.raw$Type <- factor(seu.ref.raw$Type,levels = c('Virgin','G14.5','P7NL'))

rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)
G0G14.5P7NL.tp0.1m <- do.call(cbind,by(t(as.matrix(exp(seu.ref.raw@assays$RNA@data)-1)),seu.ref.raw@meta.data$Type,colMeans))
G0G14.5P7NL.tp0.1m2 <- cbind(genes.inf.input[rownames(G0G14.5P7NL.tp0.1m),c('EnsemblGeneID','Symbol')],G0G14.5P7NL.tp0.1m)
MyWriteTable(G0G14.5P7NL.tp0.1m2,'G0G14.5P7NL.tp0.1m.20221203.tab')


ensyme.list <- list()

tp0.1M.pval.list <- list()

ensyme.list[['HDAC.list']] <- grep('Hdac',genes.inf.input$Symbol,value = T)
ensyme.list[['Sirt.list']] <- grep('Sirt',genes.inf.input$Symbol,value = T)
ensyme.list[['Kmt.list']] <- grep('Kmt',genes.inf.input$Symbol,value = T)#methyltransferase
ensyme.list[['Setd.list']] <- grep('Setd',genes.inf.input$Symbol,value = T)#methyltransferase

ensyme.list[['HDAC.list']] <- c(paste('Hdac',c(3,4,5,6,10,11),sep = ''),paste('Sirt',c(3,6,7),sep = ''))

logtpm.df <- seu.ref.raw@assays$RNA@data[ensyme.list[['HDAC.list']],]
min.list<- apply(logtpm.df,1,max)

ensyme.list[['p300']] <- c('Kat2a','Ep300')



for(tmp.list in names(ensyme.list)){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/histone/',tmp.list,'.6.tp0.1m.pdf',sep = '.'),12,15)
  par(mfrow = c(3,4))
  tp0.1M.pval.list[[tmp.list]] <- Myseuviopvalplot(seu.ref.raw,sym.list = ensyme.list[[tmp.list]],#c('Hdac3','Hdac5','Hdac10','Sirt7'),
                                                   type = 'Type',type.colors = ref.time.colors2[ c('Virgin','G14.5','P7NL')],log = T,y.cex.axis=1,ext=0,alpha.col = 1,gap = 0.5)
  dev.off()
}
MyWriteTable( tp0.1M.pval.list[[tmp.list]],row.names = T,'Sirt.pval.tab')




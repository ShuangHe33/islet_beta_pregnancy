MyGosheep <-
  function(selected.id,
           universe.id,
           pvalue=1,
           qvalue=1,
           GOanno=oaries_Entrez_ID2goID.res,
           readable=TRUE,
           ont='BP',
           ...){
    library(clusterProfiler)
    
    TERM2GENE=GOanno[GOanno$ONTOLOGY%in%ont,1:2]
    TERM2NAME=GOanno[GOanno$ONTOLOGY%in%ont,c('go_id','term')]
    
    goresult <- enricher(
      gene = selected.id,
      universe=as.character(universe.id),
      TERM2GENE=TERM2GENE,
      TERM2NAME=TERM2NAME,
      pvalueCutoff  = pvalue,
      qvalueCutoff  = qvalue
    )
    return(goresult)
  }

MyGOsheepwritetable <-
  function(genelist,
           go.out,
           universe.gene=na.omit(oaries_Entrez_ID2goID.res$entrezgene_id),
           pvalue = 1,
           qvalue = 1,
           ont=names(table(oaries_Entrez_ID2goID.res$ONTOLOGY)),
           GOanno=oaries_Entrez_ID2goID.res,
           ...){
    go.result <- MyGosheep(genelist,
                           universe.gene,
                           pvalue = pvalue,
                           qvalue = qvalue,
                           ont=ont,
                           GOanno=GOanno,
                           ...)
    
    MyNCBI2Symbol <- function(mygoid,genelist){
      goid2symbol <- 1:length(mygoid)
      for(i in 1:length(mygoid)){
        goid <- strsplit(mygoid[i],"/")
        goid <- goid[[1]]
        mygoid2symbol <- genes.inf.input[goid,'Symbol']
        mygoid2symbol <- paste(mygoid2symbol,collapse = "/")
        goid2symbol[i] <- mygoid2symbol
        
      }
      return(goid2symbol)
    }
    go.result<- go.result@result
    go.result$Symbol <- NA
    go.result$Symbol <- MyNCBI2Symbol(go.result$geneID,genelist)
    
    MyWriteTable(go.result,
                 go.out)
    return(go.result)
  }

MylogCalpcaSi <- 
  function(seu.tmp,
           pc.use=10,
           group2='Glut2H_2',
           graph = F){
    library(cluster)
    seu.tmp$tmp.cluster <- 1
    seu.tmp@meta.data[seu.tmp$beta==group2,'tmp.cluster'] <- 2 
    
    data.cor <- cor(t(seu.tmp@reductions$pca@cell.embeddings[,1:pc.use]),
                    method = "s")
    data.dist <- as.dist(abs(1 - data.cor))
    si <- silhouette(seu.tmp$tmp.cluster,data.dist)
    if(graph){
      plot(si, nmax = 80, cex.names = 0.5)
    }
    return(mean(si[,3]))
  }




MyLM.coe.parallel <- 
  function(tpm,
           variable,
           node=10){
    cl <- makeCluster(node)
    calLM.coe=function(tpm,variable){
      LM=lm(tpm~variable)
      if(nrow(summary(LM)$coefficients)==2){
        return(summary(LM)$coefficients[2,1])
      }
      else{return(NA)}
    }
    calLM.both=function(tpm,variable){
      LM=lm(tpm~variable)
      if(nrow(summary(LM)$coefficients)==2){
        return(summary(LM)$coefficients[2,4])
      }
      else{return(NA)}
    }
    
    pval=parApply(cl,X = tpm,MARGIN = 1,calLM.both,variable=variable)
    names(pval)=rownames(tpm)
    coe=parApply(cl,X = tpm,MARGIN = 1,calLM.coe,variable=variable)
    names(coe)=rownames(tpm)
    LM.res <- data.frame('gene' = rownames(tpm),
                         'pval' = pval,
                         'r'= coe
    )
    LM.res$BH.adjust <- p.adjust(LM.res$pval,'BH')
    return(LM.res)
  }


Mypairwiseseuratmarker <- function(src,type.list,gene.include=genes.keep,assay = 'rawTP01M',slot = 'data',pairtype='pairwise',...){
  KO.DEG.list <- list()
  if(pairtype=='pairwise'){
    for(type1 in type.list){
      print(type1)
      type.list <- type.list[-1]
      for(type2 in type.list){
        print(type2)
        KO.DEG.list[[paste(type1,type2,sep = '_')]] <- Myseufindmarker(src,gene.include = gene.include,ident.1 = type1,assay = assay,slot = slot,
                                                                       ident.2 = type2,c1 = type1,c2 = type2)
        KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster <- factor(KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster,levels = c(type1,type2))
        
      }
    }
  }
  if(pairtype=='adjacent'){
    for(k in 1:(length(type.list)-1)){
      type1 <- type.list[k];print(type1)
      type2 <- type.list[k+1];print(type2)
      KO.DEG.list[[paste(type1,type2,sep = '_')]] <- Myseufindmarker(seu.hepa.mnncorrect,gene.include = genes.keep,ident.1 = type1,assay = assay,slot = slot,
                                                                     ident.2 = type2,c1 = type1,c2 = type2)
      KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster <- factor(KO.DEG.list[[paste(type1,type2,sep = '_')]]$cluster,levels = c(type1,type2))
      k <- k+1
      
    }
  }
  return(KO.DEG.list)
}
#######from ypl#####
#--------
Myplsda <- function(src,tmp.group,gene.input,n.comp=3,plot=T){
  library(mixOmics)
  plsda.datatm.mg3ep = plsda(t(src@assays$RNA@data[gene.input,]),src@meta.data[,tmp.group], ncomp = n.comp)
  background <- background.predict(plsda.datatm.mg3ep, comp.predicted=2 ,dist = "max.dist") 
  #plotVar(plsda.datatm)
  if(plot==T){
    plotIndiv(plsda.datatm.mg3ep, ind.names = FALSE, 
              legend=TRUE,ellipse = TRUE, title="sPLS-DA - final result")
    plotIndiv(plsda.datatm.mg3ep, comp = 1:2, 
              ind.names = FALSE, title = "Maximum distance",
              legend = TRUE,  background = background, ellipse = TRUE)
    auc.plsda.mg3ep <- auroc(plsda.datatm.mg3ep)
  }
  src@reductions$plsda =
    src@reductions$pca
  src@reductions$plsda@cell.embeddings = 
    plsda.datatm.mg3ep $variates$X
  src@reductions$plsda@feature.loadings = 
    plsda.datatm.mg3ep $loadings$X
  colnames(src@reductions$plsda@cell.embeddings) = 
    c("Comp_1","Comp_2","Comp_3")
  src@reductions$plsda@key = "Comp_"
  return(src)
}

#####from WX###
MyManifold2Gg <- 
  function(manifold,
           sample.inf,
           x.dim = 1,
           y.dim = 2,
           color_manual = brewer.pal(9,"Set1")){
    library(ggplot2)
    print("sample.inf should have the same order to manifold matrix")
    manifold.df <- cbind(x.pos = manifold[,x.dim],
                         y.pos = manifold[,y.dim],
                         sample.inf)
    p <- ggplot(data = manifold.df, 
                mapping = aes(x = x.pos, 
                              y = y.pos,
                              label = SampleName))
    p <- p + xlab(paste("Component.",x.dim,sep = ""))
    p <- p + ylab(paste("Component.",y.dim,sep = ""))
    p <- p + guides(colour=guide_legend(title=NULL))
    p <- p + scale_color_manual(values = color_manual)
    
    p <- p + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            aspect.ratio=1)
    return(p)
  }

Myproject3Dto2D <- 
  function(coord_3d,par3d){
    coord_3d_rotated=par3d$modelMatrix%*%rbind(t(coord_3d),rep(1,nrow(coord_3d)))
    
    coord_3d_stretched=par3d$projMatrix%*%coord_3d_rotated
    return(t(coord_3d_stretched[1:2,]))
  }

MyrunMICparallel <- 
  function(tpm,
           variable,node=10)
{
  cl <- makeCluster(node)
  MyMIC <- 
  function(tpm,
           variable, node=10){
    library(minerva)
    MIC=mine(variable,as.numeric(tpm),n.cores=node)$MIC
    names(MIC)=rownames(tpm)
    return(MIC)
  }
MICres=parApply(cl,X = tpm,MARGIN = 1,MyMIC,variable=variable)
names(MICres)=rownames(tpm)
MICres=MICres[!is.na(MICres)]
MICres=sort(MICres,decreasing = T)
return(MICres)
}
########
MyPseudotimevio3 <- 
  function(si.tab,plsr.zonation,Type,time.colors  = time.colors,exp.col=pc12.heatmap.col,rep = F,size.point=4,box.lwd=2,box.width=1,shape.type=c(19,25),...){
    if(rep == F){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type), 
                           colour = Type#,
                           #shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = shape.type) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(aes(x = Type, color = plsr.zonation),groupOnX = TRUE,
                         bandwidth = 1,
                         #dodge.width = 0.6,
                         alpha = 1,
                         size = size.point,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = T#,
                         #method = "pseudorandom"
        ) +scale_color_gradientn(colours = exp.col)+
        geom_violin(aes(fill = Type),
                    alpha = 0,
                    #position = dodge,
                    color = "grey30",
                    lwd=box.lwd,
                    # outlier.shape = NA,
                    show.legend = F,
                    width = box.width
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    } 
    else if(rep == T){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type_rep), 
                           colour = Type_rep,
                           shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = c(17,19)) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(groupOnX = TRUE,
                         bandwidth = 1,
                         dodge.width = 0.6,
                         alpha = 0.8,
                         size = 4,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = 10#,
                         #method = "pseudorandom"
        ) +
        geom_boxplot(aes(fill = Type_rep),
                     alpha = 0.3,
                     #position = dodge,
                     color = "grey30",
                     #outlier.shape = NA,
                     show.legend = F,
                     lwd = box.lwd
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    }
    return(p)
  }

MyPseudotimevio2 <- 
  function(si.tab,plsr.zonation,Type,time.colors  = time.colors,zonation.col=zonation.col, rep = F,size.point=4,box.lwd=2,box.width=1,dodge.width=0.6,shape.type=c(19,25),...){
    if(rep == F){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type), 
                           colour = Type#,
                           #shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = shape.type) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(aes(x = Type, color = plsr.mergerep.zonation),groupOnX = TRUE,
                         bandwidth = 1,
                         dodge.width = dodge.width,
                         alpha = 0.8,
                         size = size.point,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = T#,
                         #method = "pseudorandom"
        ) +
        geom_violin(aes(fill = Type),
                    alpha = 0,
                    #position = dodge,
                    color = "grey30",
                    lwd=box.lwd,
                    # outlier.shape = NA,
                    show.legend = F,
                    width = box.width
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    } 
    else if(rep == T){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type_rep), 
                           colour = Type_rep,
                           shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = c(17,19)) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(groupOnX = TRUE,
                         bandwidth = 1,
                         dodge.width = 0.6,
                         alpha = 0.8,
                         size = 4,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = 10#,
                         #method = "pseudorandom"
        ) +
        geom_boxplot(aes(fill = Type_rep),
                     alpha = 0.3,
                     #position = dodge,
                     color = "grey30",
                     #outlier.shape = NA,
                     show.legend = F,
                     lwd = box.lwd
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    }
    return(p)
  }

MyPseudotimevio <- 
  function(si.tab,plsr.zonation,Type,time.colors  = time.colors,zonation.col=zonation.col, rep = F,size.point=4,box.lwd=2,box.width=1,shape.type=c(19,25),...){
    if(rep == F){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type), 
                           colour = Type#,
                           #shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = shape.type) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(#aes(x = Type, color = plsr.mergerep.zonation),
                         groupOnX = TRUE,
                         bandwidth = 1,
                         #dodge.width = 0.6,
                         alpha = 0.8,
                         size = size.point,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = T#,
                         #method = "pseudorandom"
        ) +
        geom_violin(aes(fill = Type),
                    alpha = 0,
                    #position = dodge,
                    color = "grey30",
                    lwd=box.lwd,
                    # outlier.shape = NA,
                    show.legend = F,
                    width = box.width
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    } 
    else if(rep == T){
      p.time <- ggplot(data = si.tab,
                       aes(y=plsr.zonation,
                           x=as.numeric(Type_rep), 
                           colour = Type_rep,
                           shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = c(17,19)) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -plsr.zonation,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(groupOnX = TRUE,
                         bandwidth = 1,
                         dodge.width = 0.6,
                         alpha = 0.8,
                         size = 4,
                         #show.legend = T,
                         stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = 10#,
                         #method = "pseudorandom"
        ) +
        geom_boxplot(aes(fill = Type_rep),
                     alpha = 0.3,
                     #position = dodge,
                     color = "grey30",
                     #outlier.shape = NA,
                     show.legend = F,
                     lwd = box.lwd
        )+
        scale_fill_manual(values = time.colors)+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    }
    return(p)
  }
MylimmafcDataplot <- function(bed,
                              up.name,down.name,
                              phase1,
                              phase2,
                              fc.state,
                              cutoff.K,
                              xlim,
                              ylim
){
  mcols(bed)[,fc.state] <- "-"
  mcols(bed)[,fc.state][names(bed) %in% up.name &
                          mcols(bed)[,phase2] > cutoff.K] <- "up"
  mcols(bed)[,fc.state][names(bed) %in% down.name &
                          mcols(bed)[,phase1] > cutoff.K] <- "down"
  mcols(bed)[,fc.state] <- factor(mcols(bed)[,fc.state],
                                  levels = c("-",
                                             "down",
                                             "up"))
  plot(log2(mcols(bed)[,phase1]),
       log2(mcols(bed)[,phase2]),
       xlim = xlim,
       ylim = ylim,
       col = color.fc[mcols(bed)[,fc.state]],
       pch = 16,
       cex = 0.5)
  box(lwd = 8)
  MyText(paste(fc.state,"\n","up:",sum(mcols(bed)[,fc.state] == "up"),"\n",
               "down:",sum(mcols(bed)[,fc.state] == "down"),"\n",
               "-:",sum(mcols(bed)[,fc.state] == "-"),"\n"))
}
Myseucompentropy <- function(object,feature=rownames(object@assays$RNA@data)) 
{
  probs <- t(t(exp(object@assays$RNA@data[feature,])-1)/apply(exp(object@assays$RNA@data[feature,])-1, 2, sum))
  object$entropy <- -apply(probs * log(probs + 1e-10)/log(nrow(object@assays$RNA@data[feature,])), 
                           2, sum)
  return(object)
}


Mycoord <- function(coord,dims=c('PC_1','PC_2'),class=NULL,classorder=NULL){
  coord.list=list()
  for (i in 1:length(classorder)) {
    coord.list[[i]]=coord[which(coord[,class]==classorder[i]),dims]
  }
  gravity_center.list=lapply(coord.list,colMeans)
  gravity_center=do.call("rbind",gravity_center.list)
  rownames(gravity_center) <- classorder
  return(gravity_center)
}


MyExtend <- 
function(x, upstream = 0, downstream = 0, rm.neg = T) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  if(rm.neg){
    new_start[new_start<1]=1
  }
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

myplotDeviationsTsne <- function(object, tsne, var_df = NULL, sample_column = NULL, 
                                 annotation_name = NULL, shiny = interactive()) 
{
  if (is.list(tsne) && "Y" %in% names(tsne)) {
    tsne <- tsne$Y
  }
  if (nrow(tsne) != ncol(object)) {
    stop("Number of rows of tsne do not match number of columns of object. ", 
         " plotDeviationsTsne takes result of deviationsTsne for samples")
  }
  if (shiny) 
    return(plot_deviations_tsne_shiny(object, tsne, var_df, 
                                      sample_column))
  stopifnot(sample_column %in% colnames(colData(object)))
  anno <- colData(object)[, sample_column]
  out <- list()
  if (!is.null(sample_column)) {
    for (i in sample_column) {
      anno <- colData(object)[, i]
      out[[i]] <- ggplot(data.frame(x = tsne[, 1], y = tsne[, 
                                                            2], color = anno, text = colnames(object)), 
                         aes_string(x = "x", y = "y", col = "color", 
                                    text = "text")) + geom_point(size = 2) + chromVAR_theme() + 
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid::unit(0.5, 
                                                                                     "lines"))
      if (nlevels(as.factor(anno)) <= 8) {
        out[[i]] <- out[[i]] + scale_color_brewer(palette = "Dark2", 
                                                  name = i)
      }
      else {
        out[[i]] <- out[[i]] + guides(colour = guide_legend(title = i))
      }
    }
  }
  if (!is.null(annotation_name)) {
    for (i in annotation_name) {
      if (i %in% rownames(object)) {
        ix <- match(c(i), rownames(object))
      }
      else if (i %in% rowData(object)$name) {
        ix <- which(rowData(object)$name == i)
        if (length(ix) > 1) 
          ix <- ix[which.max(row_sds(assays(object[ix, 
                                                   ])$z))]
      }
      else if (is.numeric(i)) {
        ix <- i
      }
      else {
        stop("annotation_name invalid")
      }
      if ("name" %in% colnames(rowData(object))) {
        name_val <- rowData(object)[ix, "name"]
      }
      else {
        name_val <- i
      }
      out[[i]] <- ggplot(data.frame(x = tsne[, 1], y = tsne[, 
                                                            2], color = deviationScores(object)[ix, ], text = colnames(object)), 
                         aes_string(x = "x", y = "y", col = "color", 
                                    text = "text")) + geom_point(size = 2) + scale_color_gradient2(name = i, 
                                                                                                   mid = "lightgray", low = "blue", high = "red") + 
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
        theme(legend.key.size = grid::unit(0.5, "lines"))
    }
  }
  return(out)
}
MyCalUcsccomquanNormFactor <- 
  function(input.bed,
           list.prefix,
           quantile.cutoff = 0.5){
    merge.bed <- input.bed
    for (sn.tmp in list.prefix) {
      input.bed =
        input.bed[mcols(input.bed)[,sn.tmp] > quantile(mcols(merge.bed)[,sn.tmp],quantile.cutoff)]
      
    }
    cat("filtered peak count:",
        length(input.bed),
        "\n")
    norm.factor = apply(mcols(input.bed)[,list.prefix], 2, mean)
    norm.factor = norm.factor / mean(norm.factor)
    return(norm.factor)
}


Mysfnorm <- function(raw.data,
                     norm.method = "SizeFactor",
                     ercc.row = 40825:40916,
                     graph = TRUE,
                     min.cutoff = 1,
                     ...){
  raw.data <- t(raw.data)
  # plot log2 boxplot, remove 
  LogBoxPlot <- function(data,
                         min.cutoff = min.cutoff,
                         ...){
    data.na <- data
    data.na[data.na < min.cutoff] <- NA
    boxplot(log2(data.na),
            outline = F,
            ...)
  }
  
  # normalize by size factor
  SizeFactorNorm <- function(raw.data) {
    library("DESeq")
    size.factor <- estimateSizeFactorsForMatrix(raw.data)
    norm.data <-  t(t(raw.data) / size.factor)
    return(norm.data)
  }
  
  # normalize by ERCC size factor
  ERCCSizeFactorNorm <- function(raw.data,
                                 ercc.row) {
    library("DESeq")
    size.factor <- estimateSizeFactorsForMatrix(raw.data[ercc.row,])
    norm.data <-  t(t(raw.data) / size.factor)
    return(norm.data)
  }
  
  if(graph){
    LogBoxPlot(raw.data,
               min.cutoff = min.cutoff)
  }
  
  if (norm.method == "SizeFactor") {
    norm.data <- SizeFactorNorm(raw.data)
  }else if (norm.method == "ERCCSizeFactor") {
    norm.data <- ERCCSizeFactorNorm(raw.data,
                                    ercc.row)
  }else if (norm.method == "Raw") {
    norm.data <- raw.data
  }
  
  if(graph){
    LogBoxPlot(norm.data,
               min.cutoff = min.cutoff)
  }
  return(norm.data)
}

######
Myvenn <- function(gene.venn,col.list=time.colors,lwd=5,...){
  gene.venn.them <- compute.Venn(gene.venn)
  gene.venn.them.p <- VennThemes(gene.venn.them)
  
  gene.venn.them.p$Set$Set1$col <- col.list[1]
  gene.venn.them.p$Set$Set2$col <- col.list[2]
  
  gene.venn.them.p$Set$Set1$lwd <- 5
  gene.venn.them.p$Set$Set2$lwd <- 5
  
  gene.venn.them.p$SetText$Set1$fontsize <- 0.01
  gene.venn.them.p$SetText$Set2$fontsize <- 0.01
  
  
  gene.venn.them.p$FaceText$`11`$fontsize <- 0
  gene.venn.them.p$FaceText$`10`$fontsize <- 0
  gene.venn.them.p$FaceText$`01`$fontsize <- 0
  
  gene.venn.them.p$Face$`11`$fill <- NULL
  gene.venn.them.p$Face$`10`$fill <- NULL
  gene.venn.them.p$Face$`01`$fill <- NULL
  return(gene.venn.them.p)
}

##########qiu#########
MyCalUcscNormFactor2 = function(list.peak.tmp,
                               list.bdg.tmp,
                               quantile.cutoff = 0.5){
  nor.list <- list()
  bed.merge.tmp = list.peak.tmp[[1]]
  for (sn.tmp in names(list.peak.tmp)[-1]) {
    bed.merge.tmp <- reduce(c(bed.merge.tmp,
                              list.peak.tmp[[sn.tmp]]))
  }
  
  for (sn.tmp in names(list.bdg.tmp)) {
    print(sn.tmp)
    mcols(bed.merge.tmp)[,sn.tmp] = MyCalTd(bed.merge.tmp,
                                            list.bdg.tmp[[sn.tmp]],
                                            "td")
  }
  
  bed.merge.common.tmp = bed.merge.tmp
  for (sn.tmp in names(list.peak.tmp)) {
    bed.merge.common.tmp <- MyBedGrepU(bed.merge.common.tmp,
                                       list.peak.tmp[[sn.tmp]])
  }
  
  bed.merge.common.quantile.tmp = bed.merge.common.tmp
  for (sn.tmp in names(list.bdg.tmp)) {
    bed.merge.common.quantile.tmp =
      bed.merge.common.quantile.tmp[mcols(bed.merge.common.quantile.tmp)[,sn.tmp] > quantile(mcols(bed.merge.common.tmp)[,sn.tmp],quantile.cutoff)]
    
  }
  cat("filtered peak count:",
      length(bed.merge.common.quantile.tmp),
      "\n")
  norm.factor = apply(mcols(bed.merge.common.quantile.tmp)[,names(list.bdg.tmp)], 2, mean)
  norm.factor = norm.factor / mean(norm.factor)
  nor.list[['merge.bed']] <- bed.merge.tmp
  nor.list[['norm.factor']] <- norm.factor
  return(nor.list)
}
###########
MyPseudotimebox2 <- 
  function(si.tab,pseudotime,Type,time.colors  = time.colors,rep = F,heter.group.merge = 'heter.group.merge',...){
    if(rep == F){
      p.time <- ggplot(data = si.tab,
                       aes(y=pseudotime,
                           x=as.numeric(Type)
                           # colour = Type,
                           #shape = heter.group.merge#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = c(19,25,20)) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        geom_quasirandom(aes(x = Type, color = heter.group.merge),
                         groupOnX =T,
                         bandwidth = 1,
                        # dodge.width = 1,
                         alpha = 0.8,
                         size = 2#,
                         #show.legend = T,
                         # stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = 10#,
                         #method = "pseudorandom"
        )+scale_color_manual(values = time.colors)+
        geom_boxplot(aes(x = Type,fill = Type),
                     alpha = 0.2,
                     #position = dodge,
                     color = "grey30",
                     outlier.shape = NA,
                     show.legend = F#,
                    # width = 0.5
        )+
        scale_fill_manual(values = rep('white',each=15))+
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    } 
    else if(rep == T){
      p.time <- ggplot(data = si.tab,
                       aes(y=pseudotime,
                           x=as.numeric(Type_rep)#, 
                           #colour = Type_rep#,
                           #shape = proliferative#,
                           #size = 1
                       )
      ) +
        scale_shape_manual(values = c(19,25)) +
        scale_color_manual(values = time.colors) +
        theme(axis.text = element_text(size = 30, colour = "black")) +
        theme(axis.title.x = element_text(size = 40, colour = "black")) +
        theme(axis.title.y = element_text(size = 40, colour = "black")) +
        theme(legend.text = element_text(size = 40, colour = "black")) +
        theme(legend.title = element_text(size = 30, colour = "black")) +
        xlab("") +
        ylab("") +
        guides(colour=guide_legend(title=NULL)) +
        theme_bw() +
        theme(legend.key = element_blank()) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 4,
                                          colour = "black"))
      library(ggbeeswarm)
      p <- p.time +
        # scale_x_reverse()+
        # #scale_shape_manual(values = c(16,1)) +
        # geom_point(aes(x = -Pseudotime,
        #                y = -as.numeric(Age),
        #                col = Age
        #                #shape = Ngn3Exp
        # ),
        # na.rm = TRUE,
        # size = 2) +
        geom_quasirandom(aes(x = Type_rep, color = heter.group.merge),
                         groupOnX =T,
                         bandwidth = 1,
                         # dodge.width = 1,
                         alpha = 0.8,
                         size = 2#,
                         #show.legend = T,
                         # stroke = 0#,#setting stroke to 0 removes the outline around the points
                         #nbins = 30
                         #varwidth = 10#,
                         #method = "pseudorandom"
        )+scale_color_manual(values = time.colors)+
        geom_boxplot(aes(x = Type_rep,fill = Type_rep),
                     alpha = 0.2,
                     #position = dodge,
                     color = "grey30",
                     outlier.shape = NA,
                     show.legend = F
                     # width = 0.5
        )+
        scale_fill_manual(values = rep('white',each=40))+
        
        
        guides(colour = guide_legend(title = "",
                                     keywidth = 3,
                                     keyheight = 3,
                                     override.aes = list(size=12),
                                     order = 1,
                                     ncol = 1))+
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.text = element_text(size = 10, colour = "black"))
    }
    plot(p)
  }

######
MyFDL2Gg <- 
function(dif,
         sample.inf,
         x.dim = 1,
         y.dim = 2,
         color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  dif <- as.data.frame(dif)
  dif.df <- cbind(x.pos = dif[,x.dim],
                  y.pos = dif[,y.dim],
                  sample.inf[row.names(dif),])
  
  p <- ggplot(data = dif.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("FDL",x.dim,sep = ""))
  p <- p + ylab(paste("FDL",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)+scale_color_manual(values = time.colors) +
    scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2))
  return(p)
}




#######wangxin########
MyLM.parallel=function(tpm,
                       variable,
                       direction="both",
                       node=10){
  cl <- makeCluster(node)
  calLM.pos=function(tpm,variable){
    LM=lm(tpm~variable)
    if(nrow(summary(LM)$coefficients)==2){
      if(summary(LM)$coefficients[2,1]>0){
        return(summary(LM)$coefficients[2,4])
      }else{return(NA)}
    }
    else{return(NA)}
  }
  calLM.neg=function(tpm,variable){
    LM=lm(tpm~variable)
    if(nrow(summary(LM)$coefficients)==2){
      if(summary(LM)$coefficients[2,1]<0){
        return(summary(LM)$coefficients[2,4])
      }else{return(NA)}
    }
    else{return(NA)}
  }
  calLM.both=function(tpm,variable){
    LM=lm(tpm~variable)
    if(nrow(summary(LM)$coefficients)==2){
      return(summary(LM)$coefficients[2,4])
    }
    else{return(NA)}
  }
  if (direction=="pos") {
    pval=parApply(cl,X = tpm,MARGIN = 1,calLM.pos,variable=variable)
  } else if (direction=="neg") {
    pval=parApply(cl,X = tpm,MARGIN = 1,calLM.neg,variable=variable)
  } else if (direction=="both") {
    pval=parApply(cl,X = tpm,MARGIN = 1,calLM.both,variable=variable)
  }
  names(pval)=rownames(tpm)
  pval=pval[!is.na(pval)]
  pval=sort(pval)
  return(pval)
}


Mybarplot <- function(si.tab,c1,c2,xlim = 20,cols='black',...){
  type.distri <- table(si.tab[,c1],si.tab[,c2])
  type.ratio <- apply(type.distri,1,prop.table)
  type.ratio <- reshape2::melt(type.ratio)
  
  type.ratio <- apply(type.distri,1,prop.table)*100
  
  barplot(type.ratio,
          col = add.alpha(cols,alpha = 0.8),
          cex.axis = 2,
          xlim = c(0,xlim),  
          ylim = c(0,100),
          #xaxt = "n",
          yaxt = "n",
          ...
  )
  #axis(1,lwd = 5,lwd.ticks = 0, cex.axis = 2,labels = NA)
  axis(2,lwd = 5,lwd.ticks = 3, cex.axis = 2)
  #axis(1,at = 1:6,lwd = 5,lwd.ticks = 3, cex.axis = 2,labels = names(type.ratio2))
  return(type.ratio)
}


MyGozebrafish <- 
  function(selected.id,
           universe.id,
           #organism = "mmu",
           #keytype = 'ENSEMBL',
           pvalue=0.1,
           qvalue=0.5,
           ont="BP",
           readable=TRUE,
           plot=F,
           font=0.5,
           n=100,
           #slim = FALSE,
           ...){
    library(clusterProfiler)
    library(org.Dr.eg.db)
    OrgDb <- org.Dr.eg.db
    goresult <- enrichGO(
      gene = selected.id,
      universe=universe.id,
      OrgDb = OrgDb,
      #keytype =keytype,
      ont = ont,
      #pAdjustMethod = "BH",
      pvalueCutoff  = pvalue,
      qvalueCutoff  = qvalue, 
      readable=readable)
    if(plot=="TRUE"){
      emapplot(goresult
               #vertex.label.font = font,
               #n=n
      )
    }
    
    return(goresult)
  }

MyGOZebrafishwritetable <- 
  function(genelist,
           go.out,
           universe.gene=zebrafish.gene.inf$ENTREZID,
           pvalue = 0.1,
           qvalue = 1,
           ...){
    go.result <- MyGozebrafish(genelist,
                               universe.gene,
                               pvalue = pvalue,
                               qvalue = qvalue,
                               ...)
    MyWriteTable(go.result,
                 go.out)
    return(go.result)
  }
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  zebrafish <-  useMart("ensembl", dataset = "drerio_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human,attributesL = c("zfin_id_symbol"), martL = zebrafish, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


convertMouseGene2sheepList <- function(x){

  require("biomaRt")
 # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #zebrafish <-  useMart("ensembl", dataset = "drerio_gene_ensembl")
  sheep <-  useMart("ensembl", dataset = "oaries_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse,attributesL = c('entrezgene_id'), martL = sheep, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}



##########
#https://github.com/dsancheztaltavull/Bayesian-Correlations/blob/master/BaCo.R
#############
BaCo <- function(X){
  
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}
####alexander-7 https://github.com/satijalab/seurat/issues/1260
##R-Function as workaround for the scaling issue of the 'split.by' option of Feature Plot
##########
plot_seperate_features <- function(object, gene, limit = NULL,ident=NULL,pt.size=0.5,revers.order = F,numcol=NULL, colors =c('#DCD9DE','#870052'), ...){
  
  if(is.null(ident)){
    Idents(object) <- 'combined'
  }
  else{
    Idents(object) <- object[[ident]]
  }
  p <- vector("list", length = (length(gene)*length(levels(Idents(object)))))
  for(g in 1:(length(gene))){
    if(is.null(limit)==T){
      for(i in 1:length(levels(Idents(object)))){
        j <- i+((g-1)*length(levels(Idents(object))))
        p[[j]] <- FeaturePlot(object, cells = WhichCells(object,idents = levels(Idents(object))[i]), features = gene[g],order = T,cols = colors, pt.size = pt.size, sort.cell = T,...) +
          scale_color_gradient(low = colors[1], high = colors[2], limits = c(0,ceiling(max(object@assays$RNA@data[gene[g],]))),oob = squish) +
          ggtitle(paste(levels(Idents(object))[i],gene[g],sep=", "))
      }
    }
    else{
      for(i in 1:length(levels(Idents(object)))){
        j <- i+((g-1)*length(levels(Idents(object))))
        p[[j]] <- FeaturePlot(object, cells = WhichCells(object,idents = levels(Idents(object))[i]), features = gene[g],order = T,cols = colors, pt.size = pt.size, sort.cell = T,...) +
          scale_color_gradient(low = colors[1], high = colors[2], limits = limit,oob = squish) +
          ggtitle(paste(levels(Idents(object))[i],gene[g],sep=", "))
      }
    }
  }
  #correct axis in all plots
  y_axis_max <- 1:length(p)
  y_axis_min <- 1:length(p)
  x_axis_max <- 1:length(p)
  x_axis_min <- 1:length(p)
  for(i in 1:length(p)){
    x_axis_max[i] <- max(layer_scales(p[[i]])$x$range$range)
    x_axis_min[i] <- min(layer_scales(p[[i]])$x$range$range)
    y_axis_max[i] <- max(layer_scales(p[[i]])$y$range$range)
    y_axis_min[i] <- min(layer_scales(p[[i]])$y$range$range)
  }
  for(i in 1:length(p)){
    p[[i]] <- p[[i]]+coord_equal(xlim=c(min(x_axis_min),max(x_axis_max)),ylim=c(min(y_axis_min),max(y_axis_max)))
  }
  if(revers.order==T) p <- rev(p)
  gridExtra::grid.arrange(grobs = p,ncol=numcol)
}


Mygene2pseudotime <- function(genes.inf.input = genes.inf.input ,
                              gene.list,
                              tpm.data,
                              Pseudotime,
                              Time,
                              sample.list,
                              log = T,
                              point.size=3,
                              cols=liver.col
                              
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
  p.loop <- ggplot(data = time.data) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_color_manual(values = cols) +
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
                linewidth = 2,
                data = MyGetTrend2(time.data$Pseudotime,
                                   time.data[,ens],
                                   min.rpkm = 1))
    print(p.plot)
    
    run.count = run.count + 1
  }
}

#########function######
MySeuratv3FDL2Gg <- 
  function(FDL.coord,
           seurat.meta,
           x.dim = 1,
           y.dim = 2,
           color_manual = brewer.pal(9,"Set1")){
    library(ggplot2)
    tSNE.coord <- as.data.frame(FDL.coord)
    
    
    tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                     y.pos = tSNE.coord[,y.dim],
                     seurat.meta[row.names(tSNE.coord),])
    
    p <- ggplot(data = tSNE.df, 
                mapping = aes(x = x.pos, 
                              y = y.pos,
                              label = SampleName
                ))
    p <- p + xlab(paste("FDL",x.dim,sep = ""))
    p <- p + ylab(paste("FDL",y.dim,sep = ""))
    p <- p + guides(colour=guide_legend(title=NULL))
    p <- p + scale_color_manual(values = color_manual)
    
    p <- p + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    return(p)
  }
MySeuratPca2Plotly3d <- function(seurat, 
                                 x.dim = 1, 
                                 y.dim = 2, 
                                 z.dim = 3,
                                 ...){
  pca.coord <- Embeddings(object = seurat[['pca']])
  library(plotly) 
  p <- plot_ly(x =pca.coord[,x.dim],
               y =pca.coord[,y.dim],
               z =pca.coord[,z.dim],
               type ="scatter3d",
               mode = "markers",
               ...) 
  return(p) 
  detach("package:plotly",
         unload=TRUE) 
} 


MySeuratv3UMAP10x2Gg <- function(seurat,
                                 seurat.meta,
                                 x.dim = 1,
                                 y.dim = 2,
                                 color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- Embeddings(object = seurat[['umap']])
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   seurat.meta[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName
              ))
  p <- p + xlab(paste("UMAP",x.dim,sep = ""))
  p <- p + ylab(paste("UMAP",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_color_manual(values = time.colors) +
    #scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1,
                                 ncol = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2,
                                ncol = 1))
  return(p)
}

MySeuratv3PCA10x2Gg <- function(seurat,
                                seurat.meta,
                                x.dim = 1,
                                y.dim = 2,
                                color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- Embeddings(object = seurat[['pca']])
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   seurat.meta[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName
              ))
  p <- p + xlab(paste("PC",x.dim,sep = ""))
  p <- p + ylab(paste("PC",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}


MySeuratv3TSNE10x2Gg <- function(seurat,
                                 seurat.meta,
                                 x.dim = 1,
                                 y.dim = 2,
                                 color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- Embeddings(object = seurat[['tsne']])
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   seurat.meta[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName
              ))
  p <- p + xlab(paste("TSNE",x.dim,sep = ""))
  p <- p + ylab(paste("TSNE",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_color_manual(values = time.colors) +
    #scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1,
                                 ncol = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2,
                                ncol = 1))
  return(p)
}

MyfcDataplot <- function(bed,
                         phase1,
                         phase2,
                         fc.state,
                         fc.value.up,
                         fc.value.down,
                         fc.cutoff,
                         cutoff.K,
                         xlim,
                         ylim
){
  mcols(bed)[,fc.state] <- "-"
  mcols(bed)[,fc.state][mcols(bed)[,phase1] / mcols(bed)[,phase2] < (1/fc.cutoff) &
                          mcols(bed)[,phase2] > cutoff.K] <- "up"
  mcols(bed)[,fc.state][mcols(bed)[,phase2] / mcols(bed)[,phase1] < (1/fc.cutoff) &
                          mcols(bed)[,phase1] > cutoff.K] <- "down"
  mcols(bed)[,fc.state] <- factor(mcols(bed)[,fc.state],
                                  levels = c("-",
                                             "down",
                                             "up"))
  mcols(bed)[,fc.value.up] <- ifelse(mcols(bed)[,phase2] > cutoff.K,
                                     mcols(bed)[,phase1] / mcols(bed)[,phase2],
                                     "-")
  mcols(bed)[,fc.value.down] <- ifelse(mcols(bed)[,phase1] > cutoff.K,
                                       mcols(bed)[,phase2] / mcols(bed)[,phase1],
                                       "-")
  plot(log2(mcols(bed)[,phase1]),
       log2(mcols(bed)[,phase2]),
       xlim = xlim,
       ylim = ylim,
       col = color.fc[mcols(bed)[,fc.state]],
       pch = 16,
       cex = 0.5)
  box(lwd = 4)
  MyText(paste(fc.state,"\n","up:",sum(mcols(bed)[,fc.state] == "up"),"\n",
               "down:",sum(mcols(bed)[,fc.state] == "down"),"\n",
               "-:",sum(mcols(bed)[,fc.state] == "-"),"\n"))
  return(bed)
}
Myboxplot <- function(tpm.input,
                      phase,
                      label,
                      main_text='',
                      side="two.side",
                      color=box.col,
                      notch = T,
                      paired = T,
                      cex.axis = 1.5,
                      #ylim,
                      ...
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  #phase3 <- phase[3]
  boxplot(tpm.input[,phase],
          col=color,
          #ylim=ylim,
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n",
          cex.axis = 1.3,
          ...
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=3)
  wilcox.test.res1 <- wilcox.test(tpm.input[,phase1],
                                  tpm.input[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  # wilcox.test.res2 <- wilcox.test(tpm.input[,phase2],
  #                                 tpm.input[,phase3],
  #                                 alternative = side,
  #                                 paired = paired)
  
  # MyText(paste("p-value: \n",
  #              wilcox.test.res2$p.value,"\n"),
  #        text.cex = 1)
  
}
Myboxplot3 <- function(tpm.input,
                       phase,
                       label,
                       main_text='',
                       side="two.side",
                       color=box.col,
                       notch = T,
                       paired = T,
                       cex.axis = 1.5,
                       #ylim,
                       ...
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  phase3 <- phase[3]
  boxplot(tpm.input[,phase],
          col=color,
          #ylim=ylim,
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n",
          cex.axis = 1.3
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=3)
  wilcox.test.res1 <- wilcox.test(tpm.input[,phase1],
                                  tpm.input[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res2 <- wilcox.test(tpm.input[,phase2],
                                  tpm.input[,phase3],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res2$p.value,"\n"),
         text.cex = 1)
  
}

Myboxplot4 <- function(tpm.input,
                       phase,
                       label,
                       main_text='',
                       side="two.side",
                       color=box.col,
                       notch = T,
                       paired = T,
                       cex.axis = 1.5,
                       #ylim,
                       ...
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  phase3 <- phase[3]
  phase4 <- phase[4]
  boxplot(tpm.input[,phase],
          col=color,
          #ylim=ylim,
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n",
          cex.axis = 1.3
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=3)
  wilcox.test.res1 <- wilcox.test(tpm.input[,phase1],
                                  tpm.input[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res2 <- wilcox.test(tpm.input[,phase2],
                                  tpm.input[,phase3],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res2$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res3 <- wilcox.test(tpm.input[,phase3],
                                  tpm.input[,phase4],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res3$p.value,"\n"),
         text.cex = 1)
  
}


Myboxplot3 <- function(tpm.input,
                       phase,
                       label,
                       main_text='',
                       side="two.side",
                       color=box.col,
                       notch = T,
                       paired = T,
                       cex.axis = 1.5,
                       #ylim,
                       ...
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  phase3 <- phase[3]
  #phase4 <- phase[4]
  boxplot(tpm.input[,phase],
          col=color,
          #ylim=ylim,
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n",
          cex.axis = 1.3
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=3)
  wilcox.test.res1 <- wilcox.test(tpm.input[,phase1],
                                  tpm.input[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res2 <- wilcox.test(tpm.input[,phase2],
                                  tpm.input[,phase3],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res2$p.value,"\n"),
         text.cex = 1)
  
  # wilcox.test.res3 <- wilcox.test(tpm.input[,phase3],
  #                                 tpm.input[,phase4],
  #                                 alternative = side,
  #                                 paired = paired)
  # 
  # MyText(paste("p-value: \n",
  #              wilcox.test.res3$p.value,"\n"),
  #        text.cex = 1)
  
}



Mybed2boxplot <- function(bed,
                          phase,
                          label,
                          main_text='',
                          side="two.side",
                          color=box.col2,
                          notch = T,
                          paired = T,
                          cex.axis = 1.3
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  phase3 <- phase[3]
  boxplot(as.matrix(mcols(bed)[,phase]),
          col=color,
          #ylim=c(0,60),
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n"
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=2)
  wilcox.test.res1 <- wilcox.test(mcols(bed)[,phase1],
                                  mcols(bed)[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res2 <- wilcox.test(mcols(bed)[,phase2],
                                  mcols(bed)[,phase3],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res2$p.value,"\n"),
         text.cex = 1)
  
}


MyScoreMatrixBin2 <- function(windows,
                              ip.target,
                              #input.target,
                              windows.resize = 4000,
                              bin.num = 100,
                              weight.col = "TPM",
                              score.min = NULL,
                              score.max = NULL,
                              is.log = TRUE,
                              ...){
  
  MyAddExtendChr <- function(windows.bed,
                             target.bed,
                             weight.col){
    windows.chr.len <- elementNROWS(coverage(windows.bed))
    target.chr.len <- elementNROWS(coverage(target.bed))
    
    
    add.chr <- names(windows.chr.len)[!names(windows.chr.len) %in% names(target.chr.len)]
    if(length(add.chr) > 0){
      add.chr.bed <- GRanges(seqnames = Rle(add.chr),
                             ranges = IRanges(start = 1,
                                              end = windows.chr.len[add.chr]))
    }else{
      add.chr.bed <- GRanges()
    }
    
    common.chr <- names(windows.chr.len)[names(windows.chr.len) %in% names(target.chr.len)]
    extend.chr <- common.chr[windows.chr.len[common.chr] > target.chr.len[common.chr]]
    if(length(extend.chr) > 0){
      extend.chr.bed <- GRanges(seqnames = Rle(extend.chr),
                                ranges = IRanges(start = target.chr.len[extend.chr] + 1,
                                                 end = windows.chr.len[extend.chr]))
    }else{
      extend.chr.bed <- GRanges()
    }
    add.extend.chr.bed <- c(add.chr.bed,
                            extend.chr.bed)
    mcols(add.extend.chr.bed)[,weight.col] <- 0
    
    return(add.extend.chr.bed)
  }
  
  if(!is.null(windows.resize)){
    windows = resize(windows,
                     windows.resize,
                     fix = "center")
  }
  
  ip.target.add.extend <- MyAddExtendChr(windows,
                                         ip.target,
                                         weight.col)
  
  sm.ip = ScoreMatrixBin(target = c(ip.target,
                                    ip.target.add.extend),
                         windows = windows,
                         bin.num = bin.num,
                         weight.col = weight.col,
                         ...)
  #input.target.add.extend <- MyAddExtendChr(windows,
  #  input.target,
  # weight.col)
  #sm.input = ScoreMatrixBin(target = c(input.target,
  #   #        input.target.add.extend),
  # windows = windows,
  # bin.num = bin.num,
  # weight.col = weight.col,
  #  ...)
  
  sm <- sm.ip
  # sm@.Data <- sm.ip@.Data - sm.input@.Data
  if(!is.null(score.min)){
    sm@.Data[sm@.Data < score.min] <- score.min
  }
  if(!is.null(score.max)){
    sm@.Data[sm@.Data > score.max] <- score.max
  }
  if(is.log){
    sm <- log2(sm)
  }
  return(sm)
}



Myfcvalue <- function(bed,
                      phase1,
                      phase2,
                      fc.value.up,
                      fc.value.down){
  mcols(bed)[,fc.value.up] <- ifelse(mcols(bed)[,phase2] > 0,
                                     mcols(bed)[,phase1] / mcols(bed)[,phase2],
                                     "-")
  mcols(bed)[,fc.value.down] <- ifelse(mcols(bed)[,phase1] > 0,
                                       mcols(bed)[,phase2] / mcols(bed)[,phase1],
                                       "-")
  return(bed)
}
Mycor2hexbin <- function(df,
                         x,
                         y,
                         xname = "H3K27ac"){
  
  df <- data.frame(x = df[,x],
                   y = df[,y])
  p <- ggplot(df, 
              aes(log2(x+0.1),
                  log2(y+0.1)))+
    theme_minimal()+
    xlab(paste("log2(",xname,"+0.1)"))+
    ylab("log2(ATAC-seq+0.1)")+
    geom_hex(bins=50)+
    scale_fill_gradientn(colours = color.tpm)+
    geom_smooth(color = "black",method = "lm", se = FALSE)
  
  p.sub <- p+ 
    theme(panel.border = element_rect(size = 2,
                                      colour = "black",
                                      fill=NA))+
    theme(axis.text.y = element_text(size = 15, colour = "black")) +
    theme(axis.text.x= element_text(size = 20, colour = "black"))+
    theme(axis.title.x = element_text(size = 20, colour = "black")) +
    theme(axis.title.y = element_text(size = 20, colour = "black")) +
    theme(legend.text = element_text(size = 20, colour = "black")) +
    theme(legend.title = element_text(size = 20, colour = "black"))
}
Mycor2hexbinlfc <- function(df,
                            x,
                            y,
                            xLab = "log2(fc(ATAC))",
                            yLab = "log2(fc(H3K4me1))",
                            mAin = "P9-18"){
  
  df <- data.frame(x = df[,x],
                   y = df[,y])
  p <- ggplot(df, 
              aes(log2(x+0.1),
                  log2(y+0.1)))+
    theme_minimal()+
    xlab(xLab)+
    ylab(yLab)+
    ggtitle(mAin)+
    geom_hex(bins=50)+
    scale_fill_gradientn(colours = color.tpm)+
    geom_smooth(color = "black",method = "lm", se = FALSE)
  
  p.sub <- p+ 
    theme(panel.border = element_rect(size = 2,
                                      colour = "black",
                                      fill=NA))+
    theme(axis.text.y = element_text(size = 15, colour = "black")) +
    theme(axis.text.x= element_text(size = 20, colour = "black"))+
    theme(axis.title.x = element_text(size = 20, colour = "black")) +
    theme(axis.title.y = element_text(size = 20, colour = "black")) +
    theme(legend.text = element_text(size = 20, colour = "black")) +
    theme(legend.title = element_text(size = 20, colour = "black"))+
    theme(plot.title = element_text(size=20, hjust = 0.5,colour = "black"))
}
Mycor2hexbinlfcfilter <- function(df,
                                  x,
                                  y,
                                  filter = 10,
                                  xLab = "log2(fc(ATAC))",
                                  yLab = "log2(fc(H3K4me1))",
                                  mAin = "P9-18"){
  
  df <- data.frame(x = df[,x],
                   y = df[,y])
  df <- subset(df,log2(y+0.1)<10)
  df <- subset(df,log2(x+0.1)<10)
  
  p <- ggplot(df, 
              aes(log2(x+0.1),
                  log2(y+0.1)))+
    theme_minimal()+
    xlab(xLab)+
    ylab(yLab)+
    ggtitle(mAin)+
    geom_hex(bins=50)+
    scale_fill_gradientn(colours = color.tpm)+
    geom_smooth(color = "black",method = "lm", se = FALSE)
  
  p.sub <- p+ 
    theme(panel.border = element_rect(size = 2,
                                      colour = "black",
                                      fill=NA))+
    theme(axis.text.y = element_text(size = 15, colour = "black")) +
    theme(axis.text.x= element_text(size = 20, colour = "black"))+
    theme(axis.title.x = element_text(size = 20, colour = "black")) +
    theme(axis.title.y = element_text(size = 20, colour = "black")) +
    theme(legend.text = element_text(size = 20, colour = "black")) +
    theme(legend.title = element_text(size = 20, colour = "black"))+
    theme(plot.title = element_text(size=20, hjust = 0.5,colour = "black"))
}
Myfcvalue2 <- function(df,
                       phase1,
                       phase2,
                       fc.value){
  df[,fc.value] <- (df[,phase2]+0.1)/(df[,phase1]+0.1)
  return(df)
}

myHeatmap2col <- function(x,
                          #ord,
                          xlab="",
                          ylab="",
                          main="My Heatmap",
                          col=atac.state.col,
                          cutoff.min=-1,
                          cutoff.max=1,
                          
                          ...){
  op <- par(mar=c(5,2,2,2)+0.1)
  on.exit(par(op))
  nc <- NCOL(x)
  nr <- NROW(x)
  #labCol <- colnames(x)
  x[,1] <- x[,1][order(x[,1])]
  x[,2] <- x[,2][order(x[,2])]
  #x[,3] <- x[,3][order(x[,3])]
  x <- t(x)
  x[x< cutoff.min] <- cutoff.min
  x[x > cutoff.max] <- cutoff.max
  image(1L:nc, 1L:nr, 
        x, 
        xlim = 0.5 + c(0, nc),
        ylim = 0.5 +c(0, nr),
        axes = FALSE, 
        xlab=xlab,
        ylab=ylab,
        main=main,
        col=col,
        zlim = c(cutoff.min,cutoff.max),
        ...
  )
  
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0)
  axis(2, 1L:nr, labels = NA, las = 2, line = -0.5, tick = 0)
}
MyatacfcDataplot <- function(bed,
                             phase1,
                             phase2,
                             fc.state,
                             fc.cutoff,
                             cutoff.K,
                             xlim = c(-8,6),
                             ylim = c(-8,6),
                             cex = 0.5,
                             cex.axis = 1
){
  mcols(bed)[,fc.state] <- "-"
  mcols(bed)[,fc.state][mcols(bed)[,phase1] / mcols(bed)[,phase2] < (1/fc.cutoff) &
                          mcols(bed)[,phase2] > cutoff.K] <- "up"
  mcols(bed)[,fc.state][mcols(bed)[,phase2] / mcols(bed)[,phase1] < (1/fc.cutoff) &
                          mcols(bed)[,phase1] > cutoff.K ] <- "down"
  mcols(bed)[,fc.state] <- factor(mcols(bed)[,fc.state],
                                  levels = c("-",
                                             "down",
                                             "up"))
  plot(log2(mcols(bed)[,phase1]),
       log2(mcols(bed)[,phase2]),
       xlim = xlim,
       ylim = ylim,
       col = color.fc[mcols(bed)[,fc.state]],
       pch = 16,
       cex = cex,
       cex.axis= cex.axis)
  box(lwd = 8)
  MyText(paste(fc.state,"\n","up:",sum(mcols(bed)[,fc.state] == "up"),"\n",
               "down:",sum(mcols(bed)[,fc.state] == "down"),"\n",
               "-:",sum(mcols(bed)[,fc.state] == "-"),"\n"))
  return(bed)
}
Mybed2df <- function(bed
){
  df <- data.frame("atac.P9" = bed$atac.P9,
                   "atac.P18" = bed$atac.P18,
                   "atac.P60" = bed$atac.P60,
                   "k4me1.P9" = bed$H3K4me1.P9,
                   "k4me1.P18" = bed$H3K4me1.P18,
                   "k4me1.P60" = bed$H3K4me1.P60,
                   "k27ac.P9" = bed$H3K27ac.P9,
                   "k27ac.P18" = bed$H3K27ac.P18,
                   "k27ac.P60" = bed$H3K27ac.P60
  )
  rownames(df) <- names(bed)
  return(df)
}
######MystateDataplot###########
MystateDataplot <- function(bed,
                            k4me1,
                            k27ac,
                            state,
                            cutoff.k4,
                            cutoff.k27,
                            xlim,
                            ylim){
  mcols(bed)[,state] <- NA
  mcols(bed)[,state][mcols(bed)[,k4me1] < cutoff.k4] <- "unmarked"
  mcols(bed)[,state][mcols(bed)[,k4me1] >= cutoff.k4 &
                       mcols(bed)[,k27ac] < cutoff.k27] <- "primed"
  mcols(bed)[,state][mcols(bed)[,k4me1] >= cutoff.k4 &
                       mcols(bed)[,k27ac] >= cutoff.k27] <- "active"
  mcols(bed)[,state] <- factor(mcols(bed)[,state],
                               levels = c("unmarked",
                                          "primed",
                                          "active"))
  plot(log2(mcols(bed)[,k4me1]),
       log2(mcols(bed)[,k27ac]),
       xlim = xlim,
       ylim = ylim,
       col = enhancer.colors[mcols(bed)[,state]],
       pch = 16,
       cex = 0.5)
  box(lwd = 8)
  MyText(paste(state,"\n","active:",sum(mcols(bed)[,state] == "active"),"\n",
               "primed:",sum(mcols(bed)[,state] == "primed"),"\n",
               "unmarked:",sum(mcols(bed)[,state] == "unmarked"),"\n"))
  return(bed)
  
}
Mywritegeneinf <- function(genelist,gene.output){
  MyWriteTable(cbind(genes.inf.input[unique(genelist),],
                     gene.tpm[unique(genelist),3:9]),
               gene.output
  )
}

MyWriteBed2 <- function(data,
                        output.file,
                        #output.col = NULL,
                        is.gz = F,
                        zero.based = T,
                        round.col = NULL,
                        round.digits = 4,
                        sample.name = NULL,
                        colors = "0,0,0",
                        col.names = F,
                        ...){
  data <- as.data.frame(data)
  data <- data[,]
  data[,round.col] <- round(data[,round.col],round.digits)
  if(zero.based){
    data[,"start"] <- data[,"start"] - 1
  }
  
  if(is.gz){
    output.pipe <- gzfile(output.file,
                          "w")
  }else{
    output.pipe <- output.file
  }
  
  if(!is.null(sample.name)){
    ucsc.txt <- paste("track type=bedGraph name=\"",
                      sample.name,
                      "\" description=\"",
                      sample.name,
                      "\" color=",
                      colors,
                      "\n",
                      sep = "")
    cat(ucsc.txt,
        file = output.pipe)
    options(scipen = 10)
    write.table(data,
                output.pipe,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t",
                append = T,
                col.names = col.names,
                ...)
    options(scipen = 0)
  }else{
    options(scipen = 10)
    write.table(data,
                output.pipe,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t",
                col.names = col.names,
                ...)
    options(scipen = 0)
  }
  if(is.gz){
    close(output.pipe)
  }
}

MyATACstate <- function(bed,
                        atac,
                        state,
                        atac.cutoff){
  mcols(bed)[,state] <- NA
  mcols(bed)[,state][mcols(bed)[,atac] < atac.cutoff] <- "close"
  mcols(bed)[,state][mcols(bed)[,atac] >= atac.cutoff] <- "open"
  
  MyText(paste(state,"\n","close:",sum(mcols(bed)[,state] == "close"),"\n",
               "open:",sum(mcols(bed)[,state] == "open"),"\n"))
  return(bed)
  
}

Myseuvioplot <- function(seu.ob,sym.list,type = 'Type',type.colors = time.colors,log=F,slot='data',assays='RNA',x.cex.axis=1,ylab.plot='TP0.1M',
                         ...){
  beta.preg.meta.tab <- seu.ob@meta.data
  if(log==F){
    for(j in sym.list){
      print(j)
      if(slot=='data'){beta.preg.meta.tab$tmp <- as.matrix(seu.ob[[assays]]@data)[j,beta.preg.meta.tab$SampleName]}
      if(slot=='scale.data'){beta.preg.meta.tab$tmp <- as.matrix(seu.ob[[assays]]@scale.data)[j,beta.preg.meta.tab$SampleName]}
      
      tmp.tpm <- expm1(beta.preg.meta.tab$tmp)
      names(tmp.tpm) <- beta.preg.meta.tab[,type]
      MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                          tmp.tpm,
                          color.violin = add.alpha(c(#"#1f78b4",
                            type.colors),alpha = 0.7),
                          box.lwd = 3,
                          ylab.plot  = ylab.plot,
                          j,x.cex.axis=x.cex.axis,
                          ...
                          
                          
      )
    }
  }
  
  if(log==T){
    for(j in sym.list){
      print(j)
      if(slot=='data'){beta.preg.meta.tab$tmp <- as.matrix(seu.ob[[assays]]@data)[j,beta.preg.meta.tab$SampleName]}
      
      if(slot=='scale.data'){beta.preg.meta.tab$tmp <- as.matrix(seu.ob[[assays]]@scale.data)[j,beta.preg.meta.tab$SampleName]}
      tmp.tpm <- beta.preg.meta.tab$tmp
      names(tmp.tpm) <- beta.preg.meta.tab[,type]
      MyViolinBeeSwarmMed(beta.preg.meta.tab[,type],
                          tmp.tpm,
                          color.violin = add.alpha(c(#"#1f78b4",
                            type.colors),alpha = 0.7),
                          box.lwd = 3,
                          ylab.plot  = ylab.plot,x.cex.axis=x.cex.axis,
                          j,
                          ...
                          
                          
      )
    }
  }
}

Myseuviopvalplot <- function(seu.ob,sym.list,type = 'Type',type.colors = time.colors,log=F,ylab.plot='TPO.1M',name1=1,name2=2,...){
  beta.preg.meta.tab <- seu.ob@meta.data
  tmp.median <- matrix(ncol = (length(names(table(seu.ob@meta.data[,type])))+1),nrow = length(sym.list))
  rownames(tmp.median) <- sym.list
  colnames(tmp.median) <- c('p-value',names(table(seu.ob@meta.data[,type])))
  
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
                          ylab.plot  = ylab.plot,
                          j,...
                          
                          
      )
      tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name1],'tmp'],
                           beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name2],'tmp'],
                           paired = F)
      tmp.median[j,'p-value'] <- tmp.p$p.value
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[name1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name1],'tmp'])-1),2)
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[name2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name2],'tmp'])-1),2)
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
                          ylab.plot  = ylab.plot,
                          j,...
                          
                          
      )
      tmp.p <- wilcox.test(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name1],'tmp'],
                           beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name2],'tmp'],
                           paired = F)
      tmp.median[j,'p-value'] <- tmp.p$p.value
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[name1]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name1],'tmp'])-1),2)
      tmp.median[j,names(table(seu.ob@meta.data[,type]))[name2]] <- round(median(exp(beta.preg.meta.tab[beta.preg.meta.tab$Type==names(table(seu.ob@meta.data[,type]))[name2],'tmp'])-1),2)
    }
  }
  
  return(tmp.median)
}

MyICA2Gg2 <- 
  function(cds,
           sample.inf,
           x.dim = 1,
           y.dim = 2,
           reduction.key="IC",
           color_manual = c(brewer.pal(9,"Set1"),
                            brewer.pal(12,"Set3"),
                            brewer.pal(8,"Set2"))){
    library(ggplot2)
    coord <- t(reducedDimS(fabp6.cds))
    coord <- as.data.frame(coord)
    
    df <- cbind(x.pos = coord[,x.dim],
                y.pos = coord[,y.dim],
                sample.inf[row.names(coord),])
    
    p <- ggplot(data = df,
                mapping = aes(x = x.pos,
                              y = y.pos,
                              label = rownames(sample.inf)))
    
    p <- p + xlab(paste(reduction.key,x.dim,sep = ""))
    p <- p + ylab(paste(reduction.key,y.dim,sep = ""))
    p <- p + guides(colour=guide_legend(title=NULL))
    p <- p + scale_color_manual(values = color_manual)
    
    p <- p +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1)+
      scale_color_manual(values = c(time.colors,"gray30")) +
      #scale_shape_manual(values = shape.group) +
      guides(colour = guide_legend(title = "",
                                   order = 1)) +
      theme(axis.text = element_text(size = 20, colour = "black")) +
      theme(axis.title.x = element_text(size = 40, colour = "black")) +
      theme(axis.title.y = element_text(size = 40, colour = "black")) +
      theme(legend.text = element_text(size = 45, colour = "black")) +
      theme(legend.title = element_text(size = 45, colour = "black")) +
      theme(panel.border = element_rect(size = 4,
                                        colour = "black")) +
      theme(legend.key = element_blank()) +
      guides(colour = guide_legend(title = "",
                                   keywidth = 3,
                                   keyheight = 3,
                                   override.aes = list(size = 12),
                                   order = 1,
                                   ncol = 1)) +
      guides(shape = guide_legend(title = "",
                                  keywidth = 3,
                                  keyheight = 3,
                                  override.aes = list(size = 12),
                                  order = 2,
                                  ncol = 1))
    return(p)
  }

Myseuratmarker <- function(seu.ob,marker.sym,reduction='umap',col=colors.exp,pt.size=1.5,plot.title='ln(TP0.1M+1)',...){
  
  marker.sym <- marker.sym[marker.sym %in% rownames(seu.ob)]
  length(marker.sym)#
  marker.sym <- unique(marker.sym)
  length(marker.sym)#57
  total.count <- length(marker.sym)
  run.count <- 1
  
  for (gene in c('nFeature_RNA','nCount_RNA',marker.sym)) {
    cat(paste(run.count, "/",total.count,"\n"))
    print(FeaturePlot(seu.ob,gene,pt.size = pt.size,reduction = reduction,...)+scale_color_gradientn(colours = col)+
            theme(plot.title = element_text(size = rel(3.5))) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.title = element_blank()) +
            theme(axis.text  = element_blank()) +
            theme(axis.ticks = element_blank()) +
            theme(legend.position="bottom")+
            theme(panel.border = element_rect(size = 4,
                                              colour = "black"))+
            guides(colour = guide_colorbar(title = plot.title,
                                           title.position = "top",
                                           barwidth = 23.6,
                                           title.hjust = 0.5,
                                           title.theme = element_text(angle = 0,
                                                                      size = 20),
                                           label.theme = element_text(angle = 0,
                                                                      size = 20),
                                           ticks = T))+ theme(aspect.ratio=1))
    run.count = run.count + 1
  }
}



Myseufindmarker <- function(seu.ob,gene.include= gene.include,ident.1,ident.2,c1,c2,logfc.threshold=0,diff='avg_log2FC',...){
  library(Seurat)
  seu.DEG <- FindMarkers(seu.ob,
                         features = gene.include,
                         ident.1 = ident.1,
                         ident.2 = ident.2,
                         logfc.threshold = logfc.threshold,
                         ... 
  )
  seu.DEG$gene <- rownames(seu.DEG)
  seu.DEG$cluster <- c1
  seu.DEG[seu.DEG[,diff]<0,'cluster'] <- c2
  seu.DEG$cluster <- factor(seu.DEG$cluster,levels = c(c1,c2))
  
  return(seu.DEG)
}

MyDEGfilterplot <- function(seu.ob,DEG.raw.filter,DEG.raw=NULL,prefix,
                            plotheatmap = F,display.cells = NULL,gene.text=1.2,
                            plotfireplot = F,heatmap.col1 = time.colors,heatmap.col2 = time.colors,fire.cols = time.colors[1:2],
                            heatmap.type = 'log.row.relat',group = c('Type','Rep'),plotRep=F,padj=F,pval.gap=2,line.lwd=3,y.cex=1.5,ext=0,x.ext=0,x.min.ext=0,x.gap=0.1,
                            width.input=8,height.input=8,point.size=1,rev.pc=T,col=T,heatmap.scale=F,dim=1,
                            ...){
  if(plotfireplot == T){
    DEG.raw$state <- '-'
    DEG.raw[DEG.raw.filter$gene,'state'] <- as.character(DEG.raw.filter$cluster)
    DEG.raw$state <- factor(DEG.raw$state,levels = c(names(table(DEG.raw.filter$cluster)),'-'))
    
    
    fire.cols <- c(fire.cols,'gray90')
    names(fire.cols) <- names(table(DEG.raw$state))
    # MyPlotColor(fire.cols,3)
    
    if(!padj){
    pdf(paste(prefix,"fireplot.pdf",sep = '.'),width.input,height.input)
    par(oma = c(0,1,0,0))
    #par(mar = c(0,0.002,0,0))
    x.max <- ceiling(max(abs(DEG.raw$avg_log2FC)))
    if(min(DEG.raw$p_val)==0){ymax=500}else{ymax <-  ceiling(max( -log10(DEG.raw$p_val)))}
    
    
    plot(-DEG.raw$avg_log2FC,
         -log10(DEG.raw$p_val),
         xlim = c(-x.max-x.ext,x.max+x.ext),
         ylim = c(0,(ymax+ext)),
         col = fire.cols[DEG.raw$state],
         pch = 16,
         cex = point.size,##pch size
         axes=FALSE,
         xlab = "log2(fold change)",
         ylab = "-log10(p-value)",
         cex.lab = 1.5
    )
    axis(2,pos=0,at=seq(0,(ymax+ext),by=pval.gap),lwd=line.lwd,cex.axis=y.cex)
    axis(1,at=seq(-x.max-x.ext,x.max+x.ext,by=x.gap),lwd=line.lwd,cex.axis=1.5,tck=0.02)##tck,>0 ticks inside,<0 ticks outside 
    plot(-DEG.raw$avg_log2FC,
         -log10(DEG.raw$p_val),
         xlim = c(-x.max-x.ext,x.max+x.ext),
         ylim = c(0,(ymax+ext)),
         col = fire.cols[DEG.raw$state],
         pch = 16,
         cex = point.size,##pch size
         axes=FALSE,
         xlab = "log2(fold change)",
         ylab = "-log10(p-value)",
         cex.lab = 1.5
    )
    axis(2,pos=0,at=seq(0,(ymax+ext),by=pval.gap),lwd=line.lwd,cex.axis=y.cex)
    axis(1,at=seq(-x.max-x.ext,x.max+x.ext,by=x.gap),lwd=line.lwd,cex.axis=1.5,tck=0.02)##tck,>0 ticks inside,<0 ticks outside
    
    text(-DEG.raw.filter[,'avg_log2FC'],
         -log10(DEG.raw.filter[,'p_val']),
         labels=rownames(DEG.raw.filter),
         cex=gene.text)
    
    MyText(paste(
      names(fire.cols)[1],sum(DEG.raw$state == names(fire.cols)[1]),"\n",
      names(fire.cols)[2],sum(DEG.raw$state == names(fire.cols)[2]),"\n"),...)
    
    dev.off()
  }
    if(padj){
      pdf(paste(prefix,"padj.fireplot.pdf",sep = '.'),width.input,height.input)
      par(oma = c(0,1,0,0))
      #par(mar = c(0,0.002,0,0))
      x.max <- ceiling(max(abs(DEG.raw$avg_log2FC)))
      if(min(DEG.raw$p_val)==0){ymax=500}else{ymax <-  ceiling(max( -log10(DEG.raw$p_val_adj)))}
      
      plot(-DEG.raw$avg_log2FC,
           -log10(DEG.raw$p_val_adj),
           xlim = c(-x.max-x.ext,x.max+x.ext),
           ylim = c(0,(ymax+ext)),
           col = fire.cols[DEG.raw$state],
           pch = 16,
           cex = point.size,##pch size
           axes=FALSE,
           xlab = "log2(fold change)",
           ylab = "-log10(p_value_adj)",
           cex.lab = 1.5
      )
      axis(2,pos=0,at=seq(0,(ymax+ext),by=pval.gap),lwd=line.lwd,cex.axis=y.cex)
      axis(1,at=seq(-x.max-x.ext,x.max+x.ext,by=x.gap),lwd=line.lwd,cex.axis=1.5,tck=0.02)##tck,>0 ticks inside,<0 ticks outside 
      plot(-DEG.raw$avg_log2FC,
           -log10(DEG.raw$p_val_adj),
           xlim = c(-x.max-x.ext,x.max+x.ext),
           ylim = c(0,(ymax+ext)),
           col = fire.cols[DEG.raw$state],
           pch = 16,
           cex = point.size,##pch size
           axes=FALSE,
           xlab = "log2(fold change)",
           ylab = "-log10(p_value_adj)",
           cex.lab = 1.5
      )
      axis(2,pos=0,at=seq(0,(ymax+ext),by=pval.gap),lwd=line.lwd,cex.axis=y.cex)
      axis(1,at=seq(-x.max-x.ext,x.max+x.ext,by=x.gap),lwd=line.lwd,cex.axis=1.5,tck=0.02)##tck,>0 ticks inside,<0 ticks outside
      
      text(-DEG.raw.filter[,'avg_log2FC'],
           -log10(DEG.raw.filter[,'p_val_adj']),
           labels=rownames(DEG.raw.filter),
           cex=gene.text)
      
      MyText(paste(
        names(fire.cols)[1],sum(DEG.raw$state == names(fire.cols)[1]),"\n",
        names(fire.cols)[2],sum(DEG.raw$state == names(fire.cols)[2]),"\n"),...)
      
      dev.off()
    }
  }
  else if(plotheatmap == T){
    seurat.data.tmp = subset(seu.ob,
                             cells = display.cells)
    seurat.data.tmp <- ScaleData(seurat.data.tmp,feature=DEG.raw.filter$gene)
    seurat.data.tmp = RunPCA(seurat.data.tmp,feature=DEG.raw.filter$gene,npcs = 5)
    si.tmp.sort = seurat.data.tmp@meta.data[order(seurat.data.tmp@reductions$pca@cell.embeddings[,dim],
                                                  decreasing = rev.pc),]
    
    
    
    
    heatmap.data <- seu.ob@assays$RNA@data[DEG.raw.filter$gene,
                                                rownames(si.tmp.sort)]
    heatmap.data = as.matrix(heatmap.data)
    
    c1Color <- MyName2Col(si.tmp.sort[,group[1]],
                          heatmap.col1)
    c1Color <- as.matrix(c1Color)
    c2Color <- MyName2Col(si.tmp.sort[,group[2]],
                          heatmap.col2)
    c2Color <- as.matrix(c2Color)
    cColor <- cbind(c1Color,c2Color)
    ccColor <- c1Color
    if(plotRep == T){
      ccColor <- cColor
    }
   if(col==TRUE){ png(paste(prefix,'DEG.heatmap.png',sep = ''),2000,3000)
    DEG.row.tree <-
      MyHeatmap(as.matrix(exp(heatmap.data)),
                type = heatmap.type,
                hc.c.data.type = "log.row.relat",
                hc.r.data.type = "log.row.relat",
                c.cov.method = "s",
                r.cov.method = "s",
                c.hc.method = "ward.D2",
                r.hc.method = "ward.D2",
                ColSideColors = ccColor,
                ColSideColorsSize = 1.5,
                return.tree = "row",
                graph = T,
                ...)
    dev.off()}
    if(col==FALSE){
   if(heatmap.scale==T){   
     seurat.data.tmp <- ScaleData(seurat.data.tmp,feature=rownames(heatmap.data))
     png(paste(prefix,'.scale.DEG.heatmap.png',sep = ''),2000,3000)
      DEG.row.tree <-
        MyHeatmap(as.matrix(exp(seurat.data.tmp@assays$RNA@scale.data[rownames(heatmap.data),colnames(heatmap.data)])),
                  type = heatmap.type,
                  color.palette = pc12.heatmap.col,
                  hc.c.data.type = "log.row.relat",
                  hc.r.data.type = "log.row.relat",
                  c.cov.method = "s",
                  r.cov.method = "s",
                  Colv = 'none',
                  c.hc.method = "ward.D2",
                  r.hc.method = "ward.D2",
                  ColSideColors = ccColor,
                  ColSideColorsSize = 1.5,
                  return.tree = "row",
                  graph = T,
                  ...)
      dev.off()}
      if(heatmap.scale==F){   png(paste(prefix,'.DEG.heatmap.png',sep = ''),2000,3000)
        DEG.row.tree <-
          MyHeatmap(as.matrix(exp(heatmap.data)),
                    type = heatmap.type,
                    
                    hc.c.data.type = "log.row.relat",
                    hc.r.data.type = "log.row.relat",
                    c.cov.method = "s",
                    r.cov.method = "s",
                    Colv = 'none',
                    c.hc.method = "ward.D2",
                    r.hc.method = "ward.D2",
                    ColSideColors = ccColor,
                    ColSideColorsSize = 1.5,
                    return.tree = "row",
                    graph = T,
                    ...)
        dev.off()}
      
      }
    return(DEG.row.tree)
  }
}

MyPseudotimebox <- function(si.tab,pseudotime,Type,time.colors  = time.colors,rep = F,size.point=4,box.lwd=2,shape.type=c(19,25),...){
  if(rep == F){
    p.time <- ggplot(data = si.tab,
                     aes(y=pseudotime,
                         x=as.numeric(Type), 
                         colour = Type#,
                         #shape = proliferative#,
                         #size = 1
                     )
    ) +
      scale_shape_manual(values = shape.type) +
      scale_color_manual(values = time.colors) +
      theme(axis.text = element_text(size = 30, colour = "black")) +
      theme(axis.title.x = element_text(size = 40, colour = "black")) +
      theme(axis.title.y = element_text(size = 40, colour = "black")) +
      theme(legend.text = element_text(size = 40, colour = "black")) +
      theme(legend.title = element_text(size = 30, colour = "black")) +
      xlab("") +
      ylab("") +
      guides(colour=guide_legend(title=NULL)) +
      theme_bw() +
      theme(legend.key = element_blank()) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 4,
                                        colour = "black"))
    library(ggbeeswarm)
    p <- p.time +
      # scale_x_reverse()+
      # #scale_shape_manual(values = c(16,1)) +
      # geom_point(aes(x = -Pseudotime,
      #                y = -as.numeric(Age),
      #                col = Age
      #                #shape = Ngn3Exp
      # ),
      # na.rm = TRUE,
      # size = 2) +
      geom_quasirandom(groupOnX = TRUE,
                       bandwidth = 1,
                       dodge.width = 0.6,
                       alpha = 0.8,
                       size = size.point,
                       #show.legend = T,
                       stroke = 0#,#setting stroke to 0 removes the outline around the points
                       #nbins = 30
                       #varwidth = T#,
                       #method = "pseudorandom"
      ) +
      geom_boxplot(aes(fill = Type),
                   alpha = 0.3,
                   #position = dodge,
                   color = "grey30",
                   lwd=box.lwd,
                   outlier.shape = NA,
                   show.legend = F#,
                  # width = box.width
      )+
      scale_fill_manual(values = time.colors)+
      
      
      guides(colour = guide_legend(title = "",
                                   keywidth = 3,
                                   keyheight = 3,
                                   override.aes = list(size=12),
                                   order = 1,
                                   ncol = 1))+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.text = element_text(size = 10, colour = "black"))
  } 
  else if(rep == T){
    p.time <- ggplot(data = si.tab,
                     aes(y=pseudotime,
                         x=as.numeric(Type_rep), 
                         colour = Type_rep,
                        # shape = proliferative#,
                         #size = 1
                     )
    ) +
      scale_shape_manual(values = c(17,19)) +
      scale_color_manual(values = time.colors) +
      theme(axis.text = element_text(size = 30, colour = "black")) +
      theme(axis.title.x = element_text(size = 40, colour = "black")) +
      theme(axis.title.y = element_text(size = 40, colour = "black")) +
      theme(legend.text = element_text(size = 40, colour = "black")) +
      theme(legend.title = element_text(size = 30, colour = "black")) +
      xlab("") +
      ylab("") +
      guides(colour=guide_legend(title=NULL)) +
      theme_bw() +
      theme(legend.key = element_blank()) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 4,
                                        colour = "black"))
    library(ggbeeswarm)
    p <- p.time +
      # scale_x_reverse()+
      # #scale_shape_manual(values = c(16,1)) +
      # geom_point(aes(x = -Pseudotime,
      #                y = -as.numeric(Age),
      #                col = Age
      #                #shape = Ngn3Exp
      # ),
      # na.rm = TRUE,
      # size = 2) +
      geom_quasirandom(groupOnX = TRUE,
                       bandwidth = 1,
                       dodge.width = 0.6,
                       alpha = 0.8,
                       size = 4,
                       #show.legend = T,
                       stroke = 0#,#setting stroke to 0 removes the outline around the points
                       #nbins = 30
                       #varwidth = 10#,
                       #method = "pseudorandom"
      ) +
      geom_boxplot(aes(fill = Type_rep),
                   alpha = 0.3,
                   #position = dodge,
                   color = "grey30",
                   outlier.shape = NA,
                   show.legend = F,
                   lwd = box.lwd
      )+
      scale_fill_manual(values = time.colors)+
      
      
      guides(colour = guide_legend(title = "",
                                   keywidth = 3,
                                   keyheight = 3,
                                   override.aes = list(size=12),
                                   order = 1,
                                   ncol = 1))+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.text = element_text(size = 10, colour = "black"))
  }
  return(p)
}



Mydiffusionmap <- 
  function(tpm,neigen = 50,...){
    library(diffusionMap)
    D <- parallelDist::parDist(t(tpm))
    dif <- diffuse(D,
                   neigen = neigen,
                   ...)
    return(dif)
  }


MySeuratDR2Gg2 <- 
  function(seurat.data,
           sample.inf,
           x.dim = 1,
           y.dim = 2,
           reduction.use="tsne",
           reduction.key="tSNE",
           estimate.variation.explain.percentage=F,
           color_manual = c(brewer.pal(9,"Set1"),
                            brewer.pal(12,"Set3"),
                            brewer.pal(8,"Set2"))){
    library(ggplot2)
    coord <- seurat.data@reductions[[reduction.use]]@cell.embeddings
    coord <- as.data.frame(coord)
    
    df <- cbind(x.pos = coord[,x.dim],
                y.pos = coord[,y.dim],
                sample.inf[row.names(coord),]
                )
    
    p <- ggplot(data = df,
                mapping = aes(x = x.pos,
                              y = y.pos,
                              label = SampleName
                              ))
    
    if(estimate.variation.explain.percentage){
      coord.x.VEP=(seurat.data@reductions[[reduction.use]]@stdev[x.dim])^2/
        seurat.data@reductions[[reduction.use]]@misc$total.variance*100
      coord.y.VEP=(seurat.data@reductions[[reduction.use]]@stdev[y.dim])^2/
        seurat.data@reductions[[reduction.use]]@misc$total.variance*100
      coord.x.VEP=paste("(",round(coord.x.VEP,1),"%)",sep = "")
      coord.y.VEP=paste("(",round(coord.y.VEP,1),"%)",sep = "")
    }else{
      coord.x.VEP=NULL
      coord.y.VEP=NULL
    }
    p <- p + xlab(paste(reduction.key,x.dim,coord.x.VEP,sep = ""))
    p <- p + ylab(paste(reduction.key,y.dim,coord.y.VEP,sep = ""))
    p <- p + guides(colour=guide_legend(title=NULL))
    p <- p + scale_color_manual(values = color_manual)
    
    p <- p +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1)+
      scale_color_manual(values = c(time.colors,"gray30")) +
      #scale_shape_manual(values = shape.group) +
      guides(colour = guide_legend(title = "",
                                   order = 1)) +
      theme(axis.text = element_text(size = 20, colour = "black")) +
      theme(axis.title.x = element_text(size = 40, colour = "black")) +
      theme(axis.title.y = element_text(size = 40, colour = "black")) +
      theme(legend.text = element_text(size = 45, colour = "black")) +
      theme(legend.title = element_text(size = 45, colour = "black")) +
      theme(panel.border = element_rect(size = 4,
                                        colour = "black")) +
      theme(legend.key = element_blank()) +
      guides(colour = guide_legend(title = "",
                                   keywidth = 3,
                                   keyheight = 3,
                                   override.aes = list(size = 12),
                                   order = 1,
                                   ncol = 1)) +
      guides(shape = guide_legend(title = "",
                                  keywidth = 3,
                                  keyheight = 3,
                                  override.aes = list(size = 12),
                                  order = 2,
                                  ncol = 1))
    return(p)
  }


MyLogBoxPlot2 <- 
  function(data,
           min.cutoff = min.cutoff,
           log = T,
           cex.axis = 2,cex.lab = 2,cex.main = 2,
           ...){
    data.na <- data
    data.na[data.na < min.cutoff] <- 0
    if(log ==  T){
      boxplot(log2(data.na),
              outline = F,
              ...)
    }
    else{
      boxplot(data.na,
              outline = F,
              ...)
    }
  }

MyLogVioPlot2 <- 
  function(data,
           min.cutoff = min.cutoff,
           log = T,
           cex.axis = 2,cex.lab = 2,cex.main = 2,
           ...){
    data.na <- data
    data.na[data.na < min.cutoff] <- 0
    if(log ==  T){
      vioplot::vioplot(log2(data.na),
                       outline = F,
                       ...)
    }
    else{
      vioplot::vioplot(data.na,
                       outline = F,
                       ...)
    }
  }



MyLogVioPlot <- 
function(data,
         min.cutoff = min.cutoff,
         log = T,
         cex.axis = 2,cex.lab = 2,cex.main = 2,
         ...){
  data.na <- data
  data.na[data.na < min.cutoff] <- NA
  if(log ==  T){
    vioplot::vioplot(log2(data.na),
            outline = F,
            ...)
  }
  else{
    vioplot::vioplot(data.na,
            outline = F,
            ...)
  }
}




MyCalTp0.1mExcludeGene <- 
function(read.count.data,
         gene.length,
         gene.exclude = NA){
  trans.data <- read.count.data / gene.length * 1000
  trans.sum <- colSums(trans.data[!rownames(trans.data) %in% gene.exclude,])
  tpm.data <- t(t(trans.data) / trans.sum * 100000)
  return(tpm.data)
}

######by Qiu########
MyCalUcscNormFactor = function(merge.bed,
                               list.peak.tmp,
                               max.cutoff= 200,
                               quantile.cutoff = 0.5){
  bed.merge.common.tmp = merge.bed
  for (sn.tmp in names(list.peak.tmp)) {
    bed.merge.common.tmp <- MyBedGrepU(bed.merge.common.tmp,
                                       list.peak.tmp[[sn.tmp]])
  }
  
  bed.merge.common.quantile.tmp = bed.merge.common.tmp
  for (sn.tmp in names(list.peak.tmp)) {
    bed.merge.common.quantile.tmp =
      bed.merge.common.quantile.tmp[mcols(bed.merge.common.quantile.tmp)[,sn.tmp] > quantile(mcols(bed.merge.common.tmp)[,sn.tmp],quantile.cutoff) & 
                                      mcols(bed.merge.common.quantile.tmp)[,sn.tmp] <  max.cutoff]
    
  }
  cat("filtered peak count:",
      length(bed.merge.common.quantile.tmp),
      "\n")
  norm.factor = apply(mcols(bed.merge.common.quantile.tmp)[,names(list.peak.tmp)], 2, mean)
  norm.factor = norm.factor / mean(norm.factor)
}

########
MyPca2Plotly3d.supp <-
  function(pca.res,
           x.dim = 1,
           y.dim = 2,
           z.dim = 3,
           ...){
    library(plotly)
    p <- plot_ly(x = c(pca.res$ind$coord[,x.dim],pca.res$ind.sup$coord[,x.dim]),
                 y = c(pca.res$ind$coord[,y.dim],pca.res$ind.sup$coord[,y.dim]),
                 z = c(pca.res$ind$coord[,z.dim],pca.res$ind.sup$coord[,z.dim]),
                 type = "scatter3d",
                 mode = "markers",
                 ...)
    return(p)
    detach("package:plotly",
           unload=TRUE)
  }

MyTopCountGene <- function(src.data,
                           ratio = 0.05,
                           n = 20,
                           ...){
  topcount.gene <- rownames(src.data)[colSums(t(src.data@assays$RNA@counts)/colSums(src.data@assays$RNA@counts)>ratio)>n]
  return(topcount.gene)
}


MyCalumiTpmExcludeGene <- 
  function(umi.data,gene.exclude = NA){
    trans.sum <- colSums(umi.data[!rownames(umi.data) %in% gene.exclude,])
    tpm.data <- t(t(umi.data) / trans.sum * 100000)
    tpm.data <- round(tpm.data,2)
    return(tpm.data)
}

MyCalTpmExcludeGene <- 
  function(read.count.data,
           gene.length,
           gene.exclude = NA){
    trans.data <- read.count.data / gene.length * 1000
    trans.sum <- colSums(trans.data[!rownames(trans.data) %in% gene.exclude,])
    tpm.data <- t(t(trans.data) / trans.sum * 100000)
    return(tpm.data)
}

MyTimeExp <- 
function(xpt,
         ypt,
         span =.75,
         pch = 20,
         xlab = "Pseudotime",
         ylab = "Expression (TPM)",
         col.line = "#BC2B2B",
         col.polygon = "#00529533",
         bty = "n",
         col = "black",
         ...){
  y.loess <-loess(ypt~xpt,span=span)
  ypt.predict <- predict(y.loess,sort(xpt),se=T)
  plot(xpt,
       ypt,
       bty = bty,
       xlab = xlab,
       ylab = ylab, 
       pch = pch,
       col = col,
       ...)
  lines(sort(xpt),ypt.predict$fit,lwd=1,col = col.line)
  
  y.polygon <- c((ypt.predict$fit+1.96*ypt.predict$se.fit), rev(c(ypt.predict$fit-1.96*ypt.predict$se.fit)))
  x.polygon <- c(sort(xpt), rev(sort(xpt)))
  polygon(x.polygon, 
          y.polygon, 
          col = col.polygon, 
          border=NA)
}
MyTimeExp.single <- 
  function(xpt,
           ypt.matrix,
           gene,
           xlab = "Pseudotime",
           ylab = "Relative Expression",
           col.line = "gray60",
           col.mainline = "red",
           y.lim = c(0,1),
           y.tck = -0.02,
           ...){
    ypt.matrix <- ypt.matrix[,order(xpt)]
    xpt <- sort(xpt)
    
    ypt.matrix.relat <- ypt.matrix / rowMax(as.matrix(ypt.matrix))
    
    # pdf("bc.hl.pdf")
    plot(range(xpt),
         y.lim,
         col = "white",
         bty = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab,
         axes=FALSE,
         # lwd=2,
         ...)
    axis(1,at = c(min(xpt),max(xpt)),labels = NULL,cex.axis = 1.5,lwd=3,tck=0)
    axis(2,at = seq(y.lim[1],y.lim[2],by=0.2),cex.axis = 1.5,lwd=3,tck=y.tck)
    #box(lwd=3)
      ypt <- as.numeric(ypt.matrix.relat[gene,])
      y.loess <-loess(ypt~xpt)
      ypt.predict <- predict(y.loess,
                             data.frame(x=xpt),
                             se=T)
      lines(xpt,
            ypt.predict$fit,
            lwd = 1,
            col = col.line)

    
  }

MyTimeExp.merge <- 
function(xpt,
         ypt.matrix,
         xlab = "Pseudotime",
         ylab = "Relative Expression",
         col.line = "gray60",
         col.mainline = "red",
         y.lim = c(0,1),
         y.tck = -0.02,
         x.label=NULL,
         ...){
  ypt.matrix <- ypt.matrix[,order(xpt)]
  xpt <- sort(xpt)
  
  ypt.matrix.relat <- ypt.matrix / rowMax(as.matrix(ypt.matrix))
  
  # pdf("bc.hl.pdf")
  plot(range(xpt),
       y.lim,
       
       col = "white",
       bty = "n",
       xaxt = "n",
       xlab = xlab,
       ylab = ylab,
       axes=FALSE,
       #xaxt="n",
      # lwd=2,
       ...)
  axis(1,at = seq(min(xpt),max(xpt),by=1),labels = x.label,cex.axis = 1.5,lwd=3,tck=y.tck)
  axis(2,at = seq(y.lim[1],y.lim[2],by=0.2),cex.axis = 1.5,lwd=3,tck=y.tck)
  #box(lwd=3)
  for(n.plot in 1:nrow(ypt.matrix.relat)){
    ypt <- as.numeric(ypt.matrix.relat[n.plot,])
    y.loess <-loess(ypt~xpt)
    ypt.predict <- predict(y.loess,
                           data.frame(x=xpt),
                           se=T)
    lines(xpt,
          ypt.predict$fit,
          lwd = 1,
          col = col.line)
  }
  
  ypt <- colMeans(ypt.matrix.relat)
  y.loess <-loess(ypt~xpt)
  ypt.predict <- predict(y.loess,data.frame(x=xpt),se=T)
  lines(xpt,
        ypt.predict$fit, 
        lwd = 8,
        col = col.mainline)
}

MyKeggzebrafish <- 
  function(selected.id,
           universe.id,
           #organism = "mmu",
           pval.adj=0.1){
    
    library(clusterProfiler)
    # library(org.Dr.eg.db)
    # OrgDb <- org.Dr.eg.db
    kegg=enrichKEGG(gene = selected.id,organism = 'zebrafish',
                    keyType = "kegg",universe = universe.id,
                    pvalueCutoff =pval.adj,qvalueCutoff = 4*pval.adj)
    
    detach("package:clusterProfiler", unload=TRUE)
    return(as.data.frame(kegg))
  }
MykeggZebrafishwritetable <- function(genelist,
                                      universe.gene=zebrafish.gene.inf$ENTREZID,
                                      pval.adj=0.1,
                                      kegg.out,
                                      ...){
  kegg.result <- MyKeggzebrafish(genelist,
                                 universe.gene,
                                 pval.adj = pval.adj,
                                 ...
  )
  MyNCBI2Symbol <- function(mykeggid){
    keggid2symbol <- 1:length(mykeggid)
    for(i in 1:length(mykeggid)){
      keggid <- strsplit(mykeggid[i],"/")
      keggid <- keggid[[1]]
      mykeggid2symbol <- zebrafish.gene.inf[zebrafish.gene.inf$ENTREZID %in% keggid,"Symbol"]
      mykeggid2symbol <- paste(mykeggid2symbol,collapse = "/")
      keggid2symbol[i] <- mykeggid2symbol
    }
    return(keggid2symbol)
  }
  kegg.result$Symbol <- NA
  kegg.result$Symbol <- MyNCBI2Symbol(kegg.result$geneID)
  MyWriteTable(kegg.result,
               kegg.out)
  return(kegg.result)
}


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
         ylim = c(min(y),(max(y)+ext)),
         xaxt="n",
         yaxt="n",
         main = title.plot,
         cex.main = cex.main,
         
         ...)
    box(lwd = box.lwd)
    axis(1,lwd = 0,lwd.ticks = 0, at = 1:length(group.all),cex.axis = x.cex.axis,labels = group.all)
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


MyGSEA.prepare2 <- 
  function(tpm,
           sample.class,
           gct.name,
           cls.name,
           mode="discrete",
           continuous_value_name=NULL){
    print("use EnsemblGeneID")
    gct=cbind(NAME=toupper(genes.inf.input[rownames(tpm),'Symbol']),
              Description=NA,
              tpm)
    
    gct.file=file(gct.name, open = "wt")
    writeLines("#1.2",gct.file)
    writeLines(paste(nrow(tpm),ncol(tpm),sep = "\t"),gct.file)
    write.table(gct,gct.file,append = T,quote = F,row.names = F,sep = "\t")
    close(gct.file)
    
    if(mode=="discrete"){
      cls.file=file(cls.name, open = "wt")
      writeLines(paste(ncol(tpm),
                       length(table(as.character(sample.class))),1,
                       sep = " "),
                 cls.file)
      writeLines(paste(c("#",names(table(as.character(sample.class)))),
                       collapse = " "),
                 cls.file)
      writeLines(paste(sample.class,collapse = " "),cls.file)
      close(cls.file)
    }else if(mode=="continuous"){
      cls.file=file(cls.name, open = "wt")
      writeLines("#numeric",
                 cls.file)
      if(is.null(continuous_value_name)){
        writeLines("#Trajectory",
                   cls.file)
      }else{
        writeLines(paste("#",continuous_value_name,sep = ""),
                   cls.file)
      }
      writeLines(paste(sample.class,collapse = " "),cls.file)
      close(cls.file)
    }
  }

#######



MycalnewReadcount <- 
  function(raw.read.count,
           spliced.mat,
           unspliced.mat,
           #offset.mat,
           selected.genes,
           cell.dist,
           kCells=20,fit.quantile=0.02,n.cores = 1,...){
    
    rvel.cd <- gene.relative.velocity.estimates(as(as.matrix(spliced.mat),"dgCMatrix")[selected.genes,],
                                                as(as.matrix(unspliced.mat),"dgCMatrix")[selected.genes,],
                                                # smat = as(as.matrix(offset.mat),"dgCMatrix"),
                                                cell.dist=cell.dist,
                                                kCells=kCells,fit.quantile=fit.quantile,n.cores = n.cores,...)
    
    raw.read.count.selected.genes=raw.read.count[selected.genes,]
    cellsize=colSums(raw.read.count.selected.genes)
    add=t(t(rvel.cd$deltaE)*cellsize/10000)*1
    new.read.count.selected.genes=raw.read.count.selected.genes[rownames(add),]+as.matrix(add)
    
    new.read.count=raw.read.count
    new.read.count[rownames(new.read.count.selected.genes),]=new.read.count.selected.genes
    new.read.count[new.read.count<0]=0
    return(new.read.count)
  }

#####QWL########
MyFastMNN <- function(data.log,
                      batch.factor,
                     # data.log.batch.list,
                      subset.row=NULL,
                      auto.order=FALSE,
                      d=50,
                      k=20,
                      compute.variances=FALSE, 
                      cos.norm=TRUE, 
                      ndist=3,  
                      approximate=FALSE, 
                      irlba.args=list()){
  if (!is.factor(batch.factor) || length(x = batch.factor) != ncol(data.log)) {
    stop("'batch.factor' must be a factor of length ncol(data.log)")
  }
  
  library(scran)
  
  fastMNN2 <- function (batches, k = 20, cos.norm = TRUE, ndist = 3, d = 50, approximate = FALSE, 
                        irlba.args = list(), subset.row = NULL, auto.order = FALSE, 
                        pc.input = FALSE, compute.variances = FALSE, assay.type = "logcounts", 
                        get.spikes = FALSE, BNPARAM = BiocNeighbors:::KmknnParam(), BPPARAM = BiocParallel:::SerialParam()) {
    .Deprecated("batchelor::fastMNN")
    # batches <- list(...)
    nbatches <- length(batches)
    if (nbatches < 2L) {
      stop("at least two batches must be specified")
    }
    if (!pc.input) {
      out <- scran:::.SCEs_to_matrices(batches, assay.type = assay.type, 
                                       subset.row = subset.row, get.spikes = get.spikes)
      batches <- out$batches
      subset.row <- out$subset.row
      if (!is.null(subset.row) && !identical(subset.row, seq_len(nrow(batches[[1]])))) {
        batches <- lapply(batches, "[", i = subset.row,
                          drop = FALSE)
      }
      if (cos.norm) {
        batches <- lapply(batches, FUN = cosineNorm, mode = "matrix")
      }
      pc.mat <- scran:::.multi_pca(batches, approximate = approximate, 
                                   irlba.args = irlba.args, d = d, use.crossprod = TRUE)
    }
    else {
      scran:::.check_batch_consistency(batches, byrow = FALSE)
      pc.mat <- batches
    }
    var.kept <- rep(1, nbatches)
    re.order <- NULL
    if (!is.logical(auto.order)) {
      re.order <- as.integer(auto.order)
      if (!identical(sort(re.order), seq_len(nbatches))) {
        stop("integer 'auto.order' must contain a permutation of 1:nbatches")
      }
      auto.order <- FALSE
    }
    use.order <- !is.null(re.order)
    mnn.pairings <- vector("list", nbatches - 1L)
    for (bdx in 2:nbatches) {
      if (auto.order) {
        if (bdx == 2L) {
          d.out <- scran:::.define_first_merge(pc.mat, k = k, BNPARAM = BNPARAM, 
                                               BPPARAM = BPPARAM)
          refdata <- pc.mat[[d.out$first]]
          curdata <- pc.mat[[d.out$second]]
          mnn.sets <- d.out$pairs
          precomp <- d.out$precomputed
          processed <- c(d.out$first, d.out$second)
        }
        else {
          d.out <- scran:::.define_next_merge(refdata, pc.mat, 
                                              processed, precomp, k = k, BNPARAM = BNPARAM, 
                                              BPPARAM = BPPARAM)
          curdata <- pc.mat[[d.out$other]]
          mnn.sets <- d.out$pairs
          processed <- c(processed, d.out$other)
        }
      }
      else {
        if (bdx == 2L) {
          ref.idx <- if (use.order) 
            re.order[1]
          else 1L
          refdata <- pc.mat[[ref.idx]]
          processed <- ref.idx
        }
        cur.idx <- if (use.order) 
          re.order[bdx]
        else bdx
        curdata <- pc.mat[[cur.idx]]
        mnn.sets <- scran:::find.mutual.nn(refdata, curdata, k1 = k, 
                                           k2 = k, BNPARAM = BNPARAM, BPPARAM = BPPARAM)
        processed <- c(processed, cur.idx)
      }
      ave.out <- scran:::.average_correction(refdata, mnn.sets$first, 
                                             curdata, mnn.sets$second)
      overall.batch <- colMeans(ave.out$averaged)
      if (compute.variances) {
        var.before <- scran:::.compute_intra_var(refdata, curdata, 
                                                 pc.mat, processed)
      }
      refdata <- scran:::.center_along_batch_vector(refdata, overall.batch)
      curdata <- scran:::.center_along_batch_vector(curdata, overall.batch)
      if (compute.variances) {
        var.after <- scran:::.compute_intra_var(refdata, curdata, 
                                                pc.mat, processed)
        var.kept[seq_len(bdx)] <- var.kept[seq_len(bdx)] * 
          var.after/var.before
      }
      re.ave.out <- scran:::.average_correction(refdata, mnn.sets$first, 
                                                curdata, mnn.sets$second)
      curdata <- scran:::.tricube_weighted_correction(curdata, re.ave.out$averaged, 
                                                      re.ave.out$second, k = k, ndist = ndist, BNPARAM = BNPARAM, 
                                                      BPPARAM = BPPARAM)
      mnn.pairings[[bdx - 1L]] <- DataFrame(first = mnn.sets$first, 
                                            second = mnn.sets$second + nrow(refdata))
      refdata <- rbind(refdata, curdata)
    }
    if (auto.order || use.order) {
      ordering <- vector("list", nbatches)
      last <- 0L
      for (idx in processed) {
        ncells <- nrow(pc.mat[[idx]])
        ordering[[idx]] <- last + seq_len(ncells)
        last <- last + ncells
      }
      ordering <- unlist(ordering)
      refdata <- refdata[ordering, , drop = FALSE]
      relocate <- ordering
      relocate[ordering] <- seq_along(ordering)
      for (x in seq_along(mnn.pairings)) {
        current <- mnn.pairings[[x]]
        current$first <- relocate[current$first]
        current$second <- relocate[current$second]
        mnn.pairings[[x]] <- current
      }
      var.kept[processed] <- var.kept
    }
    if (!is.null(names(batches))) {
      batch.labels <- names(batches)
    }
    else {
      batch.labels <- seq_along(batches)
    }
    ncells <- vapply(pc.mat, FUN = nrow, FUN.VALUE = 0L)
    batch.ids <- Rle(batch.labels, ncells)
    output <- list(corrected = refdata, batch = batch.ids, pairs = mnn.pairings, 
                   order = batch.labels[processed])
    if (!pc.input) {
      output$rotation <- metadata(pc.mat)$rotation
    }
    if (compute.variances) {
      output$lost.var <- 1 - var.kept
    }
    return(output)
  }

  data.log.batch.list <- list()
  sn.data <- c()
  for (n.tmp in 1:length(levels(batch.factor))){
    print(n.tmp)
    data.log.batch.list[[n.tmp]] <- data.log[,batch.factor == levels(batch.factor)[n.tmp]]
    sn.data <- c(sn.data,
                 colnames(data.log[,batch.factor == levels(batch.factor)[n.tmp]]))
  }

  fMNN.res <- fastMNN2(data.log.batch.list,
                       subset.row=subset.row,
                       auto.order=auto.order,
                       d=d,
                       k=k,
                       compute.variances=compute.variances, 
                       cos.norm=cos.norm, 
                       ndist=ndist,  
                       approximate=approximate, 
                       irlba.args=irlba.args)
  
  
  rownames(fMNN.res$corrected) <- sn.data
  colnames(fMNN.res$corrected) <- paste("PC",1:d,sep = "")
  if (is.null(subset.row)) {
    rownames(fMNN.res$rotation) <- rownames(data.log)
  } else  {
    rownames(fMNN.res$rotation) <- subset.row
  }
  colnames(fMNN.res$rotation) <- paste("PC",1:d,sep = "")
  return(fMNN.res)
}




MySeuratFDL2Gg <- function(layout.3d.fdl,
                           sample.inf,
                           x.dim = 1,
                           y.dim = 2,
                           color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  row.names(sample.inf) <- sample.inf$SampleName
  tSNE.df <- cbind(x.pos = layout.3d.fdl[,x.dim],
                   y.pos = layout.3d.fdl[,y.dim],
                   sample.inf)
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("FDL",x.dim,sep = ""))
  p <- p + ylab(paste("FDL",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}



MySeuratUMAP2Gg <- 
  function(seurat,
           sample.inf,
           x.dim = 1,
           y.dim = 2,
           color_manual = brewer.pal(9,"Set1")){
    library(ggplot2)
    tSNE.coord <- seurat@dr$umap@cell.embeddings
    tSNE.coord <- as.data.frame(tSNE.coord)
    
    row.names(sample.inf) <- sample.inf$SampleName
    tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                     y.pos = tSNE.coord[,y.dim],
                     sample.inf[row.names(tSNE.coord),])
    
    p <- ggplot(data = tSNE.df, 
                mapping = aes(x = x.pos, 
                              y = y.pos,
                              label = SampleName))
    p <- p + xlab(paste("UMAP",x.dim,sep = ""))
    p <- p + ylab(paste("UMAP",y.dim,sep = ""))
    p <- p + guides(colour=guide_legend(title=NULL))
    p <- p + scale_color_manual(values = color_manual)
    
    p <- p + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    return(p)
  }



MyPairs <- function(data,
                    ...) {
  # check replication
  MyPoints <- function(x, 
                       y){
    points(x, 
           y, 
           pch = 19, 
           cex = .2)
  }
  
  # calculate correlation
  ViewCor <- function(x, 
                      y,
                      ...){  
    cor.res <- cor(x, y, method = "pearson")
    cor.res <- round(cor.res, 
                     3)
    
    text.x <- max(x) /15
    text.y <- max(y) ^ 0.5
    
    output.text <- cor.res
    #cex.cor <- 0.8/strwidth(output.text)    
    text(text.x, 
         text.y, 
         output.text, 
         cex = 8,
         ... )
  }
  
  pairs(data, 
        log = "xy",
        lower.panel = MyPoints, 
        upper.panel = ViewCor,
        cex.labels=5.4,
        ...)
}


####by WX principlecurve######
MyPrincipalCurve <- 
function(coord,smoother = "smooth.spline",
         class=NULL,classorder=NULL,plot.true = T,...){
  library(princurve)
  if((!is.null(class))&(!is.null(classorder))){
    coord.list=list()
    for (i in 1:length(classorder)) {
      coord.list[[i]]=coord[which(class==classorder[i]),]
    }
    gravity_center.list=lapply(coord.list,colMeans)
    gravity_center=do.call("rbind",gravity_center.list)
    curve=principal_curve(x = coord,smoother = smoother,plot_iterations = plot.true,
                          start = gravity_center,...)
  }else{
    curve=principal_curve(x = coord,smoother = smoother,plot_iterations = plot.true,...)
  }
  return(curve)
}

######
MyLogBoxPlot <- function(data,
                       min.cutoff = min.cutoff,
                       log = T,
                       cex.axis = 2,cex.lab = 2,cex.main = 2,
                       ...){
  data.na <- data
  data.na[data.na < min.cutoff] <- NA
  if(log ==  T){
    boxplot(log2(data.na),
            outline = F,
            ...)
  }
  else{
    boxplot(data.na,
            outline = F,
            ...)
  }
}


###
# XuLab R function
###hs###
MyLimmaPlot <- function(de.data,
                        adj.P.Val.cutoff = 0.05,
                        fc.max,
                        fc.cutoff = 1,
                        xlim,
                        up.color = "#F25362",
                        down.color = "#487CC2",
                        undiff.color = "black",
                        xlab = "mean of normalized counts",
                        ylab = "log2 fold change",
                        count = T,
                        ...){
  de.data <- as.data.frame(de.data)
  de.color <- rep(undiff.color, nrow(de.data))
  de.color[(!is.na(de.data$adj.P.Val)) & 
             (de.data$adj.P.Val < adj.P.Val.cutoff) & 
             (de.data$logFC < -fc.cutoff)] <- down.color
  de.color[(!is.na(de.data$adj.P.Val)) & 
             (de.data$adj.P.Val < adj.P.Val.cutoff) & 
             (de.data$logFC > fc.cutoff)] <- up.color
  up.count <- sum((!is.na(de.data$adj.P.Val)) & 
                    (de.data$adj.P.Val < adj.P.Val.cutoff) & 
                    (de.data$logFC > fc.cutoff))
  down.count <- sum((!is.na(de.data$adj.P.Val)) & 
                      (de.data$adj.P.Val < adj.P.Val.cutoff) & 
                      (de.data$logFC < -fc.cutoff))
  shape.plot <- rep(16,length(de.data$logFC))
  shape.plot[abs(de.data$logFC) > fc.max] <- 17
  de.data$logFC[de.data$logFC > fc.max] <- fc.max
  de.data$logFC[de.data$logFC < -fc.max] <- -fc.max
  plot(de.data$AveExpr,
       de.data$logFC,
       xlab = xlab,
       ylab = ylab,
       ylim = c(-fc.max, fc.max),
       xlim = xlim,
       col = de.color,
       #log = "x",
       pch = shape.plot,
       cex = 1,
       ...)
  box(lwd=2)
  if(count){
    text(1,
         fc.max / 2,
         as.character(up.count),
         col = up.color,
         cex = 1)
    text(1,
         -(fc.max / 2),
         as.character(down.count),
         col = down.color,
         cex = 1)
  }
}
MyLimmaPlot2 <- function(de.data,
                         P.Val.cutoff = 0.05,
                         fc.max,
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
  de.data$state <- '-'
  de.data[(!is.na(de.data$P.Value)) & 
            (de.data$P.Value <= P.Val.cutoff) & 
            (de.data$logFC <= -fc.cutoff),'state'] <- 'down'
  de.data[(!is.na(de.data$P.Value)) & 
            (de.data$P.Value <= P.Val.cutoff) & 
            (de.data$logFC >= fc.cutoff),'state'] <- 'up'
  de.data$state <- factor(de.data$state,levels = c('-','up','down'))
  de.data <- de.data[order(de.data$state),]
  up.count <- sum(de.data$state=='up')
  down.count <- sum(de.data$state=='down')
  # shape.plot <- rep(16,length(de.data$logFC))
  # shape.plot[abs(de.data$logFC) >= fc.max] <- 17
  de.data$logFC[de.data$logFC >= fc.max] <- fc.max
  de.data$logFC[de.data$logFC <= -fc.max] <- -fc.max
  
  plot(de.data$AveExpr,
       de.data$logFC,
       xlab = xlab,
       ylab = ylab,
       ylim = c(-fc.max, fc.max),
       xlim = xlim,
       col = de.color[de.data$state],
       #log = "x",
       pch = 16,
       cex = 1,
       ...)
  box(lwd=box.lwd)
  if(count){
    text(1,
         fc.max / 2,
         as.character(up.count),
         #col = up.color,
         cex = text.cex)
    text(1,
         -(fc.max / 2),
         as.character(down.count),
         #col = down.color,
         cex = text.cex)
  }
  return(de.data)
}


MyLimmaPlot3 <- function(de.data,
                         P.Val.cutoff = 0.05,
                         fc.max,
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
  de.data$logFC[de.data$logFC >= fc.max] <- fc.max
  de.data$logFC[de.data$logFC <= -fc.max] <- -fc.max
  
  plot(de.data$AveExpr,
       de.data$logFC,
       xlab = xlab,
       ylab = ylab,
       ylim = c(-fc.max, fc.max),
       xlim = xlim,
       col = de.color[de.data$state],
       #log = "x",
       pch = 16,
       cex = 1,
       ...)
  box(lwd=box.lwd)
  if(count){
    text(1,
         fc.max / 2,
         as.character(up.count),
         #col = up.color,
         cex = text.cex)
    text(1,
         -(fc.max / 2),
         as.character(down.count),
         #col = down.color,
         cex = text.cex)
  }
}

MybedViolinMed <-  function(bed,
                            samples,
                            color.violin = "gray60",
                            title.plot = "",
                            xlab.plot = "",
                            ylab.plot = "",
                            box.lwd = 4,
                          
                            cex.main = 3,
                            cex.axis = 1.5,
                            log = TRUE,
                            ...){
  library(vioplot)
  y <- as.data.frame(bed)[,samples]
  if(log){
    y <- log2(y+0.1)
  }
  
  plot(1:2,1:2,col = "white",
       xlab = xlab.plot,
       ylab = ylab.plot,
       xlim = c(0.5,length(samples)+0.5),
       ylim = c(min(y),max(y)),
       xaxt="n",
       main = title.plot,
       cex.main = cex.main,
       ...)
  box(lwd = box.lwd)
  axis(1,lwd = 0,lwd.ticks = 0, at = 1:length(samples),cex.axis = cex.axis,labels = samples)
  for (n.plot in 1:length(samples)){
    vioplot(y[,samples[n.plot]],
            at = n.plot,
            drawRect = T,
            rectCol = "white",
            colMed = "black",
            # pchMed = 1,
            col = color.violin,
            add = T,
            ...)
  }
}

MyfcDataplot2 <- 
function(bed,
         phase1,
         phase2,
         fc.state,
         fc.value.up,
         fc.value.down,
         fc.cutoff,
         cutoff.K,
         xlim,
         ylim,
         lwd=5
){
  mcols(bed)[,fc.state] <- "-"
  mcols(bed)[,fc.state][mcols(bed)[,phase1] / mcols(bed)[,phase2] < (1/fc.cutoff) &
                          mcols(bed)[,phase2] > cutoff.K] <- "up"
  mcols(bed)[,fc.state][mcols(bed)[,phase2] / mcols(bed)[,phase1] < (1/fc.cutoff) &
                          mcols(bed)[,phase1] > cutoff.K] <- "down"
  mcols(bed)[,fc.state] <- factor(mcols(bed)[,fc.state],
                                  levels = c("-",
                                             "down",
                                             "up"))
  mcols(bed)[,fc.value.up] <- ifelse(mcols(bed)[,phase2] > cutoff.K,
                                     mcols(bed)[,phase1] / mcols(bed)[,phase2],
                                     "-")
  mcols(bed)[,fc.value.down] <- ifelse(mcols(bed)[,phase1] > cutoff.K,
                                       mcols(bed)[,phase2] / mcols(bed)[,phase1],
                                       "-")
  plot(log2(mcols(bed)[,phase1]),
       log2(mcols(bed)[,phase2]),
       xlim = xlim,
       ylim = ylim,
       col = color.fc[mcols(bed)[,fc.state]],
       pch = 16,
       cex = 0.5)
  box(lwd = lwd)
  MyText(paste(fc.state,"\n","up:",sum(mcols(bed)[,fc.state] == "up"),"\n",
               "down:",sum(mcols(bed)[,fc.state] == "down"),"\n",
               "-:",sum(mcols(bed)[,fc.state] == "-"),"\n"))
  return(bed)
}
########WX######
#limma pipeline
Mylimma <- function(exp, 
                 condition1.samples, 
                 condition2.samples){
  library(limma)
  exp=exp[,c(condition1.samples,condition2.samples)]
  group_list=c(rep("condition1",length(condition1.samples)),
               rep("condition2",length(condition2.samples)))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exp)
  contrast.matrix=makeContrasts("condition2-condition1",levels = design)
  fit=lmFit(exp,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  tempOutput = na.omit(tempOutput)
  tempOutput
}

##########
MyPlotTrend <- function(
  tpm.data,
  time.data,
  col = "black",
  line.col = "black",
  plot.line=F,
  line.width=5,
  point.width=3,
  box.col,
  ...) {
  tpm.data.max <- tpm.data / apply(tpm.data,1,quantile,probs = 0.99)
  
  xpt <- time.data
  ypt <- colMeans(tpm.data.max)
  
  y.loess <-loess(ypt~xpt)
  ypt.predict <- predict(y.loess,data.frame(x=xpt),se=T)
  # trend.merge[[tree.n]] <- ypt.predict$fit
  trend.merge <- ypt.predict$fit
  
  plot(xpt,
       ypt.predict$fit, 
       lwd = line.width,pch=19,
       cex=point.width,
       col = col,
       bty = "n",
       xaxt = "n",
       xlab = "",
       ylab = "",
       ... )
  if(plot.line){  lines(xpt,
                   ypt.predict$fit,
                   lwd = line.width,
                   col = line.col,
                   ...)}

  box(lwd=4,
      col = box.col)
  return(trend.merge)
}


Myfilter.outlier <- function(data,
                             top.quantile = 1,
                             down.quantile = "NULL",
                             return.mean = TRUE
){
  data <- as.numeric(data)
  if(length(unique(data)) > 1 ){
    if(top.quantile != "NULL"){
      top.cutoff <- quantile(data,top.quantile)
      filter.top <- sum(data >= top.cutoff)
      print(paste("filter.top",filter.top,sep = ":"))
      if(down.quantile != "NULL"){
        down.cutoff <- quantile(data,down.quantile)
        filter.down <- sum(data <= down.cutoff)
        print(paste("filter.down",filter.down,sep = ":"))
        data.filter <- data[data < top.cutoff & data > down.cutoff]
      }
      else if(down.quantile == "NULL"){
        data.filter <- data[data < top.cutoff]
      }
    }
    if(top.quantile == "NULL"){
      if(down.quantile != "NULL"){
        down.cutoff <- quantile(data,down.quantile)
        filter.down <- sum(data <= down.cutoff)
        print(paste("filter.down",filter.down,sep = ":"))
        data.filter <- data[data > down.cutoff]
      }
      else if(down.quantile == "NULL"){
        print("Hey, what are you doing?")
        print("Please input at least one quantile")
      }
    }
  }
  else{
    data.filter <- data
    print(paste("filter.top",0,sep = ":"))
    print(paste("filter.down",0, sep= ":"))
  }
  if(return.mean == TRUE){
    data.mean <- mean(data.filter)
    return(data.mean)
  }
  else {return(data.filter)}
}

Myfilter.outlier2 <- function(data,
                              top.numer = 5,
                              down.number = 5,
                              return.mean = TRUE
){
  data <- as.numeric(data)
  data.order <- sort(data,decreasing = T)
  data.filter <- data.order[(top.numer+1):(length(data)-down.number)]
  if(return.mean == TRUE){
    data.mean <- mean(data.filter)
    return(data.mean)
  }
  else {return(data.filter)}
}

MyPlotAsBitmap <- function(command,main="",tempfile="temp.png",width=500,height = 700){
  library(grid)
  library(gridExtra)
  library(png)
  png(tempfile,width,height)
  eval(parse(text = command))
  dev.off()
  raster=rasterGrob(readPNG(tempfile),
                    interpolate = T)
  plot.new()
  text(0,1,labels = main)
  grid.arrange(raster,newpage = F)
}

MySeuratTsne10x2Gg <- function(seurat,
                               seurat.meta,
                               x.dim = 1,
                               y.dim = 2,
                               color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- seurat@dr$tsne@cell.embeddings
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   seurat.meta[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName
              ))
  p <- p + xlab(paste("tSNE",x.dim,sep = ""))
  p <- p + ylab(paste("tSNE",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}

MySeuratUMAP10x2Gg <- function(seurat,
                               seurat.meta,
                               x.dim = 1,
                               y.dim = 2,
                               color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- seurat@dr$umap@cell.embeddings
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   seurat.meta[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName
              ))
  p <- p + xlab(paste("UMAP",x.dim,sep = ""))
  p <- p + ylab(paste("UMAP",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}


Mydif2Gg <- function(dif,
                     sample.inf,
                     x.dim = 1,
                     y.dim = 2,
                     color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  dif <- as.data.frame(dif)
  dif.df <- cbind(x.pos = dif[,x.dim],
                  y.pos = dif[,y.dim],
                  sample.inf[row.names(dif),])
  
  p <- ggplot(data = dif.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("DC",x.dim,sep = ""))
  p <- p + ylab(paste("DC",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)+ 
    #scale_y_reverse() +
    scale_color_manual(values = time.colors) +
    scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2))
  
  return(p)
}
MyCor <- function(exp.data,
                  var.gene = NULL,
                  exp.cutoff = 1,
                  exp.prop.whole.max = 0.9,
                  exp.prop.whole.min = 0,
                  vector.group = NULL,
                  exp.prop.group.min = 0,
                  cor.method = c("pearson","rho","phs"),
                  cor.cutoff = 0.3,
                  partner.cutoff = 5) {
  cor.method <- cor.method[1]
  if(!is.null(var.gene)){
    exp.data <- exp.data[rownames(exp.data) %in% var.gene,]
  }
  exp.data <- exp.data[rowSums(exp.data >= exp.cutoff) <= exp.prop.whole.max * ncol(exp.data) &
                         rowSums(exp.data >= exp.cutoff) >= exp.prop.whole.min * ncol(exp.data),]
  if(!is.null(vector.group)){
    gene.group <- c()
    for(group.tmp in unique(vector.group)){
      exp.data.tmp <- exp.data[,vector.group %in% exp.data.tmp]
      gene.group.tmp <- rownames(exp.data.tmp)[rowSums(exp.data.tmp >= exp.cutoff) >= exp.prop.group.min * ncol(exp.data.tmp)]
      gene.group <- c(gene.group,
                      gene.group.tmp)
    }
    exp.data <- exp.data[rownames(exp.data) %in% gene.group,]
  }
  
  if(cor.method == "pearson"){
    cor.data <- cor(t(exp.data),method = "pearson")
  } else if (cor.method == "rho") {
    cor.data <- propr:::lr2rho(t(exp.data))
    rownames(cor.data) <- rownames(exp.data)
    colnames(cor.data) <- rownames(exp.data)
  } else if (cor.method == "phs") {
    cor.data <- propr:::lr2phs(t(exp.data))
    rownames(cor.data) <- rownames(exp.data)
    colnames(cor.data) <- rownames(exp.data)
  }
  
  diag(cor.data) <- 0
  gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
  while(nrow(cor.data) > length(gene.cor)) {
    # cat(nrow(cor.data)," ",length(gene.cor),"\n")
    cor.data <- cor.data[gene.cor,gene.cor]
    gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
  }
  return(cor.data)
}
#######form WX#####
MyPlotAsBitmap=function(command,main="",tempfile="temp.png",size=500){
  library(grid)
  library(gridExtra)
  library(png)
  png(tempfile,size,size)
  eval(parse(text = command))
  dev.off()
  raster=rasterGrob(readPNG(tempfile),
                    interpolate = T)
  plot.new()
  text(0,1,labels = main)
  grid.arrange(raster,newpage = F)
}

MyGSEA.prepare=function(tpm,
                        sample.class,
                        gct.name,
                        cls.name,
                        mode="discrete",
                        continuous_value_name=NULL){
  print("use EnsemblGeneID")
  gct=cbind(NAME=rownames(tpm),
            Description=NA,
            tpm)
  
  gct.file=file(gct.name, open = "wt")
  writeLines("#1.2",gct.file)
  writeLines(paste(nrow(tpm),ncol(tpm),sep = "\t"),gct.file)
  write.table(gct,gct.file,append = T,quote = F,row.names = F,sep = "\t")
  close(gct.file)
  
  if(mode=="discrete"){
    cls.file=file(cls.name, open = "wt")
    writeLines(paste(ncol(tpm),
                     length(table(as.character(sample.class))),1,
                     sep = " "),
               cls.file)
    writeLines(paste(c("#",names(table(as.character(sample.class)))),
                     collapse = " "),
               cls.file)
    writeLines(paste(sample.class,collapse = " "),cls.file)
    close(cls.file)
  }else if(mode=="continuous"){
    cls.file=file(cls.name, open = "wt")
    writeLines("#numeric",
               cls.file)
    if(is.null(continuous_value_name)){
      writeLines("#Trajectory",
                 cls.file)
    }else{
      writeLines(paste("#",continuous_value_name,sep = ""),
                 cls.file)
    }
    writeLines(paste(sample.class,collapse = " "),cls.file)
    close(cls.file)
  }
}
MyGSEA.run=function(gct,cls,gmt,mode="discrete"){
  if(mode=="discrete"){
    mycommand=paste("-cp F:\\gsea-3.0.jar -Xmx4096m xtools.gsea.Gsea -res",gct,
                    "-cls", cls,
                    "-gmx", gmt,
                    "-collapse false", sep=" ")
    system2(command='java', args=mycommand)
  }else if(mode=="continuous"){
    print("mode:continuous")
    mycommand=paste("-cp F:\\gsea-3.0.jar -Xmx4096m xtools.gsea.Gsea -res",gct,
                    "-cls", cls,
                    "-gmx", gmt,
                    "-collapse false",
                    "-metric Pearson", sep=" ")
    print(mycommand)
    system2(command='java', args=mycommand)
  }
}



MyordergenewithPseudotime <- function(ordered.exp,genelist,graph=F,return.peak.pos=F){
  ordered.exp=ordered.exp[genelist,]
  peak.pos=c()
  for (gene in genelist) {
    state=unlist(ordered.exp[gene,])>quantile(unlist(ordered.exp[gene,]),0.8)
    peak.pos=c(peak.pos,mean(which(state)))
  }
  names(peak.pos)=genelist
  if(graph){
    hist(peak.pos,breaks = round(ncol(ordered.exp)/10,0),xlim = c(0,ncol(ordered.exp)))
  }
  if(return.peak.pos==F){
    return(names(sort(peak.pos)))
  }else{
    return(peak.pos)
  }
}
MyCo <- 
function(exp.data,
         var.gene = NULL,
         exp.cutoff = 1,
         exp.prop.whole.min = 0,
         exp.prop.whole.max = 2,
         vector.group = NULL,
         exp.prop.group.min = 0,
         exp.prop.group.max = 2,
         cor.method = c("pearson","rho","phs",'spearman','Baco'),
         cor.cutoff = 0,
         partner.cutoff = 0,
         refine.cor = F) {
  cor.method <- cor.method[1]
  if(!is.null(var.gene)){
    exp.data <- exp.data[rownames(exp.data) %in% var.gene,]
  }
  exp.data <- exp.data[rowSums(exp.data >= exp.cutoff) <= exp.prop.whole.max * ncol(exp.data) &
                         rowSums(exp.data >= exp.cutoff) >= exp.prop.whole.min * ncol(exp.data),]
  if(!is.null(vector.group)){
    gene.group.min <- c()
    for(group.tmp in unique(vector.group)){
      exp.data.tmp <- exp.data[,vector.group %in% group.tmp]
      gene.group.min.tmp <- rownames(exp.data.tmp)[rowSums(exp.data.tmp >= exp.cutoff) >= exp.prop.group.min * ncol(exp.data.tmp)]
      gene.group.min <- c(gene.group.min,
                          gene.group.min.tmp)
    }
    group.tmp = unique(vector.group)[1]
    exp.data.tmp <- exp.data[,vector.group %in% group.tmp]
    gene.group.max <- rownames(exp.data.tmp)[rowSums(exp.data.tmp >= exp.cutoff) >= exp.prop.group.max * ncol(exp.data.tmp)]
    for(group.tmp in unique(vector.group)[-1]){
      exp.data.tmp <- exp.data[,vector.group %in% group.tmp]
      gene.group.max.tmp <- rownames(exp.data.tmp)[rowSums(exp.data.tmp >= exp.cutoff) >= exp.prop.group.max * ncol(exp.data.tmp)]
      gene.group.max <- intersect(gene.group.max,
                                  gene.group.max.tmp)
    }
    
    exp.data <- exp.data[rownames(exp.data) %in% gene.group.min,]
    exp.data <- exp.data[!rownames(exp.data) %in% gene.group.max,]
  }
  
  if(cor.method == "pearson"){
    cor.data <- cor(t(exp.data),method = "pearson")
  } else if (cor.method == "rho") {
    cor.data <- propr:::lr2rho(t(exp.data))
    rownames(cor.data) <- rownames(exp.data)
    colnames(cor.data) <- rownames(exp.data)
  } else if (cor.method == "phs") {
    cor.data <- propr:::lr2phs(t(exp.data))
    rownames(cor.data) <- rownames(exp.data)
    colnames(cor.data) <- rownames(exp.data)
  }else if(cor.method == "Baco"){
    cor.data <- BaCo(t(exp.data))
    rownames(cor.data) <- rownames(exp.data)
    colnames(cor.data) <- rownames(exp.data)
  }else if (cor.method == "spearman") {
    cor.data <- cor(t(exp.data),method = "spearman")}
  
  diag(cor.data) <- 0
  gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
  while(nrow(cor.data) > length(gene.cor)) {
    # cat(nrow(cor.data)," ",length(gene.cor),"\n")
    cor.data <- cor.data[gene.cor,gene.cor]
    gene.cor <- rownames(cor.data)[rowSums(cor.data > cor.cutoff) >= partner.cutoff]
  }
  if(refine.cor){
    cor.data[cor.data < cor.cutoff] <- 0
  }
  return(cor.data)
}


## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#######
MyPrincipalCurve <- 
function(coord,smoother = "smooth.spline",
         class=NULL,classorder=NULL,plot.true = T,...){
  library(princurve)
  if((!is.null(class))&(!is.null(classorder))){
    coord.list=list()
    for (i in 1:length(classorder)) {
      coord.list[[i]]=coord[which(class==classorder[i]),]
    }
    gravity_center.list=lapply(coord.list,colMeans)
    gravity_center=do.call("rbind",gravity_center.list)
    curve=principal_curve(x = coord,smoother = smoother,plot_iterations = plot.true,
                          start = gravity_center,...)
  }else{
    curve=principal_curve(x = coord,smoother = smoother,plot_iterations = plot.true,...)
  }
  return(curve)
}
MyNetOrdGene <- 
function(graph.data,
         gene.cluster,
         classorder = NULL,
         class.start = NULL,
         smoother = "periodic_lowess",
         maxit = 10,
         plot.true = F,
         rev_order = F,
         seed.n = 1){
  set.seed(seed.n)
  gene.prin_curv <- MyPrincipalCurve(layout_with_fr(graph.data),
                                     smoother = smoother,
                                     class = gene.cluster,
                                     classorder = classorder,
                                     plot.true = plot.true,
                                     maxit = maxit) 
  
  gene.count.mid <- ceiling(length(gene.cluster)/2)
  gene.sort.1 <- names(gene.cluster)[order(gene.prin_curv$lambda,decreasing = rev_order)]
  gene.sort.2 <- c(gene.sort.1[gene.count.mid:length(gene.sort.1)],gene.sort.1[1:gene.count.mid-1])
  
  class_in_sort1 <- c()
  class_in_sort2 <- c()
  for (class.tmp in unique(gene.cluster)) {
    gene.tmp <- names(gene.cluster)[gene.cluster == class.tmp]
    if (length(gene.tmp) == 1 | sd(which(gene.sort.1 %in%  gene.tmp)) < sd(which(gene.sort.2 %in%  gene.tmp))){
      class_in_sort1 <- c(class_in_sort1,class.tmp)
    } else {
      class_in_sort2 <- c(class_in_sort2,class.tmp)
    }
  }
  
  
  if(is.null(classorder)){
    class.sort.1 <- gene.cluster[order(gene.prin_curv$lambda,decreasing = rev_order)]
    class.sort.2 <- c(class.sort.1[gene.count.mid:length(class.sort.1)],class.sort.1[1:gene.count.mid-1])
    
    class.sort.1 <- class.sort.1[class.sort.1 %in% class_in_sort1]
    class.sort.2 <- class.sort.2[class.sort.2 %in% class_in_sort2]
    
    class.sort.1.mean <- c()
    class.sort.2.mean <- c()
    for (class.tmp in class_in_sort1){
      class.sort.1.mean <- c(class.sort.1.mean,
                             mean(which(class.sort.1 == class.tmp)))
    }
    for (class.tmp in class_in_sort2){
      class.sort.2.mean <- c(class.sort.2.mean,
                             mean(which(class.sort.2 == class.tmp)))
    }
    classorder <- c(class_in_sort1[order(class.sort.1.mean)],
                    class_in_sort2[order(class.sort.2.mean)])
  }
  
  if(!is.null(class.start)){
    classorder <- c(classorder[which(classorder == class.start):length(classorder)],
                    classorder[1:which(classorder == class.start)-1])
  }
  
  cat("classorder:",classorder,"\n")
  gene.sort.out <- c()
  for (class.tmp in classorder) {
    gene.tmp <- names(gene.cluster)[gene.cluster == class.tmp]
    if (class.tmp %in% class_in_sort1){
      gene.tmp.sort <- gene.sort.1[gene.sort.1 %in% gene.tmp]
    } else if (class.tmp %in% class_in_sort2) {
      gene.tmp.sort <- gene.sort.2[gene.sort.2 %in% gene.tmp]
    }
    gene.sort.out <- c(gene.sort.out,
                       gene.tmp.sort)
  }
  return(gene.sort.out)
}



Mymonocle2gg <- function(tsne.gg,
                         #sample.inf,
                         x.dim = 1,
                         y.dim = 2,
                         color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  p <- tsne.gg + xlab(paste("tSNE",x.dim,sep = ""))
  p <- p + ylab(paste("tSNE",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}
# LDA result to ggplot2
MyLda2Gg <- function(lda,
                     tpm,
                     sample.inf,
                     x.dim = 1,
                     y.dim = 2,
                     color_manual = brewer.pal(9,"Set1"),scale=T){
  library(ggplot2)
  lda.percent <- prop.table(lda$svd^2)
  if(scale==T){
    LDAprojection <- scale(t(as.matrix(tpm)),scale=FALSE) %*% lda$scaling
  }else{LDAprojection <- t(as.matrix(tpm))%*% lda$scaling}
  row.names(sample.inf) <- sample.inf$SampleName
  lda.df <- cbind(x.pos = LDAprojection[,x.dim],
                  y.pos = LDAprojection[,y.dim],
                  sample.inf[row.names(LDAprojection),])
  
  p <- ggplot(data = lda.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos#,
                           # label = SampleName
                            ))
  p <- p + xlab(paste("LD",x.dim,"(",round(lda.percent[x.dim]*100,1),"%)",sep = ""))
  p <- p + ylab(paste("LD",y.dim,"(",round(lda.percent[y.dim]*100,1),"%)",sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)+scale_color_manual(values = c(time.colors,"gray30")) +
    #scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1,
                                 ncol = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2,
                                ncol = 1))
  return(p)
}

#run LDA
Mylda <- function(tpm,
                  class,
                  plot = T
){
  library(MASS)
  exp <- as.data.frame(t(as.matrix(tpm)))
  exp <- cbind(exp,type=class)
  LDA <- lda(type~.,data=exp)
  if(plot){
    plot(LDA)
    print(summary(LDA))
  }
  return(LDA)
}
# do Linear model corelation analysis
MyLM <- function(tpm,
              variable,
              direction="both"){
  tpm2variable=as.data.frame(cbind(t(tpm),variable=variable[colnames(tpm)]))
  p.val=c()
  if (direction=="both") {
    for (gene in colnames(tpm2variable)[-ncol(tpm2variable)]) {
      LM=lm(data=tpm2variable,variable~get(gene))
      if(nrow(summary(LM)$coefficients)==2){
        p.val=c(p.val,summary(LM)$coefficients[2,4])}
      else{p.val=c(p.val,NA)}
    }
  }else if (direction=="pos") {
    for (gene in colnames(tpm2variable)[-ncol(tpm2variable)]) {
      LM=lm(data=tpm2variable,variable~get(gene))
      if(nrow(summary(LM)$coefficients)==2){
        if(summary(LM)$coefficients[2,1]>0){
          p.val=c(p.val,summary(LM)$coefficients[2,4])
        }else{p.val=c(p.val,NA)}
      }
      else{p.val=c(p.val,NA)}
    }
  }else if (direction=="neg") {
    for (gene in colnames(tpm2variable)[-ncol(tpm2variable)]) {
      LM=lm(data=tpm2variable,variable~get(gene))
      if(nrow(summary(LM)$coefficients)==2){
        if(summary(LM)$coefficients[2,1]<0){
          p.val=c(p.val,summary(LM)$coefficients[2,4])
        }else{p.val=c(p.val,NA)}
      }
      else{p.val=c(p.val,NA)}
    }
  }else{
    stop("direction = both/pos/neg")
  }
  names(p.val)=colnames(tpm2variable)[-ncol(tpm2variable)]
  p.val=sort(p.val)
  return(p.val)
}
#filter gene with LDA
MyGeneLda <- function(lda,
                      tpm,
                      dim,
                      p.val.cut=0.001,
                      direction = "both"){
  lda.scale <- lda$scaling
  LDAprojection <- scale(t(as.matrix(tpm)),scale=FALSE) %*% lda.scale
  LD.dim <- LDAprojection[,dim]
  p.val <- MyLM(tpm,LD.dim)
  select.gene <- names(p.val[p.val<p.val.cut])
  
  lda.res <- lda.scale[,dim]
  names(lda.res) <- rownames(lda.scale)
  lda.res <- lda.res[select.gene]
  if(direction == "both"){
    filter.gene <- select.gene
  }
  else if(direction == "pos"){
    filter.gene <- names(lda.res[lda.res>0])
  }
  else if(direction == "neg"){
    filter.gene <- names(lda.res[lda.res<0])
  }
  else{
    stop("direction = both/pos/neg")
  }
  filter.gene <- unique(filter.gene)
  return(filter.gene)
}
#######

MySwne2Gg <- function(swne.embedding,
                      sample.inf,
                      x.dim = 1,
                      y.dim = 2,
                      color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  data.coord <- swne.embedding$sample.coords
  data.coord <- as.data.frame(data.coord)
  
  row.names(sample.inf) <- sample.inf$SampleName
  data.df <- cbind(x.pos = data.coord[,x.dim],
                   y.pos = data.coord[,y.dim],
                   sample.inf[row.names(data.coord),])
  
  p <- ggplot(data = data.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("SWNE_",x.dim,sep = ""))
  p <- p + ylab(paste("SWNE_",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}

MyPlotAsBitmap <- function(command,main="",tempfile="temp.png",size=500){
  library(grid)
  library(gridExtra)
  library(png)
  png(tempfile,size,size)
  eval(parse(text = command))
  dev.off()
  raster <- rasterGrob(readPNG(tempfile),
                    interpolate = T)
  plot.new()
  text(0,1,labels = main)
  grid.arrange(raster,newpage = F)
}

Myfilter <- function(exp,
                     gene,
                     pearson.threshold=0.5,
                     partner.threshold=5,
                     expression.threshold=1,
                     top.dispersion.interval=0.8,
                     bottom.dispersion.interval=0){
  pearson <- cor(t(exp)[,gene],method = "pearson")
  # hist(pearson[upper.tri(pearson)],breaks = 100)
  gene.filtered <- rownames(pearson)[rowSums(pearson>pearson.threshold)>=partner.threshold]
  gene.filtered <- gene.filtered[rowSums(exp[gene.filtered,]>=expression.threshold)<=top.dispersion.interval*ncol(exp)&
                                rowSums(exp[gene.filtered,]>=expression.threshold)>=bottom.dispersion.interval*ncol(exp)]
  return(gene.filtered)
}
############
Myscexp <- function(xpt,
         ypt.matrix,
         xlab = "Pseudotime",
         ylab = "Relative Expression",
         col.line = "gray60",
         col.mainline = "red",
         ...){
  ypt.matrix <- ypt.matrix[,order(xpt)]
  xpt <- sort(xpt)
  
  ypt.matrix.relat <- ypt.matrix / rowMax(as.matrix(ypt.matrix))
  
  # pdf("bc.hl.pdf")
  plot(range(xpt),
       c(0,1),
       col = "white",
       bty = "n",
       xaxt = "n",
       xlab = xlab,
       ylab = ylab,
       ...)
  
  for(n.plot in 1:nrow(ypt.matrix.relat)){
    ypt <- as.numeric(ypt.matrix.relat[n.plot,])
    y.loess <-loess(ypt~xpt)
    ypt.predict <- predict(y.loess,
                           data.frame(x=xpt),
                           se=T)
    lines(xpt,
          ypt.predict$fit,
          lwd = 1,
          col = col.line)
  }
  
  ypt <- colMeans(ypt.matrix.relat)
  y.loess <-loess(ypt~xpt)
  ypt.predict <- predict(y.loess,data.frame(x=xpt),se=T)
  lines(xpt,
        ypt.predict$fit, 
        lwd = 8,
        col = col.mainline)
}

# print environment information 
sessionInfo()

# import library
library(RColorBrewer)
library(Biobase)

######myheatmap2#####
myHeatmap2 <- function(x,
         #ord,
         xlab="",
         ylab="",
         main="My Heatmap",
         col=atac.state.col,
         cutoff.min=-1,
         cutoff.max=1,
         ord = TRUE,
         labCol = "",
         cex = 1.5,
         sort = TRUE,
         ...){
  op <- par(mar=c(5,2,2,2)+0.1)
  on.exit(par(op))
  nc <- NCOL(x)
  nr <- NROW(x)
  if(ord == TRUE){
    for (i in 1:ncol(x)){
    x[,i] <- x[,i][order(x[,i],decreasing = sort)]
    }
  }
  x <- t(x)
  x[x< cutoff.min] <- cutoff.min
  x[x > cutoff.max] <- cutoff.max
  image(1L:nc, 1L:nr, 
        x, 
        xlim = 0.5 + c(0, nc),
        ylim = 0.5 +c(0, nr),
        axes = FALSE, 
        xlab=xlab,
        ylab=ylab,
        main=main,
        col=col,
        zlim = c(cutoff.min,cutoff.max),
        ...
  )
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0,cex.axis = cex)
  axis(2, 1L:nr, labels = NA, las = 2, line = -0.5, tick = 0)
  return(x)
}

##########boxplot###
Myboxplot5 <- function(tpm.input,
         phase,
         label,
         main_text='',
         side="two.side",
         color=box.col,
         notch = T,
         paired = T,
         cex.axis = 1.5,
         cex.main  = 2,
         #ylim,
         ...
){
  phase1 <- phase[1]
  phase2 <- phase[2]
  phase3 <- phase[3]
  phase4 <- phase[4]
  phase5 <- phase[5]
  boxplot(tpm.input[,phase],
          col=color,
          #ylim=ylim,
          outline = F,
          notch = notch,
          main=main_text,
          xaxt="n",
          cex.axis = 1.3,
          cex.main = cex.main,
          ...
  )
  axis(1,at = 1:length(phase),labels = label,cex.axis = cex.axis)
  box(lwd=3)
  wilcox.test.res1 <- wilcox.test(tpm.input[,phase1],
                                  tpm.input[,phase2],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res1$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res2 <- wilcox.test(tpm.input[,phase2],
                                  tpm.input[,phase3],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res2$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res3 <- wilcox.test(tpm.input[,phase3],
                                  tpm.input[,phase4],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res3$p.value,"\n"),
         text.cex = 1)
  
  wilcox.test.res4 <- wilcox.test(tpm.input[,phase4],
                                  tpm.input[,phase5],
                                  alternative = side,
                                  paired = paired)
  
  MyText(paste("p-value: \n",
               wilcox.test.res4$p.value,"\n"),
         text.cex = 1)
  
}

#######correation#####
##############?Խ???################
Myspearman.cor.test <- function(df,
                                method="spearman",
                                exact=NULL){
  num <- nrow(df)
  rho.df <- matrix(rep(0,num**2),
                   nrow=num,
                   ncol=num)
  colnames(rho.df) <- rownames(df)
  rownames(rho.df) <- rownames(df)
  
  p.value.df <- matrix(rep(1,num**2),
                       nrow=num,
                       ncol=num)
  colnames(p.value.df) <- rownames(df)
  rownames(p.value.df) <- rownames(df)
  
  for(i in seq(1,num)){
    for(j in seq(i,num)){
      tmp.cor <- cor.test(as.numeric(df[i,]),
                          as.numeric(df[j,]),
                          method=method,
                          exact=exact
      )
      rho.df[i,j] <- tmp.cor$estimate
      p.value.df[j,i] <- tmp.cor$p.value
    }
  }
  result <- list("rho"=rho.df,
                 "p.value"=p.value.df)
  return(result)
}
#########????######
Myspearman.cor.test2 <- function(df,
                                 method="spearman",
                                 exact=NULL){
  num <- nrow(df)
  rho.df <- matrix(rep(0,num**2),
                   nrow=num,
                   ncol=num)
  colnames(rho.df) <- rownames(df)
  rownames(rho.df) <- rownames(df)
  
  p.value.df <- matrix(rep(1,num**2),
                       nrow=num,
                       ncol=num)
  colnames(p.value.df) <- rownames(df)
  rownames(p.value.df) <- rownames(df)
  
  for(i in seq(1,num)){
    for(j in seq(1,num)){
      tmp.cor <- cor.test(as.numeric(df[i,]),
                          as.numeric(df[j,]),
                          method=method,
                          exact=exact
      )
      rho.df[i,j] <- tmp.cor$estimate
      p.value.df[i,j] <- tmp.cor$p.value
    }
  }
  result <- list("rho"=rho.df,
                 "p.value"=p.value.df)
  return(result)
}
Myspearman.cor.test3 <- function(df,
                                 method="pearson",
                                 exact=NULL){
  num <- ncol(df)
  rho.df <- matrix(rep(0,num**2),
                   nrow=num,
                   ncol=num)
  colnames(rho.df) <- colnames(df)
  rownames(rho.df) <- colnames(df)
  
  p.value.df <- matrix(rep(1,num**2),
                       nrow=num,
                       ncol=num)
  colnames(p.value.df) <- colnames(df)
  rownames(p.value.df) <- colnames(df)
  
  for(i in seq(1,num)){
    for(j in seq(1,num)){
      tmp.cor <- cor.test(as.numeric(df[,i]),
                          as.numeric(df[,j]),
                          method=method,
                          exact=exact
      )
      rho.df[i,j] <- tmp.cor$estimate
      p.value.df[i,j] <- tmp.cor$p.value
    }
  }
  result <- list("rho"=rho.df,
                 "p.value"=p.value.df)
  return(result)
}
#######hist2norm#######
histWithNormal <- function(d,xlab="",main="",breaks=10,
                           ylim=c(0,length(d)/2),xlim=c(min(d),max(d)),
                           includeMeanLine=T,col="lightgray") {
  g = d
  h<-hist(g, ylim=ylim,xlim=xlim,breaks=breaks, density=10, col=col, xlab=xlab, main=main) 
  xfit<-seq(min(g),max(g),length=100) 
  yfit<-dnorm(xfit,mean=mean(g),sd=sd(g)) 
  yfit <- yfit*diff(h$mids[1:2])*length(g) 
  lines(xfit, yfit, col="black", lwd=2)
  if(includeMeanLine) {
    m = mean(d,na.rm=T)
    segments(m,dnorm(m,mean=m,sd=sd(g))*diff(h$mids[1:2])*length(g),m,0,lwd=2,col="blue")
  }
}
###chip-seq fold change####
MyfcDataplot <- function(bed,
                         phase1,
                         phase2,
                         fc.state,
                         fc.cutoff,
                         cutoff.K,
                         xlim,
                         ylim
){
  mcols(bed)[,fc.state] <- "-"
  mcols(bed)[,fc.state][mcols(bed)[,phase1] / mcols(bed)[,phase2] < (1/fc.cutoff) &
                          mcols(bed)[,phase2] > cutoff.K] <- "up"
  mcols(bed)[,fc.state][mcols(bed)[,phase2] / mcols(bed)[,phase1] < (1/fc.cutoff) &
                          mcols(bed)[,phase1] > cutoff.K] <- "down"
  mcols(bed)[,fc.state] <- factor(mcols(bed)[,fc.state],
                                  levels = c("-",
                                             "down",
                                             "up"))
  plot(log2(mcols(bed)[,phase1]),
       log2(mcols(bed)[,phase2]),
       xlim = xlim,
       ylim = ylim,
       col = color.fc[mcols(bed)[,fc.state]],
       pch = 16,
       cex = 0.5)
  box(lwd = 8)
  MyText(paste(fc.state,"\n","up:",sum(mcols(bed)[,fc.state] == "up"),"\n",
               "down:",sum(mcols(bed)[,fc.state] == "down"),"\n",
               "-:",sum(mcols(bed)[,fc.state] == "-"),"\n"))
  return(bed)
}

####promoter state###
MyPromoterstateDataplot <- function(promoter.bed,
                                    
                                    H3K4me3,
                                    H3K27me3,
                                    statName,
                                    cutoff.K4,
                                    cutoff.K27,
                                    xlim,
                                    ylim,
                                    cex = 0.5,
                                    ...
){
  
  mcols(promoter.bed)[,statName] <- NA
  mcols(promoter.bed)[mcols(promoter.bed)[,H3K4me3] < cutoff.K4   &
                        mcols(promoter.bed)[,H3K27me3] < cutoff.K27,statName] <- "unmarked"
  mcols(promoter.bed)[mcols(promoter.bed)[,H3K4me3] >= cutoff.K4 &
                        mcols(promoter.bed)[,H3K27me3] < cutoff.K27,statName] <- "H3K4me3"
  mcols(promoter.bed)[mcols(promoter.bed)[,H3K4me3] < cutoff.K4 &
                        mcols(promoter.bed)[,H3K27me3] >= cutoff.K27,statName] <- "H3K27me3"
  mcols(promoter.bed)[mcols(promoter.bed)[,H3K4me3] >= cutoff.K4 &
                        mcols(promoter.bed)[,H3K27me3]>= cutoff.K27,statName] <- "bivalent"
  mcols(promoter.bed)[,statName] <- factor(mcols(promoter.bed)[,statName],
                                           levels = c("unmarked",
                                                      "H3K4me3",
                                                      "H3K27me3",
                                                      "bivalent"))
  
  plot(log2(mcols(promoter.bed)[,H3K4me3]),
       log2(mcols(promoter.bed)[,H3K27me3]),
       pch = 16,
       cex = cex,
       xlim = xlim,
       ylim = ylim,
       col = promoter.colors[mcols(promoter.bed)[,statName]],
  ...)
  box(lwd = 8)
  MyText(paste(H3K4me3,":",sum(mcols(promoter.bed)[,statName] == "H3K4me3"),"\n",
               H3K27me3,":",sum(mcols(promoter.bed)[,statName] == "H3K27me3"),"\n",
               "bivalency:",sum(mcols(promoter.bed)[,statName] == "bivalent"),"\n",
               "unmarked:",sum(mcols(promoter.bed)[,statName] == "unmarked")),
         text.cex = 2)
  return(promoter.bed)
}
##########enhancer state########
MystateDataplot <- function(bed,
                            k4me1,
                            k27ac,
                            state,
                            cutoff.k4,
                            cutoff.k27,
                            xlim,
                            ylim){
  mcols(bed)[,state] <- NA
  mcols(bed)[,state][mcols(bed)[,k4me1] < cutoff.k4] <- "unmarked"
  mcols(bed)[,state][mcols(bed)[,k4me1] >= cutoff.k4 &
                       mcols(bed)[,k27ac] < cutoff.k27] <- "primed"
  mcols(bed)[,state][mcols(bed)[,k4me1] >= cutoff.k4 &
                       mcols(bed)[,k27ac] >= cutoff.k27] <- "active"
  mcols(bed)[,state] <- factor(mcols(bed)[,state],
                               levels = c("unmarked",
                                          "primed",
                                          "active"))
  plot(log2(mcols(bed)[,k4me1]),
       log2(mcols(bed)[,k27ac]),
       xlim = xlim,
       ylim = ylim,
       col = enhancer.colors[mcols(bed)[,state]],
       pch = 16,
       cex = 0.5)
  box(lwd = 8)
  MyText(paste(state,"\n","active:",sum(mcols(bed)[,state] == "active"),"\n",
               "primed:",sum(mcols(bed)[,state] == "primed"),"\n",
               "unmarked:",sum(mcols(bed)[,state] == "unmarked"),"\n"))
  return(bed)
  
}
########remove redundancy of enriched GO terms####
MyslimGO <- function(enrichGO.result,
                     cutoff=0.7,
                     by="p.adjust",
                     select_fun=min,
                     plot=TRUE,
                     font=0.5,
                     n=50){
  require(clusterProfiler)
  bp2 <- simplify(enrichGO.result,
                  cutoff=cutoff,
                  by=by, 
                  select_fun=select_fun)
  if(plot=="TRUE"){emapplot(bp2,
                             vertex.label.font = font,
                             n=n)}
  return(bp2)
}
######
MyGo2 <- function(selected.ens,
                  universe.ens,
                  organism = "mmu",
                  #keytype = 'ENSEMBL',
                  pvalue=0.1,
                  qvalue=0.5,
                  ont="BP",
                  readable=TRUE,
                  plot=F,
                  font=0.5,
                  n=100,
                  #slim = FALSE,
                  ...){
  library(clusterProfiler)
  if(organism == "mmu"){
    library(org.Mm.eg.db)
    OrgDb <- org.Mm.eg.db 
  }
  if(organism == "hsa"){
    library(org.Hs.eg.db)
    OrgDb <- org.Hs.eg.db
  } 
  ens2ent.select=AnnotationDbi::select(OrgDb,keys = selected.ens,
                                       columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
  ens2ent.universe=AnnotationDbi::select(OrgDb,keys = universe.ens,
                                         columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
  ens2ent.select=unique(ens2ent.select$ENTREZID)
  ens2ent.universe=unique(ens2ent.universe$ENTREZID)
  
  goresult <- enrichGO(
    gene = ens2ent.select,
    universe=ens2ent.universe,
    OrgDb = OrgDb,
    #keytype =keytype,
    ont = ont,
    #pAdjustMethod = "BH",
    pvalueCutoff  = pvalue,
    qvalueCutoff  = qvalue, 
    readable=readable)
  if(plot=="TRUE"){
    emapplot(goresult
            #vertex.label.font = font,
            #n=n
            )
  }

  return(goresult)
}  
MyATACstate <- function(bed,
                        atac,
                        state,
                        atac.cutoff){
  mcols(bed)[,state] <- NA
  mcols(bed)[,state][mcols(bed)[,atac] < atac.cutoff] <- "close"
  mcols(bed)[,state][mcols(bed)[,atac] >= atac.cutoff] <- "open"
  
  MyText(paste(state,"\n","close:",sum(mcols(bed)[,state] == "close"),"\n",
               "open:",sum(mcols(bed)[,state] == "open"),"\n"))
  return(bed)
  
}
#############
MyGOwritetable <- function(genelist,
                           go.out,
                           universe.gene=genes.inf.input$EnsemblGeneID,
                           pvalue = 0.1,
                           qvalue = 1,
                           organism = "mmu",
                           ...){
  go.result <- MyGo2(genelist,
                     universe.gene,
                     pvalue = pvalue,
                     qvalue = qvalue,
                     organism = organism,
                     ...)
  MyWriteTable(go.result,
               go.out)
  return(go.result)
}


Mykeggwritetable <- 
  function(genelist,
           kegg.out,
           universe.gene=genes.inf.input$EnsemblGeneID,
           pval.adj=0.1,
           ...){
    kegg.result <- MyKegg(genelist,
                          universe.gene,
                          pval.adj = pval.adj,
                          ...
    )
    MyNCBI2Symbol <- function(mykeggid,genelist){
      keggid2symbol <- 1:length(mykeggid)
      for(i in 1:length(mykeggid)){
        keggid <- strsplit(mykeggid[i],"/")
        keggid <- keggid[[1]]
        mykeggid2symbol.tab <- genes.inf.input[genes.inf.input$NCBIID %in% keggid,]
        mykeggid2symbol.tab <- mykeggid2symbol.tab[mykeggid2symbol.tab$EnsemblGeneID %in% genelist,]
        rownames(mykeggid2symbol.tab) <- mykeggid2symbol.tab$NCBIID
        mykeggid2symbol <- mykeggid2symbol.tab[keggid,'Symbol']
        mykeggid2symbol <- paste(mykeggid2symbol,collapse = "/")
        keggid2symbol[i] <- mykeggid2symbol
        
      }
      return(keggid2symbol)
    }
    kegg.result$Symbol <- NA
    kegg.result$Symbol <- MyNCBI2Symbol(kegg.result$geneID,genelist)
    MyWriteTable(kegg.result,
                 kegg.out)
    return(kegg.result)
  }



MyReactomewritetable <- function(genelist,
                             kegg.out,
                             universe.gene=genes.inf.input$EnsemblGeneID,
                             pvalue=0.1,
                             ...){
  kegg.result <- MyReactome(genelist,
                        universe.gene,
                        pvalue  = pvalue,
                        ...
  )
  kegg.result <- as.data.frame(kegg.result)
  MyNCBI2Symbol <- function(mykeggid){
    keggid2symbol <- 1:length(mykeggid)
    for(i in 1:length(mykeggid)){
      keggid <- strsplit(mykeggid[i],"/")
      keggid <- keggid[[1]]
      mykeggid2symbol <- genes.inf.input[genes.inf.input$NCBIID %in% keggid,"Symbol"]
      mykeggid2symbol <- paste(mykeggid2symbol,collapse = "/")
      keggid2symbol[i] <- mykeggid2symbol
    }
    return(keggid2symbol)
  }
  kegg.result$Symbol <- NA
  kegg.result$Symbol <- MyNCBI2Symbol(kegg.result$geneID)
  MyWriteTable(kegg.result,
               kegg.out)
  return(kegg.result)
}


MyKegg <- 
function(selected.ens,
         universe.ens,
         organism = "mmu",
         pval.adj=0.1){
  print("Default species is mm")
  library(clusterProfiler)
  if(organism == "mmu"){  
    library(org.Mm.eg.db)
    ens2ent.select=AnnotationDbi::select(org.Mm.eg.db,keys = selected.ens,
                                         columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    ens2ent.universe=AnnotationDbi::select(org.Mm.eg.db,keys = universe.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    detach("package:org.Mm.eg.db", unload=TRUE)
  }
  if(organism == "hsa"){  
    library(org.Hs.eg.db)
    ens2ent.select=AnnotationDbi::select(org.Hs.eg.db,keys = selected.ens,
                                         columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    ens2ent.universe=AnnotationDbi::select(org.Hs.eg.db,keys = universe.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    detach("package:org.Hs.eg.db", unload=TRUE)
  }
  ens2ent.select=unique(ens2ent.select[c("ENSEMBL","ENTREZID")])
  ens2ent.universe=unique(ens2ent.universe[c("ENSEMBL","ENTREZID")])
  kegg=enrichKEGG(gene = ens2ent.select$ENTREZID,organism = organism,
                  keyType = "kegg",universe = ens2ent.universe$ENTREZID,
                  pvalueCutoff =pval.adj,qvalueCutoff = 4*pval.adj)
  
  detach("package:clusterProfiler", unload=TRUE)
  return(as.data.frame(kegg))
}


MyReactome <- 
  function(selected.ens,
           universe.ens,
           organism = "mouse",
           pvalue=0.1){
    print("Default species is mm")
    library(ReactomePA)
    if(organism == "mouse"){  
      library(org.Mm.eg.db)
      ens2ent.select=AnnotationDbi::select(org.Mm.eg.db,keys = selected.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
      ens2ent.universe=AnnotationDbi::select(org.Mm.eg.db,keys = universe.ens,
                                             columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
      detach("package:org.Mm.eg.db", unload=TRUE)
    }
    if(organism == "human"){  
      library(org.Hs.eg.db)
      ens2ent.select=AnnotationDbi::select(org.Hs.eg.db,keys = selected.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
      ens2ent.universe=AnnotationDbi::select(org.Hs.eg.db,keys = universe.ens,
                                             columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
      detach("package:org.Hs.eg.db", unload=TRUE)
    }
    ens2ent.select=unique(ens2ent.select[c("ENSEMBL","ENTREZID")])
    ens2ent.universe=unique(ens2ent.universe[c("ENSEMBL","ENTREZID")])

    kegg=ReactomePA::enrichPathway(gene = na.omit(ens2ent.select$ENTREZID),
                                   organism = organism,
                   # keyType = "kegg",
                    universe = na.omit(ens2ent.universe$ENTREZID),
                    readable = T,
                    pvalueCutoff = pvalue#,
                    #qvalueCutoff = 4*pval.adj
                   )
    detach("package:ReactomePA", unload=TRUE)
    return(kegg)
  }

Mykeggpathview <- function(kegg.result,
                           keggidlist
){
  library(clusterProfiler)
  for(i in 1:length(keggidlist)){
    browseKEGG(kegg.result,
               pathID = keggidlist[i]
    )
  }
}
Mykeggpathview2 <- function(gene.data,#format1:only NCBIIDgenelist;format2:names(NCBIID genelist) with foldchange
                            pathway.idList,##mmu05142
                            species = "mouse",
                            out.suffix,
                            kegg.native = T#TRUE:png;FALSE :output pdf
){
  library(pathview)
  for(i in 1:length(pathway.idList)){
    pathview(gene.data = gene.data,
             pathway.id = pathway.idList[i],
             species = species,
             out.suffix = out.suffix, 
             kegg.native = kegg.native)
  }
}
############
MyKegg <- function(selected.ens,
                universe.ens,
                organism = "mmu",
                pval.adj=0.1){
  print("Default species is mm")
  library(clusterProfiler)
  if(organism == "mmu"){  
    library(org.Mm.eg.db)
    ens2ent.select=AnnotationDbi::select(org.Mm.eg.db,keys = selected.ens,
                                         columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    ens2ent.universe=AnnotationDbi::select(org.Mm.eg.db,keys = universe.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    detach("package:org.Mm.eg.db", unload=TRUE)
  }
  if(organism == "hsa"){  
    library(org.Hs.eg.db)
    ens2ent.select=AnnotationDbi::select(org.Hs.eg.db,keys = selected.ens,
                                         columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    ens2ent.universe=AnnotationDbi::select(org.Hs.eg.db,keys = universe.ens,
                                           columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
    detach("package:org.Hs.eg.db", unload=TRUE)
  }
  ens2ent.select=unique(ens2ent.select[c("ENSEMBL","ENTREZID")])
  ens2ent.universe=unique(ens2ent.universe[c("ENSEMBL","ENTREZID")])
  kegg=enrichKEGG(gene = ens2ent.select$ENTREZID,organism = organism,
                  keyType = "kegg",universe = ens2ent.universe$ENTREZID,
                  pvalueCutoff =pval.adj,qvalueCutoff = 4*pval.adj)
  detach("package:clusterProfiler", unload=TRUE)
  return(as.data.frame(kegg))
}

# read table, stringsAsFactors = F
MyReadDelim <- function(input.file,
                        ...){
  data <- read.delim(input.file,
                     stringsAsFactors = FALSE,
                     ...)
  return(data)
}

Myent2ens <- function(entlist,
                      org = "mmu",
                      type = "ens",
                      gene.inf = genes.inf.input,
                      symbol = TRUE){
  
  if(org=="mmu"){
    library(org.Mm.eg.db)
    ent2ens.list <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = entlist,
                                          columns = c("ENTREZID","ENSEMBL"),keytype = "ENTREZID")
  }
  
  if(org == "hsa"){
    library(org.Hs.eg.db)
    ent2ens.list <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = entlist,
                                          columns = c("ENTREZID","ENSEMBL"),keytype = "ENTREZID")
  }
  if(type=="symbol"){
    ens.list <- gene.inf[ent2ens.list[,"ENSEMBL"],"Symbol"]
  }
  else{ens.list <- ent2ens.list[,"ENSEMBL"]}
  if(symbol == "TRUE"){
    ens.list2 <- ens.list
    ens.list <- genes.inf.input[ens.list2,"Symbol"]
  }
  return(ens.list)
  
  
}

Myens2ent <- function(enslist,
                      org = "mmu"
                      ){
  
  if(org=="mmu"){
    library(org.Mm.eg.db)
    ens2ent.list <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = enslist,
                                          columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
  }
  
  if(org == "hsa"){
    library(org.Hs.eg.db)
    ens2ent.list <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = enslist,
                                          columns = c("ENSEMBL","ENTREZID"),keytype = "ENSEMBL")
  }
  return(ens2ent.list[,2])
  
}
# show text
MyText <- function(output.text,
                   x.pos = 0.5,
                   y.pos = 0.5,
                   text.cex = 3,
                   ...){
  plot(0:1, 
       0:1, 
       col = "white", 
       xlab = "", 
       ylab = "")
  text(x.pos, 
       y.pos, 
       output.text, 
       cex = text.cex,
       ...)
}

# write .tab
MyWriteTable <- function(data,
                         output.file,
                         row.names = FALSE,
                         quote = FALSE,
                         sep = "\t",
                         col.names=T,
                         ...){
  write.table(data,
              output.file,
              row.names = row.names,
              col.names =  col.names,
              quote = quote,
              sep = sep,
              ...)
}

# check replication, calculate correlation
MyPairs <- function(data,
                    method = 'pearson',
                    ...) {
  # check replication
  method <- method
  MyPoints <- function(x, 
                       y){
    points(x, 
           y, 
           pch = 19, 
           cex = .2)
  }
  
  # calculate correlation
  ViewCor <- function(x, 
                      y,
                      ...){  
    cor.res <- cor(x, y, method = 'pearson',
                   ...)
    cor.res <- round(cor.res, 
                     3)
    
    text.x <- max(x) ^ 0.5
    text.y <- max(y) ^ 0.5
    
    output.text <- cor.res
    
    text(text.x, 
         text.y, 
         output.text, 
         cex = 3)
  }
  
  pairs(data, 
        log = "xy",
        lower.panel = MyPoints, 
        upper.panel = ViewCor,
        ...)
}

# Calculate RPKM
MyCalRpkm <- function(read.count.data, 
                      gene.length) {
  total.read.count <- colSums(read.count.data)
  rpkm.data <- t(t(read.count.data) / total.read.count) / gene.length * 1e9
  return(rpkm.data)
}

# Calculate TPM using read count
MyCalTpm <- function(read.count.data,
                     gene.length,
                     gene.exclude = NA){
  trans.data <- read.count.data / gene.length * 1000
  trans.data <- trans.data[!rownames(trans.data) %in% gene.exclude,]
  trans.sum <- colSums(trans.data)
  tpm.data <- t(t(trans.data) / trans.sum * 1000000)
  return(tpm.data)
}

# Calculate TPM using UMI
MyCalUmiTpm <- function(umi.data){
  trans.sum <- colSums(umi.data[!grepl("ERCC",rownames(umi.data)),])
  tpm.data <- t(t(umi.data) / trans.sum * 100000)
  tpm.data <- round(tpm.data,2)
  return(tpm.data)
}

# Normalize RPKM
MyNorm <- function(raw.data,
                   norm.method = "SizeFactor",
                   ercc.row = 40825:40916,
                   graph = TRUE,
                   min.cutoff = 1,
                   ...){
  # plot log2 boxplot, remove 
  LogBoxPlot <- function(data,
                         min.cutoff = min.cutoff,
                         ...){
    data.na <- data
    data.na[data.na < min.cutoff] <- NA
    boxplot(log2(data.na),
            outline = F,
            ...)
  }
  
  # normalize by size factor
  SizeFactorNorm <- function(raw.data) {
    library("DESeq2")
    size.factor <- estimateSizeFactorsForMatrix(raw.data)
    norm.data <-  t(t(raw.data) / size.factor)
    return(norm.data)
  }
  
  # normalize by ERCC size factor
  ERCCSizeFactorNorm <- function(raw.data,
                                 ercc.row) {
    library("DESeq2")
    size.factor <- estimateSizeFactorsForMatrix(raw.data[ercc.row,])
    norm.data <-  t(t(raw.data) / size.factor)
    return(norm.data)
  }
  
  if(graph){
    LogBoxPlot(raw.data,
               min.cutoff = min.cutoff)
  }
  
  if (norm.method == "SizeFactor") {
    norm.data <- SizeFactorNorm(raw.data)
  }else if (norm.method == "ERCCSizeFactor") {
    norm.data <- ERCCSizeFactorNorm(raw.data,
                                    ercc.row)
  }else if (norm.method == "Raw") {
    norm.data <- raw.data
  }
  
  if(graph){
    LogBoxPlot(norm.data,
               min.cutoff = min.cutoff)
  }
  return(norm.data)
}

# filter genes by coefficient of variation
# return filter data
MyGeneCv <- function(data,
                     cv.min,
                     cv.max=Inf
                     ){
  filter.gene <- which(apply(data,
                             1,
                             function(data){sd(data)/mean(data)}) > cv.min & 
                         apply(data,
                               1,
                               function(data){sd(data)/mean(data)}) < cv.max
                         )
  filter.data <- data[filter.gene,]
  
  return(filter.data)
}

# filter genes by gene expression
# return filter data
MyGeneExp <- function(data,
                      gene.cutoff = 1,
                      sample.cutoff = 2,
                      exp = "yes"){
  if(exp == "yes"){
    filter.gene <- which(rowSums(data > gene.cutoff) >= sample.cutoff)
  }
  if(exp == "no"){
    filter.gene <- which(rowSums(data > gene.cutoff) <= (ncol(data) - sample.cutoff))
  }
  filter.data <- data[filter.gene,]
  return(filter.data)
}

# filter genes by ERCC
# input is read count
# last 92 rows are ERCC read count
# return filter gene
MyGeneErcc <- function(data,
                       graph = TRUE,
                       padjA.cutoff = 0.1){
  library(genefilter)
  
  library(statmod)
  library(DESeq2)
  
  countsMmus <- data[1:(nrow(data)-92),]#hs:first 91 
  countsERCC <- data[(nrow(data)-91):nrow(data),]#hs:92 to the end 
  
  sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
  sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
  rbind( sfMmus, sfERCC )
  
  nCountsERCC <- t( t(countsERCC) / sfERCC )
  nCountsMmus <- t( t(countsMmus) / sfMmus )
  
  meansERCC <- rowMeans( nCountsERCC )
  varsERCC <- rowVars( nCountsERCC )
  cv2ERCC <- varsERCC / meansERCC^2
  meansMmus <- rowMeans( nCountsMmus )
  varsMmus <- rowVars( nCountsMmus )
  cv2Mmus <- varsMmus / meansMmus^2
  
  minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
  useForFitA <- meansERCC >= minMeanForFitA
  minMeanForFitA
  table( useForFitA )
  
  fitA <- glmgam.fit(cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ), 
                     cv2ERCC[useForFitA])
  
  minBiolDisp <- .5^2
  xi <- mean( 1 / sfERCC )
  m <- ncol(countsMmus)
  psia1thetaA <- mean( 1 / sfERCC ) + 
    ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
  cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
  testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
  pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
  padjA <- p.adjust( pA, "BH" )
  table( padjA < padjA.cutoff )
  
  predictedNoiseCV2 <- psia1thetaA / meansMmus + coefficients(fitA)["a0"]
  useInPlot <- meansMmus>minMeanForFitA
  selectGenes <- rownames(countsMmus)[padjA < padjA.cutoff & !is.na(padjA)]
  
  if(graph){
    plot( x = meansMmus[meansMmus >0 & cv2Mmus > 0], y =  cv2Mmus[meansMmus >0 & cv2Mmus > 0], 
          col = ifelse( padjA[meansMmus >0 & cv2Mmus > 0] < padjA.cutoff, "#C0007090", "#70500040" ),pch=20, cex=.2,
          log="xy", ylim=c(.01,100),
          xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)",
          main = paste(ncol(data),"Samples,",length(selectGenes),"Genes"))
    
    
    # Plot the plant genes, use a different color if they are highly variable
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, coefficients(fitA)["a1tilde"] / xg + coefficients(fitA)["a0"], col="#FF000080", lwd=3 )
    
    points( meansERCC, cv2ERCC, pch=20, cex=1,  col="#0060B8A0"  )
  }
  return(selectGenes)
}

# # filter genes by all genes
# # input is read count
# # return filter gene
# MyGeneAll <- function(){
#   
# }

# filter genes by PCA
# return filter gene
MyGenePca <- function(pca.res,
                      dim,
                      proba = 0.05,
                      direction = "both"){
  dim.res <- dimdesc(pca.res, 
                     axes = dim,
                     proba = proba)
  
  
  if(length(dim) == 1){
    dim.res <- as.data.frame(dim.res)
    if(direction == "both"){
      filter.gene <- row.names(dim.res)
    }else if(direction == "pos"){
      filter.gene <- row.names(dim.res[dim.res[,1] > 0, ])
    }else if(direction == "neg"){
      filter.gene <- row.names(dim.res[dim.res[,1] < 0, ])
    }else{
      stop("direction = both/pos/neg")
    }
  }
  if(length(dim) > 1){
    dim.res.merge <- do.call("rbind", unlist(dim.res, recursive = FALSE))
    if(direction == "both"){
      filter.gene <- row.names(dim.res.merge)
    }else if(direction == "pos"){
      filter.gene <- row.names(dim.res.merge[dim.res.merge[,1] > 0, ])
    }else if(direction == "neg"){
      filter.gene <- row.names(dim.res.merge[dim.res.merge[,1] < 0, ])
    }else{
      stop("direction = both/pos/neg")
    }
  }
  filter.gene <- unique(filter.gene)
  return(filter.gene)
}

# use log data to do pca
MyLogPca <- function(data,
                     log.add = 0.1,
                     graph = FALSE,
                     ...){
  library(FactoMineR)
  pca.data <- t(log2(data + log.add))
  pca.res <- FactoMineR::PCA(pca.data, 
                             graph = graph,
                             ...)
  return(pca.res)
}

# PCA result to ggplot2
MyPca2Gg <- function(pca.res,
                     sample.inf,
                     x.dim = 1,
                     y.dim = 2,
                     color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  pca <- as.data.frame(pca.res$eig)
  pc.percent <- pca$`percentage of variance`
  pca.coord <- pca.res$ind$coord
  pca.coord <- as.data.frame(pca.coord)
  pca.coord.ind.sup <- pca.res$ind.sup$coord
  pca.coord.ind.sup <- as.data.frame(pca.coord.ind.sup)
  pca.coord.merge <- rbind(pca.coord,
                           pca.coord.ind.sup)
  row.names(sample.inf) <- sample.inf$SampleName
  pca.df <- cbind(x.pos = pca.coord.merge[,x.dim],
                  y.pos = pca.coord.merge[,y.dim],
                  sample.inf[row.names(pca.coord.merge),])
  
  p <- ggplot(data = pca.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("PC",x.dim,"(",round(pc.percent[x.dim],1),"%)",sep = ""))
  p <- p + ylab(paste("PC",y.dim,"(",round(pc.percent[y.dim],1),"%)",sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    
    #scale_y_reverse() +
    scale_color_manual(values = c(time.colors,"gray30")) +
    #scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1,
                                 ncol = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2,
                                ncol = 1))+
    scale_color_manual(values = c(time.colors,"gray30")) +
    #scale_shape_manual(values = shape.group) +
    guides(colour = guide_legend(title = "",
                                 order = 1)) +
    theme(axis.text = element_text(size = 20, colour = "black")) +
    theme(axis.title.x = element_text(size = 40, colour = "black")) +
    theme(axis.title.y = element_text(size = 40, colour = "black")) +
    theme(legend.text = element_text(size = 45, colour = "black")) +
    theme(legend.title = element_text(size = 45, colour = "black")) +
    theme(panel.border = element_rect(size = 4,
                                      colour = "black")) +
    theme(legend.key = element_blank()) +
    guides(colour = guide_legend(title = "",
                                 keywidth = 3,
                                 keyheight = 3,
                                 override.aes = list(size = 12),
                                 order = 1,
                                 ncol = 1)) +
    guides(shape = guide_legend(title = "",
                                keywidth = 3,
                                keyheight = 3,
                                override.aes = list(size = 12),
                                order = 2,
                                ncol = 1))
  return(p)
}

# PCA PC1 PC2 PC3 to plotly 3d scatter
# MyPca2Plotly3d <- function(pca.res,
#                            text,
#                            color = color,
#                            x.dim = 1,
#                            y.dim = 2,
#                            z.dim = 3,
#                            color_manual = brewer.pal(9,"Set1"),
#                            ...){
#   library(plotly)
#   color_manual = color_manual[order(unique(color))]
#   pca.coord <- pca.res$ind$coord
#   pca.coord <- as.data.frame(pca.coord)
#   pca.df <- data.frame(text = text,
#                        x = pca.coord[,x.dim],
#                        y = pca.coord[,y.dim],
#                        z = pca.coord[,z.dim])
# 
#   p <- plot_ly(pca.df,
#                x = x,
#                y = y,
#                z = z,
#                text = text,
#                type = "scatter3d",
#                mode = "markers",
#                colors = color_manual,
#                color = color,
#                ...)
#   return(p)
# }
MyPca2Plotly3d <- function(pca.res,
                           x.dim = 1,
                           y.dim = 2,
                           z.dim = 3,
                           ...){
  library(plotly)
  p <- plot_ly(x = pca.res$ind$coord[,x.dim],
               y = pca.res$ind$coord[,y.dim],
               z = pca.res$ind$coord[,z.dim],
               type = "scatter3d",
               mode = "markers",
               ...)
  return(p)
  detach("package:plotly",
         unload=TRUE)
}

# post plotly to web
MyPlotlyPost <- function(p,
                         filename,
                         sharing = "public",
                         fileopt = "overwrite",
                         username = "qwl",
                         key = "rk91eovugx"){
  library(plotly)
  
  Sys.setenv("plotly_username" = username)
  Sys.setenv("plotly_api_key" = key)
  plotly_POST(p, 
              filename = filename, 
              sharing = sharing,
              fileopt = fileopt)
}

# heatmap
MyHeatmap <- function(data,
                      graph = TRUE,
                      type = c("log.row.zscore",
                               "log.row.relat",
                               "row.relat",
                               "raw",
                               "log.raw",
                               "row.zscore"),
                      return.tree = "none",
                      color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), 
                                                       space="Lab"),
                      
                      hc.c.data.type = c("raw",
                                         "log.row.zscore",
                                         "log.row.relat",
                                         "row.relat",
                                         "row.zscore",
                                         "log.raw"),
                      hc.r.data.type = c("raw",
                                         "log.raw",
                                         "log.row.zscore",
                                         "row.zscore",
                                         "log.row.relat",
                                         "row.relat"),
                      c.cov.method = "spearman",
                      r.cov.method = "spearman",
                      c.hc.method = "ward.D",
                      r.hc.method = "ward.D",
                      Colv = "do",
                      Rowv = "do",
                      reoder.row = F,
                      density.info = "density",
                      labRow = NA,
                      margins = c(8,2),
                      trace = "none",
                      key.title = "",
                      show.count = T,
                      dist = "none",
                      ...){
  heatmap.3 <- function(x,
                        Rowv = TRUE, 
                        Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",
                        key.title = "",
                        ...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        # title("Color Key\nand Density Plot")
        title(key.title)
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        # title("Color Key\nand Histogram")
        title(key.title)
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  Myzscore <- function(input.data){
    data.zscore <- (input.data-mean(input.data))/sd(input.data)
    return(data.zscore)
  }
  hc.c.data.type <- match.arg(hc.c.data.type)
  if(hc.c.data.type == "raw"){
    hc.c.data <- data
  }else if(hc.c.data.type == "log.raw"){
    hc.c.data <- as.matrix(log2(data + 1))
  }else if(hc.c.data.type == "log.row.relat"){
    hc.c.data <- as.matrix(log2(data + 1) / rowMax(log2(data + 1)))
  }else if(hc.c.data.type == "row.relat"){
    hc.c.data <- as.matrix(data / rowMax(data))
  }
  else if(hc.c.data.type == "row.zscore"){
    hc.c.data <- as.matrix(apply(data,1,Myzscore))}
  else if(hc.c.data.type == "log.row.zscore"){
    hc.c.data <- as.matrix(apply(log2(data+1),1,Myzscore))}
  
  if(Colv == "do"){
    if(dist == "euclidean"){
      c.hc <- hclust(dist(t(hc.c.data)), 
                     method = c.hc.method)
    }
    else{
      c.cov <- cor(hc.c.data,
                   method = c.cov.method)  
      c.hc <- hclust(as.dist(abs(1 - c.cov)), 
                     method = c.hc.method)
    }
    
    Colv <- as.dendrogram(c.hc)
    Colv <- reorder(Colv, 
                    1:100000, 
                    mean)
  }
  
  hc.r.data.type <- match.arg(hc.r.data.type)
  if(hc.r.data.type == "raw"){
    hc.r.data <- data
  }else if(hc.r.data.type == "log.raw"){
    hc.r.data <- as.matrix(log2(data + 1))
  }else if(hc.r.data.type == "log.row.relat"){
    hc.r.data <- as.matrix(log2(data + 1) / rowMax(log2(data + 1)))
  }else if(hc.r.data.type == "row.relat"){
    hc.r.data <- as.matrix(data / rowMax(data))
  }else if(hc.r.data.type == "row.zscore"){
    hc.r.data <- as.matrix(apply(data,1,Myzscore))}
  else if(hc.r.data.type == "log.row.zscore"){
    hc.r.data <- as.matrix(apply(log2(data+1),1,Myzscore))}
  
  if(Rowv == "do"){
    r.cov <- cor(t(hc.r.data),
                 method = r.cov.method)
    r.hc <- hclust(as.dist(abs(1 - r.cov)), 
                   method = r.hc.method)  
    
    Rowv <- as.dendrogram(r.hc)
    if(reoder.row) {
      Rowv <- reorder(Rowv,
                      1:10000,
                      mean)
      Rowv <- rev(Rowv)
    }
  }
  
  if(graph){
    if(show.count){
      MyText(paste(dim(data)[2],"samples,",dim(data)[1],"genes"))
    }
    type <- match.arg(type)
    if(type == "log.row.zscore"){
      heatmap.data <- as.matrix(log2(data + 1))
      heatmap.3(heatmap.data,
                scale = "row",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }
    else if(type == "row.zscore"){
      heatmap.data <- as.matrix(data)
      heatmap.3(heatmap.data,
                scale = "row",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }else if(type == "log.row.relat"){
      heatmap.data <- as.matrix(log2(data + 1) / rowMax(log2(data + 1)))
      heatmap.3(heatmap.data,
                KeyValueName = "Row Relative Value",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }else if(type == "row.relat"){
      heatmap.data <- as.matrix(data / rowMax(data))
      heatmap.3(heatmap.data,
                KeyValueName = "Row Relative Value",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }else if(type == "raw"){
      heatmap.data <- as.matrix(data)
      heatmap.3(heatmap.data,
                # KeyValueName = "Row Relative Value",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }
    else if(type == "log.raw"){
      heatmap.data <- as.matrix(log2(data + 1))
      heatmap.3(heatmap.data,
                KeyValueName = "log2(Raw)",
                Colv = Colv,
                Rowv = Rowv,
                col = color.palette,
                density.info = density.info,
                trace = trace,
                labRow = labRow,
                margins = margins,
                key.title = key.title,
                ...)
    }
  }
  if (return.tree == "col"){
    return(c.hc)
  }else if (return.tree == "row"){
    return(r.hc)
  }
}

# cut tree and assign color
Tree2Col <- function(tree,
                     count,
                     palette,
                     graph = TRUE){
  if(graph){
    plot(1:count,
         1:count,
         col = palette[1:count],
         pch = 20,
         cex = 5,
         xlab = "",
         ylab = "")
  }
  
  color.res <- palette[cutree(tree,
                              count)]
  return(color.res)
}

# cut tree and get a branch of tree
MyTree2Name <- function(tree,
                        count,
                        branch){
  cutree.res <- cutree(tree,
                       count)
  breach.name <- cutree.res[cutree.res %in% branch]
  breach.name <- names(breach.name)
  return(breach.name)
}

# convert name to color
MyName2Col <- function(name,
                       palette,
                       is.row = FALSE){
  name.factor <- as.factor(name)
  name.color <- palette[name.factor]
  name.color <- as.matrix(name.color)
  if(is.row){
    name.color <- t(name.color)
  }
  return(name.color)
}
# Deseq2 pipeline
MyRunDeseq2 <- function(rc.data, 
                        condition1.samples, 
                        condition2.samples,
                        graph = TRUE,
                        alpha = 0.05,
                        ...){
  library(DESeq2)
  
  count.data <- as.matrix(rc.data[ ,c(condition1.samples, condition2.samples)])
  samples.name <- colnames(count.data)
  condition <- rep(NA, ncol(count.data))
  condition[samples.name %in% condition1.samples] <- "condition1"
  condition[samples.name %in% condition2.samples] <- "condition2"
  col.data <- as.data.frame(condition)
  dds <- DESeqDataSetFromMatrix(countData = count.data, colData = col.data, design = ~ condition)
  dds$condition <- factor(dds$condition, levels = c("condition1","condition2"))
  # dds <- DESeq(dds, minReplicatesForReplace = 3)
  dds <- DESeq(dds, 
               ...)
  res <- results(dds)
  if(graph){
    plotMA(res, 
           alpha = alpha)
  }
  return(res)
}

# get DE genes from Deseq2 resulut
MyFilterDeseq2 <- function(de.res,
                           fc,
                           padj.cutoff = 0.05){
  de.res <- as.data.frame(de.res)
  
  de.filter.res <- de.res[!is.na(de.res$padj),]
  de.filter.res <- de.filter.res[de.filter.res$padj < padj.cutoff,]
  if(fc == "pos"){
    de.filter.res <- de.filter.res[de.filter.res$log2FoldChange > 0,]
  }else if(fc == "neg"){
    de.filter.res <- de.filter.res[de.filter.res$log2FoldChange < 0,]
  }else{
    stop("fc = pos/neg")
  }
  return(de.filter.res)
}

# get DE genes from Deseq2 resulut (foldchange + padj)
MyFilterDeseq2Fc <- function(de.res,
                             log2fc = 0,
                             padj.cutoff = 0.05){
  de.res <- as.data.frame(de.res)
  
  de.filter.res <- de.res[!is.na(de.res$padj),]
  de.filter.res <- de.filter.res[de.filter.res$padj < padj.cutoff,]
  if(log2fc > 0){
    de.filter.res <- de.filter.res[de.filter.res$log2FoldChange > log2fc,]
  }else if(log2fc < 0){
    de.filter.res <- de.filter.res[de.filter.res$log2FoldChange < log2fc,]
  }
  return(de.filter.res)
}

# judge DE genes by fold change
MyFilterFc <- function(gene.name,
                       up.data,
                       down.data,
                       fc = 2,
                       up.min = 1){
  de.gene <- gene.name[up.data > up.min &
                         (down.data / up.data) < (1 / fc) ]
  return(de.gene)
}

# barplot of annotation
MyAnnBarPlot <- function(term,
                         value,
                         color = "gray",
                         xlab = "-log10(p-value)",
                         term.cex = 1,
                         ...){
  data <- data.frame(term = as.character(rev(term)),
                     value = as.numeric(rev(value)),
                     color = as.character(rev(color)),
                     stringsAsFactors = F)
  
  x.axis.length <- ceiling(max(data$value))
  
  op <- par(mar = c(4, 1, 1, 1))
  barplot.out <- barplot(data$value, 
                         horiz = T, 
                         axes = F,
                         col = data$color,
                         xlim = c(-1.5 * x.axis.length,
                                  x.axis.length))
  text(x = -0.2,
       y = barplot.out,
       label = data$term,
       adj = 1,
       cex = term.cex)
  
  axis(side = 1,
       at = seq(0,x.axis.length,2),
       labels = seq(0,x.axis.length,2),
       tick = T)
  
  mtext(text = xlab, 
        side = 1,
        line = 2,
        at = par('usr')[2] / 2,
        cex = 1)
  par(op)
}

# annotate gene to KEGG database using GOstats 
MyGostatsKegg <- function(selected.ens,
                          universe.ens){
  library(GOstats)
  library(org.Mm.eg.db)
  selectedIDs <- mget(selected.ens, 
                      revmap(org.Mm.egENSEMBL),
                      ifnotfound = NA)
  selectedIDs <- unlist(selectedIDs)
  selectedIDs <- selectedIDs[!is.na(selectedIDs)]
  selectedIDs <- unique(selectedIDs)
  
  universeIDs <- mget(universe.ens, 
                      revmap(org.Mm.egENSEMBL),
                      ifnotfound = NA)
  universeIDs <- unlist(universeIDs)
  universeIDs <- universeIDs[!is.na(universeIDs)]
  universeIDs <- unique(universeIDs)
  KeggParams <- new("KEGGHyperGParams", 
                    geneIds = selectedIDs, 
                    universeGeneIds = universeIDs, 
                    annotation ="org.Mm.eg")
  KeggResults <- hyperGTest(KeggParams)
  out.res <- summary(KeggResults)
  out.res <- out.res[,c("KEGGID",
                        "Term",
                        "Count",
                        "Size",
                        "Pvalue",
                        "OddsRatio",
                        "ExpCount")]
  return(out.res)
}

# annotate gene to GO database using GOstats 
MyGo <- function(selected.ens,
                 universe.ens,
                 ontology.type,
                 gene.inf){
  library(GOstats)
  library(org.Mm.eg.db)
  selectedIDs <- mget(selected.ens, 
                      revmap(org.Mm.egENSEMBL),
                      ifnotfound = NA)
  selectedIDs <- unlist(selectedIDs)
  selectedIDs <- selectedIDs[!is.na(selectedIDs)]
  selectedIDs <- unique(selectedIDs)
  
  universeIDs <- mget(universe.ens, 
                      revmap(org.Mm.egENSEMBL),
                      ifnotfound = NA)
  universeIDs <- unlist(universeIDs)
  universeIDs <- universeIDs[!is.na(universeIDs)]
  universeIDs <- unique(universeIDs)
  if(ontology.type == "KEGG"){
    KeggParams <- new("KEGGHyperGParams", 
                      geneIds = selectedIDs, 
                      universeGeneIds = universeIDs, 
                      annotation ="org.Mm.eg")
    KeggResults <- hyperGTest(KeggParams)
    out.res <- summary(KeggResults)
    if(nrow(out.res) != 0){
      out.res <- out.res[,c("KEGGID",
                            "Term",
                            "Count",
                            "Size",
                            "Pvalue",
                            "OddsRatio",
                            "ExpCount")]
      out.res$GeneID <- NA
      out.res$EnsemblID <- NA
      out.res$Symbol <- NA
      for(term.num in 1:nrow(out.res)){
        term.id <- out.res[term.num,1]
        all.gene.id <- KeggResults@catToGeneId[[term.id]]
        term.gene.id <- selectedIDs[selectedIDs %in% all.gene.id]
        term.gene.id <- sort(term.gene.id)
        term.gene.ens <- mget(term.gene.id, 
                              org.Mm.egENSEMBL,
                              ifnotfound = NA)
        term.gene.ens <- unlist(term.gene.ens)
        term.gene.ens <- term.gene.ens[!is.na(term.gene.ens)]
        term.gene.ens <- unique(term.gene.ens)
        term.gene.ens <- sort(term.gene.ens)
        term.gene.sym <- gene.inf[gene.inf$EnsemblGeneID %in% term.gene.ens,"Symbol"]
        term.gene.sym <- unique(term.gene.sym)
        term.gene.sym <- sort(term.gene.sym)
        out.res[term.num,"GeneID"] <- paste(term.gene.id, 
                                            collapse = " ")
        out.res[term.num,"EnsemblID"] <- paste(term.gene.ens, 
                                               collapse = " ")
        out.res[term.num,"Symbol"] <- paste(term.gene.sym, 
                                            collapse = " ")
      }
    }
  }
  if(ontology.type == "BP"){
    params <- new("GOHyperGParams",
                  geneIds=selectedIDs,
                  universeGeneIds=universeIDs,
                  annotation="org.Mm.eg",
                  ontology=ontology.type)
    GoResults <- hyperGTest(params)
    out.res <- summary(GoResults)
    if(nrow(out.res) != 0){
      out.res <- out.res[,c("GOBPID",
                            "Term",
                            "Count",
                            "Size",
                            "Pvalue",
                            "OddsRatio",
                            "ExpCount")]
      out.res$GeneID <- NA
      out.res$EnsemblID <- NA
      out.res$Symbol <- NA
      for(term.num in 1:nrow(out.res)){
        term.id <- out.res[term.num,1]
        all.gene.id <- GoResults@goDag@nodeData@data[[term.id]]["geneIds"]
        all.gene.id <- unlist(all.gene.id)
        term.gene.id <- selectedIDs[selectedIDs %in% all.gene.id]
        term.gene.id <- sort(term.gene.id)
        term.gene.ens <- mget(term.gene.id, 
                              org.Mm.egENSEMBL,
                              ifnotfound = NA)
        term.gene.ens <- unlist(term.gene.ens)
        term.gene.ens <- term.gene.ens[!is.na(term.gene.ens)]
        term.gene.ens <- unique(term.gene.ens)
        term.gene.ens <- sort(term.gene.ens)
        term.gene.sym <- gene.inf[gene.inf$EnsemblGeneID %in% term.gene.ens,"Symbol"]
        term.gene.sym <- unique(term.gene.sym)
        term.gene.sym <- sort(term.gene.sym)
        out.res[term.num,"GeneID"] <- paste(term.gene.id, 
                                            collapse = " ")
        out.res[term.num,"EnsemblID"] <- paste(term.gene.ens, 
                                               collapse = " ")
        out.res[term.num,"Symbol"] <- paste(term.gene.sym, 
                                            collapse = " ")
      }
    }
  }
  if(ontology.type == "CC"){
    params <- new("GOHyperGParams",
                  geneIds=selectedIDs,
                  universeGeneIds=universeIDs,
                  annotation="org.Mm.eg",
                  ontology=ontology.type)
    GoResults <- hyperGTest(params)
    out.res <- summary(GoResults)
    if(nrow(out.res) != 0){
      out.res <- out.res[,c("GOCCID",
                            "Term",
                            "Count",
                            "Size",
                            "Pvalue",
                            "OddsRatio",
                            "ExpCount")]
      out.res$GeneID <- NA
      out.res$EnsemblID <- NA
      out.res$Symbol <- NA
      for(term.num in 1:nrow(out.res)){
        term.id <- out.res[term.num,1]
        all.gene.id <- GoResults@goDag@nodeData@data[[term.id]]["geneIds"]
        all.gene.id <- unlist(all.gene.id)
        term.gene.id <- selectedIDs[selectedIDs %in% all.gene.id]
        term.gene.id <- sort(term.gene.id)
        term.gene.ens <- mget(term.gene.id, 
                              org.Mm.egENSEMBL,
                              ifnotfound = NA)
        term.gene.ens <- unlist(term.gene.ens)
        term.gene.ens <- term.gene.ens[!is.na(term.gene.ens)]
        term.gene.ens <- unique(term.gene.ens)
        term.gene.ens <- sort(term.gene.ens)
        term.gene.sym <- gene.inf[gene.inf$EnsemblGeneID %in% term.gene.ens,"Symbol"]
        term.gene.sym <- unique(term.gene.sym)
        term.gene.sym <- sort(term.gene.sym)
        out.res[term.num,"GeneID"] <- paste(term.gene.id, 
                                            collapse = " ")
        out.res[term.num,"EnsemblID"] <- paste(term.gene.ens, 
                                               collapse = " ")
        out.res[term.num,"Symbol"] <- paste(term.gene.sym, 
                                            collapse = " ")
      }
    }
  }
  if(ontology.type == "MF"){
    params <- new("GOHyperGParams",
                  geneIds=selectedIDs,
                  universeGeneIds=universeIDs,
                  annotation="org.Mm.eg",
                  ontology=ontology.type)
    GoResults <- hyperGTest(params)
    out.res <- summary(GoResults)
    if(nrow(out.res) != 0){
      out.res <- out.res[,c("GOMFID",
                            "Term",
                            "Count",
                            "Size",
                            "Pvalue",
                            "OddsRatio",
                            "ExpCount")]
      out.res$GeneID <- NA
      out.res$EnsemblID <- NA
      out.res$Symbol <- NA
      for(term.num in 1:nrow(out.res)){
        term.id <- out.res[term.num,1]
        all.gene.id <- GoResults@goDag@nodeData@data[[term.id]]["geneIds"]
        all.gene.id <- unlist(all.gene.id)
        term.gene.id <- selectedIDs[selectedIDs %in% all.gene.id]
        term.gene.id <- sort(term.gene.id)
        term.gene.ens <- mget(term.gene.id, 
                              org.Mm.egENSEMBL,
                              ifnotfound = NA)
        term.gene.ens <- unlist(term.gene.ens)
        term.gene.ens <- term.gene.ens[!is.na(term.gene.ens)]
        term.gene.ens <- unique(term.gene.ens)
        term.gene.ens <- sort(term.gene.ens)
        term.gene.sym <- gene.inf[gene.inf$EnsemblGeneID %in% term.gene.ens,"Symbol"]
        term.gene.sym <- unique(term.gene.sym)
        term.gene.sym <- sort(term.gene.sym)
        out.res[term.num,"GeneID"] <- paste(term.gene.id, 
                                            collapse = " ")
        out.res[term.num,"EnsemblID"] <- paste(term.gene.ens, 
                                               collapse = " ")
        out.res[term.num,"Symbol"] <- paste(term.gene.sym, 
                                            collapse = " ")
      }
    }
  }
  return(out.res)
}
# Venn plot
Myvenn <- 
  function(gene.venn,col.list=time.colors,lwd=5,...){
    gene.venn.them <- compute.Venn(gene.venn)
    gene.venn.them.p <- VennThemes(gene.venn.them)
    
    gene.venn.them.p$Set$Set1$col <- col.list[1]
    gene.venn.them.p$Set$Set2$col <- col.list[2]
    
    gene.venn.them.p$Set$Set1$lwd <- lwd
    gene.venn.them.p$Set$Set2$lwd <- lwd
    
    gene.venn.them.p$SetText$Set1$fontsize <- 0.01
    gene.venn.them.p$SetText$Set2$fontsize <- 0.01
    
    
    gene.venn.them.p$FaceText$`11`$fontsize <- 0
    gene.venn.them.p$FaceText$`10`$fontsize <- 0
    gene.venn.them.p$FaceText$`01`$fontsize <- 0
    
    gene.venn.them.p$Face$`11`$fill <- NULL
    gene.venn.them.p$Face$`10`$fill <- NULL
    gene.venn.them.p$Face$`01`$fill <- NULL
    return(gene.venn.them.p)
}

# plot Deseq result
# MyDePlot <- function(de.data,
#                      padj.cutoff = 0.05,
#                      up.color = "red",
#                      down.color = "blue",
#                      undiff.color = "black",
#                      xlab = "MeanExpression",
#                      ylab = "log2FoldChange",
#                      ...){
#   de.data <- as.data.frame(de.data)
#   de.color <- rep(undiff.color, nrow(de.data))
#   de.color[(!is.na(de.data$padj)) & 
#              (de.data$padj < padj.cutoff) & 
#              (de.data$log2FoldChange < 0)] <- down.color
#   de.color[(!is.na(de.data$padj)) & 
#              (de.data$padj < padj.cutoff) & 
#              (de.data$log2FoldChange > 0)] <- up.color
#   up.count <- sum((!is.na(de.data$padj)) & 
#                     (de.data$padj < padj.cutoff) & 
#                     (de.data$log2FoldChange > 0))
#   down.count <- sum((!is.na(de.data$padj)) & 
#                       (de.data$padj < padj.cutoff) & 
#                       (de.data$log2FoldChange < 0))
#   
#   fc.max <- max(abs(de.data$log2FoldChange[!is.na(de.data$log2FoldChange)]))
#   fc.max <- ceiling(fc.max)
#   
#   plot(de.data$baseMean,
#        de.data$log2FoldChange,
#        xlab = xlab,
#        ylab = ylab,
#        ylim = c(-fc.max, fc.max),
#        col = de.color,
#        log = "x",
#        pch = 20,
#        cex = 0.5
#        )
#   text(1,
#        fc.max / 2,
#        as.character(up.count),
#        col = up.color,
#        cex = 1)
#   text(1,
#        -(fc.max / 2),
#        as.character(down.count),
#        col = down.color,
#     
# cex = 1)
# }


MyDePlot <- function(de.data,
                     padj.cutoff = 0.05,
                     fc.max,
                     xlim,
                     up.color = "#F25362",
                     down.color = "#487CC2",
                     undiff.color = "black",
                     xlab = "mean of normalized counts",
                     ylab = "log2 fold change",
                     fc.cutoff=2,
                     count = T,
                     pval=F,
                     ...){
  de.data <- as.data.frame(de.data)
  de.color <- rep(undiff.color, nrow(de.data))
  
  if(pval==F){  de.color[(!is.na(de.data$padj)) & 
                (de.data$padj < padj.cutoff) & 
                (de.data$log2FoldChange < -log2(fc.cutoff))] <- down.color
    de.color[(!is.na(de.data$padj)) & 
               (de.data$padj < padj.cutoff) & 
               (de.data$log2FoldChange > log2(fc.cutoff))] <- up.color
    up.count <- sum((!is.na(de.data$padj)) & 
                      (de.data$padj < padj.cutoff) & 
                      (de.data$log2FoldChange > log2(fc.cutoff)))
    down.count <- sum((!is.na(de.data$padj)) & 
                        (de.data$padj < padj.cutoff) & 
                        (de.data$log2FoldChange < -log2(fc.cutoff)))}
   if(pval==T){
    de.color[(!is.na(de.data$pvalue)) & 
               (de.data$pvalue < padj.cutoff) & 
               (de.data$log2FoldChange < -log2(fc.cutoff))] <- down.color
    de.color[(!is.na(de.data$pvalue)) & 
               (de.data$pvalue < padj.cutoff) & 
               (de.data$log2FoldChange > log2(fc.cutoff))] <- up.color
    up.count <- sum((!is.na(de.data$pvalue)) & 
                      (de.data$pvalue < padj.cutoff) & 
                      (de.data$log2FoldChange > log2(fc.cutoff)))
    down.count <- sum((!is.na(de.data$pvalue)) & 
                        (de.data$pvalue < padj.cutoff) & 
                        (de.data$log2FoldChange < -log2(fc.cutoff)))
   }
  
  
  shape.plot <- rep(16,length(de.data$log2FoldChange))
  shape.plot[abs(de.data$log2FoldChange) > fc.max] <- 17
  de.data$log2FoldChange[de.data$log2FoldChange > fc.max] <- fc.max
  de.data$log2FoldChange[de.data$log2FoldChange < -fc.max] <- -fc.max
  plot(de.data$baseMean,
       de.data$log2FoldChange,
       xlab = xlab,
       ylab = ylab,
       ylim = c(-fc.max, fc.max),
       xlim = xlim,
       col = de.color,
       log = "x",
       pch = shape.plot,
       cex = 1)
  box(lwd=2)
  if(count){
    text(1,
         fc.max / 2,
         as.character(up.count),
         col = up.color,
         cex = 1)
    text(1,
         -(fc.max / 2),
         as.character(down.count),
         col = down.color,
         cex = 1)
  }
}

# calculate percentile rank
MyPerRank <- function(data){
  res <- rank(data) / length(data) * 100
  return(res)
}

# convert PCA result to minimum spanning tree ggplot2
MyPca2Mst <- function(
  pca.res,
  sample.inf,
  x.dim = 1,
  y.dim = 2,
  color_manual = brewer.pal(9,"Set1")){
  
  library(ggplot2)
  library(monocle)
  pc.percent <- pca.res$eig$"percentage of variance"
  pca.coord <- pca.res$ind$coord
  pca.coord <- as.data.frame(pca.coord)
  row.names(sample.inf) <- sample.inf$SampleName
  pca.df <- cbind(x.pos = pca.coord[,x.dim],
                  y.pos = pca.coord[,y.dim],
                  sample.inf[row.names(pca.coord),])
  
  
  adjusted_S <- pca.coord[,c(x.dim,y.dim)]
  dp <- as.matrix(dist(adjusted_S))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  
  edge_df <- merge(pca.df, edge_list, by.x = "SampleName", 
                   by.y = "source", all = TRUE)
  edge_df <- rename(edge_df, c(x.pos = "x.pos.source", 
                               y.pos = "y.pos.source"))
  edge_df <- merge(edge_df, pca.df[, c("SampleName", 
                                       "x.pos", "y.pos")], 
                   by.x = "target", 
                   by.y = "SampleName", 
                   all = TRUE)
  edge_df <- rename(edge_df, c(x.pos = "x.pos.target", 
                               y.pos = "y.pos.target")) 
  p <- ggplot(data = edge_df)
  p <- p + xlab(paste("PC",x.dim,"(",round(pc.percent[x.dim],1),"%)",sep = ""))
  p <- p + ylab(paste("PC",y.dim,"(",round(pc.percent[y.dim],1),"%)",sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  p.point <- p + geom_point(aes(x = x.pos.source, 
                                y = y.pos.source,
                                label = SampleName),
                            col = "white",
                            na.rm = TRUE) 
  p.segment <- p.point +  geom_segment(aes(x = x.pos.source, 
                                           y = y.pos.source, 
                                           xend = x.pos.target, 
                                           yend = y.pos.target),
                                       data = edge_df,
                                       col = "white",
                                       na.rm = TRUE)
  p.segment
}
# # convert PCA result to minimum spanning tree ggplot2(no backbone)
# MyPca2Mst <- function(pca.res,
#                       sample.inf,
#                       x.dim = 1,
#                       y.dim = 2,
#                       color_manual = brewer.pal(9,"Set1")){
#   library(ggplot2)
#   pc.percent <- pca.res$eig$"percentage of variance"
#   pca.coord <- pca.res$ind$coord
#   pca.coord <- as.data.frame(pca.coord)
#   row.names(sample.inf) <- sample.inf$SampleName
#   pca.df <- cbind(x.pos = pca.coord[,x.dim],
#                   y.pos = pca.coord[,y.dim],
#                   sample.inf[row.names(pca.coord),])
#   
#   p <- ggplot(data = pca.df, 
#               mapping = aes(x = x.pos, 
#                             y = y.pos))
#   p <- p + xlab(paste("PC",x.dim,"(",round(pc.percent[x.dim],2),"%)",sep = ""))
#   p <- p + ylab(paste("PC",y.dim,"(",round(pc.percent[y.dim],2),"%)",sep = ""))
#   p <- p + guides(colour=guide_legend(title=NULL))
#   p <- p + scale_color_manual(values = color_manual)
#   p <- p + 
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank())
#   library(ape)
#   m <- pca.df[,1:2]
#   r <- mst(dist(m))
#   n <- dim(r)[1]
#   w <- which(r!=0)
#   i <- as.vector(row(r))[w]
#   j <- as.vector(col(r))[w]
#   segment.df <- data.frame(x = m[i,1], 
#                            y = m[i,2], 
#                            xend = m[j,1], 
#                            yend =m[j,2])
#   
#   p.point <- p + geom_point() 
#   p.segment <- p.point +  geom_segment(aes(x = x, 
#                                            y = y, 
#                                            xend = xend, 
#                                            yend = yend),
#                                        data = segment.df)
#   p.segment
# }

# Calculate pseudotime using pca result
MyPca2Pseudotime <- function(pca.res,
                             x.dim = 1,
                             y.dim = 2,
                             num_paths = 1,
                             reverse = FALSE,
                             root_cell = NULL){
  library(monocle)
  pca.coord <- rbind(pca.res$ind$coord,
                     pca.res$ind.sup$coord)
  adjusted_S <- pca.coord[,c(x.dim,y.dim)]
  
  dp <- as.matrix(dist(adjusted_S))
  gp <- graph.adjacency(dp, 
                        mode = "undirected", 
                        weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  next_node <<- 0
  res <- monocle:::pq_helper(dp_mst, 
                             use_weights = FALSE, 
                             root_node = root_cell)
  cc_ordering <- monocle:::extract_good_branched_ordering(res$subtree, 
                                                          res$root, 
                                                          dp, 
                                                          num_paths, 
                                                          reverse)[,c(1,3)]
  colnames(cc_ordering) <- c("SampleName",
                             "Pseudotime")
  row.names(cc_ordering) <- cc_ordering$SampleName
  cc_ordering$SampleName <- as.character(cc_ordering$SampleName)
  cc_ordering
}

# get backbone dataframe
MyPca2Backbone <- function(
  pca.res,
  sample.inf,
  pseudotime.df,
  x.dim = 1,
  y.dim = 2){
  
  library(ggplot2)
  pc.percent <- pca.res$eig$"percentage of variance"
  pca.coord <- pca.res$ind$coord
  pca.coord <- as.data.frame(pca.coord)
  row.names(sample.inf) <- sample.inf$SampleName
  pca.df <- cbind(x.pos = pca.coord[,x.dim],
                  y.pos = pca.coord[,y.dim],
                  Pseudotime = pseudotime.df[row.names(pca.coord),
                                             "Pseudotime"],
                  sample.inf[row.names(pca.coord),])
  
  
  adjusted_S <- pca.coord[,c(x.dim,y.dim)]
  dp <- as.matrix(dist(adjusted_S))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  
  
  diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, 
                                                         weights = NA)]$name))
  colnames(diam) <- c("SampleName")
  diam <- arrange(merge(pca.df, 
                        diam, 
                        by.x = "SampleName", 
                        by.y = "SampleName"), 
                  Pseudotime)
  diam <- rename(diam, c(x.pos = "x.pos.backbone", 
                         y.pos = "y.pos.backbone"))
  diam
}


# get kegg gene list
# MyGetKeggEns("mmu04110")
MyGetKeggEns <- function(pathway_map){
  library(KEGGREST)
  library(org.Mm.eg.db)
  kegg.data <- keggGet(pathway_map)[[1]]
  gene_describe <- kegg.data$GENE
  gene <- gene_describe[seq(1,length(gene_describe),2)]
  ens <- mget(gene, 
              org.Mm.egENSEMBL,
              ifnotfound = NA)
  ens <- unlist(ens)
  ens <- ens[!is.na(ens)]
  ens <- unique(ens)
  return(ens)
}

# get boxplot outlier data and name
MyGetOutlier <- function(box.data,
                         box.name,
                         ...){
  box.res <- boxplot(box.data,
                     ...)
  outlier.data <- box.data[box.data %in% box.res$out]
  names(outlier.data) <- box.name[box.data %in% box.res$out]
  outlier.data <- sort(outlier.data,
                       decreasing = T)
  return(outlier.data)
}

MyGetTrend2 <- function(data.pseudotime,
                       data.rpkm,
                       min.rpkm = 1,
                       trend.formula = "log2rpkm ~ sm.ns(pseudotime, df=3)"){
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

# get trend data
MyGetTrend <- function(data.pseudotime,
                       data.rpkm,
                       min.rpkm = 1,
                       trend.formula = "log2rpkm ~ sm.ns(pseudotime, df=3)"){
  data <- data.frame(pseudotime = data.pseudotime,
                     rpkm = data.rpkm)
  data$log2rpkm <- log2(data$rpkm + min.rpkm)
  expectation <- tryCatch({vg <- suppressWarnings(vgam(formula = as.formula(trend.formula), 
                                                       family = tobit(), 
                                                       data = data, 
                                                       maxit = 30, 
                                                       checkwz = FALSE))
  res <- 2 ^ (predict(vg, 
                      type = "response")) - min.rpkm
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



# Box BeeSwarm
MyBoxBeeSwarm <- function(x,
                          y,
                          col,
                          pch = 20,
                          outline = F,
                          ...){
  library(beeswarm)
  boxplot(y ~ x,
          outline = outline,
          ...)
  beeswarm(y ~ x,
           col = col,
           pch = pch,
           add = TRUE,
           ...)
}


# Install Bioconductor Package
MyBioconductor <- function(package.name){
  #options(download.file.method="libcurl", url.method="libcurl")
  #source("https://www.bioconductor.org/biocLite.R")
  source("http://bioconductor.org/biocLite.R")
  biocLite(package.name)
}

# plot color
MyPlotColor <- function(palette,
                        point.count){
  plot(1:point.count,
       1:point.count,
       col = palette,
       pch = 20,
       cex = 5,
       xlab = "",
       ylab = "")
}

# order sample by PCA
MyPcaOrder <- function(pca.res,
                       order.pc = 1,
                       order.reverse = F){
  pc.coord <- pca.res$ind$coord[,order.pc]
  pc.coord <- sort(pc.coord)
  if(order.reverse){
    pc.coord <- rev(pc.coord)
  }
  sn.order <- names(pc.coord)
}

# Order gene by PCA
MyPcaDimd <- function(pca.res,
                      order.pc = 1,
                      order.reverse = F,
                      proba.cutoff = 1){
  dimd.data <- dimdesc(pca.res, 
                       axes = order.pc,
                       proba = proba.cutoff)
  dimd.data <- as.data.frame(dimd.data,
                             stringsAsFactors = F)
  dimd.data <- cbind(ID = rownames(dimd.data),
                     dimd.data,
                     stringsAsFactors = F)
  colnames(dimd.data) <- c("ID",
                           "correlation",
                           "p.value") 
  
  if(order.reverse){
    dimd.data <- dimd.data[rev(rownames(dimd.data)),]
  }
  return(dimd.data)
}



# operate GenomicRanges meta column
MyBedScoreAnd <- function(bdg1,
                          bdg2,
                          score.col = "TagCount",
                          operation = "+"){
  library(GenomicRanges)
  bdg.res <- disjoin(c(bdg1,
                       bdg2))
  mcols(bdg.res)[,score.col] <- 0
  overlap1.res <- findOverlaps(bdg.res,bdg1)
  mcols(bdg.res)[overlap1.res@from,score.col] <- mcols(bdg.res)[overlap1.res@from,score.col] + mcols(bdg1)[overlap1.res@to,score.col]
  overlap2.res <- findOverlaps(bdg.res,bdg2)
  if(operation == "+"){
    mcols(bdg.res)[overlap2.res@from,score.col] <- mcols(bdg.res)[overlap2.res@from,score.col] + mcols(bdg2)[overlap2.res@to,score.col]
  }
  if(operation == "-"){
    mcols(bdg.res)[overlap2.res@from,score.col] <- mcols(bdg.res)[overlap2.res@from,score.col] - mcols(bdg2)[overlap2.res@to,score.col]
  }
  return(bdg.res)
}


# intersectBed -v -a bed1 -b bed2
MyBedGrepV <- function(bed1,
                       bed2){
  library(GenomicRanges)
  bed1[countOverlaps(bed1,
                     bed2) == 0]
}

# intersectBed -u -a bed1 -b bed2
MyBedGrepU <- function(bed1,
                       bed2){
  library(GenomicRanges)
  bed1[countOverlaps(bed1,
                     bed2) > 0]
}

# reduce bedGraph
MyBdgReduce0 <- function(bdg.raw,
                         score.col = "TPM",
                         min.cutoff = 0){
  mcols(bdg.raw)[,score.col][mcols(bdg.raw)[,score.col] < min.cutoff] <- min.cutoff
  
  bdg.min <- bdg.raw[mcols(bdg.raw)[,score.col] == min.cutoff]
  bdg.not_min <- bdg.raw[mcols(bdg.raw)[,score.col] != min.cutoff]
  
  bdg.min.reduce <- reduce(bdg.min)
  mcols(bdg.min.reduce)[,score.col] <- min.cutoff
  
  bdg.processed <- c(bdg.min.reduce,
                     bdg.not_min)
  bdg.processed <- sort(bdg.processed)
  
  return(bdg.processed)
  
}

# setting baseline
MyBdgReduceBaseline <- function(bdg.raw,
                                score.col = "td",
                                min.cutoff = 0){
  bed.size <- range(bdg.raw)
  bdg.not_min <- bdg.raw[mcols(bdg.raw)[,score.col] > min.cutoff]
  
  bdg.min <- disjoin(c(bed.size,
                       granges(bdg.not_min)))
  bdg.min <- MyBedGrepV(bdg.min,
                        bdg.not_min)
  
  mcols(bdg.min)$value <- min.cutoff
  colnames(mcols(bdg.min)) <- score.col
  
  
  
  bdg.processed <- c(bdg.min,
                     bdg.not_min)
  bdg.processed <- sort(bdg.processed)
  
  return(bdg.processed)
}

# write bed file
MyWriteBed <- function(data,
                       output.file,
                       output.col = NULL,
                       is.gz = F,
                       zero.based = T,
                       round.col = NULL,
                       round.digits = 4,
                       sample.name = NULL,
                       colors = "0,0,0",
                       col.names = F,
                       ...){
  data <- as.data.frame(data)
  data <- data[,c("seqnames",
                  "start",
                  "end",
                  output.col)
               ]
  data[,round.col] <- round(data[,round.col],round.digits)
  if(zero.based){
    data[,"start"] <- data[,"start"] - 1
  }
  
  if(is.gz){
    output.pipe <- gzfile(output.file,
                          "w")
  }else{
    output.pipe <- output.file
  }
  
  if(!is.null(sample.name)){
    ucsc.txt <- paste('track type=bedGraph name="',sample.name,'" description="',sample.name,'" color="',colors,'" visibility=full','\t', 'yLineOnOff=on','\t', 'autoScale=on','\t', 'yLineMark="0.0"','\t', 'alwaysZero=on','\t', 'graphType=bar','\t', 'maxHeightPixels=128:75:11','\t', 'windowingFunction=maximum','\t', 'smoothingWindow=off','\t',
                      sep = "")
    cat(ucsc.txt,
        file = output.pipe)
    options(scipen = 10)
    write.table(data,
                output.pipe,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t",
                append = T,
                col.names = col.names,
                ...)
    options(scipen = 0)
  }else{
    options(scipen = 10)
    write.table(data,
                output.pipe,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t",
                col.names = col.names,
                ...)
    options(scipen = 0)
  }
  if(is.gz){
    close(output.pipe)
  }
}


# merge rep peak
MyBedMergeRep <- function(rep1.bed,
                          rep2.bed,
                          min.gapwidth = 1L){
  merge.bed <- reduce(c(rep1.bed,
                        rep2.bed),
                      min.gapwidth = min.gapwidth)
  merge.bed <- MyBedGrepU(merge.bed,
                          rep1.bed)
  merge.bed <- MyBedGrepU(merge.bed,
                          rep2.bed)
  return(merge.bed)
}

# ScoreMatrixBin: IP - input
MyScoreMatrixBin <- function(windows,
                             ip.target,
                             input.target,
                             windows.resize = 4000,
                             bin.num = 100,
                             weight.col = "TPM",
                             score.min = NULL,
                             score.max = NULL,
                             is.log = TRUE,
                             ...){
  
  MyAddExtendChr <- function(windows.bed,
                             target.bed,
                             weight.col){
    windows.chr.len <- elementNROWS(coverage(windows.bed))
    target.chr.len <- elementNROWS(coverage(target.bed))
    
    
    add.chr <- names(windows.chr.len)[!names(windows.chr.len) %in% names(target.chr.len)]
    if(length(add.chr) > 0){
      add.chr.bed <- GRanges(seqnames = Rle(add.chr),
                             ranges = IRanges(start = 1,
                                              end = windows.chr.len[add.chr]))
    }else{
      add.chr.bed <- GRanges()
    }
    
    common.chr <- names(windows.chr.len)[names(windows.chr.len) %in% names(target.chr.len)]
    extend.chr <- common.chr[windows.chr.len[common.chr] > target.chr.len[common.chr]]
    if(length(extend.chr) > 0){
      extend.chr.bed <- GRanges(seqnames = Rle(extend.chr),
                                ranges = IRanges(start = target.chr.len[extend.chr] + 1,
                                                 end = windows.chr.len[extend.chr]))
    }else{
      extend.chr.bed <- GRanges()
    }
    add.extend.chr.bed <- c(add.chr.bed,
                            extend.chr.bed)
    mcols(add.extend.chr.bed)[,weight.col] <- 0
    
    return(add.extend.chr.bed)
  }
  
  if(!is.null(windows.resize)){
    windows = resize(windows,
                     windows.resize,
                     fix = "center")
  }
  
  ip.target.add.extend <- MyAddExtendChr(windows,
                                         ip.target,
                                         weight.col)
  sm.ip = ScoreMatrixBin(target = c(ip.target,
                                    ip.target.add.extend),
                         windows = windows,
                         bin.num = bin.num,
                         weight.col = weight.col,
                         ...)
  input.target.add.extend <- MyAddExtendChr(windows,
                                            input.target,
                                            weight.col)
  sm.input = ScoreMatrixBin(target = c(input.target,
                                       input.target.add.extend),
                            windows = windows,
                            bin.num = bin.num,
                            weight.col = weight.col,
                            ...)
  
  sm <- sm.ip - sm.input
  sm@.Data <- sm.ip@.Data - sm.input@.Data
  if(!is.null(score.min)){
    sm@.Data[sm@.Data < score.min] <- score.min
  }
  if(!is.null(score.max)){
    sm@.Data[sm@.Data > score.max] <- score.max
  }
  if(is.log){
    sm <- log2(sm)
  }
  return(sm)
}



# modified function multiHeatMatrix
MyMultiHeatMatrix <- function (sml, 
                               grid = TRUE, 
                               col = NULL, 
                               xcoords = NULL, 
                               group = NULL, 
                               group.col = NULL,
                               order = FALSE, 
                               user.order = FALSE, 
                               winsorize = c(0, 100),
                               clustfun = FALSE, 
                               clust.matrix = NULL, 
                               column.scale = TRUE, 
                               matrix.main = NULL,
                               common.scale = FALSE, 
                               legend = TRUE, 
                               legend.name = NULL, 
                               cex.legend = 0.8,
                               xlab = NULL,
                               cex.lab = 1.2, 
                               cex.main = 1.5, 
                               cex.axis = 1.2,
                               newpage = TRUE,
                               lwd.heatmap = 10,
                               lwd.legend = 5,
                               col.range = NULL) {
  library(genomation)
  
  
  MyHeatLegendX <- function (min, max, cols, legend.name, main = TRUE, cex.legend = 1, 
                             cex.lab = 1, hjust = 0, vjust = 0) 
  {
    vals = seq(min, max, length.out = 100)
    rng <- range(vals, na.rm = TRUE)
    m <- (vals - min)/(max - min)
    rasta = rgb(colorRamp(cols)(m), maxColorValue = 255)
    grid.raster(matrix((rasta), nrow = 1), interpolate = FALSE, 
                height = unit(1, "npc"), width = unit(1, "npc"))
    # label = pretty(c(min, max), n = 5)
    label = c(min, max)
    at = seq(0, 1, length.out = length(label))
    # grid.xaxis(at = at, 
    #            label = label,
    #            main = main, 
    #            edits = gEdit("labels", 
    #                          hjust = hjust,
    #                          vjust = vjust), 
    #            gp = gpar(cex = cex.legend))
    grid.xaxis(at = at, 
               label = label,
               main = main, 
               edits = gEdit("labels", 
                             hjust = hjust,
                             vjust = vjust), 
               gp = gpar(cex = cex.legend,
                         lwd = 0))
    my.y = -3
    grid.text(legend.name, y = unit(my.y, "npc"), gp = gpar(cex = cex.lab,
                                                            lwd = 0))
  }
  
  MyGridHeat<- function (mat, col, xcoords, xlab, cex.lab, cex.axis, angle = 0, 
                         hjust = 0, vjust = 0, rng = NULL) 
  {
    MyConvertToColors <- function (mat, cols, rng = NULL) 
    {
      if (is.null(rng)) {
        rng <- range(mat, na.rm = TRUE)
      }
      nc = length(cols)
      if (diff(rng) == 0) 
        rng <- if (rng[1L] == 0) {
          c(-1, 1)
        }
      else {
        rng[1L] + c(-0.4, 0.4) * abs(rng[1L])
      }
      mat <- (mat - rng[1L])/diff(rng)
      zi <- floor((nc - 1e-05) * mat + 1e-07)
      # zi[zi < 0 | zi >= nc] <- NA
      zi[zi < 0 ] <- 0
      zi[zi >= nc ] <- nc
      zc <- cols[zi + 1L]
      dim(zc) <- dim(mat)
      return(zc)
    }
    # mat2 = .convertToColors(mat, col, rng)
    mat2 = MyConvertToColors(mat, col, rng)
    ras = grid.raster(mat2, interpolate = FALSE, width = unit(1, 
                                                              "npc"), height = unit(1, "npc"))
    # at = seq(0, 1, length.out = 5)
    at = c(0,1)
    # label = seq(min(xcoords), max(xcoords), length.out = 5)
    label = c(min(xcoords), max(xcoords))
    # ax = grid.xaxis(at = at, 
    #                 label = formatC(label, 
    #                                 digits = 4, 
    #                                 format = "g"), 
    #                 edits = gEdit("labels", 
    #                               rot = angle, 
    #                               hjust = hjust, 
    #                               vjust = vjust), 
    #                 gp = gpar(cex = cex.axis))
    ax = grid.xaxis(at = at, 
                    label = formatC(label, 
                                    digits = 4, 
                                    format = "g"), 
                    edits = gEdit("labels", 
                                  rot = angle, 
                                  hjust = hjust, 
                                  vjust = vjust), 
                    gp = gpar(cex = cex.axis,
                              lwd = 0))
    grid.text(xlab, 
              y = unit(-2.5, "lines"), 
              gp = gpar(cex = cex.lab))
  }
  
  if (class(sml) != "ScoreMatrixList") {
    stop("'sml' is not a ScoreMatrix object\n")
  }
  if (is.list(col) & length(col) != length(sml)) {
    stop("'col' is a list and its length does not match the length of ScoreMatrixList\n")
    col = NULL
  }
  if (is.list(xcoords) & length(xcoords) != length(sml)) {
    stop("'xcoords' is a list and its length does not match the length of ScoreMatrixList\n")
    xcoords = NULL
  }
  if (!is.null(legend.name) & length(legend.name) > 1 & length(legend.name) != 
      length(sml)) {
    stop("'legend.name' should match the length of the 'sml' ", 
         "if it is not a length 1 vector\n")
  }
  else if (length(legend.name) == 1) {
    legend.name = rep(legend.name, length(sml))
  }
  if (!is.null(xlab) & length(xlab) > 1 & length(xlab) != length(sml)) {
    stop("'xlab' should match the length of the 'sml' ", 
         "if it is not a length 1 vector\n")
  }
  else if (length(xlab) == 1) {
    xlab = rep(xlab, length(sml))
  }
  if (length(unique(sapply(sml, nrow))) > 1) {
    warning("\nThe row numbers are different\n", "attempting to get common rows and to reorder rows\n", 
            "using 'intersectScoreMatrixList()'\n")
    sml = intersectScoreMatrixList(sml, reorder = TRUE)
  }
  if (!is.null(clust.matrix)) {
    clust.matrix = unique(clust.matrix)
    if (is.numeric(clust.matrix)) {
      if (length(clust.matrix) > length(sml) | max(clust.matrix) > 
          length(sml) | 0 %in% clust.matrix) {
        warning("\n'clust.matrix' vector shouldn't be longer than 'sml'\n", 
                "and shouldn't have greater values that length of 'sml'\n", 
                "clustering all matrices\n")
        clust.matrix = NULL
      }
    }
    else {
      if (length(clust.matrix) > length(sml) | sum(clust.matrix %in% 
                                                   names(sml)) == length(clust.matrix)) {
        warning("\n'clust.matrix' vector shouldn't be longer than 'sml'\n", 
                "and should contain names of matrices of 'sml'\n", 
                "clustering all matrices\n")
        clust.matrix = NULL
      }
    }
  }
  mat.list = lapply(sml, function(x) x@.Data)
  group.vector = NULL
  group.names = NULL
  if (winsorize[2] < 100 | winsorize[1] > 0) {
    mat.list = lapply(mat.list, function(x) .winsorize(x@.Data, 
                                                       winsorize))
  }
  if (!identical(clustfun, FALSE) | order) {
    if (!is.null(clust.matrix)) {
      mat2 = do.call("cbind", mat.list[clust.matrix])
    }
    else {
      mat2 = do.call("cbind", mat.list)
    }
    if (column.scale) {
      mat2 = scale(mat2)
      mat2[is.nan(mat2)] = 0
    }
  }
  if (!identical(clustfun, FALSE)) {
    if (any(is.na(mat2))) {
      mat3 = impute.knn(mat2, k = 10, rowmax = 0.5, colmax = 0.8, 
                        maxp = 1500)$data
      clu = clustfun(mat3)
    }
    else {
      clu = clustfun(mat2)
    }
    group.vector = clu
    mat.list = lapply(mat.list, function(x) x[order(group.vector), 
                                              ])
    group.vector = group.vector[order(group.vector)]
    if (order) {
      g.factor = factor(group.vector, levels = unique(group.vector))
      my.order = order(group.vector, -rowSums(mat2, na.rm = TRUE))
      mat.list = lapply(mat.list, function(x) x[my.order, 
                                                ])
      group.vector = group.vector[my.order]
    }
    group.names = unique(group.vector)
  }
  if (!is.null(group) & identical(clustfun, FALSE)) {
    if (is.list(group)) {
      win.numbs = (lapply(group, function(x) unique(x)))
      win.vec = unlist(win.numbs)
      if (any(table(unlist(win.numbs)) > 1)) {
        stop("'group' containing a list must not have duplicated numbers\n")
      }
      row.ids = rownames(mat.list[[1]])
      group.vector = rep(0, length(row.ids))
      for (i in 1:length(win.numbs)) {
        group.vector[row.ids %in% win.numbs[[i]]] = i
      }
      if (!is.null(names(group))) {
        group.names = names(group)
      }
      if (all(group.vector == 0)) {
        stop("None of the elements in 'group' are a part of rownames(mat) \n")
      }
      if (any(group.vector == 0)) {
        warning("Number of elements in 'group' argument is less then nrow(mat) \n", 
                "Dropping rows from 'mat' that are not contained in 'group'\n")
        mat.list = lapply(mat.list, function(x) x[group.vector > 
                                                    0, ])
        group.vector = group.vector[group.vector > 0]
      }
    }
    else if (is.factor(group)) {
      if (length(group) != nrow(mat.list[[1]])) {
        stop("'group' is a factor, and its length should be equal to nrow(mat)\n")
      }
      group = factor(as.character(group), levels = as.character(unique(group)))
      group.names = levels(group)
      levels(group) = 1:length(levels(group))
      group.vector = as.numeric(group)
    }
    else {
      stop("'group' must be a factor or a list\n")
    }
    mat.list = lapply(mat.list, function(x) x[order(group.vector), 
                                              ])
    group.vector = group.vector[order(group.vector)]
    if (order) {
      mat2 = do.call("cbind", mat.list)
      if (column.scale) {
        mat2 = scale(mat2)
        mat2[is.nan(mat2)] = 0
      }
      my.order = order(group.vector, -rowSums(mat2, na.rm = TRUE))
      mat.list = lapply(mat.list, function(x) x[my.order, 
                                                ])
      group.vector = group.vector[my.order]
    }
  }
  else if (order & identical(clustfun, FALSE)) {
    mat2 = do.call("cbind", mat.list)
    if (column.scale) {
      mat2 = scale(mat2)
      mat2[is.nan(mat2)] = 0
    }
    order.vector = rep(1, nrow(mat2))
    names(order.vector) = rownames(mat2)
    my.order = order(-rowSums(mat2, na.rm = TRUE))
    mat.list = lapply(mat.list, function(x) x[my.order, ])
    order.vector = order.vector[my.order]
  }
  if (!identical(user.order, FALSE) & !is.null(group.names)) {
    if (length(user.order) != length(group.names)) {
      warning(paste0("length of 'user.order' vector (", 
                     length(user.order), ") should be the the same as number of clusters (", 
                     length(group.names), "). Skipping it..\n"))
    }
    else {
      gv.fac.ord <- sort(factor(group.vector, levels = user.order))
      group.vector <- group.vector[names(gv.fac.ord)]
      group.names <- user.order
    }
  }
  else if (!identical(user.order, FALSE) & is.null(group.names)) {
    warning("There are no groups or clusters to order. Skipping it..")
  }
  if (!grid) {
    plot.new()
    # vps <- baseViewports()
    vps <- gridBase::baseViewports()
    pushViewport(vps$figure)
  }
  else {
    if (newpage) 
      grid.newpage()
  }
  if (is.null(matrix.main) & !is.null(names(sml)) & class(sml) == 
      "ScoreMatrixList") {
    matrix.main = names(sml)
  }
  else if (!is.null(matrix.main) & length(matrix.main) != length(sml)) {
    warning("'matrix.main' length does not match to the 'sml' length\n", 
            "setting it to NULL")
    matrix.main = NULL
  }
  l.sml = length(mat.list)
  hw = 40 * (1 - (convertX(unit(4, "lines"), "npc", valueOnly = TRUE)))/(l.sml * 
                                                                           44 + 7)
  hh = 0.7
  hy = 0.6
  if (!is.null(group.vector)) {
    sideVp <- viewport(width = unit(hw/8, "npc"), height = unit(hh, 
                                                                "npc"), x = convertX(unit(4, "lines"), "npc"), y = unit(hy, 
                                                                                                                        "npc"), just = "left")
    pushViewport(sideVp)
    grid.rect()
    .rowSideCol(group.vector, group.names = group.names, 
                group.col = group.col, cex.lab = cex.lab)
    popViewport()
  }
  heat.startCoord = convertX(unit(4, "lines"), "npc", valueOnly = TRUE) + 
    (hw/8) + (hw/20)
  common.range = NULL
  if (common.scale) {
    common.range = range(mat.list, na.rm = TRUE)
  }
  for (i in 1:length(mat.list)) {
    if (is.list(xcoords)) {
      cxcoords = xcoords[[i]]
    }
    else if (is.vector(xcoords) | is.null(xcoords)) {
      cxcoords = xcoords
    }
    else if (!is.null(xcoords)) {
      warning("xcoords should be a vector or a list or NULL,nothing else!!\n", 
              "setting xcoords to NULL")
      cxcoords = NULL
    }
    if (length(cxcoords) == 2) {
      cxcoords = seq(cxcoords[1], cxcoords[2], length.out = ncol(mat.list[[i]]))
    }
    if (is.list(col)) {
      ccol = col[[i]]
    }
    else if (is.vector(col) | is.null(col)) {
      ccol = col
    }
    else if (!is.null(col)) {
      warning("col should be a vector or a list or NULL,nothing else!!\n", 
              "setting colors to default")
      ccol = NULL
    }
    if (!is.null(cxcoords)) {
      if (length(cxcoords) != ncol(mat.list[[i]])) 
        stop("xcoords has wrong length: ", length(cxcoords), 
             " \n", " it should be equal to the number of columns of ScoreMatrix\n", 
             " which is", ncol(mat.list[[i]]), "\n")
    }
    else {
      cxcoords = 1:ncol(mat.list[[i]])
    }
    if (is.null(ccol)) 
      ccol = .jets(100)
    heatVp <- viewport(width = unit(hw, "npc"), height = unit(hh, 
                                                              "npc"), x = unit(heat.startCoord, "npc"), y = unit(hy, 
                                                                                                                 "npc"), just = "left")
    if (common.scale) {
      rng = common.range
    }else if (is.list(col.range)){
      rng = col.range[[i]]
    }else if (!is.null(col.range)){
      rng = col.range
    }else {
      rng = range(mat.list[[i]], na.rm = TRUE)
    }
    pushViewport(heatVp)
    # grid.rect()
    grid.rect(gp = gpar(lwd = lwd.heatmap))
    # .gridHeat(mat.list[[i]], ccol, cxcoords, xlab[i], cex.lab, 
    #           cex.axis, angle = 60, hjust = 0.6, vjust = -0.5, 
    #           rng = common.range)
    MyGridHeat(mat.list[[i]],   
               ccol,
               cxcoords, 
               xlab[i], 
               cex.lab, 
               cex.axis,
               angle = 60, 
               hjust = 0.6, 
               vjust = -0.5, 
               # rng = common.range)
               rng = rng)
    popViewport()
    if (legend) {
      legendVp <- viewport(width = unit(hw * 0.7, "npc"), 
                           height = unit(0.5, "lines"), x = unit(heat.startCoord + 
                                                                   hw * 0.15, "npc"), y = unit(0.1, "npc"), just = "left")
      pushViewport(legendVp)
      # grid.rect()
      grid.rect(gp = gpar(lwd = lwd.legend))
      
      # .heatLegendX(min = rng[1], max = rng[2], ccol, legend.name[i], 
      #              main = TRUE, cex.legend, cex.lab, vjust = -0.5, 
      #              hjust = 0.5)
      MyHeatLegendX(min = rng[1], max = rng[2], ccol, legend.name[i], 
                    main = TRUE, cex.legend, cex.lab, vjust = -0.5, 
                    hjust = 0.5)
      popViewport()
    }
    if (!is.null(matrix.main)) {
      title.y = unit(0.96, "npc")
      grid.text(matrix.main[i], y = title.y, x = unit(heat.startCoord + 
                                                        hw * 0.5, "npc"), gp = gpar(cex = cex.main), 
                just = "bottom")
    }
    heat.startCoord = heat.startCoord + hw + (hw/20)
  }
  if (!grid) {
    popViewport()
  }
  if (!identical(clustfun, FALSE) | !is.null(group.vector)) {
    return(invisible(group.vector))
  }
  else if (order & is.null(group.vector)) {
    return(invisible(order.vector))
  }
}


# modified function readTranscriptFeatures
MyReadTranscriptFeatures <- function(location,
                                     remove.unusual = FALSE,
                                     up.flank = 1000,
                                     down.flank = 1000,
                                     unique.prom = FALSE){
  library(genomation)
  # readBed6
  message('Reading the table...\r')
  bed=genomation:::readTableFast(location,header=FALSE,skip="auto")                    
  if(remove.unusual)
    bed=bed[grep("_", as.character(bed[,1]),invert=TRUE),]
  
  # introns
  message('Calculating intron coordinates...\r')
  introns    = convertBed2Introns(bed)
  # exons
  message('Calculating exon coordinates...\r')
  exons    = convertBed2Exons(bed)
  
  # get the locations of TSSes
  message('Calculating TSS coordinates...\r')
  tss=bed
  tss[,2] = tss[,2] + 1
  #  + strand
  tss[tss$V6=="+",3] = tss[tss$V6=="+",2]
  #  - strand
  tss[tss$V6=="-",2]=tss[tss$V6=="-",3]
  
  tssg = GRanges(seqnames=as.character(tss$V1),
                 ranges=IRanges(start=tss$V2, end=tss$V3),
                 strand=as.character(tss$V6),
                 score=rep(0,nrow(tss)),
                 name=tss$V4)
  
  message('Calculating promoter coordinates...\r')
  # get the locations of promoters
  # + strand
  # bed[bed$V6=="+",3] = bed[bed$V6=="+",2] + down.flank 
  # bed[bed$V6=="+",2] = bed[bed$V6=="+",2] - up.flank
  bed[bed$V6=="+",3] = bed[bed$V6=="+",2] + down.flank + 1
  bed[bed$V6=="+",2] = bed[bed$V6=="+",2] - up.flank + 1
  
  
  #  - strand
  bed[bed$V6=="-",2] = bed[bed$V6=="-",3] - down.flank
  bed[bed$V6=="-",3] = bed[bed$V6=="-",3] + up.flank
  
  
  if(!unique.prom){
    prom.df = (bed[,c(1,2,3,4,6)])
    prom = GRanges(seqnames=as.character(prom.df$V1),
                   ranges = IRanges(start=prom.df$V2, end=prom.df$V3),
                   strand = as.character(prom.df$V6),
                   score=rep(0,nrow(prom.df)),
                   name=prom.df$V4)
  }else{
    prom.df = unique(bed[,c(1,2,3,6)])
    prom = GRanges(seqnames=as.character(prom.df$V1),
                   ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
                   strand=as.character(prom.df$V6),
                   score=rep(0,nrow(prom.df)),
                   name=rep(".",nrow(prom.df)) )
  }
  
  message('Outputting the final GRangesList...\r\n')
  GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
}


# compare peak associated gene expression
MyBoxComp <- function(bed.A,
                      bed.B,
                      mcol.gene,
                      mcol.rpkm.A,
                      mcol.rpkm.B,
                      name.A,
                      name.B,
                      col.A = "#5B9CD6",
                      col.B = "#ED7D31"){
  rpkm.A <- unique(mcols(bed.A)[,c(mcol.gene,
                                   mcol.rpkm.A)])[,mcol.rpkm.A]
  rpkm.A <- rpkm.A[!is.na(rpkm.A)]
  rpkm.B <- unique(mcols(bed.B)[,c(mcol.gene,
                                   mcol.rpkm.B)])[,mcol.rpkm.B]
  rpkm.B <- rpkm.B[!is.na(rpkm.B)]
  
  wilcox.p.value <- min(wilcox.test(rpkm.A,
                                    rpkm.B,
                                    alternative = "less")$p.value,
                        wilcox.test(rpkm.A,
                                    rpkm.B,
                                    alternative = "greater")$p.value)
  wilcox.p.value <- signif(wilcox.p.value,3)
  
  rpkm.plot <- list(rpkm.A,
                    rpkm.B)
  names(rpkm.plot) <- c(name.A,
                        name.B)
  
  boxplot(rpkm.plot,
          outline = F,
          col = c(col.A,
                  col.B),
          main = paste("p-value =",
                       wilcox.p.value))
}



# compare peak associated gene expression(bivalent)
MyBoxComp3 <- function(bed.A,
                       bed.B,
                       bed.C,
                       mcol.gene,
                       mcol.rpkm.A,
                       mcol.rpkm.B,
                       mcol.rpkm.C,
                       name.A,
                       name.B,
                       name.C,
                       col.A = "#E913E2",
                       col.B = "#EA9225",
                       col.C = "#22F34E"){
  rpkm.A <- unique(mcols(bed.A)[,c(mcol.gene,
                                   mcol.rpkm.A)])[,mcol.rpkm.A]
  rpkm.A <- rpkm.A[!is.na(rpkm.A)]
  rpkm.B <- unique(mcols(bed.B)[,c(mcol.gene,
                                   mcol.rpkm.B)])[,mcol.rpkm.B]
  rpkm.B <- rpkm.B[!is.na(rpkm.B)]
  rpkm.C <- unique(mcols(bed.C)[,c(mcol.gene,
                                   mcol.rpkm.C)])[,mcol.rpkm.C]
  rpkm.C <- rpkm.C[!is.na(rpkm.C)]
  
  wilcox.p.value.AB <- min(wilcox.test(rpkm.A,
                                       rpkm.B,
                                       alternative = "less")$p.value,
                           wilcox.test(rpkm.A,
                                       rpkm.B,
                                       alternative = "greater")$p.value)
  wilcox.p.value.AB <- signif(wilcox.p.value.AB,3)
  
  wilcox.p.value.BC <- min(wilcox.test(rpkm.B,
                                       rpkm.C,
                                       alternative = "less")$p.value,
                           wilcox.test(rpkm.B,
                                       rpkm.C,
                                       alternative = "greater")$p.value)
  wilcox.p.value.BC <- signif(wilcox.p.value.BC,3)
  
  wilcox.p.value.AC <- min(wilcox.test(rpkm.A,
                                       rpkm.C,
                                       alternative = "less")$p.value,
                           wilcox.test(rpkm.A,
                                       rpkm.C,
                                       alternative = "greater")$p.value)
  wilcox.p.value.AC <- signif(wilcox.p.value.AC,3)
  
  rpkm.plot <- list(rpkm.A,
                    rpkm.B,
                    rpkm.C)
  names(rpkm.plot) <- c(name.A,
                        name.B,
                        name.C)
  
  boxplot(rpkm.plot,
          outline = F,
          col = c(col.A,
                  col.B,
                  col.C),
          main = paste(wilcox.p.value.AB,
                       wilcox.p.value.BC,
                       wilcox.p.value.AC))
}


# Add Extend Chr
MyAddExtendChr <- function(windows.bed,
                           target.bed,
                           weight.col){
  windows.chr.len <- elementNROWS(coverage(windows.bed))
  target.chr.len <- elementNROWS(coverage(target.bed))
  
  
  add.chr <- names(windows.chr.len)[!names(windows.chr.len) %in% names(target.chr.len)]
  if(length(add.chr) > 0){
    add.chr.bed <- GRanges(seqnames = Rle(add.chr),
                           ranges = IRanges(start = 1,
                                            end = windows.chr.len[add.chr]))
  }else{
    add.chr.bed <- GRanges()
  }
  
  common.chr <- names(windows.chr.len)[names(windows.chr.len) %in% names(target.chr.len)]
  extend.chr <- common.chr[windows.chr.len[common.chr] > target.chr.len[common.chr]]
  if(length(extend.chr) > 0){
    extend.chr.bed <- GRanges(seqnames = Rle(extend.chr),
                              ranges = IRanges(start = target.chr.len[extend.chr] + 1,
                                               end = windows.chr.len[extend.chr]))
  }else{
    extend.chr.bed <- GRanges()
  }
  add.extend.chr.bed <- c(add.chr.bed,
                          extend.chr.bed)
  mcols(add.extend.chr.bed)[,weight.col] <- 0
  
  return(add.extend.chr.bed)
}

# calculate tag density
MyCalTd <- function(windows,
                    ip.target,
                    weight.col,
                    ...){
  
  MyAddExtendChr <- function(windows.bed,
                             target.bed,
                             weight.col){
    windows.chr.len <- elementNROWS(coverage(windows.bed))
    target.chr.len <- elementNROWS(coverage(target.bed))
    
    
    add.chr <- names(windows.chr.len)[!names(windows.chr.len) %in% names(target.chr.len)]
    if(length(add.chr) > 0){
      add.chr.bed <- GRanges(seqnames = Rle(add.chr),
                             ranges = IRanges(start = 1,
                                              end = windows.chr.len[add.chr]))
    }else{
      add.chr.bed <- GRanges()
    }
    
    common.chr <- names(windows.chr.len)[names(windows.chr.len) %in% names(target.chr.len)]
    extend.chr <- common.chr[windows.chr.len[common.chr] > target.chr.len[common.chr]]
    if(length(extend.chr) > 0){
      extend.chr.bed <- GRanges(seqnames = Rle(extend.chr),
                                ranges = IRanges(start = target.chr.len[extend.chr] + 1,
                                                 end = windows.chr.len[extend.chr]))
    }else{
      extend.chr.bed <- GRanges()
    }
    add.extend.chr.bed <- c(add.chr.bed,
                            extend.chr.bed)
    mcols(add.extend.chr.bed)[,weight.col] <- 0
    
    return(add.extend.chr.bed)
  }
  
  
  ip.target.add.extend <- MyAddExtendChr(windows,
                                         ip.target,
                                         weight.col)
  
  td <- as.numeric(ScoreMatrixBin(target = c(ip.target,
                                             ip.target.add.extend), 
                                  windows = windows,
                                  bin.num = 1,
                                  weight.col = weight.col))
  return(td)
}


# matrix to heatmap
MyImage <- function(data.image,
                    cutoff.min,
                    cutoff.max,
                    color.image = rgb(colorRamp(c("white","red"))(seq(0,1,0.01)), maxColorValue = 255),
                    box.lwd = 1,
                    ...){
  
  data.image.m <- data.image
  data.image.m[data.image.m < cutoff.min] <- cutoff.min
  data.image.m[data.image.m > cutoff.max] <- cutoff.max
  
  image(data.image.m,
        col  = color.image,
        zlim = c(cutoff.min,cutoff.max),
        axes = F,
        ...)
  box(lwd = box.lwd)
}

# Get PCA coordinate from seurat object
MySeuratPca2Gg <- function(seurat.data,
                           sample.inf,
                           x.dim = 1,
                           y.dim = 2,
                           color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  pca.coord <- seurat.data@dr$pca@cell.embeddings
  pca.coord <- as.data.frame(pca.coord)
  
  row.names(sample.inf) <- sample.inf$SampleName
  pca.df <- cbind(x.pos = pca.coord[,x.dim],
                  y.pos = pca.coord[,y.dim],
                  sample.inf[row.names(pca.coord),])
  
  p <- ggplot(data = pca.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("PC",x.dim,sep = ""))
  p <- p + ylab(paste("PC",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}

# Get t-SNE coordinate from seurat object
MySeuratTsne2Gg <- function(seurat,
                            sample.inf,
                            x.dim = 1,
                            y.dim = 2,
                            color_manual = brewer.pal(9,"Set1")){
  library(ggplot2)
  tSNE.coord <- seurat@dr$tsne@cell.embeddings
  tSNE.coord <- as.data.frame(tSNE.coord)
  
  row.names(sample.inf) <- sample.inf$SampleName
  tSNE.df <- cbind(x.pos = tSNE.coord[,x.dim],
                   y.pos = tSNE.coord[,y.dim],
                   sample.inf[row.names(tSNE.coord),])
  
  p <- ggplot(data = tSNE.df, 
              mapping = aes(x = x.pos, 
                            y = y.pos,
                            label = SampleName))
  p <- p + xlab(paste("tSNE",x.dim,sep = ""))
  p <- p + ylab(paste("tSNE",y.dim,sep = ""))
  p <- p + guides(colour=guide_legend(title=NULL))
  p <- p + scale_color_manual(values = color_manual)
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  return(p)
}



# Seurat to plotly 3d scatter
MySeuratPca2Plotly3d <- function(seurat.data,
                                 x.dim = 1,
                                 y.dim = 2,
                                 z.dim = 3,
                                 ...){
  library(plotly)
  
  p <- plot_ly(x = seurat.data@pca.rot[,x.dim], 
               y = seurat.data@pca.rot[,y.dim], 
               z = seurat.data@pca.rot[,z.dim], 
               type = "scatter3d", 
               mode = "markers",
               ...)
  return(p)
  detach("package:plotly", 
         unload=TRUE)
}

MyCalSi <- function(data,
                    k,
                    graph = F){
  library(cluster)
  data.log2.relat <- as.matrix(log2(data + 1) / rowMax(log2(data + 1)))
  data.cor <- cor(data.log2.relat,
                  method = "s")
  data.dist <- as.dist(abs(1 - data.cor))
  data.hc <- hclust(data.dist, 
                    method = "ward.D2") 
  si <- silhouette(cutree(data.hc, k = k), data.dist)
  if(graph){
    plot(si, nmax = 80, cex.names = 0.5)
  }
  return(mean(si[,3]))
}
  
  MySCDE <- function(data.rc,
           sn.a,
           sn.b,
           n.cores = 1){
    
    library(scde)
    library(parallel)
    
    sg <- factor(c(rep("a",length(sn.a)),
                   rep("b",length(sn.b))), 
                 levels = c("a", "b"))
    names(sg) <- c(sn.a,
                   sn.b)
    
    cd <- clean.counts(data.rc[,c(sn.a,
                                  sn.b)], 
                       min.lib.size= 1, 
                       min.reads = 1,
                       min.detected = 1)
    
    o.ifm <- scde.error.models(counts = cd,
                               groups = sg,
                               n.cores = n.cores,
                               threshold.segmentation = TRUE,
                               save.crossfit.plots = FALSE,
                               save.model.plots = FALSE,
                               verbose = 1)
    
    valid.cells <- o.ifm$corr.a > 0
    table(valid.cells)
    
    o.ifm <- o.ifm[valid.cells, ]
    
    # estimate gene expression prior
    o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
    
    # define two groups of cells
    groups <- sg[rownames(o.ifm)]
    
    # run differential expression tests on all genes.
    ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
    ediff$p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
    ediff$p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
    ediff
  }


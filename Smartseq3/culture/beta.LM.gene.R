
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  print(seu)
  si.select <- seu.vitro.nor.list[[seu]]@meta.data
  pregScore=c()
  for (time in names(table(seu.vitro.nor.list[[seu]]$Type))) {
    regress=lm(-Embeddings(seu.vitro.nor.list[[seu]],'pca')[si.select$Type==time,2]~
                 Embeddings(seu.vitro.nor.list[[seu]],'pca')[si.select$Type==time,1])$coefficients[2]
    temp=-c(1,regress)%*%rbind(Embeddings(seu.vitro.nor.list[[seu]],'pca')[si.select$Type==time,1],
                               Embeddings(seu.vitro.nor.list[[seu]],'pca')[si.select$Type==time,2])
    # temp=temp-mean(temp)
    names(temp)=si.select[si.select$Type==time,]$SampleName
    pregScore=c(pregScore,temp)
  }
  seu.vitro.nor.list[[seu]]@meta.data[names(pregScore),'KOScore'] <- pregScore
  seu.vitro.nor.list[[seu]]$pseudotime <- seu.vitro.nor.list[[seu]]$KOScore
}
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  print(MyPseudotimebox(seu.vitro.nor.list[[seu]]@meta.data,time.colors = time.colors,size.point = 2)+theme(aspect.ratio = 1.5))
}
seu.vitro.nor.list$sm3.A485$KOScore <- -seu.vitro.nor.list$sm3.A485$KOScore
seu.vitro.nor.list$Nif$KOScore <- -seu.vitro.nor.list$Nif$KOScore

###########
preg.gene <- list()
gene.input.list <- list()
cl<-makeCluster(10)
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  gene.input.list[[seu]] <- rownames(MyGeneExp(seu.vitro.nor.list[[seu]]@assays$RNA@data,1,10))
  gene.input.list[[seu]] <- setdiff(gene.input.list[[seu]],pre.ambigous.sym)
  # preg.gene[[seu]] <- MyLM.parallel(seu.vitro.nor.list[[seu]]@assays$RNA@data[gene.input.list[[seu]],colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$KOScore)]],
  #                                   variable = sort(seu.vitro.nor.list[[seu]]$KOScore))
}

###########preg##########
dir.create('direct/20221120/path/')
preg.gene.filter <- list()
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  preg.gene.filter[[seu]] <- preg.gene[[seu]][-log10(preg.gene[[seu]])>=8]
}
sapply(preg.gene.filter, length)


preg.row.tree <- list()
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$KOScore)],"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
png(paste('direct/20221120/path/',seu,'doCol.koscore.pval8.pregTop.heatmap.png',sep = ''),2000,3000)
preg.row.tree[[seu]] <-
  MyHeatmap(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[names(preg.gene.filter[[seu]]),colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$KOScore)]],
            type = "log.row.zscore",
            hc.c.data.type = "log.row.relat",
            hc.r.data.type = "log.row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D2",
            r.hc.method = "ward.D2",
            ColSideColors = c1Color,
            ColSideColorsSize = 1.5,
            return.tree = "row",
            #Colv = 'none',
            dendrogram = 'both',
            graph = T)
dev.off()
}

labels(as.dendrogram(preg.row.tree$Nif)[[1]][[2]][[1]][[2]][[2]])


sn.DEG.list[['Nif']] <- list(c("Ctrl_DMSO","Preg_DMSO"),c("Ctrl_20um_Nif","Preg_20um_Nif"))
names(sn.DEG.list[['Nif']]) <- c('DMSO','Nif')

sn.DEG.list[['sm3.A485']] <- list(c("Ctrl_DMSO","Preg_DMSO"),c("Ctrl_A485","Preg_A485"))
names(sn.DEG.list[['sm3.A485']]) <- c('DMSO','A485')

sn.DEG.list[['sm3.A485.Ctrl']] <- list(c("Ctrl_DMSO"),c("Ctrl_A485"))
names(sn.DEG.list[['sm3.A485.Ctrl']]) <- c("Ctrl_DMSO","Ctrl_A485")

sn.DEG.list[['sm3.A485.add']] <- list(c("Ctrl_DMSO","Preg_DMSO"),c("Ctrl_A485","Preg_A485"))
names(sn.DEG.list[['sm3.A485.add']]) <- c('DMSO','A485')

sn.DEG.list[['sm3.A485.add.Ctrl']] <- list(c("Ctrl_DMSO"),c("Ctrl_A485"))
names(sn.DEG.list[['sm3.A485.add.Ctrl']]) <- c("Ctrl_DMSO","Ctrl_A485")



sn.DEG.list[['Nif.Ctrl']] <- list(c("Ctrl_DMSO"),c("Ctrl_20um_Nif"))
names(sn.DEG.list[['Nif.Ctrl']]) <- c("Ctrl_DMSO","Ctrl_20um_Nif")

sn.DEG.list[['Nif.A485.Ctrl']] <- list(c('Ctrl_DMSO'),c('Ctrl_A485_Nif'))
names(sn.DEG.list[['Nif.A485.Ctrl']]) <- c("Ctrl_DMSO",'Ctrl_A485_Nif')


seu <- 'sm3.A485.add'
seu.vitro.nor.list[[seu]]$Type_batch <- paste(seu.vitro.nor.list[[seu]]$Type,seu.vitro.nor.list[[seu]]$SeqDate,sep = '_')
sn.DEG.list[['sm3.A485.add.Ctrl']] <- list(c("Ctrl_DMSO_20221120"),c("Ctrl_A485_20221120"))
names(sn.DEG.list[['sm3.A485.add.Ctrl']]) <- c("Ctrl_DMSO","Ctrl_A485")


sn.DEG.list[['sm3.A485.addpregA485']] <- list(c("Ctrl_DMSO","Preg_DMSO"),c("Ctrl_A485","Preg_A485"))
names(sn.DEG.list[['sm3.A485.addpregA485']]) <- c('DMSO','A485')

sn.DEG.list[['sm3.A485.addpregA485.Ctrl']] <- list(c("Ctrl_DMSO"),c("Ctrl_A485"))
names(sn.DEG.list[['sm3.A485.addpregA485.Ctrl']]) <- c("Ctrl_DMSO","Ctrl_A485")


KO.DEG.nor.list <- list()
seu.vitro.nor.list <- readRDS('seu.vitro.nor.list.rds')
for(seu in names(seu.vitro.nor.list)[c(1)]){
  seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$Type)
  #seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$Type_batch)
  
  for(tmp in names(sn.DEG.list)[9]){
    KO.DEG.nor.list[[paste(names(sn.DEG.list[[tmp]])[1],names(sn.DEG.list[[tmp]])[2],sep = '_')]] <- Myseufindmarker(seu.vitro.nor.list[[seu]],gene.include = gene.input.list[[seu]],
                                                                                                                     ident.1 = sn.DEG.list[[tmp]][[1]],
                                                                                                                     ident.2 = sn.DEG.list[[tmp]][[2]],c1 = names(sn.DEG.list[[tmp]])[1],c2 = names(sn.DEG.list[[tmp]])[2]
    )
  }
}

KO.DEG.nor.list.filter <- list()
for(tmp.cluster in names(KO.DEG.nor.list)[3]){
  KO.DEG.nor.list.filter[[tmp.cluster]] <- KO.DEG.nor.list[[tmp.cluster]][KO.DEG.nor.list[[tmp.cluster]]$p_val<=0.0001&
                                                                            abs( KO.DEG.nor.list[[tmp.cluster]]$avg_logFC) >= log(1.5) & 
                                                                            (KO.DEG.nor.list[[tmp.cluster]]$pct.1>=0.3 |  KO.DEG.nor.list[[tmp.cluster]]$pct.2>=0.3),]
}
table(KO.DEG.nor.list.filter$DMSO_Nif$cluster)
table(KO.DEG.nor.list.filter$DMSO_A485$cluster)
table(KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_A485$cluster)
table(KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_20um_Nif$cluster)
table(KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_A485_Nif$cluster)


saveRDS(KO.DEG.nor.list,'direct/20221120/KO.DEG.nor.list.rds')
rm(seu.vitro.nor.list)
save.image('direct/20221120/direct.20221120.RData')
load('direct/20221120/direct.20221120.RData')

KO.DEG.nor.list.filter[[paste('sm3.A485','Ctrl',sep = '_')]] <- KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_A485
KO.DEG.nor.list.filter[[paste('Nif','Ctrl',sep = '_')]] <- KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_20um_Nif
KO.DEG.nor.list.filter[[paste("sm3.Nif.A485.only",'Ctrl',sep = '_')]] <- KO.DEG.nor.list.filter$Ctrl_DMSO_Ctrl_A485_Nif

for(tmp.cluster in c("sm3.Nif.A485.only_Ctrl","Nif_Ctrl","Ctrl_DMSO_Ctrl_A485")){
  MyWriteTable(cbind(genes.inf.input[KO.DEG.nor.list.filter[[tmp.cluster]]$gene,],KO.DEG.nor.list.filter[[tmp.cluster]]),paste('direct/20221120/path/gene_inf/20230704.fc1.5.pval0.0001.',seu,'.',tmp.cluster,'.si.tab',sep = ''))
}


saveRDS(seu.vitro.nor.list,'seu.vitro.nor.list.rds')
seu.vitro.nor.list$Nif$Type <- factor(seu.vitro.nor.list$Nif$Type,levels = c('Preg_DMSO','Ctrl_DMSO','Preg_20um_Nif','Ctrl_20um_Nif'))

seu.vitro.nor.list[['CtrlPregnifA485']] <- subset(seu.vitro.nor.list$sm3.NifA485.all,
                                                  cells = colnames(seu.vitro.nor.list$sm3.NifA485.all)[seu.vitro.nor.list$sm3.NifA485.all$Type %in% names(table(seu.vitro.nor.list$sm3.Nif.A485.only$Type))]
                                                  )
seu.vitro.nor.list[['CtrlPregnifA485']]$Type <- factor(seu.vitro.nor.list[['CtrlPregnifA485']]$Type,levels = names(table(seu.vitro.nor.list$sm3.Nif.A485.only$Type)))
seu <- 'CtrlPregnifA485'
seu <- 'sm3.Nif.A485.only'
for(seu in names(seu.vitro.nor.list)[c(1,6)]){
  c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[order(seu.vitro.nor.list[[seu]]$Type) #colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$KOScore)]
                                                            ,"Type"],
                        c('skyblue1','sienna1','deeppink1','purple1'))
  c1Color <- as.matrix(c1Color)
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'doCol.CtrlDEGfc1.5.pval0.0001.heatmap.pdf',sep = ''),10,13)
  preg.row.tree[[seu]] <-
    MyHeatmap(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[KO.DEG.nor.list.filter[['Ctrl_DMSO_Ctrl_A485']]$gene,
                                                                        colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$Type)]
                                                                        ],
              type = "log.row.relat",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = c1Color,
              ColSideColorsSize = 1.5,
              return.tree = "row",
              #Colv = 'none',
              dendrogram = 'both',
              graph = T)
  rColor <- MyName2Col(cutree(preg.row.tree[[seu]],
                              2),
                       gene.tree.colors,
                       is.row = T)
  preg.row.tree[[seu]] <-
    MyHeatmap(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[KO.DEG.nor.list.filter[['Ctrl_DMSO_Ctrl_A485']]$gene,
                                                                       colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$Type)]
                                                                        ],
              type = "log.row.zscore",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",#color.palette = pc12.heatmap.col,
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = c1Color,
              ColSideColorsSize = 1.5,
              RowSideColors = rColor,
              RowSideColorsSize = 1.5,
              return.tree = "row",
              #Colv = 'none',
              dendrogram = 'both',
              graph = T)
  dev.off()
}

boxplot(-log10(KO.DEG.nor.list$Ctrl_DMSO_Ctrl_20um_Nif[c(labels(as.dendrogram(preg.row.tree[[seu]])[[2]][[2]])),'p_val']),
        -log10(KO.DEG.nor.list$Ctrl_DMSO_Ctrl_20um_Nif[c(labels(as.dendrogram(preg.row.tree[[seu]])[[2]][[1]])),'p_val'])
        )
##########
preg.gene.list <- list()
for(seu in names(seu.vitro.nor.list)[c(1,6)]
    ){
 
  preg.gene.list[[seu]][['Ctrl.DMSO']] <- c(labels(as.dendrogram(preg.row.tree[[seu]])[[2]]))
  preg.gene.list[[seu]][['Ctrl.treat']] <- labels(as.dendrogram(preg.row.tree[[seu]])[[1]])
}


sapply(preg.gene.list$sm3.A485, length)


dir.create('direct/20221120/path/go_kegg')
dir.create('direct/20221120/path/gene_inf')

LM.gene.inf <- list()

for(seu in names(seu.vitro.nor.list)[c(1,6)]){
  for(tmp.cluster in names(preg.gene.list[[seu]])){
    MyGOwritetable(genes.inf.input[preg.gene.list[[seu]][[tmp.cluster]],1],pvalue = 1,paste('direct/20221120/path/go_kegg/go.CtrlDEG.fc1.5.pval0.0001',seu,tmp.cluster,'tab',sep = '.'))
  }
}

MyGOwritetable(genes.inf.input[labels(as.dendrogram(preg.row.tree[[seu]])[[2]][[2]][[2]]),1],pvalue = 1,paste('direct/20221120/path/go_kegg/go.pval8.LM.Koscore',seu,'.2222.tab',sep = '.'))

#########
for(seu in names(seu.vitro.nor.list)[c(1,3)]){
  seu.vitro.nor.list[[seu]] <- ScaleData(seu.vitro.nor.list[[seu]],features = unlist(preg.gene.list[[seu]]))
  seu.vitro.nor.list[[seu]] <- RunPCA(seu.vitro.nor.list[[seu]],features = unlist(preg.gene.list[[seu]]))
  seu.vitro.nor.list[[seu]] <- SetIdent(seu.vitro.nor.list[[seu]],value = seu.vitro.nor.list[[seu]]$Type)
  print(DimPlot(seu.vitro.nor.list[[seu]],
                reduction = "pca",
                cols = c(time.colors,brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 2,
                label = F))

  seu.vitro.nor.list[[seu]]$heatmap.order <- Embeddings(seu.vitro.nor.list[[seu]],'pca')[,1]
  # if(seu=='Nif'){seu.vitro.nor.list[[seu]]$heatmap.order <- -Embeddings(seu.vitro.nor.list[[seu]],'pca')[,1]  }
preg.1.order <- MyordergenewithPseudotime(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[preg.gene.list[[seu]][['Ctrl.DMSO']],colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$heatmap.order,decreasing = T)]],
                                          graph = T,
                                          preg.gene.list[[seu]][['Ctrl.DMSO']])

preg.2.order <- MyordergenewithPseudotime(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[preg.gene.list[[seu]][['Ctrl.treat']],colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$heatmap.order,decreasing = T)]],
                                          graph = T,
                                          preg.gene.list[[seu]][['Ctrl.treat']])
c1Color <- MyName2Col(seu.vitro.nor.list[[seu]]@meta.data[colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$heatmap.order)],"Type"],
                      time.colors)
c1Color <- as.matrix(c1Color)

pdf(paste('direct/20221120/path/KOscore.',seu,".CtrlDEG.pval0.0001fc1.3.pdf",sep =  ''),
    10,13)
r1Color <- c(#rep(gene.tree.colors[1],length(beta.cc.order)),
  rep(gene.tree.colors[2],length(preg.1.order)),
  rep(gene.tree.colors[3],length(preg.2.order))
)
r1Color <- as.matrix(t(r1Color))
MyText(paste('c1:',length(preg.1.order),'\n',
             'c2:',length(preg.2.order),'\n',
             sep = ''))
MyHeatmap(as.matrix(as.matrix(exp(seu.vitro.nor.list[[seu]]@assays$RNA@data))[c(#beta.cc.order,
  preg.1.order,
  preg.2.order),colnames(seu.vitro.nor.list[[seu]])[order(seu.vitro.nor.list[[seu]]$heatmap.order)]]),
  type = "log.row.zscore",
  #hc.c.data.type = "log.row.relat",
  hc.r.data.type = "log.row.relat",
  #  color.palette = pc12.heatmap.col,
  #c.cov.method = "s",
  # r.cov.method = "s",
  Colv = "none",
  Rowv = "none",
  dendrogram = "none",
  #c.hc.method = "ward.D2",
  #r.hc.method = "ward.D2",
  RowSideColors = r1Color,
  RowSideColorsSize = 1.5,
  ColSideColors = c1Color,
  ColSideColorsSize = 1.5#,
  #return.tree = "none"
)
dev.off()
}
#######
rm(seu.vitro.nor.list,seu.ges.beta,seu.A485Nif,seu.vivo.vitro.list,seu.ko)
gc()
saveRDS(seu.vitro.nor.list,'seu.vitro.nor.list.rds')
save.image('direct/20221120/direct.20221120.RData')
save.image('direct/20221120/direct.20230606.RData')

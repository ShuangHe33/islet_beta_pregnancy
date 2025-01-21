dir.create('20221121')
setwd('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/KO/KOdirect/20220913/20221106/20221121/')
load('KOdirect.RData')
save.image('KOdirect.RData')

seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')

seu.p300 <- seu.ref.KO.list2$G0G14.5P300WTHMKO
ambigous.sym <- c(pre.ambigous.sym,beta.exo.sym,TAM.sym,TAM.sym2,Novaseq.batch.gene,exclude.stat3HMKO.gene)
save(seu.p300,ambigous.sym,file = 'P300KO.RData')
saveRDS(ambigous.sym,'ambigous.sym.rds')

seu <- 'G0G18.5Acss2WTHMKO'
seu.ref.KO.list <- list()
seu.ref.KO.list2[[seu]]$KOtype <- seu.ref.KO.list2[[seu]]$Type
seu.ges.beta$KOtype <- '/'
seu.ref.KO.list[[seu]] <- merge(seu.ges.beta,seu.ref.KO.list2[[seu]])
seu.ref.KO.list[[seu]]$KOtype <- factor(seu.ref.KO.list[[seu]]$KOtype,levels = c('/','G0_Acss2KOWT','G0_Acss2KOHMKO','G18.5_Acss2WT','G18.5_Acss2HMKO'))
seu.ref.KO.list2[['ref.Acss2']] <- seu.ref.KO.list[[seu]]
rm(seu.ref.KO.list)

seu <- 'ref.Acss2'
seu.ref.KO.list <- readRDS('../../seu.ref.KO.list.rds')
seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
seu <- 'G0G14.5G18.5Acss2WTHMKO'
seu.ref.KO.list2[[seu]] <- merge(subset(seu.ref.KO.list$G14.5G18.5Acss2HMKO,cells = colnames(seu.ref.KO.list$G14.5G18.5Acss2HMKO)[
  seu.ref.KO.list$G14.5G18.5Acss2HMKO$batch %in% c('92229','92932','92206')
]), seu.ref.KO.list2$G0G18.5Acss2WTHMKO)

seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type <- factor(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G14.5_Acss2WT','G14.5_Acss2HMKO','G18.5_Acss2WT','G18.5_Acss2HMKO'))
table(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type)


seu <- 'G0G14.5G18.5Acss2WTHMKO_2'
seu.ref.KO.list2[[seu]] <- merge(subset(seu.ref.KO.list$G14.5G18.5Acss2HMKO.add,cells = colnames(seu.ref.KO.list$G14.5G18.5Acss2HMKO.add)[
  seu.ref.KO.list$G14.5G18.5Acss2HMKO.add$Type %in% c('G14.5_Acss2WT','G14.5_Acss2HMKO')
  ]), seu.ref.KO.list2$G0G18.5Acss2WTHMKO)

seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type <- factor(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G14.5_Acss2WT','G14.5_Acss2HMKO','G18.5_Acss2WT','G18.5_Acss2HMKO'))
table(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO$Type)

seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2$Type <- factor(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G14.5_Acss2WT','G14.5_Acss2HMKO','G18.5_Acss2WT','G18.5_Acss2HMKO'))
table(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2$Type)

seu.ges.beta.sub <- subset(seu.ges.beta,cells = colnames(seu.ges.beta)[seu.ges.beta$Type %in% c('Virgin','G6.5','G10.5','G14.5','G18.5')])
seu.ges.beta.sub <- subset(seu.ges.beta.sub,cells = colnames(seu.ges.beta.sub)[!seu.ges.beta.sub$Type_rep %in% c('Virgin_rep1','Virgin_rep2')])

seu <- 'sub.ref.G0G14.5G18.5Acss2'

seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2$KOtype <- seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2$Type
seu.ref.KO.list2[[seu]] <- merge(seu.ges.beta.sub,seu.ref.KO.list2$G0G18.5Acss2WTHMKO)

seu.ref.KO.list2[[seu]]$Type <- factor(seu.ref.KO.list2[[seu]]$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G6.5','G10.5','G18.5_Acss2WT','G18.5_Acss2HMKO'))

seu <- 'sub.ref.G0G14.5G18.5Acss2.2'
seu.ref.KO.list2[[seu]] <- subset(seu.ref.KO.list2$sub.ref.G0G14.5G18.5Acss2,cells = colnames(seu.ref.KO.list2$sub.ref.G0G14.5G18.5Acss2)[
  !seu.ref.KO.list2$sub.ref.G0G14.5G18.5Acss2$Type_rep %in% c('G14.5_Acss2WT_rep2','G14.5_Acss2WT_rep3')
])

seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
seu.ges.beta.sub <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/cluster/seu.beta.ges.post.rds')
seu.ges.beta.sub <- subset(seu.ges.beta.sub,cells = colnames(seu.ges.beta.sub)[seu.ges.beta.sub$Type %in% c('G6.5')])

for(seu in names(seu.ref.KO.list2)){
  seu.ref.KO.list2[[paste('G6.5',seu,sep = '_')]] <- merge(seu.ges.beta.sub,seu.ref.KO.list2[[seu]])
}

names(seu.ref.KO.list2)
table(seu.ref.KO.list2$G0G18.5Acss2WTHMKO$Type,seu.ref.KO.list2$G0G18.5Acss2WTHMKO$batch)
table(seu.ref.KO.list2$G0G18.5Stat3WTHMKO$Type,seu.ref.KO.list2$G0G18.5Stat3WTHMKO$batch)
table(seu.ref.KO.list2$G0G14.5P300WTHMKO$Type,seu.ref.KO.list2$G0G14.5P300WTHMKO$batch)


seu.ref.KO.list2$G6.5_G0G18.5Acss2WTHMKO$Type <- factor(seu.ref.KO.list2$G6.5_G0G18.5Acss2WTHMKO$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G6.5','G18.5_Acss2WT','G18.5_Acss2HMKO'))
seu.ref.KO.list2$G6.5_G0G18.5Stat3WTHMKO$Type <- factor(seu.ref.KO.list2$G6.5_G0G18.5Stat3WTHMKO$Type,levels = c('G0_Stat3KOWT','Virgin_Stat3HMKO','G6.5','G18.5_Stat3WT','G18.5_Stat3HMKO'))
seu.ref.KO.list2$G6.5_G0G14.5P300WTHMKO$Type <- factor(seu.ref.KO.list2$G6.5_G0G14.5P300WTHMKO$Type,levels = c('Virgin_P300KOWT','Virgin_P300HMKO','G6.5','G14_5_P300WT','G14_5_P300HMKO'))

for(seu in names(seu.ref.KO.list2)){
seu.ref.KO.list2[[paste(seu,'.WT',sep = '')]] <- subset(seu.ref.KO.list2[[seu]],cells = colnames(seu.ref.KO.list2[[seu]])[
  seu.ref.KO.list2[[seu]]$Type %in% names(table(seu.ref.KO.list2[[seu]]$Type))[c(1,3)]
])
}

###########
vst.row.tree2 <- list()
for(seu in names(seu.ref.KO.list2)[4:5]){
  seu.ref.KO.list2[[seu]] <- FindVariableFeatures(seu.ref.KO.list2[[seu]], selection.method = "vst", nfeatures = 2000)
  VariableFeatures(seu.ref.KO.list2[[seu]]) <- setdiff(VariableFeatures(seu.ref.KO.list2[[seu]]),
                                                      c(pre.ambigous.sym,beta.exo.sym,TAM.sym,TAM.sym2,Novaseq.batch.gene,exclude.stat3HMKO.gene))
  
  c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Type"],
                        time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Rep"],
                        colors.group)
  c2Color <- as.matrix(c2Color)

  cColor <- cbind(c1Color,c2Color)

  png(paste(seu,'.rmcc.heatmap.raw.png',sep = ''),2000,3000)
  vst.row.tree2[[seu]] <-
    MyHeatmap(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[VariableFeatures(seu.ref.KO.list2[[seu]]),colnames(seu.ref.KO.list2[[seu]])],
              type = "log.row.relat",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = cColor,
              ColSideColorsSize = 2,
              return.tree = "row",
              graph = T)
  dev.off()
  #VariableFeatures(seu.ref.KO.list2[[seu]]) <- setdiff(VariableFeatures(seu.ref.KO.list2[[seu]]),c(labels(as.dendrogram(vst.row.tree2[[seu]])[[1]])))
}
VariableFeatures(seu.ref.KO.list2$G0G14.5P300WTHMKO.WT) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G14.5P300WTHMKO.WT),c(labels(as.dendrogram(vst.row.tree2$G0G14.5P300WTHMKO.WT)[[1]])))
VariableFeatures(seu.ref.KO.list2$G0G18.5Stat3WTHMKO.WT) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G18.5Stat3WTHMKO.WT),c(labels(as.dendrogram(vst.row.tree2$G0G18.5Stat3WTHMKO.WT)[[2]][[1]])))
VariableFeatures(seu.ref.KO.list2$G0G18.5Acss2WTHMKO.WT) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G18.5Acss2WTHMKO.WT),c(labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO.WT)[[2]][[1]]),
                                                                                                                               labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO.WT)[[1]][[2]][[2]][[1]])
                                                                                                                               ))

# VariableFeatures(seu.ref.KO.list2$G0G14.5P300WTHMKO) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G14.5P300WTHMKO),c(
#   labels(as.dendrogram(vst.row.tree2$G0G14.5P300WTHMKO)[[1]])))
# 
# VariableFeatures(seu.ref.KO.list2$G0G18.5Acss2WTHMKO) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G18.5Acss2WTHMKO),c(
#   labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]])))
# 
# VariableFeatures(seu.ref.KO.list2$G0G18.5Stat3WTHMKO) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G18.5Stat3WTHMKO),c(
#   labels(as.dendrogram(vst.row.tree2$G0G18.5Stat3WTHMKO)[[2]][[1]]),
#   labels(as.dendrogram(vst.row.tree2$G0G18.5Stat3WTHMKO)[[2]][[2]][[2]][[1]])
#   ))
# 
# VariableFeatures(seu.ref.KO.list2$ref.Acss2) <- setdiff(VariableFeatures(seu.ref.KO.list2$ref.Acss2),c(
#   labels(as.dendrogram(vst.row.tree2$ref.Acss2)[[1]])))
# 
# VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO),c(
#   labels(as.dendrogram(vst.row.tree2$G0G14.5G18.5Acss2WTHMKO)[[1]])))
# 
# VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO) <- c(VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO),
#                                                                 labels(as.dendrogram(vst.row.tree2$G0G14.5G18.5Acss2WTHMKO)[[1]])
#                                                                 )
# 
# VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2) <- setdiff(VariableFeatures(seu.ref.KO.list2$G0G14.5G18.5Acss2WTHMKO_2),c(
#   labels(as.dendrogram(vst.row.tree2$G0G14.5G18.5Acss2WTHMKO_2)[[1]])))
# 
# VariableFeatures(seu.ref.KO.list2$sub.ref.G0G14.5G18.5Acss2) <- setdiff(VariableFeatures(seu.ref.KO.list2$sub.ref.G0G14.5G18.5Acss2),c(
#   labels(as.dendrogram(vst.row.tree2$sub.ref.G0G14.5G18.5Acss2)[[1]])))
###########
cc.col.tree.list <- list()
for (seu in names(seu.ref.KO.list2)[c(1)]) {
  
  c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Type"],
                        ref.time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Rep"],
                        rep.colors[-1])
  c2Color <- as.matrix(c2Color)
  
  cColor <- cbind(c1Color,c2Color)
  
  png(paste(seu,".col.cc.heatmap.png",sep = ""),2000,2000)
  cc.col.tree.list[[seu]] <-
    MyHeatmap(as.matrix(MyGeneExp(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data)[labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]]),colnames(seu.ref.KO.list2[[seu]])],1,1)),
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
cc.cell.list <- list()
cc.cell.list[[seu]] <- c(labels(as.dendrogram(cc.col.tree.list$G0G18.5Acss2WTHMKO)[[1]]),
                         labels(as.dendrogram(cc.col.tree.list$G0G18.5Acss2WTHMKO)[[2]][[1]])
                         )

for (seu in names(seu.ref.KO.list2)[c(1)]) {
  seu.ref.KO.list2[[seu]]$proliferation <- 'quiescent'
  seu.ref.KO.list2[[seu]]@meta.data[cc.cell.list[[seu]],'proliferation'] <- 'proliferative'
}
#########
Vbeta.var.co15.list2 <- list()
for(seu in names(seu.ref.KO.list2)[4:6]){
  Vbeta.var.co15.list2[[seu]] <- MyCo(as.matrix(seu.ref.KO.list2[[seu]]@assays$RNA@data),
                                     var.gene = c(VariableFeatures(seu.ref.KO.list2[[seu]])
                                                  #,labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]])
                                                  ),
                                     exp.prop.whole.max = 1,
                                     exp.prop.whole.min = 0.01,
                                     # vector.group = samples.inf.qc$GroupNameFig1,
                                     # exp.prop.group.min = 0.1,
                                     # exp.prop.group.max = 0.5,
                                     cor.method = "rho",
                                     exp.cutoff = 1,
                                     cor.cutoff = 0.2,
                                     partner.cutoff = 10,
                                     refine.cor = T)
  
  c1Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Type"],
                        time.colors)
  c1Color <- as.matrix(c1Color)
  c2Color <- MyName2Col(seu.ref.KO.list2[[seu]]@meta.data[,"Rep"],
                        colors.group)
  c2Color <- as.matrix(c2Color)
  
  cColor <- cbind(c1Color,c2Color)
  png(paste(seu,'.rmcc.Vstvar.2000.cor0.2.1.10.0.02.raw.5.png',sep = ''),2000,3000)
  corVG1418.row.tree.list[[seu]] <-
    MyHeatmap(as.matrix(exp(seu.ref.KO.list2[[seu]]@assays$RNA@data))[rownames(Vbeta.var.co15.list2[[seu]]),colnames(seu.ref.KO.list2[[seu]])],
              type = "log.row.relat",
              hc.c.data.type = "log.row.relat",
              hc.r.data.type = "log.row.relat",
              c.cov.method = "s",
              r.cov.method = "s",
              c.hc.method = "ward.D2",
              r.hc.method = "ward.D2",
              ColSideColors = cColor,
              ColSideColorsSize = 2,
              return.tree = "row",
              graph = T)
  dev.off()
}
labels(as.dendrogram(corVG1418.row.tree.list[[seu]])[[2]][[2]][[1]])
Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT <- Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT[setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                      c(labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]]),
                                                        labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[2]][[1]][[2]][[1]])
                                                        )),setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                   c(labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]]),
                                                                     labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[2]][[1]][[2]][[1]])
                                                                   ))]

Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT <- Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT[setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                 c(
                                                                                                   labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[2]][[2]][[2]][[1]])
                                                                                                 )),setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                            c(
                                                                                                              labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[2]][[2]][[2]][[1]])
                                                                                                            ))]
Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT <- Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT[setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                 c(
                                                                                                   labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]][[2]][[1]][[1]])
                                                                                                 )),setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                            c(
                                                                                                              labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]][[2]][[1]][[1]])
                                                                                                            ))]
Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT <- Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT[setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                 c(
                                                                                                   labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]][[1]])
                                                                                                 )),setdiff(rownames(Vbeta.var.co15.list2$G0G18.5Stat3WTHMKO.WT),
                                                                                                            c(
                                                                                                              labels(as.dendrogram(corVG1418.row.tree.list$G0G18.5Stat3WTHMKO.WT)[[1]][[1]])
                                                                                                            ))]
############
seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
DEG.tab <- MyReadDelim('../../../../../beta/20221203/Glut2H/cluster/marker/Glut2H.DEG.tab')

for(seu in names(seu.ref.KO.list2)[4:6]){
  seu.ref.KO.list2[[seu]] <- ScaleData(seu.ref.KO.list2[[seu]],features = rownames(Vbeta.var.co15.list2[[seu]]))
 # seu.ref.KO.list2[[seu]] <- RunPCA(seu.ref.KO.list2[[seu]],features = c(rownames(Vbeta.var.co15.list2[[seu]]),labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]])))
 # seu.ref.KO.list2[[seu]] <- ScaleData(seu.ref.KO.list2[[seu]],features = rownames(Vbeta.var.co15.list2[[seu]]))
 seu.ref.KO.list2[[seu]] <- RunPCA(seu.ref.KO.list2[[seu]],features = c(rownames(Vbeta.var.co15.list2[[seu]])))
 # 
 #  
  # seu.ref.KO.list2[[seu]] <- ScaleData(seu.ref.KO.list2[[seu]],features = VariableFeatures(seu.ref.KO.list2[[seu]]))
   #seu.ref.KO.list2[[seu]] <- RunPCA(seu.ref.KO.list2[[seu]],features = VariableFeatures(seu.ref.KO.list2[[seu]]))
  # 
  
  # seu.ref.KO.list2[[seu]] <- ScaleData(seu.ref.KO.list2[[seu]],features = c(VariableFeatures(seu.ref.KO.list2[[seu]]),labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]])))
  # seu.ref.KO.list2[[seu]] <- RunPCA(seu.ref.KO.list2[[seu]],features = c(VariableFeatures(seu.ref.KO.list2[[seu]]),labels(as.dendrogram(vst.row.tree2$G0G18.5Acss2WTHMKO)[[2]][[1]])))
  # 
  pc.use <- 1:6;#cor0.15
  seu.ref.KO.list2[[seu]] <- FindNeighbors(seu.ref.KO.list2[[seu]], dims = pc.use)
  seu.ref.KO.list2[[seu]] <- FindClusters(seu.ref.KO.list2[[seu]], resolution = 1.2)
  seu.ref.KO.list2[[seu]] <- RunUMAP(seu.ref.KO.list2[[seu]],dims = pc.use)
  # 
  # # seu.ref.KO.list2[[seu]] <- RunTSNE(seu.ref.KO.list2[[seu]],
  # #                                   dims = pc.use,
  # #                                   perplexity= round((30+ncol(seu.ref.KO.list2[[seu]])/100)),
  # #                                   check_duplicates = F)
  # 
  pdf(paste(seu,'pc6.res1.2.submarker.umap.WT.pdf',sep = ''),
       6,7)
  seu.ref.KO.list2[[seu]] <-
    SetIdent(seu.ref.KO.list2[[seu]], value = seu.ref.KO.list2[[seu]]$RNA_snn_res.1.2)
  print(DimPlot(seu.ref.KO.list2[[seu]],
                reduction = "umap",
                cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 2,
                label = T))
  # 
  seu.ref.KO.list2[[seu]] <-
    SetIdent(seu.ref.KO.list2[[seu]], value = seu.ref.KO.list2[[seu]]$Type)
  
  print(DimPlot(seu.ref.KO.list2[[seu]],
                reduction = "umap",
                cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 2,
                label = F))
  Myseuratmarker(seu.ref.KO.list2[[seu]],c('Un3','Slc2a2','Prlr',DEG.tab$gene[DEG.tab$cluster=='Glut2H2_2']),reduction = 'umap')
  
  # print(DimPlot(seu.ref.KO.list2[[seu]],
  #               reduction = "pca",
  #               cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
  #               label.size = 6,
  #               sizes.highlight = 4,
  #               pt.size = 2,
  #               label = F))
  # 
  seu.ref.KO.list2[[seu]] <-
    SetIdent(seu.ref.KO.list2[[seu]], value = seu.ref.KO.list2[[seu]]$Type_rep)
  print(DimPlot(seu.ref.KO.list2[[seu]],
                reduction = "umap",
                cols = c(time.colors,'gray30','black',brewer.pal(8,"Set2"),brewer.pal(8,"Set1")),
                label.size = 6,
                sizes.highlight = 4,
                pt.size = 2,
                label = F))
 dev.off()
}

seu.ref.KO.list2$G0G14.5P300WTHMKO.WT <- SetIdent(seu.ref.KO.list2$G0G14.5P300WTHMKO.WT,value = seu.ref.KO.list2$G0G14.5P300WTHMKO.WT$RNA_snn_res.0.5)
subG14.5WT <- Myseufindmarker(seu.ref.KO.list2$G0G14.5P300WTHMKO.WT,gene.include = gene.input3,ident.1 = '1',ident.2 = '2',c1 = '1',c2 = '2')
saveRDS(subG14.5WT,'subG14.5WT.DEG.rds')
subG14.5WT.filter <- subG14.5WT[abs(subG14.5WT$avg_logFC)>=log(1.5)&subG14.5WT$p_val<=0.001,]
table(subG14.5WT.filter$cluster)

subG14.5WT.filter <- subG14.5WT.filter[order(subG14.5WT.filter$avg_logFC,decreasing = F),]
subG14.5WT['Lrrc55',]

pdf('P300.res0.2.2.marker.2.pdf',6,7)
Myseuratmarker(seu.ref.KO.list2$G0G14.5P300WTHMKO.WT,subG14.5WT.filter$gene,reduction = 'umap')
dev.off()
Myseuratmarker(seu.ref.KO.list2$G0G18.5Stat3WTHMKO,c('Tph1','Cldn8','Oxtr'),reduction = 'umap')
Myseuratmarker(seu.ref.KO.list2$G0G18.5Acss2WTHMKO,c('Tph1','Cldn8'),reduction = 'umap')


seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
seu.ref.KO.list2[[seu]]@meta.data[seu.ref.KO.list2[[seu]]$Treatment=='G0_Stat3KOWT','Rep'] <- 'rep1'
set.seed(100)

selected.G0P300KO.sn <- colnames(seu.ref.KO.list2[[seu]])[seu.ref.KO.list2[[seu]]$Type=='Virgin_P300HMKO'][sample(1:83,40)]
seu.ref.KO.list2[[seu]]@meta.data[selected.G0P300KO.sn,'Rep'] <- 'rep2'

seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')
load('KOdirect.RData')

for(seu in names(seu.ref.KO.list2)[2:3]){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/KO/',seu,'.cg.vst.pc12.circle.0.7.3.pdf',sep = ''),15,12)
  seu.ref.KO.list2[[seu]]$Type2 <- seu.ref.KO.list2[[seu]]$Type
  seu.ref.KO.list2[[seu]]$Type2 <- factor(seu.ref.KO.list2[[seu]]$Type2,levels = rownames(ko.coord.list[[seu]]))
  p.pca <- MySeuratDR2Gg2(seu.ref.KO.list2[[seu]],seu.ref.KO.list2[[seu]]@meta.data,reduction.use = 'pca',
                          reduction.key = 'PC',estimate.variation.explain.percentage = T,
                          x.dim = 1,y.dim = 2)
  
  p.pca$data_ <- p.pca$data
  p.pca$data_$Type_ <- p.pca$data_$Type
  p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(seu.ref.KO.list2[[seu]]$Type)))
  p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]
  

  # p.pca$data_ <- p.pca$data
  # p.pca$data_$Type_ <- p.pca$data_$heter.group.merge
  # p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = c('beta3','beta2','beta1'))
  # p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]
  # 
  # 
  # plot(p.pca+
  #        scale_color_manual(values = c(time.colors[c(1,2,15,12)])) +
  #        scale_shape_manual(values = c(17,19))+
  #        geom_point(aes(x = -x.pos,
  #                       y = -y.pos,
  #                       col = Type,
  #                       shape = proliferation
  #        ),
  #        size = 3
  #        )+theme(aspect.ratio = 1))
  
  # plot(p.pca+
  #        scale_color_manual(values = c(group.col,ko.time.col[names(table(seu.ref.KO.list2[[seu]]$Type))[4:5]])) +
  #        geom_point(aes(x = -x.pos,
  #                       y = y.pos,
  #                       col = heter.group#,
  #                       # shape = State
  #        ),
  #        size = 4
  #        )+theme(aspect.ratio = 1))
  
  # plot(p.pca+
  #        scale_color_manual(values = c('gray90',time.colors)) +
  #        geom_point(aes(x = -x.pos,
  #                       y = y.pos,
  #                       col = KOtype#,
  #                       # shape = State
  #        ),
  #        size = 4
  #        )+theme(aspect.ratio = 1))
  

  plot(p.pca+
         scale_color_manual(values = colors.list[[seu]][c(1,4,3,2)]) +
         scale_shape_manual(values = c(19,17,17))+
         geom_point(aes(x = x.pos,
                        y = y.pos,
                        col = Type,
                        shape = Rep
         ),
         size = 4
         
         )+theme(aspect.ratio = 1)+stat_ellipse(aes(x = x.pos,y = y.pos,col=Type2,size=I(2.5)),level=0.7,type='norm')+
         geom_point(aes(x=ko.coord.list[[seu]][1,1],y=ko.coord.list[[seu]][1,2],size=6),colour='black')+
         geom_point(aes(x=ko.coord.list[[seu]][2,1],y=ko.coord.list[[seu]][2,2],size=6),colour='black')#+
         #geom_point(aes(x=ko.coord.list[[seu]][3,1],y=ko.coord.list[[seu]][3,2],size=6),colour='black')
  )
  dev.off()
}
saveRDS(seu.ref.KO.list2,'seu.ref.KO.list2.rds')
rm(seu.ref.KO.list2)
save.image('KOdirect.RData')

########
# library(rgl)
# library(plotly)
# library(igraph)
# plot3d(x=seu.ref.KO.list2[[seu]]@reductions$pca@cell.embeddings[,1],
#        y=seu.ref.KO.list2[[seu]]@reductions$pca@cell.embeddings[,2],
#        z=seu.ref.KO.list2[[seu]]@reductions$pca@cell.embeddings[,3],
#        xlab='PC1',ylab='PC2',zlab='PC3',
#        size = 5,
#        col=MyName2Col(seu.ref.KO.list2[[seu]]$Type_rep,
#                       c(time.colors)))
# set.seed(5)
# graph.list <- list()
# graph.list[[seu]] = graph.adjacency(as.matrix(seu.ref.KO.list2[[seu]]@graphs$RNA_snn),mode = "undirected",weighted = T)
# layout.3d = layout_with_fr(graph.list[[seu]],dim = 3)
# 
# pan.data.par3d <- par3d()
# open3d(zoom = pan.data.par3d$zoom, userMatrix = pan.data.par3d$userMatrix, windowRect = pan.data.par3d$windowRect)
# 
# size.point <- 5
# plot3d(layout.3d,
#        # aspect = 'iso',
#        xlab='FDL1',ylab='FDL2',zlab='FDL3',
#        size = size.point,
#        col = time.colors[seu.ref.KO.list2[[seu]]@meta.data$Type])
######pseudotime#######
for(seu in names(seu.ref.KO.list2)[1:3]){
  print(seu)
  si.select <- seu.ref.KO.list2[[seu]]@meta.data
  pregScore=c()
  for (time in names(table(seu.ref.KO.list2[[seu]]$Type))) {
    regress=lm(-Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,2]~
                 Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,1])$coefficients[2]
    temp=-c(1,regress)%*%rbind(Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,1],
                               Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,2])
    # temp=temp-mean(temp)
    names(temp)=si.select[si.select$Type==time,]$SampleName
    pregScore=c(pregScore,temp)
  }
  seu.ref.KO.list2[[seu]]@meta.data[names(pregScore),'pregScore'] <- pregScore
  seu.ref.KO.list2[[seu]]$pseudotime <- seu.ref.KO.list2[[seu]]$pregScore
  if(grepl('Acss2',seu)){
    seu.ref.KO.list2[[seu]]@meta.data[,'pregScore'] <- -Embeddings(seu.ref.KO.list2[[seu]],'pca')[,1]
    seu.ref.KO.list2[[seu]]$pseudotime <- seu.ref.KO.list2[[seu]]$pregScore
  }
}

for(seu in names(seu.ref.KO.list2)){
  print(seu)
  si.select <- seu.ref.KO.list2[[seu]]@meta.data
  pregScore=c()
  for (time in names(table(seu.ref.KO.list2[[seu]]$Type))) {
    regress=lm(Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,1]~
                 Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,2])$coefficients[2]
    temp=-c(1,regress)%*%rbind(Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,2],
                               Embeddings(seu.ref.KO.list2[[seu]],'pca')[si.select$Type==time,1])
    # temp=temp-mean(temp)
    names(temp)=si.select[si.select$Type==time,]$SampleName
    pregScore=c(pregScore,temp)
  }
  seu.ref.KO.list2[[seu]]@meta.data[names(pregScore),'KOScore'] <- pregScore
  seu.ref.KO.list2[[seu]]$pseudotime <- seu.ref.KO.list2[[seu]]$KOScore
  # if(seu%in%c('G0G18.5Acss2WTHMKO')){
  #   seu.ref.KO.list2[[seu]]@meta.data[,'pregScore'] <- -Embeddings(seu.ref.KO.list2[[seu]],'pca')[,1]
  #   seu.ref.KO.list2[[seu]]$pseudotime <- seu.ref.KO.list2[[seu]]$pregScore
  # }
}

seu.ref.KO.list2[[seu]]$Type2 <- seu.ref.KO.list2[[seu]]


seu.ref.KO.list2[[seu]]$pseudotime <- -seu.ref.KO.list2[[seu]]$KOScore
seu.ref.KO.list2[[seu]]$pseudotime <- -Embeddings(seu.ref.KO.list2[[seu]],'pca')[,2]
for(seu in names(seu.ref.KO.list2)[1:3]){
  #seu.ref.KO.list2[[seu]]$pseudotime <- -seu.ref.KO.list2[[seu]]$pseudotime
  print(MyPseudotimebox(seu.ref.KO.list2[[seu]]@meta.data,time.colors = time.colors,size.point = 2)+theme(aspect.ratio = 1.5))
}

print(MyPseudotimebox(seu.ref.KO.list2[[seu]]@meta.data,time.colors = time.colors,size.point = 2)+theme(aspect.ratio = 1.5))
print(MyPseudotimebox(seu.ref.KO.list2[[seu]]@meta.data,time.colors = time.colors,size.point = 2,rep = T)+theme(aspect.ratio = 1.5))


print(MyPseudotimebox(seu.ref.KO.list2[[seu]]@meta.data,time.colors = time.colors,size.point = 2,rep = T)+theme(aspect.ratio = 1.5))

seu.ref.KO.list2[[seu]]@meta.data$Type_rep <- as.factor(seu.ref.KO.list2[[seu]]@meta.data$Type_rep)
tmp.selected.si.tab <- seu.ref.KO.list2[[seu]]@meta.data[!seu.ref.KO.list2[[seu]]@meta.data$Type_rep %in% c('G14.5_Acss2WT_rep2','G14.5_Acss2WT_rep3'),]
print(MyPseudotimebox(tmp.selected.si.tab,time.colors = time.colors[c(1:4,15,12)],size.point = 2)+theme(aspect.ratio = 1.5))



seu.ref.KO.list2$G0G18.5Stat3WTHMKO$KOScore <- -seu.ref.KO.list2$G0G18.5Stat3WTHMKO$KOScore
seu.ref.KO.list2$G0G14.5P300WTHMKO$KOScore <- -seu.ref.KO.list2$G0G14.5P300WTHMKO$KOScore
seu.ref.KO.list2$G0G14.5P300WTHMKO$pregScore <- -seu.ref.KO.list2$G0G14.5P300WTHMKO$pregScore
seu.ref.KO.list2$G0G18.5Stat3WTHMKO$pregScore <- -seu.ref.KO.list2$G0G18.5Stat3WTHMKO$pregScore


seu.ref.KO.list2$G0G18.5Stat3WTHMKO$pseudotime <-  -seu.ref.KO.list2$G0G18.5Stat3WTHMKO$pseudotime
seu.ref.KO.list2$G0G14.5P300WTHMKO$pseudotime <-  -seu.ref.KO.list2$G0G14.5P300WTHMKO$pseudotime

seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')

colors.list <- list()
colors.list[["G0G18.5Acss2WTHMKO"]] <- c('limegreen','darkorange1','royalblue1','magenta1')
colors.list[["G0G18.5Stat3WTHMKO"]] <- c('aquamarine4','tan1','cornflowerblue','violetred')
colors.list[["G0G14.5P300WTHMKO"]] <- c('lightcoral','seagreen3','purple3','bisque3')

select.si.list2 <- list()
for(seu in names(seu.ref.KO.list2)[1:3]){
  # seu.ref.KO.list2[[seu]]$pseudotime <- -seu.ref.KO.list2[[seu]]$pseudotime
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/KO/',seu,'.cc.vst.pseudobox.pdf',sep = ''),10,6)
  print(MyPseudotimebox(seu.ref.KO.list2[[seu]]@meta.data,time.colors = colors.list[[seu]],size.point = 2)+theme(aspect.ratio = 1.5))
  dev.off()
  select.si.list2[[seu]] <- seu.ref.KO.list2[[seu]]@meta.data
}
#saveRDS(seu.ref.KO.list2,'seu.ref.KO.list2.rds')
#########pval#######
projct.pval.list2 <- list()
for(seu in names(select.si.list2)){

  age.list <- names(table(select.si.list2[[seu]]$Type))
  projct.pval.list2[[seu]] <- matrix("-",
                                    ncol = length(age.list),
                                    nrow = length(age.list)
  )
  colnames(projct.pval.list2[[seu]]) <- age.list
  rownames(projct.pval.list2[[seu]]) <- age.list
}

for(seu in names(select.si.list2)){
  age.list <- names(table(select.si.list2[[seu]]$Type))
  rm <- age.list [1]
  for(a1 in age.list){
    rm <- unique(c(rm,a1))
    age.list2 <- age.list[!age.list %in% rm]
    for(a2 in age.list2){
      tmp.p <- wilcox.test(select.si.list2[[seu]][select.si.list2[[seu]]$Type==a1,'pseudotime'],
                           select.si.list2[[seu]][select.si.list2[[seu]]$Type==a2,'pseudotime'],
                           paired = F
      )
      # tmp.p <- var.test(select.si.list2[[seu]][select.si.list2[[seu]]$Type==a1,'pseudotime'],
      #                      select.si.list2[[seu]][select.si.list2[[seu]]$Type==a2,'pseudotime']#,
      #                      #paired = F
      # )
      projct.pval.list2[[seu]][a1,a2] <- tmp.p$p.value
    }
  }
}
##########
for(seu in names(select.si.list2)){
  MyWriteTable(projct.pval.list2[[seu]],row.names = T,paste(seu,'.cc.vst.pseudotime.pval.tab',sep = ''))
  #MyWriteTable(select.si.list2[[seu]],row.names = T,paste(seu,'si.tab',sep = ''))
  #MyWriteTable(table(select.si.list2[[seu]]$Type),row.names = T,paste(seu,'type.count.tab',sep = ''))
}
######3dPCA######
library(plotly)
library(rgl)
plot3d(x = Embeddings(seu.ref.KO.list2[[seu]],'pca')[,1],
       y = Embeddings(seu.ref.KO.list2[[seu]],'pca')[,2],
       z = Embeddings(seu.ref.KO.list2[[seu]],'pca')[,3],
       xlab='PC1',ylab='PC2',zlab='PC3',
       size = 8,
       col=MyName2Col(seu.ref.KO.list2[[seu]]$Type,
                      time.colors))
########
select.si.list[[seu]] <- seu.ref.KO.list2[[seu]]@meta.data[seu.ref.KO.list2[[seu]]$KOtype!='/',]
select.si.list[[seu]]$Type <- factor(select.si.list[[seu]]$Type,levels = c('G0_Acss2KOWT','G0_Acss2KOHMKO','G18.5_Acss2WT','G18.5_Acss2HMKO'))
print(MyPseudotimebox(select.si.list[[seu]],time.colors = time.colors[5:8],size.point = 2)+theme(aspect.ratio = 1.5))
projct.pval.list <- list()
for(seu in names(select.si.list)[5:6]){
  age.list <- names(table(select.si.list[[seu]]$Type))
  projct.pval.list[[seu]] <- matrix("-",
                                    ncol = length(age.list),
                                    nrow = length(age.list)
  )
  colnames(projct.pval.list[[seu]]) <- age.list
  rownames(projct.pval.list[[seu]]) <- age.list
}

for(seu in names(select.si.list)){
  age.list <- names(table(select.si.list[[seu]]$Type))
  rm <- age.list [1]
  for(a1 in age.list){
    rm <- unique(c(rm,a1))
    age.list2 <- age.list[!age.list %in% rm]
    for(a2 in age.list2){
      tmp.p <- wilcox.test(select.si.list[[seu]][select.si.list[[seu]]$Type==a1,'pseudotime'],
                           select.si.list[[seu]][select.si.list[[seu]]$Type==a2,'pseudotime'],
                           paired = F
      )
      projct.pval.list[[seu]][a1,a2] <- tmp.p$p.value
    }
  }
}

for(seu in names(select.si.list)){
  MyWriteTable(projct.pval.list[[seu]],row.names = T,paste(seu,'ref.DEGG0G18.5.varrmcc.pval.tab',sep = ''))
  #MyWriteTable(select.si.list[[seu]],row.names = T,paste(seu,'.si.tab',sep = ''))
  # MyWriteTable(table(select.si.list[[seu]]$Type),row.names = T,paste(seu,'.type.count.tab',sep = ''))
}
###########
saveRDS(seu.ref.KO.list2,'seu.ref.KO.list2.rds')
rm(seu.ref.KO.list2)
save.image('KOdirect.RData')
seu.ref.KO.list2 <- readRDS('seu.ref.KO.list2.rds')

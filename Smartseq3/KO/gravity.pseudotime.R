

ko.time.col <- c(time.colors[c(5:8)],time.colors[c(5:8)],time.colors[c(1,2,15,12)])
names(ko.time.col) <- c(names(table(seu.ref.KO.list2$G0G14.5P300WTHMKO$Type)),names(table(seu.ref.KO.list2$G0G18.5Stat3WTHMKO$Type)),
                        names(table(seu.ref.KO.list2$G0G18.5Acss2WTHMKO$Type))
                        )
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

ko.coord.list <- list()
pca.df.list <- list()
seu <- 'G0G14.5P300WTHMKO'
for(seu in names(seu.ref.KO.list2)[1:3]){
  pca.df.list[[seu]] <- cbind(Embeddings(seu.ref.KO.list2[[seu]],'pca'),seu.ref.KO.list2[[seu]]@meta.data[,c('SampleName','Type')])
  #pca.df.list[[seu]] <- cbind(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,],seu.ref.KO.list2[[seu]]@meta.data[,c('SampleName','Type')])
  
  pca.df.list[[seu]] <- pca.df.list[[seu]][pca.df.list[[seu]]$Type %in% names(table(seu.ref.KO.list2[[seu]]$Type))[c(1,2,3)],]

  ko.coord.list[[seu]] <- Mycoord(pca.df.list[[seu]],dims=1:4,
                                  class = 'Type',
                                  names(table(seu.ref.KO.list2[[seu]]$Type))[c(1,2,3)])
}
for(seu in names(seu.ref.KO.list2)[1:3]){
  x <- ko.coord.list[[seu]][2,]-ko.coord.list[[seu]][1,]
  x_norm <- x/sqrt(sum(x^2))
  
  seu.ref.KO.list2[[seu]]$pseudotime <- as.matrix(Embeddings(seu.ref.KO.list2[[seu]],'pca')[,1:4]) %*% as.matrix(x_norm)
  #seu.ref.KO.list2[[seu]]$pseudotime <- as.matrix(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,1:2]) %*% as.matrix(x_norm)
}

seu <- 'G0G14.5P300WTHMKO'
for(seu in names(seu.ref.KO.list2)[1:3]){
  pca.df.list[[seu]] <- cbind(Embeddings(seu.ref.KO.list2[[seu]],'pca'),seu.ref.KO.list2[[seu]]@meta.data[,c('SampleName','Type')])
  #pca.df.list[[seu]] <- cbind(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,],seu.ref.KO.list2[[seu]]@meta.data[,c('SampleName','Type')])
  
  pca.df.list[[seu]] <- pca.df.list[[seu]][pca.df.list[[seu]]$Type %in% names(table(seu.ref.KO.list2[[seu]]$Type))[c(1,2,3)],]
  
  ko.coord.list[[seu]] <- Mycoord(pca.df.list[[seu]],dims=1:4,
                                  class = 'Type',
                                  names(table(seu.ref.KO.list2[[seu]]$Type))[c(1,2,3)])
}

for(seu in names(seu.ref.KO.list2)[1:3]){
  x <- ko.coord.list[[seu]][2,]-ko.coord.list[[seu]][1,]
  x_norm <- x/sqrt(sum(x^2))
  
  seu.ref.KO.list2[[seu]]$cg.koscore <- as.matrix(Embeddings(seu.ref.KO.list2[[seu]],'pca')[,1:4]) %*% as.matrix(x_norm)
  #seu.ref.KO.list2[[seu]]$pseudotime <- as.matrix(dif.beta.ref[[seu]][seu.ref.KO.list2[[seu]]$SampleName,1:2]) %*% as.matrix(x_norm)
}


seu.ref.KO.list2 <- readRDS('20221106/20221121/seu.ref.KO.list2.rds')
saveRDS(seu.ref.KO.list2,'20221106/20221121/seu.ref.KO.list2.rds')

select.si.list2 <- list()
for(seu in names(select.si.list2)){
  # seu.ref.KO.list2[[seu]]$pseudotime <- -seu.ref.KO.list2[[seu]]$pseudotime
  select.si.list2[[seu]] <- seu.ref.KO.list2[[seu]]@meta.data
  select.si.list2[[seu]]$pseudotime <- select.si.list2[[seu]]$cg.koscore
  pdf(paste(seu,'.gravity.vst.4pcs.pseudobox.2.pdf',sep = ''),10,6)
  print(MyPseudotimebox(select.si.list2[[seu]],time.colors = time.colors[c(5,6,10,8)],size.point = 2,box.lwd = 1.5)+theme(aspect.ratio = 1.5))
  dev.off()
  
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

for(seu in names(select.si.list2)){
  MyWriteTable(projct.pval.list2[[seu]],row.names = T,paste('20221106/20221121/',seu,'.gravity.vst.4pcs.pval.tab',sep = ''))
  #MyWriteTable(select.si.list[[seu]],row.names = T,paste(seu,'si.tab',sep = ''))
  MyWriteTable(table(select.si.list[[seu]]$Type),row.names = T,paste('20221106/20221121/',seu,'type.count.tab',sep = ''))
}

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

seu <- 'beta'
ko.coord.list <- list()
pca.df.list <- list()
seu.beta <- readRDS('seu.beta.ges.post.rds')
  pca.df.list[[seu]] <- cbind(Embeddings(seu.beta,'pca'),seu.beta@meta.data[,c('SampleName','Type')])
  #pca.df.list[[seu]] <- cbind(dif.beta.ref[[seu]][seu.beta$SampleName,],seu.beta@meta.data[,c('SampleName','Type')])
  
  #pca.df.list[[seu]] <- pca.df.list[[seu]][pca.df.list[[seu]]$Type %in% names(table(seu.beta$Type))[c(1,3)],]
  
  ko.coord.list[[seu]] <- Mycoord(pca.df.list[[seu]],dims=1:2,
                                  class = 'Type',
                                  names(table(seu.beta$Type)))
  
  regress <- lm(ko.coord.list[[seu]][,2]~
                  ko.coord.list[[seu]][,1])$coefficients[2]
  saveRDS(regress,'regress.value.rds')
  
  temp=-c(1,regress)%*%rbind(Embeddings(seu.beta,'pca')[,1],
                             Embeddings(seu.beta,'pca')[,2])
  seu.beta@meta.data[colnames(temp),]$pseudotime <- temp[1,]
  # x <- ko.coord.list[[seu]][2,]-ko.coord.list[[seu]][1,]
  # x_norm <- x/sqrt(sum(x^2))
  # 
  # seu.beta$pseudotime <- as.matrix(Embeddings(seu.beta,'pca')[,1:4]) %*% as.matrix(x_norm)
  #seu.beta$pseudotime <- as.matrix(dif.beta.ref[[seu]][seu.beta$SampleName,1:2]) %*% as.matrix(x_norm)

select.si.list2 <- list()
for(seu in names(select.si.list2)[1]){
   seu.beta$pseudotime <- -seu.beta$pseudotime
  # select.si.list2[[seu]] <- seu.beta@meta.data
  pdf(paste(seu,'.gravity.vst.4pcs.pseudobox.pdf',sep = ''),10,6)
  print(MyPseudotimebox(seu.beta@meta.data,time.colors = ref.time.colors,size.point = 2,box.lwd = 1.5))
  dev.off()
}

devtools::install_github("GuangchuangYu/yyplot")
devtools::install_github('fawda123/ggord')

library(gg_ordiplot)

library(ggord)
library(ggplot2)
pca.res <- FactoMineR::PCA(t(as.matrix(seu.beta@assays$RNA@data[setdiff(all.var.co2,'6430548M08Rik'),])))
pca.res$ind$coord[,1] <- seu.beta@reductions$pca@cell.embeddings[,'PC_1']
pca.res$ind$coord[,2] <- seu.beta@reductions$pca@cell.embeddings[,'PC_2']
pdf('circle.0.6.0.7.0.8.pdf',10,10)
ggord(pca.res,seu.beta$Type,col=ref.time.colors,arrow=NULL,vec_ext =0,txt=NULL,ellipse=T,ellipse_pro=0.6,
      poly = F,alpha=0,size=I(5))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+theme(aspect.ratio = 1)+aes(size=I(2.5))

ggord(pca.res,seu.beta$Type,col=ref.time.colors,arrow=NULL,vec_ext =0,txt=NULL,ellipse=T,ellipse_pro=0.7,
      poly = F,alpha=0,size=I(5))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+theme(aspect.ratio = 1)+aes(size=I(2.5))
ggord(pca.res,seu.beta$Type,col=ref.time.colors,arrow=NULL,vec_ext =0,txt=NULL,ellipse=T,ellipse_pro=0.8,
      poly = F,alpha=0,size=I(5))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+theme(aspect.ratio = 1)+aes(size=I(2.5))
  #+ylim(c(-6,6))+xlim(c(-10,10))
  # theme(panel.border = element_rect(size = 4,
  #                                   colour = "black"))+
  # theme(axis.text = element_text(size = 20, colour = "black")) +
  # theme(axis.title.x = element_text(size = 40, colour = "black")) +
  # theme(axis.title.y = element_text(size = 40, colour = "black")) +
  # theme(legend.text = element_text(size = 45, colour = "black")) +
  # theme(legend.title = element_text(size = 45, colour = "black")) +
  # theme(panel.border = element_rect(size = 4,
  #                                   colour = "black"))
dev.off()

pdf('cg.lm.pdf',5,5)
plot(ko.coord.list[[seu]],col=ref.time.colors,pch=20,cex=2,xlim = c(-10,10),ylim = c(-6,6))
abline(lm(ko.coord.list[[seu]][,2]~
            ko.coord.list[[seu]][,1]))
dev.off()
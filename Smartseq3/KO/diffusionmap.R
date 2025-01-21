library(diffusionMap)
dif.beta.ref <- list()
for (seu in names(seu.ref.KO.list2)[1:3]) {
  dif.beta.ref[[seu]] <- Mydiffusionmap(as.matrix(seu.ref.KO.list2[[seu]]@assays$RNA@data)[VariableFeatures(seu.ref.KO.list2[[seu]]),])
  dif.beta.ref[[seu]] <- as.data.frame(dif.beta.ref[[seu]]$X)
  rownames(dif.beta.ref[[seu]]) <- colnames(seu.ref.KO.list2[[seu]])
}
df.ko.coord.list <- list()
for(seu in names(seu.ref.KO.list2)){
  df.ko.coord.list[[seu]] <- Mycoord(cbind(dif.beta.ref[[seu]],seu.ref.KO.list2[[seu]]@meta.data[rownames(dif.beta.ref[[seu]]),c('SampleName','Type')]),dims=1:2,
                                     class = 'Type',
                                     names(table(seu.ref.KO.list2[[seu]]$Type)))
}


# plot3d(x=dif.beta.ref.df[,1],
#        y=dif.beta.ref.df[,2],
#        z=dif.beta.ref.df[,3],
#        xlab='DC1',ylab='DC2',zlab='DC3',
#        size = 5,
#        col=MyName2Col(seu.ref.KO.list2[[seu]]$Type,
#                       time.colors))

p.dif <- Mydif2Gg(dif.beta.ref[[seu]],seu.ref.KO.list2[[seu]]@meta.data)
p.age <- p.dif+
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

pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/KO/',seu,'.vst.rmcc.dif12.time2.pdf',sep = ''),
    15,
    12)
print(p.age +
  scale_color_manual(values =  colors.list[[seu]]) +
    scale_shape_manual(values = c(19,17,17))+
  geom_point(aes(x =x.pos,
                 y =y.pos,
                 col = Type,
                 shape = Rep
  ),
  size = 4
  )#+theme(aspect.ratio = 1)+theme(aspect.ratio = 1)+stat_ellipse(aes(x = x.pos,y = y.pos,col=Type2,size=I(2.5)),level=0.7,type='norm')+
    #geom_point(aes(x=df.ko.coord.list[[seu]][1,1],y=df.ko.coord.list[[seu]][1,2],size=6),colour='black')+
    #geom_point(aes(x=df.ko.coord.list[[seu]][2,1],y=df.ko.coord.list[[seu]][2,2],size=6),colour='black')+
    #geom_point(aes(x=df.ko.coord.list[[seu]][3,1],y=df.ko.coord.list[[seu]][3,2],size=6),colour='black')
  )

dev.off()

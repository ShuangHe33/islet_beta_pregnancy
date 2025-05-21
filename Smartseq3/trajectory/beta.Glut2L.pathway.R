################pathway########
seu.ref.L <- readRDS('Glut2L/seu.ref.L.endo.celltype.rds')
dir.create('Glut2L/LMgene/')
all.var.coL.02.cc <- rownames(MyCo(as.matrix(seu.ref.L@assays$RNA@data),
                             var.gene = setdiff(c(VariableFeatures(seu.ref.L),pre.cc.sym),c(pre.ambigous.sym)),
                             exp.prop.whole.max = 1,
                             exp.prop.whole.min = 0.02,
                             # vector.group = samples.inf.qc$GroupNameFig1,
                             # exp.prop.group.min = 0.1,
                             # exp.prop.group.max = 0.5,
                             cor.method = "rho",
                             cor.cutoff = 0.2,
                             partner.cutoff = 8,
                             refine.cor = T))
length(all.var.coL.02.cc)#267

saveRDS(all.var.coL,'Glut2L/all.var.coL.0.2.8.1.0.05.rds')

grep('Mt2',all.var.coL)
grep('Ucn3',all.var.coL)

png('Glut2L/LMgene/cc.cor0.2.1.0.02.8.rho.png',2000,3000)
var.co0.2.row.treeLcc <-
  MyHeatmap(as.matrix(exp(seu.ref.L@assays$RNA@data))[all.var.coL.02.cc,],
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

var.co0.2.row.treeLcc.den <- as.dendrogram(var.co0.2.row.treeLcc)
labels(var.co0.2.row.treeLcc.den[[1]][[2]][[1]])
labels(var.co0.2.row.treeLcc.den[[1]][[2]][[1]][[1]])

grep('Cd9',labels(var.co0.2.row.treeLcc.den[[1]][[1]]))
grep('Cd79a',labels(var.co0.2.row.treeLcc.den[[2]][[1]]))

Glut2L.ambigous.sym <- c(labels(var.co0.2.row.treeLcc.den[[1]][[1]]),labels(var.co0.2.row.treeLcc.den[[2]][[1]]))

Glut2L.cc.pca.input.sym <- setdiff(all.var.coL.02.cc,Glut2L.ambigous.sym)

all.var.coL.02 <- all.var.coL
saveRDS(Glut2L.cc.pca.input.sym,'Glut2L/Glut2L.cc.pca.input.sym.rds')

Glut2L.cc.pca.input.sym <- readRDS('../../20220913/cluster/Glut2L/Glut2L.cc.pca.input.sym.rds')
all.var.coL.old <- readRDS('../../20220913/cluster/Glut2L/all.var.coL.02.rds')

all.var.coL <- setdiff(all.var.coL,Glut2L.ambigous.sym)

venn.co <- Venn(list('old' = all.var.coL.old,
          'new' = all.var.coL
          ))
plot(venn.co)

Glut2L.cc.pca.input.sym.new <- readRDS('Glut2L/Glut2L.cc.pca.input.sym.rds')
Glut2L.cc.pca.input.sym.new <- Glut2L.cc.pca.input.sym.new[!Glut2L.cc.pca.input.sym.new %in% Glut2L.cc.pca.input.sym]
seu.ref.L <- readRDS('Glut2L/seu.ref.L.endo.celltype.rds')

seu.ref.L <- ScaleData(seu.ref.L, features = Glut2L.cc.pca.input.sym)
seu.ref.L <- RunPCA(seu.ref.L, features = all.var.coL)#pseudotime
seu.ref.L <- RunPCA(seu.ref.L, features = Glut2L.cc.pca.input.sym)

seu.ref.L <- readRDS('seu.ref.L.endo.celltype.rds')

pdf('G:/lab/Article/heshuang/BYLW/sm3/ref/FigS3B.ref.rmcc.cc.GLUT2L.refPreggene.pc12.final.text.pdf',14.5,8)
p.pca <- MySeuratDR2Gg2(seu.ref.L,seu.ref.L@meta.data,reduction.use = 'pca',reduction.key = 'PC',estimate.variation.explain.percentage = T)
p.pca$data_ <- p.pca$data
p.pca$data_$Type_ <- p.pca$data_$Type
p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = rev(names(table(seu.ref.L$Type))))
p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

plot(p.pca+
       scale_color_manual(values = ref.time.colors2) +
       geom_point(aes(x =-x.pos,
                      y = -y.pos,
                      col = Type#,
                      # shape = State
       ),
       size = 5
       )+theme(aspect.ratio=1))

#############cc marker#######
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
  guides(colour = guide_colorbar(title = "ln(TP0.1M + 1)",
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

pdf("Glut2L/cc.marker.cellcycle.pdf",
    6,
    7)
for (plot.ens in cellcycle.marker.set){
  print(plot.ens)
  print(p.plot <- p.loop +
          geom_point(aes(x = -x.pos,
                         y = y.pos,
                         
                         color = seu.ref.L@assays$RNA@data[plot.ens,rownames(p.pca$data)]
                         #shape = Age
          ),size = 5) +
          labs(title = plot.ens)+
          theme(aspect.ratio=1))
  
}
dev.off()


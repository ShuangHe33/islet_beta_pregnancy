
library(Seurat)

seu.ges.beta <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/cluster/seu.beta.ges.post.rds')
ref.var.co.gene <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/20221203/Glut2H/cluster/all.var.co2.10.0.01.rmcc.pathway.rds')
length(ref.var.co.gene)
seu.ges.beta <- RunPCA(seu.ges.beta,features = ref.var.co.gene)
seu.ges.beta$KOtype <- '/'
seu.vitro.nor.list <- readRDS('seu.vitro.nor.list.rds')
seu.vitro.nor.list$sm3.NifA485.all@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalTp0.1mExcludeGene(as.matrix(seu.vitro.nor.list$sm3.NifA485.all@assays$RNA@counts),
                                                                                                                 genes.inf.input[rownames(seu.vitro.nor.list$sm3.NifA485.all@assays$RNA@counts),'GeneLength'],
                                                                                                                 c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list[['sm3.A485.add']] <- subset(seu.vitro.nor.list$sm3.NifA485.all,cells = colnames(seu.vitro.nor.list$sm3.NifA485.all)[
  seu.vitro.nor.list$sm3.NifA485.all$Type %in% c('Ctrl_DMSO','Ctrl_A485','Preg_DMSO','Preg_A485')])
seu.vitro.nor.list[['sm3.A485.add']] <- subset(seu.vitro.nor.list[['sm3.A485.add']],cells = colnames(seu.vitro.nor.list[['sm3.A485.add']])[!seu.vitro.nor.list[['sm3.A485.add']]$SampleName %in% seu.vitro.nor.list$Nif$SampleName ])
ncol(seu.vitro.nor.list[['sm3.A485.add']])
length(seu.vitro.nor.list[['sm3.A485.add']]$SampleName)

table(seu.vitro.nor.list[['sm3.A485.add']]$Type,seu.vitro.nor.list[['sm3.A485.add']]$Rep)
seu.vitro.nor.list[['sm3.A485.add']]@meta.data[seu.vitro.nor.list$sm3.Nif.A485$SampleName[seu.vitro.nor.list$sm3.Nif.A485$Rep=='rep1'&seu.vitro.nor.list$sm3.Nif.A485$Type %in% c('Ctrl_DMSO','Preg_DMSO','Preg_A485')],'Rep'] <- 'rep3'
seu.vitro.nor.list[['sm3.A485.add']]@meta.data[seu.vitro.nor.list$sm3.Nif.A485$SampleName[seu.vitro.nor.list$sm3.Nif.A485$Rep=='rep2'&seu.vitro.nor.list$sm3.Nif.A485$Type %in% c('Ctrl_DMSO','Preg_DMSO','Preg_A485')],'Rep'] <- 'rep4'
table(seu.vitro.nor.list[['sm3.A485.add']]$Type,seu.vitro.nor.list[['sm3.A485.add']]$Rep)


seu.vitro.nor.list[['sm3.A485.addpregA485']] <- subset(seu.vitro.nor.list[['sm3.A485.add']],cells = colnames(seu.vitro.nor.list[['sm3.A485.add']])[
  seu.vitro.nor.list[['sm3.A485.add']]$Type_batch %in% c('Ctrl_A485_20221120','Ctrl_DMSO_20221120','Preg_A485_20221120','Preg_A485_20230606','Preg_DMSO_20221120')
])


seu.vitro.nor.list$sm3.A485.addpregA485@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalTp0.1mExcludeGene(as.matrix(seu.vitro.nor.list$sm3.A485.addpregA485@assays$RNA@counts),
                                                                                                                 genes.inf.input[rownames(seu.vitro.nor.list$sm3.A485.addpregA485@assays$RNA@counts),'GeneLength'],
                                                                                                                 c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list$sm3.A485.add$Type <- factor(seu.vitro.nor.list$sm3.A485.add$Type,levels = c('Ctrl_DMSO','Ctrl_A485','Preg_DMSO','Preg_A485'))



seu.vitro.nor.list[['sm3.NifA485.all']]@meta.data[seu.vitro.nor.list$sm3.Nif.A485$SampleName[seu.vitro.nor.list$sm3.Nif.A485$Rep=='rep1'&!seu.vitro.nor.list$sm3.Nif.A485$Type%in%c('Ctrl_A485_Nif','Preg_A485_Nif')],'Rep'] <- 'rep3'
seu.vitro.nor.list[['sm3.NifA485.all']]@meta.data[seu.vitro.nor.list$sm3.Nif.A485$SampleName[seu.vitro.nor.list$sm3.Nif.A485$Rep=='rep2'&!seu.vitro.nor.list$sm3.Nif.A485$Type%in%c('Ctrl_A485_Nif','Preg_A485_Nif')],'Rep'] <- 'rep4'
seu.vitro.nor.list[['sm3.NifA485.all']]@meta.data[seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Rep=='rep1'&seu.vitro.nor.list$Nif$Type%in%c('Ctrl_DMSO','Preg_DMSO')],'Rep'] <- 'rep5'
seu.vitro.nor.list[['sm3.NifA485.all']]@meta.data[seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Rep=='rep2'&seu.vitro.nor.list$Nif$Type%in%c('Ctrl_DMSO','Preg_DMSO')],'Rep'] <- 'rep6'

table(seu.vitro.nor.list[['sm3.NifA485.all']]$Type,seu.vitro.nor.list[['sm3.NifA485.all']]$Rep)


seu.vivo.vitro.list <- list()
seu.vitro.nor.list$sm3.A485.add$KOtype <- as.character(seu.vitro.nor.list$sm3.A485.add$Type)

seu.vivo.vitro.list[['sm3.A485.add']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.A485.add)
seu.vivo.vitro.list[['sm3.A485.add']]@meta.data[seu.vivo.vitro.list[['sm3.A485.add']]$KOtype=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['sm3.A485.add']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.A485.add']]$KOtype,levels = c('/',names(table(seu.vitro.nor.list$sm3.A485.add$Type))))


seu.vitro.nor.list$sm3.A485.add$KOtype <- as.character(seu.vitro.nor.list$sm3.A485.add$Type)
seu.vivo.vitro.list[['sm3.A485.addpregA485']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.A485.addpregA485)
seu.vivo.vitro.list[['sm3.A485.addpregA485']]@meta.data[seu.vivo.vitro.list[['sm3.A485.addpregA485']]$KOtype=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['sm3.A485.addpregA485']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.A485.addpregA485']]$KOtype,levels = c('/',names(table(seu.vitro.nor.list$sm3.A485.addpregA485$Type))))


seu.vitro.nor.list[['sm3.NifA485.all.ex']] <- subset(seu.vitro.nor.list[['sm3.NifA485.all']],cells = colnames(seu.vitro.nor.list[['sm3.NifA485.all']])[
  !seu.vitro.nor.list[['sm3.NifA485.all']]$SampleName %in% seu.vitro.nor.list$sm3.A485.addpregA485$SampleName
  ])
table(seu.vitro.nor.list$sm3.NifA485.all.ex$Type,seu.vitro.nor.list$sm3.NifA485.all.ex$Rep)
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$Rep=='rep3','Rep'] <- 'rep1'
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$Rep=='rep4','Rep'] <- 'rep2'
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$Rep=='rep5','Rep'] <- 'rep3'
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$Rep=='rep6','Rep'] <- 'rep4'
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$SampleName %in% colnames(seu.vitro.nor.list$sm3.Nif.A485)[seu.vitro.nor.list$sm3.Nif.A485$Type_rep=='Preg_20um_Nif_rep1'],'Rep'] <- 'rep3'
seu.vitro.nor.list$sm3.NifA485.all.ex@meta.data[seu.vitro.nor.list$sm3.NifA485.all.ex$SampleName %in% colnames(seu.vitro.nor.list$sm3.Nif.A485)[seu.vitro.nor.list$sm3.Nif.A485$Type_rep=='Preg_20um_Nif_rep2'],'Rep'] <- 'rep4'

seu.vitro.nor.list[['sm3.NifA485.all.ex']]$Type <- factor(seu.vitro.nor.list[['sm3.NifA485.all.ex']]$Type,levels = c('Ctrl_DMSO','Ctrl_20um_Nif','Ctrl_A485_Nif','Preg_DMSO','Preg_20um_Nif','Preg_A485_Nif'))
seu.vitro.nor.list$sm3.NifA485.all.ex@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalTp0.1mExcludeGene(as.matrix(seu.vitro.nor.list$sm3.NifA485.all.ex@assays$RNA@counts),
                                                                                                                    genes.inf.input[rownames(seu.vitro.nor.list$sm3.NifA485.all.ex@assays$RNA@counts),'GeneLength'],
                                                                                                                    c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list[['sm3.NifA485.all.ex']]$KOtype <- as.character(seu.vitro.nor.list[['sm3.NifA485.all.ex']]$Type)
seu.vivo.vitro.list[['sm3.NifA485.all.ex']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.NifA485.all.ex)
seu.vivo.vitro.list[['sm3.NifA485.all.ex']]@meta.data[as.character(seu.vivo.vitro.list[['sm3.NifA485.all.ex']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['sm3.NifA485.all.ex']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.NifA485.all.ex']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Ctrl_A485_Nif','Preg_DMSO','Preg_20um_Nif','Preg_A485_Nif'))


seu.vitro.nor.list[['sm3.Nif.A485']] <- seu.NS.MS.NifA485
seu.vitro.nor.list[['sm3.Nif.A485']]$KOtype <- as.character(seu.vitro.nor.list[['sm3.Nif.A485']]$Type)

seu.vivo.vitro.list[['sm3.Nif.A485']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.Nif.A485)
seu.vivo.vitro.list[['sm3.Nif.A485']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.Nif.A485']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_A485_Nif','Preg_DMSO','Preg_A485','Preg_20um_Nif','Preg_A485_Nif'))
seu.vivo.vitro.list[['sm3.Nif.A485']]@meta.data[seu.vivo.vitro.list[['sm3.Nif.A485']]$KOtype=='/','Rep'] <- 'rep1'

seu.vitro.nor.list$sm3.A485$KOtype <- as.character(seu.vitro.nor.list$sm3.A485$Type)
seu.vivo.vitro.list <- list()
seu.vivo.vitro.list[['Nif']] <- merge(seu.ges.beta,seu.vitro.nor.list$Nif)
seu.vivo.vitro.list[['sm3.A485']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.A485)

seu.vivo.vitro.list[['Nif']]@meta.data[as.character(seu.vivo.vitro.list[['Nif']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['sm3.A485']]@meta.data[seu.vivo.vitro.list[['sm3.A485']]$KOtype=='/','Rep'] <- 'rep1'

seu.vivo.vitro.list[['Nif']]$KOtype <- factor(seu.vivo.vitro.list[['Nif']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))
seu.vivo.vitro.list[['sm3.A485']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.A485']]$KOtype,levels = c('/',names(table(seu.vitro.nor.list$sm3.A485$Type))))

seu.vitro.nor.list$Nif$KOtype <- as.character(seu.vitro.nor.list$Nif$Type)
seu.vitro.nor.list$sm3.A485$KOtype <- as.character(seu.vitro.nor.list$sm3.A485$Type)
seu.vitro.nor.list$sm3.Nif.A485$KOtype <- as.character(seu.vitro.nor.list$sm3.Nif.A485$Type)


#seu.vivo.vitro.list[['sm3.NifA485.all']] <- merge(seu.ges.beta,seu.vitro.nor.list[c('Nif','sm3.A485','sm3.Nif.A485')])
seu.vitro.nor.list[['sm3.NifA485.all']]$KOtype <- as.character(seu.vitro.nor.list[['sm3.NifA485.all']]$Type)
seu.vivo.vitro.list[['sm3.NifA485.all']] <- merge(seu.ges.beta,seu.vitro.nor.list$sm3.NifA485.all)
seu.vivo.vitro.list[['sm3.NifA485.all']]@meta.data[as.character(seu.vivo.vitro.list[['sm3.NifA485.all']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['sm3.NifA485.all']]$KOtype <- factor(seu.vivo.vitro.list[['sm3.NifA485.all']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_A485','Ctrl_20um_Nif','Ctrl_A485_Nif','Preg_DMSO','Preg_A485','Preg_20um_Nif','Preg_A485_Nif'))



seu.vitro.nor.list[['sm3.Nif.add']] <- subset(seu.vitro.nor.list$sm3.NifA485.all.ex,cells = colnames(seu.vitro.nor.list$sm3.NifA485.all.ex)[
  seu.vitro.nor.list$sm3.NifA485.all.ex$Type %in% c('Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif')])
ncol(seu.vitro.nor.list[['sm3.Nif.add']])#408
length(seu.vitro.nor.list[['sm3.Nif.add']]$SampleName)#408
table(seu.vitro.nor.list[['sm3.Nif.add']]$Type,seu.vitro.nor.list[['sm3.Nif.add']]$Rep)
seu.vitro.nor.list$sm3.Nif.add$Type <- factor(seu.vitro.nor.list$sm3.Nif.add$Type,levels = c('Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))


seu.vivo.vitro.list <- list()
seu.vivo.vitro.list[['sm3.NS.MS']] <- merge(seu.ges.beta,c(subset(seu.vitro.nor.list$Nif,Type %in% c('Ctrl_DMSO','Preg_DMSO')),
                                                           subset(seu.vitro.nor.list$sm3.A485,Type %in% c('Ctrl_DMSO','Preg_DMSO'))
                                                           ))

set.seed(10)
seu.vitro.nor.list[['Nif.sampling']] <- subset(seu.vitro.nor.list$Nif,cells = c(seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'),30)],
                                                                                seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'),30)],
                                                                                seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_DMSO'),30)],
                                                                                seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'),30)]
))
table(seu.vitro.nor.list[['Nif.sampling']]$Type)

seu.vivo.vitro.list[['Nif.sampling']] <- merge(seu.ges.beta,seu.vitro.nor.list$Nif.sampling)
seu.vivo.vitro.list[['Nif.sampling']]@meta.data[as.character(seu.vivo.vitro.list[['Nif.sampling']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['Nif.sampling']]$KOtype <- factor(seu.vivo.vitro.list[['Nif.sampling']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))


set.seed(100)
cell.1 <- seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'),40)]
set.seed(100)
cell.2 <- seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'),40)]
set.seed(100)
cell.3 <- seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_DMSO'),40)]
set.seed(100)
cell.4 <- seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'),40)]



seu.vitro.nor.list[['Nif.sampling40']] <- subset(seu.vitro.nor.list$Nif,cells =c(cell.1,cell.2,cell.3,cell.4))
table(seu.vitro.nor.list[['Nif.sampling40']]$Type)

seu.vivo.vitro.list[['Nif.sampling40']] <- merge(seu.ges.beta,seu.vitro.nor.list[['Nif.sampling40']])
seu.vivo.vitro.list[['Nif.sampling40']]@meta.data[as.character(seu.vivo.vitro.list[['Nif.sampling40']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['Nif.sampling40']]$KOtype <- factor(seu.vivo.vitro.list[['Nif.sampling40']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))

set.seed(50)
seu.vitro.nor.list[['Nif.sampling25']] <- subset(seu.vitro.nor.list$Nif,cells = c(seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'),25)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'),25)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_DMSO'),25)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'),25)]
))
table(seu.vitro.nor.list[['Nif.sampling25']]$Type)

seu.vivo.vitro.list[['Nif.sampling25']] <- merge(seu.ges.beta,seu.vitro.nor.list[['Nif.sampling25']])
seu.vivo.vitro.list[['Nif.sampling25']]@meta.data[as.character(seu.vivo.vitro.list[['Nif.sampling25']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['Nif.sampling25']]$KOtype <- factor(seu.vivo.vitro.list[['Nif.sampling25']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))


set.seed(10)
seu.vitro.nor.list[['Nif.sampling20']] <- subset(seu.vitro.nor.list$Nif,cells = c(seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_DMSO'),20)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Ctrl_20um_Nif'),20)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_DMSO'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_DMSO'),20)],
                                                                                  seu.vitro.nor.list$Nif$SampleName[seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'][sample(1:sum(seu.vitro.nor.list$Nif$Type=='Preg_20um_Nif'),20)]
))
table(seu.vitro.nor.list[['Nif.sampling20']]$Type)


seu.vivo.vitro.list[['Nif.sampling20']] <- merge(seu.ges.beta,seu.vitro.nor.list[['Nif.sampling20']])
seu.vivo.vitro.list[['Nif.sampling20']]@meta.data[as.character(seu.vivo.vitro.list[['Nif.sampling20']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['Nif.sampling20']]$KOtype <- factor(seu.vivo.vitro.list[['Nif.sampling20']]$KOtype,levels = c('/','Ctrl_DMSO','Ctrl_20um_Nif','Preg_DMSO','Preg_20um_Nif'))


saveRDS(seu.vivo.vitro.list,'seu.vitro.nor.list.rds')


seu.Hserum <- readRDS('G:/project/pregnant_mouse/summary_data/mSTRT_data/culture/src.Hserum.beta.rds')

seu.Hserum@assays$RNA@data <- Matrix::Matrix(log(as.matrix(MyCalumiTpmExcludeGene(as.matrix(seu.Hserum@assays$RNA@counts),
                                                                                  c(topgene)))+1),sparse = T)
seu.Hserum$KOtype <- seu.Hserum$Type



seu.vitro.nor.list[['barcodeHSPrl']] <- seu.Hserum
seu.vitro.nor.list[['barcodeHSPrl.all']] <- seu.Hserum

seu.vitro.nor.list$barcodeHSPrl <- subset(seu.vitro.nor.list$barcodeHSPrl,cells = colnames(seu.vitro.nor.list$barcodeHSPrl)[seu.vitro.nor.list$barcodeHSPrl$Type %in% c('human-serum-Ctrl','Prolactin-human-serum-Ctrl')])
seu.vitro.nor.list$barcodeHSPrl.all <- subset(seu.vitro.nor.list$barcodeHSPrl.all,cells = colnames(seu.vitro.nor.list$barcodeHSPrl.all)[seu.vitro.nor.list$barcodeHSPrl.all$Type %in% c('human-serum-Ctrl','human-serum-Preg','Prolactin-human-serum-Ctrl')])

saveRDS(seu.vitro.nor.list[4:5],'direct/20221120/seu.vitro.barcode.list.rds')



tmp.seu <- readRDS('seu.vitro.nor.list.rds')
seu.vitro.nor.list <- readRDS('direct/20221120/seu.vitro.barcode.list.rds')

seu.vitro.nor.list$barcodeHSPrl@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalumiTpmExcludeGene(as.matrix(seu.vitro.nor.list$barcodeHSPrl@assays$RNA@counts[rownames(tmp.seu$Nif),]),
                                                                                  c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list$barcodeHSPrl@assays$RNA@counts <- seu.vitro.nor.list$barcodeHSPrl@assays$RNA@counts[rownames(tmp.seu$Nif),]
seu.vitro.nor.list$barcodeHSPrl$Type <- factor(seu.vitro.nor.list$barcodeHSPrl$Type,levels = c('human-serum-Ctrl','Prolactin-human-serum-Ctrl'))



seu.vitro.nor.list$barcodeHSPrl.all@assays$RNA@data <- Matrix::Matrix(log(MyNorm(as.matrix(MyCalumiTpmExcludeGene(as.matrix(seu.vitro.nor.list$barcodeHSPrl.all@assays$RNA@counts[rownames(seu.vitro.nor.list$barcodeHSPrl),]),
                                                                                                                  c('Ins1','Ins2'))))+1),sparse = T)

seu.vitro.nor.list$barcodeHSPrl.all@assays$RNA@counts <- seu.vitro.nor.list$barcodeHSPrl.all@assays$RNA@counts[rownames(seu.vitro.nor.list$barcodeHSPrl),]
seu.vitro.nor.list$barcodeHSPrl.all$Type <- factor(seu.vitro.nor.list$barcodeHSPrl.all$Type,levels = c('human-serum-Ctrl','human-serum-Preg','Prolactin-human-serum-Ctrl'))

tmp.name <- colnames(seu.vitro.nor.list$barcodeHSPrl)[
seu.vitro.nor.list$barcodeHSPrl$Type=='human-serum-Ctrl'][sample(1:136,74)]
seu.vitro.nor.list$barcodeHSPrl <- subset(seu.vitro.nor.list$barcodeHSPrl,cells=colnames(seu.vitro.nor.list$barcodeHSPrl)[!colnames(seu.vitro.nor.list$barcodeHSPrl)%in%tmp.name])
seu.vitro.nor.list$barcodeHSPrl$Type <- factor(seu.vitro.nor.list$barcodeHSPrl$Type,levels = c('human-serum-Ctrl','Prolactin-human-serum-Ctrl'))

pdf('direct/20221120/HSPrl.Acss2.sampleto62.pdf',5,6.5)
par(mar=c(3,4))
Myseuvioplot(seu.vitro.nor.list$barcodeHSPrl,c('Acss2','Chgb','Ovol2'))
dev.off()


seu.vivo.vitro.list <- list()
seu.vivo.vitro.list[['barcodeHSPrl']] <- merge(seu.ges.beta,seu.vitro.nor.list$barcodeHSPrl)
seu.vivo.vitro.list[['barcodeHSPrl']]@meta.data[as.character(seu.vivo.vitro.list[['barcodeHSPrl']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['barcodeHSPrl']]$KOtype <- factor(seu.vivo.vitro.list[['barcodeHSPrl']]$KOtype,levels = c('/','human-serum-Ctrl','Prolactin-human-serum-Ctrl'))
seu.vivo.vitro.list$barcodeHSPrl@meta.data[seu.vivo.vitro.list$barcodeHSPrl$Type=='human-serum-Ctrl'&  seu.vivo.vitro.list$barcodeHSPrl$Rep %in% c('rep2','rep5'),'Rep'] <- 'rep1'
seu.vivo.vitro.list$barcodeHSPrl@meta.data[seu.vivo.vitro.list$barcodeHSPrl$Type=='human-serum-Ctrl'&seu.vivo.vitro.list$barcodeHSPrl$Rep %in% c('rep3','rep4','rep6'),'Rep'] <- 'rep2'



seu.vivo.vitro.list[['barcodeRS']] <- merge(seu.ges.beta,seu.vitro.nor.list$barcodeFM)
seu.vivo.vitro.list[['barcodeRS']]@meta.data[as.character(seu.vivo.vitro.list[['barcodeRS']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['barcodeRS']]$KOtype <- factor(seu.vivo.vitro.list[['barcodeRS']]$KOtype,levels = c('/','Female_Culture-4d-con','Male_Culture-4d-con','Female_Culture-4d-preg','Male_Culture-4d-preg'))




seu.vivo.vitro.list <- list()
seu.vivo.vitro.list[['barcodeHSPrl.all']] <- merge(seu.ges.beta,seu.vitro.nor.list$barcodeHSPrl.all)
seu.vivo.vitro.list[['barcodeHSPrl.all']]@meta.data[as.character(seu.vivo.vitro.list[['barcodeHSPrl.all']]@meta.data$KOtype)=='/','Rep'] <- 'rep1'
seu.vivo.vitro.list[['barcodeHSPrl.all']]$KOtype <- factor(seu.vivo.vitro.list[['barcodeHSPrl.all']]$KOtype,levels = c('/','human-serum-Ctrl','human-serum-Preg','Prolactin-human-serum-Ctrl'))
seu.vivo.vitro.list$barcodeHSPrl.all@meta.data[seu.vivo.vitro.list$barcodeHSPrl.all$Type=='human-serum-Ctrl'&  seu.vivo.vitro.list$barcodeHSPrl.all$Rep %in% c('rep2','rep5'),'Rep'] <- 'rep1'
seu.vivo.vitro.list$barcodeHSPrl.all@meta.data[seu.vivo.vitro.list$barcodeHSPrl.all$Type=='human-serum-Ctrl'&seu.vivo.vitro.list$barcodeHSPrl.all$Rep %in% c('rep3','rep4','rep6'),'Rep'] <- 'rep2'
seu.vivo.vitro.list$barcodeHSPrl.all@meta.data[seu.vivo.vitro.list$barcodeHSPrl.all$Type=='human-serum-Preg'&seu.vivo.vitro.list$barcodeHSPrl.all$Rep %in% c('rep3','rep4','rep5','rep6'),'Rep'] <- 'rep2'
seu.vivo.vitro.list$barcodeHSPrl.all@meta.data[seu.vivo.vitro.list$barcodeHSPrl.all$Type=='human-serum-Preg'&seu.vivo.vitro.list$barcodeHSPrl.all$Rep %in% c('rep2'),'Rep'] <- 'rep1'


saveRDS(seu.vitro.nor.list$barcodeHSPrl,'barcodeHSPrl.nor.rds')
saveRDS(seu.vitro.nor.list,'direct/20221120/seu.vitro.barcode.nor.list.rds')

seu.vitro.nor.list$barcodeHSPrl$Pseudotime <- 
Mygene2pseudotime(gene.list = c('Acss2','Chgb'),tpm.data = expm1(seu.vitro.nor.list$barcodeHSPrl@assays$RNA@data),log = T,
                  Pseudotime = select.si.list$barcodeHSPrl$pseudotime,Time = select.si.list$barcodeHSPrl$Type,
                  sample.list = select.si.list$barcodeHSPrl$SampleName,cols = time.colors
                  )

Myseuratmarker(seu.vitro.nor.list$barcodeHSPrl,'Cdb1',reduction = 'pca')

seu.vitro.nor.list <- readRDS('seu.vitro.nor.list.rds')
tmp.seu <- seu.vitro.nor.list$sm3.NifA485.all
tmp.seu <- seu.vitro.nor.list$sm3.Nif.A485.only

tmp.seu <- subset(tmp.seu,cells = colnames(seu.vitro.nor.list$sm3.NifA485.all)[seu.vitro.nor.list$sm3.NifA485.all$Type%in%c('Ctrl_DMSO','Preg_DMSO')])
tmp.seu$Type <- factor(tmp.seu$Type,levels = c('Ctrl_DMSO','Preg_DMSO'))

tmp.seu <- subset(tmp.seu,cells = c(colnames(tmp.seu)[tmp.seu$Type=='Ctrl_DMSO'][sample(1:213,70)],
colnames(tmp.seu)[tmp.seu$Type=='Preg_DMSO'][sample(1:216,70)]
))

pdf('direct/20221120/A485Nif.HM.Acss2.all.1.pdf',5,6.5)
Myseuvioplot(tmp.seu,'Acss2')
Myseuvioplot(tmp.seu,'Acss2',log = T)
dev.off()

tmp.seu <- seu.vitro.nor.list$sm3.A485
tmp.p <- wilcox.test(tmp.seu@assays$RNA@data['Acss2',colnames(tmp.seu)[tmp.seu$Type=='Ctrl_DMSO']],
                     tmp.seu@assays$RNA@data['Acss2',colnames(tmp.seu)[tmp.seu$Type=='Preg_DMSO']],paired = F
)
tmp.p$p.value#A485 0.03#Nif 0.99 #A485+Nif 0.07#sample 100 0.0035 #0.009#150 0.001


all.tab <- MyReadDelim('direct/20221120/sm3.NifA485.all.cg.ref.si.tab')
rownames(all.tab) <- all.tab$pseudotime
all.tab <- all.tab[all.tab$Type %in% c('Ctrl_DMSO','Preg_DMSO'),]
all.tab$Type2 <- paste(all.tab$Type,all.tab$SeqDate,sep = '_')
all.tab$Type2 <- as.factor(all.tab$Type2)
all.tab$Type <- all.tab$Type2

pdf('direct/20221120/NS.MS.batch.pdf',10,6)
MyPseudotimebox(all.tab,time.colors = time.colors)
dev.off()

setwd('../../')
seu <- 'barcodeHSPrl'
seu <- 'barcodeHSPrl.all'
seu <- 'sm3.NifA485.all.ex'
seu <- 'Nif.sampling40'
seu <- 'Nif.sampling25'

seu <- 'sm3.NS.MS'

width=17;height=15
for(seu in c('Nif','sm3.A485')
    #names(seu.vivo.vitro.list)
    ){
pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.rep.project.norep.pdf',sep = ''),width,height)

  p.pca <- MySeuratDR2Gg2(seu.vivo.vitro.list[[seu]],seu.vivo.vitro.list[[seu]]@meta.data,reduction.use = 'pca',
                          reduction.key = 'PC',estimate.variation.explain.percentage = T,
                          x.dim = 1,y.dim = 2)
  p.pca$data_ <- p.pca$data
  p.pca$data_$Type_ <- p.pca$data_$KOtype
  p.pca$data_$Type_ <- factor(p.pca$data_$Type_,levels = names(table(seu.vivo.vitro.list[[seu]]$KOtype)))
  p.pca$data <- p.pca$data_[order(p.pca$data_$Type_),]

  plot(p.pca+
         scale_color_manual(values = c('gray80',c('skyblue1','sienna1','deeppink1','purple1'))) +
         scale_shape_manual(values = c(19,17,8,15,1,6))+
         geom_point(aes(x = x.pos,
                        y = y.pos,
                        col = KOtype,
                        shape = Rep
         ),
         size = 6#,#alpha=.6
         )+theme(aspect.ratio = 1)
       #+stat_ellipse(aes(x = x.pos,y = y.pos,col=KOtype2,size=I(2.5)),level=0.7,type='norm')
       )

  dev.off()
}



########
ko.coord.list <- list()
pca.df.list <- list()
seu <- 'barcodeHSPrl.all'
for(seu in names(seu.vivo.vitro.list)[4:5]){
  pca.df.list[[seu]] <- cbind(pca.project.ko.list[[seu]][seu.ges.beta$SampleName,],seu.vivo.vitro.list[[seu]]@meta.data[seu.ges.beta$SampleName,c('SampleName','Type')])
  pca.df.list[[seu]]$Type <- factor(pca.df.list[[seu]]$Type,levels = names(ref.time.colors))
  ko.coord.list[[seu]] <- Mycoord(pca.df.list[[seu]],dims=1:2,
                                  class = 'Type',
                                  names(table(seu.ges.beta$Type)))
  
  regress <- lm(ko.coord.list[[seu]][,2]~
                  ko.coord.list[[seu]][,1])$coefficients[2]
  saveRDS(regress,paste(seu,'.regress.value.rds',sep = ''))
  
  temp=-c(1,regress)%*%rbind(Embeddings(seu.vivo.vitro.list[[seu]],'pca')[,1],
                             Embeddings(seu.vivo.vitro.list[[seu]],'pca')[,2])
  seu.vivo.vitro.list[[seu]]@meta.data[colnames(temp),]$pseudotime <- temp[1,]
}


pdf(paste0(seu,'.cg.lm.pdf'),5,5)
plot(ko.coord.list[[seu]],col=ref.time.colors,pch=20,cex=2,xlim = c(-15,15),ylim = c(-15,10))
abline(lm(ko.coord.list[[seu]][,2]~
            ko.coord.list[[seu]][,1]))
dev.off()


select.si.list <- list()
for(seu in names(seu.vivo.vitro.list)[4:5]){
  select.si.list[[seu]] <- seu.vivo.vitro.list[[seu]]@meta.data[seu.vivo.vitro.list[[seu]]$KOtype!='/',]
  select.si.list[[seu]]$Type <- factor(select.si.list[[seu]]$Type,levels = names(table(seu.vivo.vitro.list[[seu]]$KOtype))[names(table(seu.vivo.vitro.list[[seu]]$KOtype))!='/'])
}

select.si.list[['Nif']] <- MyReadDelim('Nif.cg.ref.si.tab')
select.si.list[['sm3.A485']] <- MyReadDelim('sm3.A485.cg.ref.si.tab')

width=10;height=6
for(seu in names(select.si.list)){
  pdf(paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.cg.ref.project.pseudobox.pdf',sep = ''),width,height)
  select.si.list[[seu]]$pseudotime <- -select.si.list[[seu]]$pseudotime
  print(MyPseudotimebox(select.si.list[[seu]],time.colors = c('skyblue1','sienna1','deeppink1','purple1'),size.point = 4,box.lwd = 2)+theme(aspect.ratio = 1.5))
  dev.off()
}

MyWriteTable(select.si.list[[seu]],'direct/20221120/Nif.A485.NS.MS.pseudotime.si.tab')

#########pval#######
projct.pval.list <- list()
for(seu in names(select.si.list)[4:5]){
  age.list <- names(table(select.si.list[[seu]]$Type))[table(select.si.list[[seu]]$Type)!=0]
  projct.pval.list[[seu]] <- matrix("-",
                                    ncol = length(age.list),
                                    nrow = length(age.list)
  )
  colnames(projct.pval.list[[seu]]) <- age.list
  rownames(projct.pval.list[[seu]]) <- age.list
}

for(seu in names(select.si.list)[4:5]){
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

for(seu in names(select.si.list)[4:5]){
  MyWriteTable(projct.pval.list[[seu]],row.names = T,paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.cg.ref.project.pval.tab',sep = ''))
  MyWriteTable(select.si.list[[seu]],row.names = T,paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.cg.ref.si.tab',sep = ''))
  MyWriteTable(table(select.si.list[[seu]]$Type),row.names = T,paste('G:/lab/Article/heshuang/BYLW/sm3/culture/',seu,'.ref.type.count.tab',sep = ''))
}

#########
saveRDS(pca.project.ko.list,'direct/20221120/cg.pca.project.ko.list.20230511.rds')
rm(seu.vivo.vitro.list,seu.vitro.nor.list)
rm(seu.ges.beta,seu.VG1418,seu.Hserum)
save.image('direct/20221120/direct.20221120.RData')
save.image('direct/20221120/direct.20230511.barcoed.RData')

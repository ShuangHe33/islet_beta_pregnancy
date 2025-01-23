zonation.col <- c('#E83412','#8D3AEB','#0026F7')
names(zonation.col) <- c('PV','Intermediate','CV')
########read file#####
color.image <- rgb(colorRamp(c("white","red"))(seq(0,1,0.01)), maxColorValue = 255)
k4me1.color.image <- rgb(colorRamp(c("white","darkgreen"))(seq(0,1,0.01)), maxColorValue = 255)
atac.color.image <-  rgb(colorRamp(c("white","orange"))(seq(0,1,0.01)), maxColorValue = 255)
p300.color.image <-  rgb(colorRamp(c("white","purple"))(seq(0,1,0.01)), maxColorValue = 255)


pc12.heatmap.col <- colorRampPalette(c( "midnightblue","dodgerblue3","white","darkorange2","red","darkred"), 
                                     space="Lab")(30)
colors.group <- c("#E08698",
                  "#999999",
                  "#A9AF62",
                  "#41BAB0",
                  "#8D9AC6",
                  "#E579E5",
                  # "black",
                  "#ff7f00",
                  "#6a3d9a",
                  "#e31a1c",
                  "#b15928",
                  "#33a02c",
                  "#a6cee3",
                  "#67001f",
                  "#969696"
                  
)
shape.group <- c(16,8,2)

colors.exp <- colorRampPalette(c("#0000FF","white","#FF0000"),
                               space="Lab")(50)

gene.tree.colors <- brewer.pal(8,"Set2")

time.colors <- c("#4DAF4A",
                 # "black",
                 # "#F01414",
                 # "#5F9EA0",#P21-2-B-con
                 "#FF7F00",
                 "#7570B3", 
                 "#e31a1c",
                 "#452cff",#G5.5
                 "#A65628",
                 "#999999",
                 "#F781BF",
                 "#234682",
                 # "#ffe478",#P2-B
                 "#5F9EA0",
                 #"#ff8e9d",
                 #G14.5-HFD,
                 #"#ff723f",
                 "#cab2d6",
                 "#ba984d",
                 "#1f78b4",
                 "#ff4a21",#P7-14-21,P21-B
                 "#752a78",
                 "#a6cee3",
                 "#53751c",
                 "#C71585",
                 "#F01414",
                 "#ff8e9d",
                 "#234682",#G14.5-HFD,
                 "#ff723f",#control culture for 4 days
                 "#d96ca8",#preg culture for 4 days
                 "#0f8575",#G18.5-KO
                 "#5277cf",#Nif-unpreg-Ctrl
                 "#de562c",#Nif-preg-Ctrl
                 "#5a7368",#Nif-preg
                 "#c77fc1",#C646-unpreg-Ctrl,
                 "#3fa0a3",#C646-preg-Ctrl
                 "#c7a96c",#C646-perg
                 "#805e2b",#human serum ctrl
                 "#F01414",
                 "#d15126",
                 #"#6bff9c",#G14.5-2nd
                 #"#5ce6ba"#,
                 # "#39a866",#human serum preg
                 "#6343c4"#human serum GDM
                 
                 
)

ref.time.colors <- time.colors[1:16]
names(ref.time.colors) <- c('Virgin','G3.5','G5.5','G6.5','G8.5','G10.5','G12.5','G14.5','G16.5','G18.5',
                            'P2L','P7L','P14L','P2NL','P7NL','P14NL'
)
MyPlotColor(ref.time.colors,15)

rep.colors <- c('gray90',time.colors[1:4])
names(rep.colors) <- c("*",
                       "rep1",
                       "rep2",
                       "rep3",
                       "rep4"
)
atac.state.col <- rgb(colorRamp(c("#2171b5","black","yellow"))(seq(0,1,0.01)), maxColorValue = 255)
color.tpm <- colorRampPalette(c("#0000FF","white","#FF0000"),
                              space="Lab")(50)
color.atac <- colorRampPalette(c("#3288BD","#EAF79B","#FA834B","#C00000"),
                               space="Lab")(50)
color.image2 <- colorRampPalette(c("#F9E3D2","#F58F66","#DE2E44","#901D5B","#06071C"),
                                 space="Lab")(50)
color.fc <- c("gray60",
              "#6E369E",
              "#FFBE25")
meta.col2 <- c("#33a02c",
               "#e31a1c")
names(meta.col2) <- c("P9",
                      "P18")
meta.col <- c("#e31a1c",
              "#1f78b4")
names(meta.col) <- c("P18",
                     "P60")
box.col2 <- c("#fed976",
              "#fd8d3c",
              "#bd0026")
names(box.col2) <- c("P9",
                     "P18",
                     "P60")
k27ac.col <- c("#a6cee3",
               "#1f78b4",
               "#6495ED"
)
promoter.colors <- c("gray30",
                     "#66c2a5",
                     "#e78ac3",
                     "#fc8d62")

names(promoter.colors) <- c("unmarked",
                            "H3K4me3",
                            "H3K27me3",
                            "bivalent")
library(genomation)
library(GenomicRanges)
tss.bed <- readGeneric("/lustre1/heshuang/work/genome/mm10/mm10.tss.bed",
                       strand = 4,
                       meta.cols = list(gene.id = 5,
                                        trans.id = 6),
                        zero.based = T)
names(tss.bed) <- mcols(tss.bed)$trans.id
# gene.tpm <- read.table("H:/system/E/beta-chip/promoter/2017-12-22 promoter/Bulk_Cell_RNA-seq_Pooled_Normalized_Gene_TPM.txt",
#                        # row.names = 'ID',
#                        header = T,
#                        stringsAsFactors = FALSE
# )
# gene.tpm <- read.delim("G:/lab/gene inf/mm10/Bulk_Cell_RNA-seq_Gene_Read_Count-4.tpm.txt",
#                        header = T,
#                        stringsAsFactors = F)
# 
# row.names(gene.tpm) <- gene.tpm$ID
# sum(is.na((gene.tpm)))
# 
# gene.split.tpm <- read.delim("G:/lab/gene inf/mm10/Bulk_Cell_RNA-seq_Gene_Read_Count-4.tpm.split.txt",
#                        header = T,
#                        stringsAsFactors = F)
# 
# row.names(gene.split.tpm) <- gene.split.tpm$ID

# gene.tpm <- read.delim("G:/project/endocrine_chromatin/beta/analysis/bulk-RNA-seq/merge/beta.tpm.merge.norm.rmP12.1.tab",
#                        header = T,
#                        stringsAsFactors = F)
# row.names(gene.tpm) <- gene.tpm$ID


# bulk_inf <- MyReadDelim("G:/lab/gene inf/mm10/bulk_inf.20161209.tab")
# rownames(bulk_inf) <- bulk_inf$EnsemblGeneID

# gene.parts <-  MyReadTranscriptFeatures("G:/lab/genome/mm10/mm10_convert_ens_filter.refFlat.bed",
#                                         up.flank = 2000,
#                                         down.flank = 2000)

genes.inf.input <- MyReadDelim("/lustre1/heshuang/work/gene_inf/mm10.gene.inf.merge.v1.2.tab")

rownames(genes.inf.input) <- genes.inf.input$SymbolDedu


#hg38.genes.inf.input <- read.csv('G:/lab/gene inf/hg38.gene.inf.merge.v3.csv',stringsAsFactors = F)


# cellcycle.inf <- read.csv("G:/lab/alpha_beta supplyment//cellCycle betacell.csv",
#                           stringsAsFactors = F,
#                           header=T)
# rownames(cellcycle.inf) <- cellcycle.inf$ID
# 
# beta.ss2.immature.gene <- read.csv("G:/lab/alpha_beta supplyment/Cluster 2 genes highly expressed in immature beta-cells.csv",
#                           stringsAsFactors = F,
#                           header = T)
# row.names(beta.ss2.immature.gene) <- beta.ss2.immature.gene$ID
# 
# beta.ss2.mature.gene <- read.csv("G:/lab/alpha_beta supplyment/Cluster 3 genes highly expressed in mature beta-cells.csv",
#                         stringsAsFactors = F,
#                         header = T)
# row.names(beta.ss2.mature.gene) <- beta.ss2.mature.gene$ID
# 
# 
# sc.mature.alpha.gene <- read.csv("G:/lab/alpha_beta supplyment/alpha_ma_high_ens.csv",
#                                  stringsAsFactors = F)
# rownames(sc.mature.alpha.gene) <- sc.mature.alpha.gene$ID
# 
# sc.immature.alpha.gene <- read.csv("G:/lab/alpha_beta supplyment/alpha_imm_high_ens.csv",
#                                  stringsAsFactors = F)
# rownames(sc.immature.alpha.gene) <- sc.immature.alpha.gene$ID


# alpha.sc.tpm.df <- MyReadDelim("G:/lab/genome/mm10/ss2_alpha.sc_tpm/ss2.alpha.sc.tpm.tab")
# rownames(alpha.sc.tpm.df) <- alpha.sc.tpm.df$EnsemblGeneID
# alpha.sc.tpm.norm.df <- MyReadDelim("G:/lab/genome/mm10/ss2_alpha.sc_tpm/ss2.alpha.sc.tpm.norm.tab")
# rownames(alpha.sc.tpm.norm.df) <- alpha.sc.tpm.norm.df$EnsemblGeneID


# alpha.high.tab <- read.csv("G:/lab/alpha and beta supp/gene highly expressed in alpha cells.csv",
#                            stringsAsFactors = F)
# rownames(alpha.high.tab) <- alpha.high.tab$ID
# nrow(alpha.high.tab)#611
# 
# beta.high.tab <- read.csv("G:/lab/alpha and beta supp/gene highly expressed in beta cells.csv",
#                           stringsAsFactors = F)
# rownames(beta.high.tab) <- beta.high.tab$ID
# nrow(beta.high.tab)#762




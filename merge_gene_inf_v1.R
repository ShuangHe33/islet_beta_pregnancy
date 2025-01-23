setwd("~/Volumes/cls_home/lustre1/gene_inf/mouse")

old_inf <- read.delim("mm10.gene.inf.merge.v0.tab",
                      sep = "\t",
                      # skip = 1,
                      # header = F,
                      stringsAsFactors = F)

biomart2 <- read.csv("biomart2.csv",
                     stringsAsFactors = F)
biomart2 <- biomart2[biomart2$Gene.stable.ID %in% old_inf$EnsemblGeneID,]
biomart2 <- biomart2[!duplicated(biomart2$Gene.stable.ID),]
rownames(biomart2) <- biomart2$Gene.stable.ID
biomart2 <- biomart2[!is.na(biomart2$NCBI.gene.ID),]

merge_inf <- data.frame(old_inf[,1:6],
                        NCBIID = "-",
                        old_inf[,7:29],
                        stringsAsFactors = F)
rownames(merge_inf) <- merge_inf$EnsemblGeneID
merge_inf[biomart2$Gene.stable.ID,"NCBIID"] <- biomart2$NCBI.gene.ID

write.table(merge_inf,
            "mm10.gene.inf.merge.v1.tab",
            row.names = F,
            quote = F,
            sep = "\t")


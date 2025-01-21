setwd('../')
dir.create('DEG')
DEG.list <- list()
gene.input <- readRDS('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/beta/Glut2H/ges/gene.input.exp2.10.rds')

seu.beta <- SetIdent(seu.beta,value = seu.beta$Type)
DEG.list[['P7NL']] <- Myseufindmarker(seu.beta,gene.include = gene.input,ident.1 = 'P7NL',ident.2 = 'P7NL_HFD',c1 = 'P7NL',c2 = 'P7NL_HFD')
DEG.list[['Virgin_30HFD']] <- Myseufindmarker(seu.beta,gene.include = gene.input,ident.1 = 'Virgin',ident.2 = 'Virgin_30HFD',c1 = 'Virgin',c2 = 'Virgin_30HFD')


seu <- 'P7NL'
seu <- 'Virgin_30HFD'
#DEG.list.filter <- list()
for(seu in names(DEG.list)){
  DEG.list.filter[[seu]] <- DEG.list[[seu]][DEG.list[[seu]]$p_val_adj<=0.01&
                                                    abs( DEG.list[[seu]]$avg_logFC) >= log(1.5) & 
                                                    ( DEG.list[[seu]]$pct.1>=0.3 |  DEG.list[[seu]]$pct.2>=0.3),]
}



load('../../../../20220726/beta/CD_HFD/seu.beta.ges.CD.HFD.RData')

MyDEGfilterplot(seu.ref.beta,DEG.list.filter[[seu]],DEG.list[[seu]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/P7NL','padj',0.01,'fc',1.5,sep = ''),plotheatmap = F,plotfireplot = T,padj = F,text.cex=1,pval.gap = 10,x.gap = 0.5,point.size = 1.5,fire.cols = ref.time.colors2[c('P7NL','P7L')],ext = 15)
MyDEGfilterplot(seu.ref.beta,DEG.list.filter[[seu]],DEG.list[[seu]],prefix = paste('G:/lab/Article/heshuang/BYLW/sm3/text.P7NL','padj',0.01,'fc',1.5,sep = ''),plotheatmap = F,plotfireplot = T,padj = F,text.cex=1,pval.gap = 10,x.gap = 0.5,point.size = 1.5,fire.cols = ref.time.colors2[c('P7NL','P7L')],ext = 15,y.cex = 0.01)

table(DEG.list.filter[[seu]]$cluster)
save.image('seu.beta.ges.CD.HFD.RData')
P7NL.venn <- Venn(list('P7NL_CD' = DEG.list.filter[[seu]]$gene[DEG.list.filter[[seu]]$cluster=='P7NL'],
          'P7NL_HFD' = DEG.list.filter[[seu]]$gene[DEG.list.filter[[seu]]$cluster=='P7NL_HFD'],
          'preg.high' = preg.gene$SymbolDedu[preg.gene$cluster=='preg.up'],
          'preg.low' = preg.gene$SymbolDedu[preg.gene$cluster=='preg.down']
          ))
plot(P7NL.venn,doWeight = F,type='ellipses')
P7NL.venn@IntersectionSets$`1010`
P7NL.venn@IntersectionSets$`0101`


G14.5HFD.gene.tab <- MyReadDelim('G:/project/pregnant_mouse/beta/sm3/ref/new_sm3/ref_final/HFD/KO/20221203/DEG/gene_inf/padj0.01.fc1.5.G14.5_G14.5_HFD.si.tab')
tmp.venn <- Venn(list('G14.5HFD' = G14.5HFD.gene.tab$gene[G14.5HFD.gene.tab$cluster=='G14.5_HFD'],
          'P7NLHFD' = DEG.list.filter$P7NL$gene[DEG.list.filter$P7NL$cluster=='P7NL_HFD']
          ))
plot(tmp.venn)
boxplot(-G14.5HFD.gene.tab$avg_logFC[G14.5HFD.gene.tab$cluster=='G14.5_HFD'],
        -DEG.list$P7NL[G14.5HFD.gene.tab$gene[G14.5HFD.gene.tab$cluster=='G14.5_HFD'],]$avg_logFC,
        ylab='logFC(HFD/SD)',main = '412genes'
        )
axis(1,1:2,c('G14.5','P7NL'))


Virgin30HFD.venn <- Venn(list('Virgin_30HFD' = DEG.list.filter[[seu]]$gene[DEG.list.filter[[seu]]$cluster=='Virgin_30HFD'],
                       'Virgin_CD' = DEG.list.filter[[seu]]$gene[DEG.list.filter[[seu]]$cluster=='Virgin'],
                       'preg.high' = preg.gene$SymbolDedu[preg.gene$cluster=='preg.up'],
                       'preg.low' = preg.gene$SymbolDedu[preg.gene$cluster=='preg.down']
                       
))
plot(Virgin30HFD.venn,doWeight = F,type='ellipses')


P7NL.row.tree <- MyDEGfilterplot(seu.beta,DEG.raw.filter = DEG.list.filter[[seu]],DEG.raw = DEG.list[[seu]],prefix = 'DEG/P7NLCD_HFD',plotheatmap = T,display.cells = colnames(seu.beta)[seu.beta$Type %in% c('P7NL','P7NL_HFD')],plotRep = T)
Virgin30HFD.row.tree <- MyDEGfilterplot(seu.beta,DEG.raw.filter = DEG.list.filter[[seu]],DEG.raw = DEG.list[[seu]],prefix = 'DEG/Virgin_30HFD',plotheatmap = T,display.cells = colnames(seu.beta)[seu.beta$Type %in% c('Virgin','Virgin_30HFD')],plotRep = T)


dir.create('DEG/go_kegg')
dir.create('DEG/gene_inf')
parame.list <- list()
parame.list[['P7NL']] <- c(0.01,1.5,1.5)
parame.list[["Virgin_30HFD"]] <- c(0.05,1.4,1.4)
sn.DEG.list <- list()
sn.DEG.list[['P7NL']] <- c('P7NL','P7NL_HFD')
sn.DEG.list[["Virgin_30HFD"]] <- c('Virgin',"Virgin_30HFD")
k <- "Virgin_30HFD"
rownames(genes.inf.input) <- sub('_','-',genes.inf.input$SymbolDedu)

for(k in names(DEG.list)){
  # if(seu=="Stat3HMKO"){j <- names(DEG.list.filter)[6];k <- names(DEG.list.filter)[7]}
  # if(seu=="Acss2HMKO"){j <- names(DEG.list.filter)[3];k <- names(DEG.list.filter)[4]}
  MyWriteTable(cbind(genes.inf.input[DEG.list.filter[[k]]$gene,],DEG.list.filter[[k]]),paste('DEG/gene_inf/',seu,'.padj',parame.list[[k]][1],'fc',parame.list[[k]][2],'.si.tab',sep = ''))
  
  for(i in 1:2){
    MyGOwritetable(genes.inf.input[DEG.list.filter[[k]]$gene[DEG.list.filter[[k]]$cluster==sn.DEG.list[[seu]][[i]]],1],pvalue = 1,paste('DEG/go_kegg/go.padj.0.01.fc1.5.',seu,'.',sn.DEG.list[[seu]][[i]],'.high.tab',sep = ''))
    
  }
}
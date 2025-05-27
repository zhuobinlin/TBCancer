library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyverse)

library(survival)
library(survminer)
library(forestplot)
library(forestploter)

main.path='/hwdata/home/linzhuobin/TBCancer'
save.data='/hwdata/home/linzhuobin/TBCancer/TBG_Tumor'

####### Phenotypes of tissue-biased genes #################

## HALLMARK
hgene <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

gene_symbol <- bitr(hgene$entrez_gene, 
                    fromType = "ENTREZID", 
                    toType = c("SYMBOL"), 
                    OrgDb = org.Hs.eg.db)
colnames(gene_symbol) <- c('entrez_gene','gene_symbol')
hgene <- merge(hgene,gene_symbol)

tis_tum<-read.table('/home/linzhuobin/TBCancer/TPM/tissue.list',sep='\t',header=FALSE)
colnames(tis_tum)<-c('tissue','tumor_full','tumor')
tis_tum

tbgene <- '/home/linzhuobin/TBCancer/TPM/TBG_tumor.tsv'
tbgene <- read.table(tbgene, sep='\t',header=TRUE)
rownames(tbgene) <- tbgene$gene

deg<-data.frame()
for (i in 1:nrow(tis_tum)){
  t=tis_tum$tissue[i]
  tt=tis_tum$tumor[i]
  df<-tbgene[tbgene$tissue==t,c('gene','tissue',tt)]
  df<-df[df[3]<0,c('gene','tissue')]
  deg<-rbind(deg,df)
}

hpath <- unique(hgene$gs_name)#50
adf<-data.frame()
for (p in hpath){
  #p='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'
  hg = hgene[hgene$gs_name ==p, 3]
  g<-deg$gene[deg$gene %in% hg]
  df<-data.frame(hallmark=c(p),hgene=c(length(hg)),dtgene=c(length(g)),gene=c(toString(g)))
  adf<-rbind(adf,df)
}

adf$ratio<-adf$dtgene/adf$hgene
adf[order(adf$dtgene,decreasing = T),]
df<-adf[adf$dtgene>=10,]
df$hallmark<-factor(df$hallmark, df[order(df$dtgene),1])

p1<-ggplot(df,aes(hallmark,dtgene,fill=ratio))+
  geom_col(linewidth=0.5,width=0.6)+coord_flip()+#color="black",
  labs(x="",y = "Tissue-biased gene count")+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white", na.value = NA,limits=c(0,0.6),name="Ratio")+
  theme_classic()+
  theme(axis.text = element_text(color='black'))
p1

## SURVIVAL
ssgsea<-read.table('/home/linzhuobin/TBCancer/TBG_Tumor/TBG_ssGSEA.txt',sep='\t',header=T)

surv = read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.survival.tsv'))
surv$sample0 = surv$sample
surv$sample = str_replace_all(surv$sample,'-','.')
surv = surv[grepl('01A',surv$sample0),] # remove Solid Tissue Normal (-11A), duplicated

nn <- surv$sample
nn <- substr(nn, 1, nchar(nn) - 1)
surv$sample <- nn

sg = as.data.frame(t(ssgsea[rownames(ssgsea)=='Liver',]))
colnames(sg) = 'ssGSEA'
sg$sample = rownames(sg)

surv = merge(surv, sg)
cutoff=surv$ssGSEA <= median(surv$ssGSEA)
fit <- survfit(Surv(OS.time, OS) ~cutoff, data = surv) 

ggsurvplot(fit,
           pval = TRUE, pval.method = F, pval.coord= c(0.05, 0.05), pval.size = 5,
           conf.int = TRUE, 
           risk.table = F, 
           legend = c(0.8, 0.9), 
           legend.title = 'Group',legend.labs = c("High", "Low"),font.legend = 14,
           font.main = c(16, "bold", "darkblue"),font.tickslab = 12,
           ylab= "Overall Survival",font.x = 16, font.y =16, 
           ggtheme = theme_survminer(), 
           #palette = c("red", "cyan3") 
           palette = c("#E41A1C","#377EB8"),xlim=c(0,2500)
)

## CNV, SNV, Methylation
cnv = read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.gistic.tsv.gz'))
snv = read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.mutect2_snv.tsv.gz'))
meth= read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.methylation450.tsv.gz'))

alt = read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.alt.tsv')) # gene × sample
alt$sample = str_replace_all(alt$Sample_ID,'-','.')

## for specific gene, such as PCK1 in LIHC
gene = 'PCK1'
gg <- as.data.frame(t(TCGA_TPM[rownames(TCGA_TPM) == gene, ]))
colnames(gg) <- 'gene'
gg$sample <- rownames(gg)

gg = merge(gg, alt[alt$gene == gene, 3:6])

comparisons = list(c('Amplification','Deletion'),c('Deletion','No_CNV'),c('Amplification','No_CNV'))
ggplot(gg, aes(x=CNV, y=gene,fill=CNV))+
  geom_boxplot()+
  scale_fill_manual(values = c("Amplification" = "#E41A1C", "Deletion" = "#377EB8", "No_CNV" = "#4DAF4A"))+
  theme_classic()+ stat_compare_means(method="t.test",comparisons = comparisons,label = "p.signif")+
  ylab(paste0(gene,' expression'))+xlab('')+
  theme(legend.position='none',axis.text = element_text(color='black'))

comparisons = list(c('Hypermethylation','Hypomethylation'),c('Hypermethylation','No_METH'),c('Hypomethylation','No_METH'))
ggplot(gg, aes(x=METH, y=gene,fill=METH))+
  geom_boxplot()+
  scale_fill_manual(values = c("Hypomethylation" = "#E41A1C", "Hypermethylation" = "#377EB8", "No_METH" = "#4DAF4A"))+
  theme_classic()+ stat_compare_means(method="t.test",comparisons = comparisons,label = "p.signif")+
  ylab(paste0(gene,' expression'))+xlab('')+
  theme(legend.position='none',axis.text = element_text(color='black'))

gg1 <- gg
gg1$SNV <- ifelse(gg1$SNV!='No_SNV', 'SNV', 'No_SNV')
ggplot(gg1, aes(x=SNV, y=gene,fill=SNV))+
  geom_boxplot()+
  #scale_fill_brewer(palette="Set1")+
  scale_fill_manual(values = c("SNV" = "#377EB8", "No_SNV" = "#4DAF4A"))+
  theme_classic()+ stat_compare_means(method="t.test",label = "p.signif",label.x = c(1.5))+
  ylab(paste0(gene,' expression'))+xlab('')+
  theme(legend.position='none',axis.text = element_text(color='black'))

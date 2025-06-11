library(data.table)
library(limma)
library(stringr)
library(ggplot2)

library(GSVA)
library(GSEABase)

library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

library(tidyr)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyverse)

library(survival)
library(survminer)
library(forestplot)
library(forestploter)

library(Seurat)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(CellChat)

library(ggVennDiagram)
library(ggvenn)
library(ggpubr)

## Part1
main.path='./'
save.data='./TBG'

####### Identification of tissue-biased genes #################
## Figure 1A, in processing.py

## Load RNA expression matrix, log2(norm_count+1), from https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
TCGA_GTEX_TPM=fread(paste0(main.path,'/Toil/TcgaTargetGtex_RSEM_Hugo_norm_count.gz'),header=T) %>% data.frame()
rownames(TCGA_GTEX_TPM)=TCGA_GTEX_TPM$sample
TCGA_GTEX_TPM=TCGA_GTEX_TPM[,-1]

## Load Phenotype
pheno=read.delim(paste0(main.path,'/Toil/Pheno.txt'),header=T,row.names = 1)

## DEG between tissues
#normal_tissue <- as.data.frame(table(pheno[pheno$X_study == 'GTEX', 'detailed_category']))
#colnames(normal_tissue) <- c('detailed_category','number')
normal_tissue <- sort(unique(pheno[pheno$X_study == 'GTEX', 'detailed_category']))

valid_name <- data.frame(normal_tissue = normal_tissue)
vn <- str_split_i(valid_name$normal_tissue,' \\(',1)
vn <- str_replace_all(vn, ' - ', '_')
vn <- str_replace_all(vn, ' ', '_')
vn <- str_replace_all(vn, '-', '_')
valid_name$valid_name <- vn
write.table(valid_name, paste0(save.data, '/', 'tissue.txt'),quote = F,sep = '\t', row.names =F)

GTEX_TPM <- TCGA_GTEX_TPM[,colnames(TCGA_GTEX_TPM) %in% pheno[pheno$X_study=='GTEX','sample']]

#combinations <- combn(as.character(normal_tissue$detailed_category), 2)
#combinations <- as.data.frame(t(combinations))
#colnames(combinations) <- c('tissue1','tissue2')

n = 1:nrow(valid_name)
for (i in n){
  t1 = valid_name[i,1]
  s1 = pheno[pheno$detailed_category == t1 & pheno$X_study == 'GTEX', 'sample']
  df1= GTEX_TPM[, colnames(GTEX_TPM) %in% s1 ]
  t1 = valid_name[i,2]
  
  n2 = n[n!=i]
  for (j in n2){
    t2 = valid_name[j,1]
    s2 = pheno[pheno$detailed_category == t2 & pheno$X_study == 'GTEX', 'sample']
    df2= GTEX_TPM[, colnames(GTEX_TPM) %in% s2 ]
    t2 = valid_name[j,2]
    
    print(paste0(t1,' vs ',t2))
    
    df = cbind(df1, df2)
    group = c(rep(t1, ncol(df1)),rep(t2,ncol(df2))) %>% factor(., levels = c(t1, t2), ordered = F)
    design <- model.matrix(~0+group)
    colnames(design) <- c(t1, t2)
    rownames(design) <- colnames(df)
    aaa = paste0(t1,'-',t2)
    
    contrast.matrix<-makeContrasts(aaa,levels=design) 
    fit <- lmFit(df,design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2) 
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput)
    write.table(nrDEG,paste0(save.data, '/', paste0(t1,'-VS-',t2), '.result.xls'),
                quote = F,sep = '\t')
    
  }
}

## Find TBGene, 1 vs others
allf <- dir(save.data, pattern = '*.result.xls')
allf <- paste0(save.data, '/', allf)

TBG <- data.frame()
for (i in n){
  t1 = valid_name[i,2]
  tf = allf[grepl(paste0('/',t1), allf)]
  tdf <- data.frame()
  print(t1)
  
  for (f in tf){
    t2 = str_split_i(f, '-VS-', 2)
    t2 = str_replace_all(t2, '.result.xls', '')
    
    df <- read.table(f)
    df <- data.frame( logFC = ifelse(df$adj.P.Val < 0.05, df$logFC, 0),
                      gene  = rownames(df) )
    df$tissue <- t2
    tdf <- rbind(tdf, df)
  }
  
  tdf <- dcast(tdf, gene ~ tissue, value.var = 'logFC')
  rownames(tdf) <- tdf$gene
  tdf <- tdf[,-1]
  
  tbg <- apply(tdf, 1, function(x) sum(x > 1 ))
  tbg <- rownames(tdf)[tbg == ncol(tdf)]
  
  tbg <- data.frame( gene = tbg, tissue = t1)
  TBG <- rbind(TBG, tbg)
}

write.table(TBG,paste0(save.data, '/', 'TBG_all.tsv'),quote = F,sep = '\t')

## TBG list used in paper and database, from http://tsbcancer.canceromics.org/download/TBG_all.tsv
TBG <- read.table(paste0(save.data, '/', 'TBG_all.tsv'), sep='\t',header=T)
geneinfo <- read.table(paste0(save.data, '/', 'genetable.csv'), sep='\t',header=T)
colnames(geneinfo)[3] <- 'gene_symbol'
df <- merge(geneinfo[, 3:4], TBG[, c(1,5)])
df <- as.data.frame(table(df$smple_category, df$gene_type)) # number and gene type
colnames(df) <- c('tissue','gene_type','count')

ggplot(df, aes(fill=gene_type, y=count, x=tissue)) + 
  geom_bar(position='stack', stat='identity') + coord_flip()

####### Expression of tissue-biased genes #################
# from https://xenabrowser.net/datapages/?host=https%3A%2F%2Fgdc.xenahubs.net
save.data='./TBG_Tumor'

## Load Phenotype
pheno=read.delim(paste0(main.path,'/Toil/TcgaTargetGTEX_phenotype.txt.gz'),header=T,row.names = 1)
rownames(pheno)=rownames(pheno) %>% str_replace_all('-','.')
pheno$sample = rownames(pheno)

## Tumor vs NAT, 33 pairs
TCGA_pheno <- pheno[pheno$X_study=='TCGA',]
TCGA_TPM <- TCGA_GTEX_TPM[,colnames(TCGA_GTEX_TPM) %in% TCGA_pheno$sample]

stype <- unique(TCGA_pheno$X_sample_type)
stype <- stype[grepl('Primary|Normal', stype)]
stype <- data.frame(X_sample_type = stype)
stype$sample_group <- ifelse(grepl('Normal', stype$X_sample_type), 'NAT', 'Tumor')

TCGA_pheno <- merge(TCGA_pheno, stype) # exclude Metastatic, Recurrent Tumor...
TCGA_TPM <- TCGA_GTEX_TPM[,colnames(TCGA_GTEX_TPM) %in% TCGA_pheno$sample]
gn <- rownames(TCGA_TPM)
#table(TCGA_pheno$sample_group, TCGA_pheno$detailed_category)

normal_tissue <- sort(unique(TCGA_pheno$detailed_category))

valid_name <- data.frame(normal_tissue = normal_tissue)
vn <- str_split_i(valid_name$normal_tissue,' \\(',1)
vn <- str_replace_all(vn, ' - ', '_')
vn <- str_replace_all(vn, ' ', '_')
vn <- str_replace_all(vn, '-', '_')
valid_name$valid_name <- vn
write.table(valid_name, paste0(save.data, '/', 'tissue.txt'),quote = F,sep = '\t', row.names =F)

## for each tumor 
tumor = 'LIHC'
detailed_pheno=read.delim(paste0(main.path, '/GDC/LIHC/', 'TCGA-LIHC.GDC_phenotype.tsv.gz')) 
nn <- detailed_pheno$submitter_id.samples %>% str_replace_all('-','.')
nn <- substr(nn, 1, nchar(nn) - 1)
detailed_pheno$sample = nn
detailed_pheno = detailed_pheno[!duplicated(nn),]

detailed_stype = stype
colnames(detailed_stype)[1] = 'sample_type.samples'

detailed_pheno <- merge(detailed_pheno, detailed_stype)

#colnames(detailed_pheno)
# 01=Tumor, 11=NAT
nn <- detailed_pheno$sample
nn <- substr(nn, 1, nchar(nn) - 3)
detailed_pheno$sample_pair <- nn

ispair <- as.data.frame(table(nn))
ispair <- ispair[ispair$Freq==2, 'nn']

pair_pheno <- detailed_pheno[detailed_pheno$sample_pair %in% ispair,]
tTCGA_TPM <- t(TCGA_TPM[,colnames(TCGA_TPM) %in% pair_pheno$sample])
tTCGA_TPM <- data.frame(tTCGA_TPM)
tTCGA_TPM$sample <- rownames(tTCGA_TPM)
tTCGA_TPM <- merge(pair_pheno, tTCGA_TPM)

tTCGA_TPM <- tTCGA_TPM[order(tTCGA_TPM$sample_pair, tTCGA_TPM$sample_group),]
rownames(tTCGA_TPM) <- paste0(tTCGA_TPM$sample_pair,':',tTCGA_TPM$sample_group)

df = data.frame(t(tTCGA_TPM[,colnames(tTCGA_TPM) %in% str_replace_all(gn,'-','.')]))
group = ifelse(grepl('.NAT',colnames(df)), 'NAT', 'Tumor') %>% factor(., levels = c('NAT','Tumor'), ordered = T)
pair  = str_replace_all(colnames(df), '.NAT|.Tumor', '')

design <- model.matrix(~0 + group+pair)
colnames(design)[1:length(levels(group))] <- levels(group)
contrast.matrix <- makeContrasts(Tumor - NAT, levels=design)

fit <- lmFit(df, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)

tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.table(nrDEG,paste0(save.data, '/', paste0(tumor, '.result.xls')),quote = F,sep = '\t')

## for each sample
TBG <- read.table(paste0(main.path, '/TBG/', 'TBG_all.tsv'), sep='\t',header=T)
TBG <- TBG[, c(1,5)]
colnames(TBG) <- c('gene','tissue')

gene_sets <- list()
for (t in unique(TBG$tissue)){
  g = TBG[TBG$tissue==t, 'gene']
  gene_sets[[t]] = g
}

gsvaPar <- ssgseaParam(exprData = as.matrix(TCGA_TPM), geneSets = gene_sets, normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose=FALSE)

write.table(ssgsea, paste0(save.data, '/', 'TBG_ssGSEA.txt'),quote = F,sep = '\t', row.names =T)


## Part2
main.path='./'
save.data='./TBG_Tumor'

####### Phenotypes of tissue-biased genes #################

## Figure 1B, HALLMARK
hgene <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

gene_symbol <- bitr(hgene$entrez_gene, 
                    fromType = "ENTREZID", 
                    toType = c("SYMBOL"), 
                    OrgDb = org.Hs.eg.db)
colnames(gene_symbol) <- c('entrez_gene','gene_symbol')
hgene <- merge(hgene,gene_symbol)

tis_tum<-read.table('./TPM/tissue.list',sep='\t',header=FALSE)
colnames(tis_tum)<-c('tissue','tumor_full','tumor')
tis_tum

tbgene <- './TPM/TBG_tumor.tsv'
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

## Figure 1C, SURVIVAL
ssgsea<-read.table('./TBG_Tumor/TBG_ssGSEA.txt',sep='\t',header=T)

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

## Figure 1K, CNV, SNV, Methylation
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


## Part3

main.path='./TBCancer'
save.data='./SC'

####### Single-cell analysis of tissue-biased genes #################
## Figure 1D, clustering and annotation
sce = readRDS('./SC/LIHC/combined.TB.gene.rds')
sce = subset(sce,recluster %in% c('Tumor','Hepatocyte'))
embed <- data.frame(Embeddings(sce, reduction = "umap"))
embed <- subset(embed,UMAP_2 > 0.6)
sce = subset(sce,Cell_ID %in% rownames(embed))

DimPlot(sce, reduction = "umap", label = TRUE)

## Figure 1E, pseudotime
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
new_pd = sce@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(expression_data=data,cell_metadata=new_pd,gene_metadata= fData)
cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="group")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="group")

cds <- cluster_cells(cds)

plot_cells(cds, show_trajectory_graph = FALSE)

cds <- learn_graph(cds)

#saveRDS(cds,'./SC/LIHC/LIHC_monocle3.rds')

p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)

cds <- order_cells(cds)

p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(15,16,0.25))

embed <- data.frame(Embeddings(sce, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.3 & UMAP_1 < -6.1 & UMAP_2 > 15.1 & UMAP_2 < 15.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)

p = plot_cells(cds, color_cells_by = "pseudotime",label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, label_branch_points = F)
p = plot_cells(cds, color_cells_by = "group",label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, label_branch_points = F)
p = plot_cells(cds, color_cells_by = "pseudotime",label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, label_branch_points = F)

a = c('A1BG','F9','AFM')

plot_cells(cds,
           genes='A1BG',
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

plot_genes_branched_heatmap(cds[a,], cluster_rows = T,branch_point = 1,
                            show_rownames = T, use_gene_short_name = T,
                            cores = 1)

## Figure 1F, tissue-biased gene score and stemness score
exp = fread("./SC/LIHC/down.gene.txt")
a = as.data.frame(cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])
a$cell = rownames(a)
names(a)[1] = 'pseudotime'
a = a[order(a$pseudotime),]
a$type = round(a$pseudotime)
a[which(a$type == 0),3] = 1 
b = new_pd[,c(2,65)]
names(b)[1] = 'cell'
a = merge(a,b,by = 'cell')
a = a[order(a$pseudotime),]
anno = data.frame(row.names = a$cell, pseudotime = a$type, Group = a$group)

ct = as.data.frame(sce@assays[["RNA"]]@data)
ct = ct[which(rownames(ct) %in%  exp$gene_symbol),]
ct = ct[exp$gene_symbol,a$cell]

pheatmap(ct,cluster_cols=F,cluster_rows=F,annotation_col = anno,
         scale = 'row',show_colnames = F,show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         breaks = seq(-2,2,by = 0.04),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0)

TB_data = as.data.frame(sce@assays[["AUCell"]]@counts)
TB_data = TB_data[,a$cell]
TB_data = as.data.frame(t(TB_data))
TB_data$cell = rownames(TB_data)

use_plot = merge(TB_data,a,by = 'cell')

plot(use_plot$pseudotime, use_plot$`Liver-TB-gene`)

table(use_plot$type)

use_end = data.frame()
i = 1
for (i in 1:29) {
  use_data = use_plot[which(use_plot$type == i),]
  use_TB = data.frame(pseudotime = i, TB_exp = mean(use_data$`Liver-TB-gene`))
  use_end = rbind(use_end,use_TB)
}

use_end$type = 'TB_gene'
use_end$pseudotime = as.factor(use_end$pseudotime)
use_end$pseudotime <- factor(use_end$pseudotime, levels=c('1','2','3','4','5','6','7','8','9','10','11',
                                                          '12','13','14','15','16','17','18','19','20',
                                                          '21','22','23','24','25','26','27','28','29'))
p = ggplot(use_end) +  theme_classic()+ 
  geom_line(aes(x=pseudotime, y=TB_exp,colour=type,group=type))+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

## Figure 1G, gene intersection
markers <- FindAllMarkers(sce, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
markers3 = markers[which(markers$cluster == 1),]
markers3 = markers3[which(markers3$avg_log2FC < -0),]

aaa = exp[which(exp$gene_symbol %in% markers3$gene),]

TBG <- read.table(paste0(main.path, '/TBG/', 'TBG_all.tsv'), sep='\t',header=T)
TBG <- TBG[, c(1,5)]
colnames(TBG) <- c('gene','tissue')

gene_sets <- list()
for (t in unique(TBG$tissue)){
  g = TBG[TBG$tissue==t, 'gene']
  gene_sets[[t]] = g
}

liver_genes <- gene_sets[['Liver']]
stem_genes <- read.table('./SC/stemnes_gene.list',sep='\t',header=T)
stem_genes <- stem_genes$Gene

sc_data <- AddModuleScore(sce,features = list(Liver = liver_genes),name = "Liver_")
FeaturePlot(sc_data, "Liver_1") + scale_color_gradient(low='#001a49', high='#ecd854')

sc_data <- AddModuleScore(sce,features = list(Stem = stem_genes),name = "Stem_")
FeaturePlot(sc_data, "Stem_1") + scale_color_gradient(low='#001a49', high='#ecd854')

## Figure 1I, cell - cell interaction
sc_data = readRDS('./SC/LIHC/combined.TB.gene.rds')
#table(sc_data@meta.data$recluster2)
cellchat <- createCellChat(object = sc_data, group.by = "recluster2",assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = 'Cell-Cell Contact')

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

netVisual_circle(cellchat@net$count, weight.scale = TRUE, label.edge = FALSE)
netVisual_heatmap(cellchat, color.heatmap = "Reds")

allc = unique(sc_data@meta.data$recluster2)
tumc = allc[grepl('Tumor',allc)]
netVisual_bubble(cellchat, sources.use = tumc,targets.use = c("CD4_T"))


## Figure 1J, screen candidate genes
df <- read.table('./SC/LIHC/screen.txt',sep='\t',header=TRUE)
df$screen_score <- df$screen_score / max(df$screen_score)
df <- df[order(df$screen_score,decreasing = T),]

ggplot(df, aes(x = cor_score, y = screen_score)) +  
  geom_point() + xlim(-15,15) + theme_classic() +
  geom_text(aes(label = ifelse(df$Gene %in% head(df$Gene), df$Gene, '' ),  
    nudge_x = 0.2, check_overlap = TRUE  ))


## Figure 1H, immune infiltration
cibersort_data = fread("./TBG_Tumor/TCGA_CIBERSORT.txt")

df <- as.data.frame(t(ssgsea[rownames(ssgsea) == 'Liver',]))
df$Sample_ID <- rownames(df)
#df$Sample_ID <- str_replace_all(rownames(df),'[.]','-')
df <- merge(df, cibersort_data)

df <- df[df$CancerType=='LIHC',]
df <- df %>% mutate(group = ifelse(Liver > median(Liver),  "High","Low"))
table(df$group)

ggplot(df, aes(x=group, y=B.cells.naive,fill=group))+
  geom_boxplot()+
  scale_fill_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C"))+  
  theme_classic()

df1 <- melt(df[,5:27], id.vars = 'group')
colnames(df1) <- c('group','immune_cell','infiltration')
ggplot(df1, aes(x=immune_cell, y=infiltration,fill=group))+
  geom_boxplot()+
  scale_fill_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C"))+  
  theme_classic()

## Part4

main.path='./TBCancer'
save.data='./experiment'

####### Experiments of tissue-biased genes #################
## Figure 1L, CCK-8
CCK8 <- './experiment/CCK-8.txt'
CCK8 <- read.table( CCK8, header = T)
order_Gene <- c('Ctrl','PCK1','ADH4','HSD17B13','CYP2C8')
CCK8$Gene <- factor(CCK8$Gene, order_Gene)

ggplot(CCK8, aes(Time, OD))+
  geom_point(aes(color = Gene), size = 3)+#, shape = Gene
  geom_line(aes(color = Gene, group = Gene), linewidth = 0.8)+
  scale_x_continuous(breaks=c(0,12,24,36,48,60,72))+labs(x= "Time (h)", y= "OD value")+
  scale_color_manual(values=c('#999999','#57C3F3','#53A85F','#D6E7A3','#F1BB72'))+
  theme_bw()+theme(panel.grid = element_blank())

CCK8p <- './experiment/CCK-8.pvalue.txt'
CCK8p <- read.table( CCK8p, header = T)
CCK8p <- melt(CCK8p, id.vars = c('Time'))
colnames(CCK8p) <- c('Time','Gene','pvalue')
CCK8p$log10pvalue <- -log(CCK8p$pvalue,10)

CCK8f <- './experiment/CCK-8.foldchange.txt'
CCK8f <- read.table( CCK8f, header = T)
CCK8f <- melt(CCK8f, id.vars = c('Time'))
colnames(CCK8f) <- c('Time','Gene','foldchange')

CCK8pf <- merge(CCK8p,CCK8f)
CCK8pf$foldchange <- round(CCK8pf$foldchange,2)
CCK8pf$Time <- factor(CCK8pf$Time,c(72,48,24,12))

ggplot(CCK8pf, aes(Gene, Time, fill = log10pvalue)) + 
  geom_tile(color = "black") +
  geom_text(aes(label = foldchange) , fontface = "bold" ) + 
  scale_fill_gradient(limits=c(-log10(0.05),8),low="#FFC0CB",high ="#CB3B2F",na.value = "grey95")+
  #scale_fill_distiller(palette = "RdBu", limits = c(0, 5), na.value = "grey95") + 
  guides(size = "none") + 
  theme_minimal(base_size = 14) +  
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"))+
  coord_fixed()

## Figure 1M, Transwell
twell <- './experiment/Transwell.txt'
twell <- read.table( twell, header = T)
order_Gene <- c('Ctrl','PCK1','ADH4','HSD17B13','CYP2C8')
twell$Gene <- factor(twell$Gene, order_Gene)

ggplot(twell, aes(x=Gene, y=Value,fill=Gene))+
  geom_bar(position = "dodge",stat = "summary",fun = "mean",color="black")+
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.8,position=position_dodge(0.8))+
  stat_summary(geom='errorbar', width=0.5, fun.min=min, fun.max=max,position=position_dodge(0.8))+
  labs(x= "", y= "Migration cells per field")+ylim(0,2000)+
  scale_fill_manual(values=c('#999999','#57C3F3','#53A85F','#D6E7A3','#F1BB72'))+
  theme_bw()+theme(panel.grid = element_blank(), legend.position = "none")

## Figure 1N, Wound
wound <- './experiment/Wound.txt'
wound <- read.table( wound, header = T)
order_Gene <- c('Ctrl','PCK1','ADH4','HSD17B13','CYP2C8')
wound$Gene <- factor(wound$Gene, order_Gene)
wound$Time <- factor(wound$Time, c(24,48))

ggplot(wound, aes(x=Time, y=Value,fill=Gene))+
  geom_bar(position = "dodge",stat = "summary",fun = "mean",color="black")+
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.8,position=position_dodge(0.9))+
  stat_summary(geom='errorbar', width=0.5, fun.min=min, fun.max=max,position=position_dodge(0.9))+
  labs(x= "Time (h)", y= "Migration distance (μm)")+ylim(0,500)+
  scale_fill_manual(values=c('#999999','#57C3F3','#53A85F','#D6E7A3','#F1BB72'))+
  theme_bw()+theme(panel.grid = element_blank())

## Figure 1O, WB
PCK1 <- './experiment/PCK1_WB.txt'
PCK1 <- read.table( PCK1, header = T)
order_Gene <- c('CD44','CD326','CD133','β-catenin','SOX2','KRT19')
PCK1$Gene <- factor(PCK1$Gene, order_Gene)

ggplot(PCK1, aes(x=Gene, y=Expression,fill=Group))+
  geom_bar(position = "dodge",stat = "summary",fun = "mean",color="black")+
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.8,position=position_dodge(0.8))+
  stat_summary(geom='errorbar', width=0.5, fun.min=min, fun.max=max,position=position_dodge(0.8))+
  labs(x= "", y= "Relative expression")+ylim(0,1.8)+
  scale_fill_manual(values=c('#999999','#59A7D2'))+
  theme_bw()+theme(panel.grid = element_blank())

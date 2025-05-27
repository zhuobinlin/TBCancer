library(data.table)
library(limma)
library(stringr)
library(ggplot2)

library(GSVA)
library(GSEABase)

main.path='/hwdata/home/linzhuobin/TBCancer'
save.data='/hwdata/home/linzhuobin/TBCancer/TBG'

####### Identification of tissue-biased genes #################

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
save.data='/hwdata/home/linzhuobin/TBCancer/TBG_Tumor'

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

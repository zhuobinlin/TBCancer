library(Seurat)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(CellChat)

main.path='/hwdata/home/linzhuobin/TBCancer'
save.data='/hwdata/home/linzhuobin/TBCancer/SC'

####### Single-cell analysis of tissue-biased genes #################
## clustering and annotation
sce = readRDS('/hwdata/home/linzhuobin/TBCancer/SC/LIHC/combined.TB.gene.rds')
sce = subset(sce,recluster %in% c('Tumor','Hepatocyte'))
embed <- data.frame(Embeddings(sce, reduction = "umap"))
embed <- subset(embed,UMAP_2 > 0.6)
sce = subset(sce,Cell_ID %in% rownames(embed))

DimPlot(sce, reduction = "umap", label = TRUE)

## pseudotime
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

#saveRDS(cds,'/hwdata/home/linzhuobin/TBCancer/SC/LIHC/LIHC_monocle3.rds')

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

## tissue-biased gene score and stemness score
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
stem_genes <- read.table('/hwdata/home/linzhuobin/TBCancer/SC/stemnes_gene.list',sep='\t',header=T)
stem_genes <- stem_genes$Gene

sc_data <- AddModuleScore(sce,features = list(Liver = liver_genes),name = "Liver_")
FeaturePlot(sc_data, "Liver_1") + scale_color_gradient(low='#001a49', high='#ecd854')

sc_data <- AddModuleScore(sce,features = list(Stem = stem_genes),name = "Stem_")
FeaturePlot(sc_data, "Stem_1") + scale_color_gradient(low='#001a49', high='#ecd854')

## cell - cell interaction
sc_data = readRDS('/hwdata/home/linzhuobin/TBCancer/SC/LIHC/combined.TB.gene.rds')
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
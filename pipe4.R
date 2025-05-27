library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
#library(ggVennDiagram)
#library(ggvenn)
library(reshape2)
library(ggpubr)
library(tidyverse)

main.path='/hwdata/home/linzhuobin/TBCancer'
save.data='/hwdata/home/linzhuobin/TBCancer/experiment'

####### Experiments of tissue-biased genes #################
## WB
PCK1 <- '/hwdata/home/linzhuobin/TBCancer/experiment/PCK1_WB.txt'
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

## CCK-8
CCK8 <- '/home/linzhuobin/jupyter/TBCancer/experiment/CCK-8.txt'
CCK8 <- read.table( CCK8, header = T)
order_Gene <- c('Ctrl','PCK1','ADH4','HSD17B13','CYP2C8')
CCK8$Gene <- factor(CCK8$Gene, order_Gene)

ggplot(CCK8, aes(Time, OD))+
  geom_point(aes(color = Gene), size = 3)+#, shape = Gene
  geom_line(aes(color = Gene, group = Gene), linewidth = 0.8)+
  scale_x_continuous(breaks=c(0,12,24,36,48,60,72))+labs(x= "Time (h)", y= "OD value")+
  scale_color_manual(values=c('#999999','#57C3F3','#53A85F','#D6E7A3','#F1BB72'))+
  theme_bw()+theme(panel.grid = element_blank())

CCK8p <- '/home/linzhuobin/jupyter/TBCancer/experiment/CCK-8.pvalue.txt'
CCK8p <- read.table( CCK8p, header = T)
CCK8p <- melt(CCK8p, id.vars = c('Time'))
colnames(CCK8p) <- c('Time','Gene','pvalue')
CCK8p$log10pvalue <- -log(CCK8p$pvalue,10)

CCK8f <- '/home/linzhuobin/jupyter/TBCancer/experiment/CCK-8.foldchange.txt'
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

## Transwell
twell <- '/home/linzhuobin/jupyter/TBCancer/experiment/Transwell.txt'
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

## Wound
wound <- '/home/linzhuobin/jupyter/TBCancer/experiment/Wound.txt'
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

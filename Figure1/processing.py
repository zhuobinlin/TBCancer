#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

import os
import itertools
from scipy import stats
import math

import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

indir='/home/linzhuobin/TBCancer/Toil'
outdir='/home/linzhuobin/TBCancer/TBG_Tumor'


mdata='/home/linzhuobin/TBCancer/Toil/TcgaTargetGTEX_phenotype.txt'
mdata=pd.read_csv(mdata, sep='\t', encoding='ISO8859-1')#'gb18030', encoding='ISO8859-1'
mdata.head()

mdata._primary_site.value_counts()

## unified name
mdata['_primary_site']=mdata['_primary_site'].str.replace('Adrenal_gland','Adrenal_Gland')
mdata['_primary_site']=mdata['_primary_site'].str.replace('Cervix_Uteri','Cervix')
mdata['_primary_site']=mdata['_primary_site'].str.replace('Thyroid_Gland','Thyroid')
mdata['_primary_site']=mdata['_primary_site'].str.replace('Breast','Breast_Mammary_Tissue')

## short name
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Adrenocortical_Cancer","ACC")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Bladder_Urothelial_Carcinoma","BLCA")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Brain_Lower_Grade_Glioma","LGG")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Glioblastoma_Multiforme","GBM")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Breast_Invasive_Carcinoma","BRCA")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Cervical_and_Endocervical_Cancer","CESC")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Colon_Adenocarcinoma","COAD")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Esophageal_Carcinoma","ESCA")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Kidney_Chromophobe","KICH")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Kidney_Clear_Cell_Carcinoma","KIRC")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Kidney_Papillary_Cell_Carcinoma","KIRP")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Liver_Hepatocellular_Carcinoma","LIHC")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Lung_Adenocarcinoma","LUAD")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Lung_Squamous_Cell_Carcinoma","LUSC")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Ovarian_Serous_Cystadenocarcinoma","OV")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Pancreatic_Adenocarcinoma","PAAD")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Prostate_Adenocarcinoma","PRAD")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Skin_Cutaneous_Melanoma","SKCM")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Stomach_Adenocarcinoma","STAD")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Testicular_Germ_Cell_Tumor","TGCT")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Thyroid_Carcinoma","THCA")
mdata['primary_disease_or_tissue']=mdata['primary_disease_or_tissue'].str.replace("Uterine_Carcinosarcoma","UCS")


## group GTEx(7429) normal and TCGA normal/tumor (727/9369)
GTEX_mdata=mdata.loc[(mdata['_sample_type']=='Normal_Tissue'),]
TCGA_mdata=mdata.loc[(mdata['_study']=='TCGA'),]
TCGA_tumor_mdata=TCGA_mdata.loc[TCGA_mdata['_sample_type'].str.contains('Primary'),]
TCGA_NAT_mdata=TCGA_mdata.loc[TCGA_mdata['_sample_type'].str.contains('Normal'),]

## previously log 2 (x) transformed
tpm='/home/linzhuobin/TBCancer/Toil/TcgaTargetGtex_RSEM_Hugo_norm_count.gz'
tpm=pd.read_csv(tpm, compression="gzip", sep='\t')
tpm.head()

colnames=tpm.columns.to_list()
tpm.columns=['gene']+colnames[1:]
tpm.head()

colnames=tpm.columns.to_list()
tpm[colnames[1:]]=tpm[colnames[1:]].apply(lambda x: 2**x)
## 1 symbol with N gene_ID
tpm=tpm.groupby(['gene'])[colnames[1:]].mean().reset_index() # 58581
tpm.head()

tbgene='/home/linzhuobin/TBCancer/TBG/TBG_all.tsv'
tbgene=pd.read_csv(tbgene, sep='\t')
tbgene.head()

gidx=tbgene.smple_category.str.contains('Adipose')
adipose=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Artery')
artery=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Brain')
brain=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Cervix')
cervix=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Colon')
colon=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Esophagus')
esophagus=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Heart')
heart=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Kidney')
kidney=tbgene.loc[gidx]['gene_symbol'].to_list()

gidx=tbgene.smple_category.str.contains('Skin')
skin=tbgene.loc[gidx]['gene_symbol'].to_list()

tissue=set(tbgene.smple_category)
#mul_tissue=['Adipose','Artery','Brain','Cervix','Colon','Esophagus','Heart','Kidney','Skin']
mul_tissue=['Adipose','Artery','Brain','Colon','Esophagus','Heart','Kidney','Skin']
len(mul_tissue)

sin_tissue=[]
for t in tissue:
    if not any(mt in t for mt in mul_tissue):
        sin_tissue.append(t)
len(sin_tissue)

tbgene_dict={}
for st in sin_tissue:
    gidx=tbgene.smple_category.str.contains(st)
    tbgene_dict[st]=tbgene.loc[gidx]['gene_symbol'].to_list()
for mt in mul_tissue:
    gidx=tbgene.smple_category.str.contains(mt)
    tbgene_dict[mt]=tbgene.loc[gidx]['gene_symbol'].to_list()
#tbgene_dict

tbgene=pd.DataFrame(columns=['gene','tissue'])
for t,tl in tbgene_dict.items():
    #print(t,len(tl))
    df=pd.DataFrame(tl,columns=['gene'])
    df['tissue']=t
    tbgene=pd.concat([tbgene,df])
tbgene=tbgene.reset_index(drop=True)
tbgene.head()

tbgene.to_csv('/home/linzhuobin/TBCancer/TBG/TBG.tsv',sep='\t',index=None)

tbgene='/home/linzhuobin/TBCancer/TBG/TBG.tsv'
tbgene=pd.read_csv(tbgene, sep='\t')
tbgene.head()

def get_sid(tissue):
    tissue_dict={}
    HEA=GTEX_mdata.loc[GTEX_mdata['_primary_site']==tissue,'sample'].to_list()
    HEA=list(set(HEA)  & set(tpm.columns))
    NAT=TCGA_NAT_mdata.loc[TCGA_NAT_mdata['_primary_site']==tissue,'sample'].to_list()
    TUM=TCGA_tumor_mdata.loc[TCGA_tumor_mdata['_primary_site']==tissue,'sample'].to_list()
    
    print(len(HEA),len(NAT),len(TUM))
    HEA=tpm[['gene']+HEA]
    NAT=tpm[['gene']+NAT]
    TUM=tpm[['gene']+TUM]
    tissue_dict['HEA']=HEA
    tissue_dict['NAT']=NAT
    tissue_dict['TUM']=TUM

    ## TUM types
    tum_df=TCGA_tumor_mdata.loc[TCGA_tumor_mdata['_primary_site']==tissue].groupby('primary_disease_or_tissue')['sample'].apply(list).reset_index()
    for i in list(range(len(tum_df))):
        tt=tum_df.loc[i,'primary_disease_or_tissue']
        sid=tum_df.loc[i,'sample']
        tissue_dict[tt]=tpm[['gene']+sid]
    return tissue_dict


## 18 pair tissue: Healthy, NAT, Tumor
tissue_dict={}
tissue_dict['Adrenal_Gland']=get_sid('Adrenal_Gland')# 128 0 77
tissue_dict['Bladder']=get_sid('Bladder')#9 19 407
tissue_dict['Brain']=get_sid('Brain')#1152 5 662
tissue_dict['Breast_Mammary_Tissue']=get_sid('Breast_Mammary_Tissue')#179 113 1092
tissue_dict['Cervix']=get_sid('Cervix')#10 3 304
tissue_dict['Colon']=get_sid('Colon')#308 41 288
tissue_dict['Esophagus']=get_sid('Esophagus')#653 13 181
tissue_dict['Kidney']=get_sid('Kidney')#28 129 886
tissue_dict['Liver']=get_sid('Liver')#110 50 369
tissue_dict['Lung']=get_sid('Lung')#288 109 1011
tissue_dict['Ovary']=get_sid('Ovary')#88 0 419
tissue_dict['Pancreas']=get_sid('Pancreas')#167 4 178
tissue_dict['Prostate']=get_sid('Prostate')#100 52 495
tissue_dict['Skin']=get_sid('Skin')#556 1 102
tissue_dict['Stomach']=get_sid('Stomach')#174 36 414
tissue_dict['Testis']=get_sid('Testis')#165 0 154
tissue_dict['Thyroid']=get_sid('Thyroid')#279 59 504
tissue_dict['Uterus']=get_sid('Uterus')#78 0 57


for t,tdict in tissue_dict.items():
    os.mkdir('/home/linzhuobin/TBCancer/TPM/'+t)
    tdict['HEA'].to_csv('/home/linzhuobin/TBCancer/TPM/'+t+'/HEA.tsv',sep='\t',index=None)
    tdict['NAT'].to_csv('/home/linzhuobin/TBCancer/TPM/'+t+'/NAT.tsv',sep='\t',index=None)
    tdict['TUM'].to_csv('/home/linzhuobin/TBCancer/TPM/'+t+'/TUM.tsv',sep='\t',index=None)
    for tt in list(tdict.keys())[3:]:
        tdict[tt].to_csv('/home/linzhuobin/TBCancer/TPM/'+t+'/'+tt+'.tsv',sep='\t',index=None)

outdir='/home/linzhuobin/TBCancer/Toil/TPM'

tissue=['Adrenal_Gland','Bladder','Brain','Breast_Mammary_Tissue','Cervix','Colon','Esophagus','Kidney','Liver','Lung',\
        'Ovary','Pancreas','Prostate','Skin','Stomach','Testis','Thyroid','Uterus']#
#tbgene['tissue']=tbgene['tissue'].str.replace('Breast_Mammary_Tissue','Breast')

## sample
def weighted_log2FC(tdf,gene):
    ## TUM fold changes
    #tdf=pd.merge(hdf,tissue_dict[t2]['TUM'],how='left').fillna(0)
    tdf['HEA1']=tdf['HEA']/tdf['HEA'].sum() # -1:last col
    sample=tdf.columns.to_list()[2:-1]
    tdf[['HEA']+sample]=tdf[['HEA']+sample]+0.0000001
    tdf=pd.merge(gene,tdf,how='left').fillna(0)
    tdf[sample]=tdf.iloc[:,2:-1].apply(lambda x: np.log2(x/x[0]), axis=1).iloc[:,1:]
    ## weighted average of local fold changes
    tdf['HEA']=tdf['HEA1']
    tdf[sample]=tdf.iloc[:,2:-1].apply(lambda x: x*x[0], axis=1).iloc[:,1:]*100
    tdf=tdf.iloc[:,3:-1].sum().reset_index()
    tdf.columns=['sample','log2FC']
    return tdf

log2FC_dict={}
for t in tissue:
    print(t)
    tis_df=pd.DataFrame()
    for t2 in tissue:
        ## HEA as baseline
        gene=tbgene.loc[tbgene.tissue==t]
        #df=pd.merge(tbgene,tissue_dict[t2]['HEA'],how='right').fillna('NS')
        df=tissue_dict[t2]['HEA']
        hdf=df.iloc[:,:1]
        hdf['HEA']=df.iloc[:,1:].mean(1)

        ## TUM fold changes
        for tt in list(tissue_dict[t2].keys())[3:]:
            tdf=pd.merge(hdf,tissue_dict[t2][tt],how='left').fillna(0)
            tdf=weighted_log2FC(tdf,gene)
            tdf['tumor']=tt
            tdf['tissue']=t2
            tis_df=pd.concat([tis_df,tdf])
    
    log2FC_dict[t]=tis_df

log2fc=pd.DataFrame()
for t,df in log2FC_dict.items():
    log2fc[t]=df.groupby(['tumor'])['log2FC'].mean().reset_index()['log2FC'].to_list()
log2fc.index=df.groupby(['tumor'])['log2FC'].mean().reset_index()['tumor'].to_list()
log2fc

#log2fc.index.to_list()
tum_rank=['ACC', 'BLCA', 'LGG', 'GBM', 'BRCA', 'CESC', 'COAD', 'ESCA',  'KICH', 'KIRC', 'KIRP',  'LIHC', 
          'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'SKCM', 'STAD', 'TGCT', 'THCA', 'UCS']


df = log2fc.T[tum_rank]
#df.to_csv('/home/linzhuobin/TBCancer/Toil/TPM/tumor.Wlog2FC.tsv',sep='\t',index=True,header=True)
df

## Figure 1A
plt.figure(figsize=(6, 6))
ax=sns.heatmap(data=df,square=True,cmap="RdBu_r",
               center=0,vmin=-1.5,vmax=1,linewidths=0.3,
               cbar_kws={"label":"Weighted Avg. log2FC","shrink":0.5})#,cbar_kws={'shrink':0.5},cbar_kws={"orientation":"horizontal"}
ax.set(xlabel="", ylabel="")
ax.xaxis.tick_top()
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)


tis_tum='/home/linzhuobin/TBCancer/TPM/tissue.list'
tis_tum=pd.read_csv(tis_tum,sep='\t',names=['tissue','tumor_full','tumor'])

tbgene='/home/linzhuobin/TBCancer/TBG/TBG_all.tsv'
tbgene=pd.read_csv(tbgene, sep='\t')
tbgene=tbgene[['gene_symbol','smple_category']]
tbgene.columns=['gene','tissue']
tbgene=pd.merge(tis_tum[['tissue']],tbgene)

tissue=set(tis_tum.tissue)
tumor=['ACC', 'BLCA', 'LGG', 'GBM', 'BRCA', 'CESC', 'COAD', 'ESCA', 'KICH',
       'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'SKCM', 'STAD', 'TGCT', 'THCA', 'UCS']


deg_tum=tbgene.copy()
for i in list(range(len(tis_tum))):
    t=tis_tum.loc[i,'tissue']
    tt=tis_tum.loc[i,'tumor']
    deg=indir+'/'+tt+'.result.xls'
    deg=pd.read_csv(deg,sep='\t')
    deg['type']='NA'
    deg.loc[(deg['logFC'] >= 1 )&(deg['adj.P.Val'] <= 0.05),'type']='Up'
    deg.loc[(deg['logFC'] <=-1 )&(deg['adj.P.Val'] <= 0.05),'type']='Down'
    deg=deg.loc[deg['type']!='NA']
    deg['gene']=deg.index.to_list()
    deg=pd.merge(tbgene,deg,how='left').fillna(0)
    deg_tum[tt]=deg['logFC'].to_list()

deg_tum.to_csv('/home/linzhuobin/TBCancer/TBG_Tumor/TBG_tumor.tsv',sep='\t',index=None)
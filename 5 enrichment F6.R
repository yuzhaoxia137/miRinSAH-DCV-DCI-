#GOKEGG enrichment
#Target genes of significant miRNAs were downloaded from miRWalk (http://129.206.7.150/)
library(tidyverse)


#读取miRWalk下载文件
miR15b5p <- read.csv("miRWalk_miR-15b-5p_Targets.csv",header = T,na.strings = "")
miR215p <- read.csv("miRWalk_miR-21-5p_Targets.csv",header = T,na.strings = "")
miR243p <- read.csv("miRWalk_miR-24-3p_Targets.csv",header = T,na.strings = "")
miR2213p <- read.csv("miRWalk_miR-221-3p_Targets.csv",header = T,na.strings = "")
miR2233p <- read.csv("miRWalk_miR-223-3p_Targets.csv",header = T,na.strings = "")


#筛选在targetscan miRDB miRTarBase三个数据库中均包含的靶基因,并且score>=0.95(miRWalk中filter所推荐)
miR15b5p <- miR15b5p %>% 
  filter((!is.na(validated))&(miRDB==1)&(TargetScan==1)&(bindingp>=0.95))
miR215p <- miR215p %>% 
  filter((!is.na(validated))&(miRDB==1)&(TargetScan==1)&(bindingp>=0.95))
miR243p <- miR243p %>% 
  filter((!is.na(validated))&(miRDB==1)&(TargetScan==1)&(bindingp>=0.95))
miR2213p <- miR2213p %>% 
  filter((!is.na(validated))&(miRDB==1)&(TargetScan==1)&(bindingp>=0.95))
miR2233p <- miR2233p %>% 
  filter((!is.na(validated))&(miRDB==1)&(TargetScan==1)&(bindingp>=0.95))

#合并靶基因
targetgenes <- unique(c(miR15b5p$genesymbol,miR215p$genesymbol,miR243p$genesymbol,miR2213p$genesymbol,miR2233p$genesymbol))

#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
TANSid <- bitr(targetgenes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

#plot of KEGG enrichment
dotplot(enrichKEGG(TANSid$ENTREZID))

#plots of GO enrichment
dotplot(enrichGO(TANSid$ENTREZID, OrgDb = org.Hs.eg.db,ont = "BP"))
dotplot(enrichGO(TANSid$ENTREZID, OrgDb = org.Hs.eg.db,ont = "MF"))
dotplot(enrichGO(TANSid$ENTREZID, OrgDb = org.Hs.eg.db,ont = "CC"))



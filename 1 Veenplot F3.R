#环境信息
# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8
# [4] LC_NUMERIC=C                                LC_TIME=Chinese (Simplified)_China.utf8    
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] meta_5.2-0          GEOquery_2.64.2     Biobase_2.56.0      BiocGenerics_0.42.0 readxl_1.4.0        forcats_0.5.1       stringr_1.4.0      
#  [8] dplyr_1.0.9         purrr_0.3.4         readr_2.1.2         tidyr_1.2.0         tibble_3.1.7        ggplot2_3.3.6       tidyverse_1.3.1    



#加载包
library(tidyverse)
library(VennDiagram)

#加载整理好的miRNAs
#所有miRNA名称去除hsa- - 
#Wang WX2021数据来自于GSE165608
load("miRNAsInEachStudy.Rdata")


#数据变换
data_list <- data_list %>% 
  map(as.vector) %>% 
  map(unlist)

#确定重叠miRNA数量
data_12 <- intersect(data_list[[1]],data_list[[2]])
data_13 <- intersect(data_list[[1]],data_list[[3]])
data_23 <- intersect(data_list[[2]],data_list[[3]])
data_123 <- intersect(data_12,data_list[[3]])

#保存共有miRNA
save(data_123,file = "ComMir.Rdata")

#Venn plot
venn.plot <- draw.triple.venn(
  area1 = length(data_list[[1]]),
  area2 = length(data_list[[2]]),
  area3 = length(data_list[[3]]),
  n12 = length(data_12),
  n13 = length(data_13),
  n23 = length(data_23),
  n123 = length(data_123),
  category = c("Bache S	2017", "Bache S	2020", "Wang WX	2021"),
  fill = c("orange", "red", "green"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green")
)

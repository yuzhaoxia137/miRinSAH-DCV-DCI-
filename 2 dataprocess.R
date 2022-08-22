###加载包
library(tidyverse)
library(readxl)

#读取数据
BS2017_screen <- read_xlsx(path =  "28768799.xlsx", sheet = 1,col_names = T)
BS2017_validation <-  read_xlsx(path =  "28768799.xlsx", sheet = 2,col_names = T)


# 
#加载共有基因
load("ComMir.Rdata")
# 


#加载数据处理流 函数  
#所需格式需要注意：列名为miRNA，分组名称为DCV，内容包括（Yes, No）
dataflow <- function(rawdata,commir = data_123,scale=FALSE){
  #长数据
  rawdata <-  rawdata %>% 
    pivot_longer(cols = !DCV, names_to = "miRNA")
  #是否需要均一化处理
  if (scale) {
    rawdata$value <- scale(rawdata$value,center = TRUE,scale = TRUE)
  }
  
  #改变名称
  rawdata$miRNA <-  rawdata$miRNA %>% 
    str_replace_all("hsa-","") %>% 
    str_replace_all("-","") %>% 
    str_replace_all("\\.\\.\\.[0-9]+","")
  #基因筛选
  rawdata <- rawdata %>% 
    filter(.$miRNA %in% commir)
  #分组计算
  resultname <- rawdata %>% 
    group_by(DCV,miRNA) %>% 
    summarise(m1 = mean(value, na.rm = T),n1=n(), sd1 = sd(value, na.rm = T))
  #设置宽数据
  metaname <- resultname %>% 
    pivot_wider(names_from = DCV, values_from = c(m1,n1,sd1)) %>% 
    rename(mean.e = m1_Yes, sd.e = sd1_Yes, n.e = n1_Yes, mean.c = m1_No, sd.c = sd1_No, n.c = n1_No)
  return(metaname)
}

#函数处理
BS2017_meta_scr <- dataflow(BS2017_screen)
BS2017_meta_vali <- dataflow(BS2017_validation)


#####
#处理32248435数据
BS2020_screen <- read_xlsx(path = "32248435.xlsx",sheet = 1,col_names = T)
BS2020_validation <- read_xlsx(path = "32248435.xlsx",sheet = 2,col_names = T)

#####
#去除Possible人群
BS2020_screen <- BS2020_screen %>% 
  filter(DCV!="Possible")

BS2020_validation <- BS2020_validation %>% 
  filter(DCV!="Possible")

BS2020_meta_scr <- dataflow(BS2020_screen)
BS2020_meta_vali <- dataflow(BS2020_validation,scale = T)

#####
#处理Wang WX 2021数据，选择关联GSE165608数据
#加载包
library(GEOquery)
#####下载数据
a <- getGEO("GSE165608")
#####获取表达数据
exprdata <- exprs(a[[1]])
####提取行名列名
miRNA <- rownames(exprdata)
GSM_ <- colnames(exprdata) 

#获取分组资料
pdata <- pData(a[[1]])
groupdcv <- tibble( GSM = pdata$geo_accession, Group =pdata$`analysis group:ch1`, DCV = pdata$`dcv diagnosis:ch1`)
#转置数据补充行名
dfexprdata <- as.data.frame(t(exprdata))
data4_1 <- rownames_to_column(dfexprdata,var = "GSM") %>% 
  as_tibble()
#合并组和列
data4 <- full_join(data4_1,groupdcv,by ="GSM" )
#按GROUP进行分组  去除多余信息  去除HC健康人群
WWX2021_G1 <- data4 %>% 
  filter(Group=="GroupA") %>% 
  select(-GSM, -Group) %>% 
  select(DCV, everything()) %>% 
  filter(DCV!="HC")

WWX2021_G2 <- data4 %>% 
  filter(Group=="GroupB")%>% 
  select(-GSM, -Group)%>% 
  select(DCV, everything()) %>% 
  filter(DCV!="HC")

WWX2021_G3 <- data4 %>% 
  filter(Group=="GroupC")%>% 
  select(-GSM, -Group)%>% 
  select(DCV, everything()) %>% 
  filter(DCV!="HC")

#数据流处理
WWX2021_meta_G1 <- dataflow(WWX2021_G1)
WWX2021_meta_G2 <- dataflow(WWX2021_G2)
WWX2021_meta_G3 <- dataflow(WWX2021_G3)

#保存数据
save(BS2017_meta_scr,BS2017_meta_vali,BS2020_meta_scr,BS2020_meta_vali,WWX2021_meta_G1,WWX2021_meta_G2,WWX2021_meta_G3,file = "meta_data.Rdata")

save(BS2017_screen,BS2017_validation,BS2020_screen,BS2020_validation,WWX2021_G1,WWX2021_G2,WWX2021_G3,file = "orig_data.Rdata")


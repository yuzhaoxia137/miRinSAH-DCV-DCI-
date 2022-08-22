#加载数据
load(file = "meta_data.Rdata")

#加载包
library(meta)

#合并数据
meta_data <- rbind(BS2017_meta_scr,BS2017_meta_vali,BS2020_meta_scr,BS2020_meta_vali,WWX2021_meta_G1,WWX2021_meta_G2,WWX2021_meta_G3)

#合并分析
m1 <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, n.c = n.c, mean.c = mean.c,sd.c = sd.c,data = meta_data, subgroup = meta_data$miRNA)

#总体亚组分析森林图
forest(m1, overall = F,overall.hetstat = F, subgroup = T, layout ="subgroup",fixed = F, label.e = "DCI",label.c = "non DCI", test.effect.subgroup=T)

#Table 2
print(m1, test.effect.subgroup=T)

#Table 2 计算pvalue

mn <- data.frame()
for (i in data_123) {
  miRpmeta <- meta_data[meta_data$miRNA==i,]
  m5 <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, n.c = n.c, mean.c = mean.c,sd.c = sd.c,data = miRpmeta)
  mn <- append(mn,m5$pval.random)
}

#Table 2 计算publication bias
for (i in data_123) {
  miRpmeta <- meta_data[meta_data$miRNA==i,]
  m5 <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, n.c = n.c, mean.c = mean.c,sd.c = sd.c,data = miRpmeta)
  egg <- metabias(m5,k.min = length(miRpmeta))
  begg <- metabias(m5,k.min = length(miRpmeta),method.bias = "begg")
  print(i)
  print(list(egg$p.value,begg$p.value))
}


#Figure 4 forest plot for miR-21-5p
miR215pmeta <- meta_data[meta_data$miRNA=="miR215p",]
m2 <- metacont(n.e = n.e, mean.e = mean.e, sd.e = sd.e, n.c = n.c, mean.c = mean.c,sd.c = sd.c,data = miR215pmeta)
forest(m2)



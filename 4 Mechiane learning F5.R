#Mechiane learning
#加载数据
library(tidyverse)
load("orig_data.Rdata")

#从meta分析结果中筛选感兴趣基因

sigmiRNA <- c("miR15b5p",
              "miR215p",
              "miR2213p", 
              "miR2233p",
              "miR243p"
)

#处理BS2017数据
#变换miRNA名称
names(BS2017_screen) <- names(BS2017_screen) %>% 
  str_replace_all("(hsa)|-","") 

names(BS2017_validation) <- names(BS2017_validation) %>% 
  str_replace_all("(hsa)|-","") 

#筛选感兴趣miRNAs
BS2017_screen_s <- BS2017_screen %>% 
  select(all_of(sigmiRNA),"DCV")

BS2017_validation_s <- BS2017_validation %>% 
  select(all_of(sigmiRNA),"DCV")

save(BS2017_screen_s,BS2017_validation_s,file = "28768799_ML.Rdata")

#处理BS2020数据，筛选感兴趣miRNAs
BS2020_screen_s <- BS2020_screen %>% 
  select(all_of(sigmiRNA),"DCV")

BS2020_validation_s <- BS2020_validation %>% 
  select(all_of(sigmiRNA),"DCV")

save(BS2020_screen_s,BS2020_validation_s,file = "32248435_ML.Rdata")

#处理WWX2021数据，筛选感兴趣miRNAs
names(WWX2021_G1) <- names(WWX2021_G1) %>% 
  str_replace_all("(hsa)|-","") 
names(WWX2021_G2) <- names(WWX2021_G2) %>% 
  str_replace_all("(hsa)|-","")
names(WWX2021_G3) <- names(WWX2021_G3) %>% 
  str_replace_all("(hsa)|-","") 

#筛选感兴趣miRNA
WWX2021_G1_s <- WWX2021_G1 %>% 
  select(all_of(sigmiRNA),"DCV")
WWX2021_G2_S <- WWX2021_G2 %>% 
  select(all_of(sigmiRNA),"DCV")
WWX2021_G3_s <- WWX2021_G3 %>% 
  select(all_of(sigmiRNA),"DCV")

save(WWX2021_G1_s,WWX2021_G2_S,WWX2021_G3_s,file = "GSE165608_ML.Rdata")

#####准备数据
load("28768799_ML.Rdata")
load("32248435_ML.Rdata")
load("GSE165608_ML.Rdata")

ml_data <- bind_rows(BS2017_screen_s,BS2017_screen_s,BS2020_screen_s,BS2020_validation_s,WWX2021_G1_s,WWX2021_G2_S,WWX2021_G3_s)

ml_data$DCV <- as.factor(ml_data$DCV)

#抽样训练集和验证集
#set.seed("4321")
DCVyes <- ml_data %>% 
  filter(DCV=="Yes")
DCVno <- ml_data %>% 
  filter(DCV=="No")

subtrainyes <- sample(1:dim(DCVyes)[1],dim(DCVyes)[1]*0.7)
trainyes <- DCVyes[subtrainyes,]
testyes <- DCVyes[-subtrainyes,]

subtrainno <- sample(1:nrow(DCVno),nrow(DCVno)*0.7)
trainno <- DCVno[subtrainno,]
testno <- DCVno[-subtrainno,]

train_set <- bind_rows(trainyes,trainno)
test_set <- bind_rows(testyes,testno)

########Decision tree
library(maptree)
library(pROC)
model_1 <- rpart(DCV~., data = train_set,method = "class")
draw.tree(model_1)
model_2 <- rpart(DCV~., data = train_set,method = "class",control = rpart.control(cp=0.011))
draw.tree(model_2)
model_3 <- rpart(DCV~., data = train_set,method = "class",control = rpart.control(cp=0.02))
draw.tree(model_3)
model_4<- rpart(DCV~., data = train_set,method = "class",control = rpart.control(cp=0.031))
draw.tree(model_4)
testresult <- data.frame(DCV= test_set$DCV, predict= predict(model_1,test_set,type = "class"))
table(testresult)
testresult2 <- data.frame(DCV= test_set$DCV, predict= predict(model_2,test_set,type = "class"))
table(testresult2)
testresult3 <- data.frame(DCV= test_set$DCV, predict= predict(model_3,test_set,type = "class"))
table(testresult3)
testresult4 <- data.frame(DCV= test_set$DCV, predict= predict(model_4,test_set,type = "class"))
table(testresult4)

#ROC曲线绘制
ran_roc <- roc(testresult3$DCV,as.numeric(testresult3$predict),ci=T)
plot(ran_roc, print.auc=T, auc.polygon=T, grid=c(0.1,0.2),grid.col = c("green","red"),max.auc.polygon = T, auc.polygon.col="skyblue",print.thres=T)

#######random forest model
library(randomForest)
train_set <- na.omit(train_set)
model_F_1 <- randomForest(train_set$DCV ~ ., data = train_set, importance=TRUE,
                          proximity=TRUE)
pred <- predict(model_F_1,test_set,type = "class")
pred2 <- as.data.frame(pred)
rownames(pred2)[109] <- "109"
testresult5 <- data.frame(DCV= test_set$DCV, predict=pred2)
table(testresult5)

getTree(model_F_1,3,labelVar = T)
plot(model_F_1)
varImpPlot(model_F_1)

#ROC曲线绘制
ran_roc <- roc(testresult5$DCV,as.numeric(testresult5$pred),ci=T)
plot(ran_roc, print.auc=T, auc.polygon=T, grid=c(0.1,0.2),grid.col = c("green","red"),max.auc.polygon = T, auc.polygon.col="skyblue",print.thres=T)


#####Logistic regression model
model_F_2 <- glm(train_set$DCV ~ ., data = train_set,family = binomial())

testresult6 <- data.frame(DCV= test_set$DCV, predict= predict(model_F_2,test_set,type = "link"))
table(testresult6)
logit.pred <- factor(testresult6$predict>.5, levels = c(F,T),labels = c("No","Yes"))
logit.pref <- table(test_set$DCV, logit.pred)
logit.pref

#ROC曲线绘制
ran_roc2 <- roc(testresult6$DCV,as.numeric(logit.pred),ci=T)
plot(ran_roc2, print.auc=T, auc.polygon=T, grid=c(0.1,0.2),grid.col = c("green","red"),max.auc.polygon = T, auc.polygon.col="skyblue",print.thres=T)


####### SVM model
#需要重新加载数列并分组

library(e1071)

model_F_3 <- svm(train_set$DCV ~ ., data = train_set)
testresult7 <- data.frame(DCV= test_set$DCV, predict= predict(model_F_3,test_set,type = "class"))
table(testresult7)

ran_roc3 <- roc(testresult7$DCV,as.numeric(testresult7$predict),ci=T)
plot(ran_roc3, print.auc=T, auc.polygon=T, grid=c(0.1,0.2),grid.col = c("green","red"),max.auc.polygon = T, auc.polygon.col="skyblue",print.thres=T)


#  单层神经网络model
library(nnet)
model_F_4 <- nnet(train_set$DCV ~ ., data = train_set,size=5)
testresult8 <- data.frame(DCV= test_set$DCV, predict= predict(model_F_4,test_set,type = "class"))
table(testresult8)
ran_roc4 <- roc(testresult8$DCV,as.numeric(as.factor(testresult8$predict)),ci=T)
plot(ran_roc4, print.auc=T, auc.polygon=T, grid=c(0.1,0.2),grid.col = c("green","red"),max.auc.polygon = T, auc.polygon.col="skyblue",print.thres=T)

#nnet plot做图
library(NeuralNetTools)
plotnet(model_F_4)


###加载R包
library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE2378', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]
gset2 = getGEO('GSE9944', destdir=".", AnnotGPL = T, getGPL = T)
class(gset2)
gset2[[1]]

#提取子集
plf1<-gset[[1]]@annotation
plf2<-gset2[[1]]@annotation
#提取平台文件
GPL_data<- getGEO(filename ="GPL8300.soft.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)
GPL_data1<- getGEO(filename ="GPL8300.annot.gz", AnnotGPL = T)
GPL_data_22 <- Table(GPL_data1)
#提取表达量
exp <- exprs(gset[[1]])
probe_name<-rownames(exp)
exp2 <- exprs(gset2[[1]])
probe_name2<-rownames(exp2)

###############################################
###########                       #############
###########       数据1转ID       #############
###########                       #############
###############################################
loc<-match(GPL_data_11[,1],probe_name)
probe_exp<-exp[loc,]
raw_geneid<-(as.matrix(GPL_data_11[,"GENE_SYMBOL"]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
gene_exp_matrix=na.omit(gene_exp_matrix)

###############################################
###########                       #############
###########       数据2转ID       #############
###########                       #############
###############################################
loc2<-match(GPL_data_22[,1],probe_name2)
probe_exp2<-exp2[loc2,]
raw_geneid2<-(as.matrix(GPL_data_22[,"Gene symbol"]))
index2<-which(!is.na(raw_geneid2))
geneid2<-raw_geneid2[index2]
exp_matrix2<-probe_exp2[index2,]
geneidfactor2<-factor(geneid2)
gene_exp_matrix2<-apply(exp_matrix2,2,function(x) tapply(x,geneidfactor2,mean))
rownames(gene_exp_matrix2)<-levels(geneidfactor2)
gene_exp_matrix2=na.omit(gene_exp_matrix2)


#数据结合
geo_exp_1=as.data.frame(gene_exp_matrix)
geo_exp_2=as.data.frame(gene_exp_matrix2)
sameSample=intersect(rownames(geo_exp_1), rownames(geo_exp_2))
sameSample=as.data.frame(sameSample)
sameSample=intersect(sameSample$sameSample)
gene_exp1=geo_exp_1[sameSample,,drop=F]
gene_exp2=geo_exp_2[sameSample,,drop=F]
bindgeo=cbind(gene_exp1,gene_exp2)


#####读取分组信息#####
##数据1
pdata <- pData(gset[[1]])
group_list <- ifelse(str_detect(pdata$source_name_ch1,"Primary Optic nerve head astrocytes"), "GLC",
                     "Con")
group_list
group_list = factor(group_list,
                    levels = c("GLC","Con"))
group_list
pdata$group=group_list
##数据2
pdata2 <- pData(gset2[[1]])
group_list2 <- ifelse(str_detect(pdata2$source_name_ch1,"Human optic nerve head astrocytes"), "GLC",
                     "Con")
group_list2
group_list2 = factor(group_list2,
                    levels = c("Con","GLC"))
group_list2
pdata2$group=group_list2


#####分组信息合并#####
group1<-(as.matrix(pdata[,"group"]))
row.names(group1)=rownames(pdata)
colnames(group1)="group"
group2<-(as.matrix(pdata2[,"group"]))
row.names(group2)=rownames(pdata2)
colnames(group2)="group"
talgroup=as.data.frame(rbind(group1,group2))
talgroup_list=factor(talgroup$group,levels = c("Con","GLC"))
write.csv(talgroup,file = "group.csv")

#####进行数据矫正#####
boxplot(bindgeo,outline=T, notch=T,col=talgroup_list, las=2)
dev.off()
bindgeo_normal=normalizeBetweenArrays(bindgeo)
boxplot(bindgeo_normal,outline=T, notch=T,col=talgroup_list, las=2)
range(bindgeo_normal)
bindgeo_normal <- log2(bindgeo_normal+1)
bindgeo_normal <-as.data.frame(bindgeo_normal)
bindgeo_normal=na.omit(bindgeo_normal)
write.csv(bindgeo_normal,file = "bindgeo_exp.csv")
range(bindgeo_normal)
dev.off()

#####进行差异分析#####
design=model.matrix(~talgroup_list)
fit=lmFit(bindgeo_normal,design)#这里要注意，分组的样本与矩阵样本是否相符，不相符则去文件中调整然后读入
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
logFC=1
adj.P.Val = 0.05
k1 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="upanddown.csv")

# 安装和加载所需的包 火山图
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
library(ggplot2)
library(dplyr)

# 导入数据
data <- data.frame(
  logFC = rnorm(100, 0, 1), # 对数变化，随机生成
  pvalue = runif(100, 0, 1) # p值，随机生成
)
data$minusLog10PValue <- -log10(data$pvalue)

# 使用ggplot2绘制火山图
ggplot(data, aes(x = logFC, y = minusLog10PValue)) +  
  geom_point(alpha = 0.5) + 
  theme_minimal() + 
  labs(title = "Volcano Plot",       
       x = "log2(FoldChange)",       
       y = "-log10(P.Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.01), col = "blue", linetype = "dashed")

#绘制top40基因热图
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)  

setwd("D:/biowolf/bioR/17.heatmap")  # 注意路径分隔符在R中通常使用正斜杠（/）
inputFile <- "input.txt"      
groupFile <- "group.txt"  # 分组文件   
outFile <- "heatmap.pdf"     

rt <- read.table(inputFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)     # 读取文件
ann <- read.table(groupFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)    # 读取样本属性文件

pdf(file = outFile, width = 6, height = 5.5)
pheatmap(rt,
         annotation = ann,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_colnames = TRUE,
         scale = "row",  # 矫正
         border_color = NA,  # 此处如果不想显示边框颜色，应该是 NA 而不是 "NA"
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6)
dev.off()

#PCA主成分分析
if (!require("FactoMineR")) install.packages("FactoMineR")
if (!require("factoextra")) install.packages("factoextra")
if (!require("ggplot2")) install.packages("ggplot2")

library(FactoMineR)
library(factoextra)
library(ggplot2)

set.seed(123)  # 设置随机种子以便结果可复现
df <- data.frame(
  Var1 = rnorm(100),
  Var2 = rnorm(100),
  Var3 = rnorm(100),
  Var4 = rnorm(100),
  Group = factor(rep(c("GSE2378", "GSE9944"), each = 50))  #两个组GSE2378和GSE9944
)

pca_result <- PCA(df[, -5], graph = FALSE)  # 排除组别列

# 画出PCA图。使用ggplot2的语法
pca_plot <- fviz_pca_ind(pca_result,
                          col.ind = df$Group,  # 使用组别信息着色
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE,  # 添加置信椭圆
                          ellipse.level = 0.95,  # 置信度为95%
                          legend.title = "Groups"
                          ) + 
             theme_minimal()
             
print(pca_plot)             
             

#样本聚类图（树状图）
library(gplots)
library(stats)

set.seed(123) # 保证可重现性
sample_data <- as.dist(matrix(runif(100), nrow=32))

hc <- hclust(sample_data)     

traits <- matrix(sample.int(2, 32, replace = TRUE) - 1, nrow=2)
rownames(traits) <- c("Con", "GLC")
colnames(traits) <- hc$labels      

heatmap.2(as.matrix(traits),
          Rowv=as.dendrogram(hc),
          Colv=NA, # 不对列进行聚类
          dendrogram="row",
          trace="none", # 不绘制热图的色彩跟踪线
          col=redgreen(2), # 选择颜色
          margin=c(8,4), # 设置边距
          main="Sample dendrogram and trait heatmap")

         
# 尺度无关性图和平均连通度图 
install.packages("WGCNA")
library(WGCNA)
datExpr <- read.csv("path_to_expression_data.csv", row.names = 1)   
options(stringsAsFactors = FALSE)        
powers = c(1:20, seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))   
# 绘制尺度无关性图     
plot(sft$fitIndices[,1], -log(sft$fitIndices[,3]), xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, -log10(p)", type="n", main = "Scale independence")
text(sft$fitIndices[,1], -log(sft$fitIndices[,3]), labels=powers, col="red")
abline(h=-log(0.05), col="red")         

# 绘制平均连通度图
plot(sft$fitIndices[,1], sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,2], labels=powers, col="red")

#基因树状图和模块颜色
#install.packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
    
powers = c(1:20, seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate

adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
                            
moduleColors = labels2colors(dynamicMods)                            
                            
plot(geneTree, main = "Gene dendrogram and module colors",
     xlab = "", sub = "", cex = 0.6)
abline(h = cutHeightForDynamicTree, col = "red")

plot(1:ncol(datExpr), rep(-100, ncol(datExpr)), col = moduleColors, pch = 19, cex = 0.9, xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "Dynamic Tree Cut")

#模块-性状关系热图
#install.packages("WGCNA")
library(WGCNA)

exprData <- read.csv("path_to_gene_expression_data.csv", row.names = 1)
traitData <- read.csv("path_to_trait_data.csv", row.names = 1) 

goodSamples = goodGenes = NULL
goodSamples <- rowSums(is.na(exprData)) == 0
goodGenes <- colSums(is.na(exprData)) == 0
exprData = exprData[goodSamples, goodGenes]
traitData = traitData[goodSamples, ]

powers = c(1:20, seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(exprData, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate

Turn expression data into a TOM (topological overlap matrix)
dissTOM = 1-TOMsimilarityFromExpr(exprData, power = softPower)

geneTree = hclust(as.dist(dissTOM), method = "average")

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE)
dynamicColors = labels2colors(dynamicMods)

moduleTraitCor = cor(dynamicColors, traitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)        

labeledHeatmap(Matrix = moduleTraitCor,
                xLabels = names(traitData),
                yLabels = names(dynamicColors),
                ySymbols = names(dynamicColors),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = ""),
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1,1),
                main = "Module-trait relationships")

#特定模块成员资格（module membership，MM）和基因重要性（gene significance，GS）
library(ggplot2)

et.seed(123) # 确保可重现性
module_membership <- runif(1000, 0, 1)  # 模块成员资格的随机数据
gene_significance <- rbeta(1000, 2, 5)  # 基因重要性的随机数据

correlation <- cor(module_membership, gene_significance)

plot(module_membership, gene_significance, xlab="Module Membership in salmon module", ylab="Gene significance for Treat", 
     main="Module membership vs. gene significance", pch=19, col="red")

text(0.8, 0.8, labels=paste("cor=", round(correlation, 2)), col="black", cex=1.2)

cor.test_result <- cor.test(module_membership, gene_significance)
p_value <- cor.test_result$p.value
text(0.8, 0.75, labels=paste("p<", format.pval(p_value, digits=2)), col="black", cex=1.2)


#VennDiagram WGCNA、Immune、DEGs
install.packages("VennDiagram")
library(VennDiagram)

WGCNA_genes <- sample(1:2000, 2038, replace = FALSE)
DEGs <- sample(1:2000, 409, replace = FALSE)
Immune_genes <- sample(1:2000, 1793, replace = FALSE)

input <- list(
  WGCNA = WGCNA_genes,
  DEGs = DEGs,
  Immune = Immune_genes
)

venn.plot <- venn.diagram(
  x = input,
  category.names = c("WGCNA", "DEGs", "Immune"),
  output = TRUE,
  filename = "venn_diagram.png",
  imagetype = "png",
  height = 768,
  width = 768,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = "black",
  fill = c("skyblue", "pink1", "mediumorchid"),
  alpha = 0.50,
  label.col = "black",
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = c("skyblue", "pink1", "mediumorchid"),
  cat.cex = 2.5,
  cat.pos = 0,
  cat.dist = 0.05,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

grid.draw(venn.plot)

#GO富集分析
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

go_data <- data.frame(
  Term = c("positive regulation of immune system process", "cell activation involved in immune response", ...),
  Count = c(20, 15, ...),
  Ontology = c("BP", "CC", "MF", ...),
  PValue = c(1e-10, 1e-5, ...)
)

go_data$neg_log_pvalue <- -log10(go_data$PValue)
shape_mapping <- c("BP" = 21, "CC" = 24, "MF" = 22)

ggplot(go_data, aes(x = Count, y = Term, size = Count, color = neg_log_pvalue)) +
  geom_point(aes(shape = Ontology), stroke = 2) +  # 根据本体调整形状
  scale_shape_manual(values = shape_mapping) +      # 定义本体形状
  scale_color_gradient(low = "green", high = "red") + # 颜色渐变
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "The Most Enriched GO Terms",
       color = "-log10(pvalue)",
       shape = "ONTOLOGY")

#KEGG分析
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

kegg_data <- data.frame(
  Pathway = c("Central carbon metabolism in cancer", "Epithelial cell signaling in Helicobacter pylori infection", ...),
  GeneNumber = c(10.5, 7.5, ...),
  Count = c(12, 8, ...), # 在这里Count是基因数量
  PValue = c(0.01, 0.02, ...)
)

kegg_data$neg_log_pvalue <- -log10(kegg_data$PValue)
ggplot(kegg_data, aes(x = GeneNumber, y = Pathway, size = Count, color = neg_log_pvalue)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "KEGG Pathway Enrichment",
       x = "Gene Number",
       color = "pvalue")
       
 #基因表达水平的箱线图      
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggsignif", quietly = TRUE)) {
  install.packages("ggsignif")
}

library(ggplot2)
library(ggsignif)

# Gene - 基因名
# Expression - 基因表达值
# Group - 样本分组 (CON 或 GLC)
df <- read.csv("path_to_your_data.csv")

p <- ggplot(df, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = 1) +  # 以点的形式显示异常值
  scale_fill_manual(values = c("blue", "red")) +  # 为不同组手动设置颜色
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + # X轴标签垂直显示
  labs(fill = "Type") + # 设置图例标题

annotations <- data.frame(Gene = c("ACKR3", "ADRB2", ...), pvalue = c(0.01, 0.02, ...))
p + geom_signif(
  data = annotations,
  aes(xmin = Gene, xmax = Gene, annotations = ifelse(pvalue < 0.05, "*", ifelse(pvalue < 0.01, "**", "***")),
  y_position = 8,
  tip_length = 0.5,
  vjust = 0.5
)

print(p)

#模型的偏差与正则化参数（λ）的关系
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}
library(glmnet)

x <- matrix(rnorm(100*20), 100, 20)
y <- sample(c(0,1), 100, replace = TRUE)

cv_fit <- cv.glmnet(x, y, family="binomial")

plot(cv_fit)

#随机森林模型
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}
library(randomForest)

rf <- randomForest(x, y)
importance <- importance(rf)
ordered_importance <- importance[order(importance, decreasing = TRUE),]
barplot(ordered_importance, 
        main="Variable Importance",
        horiz=TRUE, cex.names=0.7)
        
# Venn diagram LASSO and RF      
   if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)     

list_of_sets <- list(
  LASSO = 1:11,
  "Random forest" = 9:15
)

venn.plot <- venn.diagram(
  x = list_of_sets,
  category.names = c("LASSO", "Random forest"),
  output = TRUE,
  filename = "venn_diagram.png",
  imagetype = "png",
  height = 768,
  width = 768,
  resolution = 300,
  compression = "lzw",
  col = "black",
  fill = c("yellow", "green"),
  alpha = 0.50,
  label.col = c("black", "black"),
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = c("yellow", "green"),
  cat.cex = 2.5,
  cat.fontface = "bold"
)

grid.draw(venn.plot)

#散点图矩阵 CD40LG TEK MDK
if (!requireNamespace("GGally", quietly = TRUE)) {
  install.packages("GGally")
}
library(GGally)

df <- data.frame(CD40LG = rnorm(50, 5, 1), MDK = rnorm(50, 5, 1), TEK = rnorm(50, 5, 1))
ggpairs(df, upper = list(continuous = wrap("cor", size = 5)), 
        lower = list(continuous = "points"), 
        diag = list(continuous = "barDiag"))

# CD40LG 箱线图
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

df <- data.frame(Group = rep(c("Con", "GLC"), each = 20), CD40LG = rnorm(40, mean = 6.5, sd = 0.2))
p <- ggplot(df, aes(x = Group, y = CD40LG, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = c("Con" = "blue", "GLC" = "red")) +
  theme_minimal()
  
p_value <- t.test(CD40LG ~ Group, data = df)$p.value  
  
p + annotate("text", x = 1.5, y = max(df$CD40LG) + 0.1, label = paste("p =", format.pval(p_value, digits = 3)))  
print(p)

# TEK 箱线图
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

df <- data.frame(Group = rep(c("Con", "GLC"), each = 20), TEK = rnorm(40, mean = 6.5, sd = 0.2))
p <- ggplot(df, aes(x = Group, y = TEK, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = c("Con" = "blue", "GLC" = "red")) +
  theme_minimal()
  
p_value <- t.test(TEK ~ Group, data = df)$p.value  
  
p + annotate("text", x = 1.5, y = max(df$TEK) + 0.1, label = paste("p =", format.pval(p_value, digits = 3)))  
print(p)

# MDK 箱线图
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

df <- data.frame(Group = rep(c("Con", "GLC"), each = 20), MDK = rnorm(40, mean = 6.5, sd = 0.2))
p <- ggplot(df, aes(x = Group, y = MDK, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = c("Con" = "blue", "GLC" = "red")) +
  theme_minimal()
  
p_value <- t.test(MDK ~ Group, data = df)$p.value  
  
p + annotate("text", x = 1.5, y = max(df$MDK) + 0.1, label = paste("p =", format.pval(p_value, digits = 3)))  
print(p)

#绘制pROC曲线图
install.packages("pROC")

library(pROC)              
inputFile="input.txt"      
outFile="ROC.pdf"         
setwd("D:\\biowolf\\bioR\\41.multiVarROC")              
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)  
y=colnames(rt)[1]

#定义颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
	bioCol=rainbow(ncol(rt))}

#绘制
pdf(file=outFile,width=5,height=5)
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(rt)){
	roc1=roc(rt[,y], as.vector(rt[,i]))
	lines(roc1, col=bioCol[i-1])
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

#绘制列线图
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10), xaxt = 'n', yaxt = 'n')

segments(1, 9, 9, 9)
segments(1, 7, 9, 7)
segments(1, 5, 9, 5)
segments(1, 3, 9, 3)

text(0, 9, "Points", pos = 4)
text(0, 7, "CD40LG", pos = 4)
text(0, 5, "MDK", pos = 4)
text(0, 3, "TEK", pos = 4)
text(0, 1, "Total Points", pos = 4)
text(5, 1, "Risk of Disease", pos = 3)

points(5, 9, pch = 19)
text(5, 9, "0-100", pos = 3)

points(3, 7, pch = 19)
text(3, 7, "2.6-4.6", pos = 3)

points(7, 5, pch = 19)
text(7, 5, "4.4-6.4", pos = 3)

points(2, 3, pch = 19)
text(2, 3, "6.8-4.8", pos = 3)

segments(1, 1, 9, 1)
text(c(2, 4, 6, 8), 1, c("0.1", "0.3", "0.7", "0.99"), pos = 1)

#绘制校准差异图
if (!requireNamespace("rms", quietly = TRUE)) {
  install.packages("rms")
}
library(rms)

fit <- lrm(outcome ~ predictors, data=mydata)
cal <- calibrate(fit, method="boot", B=1000)
 plot(cal)

#绘制分组箱线图
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}
library(ggplot2)
library(ggpubr)

df <- data.frame(Group = rep(c("Con", "GLC"), each = 50), B_cells = rnorm(100), T_cells = rnorm(100), ...)
# 用ggplot2创建箱线图
p <- ggplot(df, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1.5, aes(color = Group)) +
  scale_fill_manual(values = c("Con" = "blue", "GLC" = "red")) +
  theme_minimal() +
  facet_wrap(~variable, scales = 'free', ncol = 4) +
  labs(x = "Group", y = "Estimated proportion") +
  theme(legend.position = "none") # 如果你不想显示图例

p <- p + stat_compare_means(method = "t.test", label = "p.signif")

print(p)

#绘制三个变量（CD40LG、MDK、TEK）和一系列其他变量之间的相关性
if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
}
library(corrplot)

cor_matrix <- cor(df)

# 计算p值矩阵
p_matrix <- matrix(nrow = ncol(df), ncol = ncol(df))
for(i in 1:ncol(df)) {
   for(j in 1:ncol(df)) {
     p_matrix[i,j] <- cor.test(df[,i], df[,j])$p.value
   }
 }
# 创建一个自定义的颜色组
col <- colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200)
# 绘制相关性热图
corrplot(cor_matrix, method = "color", col = col, 
         type = "upper", order = "hclust", 
         addCoef.col = "black", # 添加相关性系数
         # 根据p值添加显著性星号
         p.mat = p_matrix, sig.level = c(.001, .01, .05), insig = "label_sig", 
         tl.col="black", tl.srt=45, # 调整文字颜色和角度
         # 自定义显著性标记
         pch.cex = 0.6, pch.col = "black", cl.lim = c(-0.6, 0.3))

#分组箱线图，用于展示不同细胞类型的估计比例
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

df <- data.frame(cell_type = rep(c("aDCs", "T_cells", "B_cells", ...), times = 20),
estimated_proportion = runif(20 * number_of_cell_types))
# 使用ggplot2来创建箱线图
ggplot(df, aes(x = cell_type, y = estimated_proportion, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +  # 这里选择了一个彩色调色板
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 使X轴标签垂直于轴
        legend.position = "none") +  # 移除图例，因为颜色与X轴标签一一对应
  labs(x = "Cell Type", y = "Estimated Proportion")

#比较不同分组（例如控制组和处理组）之间的数据差异
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggsignif", quietly = TRUE)) {
  install.packages("ggsignif")
}
library(ggplot2)
library(ggsignif)

df <- data.frame(Group = rep(c("Con", "GLC"), each = 100), Variable1 = rnorm(200), Variable2 = rnorm(200), ...)

df_long <- reshape2::melt(df, id.var = "Group")

p <- ggplot(df_long, aes(x = variable, y = value, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("Con" = "blue", "GLC" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转X轴标签以适应

p <- p + geom_signif(comparisons = list(c("Con", "GLC")), 
                     map_signif_level= TRUE,
                     step_increase = 0.05)
print(p)

#未分组的箱线图，用于展示不同细胞类型的估计比例
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
df <- data.frame(cell_type = rep(c("Activated.B.cell", "CD56dim.natural.killer.cell", ...), each = 20),
estimated_proportion = runif(20 * number_of_cell_types))

ggplot(df, aes(x = cell_type, y = estimated_proportion, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +  # 这里选择了一个彩色调色板
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  # 使X轴标签垂直于轴
  labs(x = "Cell Type", y = "Estimated Proportion")

#变量之间相关性的热图
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
}

library(ggplot2)
library(reshape2)
library(corrplot)

cor_matrix <- cor(df)

cor_melted <- reshape2::melt(cor_matrix)

get_significance_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

cor_melted$signif <- mapply(get_significance_annotation, sig_matrix[cbind(cor_melted$Var1, cor_melted$Var2)])

ggplot(cor_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = signif), vjust = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

#GSEA TEK分析
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

set.seed(123) # 设置随机种子，确保结果可重现
gene_list <- sort(rnorm(10000), decreasing = TRUE)
names(gene_list) <- paste0("gene", 1:10000)

selected_genes <- names(gene_list)[1:500]

sample_pathways <- data.frame(
  ID = paste0("KEGG_", sample(LETTERS, 6)),
  Description = c("Alanine, aspartate and glutamate metabolism", "Ascorbate and aldarate metabolism", 
                  "Drug metabolism cytochrome P450", "Focal adhesion", "Hypertrophic cardiomyopathy HCM", 
                  "Metabolism of xenobiotics by cytochrome P450"),
  GeneRatio = sample(1:100, 6),
  BgRatio = sample(1:100, 6),
  pvalue = runif(6, 0, 0.05),
  p.adjust = p.adjust(runif(6, 0, 0.05)),
  qvalue = runif(6, 0, 0.05),
  GeneID = I(replicate(6, sample(selected_genes, sample(20:100, 1))))
)

gsea_plot <- enrichplot::gseaplot2(
  gseaResult = sample_pathways, 
  geneSetID = sample_pathways$ID[5], 
  title = "TEK"
)

print(gsea_plot)

#GSEA CD40LG分析
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

set.seed(123)
gene_list <- sort(rnorm(8000), decreasing = TRUE)
names(gene_list) <- paste0("gene", 1:8000)

selected_genes <- names(gene_list)[1:200]

sample_pathways <- data.frame(
  ID = paste0("KEGG_", sample(LETTERS, 6)),
  Description = c("Autoimmune thyroid disease", "Cytokine-cytokine receptor interaction", 
                  "Hematopoietic cell lineage", "Metabolism of xenobiotics by cytochrome P450", 
                  "Neuroactive ligand-receptor interaction", "Steroid hormone biosynthesis"),
  GeneRatio = sample(1:100, 6),
  BgRatio = sample(1:100, 6),
  pvalue = runif(6, 0, 0.05),
  p.adjust = p.adjust(runif(6, 0, 0.05)),
  qvalue = runif(6, 0, 0.05),
  GeneID = I(replicate(6, sample(selected_genes, sample(20:50, 1))))
)

gsea_plot <- enrichplot::gseaplot(
  gseaResult = sample_pathways, 
  geneSetID = sample_pathways$ID[1], 
  title = "CD40LG"
)

print(gsea_plot)

#GSEA MDK分析
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

set.seed(123)
gene_list <- sort(rnorm(8000), decreasing = TRUE)
names(gene_list) <- paste0("gene", 1:8000)

selected_genes <- names(gene_list)[1:200]

sample_pathways <- data.frame(
  ID = paste0("KEGG_", sample(LETTERS, 6)),
  Description = c("Allograft rejection", "Antigen processing and presentation", 
                  "Autoimmune thyroid disease", "Fatty acid metabolism", 
                  "Graft-versus-host disease", "Neuroactive ligand-receptor interaction"),
  GeneRatio = sample(1:100, 6),
  BgRatio = sample(1:100, 6),
  pvalue = runif(6, 0, 0.05),
  p.adjust = p.adjust(runif(6, 0, 0.05)),
  qvalue = runif(6, 0, 0.05),
  GeneID = I(replicate(6, sample(selected_genes, sample(20:50, 1))))
)

gsea_plot <- enrichplot::gseaplot2(
  gseaResult = sample_pathways, 
  geneSetID = sample_pathways$ID[1], 
  title = "MDK"
)

print(gsea_plot)

#k均值聚类分析
install.packages("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)

set.seed(123) # 保证结果可重现
data <- matrix(rnorm(100), nrow = 50, ncol = 2)

n_runs <- 50
consensus_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))

for (i in 1:n_runs) {
  set.seed(i)
  km <- kmeans(data, centers = 2, nstart = 25)

 for (j in 1:(nrow(data)-1)) {
    for (k in (j+1):nrow(data)) {
      consensus_matrix[j, k] <- consensus_matrix[j, k] + as.integer(km$cluster[j] == km$cluster[k])
      consensus_matrix[k, j] <- consensus_matrix[j, k] # 因为是对称矩阵
    }
  }
}

consensus_matrix <- consensus_matrix / n_runs

Heatmap(consensus_matrix, name = "consensus matrix k=2", col = colorRampPalette(c("blue", "white"))(10))

#累积分布函数（CDF）
install.packages("NbClust")
install.packages("ggplot2")

library(NbClust)
library(ggplot2)

set.seed(123) # 保证结果可重现
data <- matrix(rnorm(100), nrow = 50, ncol = 2)

consensus_index_list <- lapply(2:9, function(k) {
  nb <- NbClust(data, distance = "euclidean", min.nc = k, max.nc = k, method = "kmeans")
  return(nb$consensus)
})

cdf_data <- data.frame()
for (k in 2:9) {
  cdf <- ecdf(consensus_index_list[[k-1]])
  cdf_data <- rbind(cdf_data, data.frame(cdf = cdf(seq(0, 1, length.out = 100)), k = factor(k)))
}

ggplot(cdf_data, aes(x = seq(0, 1, length.out = 100), y = cdf, color = k)) +
  geom_line() +
  labs(title = "consensus CDF", x = "consensus index", y = "CDF") +
  scale_color_discrete(name = "k")

dev.new()
print(ggplot_object)

#比较不同类别（这里是不同的免疫细胞类型）在两个不同集群（A和B）中的分布
library(ggplot2)

set.seed(123) # 确保可重现性
categories <- c('aDCs', 'B_cells', 'CD8_T_cells', 'Cytotoxic_activity', 'DCs', 'HLA', 'Inflammation', 
                'Macrophages', 'Mast_cells', 'Neutrophils', 'NK_CD56bright_cells', 'NK_CD56dim_cells', 
                'pDCs', 'T_cell_co-inhibition', 'T_cell_co-stimulation', 'T_helper_cells', 'TIL', 
                'T_regs', 'Type_I_IFN_Response', 'Type_II_IFN_Response')
cluster <- factor(rep(c('A', 'B'), each = length(categories)))
fraction <- runif(2 * length(categories), min = 0, max = 1)
cell_type <- rep(categories, 2)

data <- data.frame(cluster, cell_type, fraction)

p <- ggplot(data, aes(x = cell_type, y = fraction, fill = cluster)) +
  geom_boxplot(outlier.shape = NA) + # 隐藏异常值
  geom_jitter(shape = 16, position=position_jitter(0.2), aes(color = cluster)) + # 添加数据点
  scale_color_manual(values = c("blue", "red")) + # 设置数据点颜色
  scale_fill_manual(values = c("lightblue", "pink")) + # 设置箱线图填充颜色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # x轴标签竖直排列
  labs(title = "Cluster", y = "Fraction", x = "") + # 设置标签
  theme(legend.position = "none") # 隐藏图例

print(p)

#药物双向网络图
library(igraph)

edges <- data.frame(
  from = c(rep("TEK", 17), rep("CD40LG", 11)),
  to = c("REGORAFENIB", "CABOZANTINIB", "RAZUPROTAFIB", "MGCD-265", "GLESATINIB", "GW559768X", "VANDETANIB", "PEXMETINIB", "FORETINIB", "LOPERAMIDE", "LINIFANIB", "REBASTINIB", "AMPICILLIN", "CETIRIZINE", "ALTIRATINIB", "CE-245677", "CEP-11981",
         "FENOFIBRATE", "TORALIZUMAB", "ROSIGLITAZONE", "PG-102", "ALDESLEUKIN", "SELICRELUMAB", "IMATINIB", "DAPIROLIZUMAB PEGOL", "ATORVASTATIN", "RUPLIZUMAB", "LUCATUMUMAB")
)

graph <- graph_from_data_frame(edges, directed = TRUE)

V(graph)$color <- ifelse(V(graph)$name == "TEK" || V(graph)$name == "CD40LG", "red", "green")
V(graph)$size <- ifelse(V(graph)$name == "TEK" || V(graph)$name == "CD40LG", 30, 15)

plot(graph, 
     edge.arrow.size = 0.5, 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     vertex.frame.color = NA) # 去除节点边框

#microRNAs对这些基因的调控关系
library(igraph)

edges <- data.frame(
  from = c(rep("TEK", 34), rep("MDK", 10), rep("CD40LG", 27)),
  to = paste0("hsa-miR-", sample(100:999, 71, replace = FALSE))
)

graph <- graph_from_data_frame(edges, directed = TRUE)

V(graph)$color <- ifelse(V(graph)$name %in% c("TEK", "MDK", "CD40LG"), "red", "green")
V(graph)$size <- ifelse(V(graph)$name %in% c("TEK", "MDK", "CD40LG"), 30, 10)

plot(graph, 
     edge.arrow.size = 0.5, 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     vertex.frame.color = NA, # 去除节点边框
     vertex.label.dist = 1.5, # 标签与节点的距离
     main = "miRNA-gene regulatory network") # 添加标题
     
#基因相互作用的其他基因或转录因子     
library(igraph)

interaction_data <- data.frame(
  regulator = c(rep("TEK", 6), rep("MDK", 11), rep("CD40LG", 9)),
  target = c("FOXC1", "ARID3A", "FEV", "YY1", "FOXL1", "PRRX2", 
             "HINFP", "E2F6", "ELK4", "ZFX", "RELA", "TFAP2A", "SOX17", "PAX2","SP1", "GATA2",  
             "MAX", "USF1", "STAT1", "NFATC2", "E2F1", "USF2", "GATA3", "FOXL1", "FOXC1")
)

network <- graph_from_data_frame(interaction_data, directed = TRUE)

V(network)$color <- ifelse(V(network)$name %in% c("TEK", "MDK", "CD40LG"), "red", "lightblue")
V(network)$size <- ifelse(V(network)$name %in% c("TEK", "MDK", "CD40LG"), 20, 10)
E(network)$arrow.mode <- 0.5

plot(network, 
     edge.arrow.size = 0.5, 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     vertex.frame.color = NA, # 去除节点边框
     main = "Gene Regulatory Network") # 添加标题
     
layout <- layout_with_fr(network)
plot(network, layout = layout, edge.arrow.size = 0.5)   
     





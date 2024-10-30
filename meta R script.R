setwd("D:/MS")




library(pheatmap)
library(factoextra)
library(ggforce)
library(gridExtra)
library(ggplot2)
library(VennDiagram)
library(dplyr)
library(reshape2)
library(meta)
library(metafor)
library(piano)
library(ggpubr)

library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

library(AnnotationDbi)
library(readxl)
library(GSA)
library(visNetwork)
library(igraph)
library(ggraph)
library(openxlsx)
## df_NAs - comparison across studies 
#NA information saved as a new txt file. 
#write.table(data_NAs, file = "data_NAs.txt", row.names = TRUE, quote = FALSE, sep = "\t" )
data_NA1 = read.delim("D:/MS/(PXD021788)WK001_TMTtxt/data1_PXD021788_NAs.txt")
data_NA2 = read.delim("D:/MS/(PXD038285)MaxQuant_txt_default_2uniquepepts（285）/data2_PXD038285_NAs.txt")
data_NA3 = read.delim("D:/MS/PXD027556/data3_PXD27556_NAs.txt")
data_NA4 = read.delim("D:/MS/(PXD008354)/data4_PXD008354_NAs.txt")
data_NA5 = read.delim("D:/MS/PXD025986/data5_PXD025986_NAs.txt")
data_NA6 = read.delim("D:/MS/(PXD033101)/data6_PXD033101_NAs.txt")
data_NA7 = read.delim("D:/MS/PXD036033/data7_PXD036033_NAs.txt")

data_NAs <- rbind(data_NA1, data_NA2, data_NA3, data_NA4, data_NA5, data_NA6, data_NA7)
melted_NAs=melt(data_NAs, value.name = 'values', variable.name = 'NA_found')

#plotting - melted NAs. x-axis NA_found and fill dataset
pdf(file = "dfFULL_NAfound.pdf", width=10, height = 7)
ggplot(data = melted_NAs, 
       aes(x=NA_found,y=values,fill = Dataset)) +
  geom_bar(aes(fill = Dataset), stat = "identity", position = position_dodge())
dev.off()

#bar graph with NA=1
pdf(file = "dfFULL_1NA_found.pdf", width=10, height = 7)
ggplot(data = data_NAs, 
       aes(x=Dataset,y=NA1,fill = Dataset)) +
  geom_bar(aes(fill = Dataset), stat = "identity", position = position_dodge())
dev.off()


data1_pca = read.delim("D:/MS/(PXD021788)WK001_TMTtxt/data1_melted.txt")
data2_pca = read.delim("D:/MS/(PXD038285)MaxQuant_txt_default_2uniquepepts（285）/data2_melted.txt")
data3_pca = read.delim("D:/MS/PXD027556/data3_melted.txt")
data4_pca = read.delim("D:/MS/(PXD008354)/data4_melted.txt")
data5_pca = read.delim("D:/MS/PXD025986/data5_melted.txt")
data6_pca = read.delim("D:/MS/(PXD033101)/data6_melted.txt")
data7_pca = read.delim("D:/MS/PXD036033/data7_melted.txt")


df_full = data.frame()
# 要合并的数据框列表
data_frames <- list(data1_pca, data2_pca, data3_pca, data4_pca, data5_pca, data6_pca, data7_pca)

# 使用 data2_pca 的列名作为标准
standard_column_names <- colnames(data2_pca)

# 标准化所有数据框的列名
for(i in 1:length(data_frames)){
  colnames(data_frames[[i]]) <- standard_column_names
}

# 合并所有数据框
df_full <- do.call(rbind, data_frames)
write.table(df_full, file = "D:/MS/dfFULL_all.txt", row.names = TRUE, quote = FALSE, sep = "\t" )

#PCA metadata and numeric data prep. 
#Sample distributions.
df_full_wide = dcast(df_full, dataset + sample + group + ocular.region + labeling + instrument ~ protein, value.var="intensities")  #rows are samples, columns are proteins + metadata

metadata_columns=c('dataset','sample','group','ocular.region', 'labeling', 'instrument')   #for PCA, you need to separate metadata from columns that are numeric
temp_meta=df_full_wide[,metadata_columns]
temp_numeric=df_full_wide[,
                          !colnames(df_full_wide) %in% metadata_columns   #columns (prots) of temp_DF NOT IN metadata_columns
]
which(is.na(temp_numeric))

# 用列的均值填充NA值
temp_numeric <- temp_numeric %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

#write.table(temp_numeric, file = "../data/formatted_melted_log2/7Dec_melted_individualsampleinfo/dfFULL_tempnumeric_533variableswithnoNA_exceptCGdata_log2int_filtered_Median_scaled_noAsymAD-EBB-MCI.txt", row.names = TRUE, quote = FALSE, sep = "\t")

# PCA via prcomp 
## get_xxxx - summary 
## fviz_xxxx - visualize the samples and variables on the coordinates. 

PCAdf_full = prcomp (temp_numeric, scale = TRUE)

#PCA:group
pdf(file = "D:/MS/results/PCA_fulldata/PCAgeneral.pdf", width =8, height=6)
fviz_pca_ind(PCAdf_full, 
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-40, 40),
             col.ind = temp_meta$group,
             title = "globalPCA_dfFULL")
dev.off()

#PCA:dataset
pdf(file = "D:/MS/results/PCA_fulldata/PCAdataset.pdf", width =8, height=6)
fviz_pca_ind(PCAdf_full,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-250, 250),
             ylim = c(-100, 100),
             col.ind = temp_meta$dataset,
             title = "globalPCAdataset") + 
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
dev.off()

#PCA:labeling
pdf(file = "D:/MS/results/PCA_fulldata/PCAlabeling.pdf", width =8, height=6)
fviz_pca_ind(PCAdf_full,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-75, 75),
             ylim = c(-40, 40),
             col.ind = temp_meta$labeling,
             title = "globalPCAlabeling")
dev.off()
#
#PCA:ocular.region
pdf(file = "D:/MS/results/PCA_fulldata/PCAocular.region.pdf", width =8, height=6)
fviz_pca_ind(PCAdf_full,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-30, 30),
             col.ind = temp_meta$ocular.region,
             title = "globalPCAocular.region")
dev.off()

#PCA:instrument
pdf(file = "D:/MS/results/PCA_fulldata/PCAinstrument.pdf", width =8, height=6)
fviz_pca_ind(PCAdf_full,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-30, 30),
             col.ind = temp_meta$instrument,
             title = "globalPCAinstrument")
dev.off()




fviz_pca_var(PCAdf_full,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
)
#Results for variables
res.var = get_pca_var(PCAdf_full)
contributions = as.data.frame(res.var$contrib)


# 剔除 data1 并重新进行主成分分析
#----------------------------------------

# Combine data frames excluding data1
data_frames_no_data1 <- list(data2_pca, data3_pca, data4_pca, data5_pca, data6_pca, data7_pca)
standard_column_names <- colnames(data2_pca)

# Standardize column names across all data frames
for(i in 1:length(data_frames_no_data1)){
  colnames(data_frames_no_data1[[i]]) <- standard_column_names
}

# Merge all data frames into a single data frame
df_full_no_data1 <- do.call(rbind, data_frames_no_data1)
write.table(df_full_no_data1, file = "D:/MS/dfFULL_all_no_data1.txt", row.names = TRUE, quote = FALSE, sep = "\t")

# Prepare data for PCA
df_full_wide_no_data1 = dcast(df_full_no_data1, dataset + sample + group + ocular.region + labeling + instrument ~ protein, value.var="intensities")

metadata_columns=c('dataset','sample','group','ocular.region', 'labeling', 'instrument')
temp_meta_no_data1 = df_full_wide_no_data1[, metadata_columns]
temp_numeric_no_data1 = df_full_wide_no_data1[, !colnames(df_full_wide_no_data1) %in% metadata_columns]

# Fill NA values with column means
temp_numeric_no_data1 <- temp_numeric_no_data1 %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Perform PCA
PCAdf_full_no_data1 = prcomp(temp_numeric_no_data1, scale = TRUE)

# Visualize PCA results
pdf(file = "D:/MS/results/PCA_fulldata/PCAgeneral_no_data1.pdf", width = 8, height = 6)
fviz_pca_ind(PCAdf_full_no_data1, 
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-40, 40),
             col.ind = temp_meta_no_data1$group,
             title = "globalPCA_dfFULL_no_data1")
dev.off()

pdf(file = "D:/MS/results/PCA_fulldata/PCAdataset_no_data1.pdf", width = 8, height = 6)
fviz_pca_ind(PCAdf_full_no_data1,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-60, 60),
             ylim = c(-60, 60),
             col.ind = temp_meta_no_data1$dataset,
             title = "globalPCAdataset_no_data1") + 
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
dev.off()

pdf(file = "D:/MS/results/PCA_fulldata/PCAlabeling_no_data1.pdf", width = 8, height = 6)
fviz_pca_ind(PCAdf_full_no_data1,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-75, 75),
             ylim = c(-40, 40),
             col.ind = temp_meta_no_data1$labeling,
             title = "globalPCAlabeling_no_data1")
dev.off()

pdf(file = "D:/MS/results/PCA_fulldata/PCAocular.region_no_data1.pdf", width = 8, height = 6)
fviz_pca_ind(PCAdf_full_no_data1,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-30, 30),
             col.ind = temp_meta_no_data1$ocular.region,
             title = "globalPCAocular.region_no_data1")
dev.off()

pdf(file = "D:/MS/results/PCA_fulldata/PCAinstrument_no_data1.pdf", width = 8, height = 6)
fviz_pca_ind(PCAdf_full_no_data1,
             geom.ind = c("point"),
             pointsize = 1.5,
             addEllipses = TRUE,
             xlim = c(-90, 90),
             ylim = c(-30, 30),
             col.ind = temp_meta_no_data1$instrument,
             title = "globalPCAinstrument_no_data1")
dev.off()

# Visualize variable contributions to the PCA
fviz_pca_var(PCAdf_full_no_data1,
             col.var = "contrib", # Color by contributions to the PC 按照变量对主成分的贡献着色
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping

# Results for variables
res.var_no_data1 = get_pca_var(PCAdf_full_no_data1)
contributions_no_data1 = as.data.frame(res.var_no_data1$contrib)

#PCA bioplot.
#```{r}
### 10 highest contributions on dim1
dim1_10highest=row.names(contributions[order(contributions$Dim.1, decreasing = T),])[1:10]
### 10 highest on dim2
dim2_10highest=row.names(contributions[order(contributions$Dim.2, decreasing = T),])[1:10]

pdf(file = "D:/MS/results/PCA_fulldata/variableplot.pdf", width =8, height=6)
fviz_pca_biplot(
  PCAdf_full,
  col.var = "contrib", # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = T,     # Avoid text overlapping
  geom.ind = c("point"),
  pointsize = 1.5,
  col.ind = "grey30",
  addEllipses = F,
  xlim = c(-30, 30),
  ylim = c(-40, 40),
  select.var = list(name=c(dim1_10highest, dim2_10highest))
)
dev.off()


# Random-effects-model  
## Data prep
#use the formatted data called dataX_log2: 20% NA is allowed per group, log2 conversion, median and scale. 
# 读取归一化后的txt文件并存储为data_filtered数据框
data1_input <- read.table("D:/MS/(PXD021788)WK001_TMTtxt/data1_filtered.txt", header = TRUE)  
data2_input <- read.table("D:/MS/(PXD038285)MaxQuant_txt_default_2uniquepepts（285）/data2_filtered.txt", header = TRUE) 
data3_input <- read.table("D:/MS/PXD027556/data3_filtered.txt", header = TRUE) 
data4_input <- read.table("D:/MS/(PXD008354)/data4_PXD008354_filtered.txt", header = TRUE)
data5_input <- read.table("D:/MS/PXD025986/data5_PXD025986_filtered.txt", header = TRUE)
data6_input <- read.table("D:/MS/(PXD033101)/data6_PXD033101_filtered.txt", header = TRUE)
data7_input <- read.table("D:/MS/PXD036033/data7_PXD036033_filtered.txt", header = TRUE)

# data1
data1_input_NC <- data1_input[, 1:3]
data1_mean_NC <- rowMeans(data1_input_NC, na.rm = TRUE)
data1_sd_NC <- apply(data1_input_NC, 1, FUN = sd, na.rm = TRUE)
data1_n_NC <- rowSums(!is.na(data1_input_NC))
data1input <- data.frame(
  Mean_NC = data1_mean_NC,
  SD_NC = data1_sd_NC,
  N_NC = data1_n_NC
)

data1_input_PDR <- data1_input[, 4:6]
mean_PDR <- rowMeans(data1_input_PDR, na.rm = TRUE)
sd_PDR <- apply(data1_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR <- rowSums(!is.na(data1_input_PDR))
data1input_PDR <- data.frame(
  Mean_PDR = mean_PDR,
  SD_PDR = sd_PDR,
  N_PDR = n_PDR
)
# 合并NC组和PDR组的数据
data1_input <- cbind(data1input, data1input_PDR)
data1_input$reference <- "data1"


# data2
data2_input_NC <- data2_input[, 1:4]
data2_input_PDR <- data2_input[, 5:8]
mean_NC_data2 <- rowMeans(data2_input_NC, na.rm = TRUE)
sd_NC_data2 <- apply(data2_input_NC, 1, FUN = sd, na.rm = TRUE)
n_NC_data2 <- rowSums(!is.na(data2_input_NC))
mean_PDR_data2 <- rowMeans(data2_input_PDR, na.rm = TRUE)
sd_PDR_data2 <- apply(data2_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR_data2 <- rowSums(!is.na(data2_input_PDR))
data2_input <- data.frame(
  Mean_NC = mean_NC_data2,
  SD_NC = sd_NC_data2,
  N_NC = n_NC_data2,
  Mean_PDR = mean_PDR_data2,
  SD_PDR = sd_PDR_data2,
  N_PDR = n_PDR_data2
)
data2_input$reference <- "data2"

# data3
data3_input_NC <- data3_input[, 1:3]
data3_input_PDR <- data3_input[, 4:8]
mean_NC_data3 <- rowMeans(data3_input_NC, na.rm = TRUE)
sd_NC_data3 <- apply(data3_input_NC, 1, FUN = sd, na.rm = TRUE)
n_NC_data3 <- rowSums(!is.na(data3_input_NC))
mean_PDR_data3 <- rowMeans(data3_input_PDR, na.rm = TRUE)
sd_PDR_data3 <- apply(data3_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR_data3 <- rowSums(!is.na(data3_input_PDR))
data3_input <- data.frame(
  Mean_NC = mean_NC_data3,
  SD_NC = sd_NC_data3,
  N_NC = n_NC_data3,
  Mean_PDR = mean_PDR_data3,
  SD_PDR = sd_PDR_data3,
  N_PDR = n_PDR_data3
)
data3_input$reference <- "data3"


#data4
data4_input_NC <- data4_input[, 1:9]
data4_mean_NC <- rowMeans(data4_input_NC, na.rm = TRUE)
data4_sd_NC <- apply(data4_input_NC, 1, FUN = sd, na.rm = TRUE)
data4_n_NC <- rowSums(!is.na(data4_input_NC))
data4input <- data.frame(
  Mean_NC = data4_mean_NC,
  SD_NC = data4_sd_NC,
  N_NC = data4_n_NC
)

data4_input_PDR <- data4_input[, 10:19]
mean_PDR <- rowMeans(data4_input_PDR, na.rm = TRUE)
sd_PDR <- apply(data4_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR <- rowSums(!is.na(data4_input_PDR))
data4input_PDR <- data.frame(
  Mean_PDR = mean_PDR,
  SD_PDR = sd_PDR,
  N_PDR = n_PDR
)

data4_input <- cbind(data4input, data4input_PDR)
data4_input$reference <- "data4"

#data5
data5_input_NC <- data5_input[, 1:10]
data5_mean_NC <- rowMeans(data5_input_NC, na.rm = TRUE)
data5_sd_NC <- apply(data5_input_NC, 1, FUN = sd, na.rm = TRUE)
data5_n_NC <- rowSums(!is.na(data5_input_NC))
data5input <- data.frame(
  Mean_NC = data5_mean_NC,
  SD_NC = data5_sd_NC,
  N_NC = data5_n_NC
)

data5_input_PDR <- data5_input[, 11:32]
mean_PDR <- rowMeans(data5_input_PDR, na.rm = TRUE)
sd_PDR <- apply(data5_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR <- rowSums(!is.na(data5_input_PDR))
data5input_PDR <- data.frame(
  Mean_PDR = mean_PDR,
  SD_PDR = sd_PDR,
  N_PDR = n_PDR
)


data5_input <- cbind(data5input, data5input_PDR)
data5_input$reference <- "data5"

#data6
data6_input_NC <- data6_input[, 1:8]
data6_mean_NC <- rowMeans(data6_input_NC, na.rm = TRUE)
data6_sd_NC <- apply(data6_input_NC, 1, FUN = sd, na.rm = TRUE)
data6_n_NC <- rowSums(!is.na(data6_input_NC))
data6input <- data.frame(
  Mean_NC = data6_mean_NC,
  SD_NC = data6_sd_NC,
  N_NC = data6_n_NC
)

data6_input_PDR <- data6_input[, 9:16]
mean_PDR <- rowMeans(data6_input_PDR, na.rm = TRUE)
sd_PDR <- apply(data6_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR <- rowSums(!is.na(data6_input_PDR))
data6input_PDR <- data.frame(
  Mean_PDR = mean_PDR,
  SD_PDR = sd_PDR,
  N_PDR = n_PDR
)

data6_input <- cbind(data6input, data6input_PDR)
data6_input$reference <- "data6"

#data7
data7_input_NC <- data7_input[, 1:18]
data7_mean_NC <- rowMeans(data7_input_NC, na.rm = TRUE)
data7_sd_NC <- apply(data7_input_NC, 1, FUN = sd, na.rm = TRUE)
data7_n_NC <- rowSums(!is.na(data7_input_NC))
data7input <- data.frame(
  Mean_NC = data7_mean_NC,
  SD_NC = data7_sd_NC,
  N_NC = data7_n_NC
)

data7_input_PDR <- data7_input[, 19:38]
mean_PDR <- rowMeans(data7_input_PDR, na.rm = TRUE)
sd_PDR <- apply(data7_input_PDR, 1, FUN = sd, na.rm = TRUE)
n_PDR <- rowSums(!is.na(data7_input_PDR))
data7input_PDR <- data.frame(
  Mean_PDR = mean_PDR,
  SD_PDR = sd_PDR,
  N_PDR = n_PDR
)

# 合并NC组和PDR组的数据
data7_input <- cbind(data7input, data7input_PDR)
data7_input$reference <- "data7"


list_dfs=list(data2_input, data3_input, data4_input, data5_input, data6_input, data7_input)
list_names=c("data2", "data3", "data4", "data5", "data6", "data7")

compute_protein_stats <- function(dataframe_in, dataset_name_in) {
  # 设置输入数据框
  data_input <- dataframe_in
  
  # 获取列名
  my_cols <- colnames(data_input)
  
  # 删除列名中的数字以找到唯一组
  groups_in_cols <- unique(gsub('\\d+', '', my_cols))
  
  # 提取蛋白质列表
  protein_list <- row.names(data_input)
  
  # 初始化空数据框用于存储结果
  df_to_output <- data.frame()
  
  # 循环处理每个蛋白质
  for (protein_id in protein_list) {
    data_protein <- data_input[protein_id,]
    
    # 如果列名中包含'NC'，则进行处理
    if(any(grepl('NC', my_cols))) {
      cols_NC <- my_cols[grep('NC', my_cols)]
      NC_df <- data_protein[,cols_NC]
      
      # 计算NC的均值、数量和标准差
      Mean_NC <- mean(as.matrix(NC_df), na.rm = TRUE)
      N_NC <- length(cols_NC)
      SD_NC <- sd(NC_df, na.rm = TRUE)
      
      vec_NC <- c(Mean_NC, N_NC, SD_NC)
    } else {
      vec_NC <- c(NA, NA, NA)
    }
    
    
    
    # 如果列名中包含'PDR'，则进行处理
    if(any(grepl('PDR', my_cols))) {
      cols_PDR <- my_cols[grep('PDR', my_cols)]
      PDR_df <- data_protein[,cols_PDR]
      
      # 计算PDR的均值、数量和标准差
      Mean_PDR <- mean(as.matrix(PDR_df), na.rm = TRUE)
      N_PDR <- length(cols_PDR)
      SD_PDR <- sd(PDR_df, na.rm = TRUE)
      
      vec_PDR <- c(Mean_PDR, N_PDR, SD_PDR)
    } else {
      vec_PDR <- c(NA, NA, NA)
    }
    
    # 合并NC和PDR的统计数据，并设置蛋白质ID为行名
    out <- t(data.frame(c(vec_NC, vec_PDR)))
    row.names(out) <- protein_id
    df_to_output <- rbind(df_to_output, out)
  }
  
  # 设置输出数据框的列名
  colnames(df_to_output) <- c("Mean_NC", "N_NC", "SD_NC", "Mean_PDR", "N_PDR", "SD_PDR")
  
  # 将结果写入到文件中
  write.table(df_to_output, file = paste0("D:/MS/results/", dataset_name_in, ".txt"), row.names = TRUE, quote = FALSE, sep = "\t")
}

      


compute_protein_stats(list_dfs[[1]], list_names[1])
for(i in seq(1, length(list_dfs))){
  print(head(list_dfs[i]))
  print(list_names[i])
  compute_protein_stats(list_dfs[[i]], list_names[i])
}

#list_dfs[[1]]
#```

data2_input$protein = row.names(data2_input)
row.names(data2_input) = NULL

data3_input$protein = row.names(data3_input)
row.names(data3_input) = NULL

data4_input$protein = row.names(data4_input)
row.names(data4_input) = NULL

data5_input$protein = row.names(data5_input)
row.names(data5_input) = NULL

data6_input$protein = row.names(data6_input)
row.names(data6_input) = NULL

data7_input$protein = row.names(data7_input)
row.names(data7_input) = NULL


#Concatenate input dfs - final tables 
#```{r}
data_input = data.frame()
for( df in list(data2_input, data3_input, data4_input, data5_input, data6_input, data7_input
)){
  data_input = rbind(data_input, df)
}

write.table(data_input, file = "D:/MS/results/data_Concatenate.txt", row.names = TRUE, quote = FALSE, sep = "\t" )

## Random effects model - Run analysis
#estimator: DL
#sm: SMD 
#input tables are transferred via filezilla and located in the folder called data. 
###  metaanalysis(6 data)
#```{r}
metaanalysis = read.table("D:/MS/results/data_Concatenate.txt")

all_proteins=as.character(unique(metaanalysis$protein))
all_random_models=data.frame()
for(protein_id in all_proteins){
  temp = metaanalysis[metaanalysis$protein==protein_id,]
  rand_eff_mod=metacont(
    N_PDR, Mean_PDR, SD_PDR,
    N_NC, Mean_NC, SD_NC,
    data = temp,
    comb.fixed = FALSE,
    comb.random = TRUE,
    method.tau = "DL",
    hakn = FALSE,
    prediction = TRUE,
    sm = "SMD"
  )
  
  out=as.data.frame(rand_eff_mod[c("TE.random", "lower.random","upper.random","pval.random","seTE.random", "zval.random", "k", "Q", "tau2", "tau", "I2")])
  out$protein=protein_id
  
  all_random_models=rbind(all_random_models, out)
}

all_random_models$FDR=p.adjust(all_random_models$pval.random, method='BH')

write.table(all_random_models, file= "D:/MS/results/meta_results.txt", row.names = TRUE, quote = FALSE, sep = "\t" )
#```

# Meta-analysis results: Find significantly altered proteins (FDR < 10%) 
#Venn diagram: number of ALL distinct protein IDs and FDR<10%. 
#Find proteins with consistent alterations across datasets. From output df find proteins with FDR<10% and extract their FC from the meta-analysis input file. 
#Heatmap
#Enrichment analysis

#```{r}筛选FDR小于0.1
output = read.table("D:/MS/results/meta_results.txt")
input = read.table("D:/MS/results/data_Concatenate.txt")
FDR10 = output[abs(output$FDR) < 0.10, ]
# 筛选表达升高的蛋白 (TE.random > 0) 
increased_proteins = FDR10[FDR10$TE.random > 0, ]

# 筛选表达降低的蛋白 (TE.random < 0) 
decreased_proteins = FDR10[FDR10$TE.random < 0, ]


# 筛选 FDR10 排名前 50 的蛋白质
FDR10_sorted <- FDR10[order(FDR10$FDR), ]
top50_proteins <- head(FDR10_sorted, 50)

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
library(biomaRt)

# 使用 biomaRt 进行蛋白质 accession 转换
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 将 protein 列中的蛋白质 accession 转换为基因符号
protein_accessions <- top50_proteins$protein

# 从 UniProt ID 获取基因符号
protein_to_gene <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = protein_accessions,
  mart = mart
)

# 重命名列名
colnames(protein_to_gene) <- c("protein", "gene_symbol")

# 合并 gene_symbol 到 top50_proteins
top50_proteins <- merge(top50_proteins, protein_to_gene, by = "protein", all.x = TRUE)


# 过滤数据，只保留共同蛋白质的数据，并使用 gene_symbol 作为标注
shared_proteins <- top50_proteins$gene_symbol  # 从 top50_proteins 中提取基因符号列表

# 过滤数据，只保留共同蛋白质的数据
input_shared <- input[input$protein %in% top50_proteins$protein, c("reference", "protein", "Mean_NC", "Mean_PDR")]

# 添加基因符号
input_shared <- merge(input_shared, protein_to_gene, by = "protein", all.x = TRUE)

# 计算 FC (Fold Change)
input_shared$FC <- input_shared$Mean_PDR - input_shared$Mean_NC

# 手动替换 NA 的 gene_symbol 值
input_shared$gene_symbol[input_shared$protein == "B0YIW2"] <- "APOC3"
input_shared$gene_symbol[input_shared$protein == "C9J2H1"] <- "ITIH5"
input_shared$gene_symbol[input_shared$protein == "J3KNF6"] <- "RGMB"
input_shared$gene_symbol[input_shared$protein == "J3KQ18"] <- "DDT"
input_shared$gene_symbol[input_shared$protein == "O94985-2"] <- "Alcalpha1"
input_shared$gene_symbol[input_shared$protein == "P43251-2"] <- "BTD"
input_shared$gene_symbol[input_shared$protein == "Q13510-2"] <- "ASAH1"
input_shared$gene_symbol[input_shared$protein == "Q6P2Q0"] <- "ABCB9"



# 将数据转换为宽格式，行是基因符号，列是研究
input_shared_wide <- reshape2::dcast(input_shared, gene_symbol ~ reference, value.var="FC")

# 设置行名为gene symbol
input_shared_wide <- input_shared_wide[-50, ]
row.names(input_shared_wide) <- input_shared_wide$gene_symbol
input_shared_wide <- input_shared_wide[, -1]
# 替换 NA 为 0
input_shared_wide[is.na(input_shared_wide)] <- 0

# 保存结果
write.table(input_shared_wide, file="D:/MS/results/shared_proteins_FC_gene_symbols.txt", row.names=TRUE, quote=FALSE, sep="\t")
input_shared_wide = read.table("D:/MS/results/shared_proteins_FC_gene_symbols.txt")
proteins = rownames(input_shared_wide)
list(proteins)
# 使用biomaRt查询GO注释
library(biomaRt)

# 连接到Ensembl数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 查询基因的GO功能注释
go_data <- getBM(attributes = c('external_gene_name', 'go_id', 'name_1006'), 
                 filters = 'external_gene_name', 
                 values = proteins, 
                 mart = ensembl)

######### HEATMAP
library(RColorBrewer)
library(pheatmap)

paletteLength <- 100
col.pal <- colorRampPalette(c('blue3', 'white', 'red2'))(paletteLength)

myBreaks <- c(
  seq(-1, -0.05, length.out=50),
  seq(0.05, 1, length.out=50)
)

### Annotations for the heatmap
metadata <- data.frame(sample.size = c( "= 8", "= 8",  "= 19", "= 32", "= 16" , "= 38"))
row.names(metadata) <- colnames(input_shared_wide)

metadata$instrument <- c( "TripleTOF 5600", "Orbitrap Fusion",  "Orbitrap Fusion", "Orbitrap Fusion", "Q Exactive HF", "Orbitrap Fusion")

metadata$ocular.region <- c( "vitreous", "vitreous",  "vitreous", "vitreous", "tear", "vitreous")

metadata$labeling <- c( "label-free", "label-free", "label-free", "labeled", "label-free", "labeled")

my_color <- list(
  sample.size = c("= 8" = "#ccffcc",  "= 16" = "#99ff33", "= 19" = "#66ff66",  "= 27" = "#339900", "= 32" = "#336600", "= 38" = "#003300"), 
  ocular.region = c("vitreous" = "#ffcc99" , "tear" = "#cc6600"), 
  instrument = c("Orbitrap Fusion" = "#ffccff" , "Q Exactive HF" = "#cc33cc" , "TripleTOF 5600" = "#993399"), 
  labeling = c("labeled" = "#6666ff", "label-free" = "#ffff66")
)

pheatmap(input_shared_wide,
         color = col.pal,
         border_color = "grey95",
         breaks = myBreaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2",
         annotation_colors = my_color,
         annotation_col = metadata,
         fontsize_row = 7,
         fontsize_col = 7,
         cellwidth = 10, cellheight = 10,
         main = "Top 50 Proteins Heatmap (Gene Symbols)",
         filename = "D:/MS/results/Top50_proteins_heatmap_gene_symbols.pdf"
)





##enrichment analysis
#Prepare the input df. the meta-analysis results include protein-centric info, meaning isoforms. In addition different uniprot accesssion numbers can belong to the same gene. 
#columns needed: protein, TE.random, pval, FDR, k, lower.random, upper.random, Q 
#1. new column remove the isoform info "-XX"
#2. find duplicates and check whether the direction of change is same. 
#3. sort by pval and keep the most significant pval (remove the least significant one). 
#4. remove the IDs with inconsistent mean dif. 
#5. Uniprot ID conversion - get unique (primary) gene name 
#6. If exists, remove the obsolete entries. 


#data prep 
#```{r}
# 读取数据
meta_results <- read.table("D:/MS/results/meta_results.txt", header = TRUE)
meta_results$uniprot.input <- gsub("-.*", "", meta_results$protein)  # 移除 isoform 信息
meta_results <- meta_results[, c("protein", "uniprot.input","TE.random", "pval.random", "FDR", "k", "lower.random", "upper.random", "Q")]

# 查找重复项
data_dup <- meta_results[duplicated(meta_results$uniprot.input),]
data_dup <- meta_results[meta_results$uniprot.input %in% data_dup$uniprot.input, ]

total_dup <- c(data_dup$uniprot.input)
nodup <- meta_results[!meta_results$uniprot.input %in% total_dup, ]  # 移除所有重复项

dup_up <- data_dup[data_dup$TE.random > 0,]
dup_up <- dup_up[order(dup_up$pval.random),]
dup_up_unique <- dup_up[!duplicated(dup_up$uniprot.input, fromLast = FALSE),]

dup_dn <- data_dup[data_dup$TE.random < 0,]
dup_dn <- dup_dn[order(dup_dn$pval.random),]
dup_dn_unique <- dup_dn[!duplicated(dup_dn$uniprot.input, fromLast = FALSE),]

# 合并数据
unique <- data.frame()
for(df in list(dup_up_unique, dup_dn_unique)) {
  unique <- rbind(unique, df)
}

# 查找不一致的表达
unique <- unique[order(unique$uniprot.input),]
unique_incons <- unique[duplicated(unique$uniprot.input),]

# 移除不一致的 ID
unique_cleaned <- unique[!unique$uniprot.input %in% unique_incons$uniprot.input, ]

# 合并 cleaned 和 nodup
cleaned <- data.frame()
for(df in list(unique_cleaned, nodup)) {
  cleaned <- rbind(cleaned, df)
}

# 获取 Entrez ID
unique_uniprot_ids <- unique(cleaned$uniprot.input)
mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = unique_uniprot_ids,
                                 columns = c("ENTREZID"),
                                 keytype = "UNIPROT")

# 处理多对一映射
mapping <- mapping[!duplicated(mapping$UNIPROT), ]

# 确保数据一致性并合并
if (all(cleaned$uniprot.input %in% mapping$UNIPROT)) {
  cleaned <- left_join(cleaned, mapping, by = c("uniprot.input" = "UNIPROT"))
  print(head(cleaned))
} else {
  print("Error: Not all uniprot.input values found in mapping")
}

# 去除 NA 值
cleaned <- cleaned[!is.na(cleaned$ENTREZID) & !is.na(cleaned$pval.random) & !is.na(cleaned$TE.random), ]
# 去除 ENTREZID 中的重复项
cleaned <- cleaned[!duplicated(cleaned$ENTREZID), ]
write.xlsx(cleaned, file = "D:/MS/results/enrichment_prepare.xlsx", rowNames = FALSE)

# 创建函数来保存数据到 Excel
save_to_excel <- function(df, file_name) {
  write.xlsx(df, file = file_name, rowNames = FALSE)
}

# 进行 KEGG 富集分析
kegg_enrich <- enrichKEGG(
  gene = cleaned$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 检查 KEGG 结果是否为空
if (is.null(kegg_enrich) || nrow(as.data.frame(kegg_enrich)) == 0) {
  print("No KEGG enrichment results found.")
} else {
  print("KEGG enrichment results found.")
  kegg_enrich_df <- as.data.frame(kegg_enrich)
  
  # 计算 enrichment score
  kegg_enrich_df$enrichment_score <- as.numeric(sapply(strsplit(kegg_enrich_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # 添加被富集的基因
  kegg_enrich_df$genes <- sapply(kegg_enrich_df$geneID, function(x) paste(cleaned$ENTREZID[cleaned$ENTREZID %in% unlist(strsplit(x, "/"))], collapse = ", "))
  
  # 按 p.adjust 排序，展示前 10 个通路
  kegg_enrich_df <- kegg_enrich_df %>% arrange(p.adjust) %>% head(10)
  
  # 保存 KEGG 富集结果到 Excel
  save_to_excel(kegg_enrich_df, "D:/MS/results/kegg_enrichment_results_top10.xlsx")
}
  
  # 可视化 KEGG 富集结果：条形图
  ggplot(kegg_enrich_df, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = "Top 10 KEGG Pathway Enrichment Analysis", x = "Pathway Name", y = "Gene Count") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 20))  # 调整条形图中的字体大小
  
  ggsave("D:/MS/results/kegg_enrichment_barplot_top10.pdf", width = 10, height = 10)  # 调整保存的图像大小
  
  # 可视化 KEGG 富集结果：气泡图
  ggplot(kegg_enrich_df, aes(x = enrichment_score, y = reorder(Description, Count), size = Count, color = p.adjust)) +
    geom_point(alpha = 0.6) +
    scale_size_area(max_size = 10) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "Top 10 KEGG Pathway Enrichment Analysis", x = "Enrichment Score", y = "Pathway Name") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 20))  # 调整气泡图中的字体大小
  
  ggsave("D:/MS/results/kegg_enrichment_bubbleplot_top10.pdf", width = 10, height = 10)  # 调整保存的图像大小
}



# 进行 GO 富集分析（BP, CC, MF）
go_categories <- c("BP", "CC", "MF")
go_results <- data.frame()
for (category in go_categories) {
  go_enrich <- enrichGO(
    gene = cleaned$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = category,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
    go_enrich_df <- as.data.frame(go_enrich)
    go_enrich_df$Category <- category  # 添加类别信息
    
    # 添加被富集的基因
    go_enrich_df$genes <- sapply(go_enrich_df$geneID, function(x) paste(cleaned$ENTREZID[cleaned$ENTREZID %in% unlist(strsplit(x, "/"))], collapse = ", "))
    
    go_results <- rbind(go_results, go_enrich_df)
  }
}    


# 检查 GO 富集结果
if (nrow(go_results) == 0) {
  print("No GO enrichment results found.")
} else {
  # 计算 enrichment score
  go_results$enrichment_score <- as.numeric(sapply(strsplit(go_results$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  
  # 只显示前 10 个结果
  go_results_top10 <- go_results %>% group_by(Category) %>% top_n(10, wt = -p.adjust)
  
  # 保存 GO 富集结果到 Excel
  save_to_excel(go_results_top10, "D:/MS/results/go_enrichment_results.xlsx")
  
  # 创建函数保存每个类别的图
  plot_go_enrichment <- function(go_data, category, filename_prefix) {
    # 提取特定类别的数据
    go_data_filtered <- go_data %>% filter(Category == category)
    
    # 可视化 GO 富集结果：条形图
    ggplot(go_data_filtered, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = paste("GO Enrichment Analysis -", category), x = "GO Term", y = "Gene Count") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, face = "bold"),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text.x = element_text(size = 14),   
        axis.text.y = element_text(size = 14)   
      )
    
    ggsave(paste0("D:/MS/results/", filename_prefix, "_barplot.pdf"), width = 15, height = 10)  # 调整保存的图像大小
    
    # 可视化 GO 富集结果：气泡图
    ggplot(go_data_filtered, aes(x = enrichment_score, y = reorder(Description, Count), size = Count, color = p.adjust)) +
      geom_point(alpha = 0.6) +
      scale_size_area(max_size = 10) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = paste("GO Enrichment Analysis -", category), x = "Enrichment Score", y = "GO Term") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, face = "bold"),  # 调整标题字体大小
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14)   
      )
    
    ggsave(paste0("D:/MS/results/", filename_prefix, "_bubbleplot.pdf"), width = 15, height = 10)  # 调整保存的图像大小
  }
  
  # 依次生成 BP、CC、MF 的图
  plot_go_enrichment(go_results_top10, "BP", "go_enrichment_BP")
  plot_go_enrichment(go_results_top10, "CC", "go_enrichment_CC")
  plot_go_enrichment(go_results_top10, "MF", "go_enrichment_MF")
}

  
input_shared_wide <- read.table("D:/MS/results/shared_proteins_FC_gene_symbols.txt", 
                                header = TRUE, sep = "\t", row.names = 1)


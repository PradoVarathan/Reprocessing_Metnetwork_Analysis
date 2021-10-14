
whole <- read.csv('~/Documents/Test_Folder/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv')
#Randomize the columns

T_1 <- whole[,1:200]
T_2 <- whole[,201:ncol(whole)]
row.names(T_1) = T_1[,1]
row.names(T_2) = T_1[,1]
T_1[,1] = NULL
T_1[is.na(T_1)] = 0
T_2[is.na(T_2)] = 0


#MF
#install.packages("fastICA")
library(fastICA)
mf_1 <- fastICA(T_1, 100)
mf_2 <- fastICA(T_2,100)

Metagenes_1 <- mf_1$S
Metagenes_2 <- mf_2$S
MetaSamples_1 <- mf_1$A
MetaSamples_2 <- mf_2$A


cor_list <- function(){
  i=1
  res <- data.frame("Genes" = 1:8817,col2 = 0,col3 = 0)
  col2 = paste0("Index_",i)
  col3 = paste0("Correlation_",i)
  colnames(res) <- c("Genes",col2,col3)
  for (i in 1:nrow(res)){
    old_cor = 0
    idx = 0
    print(i)
    for (j in 1:nrow(res)){
      temp_cor = cor(Metagenes_1[i,], Metagenes_2[j,], method = "pearson")
      if (Mod(temp_cor) > old_cor){
        old_cor = temp_cor
        idx = j
      }
    }
    res$Index_1[i] = idx
    res$Correlation_1[i] = old_cor
  }
  
  i = 2
  res_2 <- data.frame("Genes" = 1:8817,col2 = 0,col3 = 0)
  col2 = paste0("Index_",i)
  col3 = paste0("Correlation_",i)
  colnames(res_2) <- c("Genes",col2,col3)
  for (i in 1:nrow(res_2)){
    old_cor = 0
    idx = 0
    print(i)
    for (j in 1:nrow(res_2)){
      temp_cor = cor(Metagenes_2[i,], Metagenes_1[j,], method = "pearson")
      if (Mod(temp_cor) > old_cor){
        old_cor = temp_cor
        idx = j
      }
    }
    res_2$Index_2[i] = idx
    res_2$Correlation_2[i] = old_cor
  }
  
  
}
res_2_mod = res_2[res$Index_1,]
res_final <- cbind(res,res_2_mod)

idx_rbh = which(res_final$Genes == res_final$Index_2)
res_rbh = res_final[idx_rbh,]

max_rdx = which(res_rbh$Correlation_1 == max(res_rbh$Correlation_2))
res_rbh[max_rdx,]

gene_names = whole$X
head(res_rbh)
gene_names[5]

result_rbh = res_rbh[,1:3]
colnames(result_rbh) = c("Genes_from_Set1", "Genes_from_Set2","Correlation")
for ( k in 1:nrow(result_rbh)){
  temp = as.integer(result_rbh$Genes_from_Set1[k])
  result_rbh$Genes_from_Set1[k] = gene_names[temp]
  temp = as.integer(result_rbh$Genes_from_Set2[k])
  result_rbh$Genes_from_Set2[k] = gene_names[temp]
}
library(dplyr)
result_rbh = result_rbh %>% arrange(desc(Correlation))
head(result_rbh)


install.packages("MCL")
library(MCL)

length(unique(result_rbh$Genes_from_Set1))
length(unique(result_rbh$Genes_from_Set2))
dim(result_rbh)

mcl(result_rbh)












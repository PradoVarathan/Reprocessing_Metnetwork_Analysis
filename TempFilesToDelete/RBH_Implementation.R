library(dplyr)
synapser::synLogin()
whole <- read.csv(synapser::synGet('syn21266454')$path, row.names = 1)
whole <- read.csv('~/Documents/Test_Folder/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv')

#Randomize the columnsa
set.seed(42)
whole_rand <- whole[ ,
                     sample(colnames(whole),
                            size = dim(whole)[2],
                            replace = FALSE
                           )
]

T_1 <- whole_rand[,
                  1:floor(ncol(whole_rand)/2)
                 ]
T_2 <- whole_rand[,
                  (floor(ncol(whole_rand)/2)+1):ncol(whole_rand)
                  ]
# @Pradeep Data is log2 so zero should essentially be the mean/median 
# Definitely what I would do, but the mean is not exactly zero but there are 
# issues re:
# mean(na.omit(as.numeric(T_1[2,])))   = 0.2278513
# median(na.omit(as.numeric(T_1[2,]))) = 0.2318029
# sd(na.omit(as.numeric(T_1[2,])))     = 0.1314889
# table(is.na(T_1[2,])) : FALSE = 102, True = 98
comp <- T_1[2,]
comp[is.na(comp)] = 0
# mean(as.numeric(comp))   = 0.1162042 - mean changes quite a bit
# median(as.numeric(comp)) = 0
# sd(as.numeric(comp))     = 0.1476958 - Variance changes very little
# If we want to set it to a functional
comp <- T_1[2,]
comp[is.na(comp)] = 0.2278513
# mean(as.numeric(comp))   = 0.2278513 - mean is the same
# median(as.numeric(comp)) = 0.2278513 
# sd(as.numeric(comp))     = 0.0936749 - Variance changes quite a bit 
# In terms of compression zero seems great, but it does change our  overall  distributions

T_1[is.na(T_1)] = 0
T_2[is.na(T_2)] = 0


#MF
#install.packages("fastICA")
library(fastICA)
mf_1 <- fastICA::fastICA(T_1, 100)
mf_2 <- fastICA::fastICA(T_2,100)

Metagenes_1 <- mf_1$S
Metagenes_2 <- mf_2$S
MetaSamples_1 <- mf_1$A
MetaSamples_2 <- mf_2$A


cor_list <- function(){
  i = 1
  
  #res <- data.frame("Genes" = 1:8817, col2 = 0, col3 = 0)
  res <- data.frame("Genes" = 1:nrow(Metagenes_2), col2 = 0, col3 = 0)
  col2 = paste0("Index_", i)
  col3 = paste0("Correlation_", i)
  colnames(res) <- c("Genes", col2, col3)
  # These inner Loop could a stand alone function implemented with parApply
  # Returns a list of each run, and convert to df with do.call(rbind,<list>)
  #' @param i A Vector to compare to each vector in j
  #' @param j A data frame of vectors to compare to i
  #' @param method  a correlation method ie "pearson", "kendall", or "spearman"
  cor_func <- function( i, j, method) {
    old_cor = 0
    idx = 0
    #print(i)
    res <- c(NA, NA)
    names(res) <- c('Index_1','Correlation_1')
    for (comp_row in 1:nrow(j)){
      temp_cor = cor(i, j[comp_row,], method = method)
      if (Mod(temp_cor) > old_cor){
        old_cor = temp_cor
        idx = comp_row
      }
    }
    #res$Index_1[i] = idx
    res$Index_1 = idx
    #res$Correlation_1[i] = old_cor
    res$Correlation_1 = old_cor
    
    return(res)
  }
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  res <- parallel::parApply(cl = cl,
                     Metagenes_1,
                     1,
                     cor_func,
                     j = Metagenes_2,
                     method = "pearson")
  # res
  res <- as.data.frame(do.call(rbind,res))
  
  # Run inverse:
  res2 <- parallel::parApply(cl = cl,
                            Metagenes_2,
                            1,
                            cor_func,
                            j = Metagenes_1,
                            method = "pearson")
  # res
  res2 <- as.data.frame(do.call(rbind,res2))
  ######  ######  ######  ######  ######  ######  ######  ######  ######
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

cor_list()

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

result_rbh = result_rbh %>% arrange(desc(Correlation))
head(result_rbh)


install.packages("MCL")
library(MCL)

length(unique(result_rbh$Genes_from_Set1))
length(unique(result_rbh$Genes_from_Set2))
dim(result_rbh)

mcl(result_rbh)












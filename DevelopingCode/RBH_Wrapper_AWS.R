

library(dplyr)

if(!("fastICA" %in% rownames(installed.packages()))){
  install.packages("ica")
  install.packages("varhandle")
}
library(ica)
library(varhandle)




# Function to find the correlation matrix between two datasets via --------

cor_func <- function( i, j, method, name_ind = 1) {
  old_cor = 0
  idx = 0
  p_vals <- NULL
  cors <- NULL
  res <- c(NA, NA, NA ,NA)
  
  for (comp_row in 1:nrow(j)){
    temp_cor = cor(i, j[comp_row,], method = method)
    cors <- c(cors,temp_cor)
    
    temp_cortest = cor.test(i, j[comp_row,], method = method)
    p_vals <- c(p_vals, temp_cortest$p.value)
    #_# Need to take Mod of both vals to ensure absolute value is taken 
    if (Mod(temp_cor) > Mod(old_cor) ){
      old_cor = temp_cor
      idx = comp_row
    }
  }
  
  adj_p <- p.adjust(p_vals, method = 'fdr' )
  
  if( length(cors[ which(adj_p < 0.05) ]) > 0 ){
    keep_ind <- which(abs(cors) == max(abs(cors[ which(adj_p < 0.05) ])))
    keep <- cors[keep_ind]
  }else{
    keep <- NA
    keep_ind <- NA
  }
  
  res[1] = idx
  res[2] = old_cor
  res[3] = keep_ind
  res[4] = keep
  names(res) <- c(paste0('Index_',name_ind),
                  paste0('Correlation_',name_ind),
                  paste0('Index_adj_',name_ind),
                  paste0('Correlation_adj_',name_ind)
  )
  return(res)
}


# To pairwise process RBH pair cor matrix and gene list  ------------------
# Returns a dataframe of RBH results and their corelation
#' @param dataset1 A data frame of expression profile 
#' @param dataset2 A data frame of expression profile to compare against dataset1 with same genes as dataset1
#' @param method  a correlation method ie "pearson", "kendall", or "spearman"
#' @param geneList a list of gene name from daataset1 and dataset2 in same order 
#' @param cl cluster object for parallel processing 

rbh_pair_processing = function(dataset1, dataset2, method, geneList, cl){
  
  #Spliting the dataset into components using ICA analysis
  mf_1 <- ica::icafast(dataset1, 100)
  mf_2 <- ica::icafast(dataset2,100)
  
  Metagenes_1 <- mf_1$S
  Metagenes_2 <- mf_2$S
  MetaSamples_1 <- mf_1$A
  MetaSamples_2 <- mf_2$A
  print("Starting first corelation matrix build")
  res_1 <-as.data.frame(t(as.matrix(parallel::parApply(
    cl = cl, Metagenes_1, 1, cor_func,
    j = Metagenes_2, method = "pearson", name_ind = 1
  ))))
  print("Starting second corelation matrix build")
  
  res_2 <- as.data.frame(t(as.matrix( parallel::parApply(
    cl = cl, Metagenes_2, 1, cor_func,
    j = Metagenes_1, method = "pearson", name_ind = 2
  )))) 
  
  res_final <- cbind(res_1,res_2)
  
  idx_rbh = which(res_final$Index_1 == res_final$Index_2)
  res_rbh = res_final[idx_rbh,]
  
  print("Starting assesing gene names")
  
  result_rbh = res_rbh[,c('Index_1','Index_2','Correlation_1')]
  colnames(result_rbh) = c("Genes_from_Set1", "Genes_from_Set2","Correlation")
  for ( k in 1:nrow(result_rbh)){
    temp = as.integer(result_rbh$Genes_from_Set1[k])
    result_rbh$Genes_from_Set1[k] = geneList[temp]
    temp = as.integer(result_rbh$Genes_from_Set2[k])
    result_rbh$Genes_from_Set2[k] = geneList[temp]
  }
  
  return(result_rbh)
}

winsorizeData <- function(x){
  library(dplyr)
  
  winsorize <- function(x,per=.99){
    up <- quantile(x,per,na.rm=T)
    low <- quantile(x,1-per,na.rm=T)
    x[x>=up] <- up
    x[x<=low] <- low
    return(x)
  }
  
  
  replaceNaMean <- function(x){
    if(sum(is.na(x))>0){
      y <- x
      y[is.na(x)] <- mean(x,na.rm=T)
      return(y)
    }else{
      return(x)
    }
  }
  
  x <- t(x)
  x <- apply(x,2,winsorize)
  x <- apply(x,2,replaceNaMean)
  x <- scale(x)
  return(x)
}



# TESTING SECTION ---------------------------------------------------------
synapser::synLogin("Pradeep_Varathan_98","pRADEEP@01")
install.packages('reader')
library(reader)

data_1 <- reader::reader(synapser::synGet('syn26387667')$path)
data_2 <- reader::reader(synapser::synGet('syn26340641')$path)

# devtools::install_github('PradoVarathan/metanetwork')
# BiocManager::install('WGCNA')
# library(metanetwork)
data_1 = winsorizeData(data_1)
data_2 = winsorizeData(data_2)
data_1 = t(data_1)
data_2 = t(data_2)

geneList1 = rownames(data_1) 
geneList2 = rownames(data_2)

geneList = intersect(geneList1,geneList2)
data_1 = data_1[geneList,]
data_2= data_2[geneList,]


#Starting cluster 
  
cl <- parallel::makeCluster(parallel::detectCores()-1)

res_1_2 = rbh_pair_processing(data_1, data_2, 'pearson', geneList = geneList, cl)

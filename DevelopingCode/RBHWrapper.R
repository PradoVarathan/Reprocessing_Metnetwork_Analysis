

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
  mf_1 <- fastICA::fastICA(dataset1, 10)
  mf_2 <- fastICA::fastICA(dataset2,10)
  
  Metagenes_1 <- mf_1$S
  Metagenes_2 <- mf_2$S
  MetaSamples_1 <- mf_1$A
  MetaSamples_2 <- mf_2$A
  
  res_1 <-as.data.frame(t(as.matrix(parallel::parApply(
    cl = cl, Metagenes_1, 1, cor_func,
    j = Metagenes_2, method = "pearson", name_ind = 1
  ))))
  res_2 <- as.data.frame(t(as.matrix( parallel::parApply(
    cl = cl, Metagenes_2, 1, cor_func,
    j = Metagenes_1, method = "pearson", name_ind = 2
  )))) 
  
  res_final <- cbind(res_1,res_2)
  
  idx_rbh = which(res_final$Index_1 == res_final$Index_2)
  res_rbh = res_final[idx_rbh,]
  
  
  result_rbh = res_rbh[,1:3]
  colnames(result_rbh) = c("Genes_from_Set1", "Genes_from_Set2","Correlation")
  for ( k in 1:nrow(result_rbh)){
    temp = as.integer(result_rbh$Genes_from_Set1[k])
    result_rbh$Genes_from_Set1[k] = geneList[temp]
    temp = as.integer(result_rbh$Genes_from_Set2[k])
    result_rbh$Genes_from_Set2[k] = geneList[temp]
  }
  
  return(result_rbh)
}





# TESTING SECTION ---------------------------------------------------------

#Obtaining the geneList and removing it from the dataframe
geneList = whole$X
whole$X = NULL


#Randomize the columns
set.seed(42)
whole_rand <- whole[ ,
                     sample(colnames(whole),
                            size = dim(whole)[2],
                            replace = FALSE
                     )
]

# Dividing the dataset into six testing parts

T_1 <- whole_rand[,
                  1:floor(ncol(whole_rand)/6)
]
T_2 <- whole_rand[,
                  (floor(ncol(whole_rand)/6)+1):(floor(ncol(whole_rand)/6)*2)
]
T_3 <- whole_rand[,
                  ((floor(ncol(whole_rand)/6)*2)+1):(floor(ncol(whole_rand)/6)*3)
]
T_4 <- whole_rand[,
                  ((floor(ncol(whole_rand)/6)*3)+1):(floor(ncol(whole_rand)/6)*4)
]
T_5 <- whole_rand[,
                  ((floor(ncol(whole_rand)/6)*4)+1):(floor(ncol(whole_rand)/6)*5)
]
T_6 <- whole_rand[,
                  ((floor(ncol(whole_rand)/6)*5)+1):ncol(whole_rand)
]


# Dividing into metagenes and metasamples via ICA
T_1[is.na(T_1)] = 0
T_2[is.na(T_2)] = 0
T_3[is.na(T_3)] = 0
T_4[is.na(T_4)] = 0
T_5[is.na(T_5)] = 0
T_6[is.na(T_6)] = 0
#MF
mf_1 <- fastICA::fastICA(T_1, 10)
mf_2 <- fastICA::fastICA(T_2,10)
mf_3 <- fastICA::fastICA(T_3,10)
mf_4 <- fastICA::fastICA(T_4,10)
mf_5 <- fastICA::fastICA(T_5,10)
mf_6 <- fastICA::fastICA(T_6,10)

Metagenes_1 <- mf_1$S
Metagenes_2 <- mf_2$S
Metagenes_3 <- mf_3$S
Metagenes_4 <- mf_4$S
Metagenes_5 <- mf_5$S
Metagenes_6 <- mf_6$S



#Starting cluster 
cl <- parallel::makeCluster(parallel::detectCores()-1)

res_1_2 = rbh_pair_processing(Metagenes_1, Metagenes_2, 'pearson', geneList = geneList, cl)

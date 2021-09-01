library(pracma)
library(Matrix)
for (i in 1:size(partitions,2)){
  # Dev
  # print(i)
  for (j in 1:size(partitions,2)){
    
    if (j<i){#metric is symetric, so save time by computing only half
      adjustedrand_pairwise[i,j] <- adjustedrand(partitionID[i,],partitionID[j,])  
      
      # Dev
      jj <- paste0(str(i),'_',str(j))
      print(jj)
    }
  }
}
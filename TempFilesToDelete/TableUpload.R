
library(dplyr, quietly = TRUE)
library(glmnet, quietly = TRUE)
library(randomForest, quietly = TRUE)
library(Hmisc, quietly = TRUE)
library(lars, quietly = TRUE)
library(WGCNA, quietly = TRUE)
library(synapser, quietly = TRUE)
library(metanetwork, quietly = TRUE)
library(githubr, quietly = TRUE)
library(c3net, quietly = TRUE)
library(config, quietly = TRUE)
library(optparse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parmigene, quietly = TRUE)
library(WGCNA, quietly = TRUE)
library(reader, quietly = TRUE)
library(tidyr)
library(plyr)
library(igraph)
library(parallel)
library(doParallel)
library(foreach)
library(linkcomm)
library(R.matlab) #to read mat file
#Obtaining the required data
data = readMat('~/Downloads/SpeakeasyRes.mat')
A = data$partition.codes[[1]][[1]]
head(A)
GeneNames = colnames(ADJ)

gene_Modules <- as.data.frame(A)
colnames(gene_Modules) <- c('Gene.ID','moduleNumber')
if(!(is.null(GeneNames))){
  gene_Modules['Gene.ID'] = GeneNames
}
gene_Modules['moduleSize'] <- 0
count_mtx <- table(gene_Modules$moduleNumber)

for (i in 1:nrow(gene_Modules)){
  temp = gene_Modules$moduleNumber[i]
  gene_Modules$moduleSize[i] = as.integer(count_mtx[temp])
}

gene_Modules = gene_Modules %>%
  group_by(Gene.ID) %>%
  dplyr::mutate(moduleNumber = factor(moduleNumber),
                moduleNumber = as.numeric(moduleNumber)) %>% filter(moduleSize > 30)


# Change cluster number to color labels
gene_Modules$moduleNumber = as.numeric(factor(gene_Modules$moduleNumber))
gene_Modules$moduleLabel = WGCNA::labels2colors(gene_Modules$moduleNumber)
mod = gene_Modules[c('Gene.ID','moduleNumber','moduleLabel')]
head(mod)

mod['algorithms'] = 'SpeakEasy (M)'

synStore(Table('syn26276910', mod))

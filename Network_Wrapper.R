
# Calling in Libraries ----------------------------------------------------

library(dplyr)
library(glmnet)
library(randomForest)
library(Hmisc)
library(lars)
library(WGCNA)
library(synapser)
library(metanetwork)
library(githubr)
library(c3net)
library(config)
library(optparse)

# Obtaining the data - From User --------------------------------------------

option_list <- list(make_option(c("-u","--synapse_user"), type="character", action = "store",
                               help = "Synapse User name"),
                   make_option(c("-p","--synapse_pass"), type="character", action = "store",
                               help = "Synapse User Password"),
                   make_option(c("-p","--config_file"), type="character", action = "store",
                               help = "Path to the complete config file"))          
req_args <- parse_args(OptionParser(option_list=option_list))

# Obtaining the data - From Synapse --------------------------------------------

#Setting up the cofig file 
Sys.setwnd(R_CONFIG_ACTIVE = "default")
config <- config::get(file = option_list$config_file)

#Linking with Project
synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
project = Project(config$project_id)
project <- synStore(project)

# Data
synID_input = config$input_synid
data = synGet(synID_input, downloadLocation = config$temp_storage_loc)

# Performing the analysis -------------------------------------------------

net_methods = config$netowrk_method
net_methods = strsplit(net_methods,",")[[1]]

'''
    Working on only the light networks for now - C3 net, mrnet and WGCNA
    Just storing the function for reference to delete later

    c3netWrapper  = function(data,path=NULL,pval=1,outputpath){
  library(c3net)
  network <- c3net(t(data))
  #save(network,file=paste0(outputpath,'result_c3net.rda'))
  network <- network*upper.tri(network)
  write.csv(network,file=paste0(outputpath,'c3netNetwork.csv'),quote=F)
}

    mrnetWrapper = function(data,path=NULL,pval=1,outputpath){
  library(parmigene)                               ########## Need installation in docker
  metanetwork::aracne(data,path,pval,outputpath)
  library(data.table)
  library(dplyr)
  if(pval==1){
    fileName <- paste0(outputpath,'aracneNetwork.csv')
  }else{
    fileName <- paste0(outputpath,'aracneThresholdNetwork.csv')
  }
  cat('fileName:',fileName,'\n')
  #load(fileName)
  #data.matrix(data.frame(data.table::fread('~/Desktop/sparrowZNetwork.csv',data.table=F),row.names=1))
  network <- data.table::fread(fileName,data.table=F) %>%
    data.frame(row.names=1) %>%
    data.matrix
  network <- network+t(network)
  gc()
  network <- parmigene::mrnet(data.matrix(network))
  #save(network,file=paste0(outputpath,'result_mrnet.rda'))
  network <- network*upper.tri(network)
  write.csv(network,file=paste0(outputpath,'mrnetNetwork.csv'),quote=F)
}
'''

for (method in net_methods){
    switch(method,
            "c3net" = )
}



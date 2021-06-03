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
library(data.table)
library(parmigene)
library(WGCNA)
library(Rmpi)
library(utilityFunctions)

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
Sys.setenv(R_CONFIG_ACTIVE = "default")
config <- config::get(file = option_list$config_file)

#Linking with Project
synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
project = Project(config$input_profile$project_id)
project <- synStore(project)

buildConsensus = function(outputpath,networkFolderId, fileName){

  #get all networks from Synapse
  foo <- synGet(paste0('select name,id from file where parentId==\'',networkFolderId,'\''))
  foo <- dplyr::filter(foo,file.name!='bicNetworks.rda' & file.name!='rankConsensusNetwork.csv')
  bar <- lapply(foo$file.id,synGet,downloadLocation=outputpath)

  loadNetwork <- function(file){
    sparrowNetwork <- data.table::fread(file,stringsAsFactors=FALSE,data.table=F)
    rownames(sparrowNetwork) <- sparrowNetwork$V1
    sparrowNetwork <- sparrowNetwork[,-1]
    gc()
    return(sparrowNetwork)
  }

  networkFiles <- sapply(bar,function(x){return(x@filePath)})
  networks <- lapply(networkFiles,loadNetwork)
  networks <- lapply(networks,data.matrix)

  networks$rankConsensus <- metanetwork::rankConsensus(networks)
  cat('built rank consensus\n')
  cat('write rank consensus\n')
  write.csv(networks$rankConsensus,file=paste0(outputpath,'rankConsensusNetwork.csv'),quote=F)
  if(!is.null(fileName)){
    library(Matrix)
    networkMethods <- sapply(bar,synGetAnnotation,which='method')
    cat('grabbed methods\n')
    #build rank consensus
    cat('updated methods\n')
    networkMethods <- c(networkMethods,'rankConsensus')
    cat('reading in data\n')
    options(stringsAsFactors = F)
    dataSet <- reader(fileName, row.names=1)
    cat('turning data into data matrix\n')
    dataSet <- data.matrix(dataSet)
    cat('build bicNetworks\n')
    #bicNetworks <- lapply(networks,metanetwork::computeBICcurve,dataSet,maxEdges=1e5)
    bicNetworks <- metanetwork::computeBICcurve(networks$rankConsensus,dataSet,maxEdges=2e5)
    #cat('make names of bicNetworks\n')
    #names(bicNetworks) <- 'rankConsensus'
    cat('save bicNetworks\n')
    save(bicNetworks,file=paste0(outputpath,'bicNetworks.rda'))
  }

}
buildConsensus(outputpath,networkFolderId,fileName)

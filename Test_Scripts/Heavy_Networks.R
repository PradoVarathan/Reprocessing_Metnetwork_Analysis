# Calling in Libraries ----------------------------------------------------

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
library(reader, quietly = TRUE)
library(Rmpi)
library(bench)
library(parallel)
library(doParallel)
#library(utilityFunctions) -->installation error

# Obtaining the data - From User --------------------------------------------

option_list <- list(make_option(c("-u","--synapse_user"), type="character", action = "store",
                                help = "Synapse User name"),
                    make_option(c("-p","--synapse_pass"), type="character", action = "store",
                                help = "Synapse User Password"),
                    make_option(c("-c","--config_file"), type="character", action = "store",
                                help = "Path to the complete config file"),
                    make_option(c("-s","--percentage_data"), type="numeric", action = "store",
                                help = "Section/Percentage of the data to be used"))          
req_args <- parse_args(OptionParser(option_list=option_list))

# Obtaining the data - From Synapse --------------------------------------------

#Setting up the cofig file 
Sys.setenv(R_CONFIG_ACTIVE = "default")
config <- config::get(file = req_args$config_file)
setwd(config$input_profile$temp_storage_loc)
  
#Linking with Project
synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
project = Project(config$input_profile$project_id)
project <- synStore(project)
  
# Data
synID_input = config$input_profile$input_synid
data = synGet(synID_input, downloadLocation = config$input_profile$temp_storage_loc)

dataFolder <- Folder('Heavy',parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)

# Performing the analysis -------------------------------------------------
  
net_methods = c('genie3','tigress')
data = reader::reader(data$path)
rows_to_use = (nrow(data)*req_args$percentage_data)/100
cols_to_use = (ncol(data)*req_args$percentage_data)/100
data = data[1:rows_to_use,1:cols_to_use]
nslaves = config$computing_specs$medium_ncores
mpi.spawn.Rslaves(nslaves=nslaves,hosts=NULL);


  
for (method in net_methods){
    registerDoParallel(cl)

  benchMarkRes = mark(switch(method,
         "genie3" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                               outputpath = config$output_profile$output_path),
         "tigress" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                                outputpath = config$output_profile$output_path))
  
,check=FALSE)
  setwd(config$input_profile$temp_storage_loc)

  benchmark_filename = paste0(method,"_",as.character(req_args$percentage_data),"_Performance.rds")
  saveRDS(benchMarkRes,benchmark_filename)

  file <- File(path = benchmark_filename, parent = dataFolder)
  file <- synStore(file)
}

 #mpi.bcast.cmd(q("no"));
  mpi.close.Rslaves()
  mpi.quit(save = "no")
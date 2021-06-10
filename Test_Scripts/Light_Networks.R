
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

#Iterating over every test config file in folder 
Sys.setenv(R_CONFIG_ACTIVE = "default")
config_files = list.files(path = req_args$config_file)

for (file_config in config_files){
  
  #Setting up the cofig file 
  file_config = paste0(req_args$config_file,"/",file_config)
  config <- config::get(file = file_config)
  setwd(config$input_profile$temp_storage_loc)
  
  #Linking with Project
  synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
  project = Project(config$input_profile$project_id)
  project <- synStore(project)
  
  # Data
  synID_input = config$input_profile$input_synid
  data = synGet(synID_input, downloadLocation = config$input_profile$temp_storage_loc)
  
  # Performing the analysis -------------------------------------------------
  
  net_methods = config$input_profile$net_methods
  data = reader::reader(data$path)
  rows_to_use = nrow(data)*req_args$percentage_data
  cols_to_use = ncol(data)*req_args$percentage_data
  data = data[1:rows_to_use,1:cols_to_use]
  
  if(is.null(config$input_profile$na_fill)){
    data <- metanetwork::winsorizeData(data)
  }else{
    data[is.na(data)] <- config$input_profile$na_fill 
  }
  
  for (method in net_methods){
    
    benchMarkRes = bench(switch(method,
           "c3net" = c3netWrapper(data, path = NULL, pval = config$input_profile$p_val_c3net, config$output_profile$output_path),# What does this path define in main function?
           "mrnet" = mrnetWrapper(data, path = NULL, pval = config$input_profile$p_val_mrnet, config$output_profile$output_path),
           "wgcna" = wgcnaTOM(data, path = NULL, pval = config$input_profile$p_val_wgcna, config$output_profile$output_path, 
                              config$input_profile$rsquaredcut, config$input_profile$defaultnaPower),
           "lassoAIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
           "lassoBIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
           "lassoCV1se" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                     outputpath = config$output_profile$output_path),
           "lassoCVmin" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                     outputpath = config$output_profile$output_path),
           "ridgeAIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
           "ridgeBIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
           "ridgeCV1se" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                     outputpath = config$output_profile$output_path),
           "ridgeCVmin" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                     outputpath = config$output_profile$output_path),
           "sparrowZ" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
           "sparrow2Z" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                    outputpath = config$output_profile$output_path),
           "genie3" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
           "tigress" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                                  outputpath = config$output_profile$output_path))
    
  ,check=FALSE)
    benchmark_filename = paste0(config$output_profile$output_path,"/",method,"_Performance_",as.character(reqs_args$percentage_data),".csv")
    write.csv(benchMarkRes,benchmark_filename,quote=F,row.names = T)
  }
}


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
library(parallel)
library(doParallel)
#library(utilityFunctions) -->installation error

# Obtaining the data - From User --------------------------------------------

option_list <- list(make_option(c("-u","--synapse_user"), type="character", action = "store",
                                help = "Synapse User name"),
                    make_option(c("-p","--synapse_pass"), type="character", action = "store",
                                help = "Synapse User Password"),
                    make_option(c("-c","--config_file"), type="character", action = "store",
                                help = "Path to the complete config file"))          
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

# Registering the parallel clusters
if(!is.na(config$computing_specs$medium_ncores)){
  nslaves = config$computing_specs$medium_ncores
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=NULL);
}
if(!is.na(config$computing_specs$heavy_ncores)){
  nslaves = config$computing_specs$heavy_ncores
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=NULL);
}

# Performing the analysis -------------------------------------------------

net_methods = config$input_profile$network_method
data = reader::reader(data$path)

if(is.null(config$input_profile$na_fill)){
  data <- metanetwork::winsorizeData(data)
}else{
  data[is.na(data)] <- config$input_profile$na_fill # Ideally, user could insert a large negative number or use min(data)
}

for (method in net_methods){# Assuming we have more methods - not developing for now

  switch(method,
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
  
  output_filename <- list.files(pattern = method )
  
  
}
  mpi.close.Rslaves()
  mpi.quit(save = "no")

# Obtaining the data - For provenance --------------------------------------------

activity <- synapser::synGet(config$input_profile$project_id)

dataFolder <- Folder(method,parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)
for (filePath in output_filename){
  file <- File(path = filePath, parent = dataFolder)
  file <- synStore(file)
}

all.annotations <- list(
  dataType = config$provenance$annotations$data_type,
  resourceType = config$provenance$annotations$resuorce_type,
  metadataType = config$provenance$annotations$metadata_type,
  isModelSystem = config$provenance$annotations$ismodelsystem,
  isMultiSpecimen = config$provenance$annotations$ismultispecimen,
  fileFormat = config$provenance$annotations$fileformat,
  grant = config$provenance$annotations$grant,
  species = config$provenance$annotations$species,
  organ = config$provenance$annotations$organ,
  tissue = config$provenance$annotations$tissue,
  study = config$provenance$annotations$study, 
  consortium = config$provenance$annotations$consortium,
  assay = config$provenance$annotations$assay
)

thisRepo = NULL
thisFile = NULL

try(
  thisRepo <- githubr::getRepo(
    repository = config$provenance$code_annotations$repository,
    ref = config$provenance$code_annotations$ref,
    refName = config$provenance$code_annotations$ref_name
  ), silent = TRUE
)
try(
  thisFile <- githubr::getPermlink(
    repository = thisRepo,
    repositoryPath = config$provenance$code_annotations$repository_path
  ), silent = TRUE
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path = output_filename[length(output_filename)],
  name = config$output_profile$output_name,
  parentId = activity$properties$id),
  used = config$input_profile$input_synid,
  activityName = config$provenance$activity_name,
  executed = thisFile,
  activityDescription = config$provenance$activity_description
)

synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)

# Formatting the network to md5 format --------------------------------------------

md5Command <- paste0('md5sum ', output_filename[length(output_filename)])
md5 <- strsplit(system(md5Command, intern = TRUE), '  ')[[1]][1]
cat(md5, '\n', file = config$output_profile$md5_output_path, sep = '')
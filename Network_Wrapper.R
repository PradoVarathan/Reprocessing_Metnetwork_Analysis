
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

# Data
synID_input = config$input_profile$input_synid
data = synGet(synID_input, downloadLocation = config$input_profile$temp_storage_loc)

# Performing the analysis -------------------------------------------------

net_methods = config$input_profile$network_method
data = read.csv(data$path)

for (method in net_methods){# Assuming we have more methods - not developing for now
  switch(method,
         "c3net" = c3netWrapper(data, path = NULL, pval = config$input_profile$p_val_c3net, config$output_profile$output_path),# What does this path define in main function?
         "mrnet" = mrnetWrapper(data, path = NULL, pval = config$input_profile$p_val_mrnet, config$output_profile$output_path),
         "wgcna" = wgcnaTOM(data, path = NULL, pval = config$input_profile$p_val_wgcna, config$output_profile$output_path, 
                            config$input_profile$rsquaredcut, config$input_profile$defaultnaPower))
  output_filename <- paste0(config$output_profile$output_path,method,'Network.csv')
  
  
}

# Obtaining the data - For provenance --------------------------------------------

activity <- synapser::synGet(config$input_profile$project_id)

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
  path = output_filename,
  name = config$output_profile$output_name,
  parentId = activity$properties$id),
  used = config$input_profile$input_synid,
  activityName = config$provenance$activity_name,
  executed = thisFile,
  activityDescription = config$provenance$activity_description
)

synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)


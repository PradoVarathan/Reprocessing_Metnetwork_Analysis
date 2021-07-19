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
library(WGCNA, quietly = TRUE)

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

#Linking with Project
synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
input_file = synGet(config$input_profile$input_proj_id,downloadLocation = temp_input_loc)
fileName = input_file[[1]]$path
project = Project(config$input_profile$project_id)
project <- synStore(project)


outputpath = config$input_profile$temp_storage_loc
networkFolderId = config$input_profile$input_folderid
pattern_id = config$input_profile$input_proj_id
buildConsensus(outputpath = outputpath,networkFolderId = networkFolderId,pattern_id = pattern_id, fileName = fileName)


# Obtaining the data - For provenance --------------------------------------------

activity <- synapser::synGet(config$input_profile$project_id)

dataFolder <- Folder(method,parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)

file <- File(path = paste0(outputpath,'rankConsensusNetwork.csv'), parent = dataFolder)
file <- synStore(file)

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

if((!is.na(config$computing_specs$heavy_ncores)) || (!is.na(config$computing_specs$medium_ncores))){
  mpi.quit(save = "no")
}
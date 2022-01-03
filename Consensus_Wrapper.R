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
library(reader, quietly = TRUE)
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
input_file = synGet(config$input_profile$input_proj_id,downloadLocation = config$input_profile$temp_input_loc)
fileName = input_file$path
project = Project(config$input_profile$project_id)
project <- synStore(project)


outputpath = config$input_profile$temp_storage_loc
networkFolderId = config$input_profile$input_folderid
pattern_id = config$input_profile$pattern_id

buildConsensus(outputpath = outputpath,networkFolderId = networkFolderId,pattern_id = pattern_id, fileName = fileName)


# Obtaining the data - For provenance --------------------------------------------

activity <- synapser::synGet(config$input_profile$project_id)

dataFolder <- Folder('Consensus',parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)

file <- File(path = paste0(outputpath,'rankConsensusNetwork.csv'), parent = dataFolder)
file <- synStore(file)
file2 <- File(path = paste0(outputpath,'bicNetworks.rda'), parent = dataFolder)
file2 <- synStore(file2)
all.annotations <- synGetAnnotations(config$input_profile$input_synid)

checkAnnotations <- function(annotations, config){
  annot_default <- list(
    dataType = NULL,
    resourceType = NULL,
    metadataType = NULL,
    isModelSystem = NULL,
    isMultiSpecimen = NULL,
    fileFormat = NULL,
    grant = NULL,
    species = NULL,
    organ = NULL,
    tissue = NULL,
    study = NULL, 
    consortium = NULL,
    assay = NULL
  )
  for (item in names(annot_default)){
    if (!is.null(config$provenance$annotations$item)){
      annot_default$item = config$provenance$annotations$item[[1]]
    }
    else if (!is.null(annotations$item)){
      annot_default$item = annotations$item[[1]]
    }
  }
}

all.annotations <- checkAnnotations(all.annotations,config)

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
  path = paste0(outputpath,'rankConsensusNetwork.csv'),
  name = 'RankConsensusNetwork',
  parentId = activity$properties$id),
  used = config$input_profile$input_synid,
  activityName = config$provenance$activity_name,
  executed = thisFile,
  activityDescription = config$provenance$activity_description
)

synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)

# Formatting the network to md5 format --------------------------------------------

md5Command <- paste0('md5sum ', paste0(outputpath,'rankConsensusNetwork.csv'))
md5 <- strsplit(system(md5Command, intern = TRUE), '  ')[[1]][1]
cat(md5, '\n', file = config$output_profile$md5_output_path, sep = '')

if((!is.na(config$computing_specs$heavy_ncores)) || (!is.na(config$computing_specs$medium_ncores))){
  mpi.quit(save = "no")
}

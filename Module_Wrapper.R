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
library(tidyr)
library(plyr)
library(igraph)
library(synapseClient)
library(parallel)
library(doParallel)
library(foreach)
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
bic_file = synGet(config$input_profile$bic_file, downloadLocation = config$input_profile$temp_input_loc)
fileName = input_file$path
project = Project(config$input_profile$project_id)
project <- synStore(project)

#Creating parallel cores
nc = detectCores()
if (nc > 2){
  cl = makeCluster(nc - 2)
} else {
  cl = makeCluster(1)
}
registerDoParallel(cl)

#### Get input data from synapse and formulate adjacency matrix ####
# Get bicNetworks.rda
bic.obj = bic_file
load(bic.obj$path) # this will load an R object nameds bicNetworks
all.used.ids = config$input_profile$bic_file # for provenance
writeLines(paste('Total number of edges', sum(as.matrix(bicNetworks$network))))

# Get rankconsensus network for weights
rank.cons = data.table::fread(fileName, data.table = F, header = T)
rownames(rank.cons) = rank.cons$V1
rank.cons$V1 = NULL
all.used.ids = c(all.used.ids, rankConsNet.id)

# Formulate adjacency matrix
adj = rank.cons
adj[!as.matrix(bicNetworks$network)] = 0
adj = data.matrix(adj)

rm(list = c('bicNetworks', 'rank.cons'))
gc()

#### Compute modules using specified algorithm ####
# Get a specific algorithm

# Running CF Finder - Under comments since it is not working wiht lisencing 
#cf_loc = synGet('syn7806853',downloadLocation = '/home/sage')
#temp_command <- paste0("unzip ",cf_loc$path," -d /home/sage/")
#system(temp_command)
#CFinder = metanetwork::findModules.CFinder(adj, '/home/sage/CFinder-2.0.6--1448/', nperm = 3, min.module.size = 30)

# GANXIS
ga_loc = synGet('syn7806859',downloadLocation = '/home/sage')
temp_command <- paste0("unzip ",ga_loc$path," -d /home/sage/")
system(temp_command)
GANXiS = metanetwork::findModules.GANXiS(adj, '/home/sage/GANXiS_v3.0.2/', nperm = 3, min.module.size = 30)

#Fast Greedy Algorithm
fast_greedy = metanetwork::findModules.fast_greedy(adj, nperm = 3, min.module.size = 30)

#Label_Prop
label_prop = metanetwork::findModules.label_prop(adj, nperm = 3, min.module.size = 30)

#Louvain
louvain = metanetwork::findModules.louvain(adj, nperm = 3, min.module.size = 30)

#Walktrap
walktrap = metanetwork::findModules.walktrap(adj, nperm = 3, min.module.size = 30)

#Infomap
infomap = metanetwork::findModules.infomap(adj, nperm = 3, min.module.size = 30)

#Link Communities
linkcommunities = metanetwork::findModules.linkcommunities(adj, nperm = 3, min.module.size = 30)

#Spinglass
spinglass = metanetwork::findModules.spinglass(adj, nperm = 3, min.module.size = 30)

#megena
megena = metanetwork::findModules.megena(adj, method = "pearson", FDR.cutoff = 0.05, module.pval = 0.05, hub.pval = 0.05,doPar = TRUE)

#Speakeasy
#speakeasy_r = metanetwork::findModules.speakeasy(adj)

algorithms = c('GANXis','fast_greedy','label_prop','louvain','walktrap','infomap','linkcommunities','spinglass','megena')
# algorithms = c('CFinder', algorithms, 'speakeasy') add after correction of CF issue and rechecking the speakeasy implementation in R

# Obtaining the data - For provenance --------------------------------------------

activity <- synapser::synGet(config$input_profile$project_id)

dataFolder <- Folder(method,parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)
for (filePath in algorithms){
  temp_out <- paste0(outputpath,filePath,'.csv')
  write.csv(filePath,temp_out)
  file <- File(path = temp_out, parent = dataFolder)
  file <- synStore(file)
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
  path = temp_out,
  name = algorithm,
  parentId = activity$properties$id),
  used = config$input_profile$input_synid,
  activityName = config$provenance$activity_name,
  executed = thisFile,
  activityDescription = config$provenance$activity_description
)

synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)

md5Command <- paste0('md5sum ', temp_out)
md5 <- strsplit(system(md5Command, intern = TRUE), '  ')[[1]][1]
cat(md5, '\n', file = config$output_profile$md5_output_path, sep = '')


}


# Formatting the network to md5 format --------------------------------------------


if((!is.na(config$computing_specs$heavy_ncores)) || (!is.na(config$computing_specs$medium_ncores))){
  mpi.quit(save = "no")
}
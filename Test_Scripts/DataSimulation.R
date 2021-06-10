#Calling in libraries
library(metanetwork)
library(synapser)

#Simulation parameters
n = 1000 # number of subjects
p = 20000 #number of genes
prop = 0.5 #probability
path = "/home/sage/SimulatedData_1000x20000.csv"

#Preparing data
data = simulateNetworkData(n,p,prop)
dataMatrix = data$data
write.table(dataMatrix, path,sep = ",",row.names = T,col.names = T,quote = F)

#Saving data to Project

project_id = 'syn25872026' # Change to your project ID
activity <- synapser::synGet(project_id)

all.annotations <- list(
  dataType = 'simulation data',
  resourceType = 'metadata',
  metadataType = 'analytical',
  isModelSystem = TRUE,
  isMultiSpecimen = FALSE,
  fileFormat = 'csv',
  grant = '',
  species = '',
  organ = '',
  tissue = '',
  study = '', 
  consortium = '',
  assay = '')

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
  path = path,
  name = 'Simulation Data',
  parentId = project_id),
  used = '',
  activityName = 'Data Simulation for metanetwork',
  executed = thisFile,
  activityDescription = 'Data simulation with 1000 subjects and 20000 genes with 0.5 probability'
)

synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)

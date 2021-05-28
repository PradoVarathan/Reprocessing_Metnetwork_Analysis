
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

for (method in net_methods){
    
}

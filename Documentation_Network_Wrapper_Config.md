   
## Preparing Config File for Network Creation 

For the current development, we use a single (default) config profile per file
 It is divided into four sections :
 
  1. Input Profile
  2. Provenanace
  3. Output Profile
  4. Computing Specs
  
 Each of these sections are further divided into subsections and tags that has to be added and required for most cases, else filled as NULL

***
#### `Input Profile`   



  
**input_synid**     
*Required*. The Synapse ID of the input file is to provided; User has to make sure that their synapse profile has access to this particular profile
    
**project_id**    
*Required*. The Synapse ID of the Project folder currently the user owns to build provenance and upload files 
    
**na_fill**    
*Required*. Value to be filled in case of missing values, else will be omitted from original dataset.  


**network_method**    
*Required*. The type of network methods to be used in a list format with quotes. 

  - The options in light network include c3net, mrnet, wgcna
  - The options in light network include lassoAIC, lassoBIC, lassoCV1se, lassoCVmin, ridgeAIC, ridgeBIC, ridgeCV1se, ridgeCVmin, sparrowZ, sparrow2Z
  - The options in light network include genie3 and tigress
     
**p_val_c3net**  
Optional. P - value threshold for the C3Net method, if specified in network_method, else leave blank or NULL
    
**p_val_mrnet**   
Optional. P - value threshold for the mrnet method, if specified in network_method, else leave blank or NULL
    
**p_val_wgcna**   
Optional. P - value threshold for the WGCNA method, if specified in network_method, else leave blank or NULL
    
**temp_storage_loc**: *Required*.Temporary file storage location for the input file which will be deleted after completion of the process
    
**rsquaredCut**   
Optional. R squared cutt-off for the WGCNA network method and if not used, leave blank or NULL
    
**defaultnaPower**   
Optional. Default NA power value for WGCNA network method and if not used, leave blank or NULL
    
***

#### `Provenance`      
    
**Annotations**   
Provenance annoation subsections - can be obtained from the parent syn ID   
More details about obtaining annotations from synapse can be found [here](https://help.synapse.org/docs/Annotations-and-Queries.2011070649.html)  

More details on the functioning of Provenance can be found [here](https://help.synapse.org/docs/Provenance.1972470373.html)  

A sample annotation profile includes the following required information and can be replaced with NULL if unknown


data_type: ['clinical','geneExpression']   
resource_type: metadata   
metadata_type: 'analytical covariates' 
ismodelsystem: FALSE   
ismultispecimen: TRUE  
fileformat: csv  
grant: U01AG046152  
species: Human  
organ: brain  
tissue: ['dorsolateral prefrontal cortex', 'Head of caudate nucleus', 'posterior cingulate cortex']  
study: ["ROSMAP","rnaSeqReprocessing"]  
consortium: AMP-AD  
assay: rnaSeq  
      
`code_annotations`  For development purpose only. **DO NOT CHANGE** the code annotations unless user determines to use another analysis pipeline

  - repository: PradoVarathan/Reprocessing_Metanetwork_Analysis # The current repository from where the analysis code is obtained from
  - ref: branch  
  - ref_name: main  
  - repository_path: Network_Wrapper.R # The code block path in the repository used for analysis   

**activity_name**    
*Required*. The activity name for Synapse Provenance path creation  

**activity_description**  
*Required*. Detailed explanation of the activity in quotes to be updated on Synapse  


***
#### `Output Profile`  

**output_path**   
*Required*. The full path for the final output file  

**output_name**  
*Required*. The output name as to be updated in Synapse Project  

**md5_output_path**  
*Required*. The output path for the md5 formatted output file  

**error_path**   
*Required*. The full path for the error file from the process  

***    
#### `Computing Specs`

**light_ncores**  
Total number of cores available for the light network processing (c3net, mrnet and wgcna) and if not being used, leave blank or NULL   


**medium_ncores**   
Total number of cores available for the medium network processing (lassoAIC, lassoBIC, lassoCV1se, lassoCVmin, ridgeAIC, ridgeBIC, ridgeCV1se, ridgeCVmin, sparrowZand sparrow2Z) and if not being used, leave blank or NULL   


**heavy_ncores**   
Total number of cores available for the heavy network processing (genie3 and tigress) and if not being used, leave blank or NULL   

    

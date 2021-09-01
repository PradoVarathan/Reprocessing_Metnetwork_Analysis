# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#calls SpeakEasy number of times=="iter" and the runs consensus clustering
#reason for multiple runs is to obtain representative clustering not driven
#by outliers and to facilitate overlapping output (if desired)

#also call "genIC" which sets up the initial conditions for all runs

#Description of inputs: (these are typically set in Layer_SpeakEasy
#"iter" is number of replicate runs - depends on how finely you want to estimate cluster confidence, but 1 is min value and 10+ are useful especially when overlapping output is requested
#"ADJ" is an adjacency matrix (weighted, unweighted, sparse or full, negative weights are fine too)
#"timesteps" is number of times to repeat the label selection - 30 tends to handle even huge matrices, as the number of labels collapses rapidly in real networks after the first few timesteps
#"multi_community_switch" should be in the range of [0 1].
#"nback" is hardcoded in layer_SpeakEasy because practical evidence suggests it never needs to be changed.  It is the number of previous time-steps overwhich labels are counted up to select the updated label.

bootstrap_SpeakEasy <- function(iter,ADJ,timesteps,nback,is_ADJ_weighted, force_efficient, main_iter, multi_community_switch,varargin){     

    cat('starting to setup ICs for all runs',main_iter)
    
    
    genIC_full <- function(ADJ,how_many_runs,nback,varargin){
        
        IC_store <- matrix(0,how_many_runs*(nback+1),max(dim(ADJ)))
        
        for (m in 1:max(dim(ADJ))){#setup initial listener histories with randomly selected neighbor labels
            
            contacts=as.list(which(ADJ[,m]==0))
            contacts = unlist(unname(contacts))
            if (length(contacts)==0){#particularly in very small modules, we might have no connections... usually diag(ADJ) is set to all ones, which also solves this, but in case it is not, this avoids an error message
                IC_store[,m] <- m
                
                
            } else {
                
                contacts_values <- ADJ[contacts,m]
                
                ## Made an edit here to change the contacts_values  to 0 rahter than 1 ##
                
                if (length(which(contacts_values==1))!=length(contacts_values)){#this should generally be the case for weighted networks
                    contacts_weights <- contacts_values-min(contacts_values)
                    # Commented the below sections since it was popping up NaN values -- Pradeep
                    contacts_weights <- contacts_weights/max(contacts_weights)#rescale weights for to [0 1]... not totally necessary, but seems to help in some cases
                    #should note that this is just for the selection of IC's - the
                    #links with negative weights in ADJ should continue to have
                    #different (not just rescaled) behavior from positive links -
                    #practically speaking tests bear this out as nodes with negative
                    #correlations to some other group of nodes end up outside of that
                    #given cluster
                    
                } else {
                    contacts_weights <- matrix(1,length(contacts),1)
                    
                }
                
                bins <- c(0,cumsum(contacts_weights/sum(contacts_weights)),1)
                bins <- as.matrix(bins)
                bins <- as.matrix(apply(bins, 1, FUN = min))
                
                library(pracma)
                bins <- sort(bins, decreasing = FALSE)
                ret_hist <- histc(rand(how_many_runs*(nback+1),1),bins)
                idxs <- c()
                for (i in 1:length(ret_hist$bin)){
                    temp <- ret_hist$bin[i]
                    idxs <- c(idxs,contacts[temp])
                }
                IC_store[,m] <- idxs
                #using cellfun (@randsample or bsxfn isn't any faster
                
            }
            
        }
        return(IC_store)
    }
    
    
    genIC_sparse_unweighted <- function(ADJ,how_many_runs,nback,varargin){
        library(pracma)
        k <- colSums(ADJ)
        kcumsum <- cumsum(k)
        basicrand <- rand((nback+1)*how_many_runs,max(dim(ADJ)))
        indices <- which(ADJ != 0, arr.ind = TRUE)
        row_ind <- indices[,1]
        basicrand_scaled <- ceil(basicrand*drop(repmat(k,size(basicrand,1),1)))
        basicrand_scaled_lifted <- basicrand_scaled+repmat(c(0,kcumsum[-length(kcumsum)]),size(basicrand_scaled,1),1)
        IC_store <- Reshape(row(basicrand_scaled_lifted),(nback+1)*how_many_runs,max(dim(ADJ)))
        return(IC_store)
    }
    
    if (is_ADJ_weighted!=0){#use the appropriate routine to generate the initial conditions for the run (we can do this more quickly for sparse unweighted networks
        IC_store <- genIC_full(ADJ,iter,nback)
    } else {
        IC_store <- genIC_sparse_unweighted(ADJ,iter,nback)
    }
    
    partitionID_lst <- matrix(nrow=max(dim(IC_store)),ncol=iter)
    for (i in 1:iter){
    
        cat('iteration #',str(i), ' of ',str(iter),main_iter)
    
        
        IC_store_relevant <- IC_store[1:nback+1,]
        IC_store <- IC_store[1,] #set in SSLPA
        speakeasy_res <- SpeakEasycore(ADJ,timesteps,IC_store_relevant,nback,force_efficient)
        listener_history = speakeasy_res$listener_history
        partitionID_lst[,i] = speakeasy_res$partitionID  #i.e.check sparse and suppress graphics
    }
    
    partition_columns <- c()
    for( i in ncol(partitionID_lst)){
        partition_columns <- c(partition_columns,as.vector(partitionID_lst[,i]))
    }
    
    ### Take a look into with @Jake
    # save_results <- 0
    # if (main_iter==1){
    #     if (save_results==1){
    #         disp('saving partitions')
    #         save record_of_all_partitions.mat ADJ partition_columns partitionID multi_community_switch#sometimes useful for very large networks if you want to save results before consensus clustering, especially if you want to try several multi-community cutoffs
    #     }
    # }
    
    cat('started consensus clustering', as.character(main_iter))
    
    virt_res <- virtual_cooccurrence(ADJ,partition_columns,partitionID, main_iter, multi_community_switch)
    
    [partition_codes partition_codes_overlapping cell_partition cell_partition_overlapping ]

}


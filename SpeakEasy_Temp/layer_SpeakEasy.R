# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




##This is the suggested fuction through which to interact with SpeakEasy
#clustering.  The loops in this funtion are related to calling SpeakEasy on
#subsets of data.  This maybe useful when there are genuine hierarchies of clusters in the
#data, although you cannot always force SpeakEasy to split clusters (unlike hierarchical clustering).
#
#Description of inputs:
#"layers" is number of times to subcluster: should be integer >=1 (1 is no subclustering)  If you want to avoid subclustering clusters smaller than a given size, enter that as the 2nd entry i.e. [3,15] does three rounds of clustering and won't touch a cluster unless it has at least 15 members
#"iter" is number of replicate runs - depends on how finely you want to estimate cluster confidence, but 1 is min value and 10+ are useful especially when overlapping output is requested
#"ADJ" is an double, sparse or full adjacency matrix (weighted, unweighted; negative weights are fine too)
#"timesteps" is number of times to repeat the label selection - 30 tends to handle even huge matrices, as the number of labels collapses rapidly in real networks after the first few timesteps
#varargin - optional input is "multi_community_switch" in the range of [0 1], which if <1 allows overlapping clusters - closer to 0 is more lenient definition of overlapping nodes.
#We suggest a value of 1/(max# of overlapping communities per node desired).
#For instance if you want nodes to be in no more than 3 clusters, multi_community_switch=.33333
#
#Description of outputs:
#"partition_codes_overlapping" is a two column matrix with node ID's in the first column and numeric label ID's in the 2nd
#"cell_partition_overlapping" is a sequence of cells, each containing a list of nodes in a given cluster.  The clusters are ordered in size.
#convenient_node ordering{i} is imply vertcat(partition_n_idx_overlapping_cell{i}{:})
#even though the output variables contains the word "overlapping", if disjoint output is requested, that's what these will contain
#
#to check that your ADJ was clustered reasonably, try running: imagesc(ADJ(convenient_node_ordering{1},convenient_node_ordering{1})

layer_SpeakEasy <- function(layers,iter,ADJ,timesteps,varargin=NULL){ ##ok<NCOMMA>
## Get ADJ properties to optimize runtime
# values in this section do not change the SpeakEasy method, just help to
# run it efficiently and determine how it is applied.  You likely don't need to
# adjust anything here, ever.

    nback <- 5#min recommended value is 5 - if you increase communities evolve more slowly... never seen a need to change this
    force_efficient <- 0#recommended to leave==0, which adaptively selects different forms of computation based on ADJ density, setting to 1 forces slower alternative
    max_ADJ_size <- 50000
    
    #set "max_ADJ_size" to the length of the largest ADJ your machine can hold in available
    #RAM 2x.  You can of course cluster sparse matrices that are larger than
    #this.  The reason for this value is that we use a different (~60# faster)
    #way of estimated expected label frequncy for relatively dense matrices if their dimension is less than this value.  The values computed are identical in either case; it's only a matter of efficiency.
    #On my laptop with 16GB, I can use a 10Kx10K full ADJ, so I'd set this value to 10000.
    #for my 128gb machine, I set "max_ADJ_size" to x with no trouble.
    #if you machine has 8GB of ram, I recommend max_ADJ_size=8000, as a comfortable value if you're not running a bunch of other stuff on the machine
    
    #deal with optional inputs for overlapping clustering and
    if (is.null(varargin)){
        multi_community_switch <- 1#default is disjoint clusters
    } else {
        if (length(varargin)>=1){
            multi_community_switch <- varargin[1]
            if (is.null(multi_community_switch)){
                multi_community_switch <- 1
            }
    
            if (multi_community_switch>1){
                cat('error - multi_community_switch should be in range of [0 1]')
                stop()
            }
    
            cat('optional input multi_community_switch set to ',str(multi_community_switch))
            if (multi_community_switch==1){
                cat('indicating disjoint output')
            }
        }
    }
    
    
    if (length(layers)==2){
        subclustersize <- layers[2]
        cat('optional input subclustersize set to ', str(subclustersize))
    } else {
        subclustersize <- 5#default...clusters smaller than this will be left alone (NOT subclustered)
    }
    
    noinputs=which(colSums(abs(ADJ))==0)#because if a node has no inputs it can't receive labels, we set it to receive it's own label, i.e. remain isolated
    ADJ[(noinputs-1)*(max(dim(ADJ)))+noinputs] <- 1
    
    ADJ_characteristics <- function(ADJ,max_ADJ_size, force_efficient){
        
        if (dim(ADJ)[1] != dim(ADJ)[2]){
            cat('your ADJ is not square, please fix')
            stop()
        }
        
        if ((max(ADJ))>1 | min(ADJ)<(-1)){
            cat('your connection strength is outside [-1 1] please fix')
            stop()
        }
        
        is_ADJ_weighted <- c()
        fraction_to_be_full <- 0.3#should be around .2 or less - if the matrix is more dense than this it will be treated as full
        
        sample_of_links = ADJ
        ADJ_density=length(which(sample_of_links!=0))/length(sample_of_links)
        cat('approximate edge density is ',str(ADJ_density))
        
        if (length(which(sample_of_links!=0))/length(sample_of_links)>fraction_to_be_full){
            ADJ_is_full <- 1
            
            if (length(which(sample_of_links>.9999))+length(which(abs(sample_of_links)<.0001))!=length(sample_of_links)){#in case there are rounding errors
                cat('weighted full')#ADJ with links of only +1 and -1 will still be considered weighted
                is_ADJ_weighted <- 1
            } else {
                cat('unweighted full')
                is_ADJ_weighted <- 0
            }
            
        } else {
            ADJ_is_full <- 0
            
            #in this case we can afford to test all links
            c <- ADJ[which(ADJ != 0)]
            if (length(which(c>.9999))!=length(c)){#in case there are rounding errors
                cat('weighted sparse')
                is_ADJ_weighted <- 1
            } else {
                cat('unweighted_sparse')
                is_ADJ_weighted <- 0
            }
            
        }
        
        
        #force_efficient--1 if the ADJ is huge or if it is small&dense
        if (force_efficient == 1){
            
            if (dim(ADJ)[1]<max_ADJ_size){#if matrix is small enough that we could use full multiplicationn
                
                if (ADJ_density>0.03){             #do not change this value; this value was selected by observing the two multiplication strategies it switches between require equal time at .03 density){
                    
                    force_efficient = 0}#ADJ dense enough that full mult is faster than pairwise (and small enough it will fit in memory)
                else {
                    
                    force_efficient = 1#'ultrasparse and small ADJ, so pairwise multi is faster'
                }
                
            } else {
                if (ADJ_density<=.03){#you still might want ADJ sparse if it's huge and RAM usage is an issue, but hoping the user necessarily does that
                    
                    force_efficient = 1}
                
                else {
                    force_efficient = 0}#if you have a huge matrix is still may be wise to put it as sparse even if density >.03.  .03 is the point at which it becomes faster to compute with a full, but you may need to convert to sparse simply for memory reasons
                
            }
            
        }
        values <- c(is_ADJ_weighted, force_efficient)
        return(values)
    }
    
    values <- ADJ_characteristics(ADJ, max_ADJ_size, force_efficient)#this is used because there are different routines for setting up the initial conditions for weighted vs unweighted networks
    is_ADJ_weighted <- values[1] 
    force_efficient <- values[2]
    
    ##  Call SpeakEasy on primary clusters and again on each cluster (if layers>1)
    cluster_stats_store <- list()
    for (i in 1:layers){
        main_iter <- i#for clarity when passing
        if (i==1){
            cat('calling main routine at level 1')
            [partition_codes{i} partition_codes_overlapping{i} cell_partition{i} cell_partition_overlapping{i}] in bootstrap_SpeakEasy(iter,ADJ,timesteps,nback,is_ADJ_weighted, force_efficient, main_iter,multi_community_switch)){
    
        } else {
            disp(['doing subclustering at level ' num2str(i)])
            printmessage('terminal output reduced',main_iter-1)
    
            for (j in 1:length(cell_partition_overlapping{i-1})){
                current_nodes <- cell_partition_overlapping{i-1}{j}#needed becase first node in sub cluster may be 99th overal etc
                if (length(current_nodes)>subclustersize){#don't even try to subcluster communites that are smaller than a certain size
    
                    [partition_codes_temp{j,1} partition_codes_overlapping_temp{j,1} cell_partition{j,1} cell_partition_overlapping_temp{j,1}] in bootstrap_SpeakEasy(iter,ADJ(current_nodes,current_nodes),timesteps,nback,is_ADJ_weighted,force_efficient,main_iter,multi_community_switch)){
    
                    partition_codes_overlapping_temp{j,1}(:,1) <- current_nodes(partition_codes_overlapping_temp{j,1}(:,1))
                    for (k in 1:length(cell_partition_overlapping_temp{j,1})){#update to main node ID's
                        cell_partition_overlapping_temp{j,1}{k} <- current_nodes( cell_partition_overlapping_temp{j,1}{k})
                    }
                } else {#if module is too small for subclustering
                    cell_partition_overlapping_temp{j,1} <- {current_nodes}
                }
            }
            cell_partition_overlapping{i} <- vertcat(cell_partition_overlapping_temp{:})
            cell_partition_overlapping_temp <- []
    
    
            for (m in 1:length(cell_partition_overlapping{i})){
            temp{m,1} <- [cell_partition_overlapping{i}{m}, repmat(m,length(cell_partition_overlapping{i}{m}),1)]
            }
            partition_codes_overlapping{i} <- vertcat(temp{:})
    
        }
        convenient_node_ordering{i} <- vertcat(cell_partition_overlapping{i}{:})
        return(partition_codes_overlapping,cell_partition_overlapping,convenient_node_ordering)
    }
}

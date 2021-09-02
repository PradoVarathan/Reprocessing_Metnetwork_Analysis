# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#If you only want disjoint clusters, you can spare yourself the pain of
#this function.  For disjoint output, SpeakEasy will only the first section to select one of several
#cluster partitions as the representative partition, based on the ARI betwen all pairs of clusters
#
#Outputs:
#"nodes_and_partition_identifiers_hard" two column output with first column as nodeID (1:length(ADJ)) and 2nd column as numeric cluster ID
#"nodes_and_partition_identifiers_overlapping"
#"cell_partition"
#"cell_partition_overlapping
# is list of nodeID's for nodes that show up in more than one community
#
#Inputs:
#"partitions" has columns which are alternate clustering results
#"partitionID" are alternate clustering results with clusters stored in cells, which make adjusted rand index faster to compute
#"accept_multi" higher values are more stringent criterion for overlapping clusters, 1==disjoint clusters

#If you desire overlapping output AND want to understand this function...
#The goal of this function is to do something that is really easy to do if
#we had access to a traditional co-occurrence matrix, but which is rather
#complicated when we cannot create one.  It is often impractical to create
#co-occurrence matrices when clustering 20K+ nodes, particularly if there
#are large clusters, because it will require huge amounts of ram - far more
#than required to store the ADJ, especially if the original ADJ is sparse.
#
#What we want to do for each node is to sum up the number of co-occurrences
#with the members of all other clusters, using the optimal cluster partition.  This is
#equivalent to summing the rows the the co-ocurrence matrix that correspond
#to the members of a given cluster (then doing this for all clusters and all nodes).
#
#The reason we want to do this is because it's how we assign
#multi-community nodes: they are nodes that co-occur with more than one
#community at a frequency greater than the input "accept_multi"

#Time saving idea that allows you to operate on just the set of "winning"
#clusters (i.e. clusters in the representative partition):
#if you look at all the partitions that have ever
#occured for a given set of nodes and you notice that at some other time
#these nodes have the same label, then that contributes to their co-occurence with that other label.
#you then need to count the occurences of all other labels for all nodes in
#a given winning cluster, which is what "fraction_counts_labels" does


function [nodes_and_partition_identifiers_hard nodes_and_partition_identifiers_overlapping cell_partition cell_partition_overlapping multicom_nodes_all] <- virtual_cooccurrence(ADJ,partitions,partitionID, main_iter, accept_multi)){
    library(pracma)
    
    partitions <- t(as.matrix(partitions))
    
    adjustedrand <- function(partitionA,partitionB){
        
        ng1 <- max(partitionA)
        ng2 <- max(partitionB)
        #ctabmat=full(sparse(partitionA,partitionB,1,ng1,ng2)); slightly faster, but memory intense for large matrices
        ctabmat <- sparseMatrix(i=as.vector(partitionA),j=as.vector(t(partitionB)),x=1L,
                                symmetric = FALSE, repr='T',dims = c(ng1,ng2))#same format as lxn
        #SpeakEasy passes in sequentially labeled clusters, so this shouldn't be an issue
        
        n <- sum(colSums(ctabmat))
        nis <- sum(rowSums(ctabmat)^2)#sum of squares of sums of rows
        njs <- sum(colSums(ctabmat)^2)#sum of squares of sums of columns
        
        t1 <- combs(n,2)#total number of pairs of entities
        t2 <- sum(colSums(ctabmat^2))#sum over rows & columnns of nij^2
        t3 <- 0.5*(nis+njs)
        
        #Expected index (for adjustment)
        nc <- (n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1))
        
        A <- (t1+t2-t3)#no. agreements
        D <-   (-t2+t3)#no. disagreements
        
        if (is.na(nc)){
            ARI <- 0
        }
        else if (t1==nc){
            ARI <- 0#avoid division by zero; if k=1, define Rand = 0
        } else {
            ARI <- (A-nc)/(t1-nc)#adjusted Rand - Hubert & Arabie 1985
        }
        return(ARI)
    }
    
    fast_ind2sub <- function(matsize,idx){
        nrs <- nrow(matsize)
        ncs <- ncol(matsize)
        r = (idx-1%nrows)+1
        c = ((idx-r)/nrows) + 1
        res <- c(r,c)
        return(res)
        
    }
    
    for (i in 2:size(partitions,2)){           #make sure partition codes are distinct across partitions
        partitions[,i] <- partitions[,i]+ max(partitions[,i-1])
    }
    
    adjustedrand_pairwise <- zeros(size(partitions,2))#holds all possible adjusted rand index comparisons of partitions
    for (i in 1:size(partitions,2)){
        # Dev
        # print(i)
        for (j in 1:size(partitions,2)){
    
            if (j<i){#metric is symetric, so save time by computing only half
                adjustedrand_pairwise[i,j] <- adjustedrand(partitionID[i,],partitionID[j,])  
            
                # Dev
                # print(j)
                }
        }
    }
    
    adjustedrand_pairwise <- adjustedrand_pairwise+ t(adjustedrand_pairwise)
    most_similar_var <-  max(colSums(adjustedrand_pairwise))
    most_similar_idx <- which.max(colSums(adjustedrand_pairwise))#select representative partition
    winning_partition <- partitions[,most_similar_idx]
    winning_members_unq <- unique(partitions[,most_similar_idx])
    
    
    cell_partition <- list()#get the indices of nodes which end up in the same cluster
    for (i in 1:length(winning_members_unq)){
        cell_partition[[i]]= c(which(winning_partition==winning_members_unq[i]))#seems like same contest as partitionID(most_similar_idx) but reordered
    }
    
    if (accept_multi<1){
        cat('getting overlapping clusters from virtual cooccurrence matrix',main_iter)
        all_partitions_unq <- unique(partitions)
    
        #markertitle is a list of clusterID's from all partitions
        #markerfound contains the positions of all the clusterID's in the ADJ
        #also, if clusterID==3, all the positions of nodes taked with clusterID==3 are in markerfound{3}
        sorted_labels <- sortrows(cbind(partitions,t(as.matrix(1:numel(partitions)))))
        sorted_labels <- as.matrix(sorted_labels)
        transitions=matrix(c(0, which(sorted_labels[2:nrow(sorted_labels),1]-sorted_labels[1:(nrow(sorted_labels)-1),1]!=0),size(sorted_labels,1)))#find sets of locations of the end of the previous label
        markertitle <- zeros(length(all_partitions_unq),1)#for all clusters, get their locations in the partitions and their numeric identifier
        markerfound = list()
        for (i in 1:(size(transitions,1)-1)){#puts the ith label in the ith cell of markerfound - remember that labels are unique across partitions
            markerfound[[i]] <- sorted_labels[(transitions[i]+1:transitions[i+1]),1]#this actually gets very big too
            markertitle[i,1] <- sorted_labels[(transitions[i]+1),1]
        }
    
        nodes_and_partition_identifiers_overlapping <- list()
        multi_store <- list()
        mutlicom_node_cell <- list() 
        cell_partition_overlapping <- list() 
        ban_module_size <- 3
        
        banned_modules <- which(lengths(cell_partition) < ban_module_size)#sometimes a node may be half in it's "own cluster" (size one) and in some other, so not really multi-comm
        for (i in 1:length(winning_members_unq)){#for each cluster in winning partition
    
            if (!(i %in% (c(0, banned_modules)))){
                winning_partition_idxs=which(winning_partition==winning_members_unq[i],arr.ind = T)
                winning_partition_chunk <- partitions[winning_partition_idxs,]#all labels ever used for a cluster in all other partitions
                all_labels_for_winning_cluster_unq <- unique(winning_partition_chunk)
                countseach <- pracma::histc(as.matrix(winning_partition_chunk),all_labels_for_winning_cluster_unq)
                countseach <- countseach$cnt
                
                fraction_counts_labels <- list()
                for (j in 1:length(all_labels_for_winning_cluster_unq)){#for all labels that have ever characterized any node in a winning cluster
    
                    linear_idx = markerfound[[all_labels_for_winning_cluster_unq[j]]]
                    fraction_counts_labels[[j]] <- (as.matrix(countseach[,j])/size(winning_partition_chunk,1))*(ones(length(linear_idx) ,1))
                }
                
                for (k in 1:length(fraction_counts_labels)){
                    if(k == 1){
                        fraction_counts_labels_db <- as.data.frame(fraction_counts_labels[[1]])
                    }else {
                        fraction_counts_labels_db <- rbind(fraction_counts_labels_db, as.data.frame(fraction_counts_labels[[k]]))
                    }
                }
                
                
               #fraction_counts_labels <- vertcat(fraction_counts_labels{:})
                [row col] in fast_ind2sub(size(partitions),vertcat(markerfound{all_labels_for_winning_cluster_unq(:)}))){
                nodes_x_partitions_coocvals <- sparse(row, col, fraction_counts_labels, size(partitions,1), size(partitions,2))# in nodes_x_partitions_coocvals rows are nodes and columns are partitions, values are mean of co-occurence with all members of a (winning) cluster using data from a non-winning partition
                nodes_x_partitions_coocvals(winning_partition_idxs,:) <- 0
                nodes_x_partitions_coocvals_holder{i,1} <- nodes_x_partitions_coocvals
    
                mutlicom_node_cell{i}       <- find(mean(nodes_x_partitions_coocvals,2)>accept_multi)#it's a bit slow to take these means sequentially, but I've tried doing them all at once and it ends up taking up too much ram for networks of 100K+
                cell_partition_overlapping{i,1} <- cat(1,mutlicom_node_cell{i},cell_partition{i})#i is the target cluster and the source is whatever cluster the node was originally in
                nodes_and_partition_identifiers_overlapping(cat(1,mutlicom_node_cell{i},cell_partition{i})) <- i){
    
            } else {#for small/banned modules
                mutlicom_node_cell{i} <- []
                cell_partition_overlapping{i,1} <- cell_partition{i}
    
            }
        }
    multicom_nodes_all <- vertcat(mutlicom_node_cell{:})
    save multi_com_list.mat multicom_nodes_all
    } else {
        cell_partition_overlapping <- cell_partition
        multicom_nodes_all <- []
        save multi_com_list.mat multicom_nodes_all#just so not confused with old result
    }
    
    if (accept_multi<1){
    printmessage(['overlapping vs discrete length: ' num2str(sum(cellfun(@length,cell_partition_overlapping))) ' vs ' num2str(sum(cellfun(@length,cell_partition)))],main_iter)
    }
    
    
    #from here on it's just arranging the output
    [trash, idx_large_partition] <- sort(cellfun('size', cell_partition, 1), 'desc}')
    cell_partition <- cell_partition(idx_large_partition)# reorder partitions by size, with no change to contents
    
    [trash, idx_large_partition] <- sort(cellfun('size', cell_partition_overlapping, 1), 'desc}')
    cell_partition_overlapping <- cell_partition_overlapping(idx_large_partition)# reorder partitions by size,
    
    
    
    cluster_density <- cell(length(cell_partition),1)#sort order of nodes within each cluster for display purposes
    for (i in 1:length(cell_partition)){
        cluster_density{i} <- mean(ADJ(cell_partition{i},cell_partition{i}))
        [tmp,ind] <- sort(cluster_density{i}, 'desc}')
        cell_partition{i} <- cell_partition{i}(ind)
    }
    #best_nodeorder_hard=vertcat(cell_partition{:});
    
    
    cluster_density_overlapping <- cell(length(cell_partition),1)#sort order of nodes within each cluster for display purposes
    for (i in 1:length(cell_partition_overlapping)){
        cluster_density_overlapping{i} <- mean(ADJ(cell_partition_overlapping{i},cell_partition_overlapping{i}))
        [tmp,ind] <- sort(cluster_density_overlapping{i}, 'desc}')
        cell_partition_overlapping{i} <- cell_partition_overlapping{i}(ind)
    }
    #best_nodeorder_hard=vertcat(cell_partition_overlapping{:});
    
    
    partition_marker_sorted_hard <- []
    partition_marker_sorted_overlapping <- []
    for (i in 1:length(cell_partition)){
        partition_marker_sorted_hard(}+1:}+length(cell_partition{i})) <- i
        partition_marker_sorted_overlapping(}+1:}+length(cell_partition_overlapping{i})) <- i
    }
    nodes_and_partition_identifiers_hard <- sortrows([vertcat(cell_partition{:}) partition_marker_sorted_hard'])){
    
    
    partition_marker_sorted_hard <- []
    for (i in 1:length(cell_partition)){
        partition_marker_sorted_hard(}+1:}+length(cell_partition_overlapping{i})) <- i
    }
    nodes_and_partition_identifiers_overlapping <- sortrows([vertcat(cell_partition_overlapping{:}) partition_marker_sorted_hard'])){
    
    
    record_stuff <- 0#cluster cood density stats
    if (record_stuff==1){
    
        partition_rows <- zeros(length(cell_partition),length(cooc))
        cooc_temp <- cooc
        for (i in 1:length(cell_partition)){
    
            partition_rows(i,:) <- sum(cooc_temp(cell_partition{i},:))./(length(cell_partition{i})-1)
    
        }
    
        [cluster row]=find(partition_rows!=0)
        value <- partition_rows(find(partition_rows>0))
    
        if (isempty(accept_multi)){
            csvwrite('SpeakEasy_cluster_assignment.csv',[row cluster value])
        }

}
#cluster_stats=cluster_stability(ADJ,cooc,cell_partition_overlapping);

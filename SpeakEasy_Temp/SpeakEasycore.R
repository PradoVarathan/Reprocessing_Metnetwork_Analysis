# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#core SpeakEasy routine that counts up actual labels, compares to
#predicted, and then selects most unexpected label for node identifier
#
#expected to be called bye bootstarp_sepak
#to enable clustering of larger sparse matrices, becuase we always need to
#get a predicted value for each label, we only compute predicted values for
#the labels a node actually receives, otherwise we'd need a full matrix of
#expectations this size of ADJ.  Even if the ADJ is small enough that we
#could do this, it's faster to do it this way if the matrix is very sparse.


SpeakEasycore <- function(ADJ,total_time,IC_store_relevant,nback,force_efficient){
    library(pracma)
    kin <- as.matrix(colSums(ADJ))
    
    
    
    aggregate_labels <- function(ADJ,labels,nback){
        library(pracma)
        labels <- as.vector(labels)
        labels_unq <- unique(labels)
        
        #sorted_labels=sortrows([labels (1:length(labels))']);  #two columns of labels and label locations
        indices <- sort(labels, index.return=TRUE)$ix#sort faster than sortrows
        temp <- rbind(labels,(1:length(labels)))
        temp <- t(temp)
        sorted_labels <- temp[indices,]
        transitions=matrix(c(0, pracma::finds(sorted_labels[2:nrow(sorted_labels),1]-sorted_labels[1:(nrow(sorted_labels)-1),1] != 0),dim(sorted_labels)[1]))#find sets of locations of the end of the previous label
        # node_identifiers=zeros(length(labels),1);
        # for i=1:size(transitions,1)-1 #relabel the... labels with sequential integers, starting with 1, which becomes a time-saver later
        #      node_identifiers(sorted_labels(transitions(i)+1:transitions(i+1),2))=i;  #maybe replace with cumsum on sprasemat based on transitions, if that's fast
        # end
        
        future_markers <- zeros(size(sorted_labels,1),1)
        future_markers[1] <- 1
        for (k in (transitions[2:length(transitions)-1])){
            if (k != 0){
                future_markers[k] = 1
            }
        }
        future_markers <- cumsum(future_markers)
        node_identifiers <- zeros(size(sorted_labels,1),1)
        node_identifiers[sorted_labels[,2]] = future_markers
        
        
        #idea is to consider the labels from different time-steps to be different
        #(even though they are not) then add the row of temp created with these pseudo-different labels, to add up the rows that are infact related to the #same core label
        #the obvious way to do this is:
        #temp=zeros(length(labels_unq), length(node_identifiers));
        #temp=bsxfun(@eq, sparse(1:length(labels_unq)).', node_identifiers');
        #but a faster way is:
        j = t(as.matrix(1:length(node_identifiers)))
        nodes_by_labels_all_times <- sparseMatrix(i= as.vector(node_identifiers), j=t(c(1:length(node_identifiers))), 
                                                  x = as.vector(t(ones(length(node_identifiers),1))),
                                                  symmetric = FALSE,repr="T")
        #we can do that becaue node identifiers are numbers sequentially starting at 1
        #in each row of "nodes_by_labels_all_times" we tick off positions (using a 1) where that label occurs in the full list of labels
        
        section_length <- dim(nodes_by_labels_all_times)[2]/nback
        #we've arranged the labels from several timesteps chunks arranged horizontally, and we want to sum each of those chunks independently and add them all back up (because the rows in each chunch are in fact synchronized)
        nodes_by_labels_all_times = as.matrix(nodes_by_labels_all_times)
        for (i in 1:nback){
            if (i==1){
                running_sum <- nodes_by_labels_all_times[,1:section_length]
            } else {
                start = pracma::ceil((section_length)*(i-1))
                end = floor(start+section_length-1)
                running_sum <- running_sum+nodes_by_labels_all_times[,start:end]
            }
        }
        
        lxn <- running_sum%*%ADJ
        return(lxn)
        #size(runing_sum,2)==length(ADJ), lnx will be sparse if ADJ is sparse
        #lxn has counts of each label (each row of lxn represents a different label) for each node (columns of lxn)
    }
    
    
    
    #initialize matrix to store chosen labels
    listener_history <- matrix(0,nback+total_time,max(dim(ADJ)))#each column is history for a single node, starting at row 1
    listener_history[1:nback+1,1:max(dim(ADJ))] <- IC_store_relevant
    #nback <- 
    for (i in (2+nback):(nback+total_time)){
        
        current_listener_history <- listener_history[(i-nback):(i-1),]
        temp_agr <- aggregate_labels(ADJ,current_listener_history,nback)#actual_counts is sparse for sparse input
        actual_counts <- temp_agr 
        active_labels <- unique(as.vector(current_listener_history))
        counts_normk <- colSums(actual_counts)  #',2)' needed for case of size-1 clusters
        count_normk_norm1 <- as.matrix(counts_normk/sum(counts_normk))#proportions of various labels normalized to 1
    
        #if matrix is very sparse, or too large to store a full ADJ, we only care about generating expected counts of labels if a node actually receives some of that label
        x_y <- which(actual_counts != 0, arr.ind = T)#x will be labels and y will be nodeID
        x <- x_y[,1]
        y <- x_y[,2]
            #two lines below are easier to understand but slightly slower separately
            #scaled_kin=nback*([full(count_normk_norm1(x))]'.* kin(y)); #scales normalized counts by total input (some nodes have more inputs and thus you would expect more of all labels)
            #expected=sparse(x,y,scaled_kin, size(actual_counts,1),size(actual_counts,2));  #same format as lxn
        expected <- sparseMatrix(i=as.vector(x),j=as.vector(t(y)),x=as.vector(nback*(as.matrix(count_normk_norm1[x])* kin[y])),
                                     symmetric = FALSE, repr='T',dims = c(length(unique(x)), dim(actual_counts)[2]))#same format as lxn
    
    
        # } else {
        #     # sum(expected)==kin*nback
        #     expected <- nback*t(as.matrix(count_normk_norm1))*(as.matrix(kin))#a bit slower for ADJ's with less than 3# density compared to option above (and more mem usage) but 10x faster on sparse
        # 
            #diagnostic
    #                 [x y]=find(actual_counts);   #x will be labels and y will be nodeID
    #                 expectedsparse=sparse(x,y,nback*([full(count_normk_norm1(x))]'.* kin(y)), size(actual_counts,1),size(actual_counts,2));  #same format as lxn
    #          full(actual_counts)
    #          full(expectedsparse)
    #          full(expected)
    #          full( expectedsparse-actual_counts)
    #          full( expected-actual_counts)
    #         length(find(min(( expectedsparse-actual_counts))>0))
       
        tem <- actual_counts-expected
        tem <- t(tem)
        max_vals <- apply(tem,1,max)
        max_idx <- apply(tem,1,which.max)
       # length(find(maxvals==0))   #might worry that for force_efficient==1 case max of some column will be zero, indicating a non-elegible label, but that never happens... could add 100 to each non-zero entry of actual_counts if worried, but really not necessary
       listener_history[i,] <-  active_labels[max_idx]
    
    }#time-step loop
    
    #identify nodes in same clusters
    sorted_labels <- t(sortrows(rbind(t(listener_history[(dim(listener_history)[1]-2),]),t(1:length(listener_history[(dim(listener_history)[1]-2),])))))
    transitions=matrix(c(0, pracma::finds(sorted_labels[(2:dim(sorted_labels)[1]),1]-sorted_labels[(1:dim(sorted_labels)[1]-1),1]!=0), size(sorted_labels,1)))#find sets of locations of the end of the previous label
    partitionID <- zeros(max(dim(ADJ)),1)#for all clusters, get their locations in the partitions and their numeric identifier
    label_assignment <-list()
    for (i in 1:(size(transitions,1)-1)){
        ids <- c((transitions[i]+1):transitions[i+1])
        partitionID[sorted_labels[ids,2]] <- i
        
         #this actually gets very big too
        label_assignment[[i]] <- sorted_labels[((transitions[i]+1):transitions[i+1]),2]
    }
    final_res <- list("label_assignment" = label_assignment, "partitionID" = partitionID, "listener_history" = listener_history)
    return(final_res)
}

# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#reports the number of each label that each cell "heard" over the past several time-steps

#in the returned matrix (lxn) rows represent labels and columns are all nodes in ADJ
#the values in lxn are the number of times a given label is an input for a given node

#the input called labels is actually the last nback labels for all nodes

#goal is to add up certain rows (inputs) of ADJ according to "labels"
#(corresponding to each row of ADJ)
#basically, add up rows that have same label identifier

#the main reason this get a little tricky is that A) it has to be done very
#quickly and B) some labels may be present at some time-steps but not at
#others

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

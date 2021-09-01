# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#SpeakEasy generate a history of labels for eah node.  Initially each node
#has no label history and this can lead to getting trapped in sub-optimal
#cluster states.  Therefore, we generate an estimated initial label history for each
#node at the start, much like setting up reasonable (but still somewhat random)
#initial conditions in an system of ODE's
#
#A more general point is the in label propagation algorithm SLPA by Xie et al. the randomness comes in the operation
#(labels are selected stochastically)
#We could of course do that in SpeakEasy, but at the moment, we put the
#randomness in the initial conditions because semi-randomly updating labels
#requires longer time to converge.
#
#We can generate initial conditions for all runs of SpeakEasy at a single
#time, which is much more efficient than creating them every separatly for
#each replicate run.


genIC_full <- function(ADJ,how_many_runs,nback,varargin){
    
    IC_store <- matrix(0,how_many_runs*(nback+1),max(dim(ADJ)))
    
    for (m in 1:max(dim(ADJ))){#setup initial listener histories with randomly selected neighbor labels
    
        contacts=which(ADJ[,m]!=0)
        if (length(contacts)==0){#particularly in very small modules, we might have no connections... usually diag(ADJ) is set to all ones, which also solves this, but in case it is not, this avoids an error message
            IC_store[,m] <- m
    
        } else {
            contacts_values <- ADJ[contacts,m]
    
            if (length(which(contacts_values==1))!=length(contacts_values)){#this should generally be the case for weighted networks
                contacts_weights <- contacts_values-min(contacts_values)
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
    
            bins <- min(c(0,cumsum(contacts_weights/sum(contacts_weights)),1))
            library(pracma)
            ret_hist <- histc(rand(how_many_runs*(nback+1),1),bins)
            IC_store[,m] <- contacts(ret_hist$bin)
            #using cellfun (@randsample or bsxfn isn't any faster
    
        }
    
    }
}


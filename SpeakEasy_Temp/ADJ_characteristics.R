# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#By checking the network density and weighted/unweighted format, we can
#potentially save time on later operations.
#For larger matrices we only sample their properties as it would be slow to
#check every link.

#Outputs:
#is_ADJ_weighted==1 for case when entries are not entirely 0 or 1
#force_efficient==1 if matrix is larger than max_ADJ_size or if it is smaller than max_ADJ_size and low density
ADJcharacteristics <- function(ADJ,max_ADJ_size, force_efficient){
    
    if (dim(ADJ)[1] != dim(ADJ)[2]){
        cat('your ADJ is not square, please fix')
        stop()
    }
    
    if ((max(ADJ))>1 || min(ADJ)<-1){
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
}
    

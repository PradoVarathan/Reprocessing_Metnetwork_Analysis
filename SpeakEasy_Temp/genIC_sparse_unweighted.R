# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#generates initial conditions for SpeakEasey
#only intended for unweighted sparse networks
#since the largest graphs with millions of nodes are likely not full or weighted, this
#is helpful as genIC_full is much slower
#
#this routine works by creating an equal number of random values for each
#node and scaling them by k, to generate the indices of neighbors

genIC_sparse_unweighted <- function(ADJ,how_many_runs,nback,varargin){
  library(pracma)
  k <- sum(ADJ)
  kcumsum <- cumsum(k)
  basicrand <- rand((nback+1)*how_many_runs,length(ADJ))
  indices <- which(ADJ != 0, arr.ind = TRUE)
  row_ind <- indices[,1]
  basicrand_scaled <- ceil(basicrand.*repmat(k,size(basicrand,1),1))
  basicrand_scaled_lifted <- basicrand_scaled+repmat(c(0,kcumsum(1:-1)),size(basicrand_scaled,1),1)
  IC_store <- Reshape(row(basicrand_scaled_lifted),(nback+1)*how_many_runs,length(ADJ))
  return(IC_store)
}




# Copyright 2015 Chris Gaiteri, Rush University, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin and Rensselaer Polytechnic Institute. All worldwide rights reserved. A license to use, copy, and modify and distribute this software for non-commercial research purposes only is hereby granted, provided that this copyright notice and accompanying disclaimer is not modified or removed from the software.

# Any paper published with results obtained using this software should  provide a statement identifying the algorithm and providing a citation to:

# Identifying robust communities and multi-community nodes by combining topdown and bottom-up approaches to clustering, Chris Gaiteri, Mingming Chen, Boleslaw Szymanski, Konstantin Kuzmin, Jierui Xie, Changkyu Lee, Timothy Blanche, Elias Chaibub Neto, Su-Chun Huang, Thomas Grabowski, Tara Madhyastha and Vitalina Komashko,Scientific Reports 5 Article number: 16361, 2015.

# DISCLAIMER: The software is distributed "AS IS" without any express or implied warranty, including but not limited to, any implied warranties of merchantability or fitness for a particular purpose or any warranty of non-infringement of any current or pending patent rights. The authors of the software make no representations about the suitability of this software for any particular purpose. The entire risk as to the quality and performance of the software is with the user. Should the software prove defective, the user assumes the cost of all necessary servicing, repair or correction. In particular, neither Rush University, Rensselaer Polytechnic Institute, nor the authors of the software are liable for any indirect, special, consequential, or incidental damages related to the software, to the maximum extent the law permits.




#posted here http://www.mathworks.com/matlabcentral/fileexchange/13916--simple--tool-for-estimating-the-number-of-clusters/content/valid_RandIndex.m
#by by Kaijun Wang
#10 Feb 2007 (Updated 08 Jul 2009)
#patched to make fast enough for huge cluster results by C.Gaiteri 2014

#RANDINDEX - calculates Rand Indices to compare two partitions
# ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the
# class membership, returns the "Hubert & Arabie adjusted Rand index".
# See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of
# Classification 2:193-218
#(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk
#Copyright (c) 2006-2007, Kaijun WANG
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:#
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the distribution

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.

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
  
  if (t1==nc){
     ARI <- 0#avoid division by zero; if k=1, define Rand = 0
  } else {
     ARI <- (A-nc)/(t1-nc)#adjusted Rand - Hubert & Arabie 1985
  }
  return(ARI)
}

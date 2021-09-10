ARG BASE_IMAGE=rstudio/r-base
FROM ${BASE_IMAGE}:centos7

ARG R_VERSION=3.6.3
ARG OS_IDENTIFIER=centos-7

# Install R
RUN wget https://cdn.rstudio.com/r/${OS_IDENTIFIER}/pkgs/R-${R_VERSION}-1-1.x86_64.rpm && \
    yum -y install epel-release && \
    yum -y install ./R-${R_VERSION}-1-1.x86_64.rpm && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm R-${R_VERSION}-1-1.x86_64.rpm && \
    yum -y remove epel-release && \
    yum clean all

RUN cd ~ && \
    sudo yum update -y && \
    echo "export PS1=\"\n\[$(tput bold)\]\[$(tput sgr0)\]\[\033[38;5;38m\]\u\[$(tput sgr0)\]\[$(tput sgr0)\]\[\033[38;5;15m\]@\h[\[$(tput sgr0)\]\[\033[38;5;214m\]\W\[$(tput sgr0)\]\[\033[38;5;15m\]]: \[$(tput sgr0)\]\"" >> ~/.bashrc && \
    . ~/.bashrc && \
    mkdir /shared && \
    mkdir /shared/rlibs && \
    sudo chown root /shared/rlibs && \
    sudo yum install -y openmpi && \
    sudo yum install -y openmpi-devel && \
    sudo yum install -y git && \
    sudo yum install -y nano && \
    sudo yum -y install libcurl libcurl-devel && \
    sudo yum -y install openssl-devel && \
    sudo yum -y install libxml2-devel && \
    sudo yum install -y libpng-devel && \
    sudo yum install -y libjpeg-turbo-devel && \
    sudo yum -y groupinstall 'Development Tools' && \
    sudo yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel

RUN echo "local({r <- getOption('repos') \nr['CRAN'] <- 'http://cran.r-project.org' \noptions(repos=r)})" >> /opt/R/3.6.3/lib/R/etc/Rprofile.site \
    && R -e "install.packages('BiocManager')"\
    && R -e "install.packages('dplyr')"\
    && R -e "install.packages('glmnet')"\
    && R -e "install.packages('randomForest')"\
    && R -e "install.packages('Hmisc')"\
    && R -e "install.packages('lars')"\
    && R -e "BiocManager::install('WGCNA')"\
    && R -e "install.packages('devtools')\
    && R -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
    && R -e "devtools::install_github('PradoVarathan/metanetwork')"\
    && R -e "devtools::install_github('brian-bot/githubr')"\
    && R -e "install.packages('c3net')"\
    && R -e "install.packages('config')"\
    && R -e "install.packages('optparse')"\
    && R -e "install.packages('parmigene')"\
    && R -e "install.packages('data.table')"\
    && R -e "install.packages('Rmpi')"\
    && R -e "install.packages('Matrix')" \
    && R -e "install.packages('reader')"\
    && R -e "install.packages('caroline')"\
    && R -e "install.packages('bench')" \
    && R -e "install.packages('doParallel')"\
    && R -e "install.packages('parallel')"\
    && R -e "install.packages('snow')"\
    && R -e "devtools::install_github('blogsdon/utilityFunctions')" \
    && R -e "install.packages('bit64')" \
    && module load mpi/openmpi-x86_64  \
    && mkdir /shared/rlibs \
    && Rscript -e 'install.packages(c("docopt", "devtools", "dplyr", "tidyr", "ggplot2", "reshape2", "knitr", "stringr", "readr", "plyr", "data.table", "rJava", "doParallel", "snow", "igraph", "Rcpp", "RcppEigen", "Rclusterpp", "RColorBrewer", "MRCE", "vbsr", "ctv", "psych", "reshape2", "vcd", "erer", "fpc", "pacman"))'\
    && Rscript -e 'install.packages(c("RJSONIO", "RCurl", "digest")); source("http://depot.sagebase.org/CRAN.R"); pkgInstall(c("synapseClient"))' \
    && Rscript -e 'source("http://www.bioconductor.org/biocLite.R") ; biocLite(c("limma", "biovizBase", "e1071", "org.Hs.eg.db", "edgeR", "AnnotationDbi", "biomaRt", "ComplexHeatmap", "FDb.InfiniumMethylation.hg19", "RDAVIDWebService", "topGO", "goseq", "GO.db", "GSVA", "org.Hs.eg.db", "edgeR", "limma", "CePa", "Biobase", "pracma", "annotate", "BiocInstaller", "Biostrings", "GEOquery", "GOstats", "graph", "GSEABase", "impute", "preprocessCore"))'\
    && Rscript -e 'library(devtools); install_github("th1vairam/CovariateAnalysis@dev")'\
    && export LD_LIBRARY_PATH=/usr/lib64/openmpi-1.10/lib:$LD_LIBRARY_PATH \
    && Rscript -e 'install.packages("Rmpi", configure.args = paste("--with-Rmpi-include=/usr/include/openmpi-x86_64","--with-Rmpi-libpath=/usr/lib64/openmpi/lib","--with-Rmpi-type=OPENMPI"))'\
    && Rscript -e 'install.packages("idr"))'

RUN mkdir ~/metanetwork \
    && mkdir ~/Analysis_Scripts \
    && mkdir /home/test_res/ \
    && git clone https://github.com/PradoVarathan/metanetwork ~/metanetwork \
    && mkdir /home/test_res/.openmpi \
    && touch /home/test_res/.openmpi/mca-params.conf \
    && echo "orte_default_hostfile=/home/test_res/default_host" >> /home/test_res/.openmpi/mca-params.conf \
    && echo "localhost slots=32" >> /home/test_res/default_host


FROM rocker/tidyverse:4.0.0

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libglpk40\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y libopenmpi-dev\
&& apt-get install -y git\
&& R -e "install.packages('BiocManager')"\
&& R -e "install.packages('dplyr')"\
&& R -e "install.packages('glmnet')"\
&& R -e "install.packages('randomForest')"\
&& R -e "install.packages('Hmisc')"\
&& R -e "install.packages('lars')"\
&& R -e "BiocManager::install('WGCNA')"\
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
&& R -e "install.packages('parallel')"

RUN mkdir ~/metanetwork \
&& mkdir ~/Analysis_Scripts \
&& git clone https://github.com/PradoVarathan/metanetwork ~/metanetwork
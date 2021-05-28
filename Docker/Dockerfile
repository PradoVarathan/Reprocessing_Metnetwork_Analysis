FROM rocker/tidyverse:4.0.0

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libglpk40\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y git\
&& R -e "install.packages('BiocManager')"\
&& R -e "install.packages('dplyr')"\
&& R -e "install.packages('glmnet')"\
&& R -e "install.packages('randomForest')"\
&& R -e "install.packages('Hmisc')"\
&& R -e "install.packages('lars')"\
&& R -e "BiocManager::install('WGCNA')"\
&& R -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
&& R -e "devtools::install_github('blogsdon/metanetwork')"\
&& R -e "devtools::install_github('brian-bot/githubr')"\
&& R -e "install.packages('c3net')"\
&& R -e "install.packages('config')"\
&& R -e "install.packages('optparse')"

RUN mkdir ~/metasample \
&& mkdir ~/Analysis_Scripts \
&& git clone https://github.com/Sage-Bionetworks/metanetwork ~/metasample\
&& git clone https://github.com/PradoVarathan/Reprocessing_Metnetwork_Analysis ~/Analysis_Scripts
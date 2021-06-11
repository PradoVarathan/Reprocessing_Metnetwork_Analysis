# Reprocessing_Metnetwork_Analysis

## This repository is for preparing code and pipeline for the network analysis prior to uploading to forked Metanetwork analysis on Sage Bionetworks.

### Getting Started
Create your docker container run environment with port forwarded IP address

```
mkdir /home/<USER>/
cd /home/<USER>/

git clone git@github.com:PradoVarathan/Reprocessing_Metnetwork_Analysis.git
```

Build your docker image and container

```
docker image build -t metanet Reprocessing_Metnetwork_Analysis/Docker/

docker run -v "/home/<USER>/Reprocessing_Metnetwork_Analysis/:/home/<USER>/Reprocessing_Metnetwork_Analysis/" -e USER=<USER> -e PASSWORD=<PASSWORD> -d -p 8787:8787 metanet
```
 
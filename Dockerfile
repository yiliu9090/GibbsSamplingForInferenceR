FROM r-base

RUN mkdir /Workspace
COPY Rprofile /Workspace
COPY RunMCMC.R /Workspace 

WORKDIR  /Workspace

RUN Rscript Initialpackage.R



ENTRYPOINT [ "Rscript" , "RunMCMC.R" ]
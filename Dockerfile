FROM r-base



RUN mkdir /Workspace
RUN apt-get update
RUN apt-get install -y libgsl0-dev

COPY Rprofile /Workspace
COPY RunMCMC.R /Workspace 

WORKDIR  /Workspace

RUN Rscript Initialpackage.R



ENTRYPOINT [ "Rscript" , "RunMCMC.R" ]
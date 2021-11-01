FROM r-base

RUN mkdir /Workspace
COPY Rprofile /Workspace
COPY RunMCMC.R /Workspace

RUN Rscript /Workspace/Initialpackage.R

WORKDIR   /Workspace

ENTRYPOINT [ "Rscript" , "RunMCMC.R" ]
FROM r-base

RUN mkdir /Workspace
COPY . /Workspace

RUN Rscript Rprofile/Initialpackage.R

WORKDIR   /Workspace

ENTRYPOINT [ "Rscript" , "RunMCMC.R" ]
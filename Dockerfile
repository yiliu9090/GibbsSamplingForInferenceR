FROM r-base



RUN mkdir /Workspace
RUN apt-get update
RUN apt-get install -y libgsl0-dev

COPY Rprofile /Workspace
COPY RunMCMC.R /Workspace 
COPY gsl.cpp /Workspace

COPY DirichletReg_0.7-1.tar.gz /Workspace
COPY jsonlite_1.7.2.tar.gz /Workspace
COPY matrixStats_0.61.0.tar.gz /Workspace
COPY maxLik_1.5-2.tar.gz /Workspace
COPY Rcpp_1.0.7.tar.gz /Workspace
COPY Rcpp11_3.1.2.0.1.tar.gz /Workspace
COPY RcppGSL_0.3.10.tar.gz /Workspace
COPY stats19_2.0.0.tar.gz /Workspace
COPY Formula_1.2-4.tar.gz /Workspace
COPY generics_0.1.1.tar.gz /Workspace
COPY miscTools_0.6-26.tar.gz /Workspace
COPY sandwich_3.0-1.tar.gz /Workspace
COPY digest_0.6.28.tar.gz /Workspace
COPY zoo_1.8-9.tar.gz /Workspace

WORKDIR  /Workspace

RUN Rscript Initialpackage.R



ENTRYPOINT [ "Rscript" , "RunMCMC.R" ]
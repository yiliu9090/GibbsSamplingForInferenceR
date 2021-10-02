

FROM r-base

RUN mkdir /Workspace
COPY Initialpackage.R /Workspace

WORKDIR   /Workspace

ENTRYPOINT [ "R" , "Entry.py" ]
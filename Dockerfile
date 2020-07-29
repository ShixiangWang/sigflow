FROM r-base:4.0.2

LABEL \
    author="Shixiang Wang" \
    email="w_shixiang@163.com" \
    description="Docker Image for Sigflow" \
    version="SigFlow v0.1 based on Sigminer v1.0.9 (platform R v4.0.2)"

## Install system dependencies
RUN apt update -y && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev &&  \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
## Install R packages which are easy to install
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
## Install R packages which are not easy to install
RUN R -e "BiocManager::install('BSgenome')" && \
    R -e "BiocManager::install('ShixiangWang/sigminer@v1.0.9', dependencies = TRUE)" && \
    R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.mm10'))" && \
    rm -rf /tmp/* /var/tmp/*
## Copy sigflow program and run test
## It is strange that the docopt cannot be installed to the first location
COPY sigflow.R /opt/
COPY ./test/ /opt/test/
RUN R -e "install.packages('docopt', lib = .libPaths()[2])" && \
    chmod u+x /opt/sigflow.R && ln -s /opt/sigflow.R /usr/bin/sigflow && \
    cd /opt/test && chmod u+x test.sh && ./test.sh && rm -rf test_results && cd /root
WORKDIR /root
## Deploy
## When ENTRYPOINT is used, the docker can be only run as a command
#ENTRYPOINT [ "sigflow" ]
#CMD [ "--help" ]
CMD [ "sigflow", "--help" ] 
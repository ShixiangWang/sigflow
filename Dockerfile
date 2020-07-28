FROM r-base:4.0.2

LABEL \
    author="Shixiang Wang" \
    email="w_shixiang@163.com" \
    description="Docker Image for Sigflow" \
    version="v0.1.0 build on R v4.0.2"

## Install system dependencies
RUN apt update -y && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev &&  \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
## Install R packages which are easy to install
RUN R -e "install.packages('docopt', repos = 'https://cloud.r-project.org')" && \
    R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
## Install R packages which are not easy to install
RUN R -e "BiocManager::install('BSgenome', dependencies = TRUE)" && \
    R -e "BiocManager::install('sigminer', dependencies = TRUE)" && \
    rm -rf /tmp/* /var/tmp/*
    ## Support all genomes directly
    #R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.mm10'))"
COPY sigflow.R /opt/
RUN chmod u+x /opt/sigflow.R && ln -s /opt/sigflow.R /usr/bin/sigflow
# Run test
COPY ./test/ /opt/
RUN ls -l /opt && cd /opt/test && chmod u+x test.sh && ./test.sh && rm -rf test_results && cd /root
WORKDIR /root
ENTRYPOINT [ "sigflow" ]
CMD [ "--help" ]

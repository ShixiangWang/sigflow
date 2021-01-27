FROM r-base:4.0.3

LABEL \
    author="Shixiang Wang" \
    maintainer="Shixiang Wang" \
    email="wangshx@shanghaitech.edu.cn" \
    description="Docker Image for Sigflow" \
    org.label-schema.license="Academic Free License v.3.0" \
    org.label-schema.vcs-url="https://github.com/ShixiangWang/sigflow/" \
    org.label-schema.vendor="XSLiu Lab Project"

## Install system dependencies
RUN apt update -y && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev &&  \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*
## Install R packages which are easy to install
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))"
## Install reference genome packages which are big
RUN R -e "BiocManager::install('BSgenome')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" && \
    R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
## Install sigminer & SigProfiler & test SigProfiler
RUN R -e "BiocManager::install('ShixiangWang/sigminer@v1.2.3', dependencies = TRUE)" && \
    rm -rf /tmp/* /var/tmp/* && \
    R -e "library('sigminer'); load(system.file('extdata', 'toy_copynumber_tally_M.RData', package = 'sigminer', mustWork = TRUE)); mat = cn_tally_M[['nmf_matrix']]; print(mat); sigprofiler_extract(mat, '/opt/test_sp_install', range = 3:4, nrun = 2L, use_conda = TRUE); cat(paste(list.files('/opt/test_sp_install', recursive = TRUE), '\n')); sigprofiler_import('/opt/test_sp_install')" && \
    rm -rf /opt/test_sp_install && \
    /root/.local/share/r-miniconda/bin/conda clean -afy
## Copy sigflow program and run test
## It is strange that the docopt cannot be installed to the first location
ENV PATH /root/.local/share/r-miniconda/bin:$PATH
COPY sigflow.R pkg_check.R /opt/
COPY ./test/ /opt/test/
RUN R --vanilla -f /opt/pkg_check.R && \
    R -e "install.packages('docopt', lib = .libPaths()[2])" && \
    chmod u+x /opt/sigflow.R && ln -s /opt/sigflow.R /usr/bin/sigflow && \
    cd /opt/test && chmod u+x test.sh && ./test.sh && rm -rf test_results && cd /root
WORKDIR /root
## Deploy
## When ENTRYPOINT is used, the docker can be only run as a command
## unless specify --entrypoint /bin/bash to access docker terminal
## see: https://phoenixnap.com/kb/docker-run-override-entrypoint
## sudo docker run -it --rm --entrypoint bash shixiangwang/sigflow
ENTRYPOINT [ "sigflow" ]
CMD [ "--help" ]
## Can entrypoint directly input outside files?
## CMD ["sigflow", "--help"]

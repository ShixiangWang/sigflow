FROM r-base:4.0.2

LABEL \
    author="Shixiang Wang" \
    email="w_shixiang@163.com" \
    description="Docker Image for Sigflow" \
    version="v0.1.0 build on R v4.0.2"

RUN apt install libcurl4-openssl-dev libxml2-dev libssl-dev &&  \
    R -e "install.packages('docopt', repos = 'https://cloud.r-project.org')" && \
    R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install('sigminer', dependencies = TRUE)"
ADD sigflow.R /opt
RUN chmod u+x /opt/sigflow.R && ln -s /opt/sigflow.R /usr/bin/sigflow
WORKDIR /root
CMD [ "sigflow" ]


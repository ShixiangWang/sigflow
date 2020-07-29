#!/usr/bin/env bash
CURRENT_DIR=$(cd `dirname $0`; pwd)
sudo docker run --rm shixiangwang/sigflow extract -i /opt/test/tcga_laml.maf.gz -o /opt/test/test_results/test_docker -m SBS -r 3
sudo docker run --rm -v $CURRENT_DIR:/root/test/ shixiangwang/sigflow extract -i /root/test/tcga_laml.maf.gz -o /root/test/test_results/test_docker -m SBS -r 3

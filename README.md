# Sigflow: Streamline Analysis Workflows for Mutational Signatures

<img alt="GitHub repo size" src="https://img.shields.io/github/repo-size/shixiangwang/sigminer.workflow"> <img alt="Docker Automated build" src="https://img.shields.io/docker/automated/shixiangwang/sigflow"> <img alt="Docker Image Version (latest by date)" src="https://img.shields.io/docker/v/shixiangwang/sigflow?color=blue"> <img alt="Docker Image Size (latest by date)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow"> <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/shixiangwang/sigflow">

```
       .__        _____.__                 
  _____|__| _____/ ____\  |   ______  _  __
 /  ___/  |/ ___\   __\|  |  /  _ \ \/ \/ /
 \___ \|  / /_/  >  |  |  |_(  <_> )     / 
/____  >__\___  /|__|  |____/\____/ \/\_/  
     \/  /_____/  
```

**Sigflow** provides useful mutational signature analysis workflows based on R package [sigminer](https://github.com/ShixiangWang/sigminer). It can auto-extract mutational signatures,
fit mutation data to COSMIC reference signatures (SBS/DBS/INDEL) and run bootstrapping analysis for
signature fitting.

Skip the following installation step if you would like to use [Docker](https://hub.docker.com/r/shixiangwang/sigflow).

<details>
<summary>Table of content</summary>
       
- [Installation](#installation)
- [Use Sigflow docker image](#use-sigflow-docker-image)
- [Usage](#usage)
- [Use cases](#use-cases)
- [Updates](#updates)
- [Test](#test)
- [Citation](#citation)
- [License](#license)

</details>

## Installation

Requirements:

- R (>3.6).
- R package [sigminer](https://github.com/ShixiangWang/sigminer).
- R package [docopt](https://cran.r-project.org/package=docopt).

Steps:

1. Install R - follow the instructions at <https://cran.r-project.org/>.
2. Install R packages, run

```r
install.packages("docopt")
install.packages("BiocManager")
BiocManager::install("sigminer", dependencies = TRUE)
# Update Sigminer version
install.packages("remotes")
remotes::install_github("ShixiangWang/sigminer")
# Install specific version by
# remotes::install_github("ShixiangWang/sigminer@v1.0.10")
```

3. Clone this repository, run

```bash
$ git clone https://github.com/ShixiangWang/sigminer.workflow
```

4. Link the R script as a executable file (command)

```bash
$ cd sigminer.workflow
$ ln -s $PWD/sigflow.R /usr/bin/sigflow  # You can choose another place instead of /usr/bin/sigflow
```

5. Try calling `sigflow`

```bash
sigflow -h
```

> Maybe you need to restart your terminal.

## Use Sigflow docker image

Use specified version (recommended way):

```bash
# docker pull shixiangwang/sigflow:version, e.g.
$ docker pull shixiangwang/sigflow:0.1
```

> NOTE: Sigflow version has no prefix `v`.

Current available tag versions:

- <img alt="Docker Image Version (tag latest semver)" src="https://img.shields.io/docker/v/shixiangwang/sigflow/0.1?color=blue"> <img alt="MicroBadger Layers (tag)" src="https://img.shields.io/microbadger/layers/shixiangwang/sigflow/0.1"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/0.1">
- <img alt="Docker Image Version (tag latest semver)" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.0?color=blue"> <img alt="MicroBadger Layers (tag)" src="https://img.shields.io/microbadger/layers/shixiangwang/sigflow/1.0"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.0">

Use latest version:

```bash
$ docker pull shixiangwang/sigflow:latest
```

> The latest version uses the latest commit from GitHub, so it may have not been
> prepared or tested. **So, be careful!**

Run the docker by:

```bash
$ docker run shixiangwang/sigflow
```

See [test/test_docker.sh](test/test_docker.sh) for examples.

If you want to go into the docker terminal, run

```bash
$ docker run --rm --entrypoint /bin/bash -it shixiangwang/sigflow
```

## Usage

### Commands and options

All Sigflow commands and options are described as the following.

```
=================================================================
sigflow: Streamline Analysis Workflows for Mutational Signatures.

Author: Shixiang Wang (wangshx@shanghaitech.edu.cn)
Copyright: AFL@2020 [https://opensource.org/licenses/AFL-3.0]

Desc:
  There are several subcommands.
  ==
  extract - extract signatures by either automated or manual way.
            Of note, when you use manual way, you need to run 2 times, 
            firstly you should set --manual to get signature estimation results,
            and secondly you should set --manual --number N to get N signatures.
  ==
  fit     - fit signatures in >=1 samples based on COSMIC reference signatures.
  ==
  bt      - run bootstrap signature fitting analysis in >=1 samples.

Usage:
  sigflow extract --input=<file> [--output=<outdir>] [--mode=<class>] [--manual --number <sigs>] [--max <max>] [--genome=<genome>] [--nrun=<runs>] [--cores=<cores>] [--sigprofiler] [--refit] [--hyper] [--verbose]
  sigflow fit --input=<file> [--output=<outdir>] [--index=<index>] [--mode=<class>] [--genome=<genome>] [--verbose]
  sigflow bt  --input=<file> [--output=<outdir>] [--index=<index>] [--mode=<class>] [--genome=<genome>] [--nrun=<runs>] [--verbose]
  sigflow (-h | --help)
  sigflow --version

Options:
  -h --help     Show help message.
  --version     Show version.
  -i <file>, --input <file>       input CSV/EXCEL/MAF file or VCF directory path.
  -o <outdir>, --output <outdir>  output directory path [default: ./sigflow_result/].
  --index <index>                 reference signature index separated by comma, e.g. '1,2,3' [default: ALL].
  -m <class>, --mode <class>      extract/fit mode, can be one of SBS, DBS, ID, MAF (for three types), CN (not supported in fit subcommand) [default: SBS].
  --manual                        enable manual extraction, set -N=0 for outputing signature estimation firstly.
  -N <sigs>, --number <sigs>      extract specified number of signatures [default: 0].
  --max <max>                     maximum signature number, default is auto-configured, should >2 [default: -1].
  -g <genome>, --genome <genome>  genome build, can be hg19, hg38 or mm10, [default: hg19].
  -r <runs>, --nrun <runs>        run times of NMF (extract) or bootstrapping (bt) to get results [default: 30].
  -T <cores>, --cores <cores>     cores to run the program, large dataset will benefit from it [default: 1].
  --refit                         refit the denovo signatures with quadratic programming or nnls.
  --hyper                         enable hyper mutation handling in COSMIC signatures (not used by SigProfiler approach).
  --sigprofiler                   enable auto-extraction by SigProfiler software.
  -v, --verbose                   print verbose message.

=================================================================
```

### Input

Sigflow supports input data in VCF/MAF/CSV/EXCEL format. The file format is auto-detected by Sigflow.

For SBS/DBS/INDEL data in CSV (including TSV) or EXCEL format, the following columns typically described in MAF format are necessary:

- `Hugo_Symbol`: gene symbol
- `Chromosome`: chromosome name, e.g. "chr1"
- `Start_Position`: start positionof the variant (1-based)
- `End_Position`: end position of the variant (1-based) 
- `Reference_Allele`: reference allele of the variant, e.g. "C"
- `Tumor_Seq_Allele2`: tumor sequence allele, e.g. "T"
- `Variant_Classification`: variant classification, e.g. "Missense_Mutation"
- `Variant_Type`: variant type, e.g. "SNP"
- `Tumor_Sample_Barcode`: sample identifier

For copy number segment data in in CSV (including TSV) or EXCEL format, the following columns are necessary:

- `Chromosome`: chromosome name, e.g. "chr1"
- `Start.bp`: start breakpoint position of segment
- `End.bp`: end breakpoint position of segment
- `modal_cn`: integer copy number value
- `sample`: sample identifier

## Use cases

Example datasets along with many example code are available in clone repository above (you can read it online at [here](https://github.com/ShixiangWang/sigminer.workflow/tree/master/test)).

The following parts give an example for each command.

Result directory of any command has the following structure.

- Files with extension`.RData` and `.rds`  are R related files to reproduce the results, and can be imported into R for further analysis and visualization.
- Files with extension`.pdf` are common visualization results used for communication.
- Files with extension `.csv` are formated data tables used for inspection, communication or further analysis.

#### `extract` command

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigminer.workflow
$ cd sigminer.workflow/test
$ sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf -m MAF -r 10 -T 4 --max 10
```

This will auto-extract SBS/DBS/INDEL signatures from data `toga_laml.maf.gz` by 10 Bayesian NMF runs with 4 computer cores (4 threads) from signature number ranges from 1 to 10, output results to directory `test_results/test_maf`.

> NOTE, in practice, set `-r` to a value  `>=10` is recommended for auto-extraction with Bayesian NMF, `>=100` for semi-automatic extraction and automatic extraction with SigProfiler (enabled by `--sigprofiler`).

#### `fit` command

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigminer.workflow
$ cd sigminer.workflow/test
$ sigflow fit -i tcga_laml.maf.gz -o test_results/test_fitting -m MAF
```

This will auto-fit input data `tcga_laml.maf.gz` to COSMIC SBS/DBS/INDEL signatures. Signature exposure data tables and plots are outputed.

#### `bt` command

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigminer.workflow
$ cd sigminer.workflow/test
$ sigflow bt -i tcga_laml.maf.gz -o test_results/test_bt -m SBS -r 5
```

This will auto-fit the random resample of input mutation profile to COSMIC SBS/DBS/INDEL signatures for specified times (here is 5). Data tables and plots of bootstrap signature exposures, errors and p values under different exposure cutoff are outputed.

> NOTE, in practice, set `-r` to a value  `>=100` is recommended.

#### How to use Docker to run Sigflow

If you use Docker to run Sigflow, you cannot directly call `sigflow` command. Instead, you should use `sudo docker run --rm -v /your_local_path:/docker_path shixiangwang/sigflow` to start a Docker container.

For example, if you want to accomplish the same task shown in `extract` command above, you need to run:

```sh
$ sudo docker run --rm -v /your_local_path:/docker_path shixiangwang/sigflow \
  extract -i /docker_path/tcga_laml.maf.gz \
          -o /docker_path/test_maf \
          -m MAF -r 10 -T 4 --max 10
```

Here,

- `--rm` will delete this container when this task is finished.

- `-v` is used for mounting your local directory `/your_local_path` as `/docker_path` in Docker image. **This is important**. You need to use the Docker container path in Sigflow arguments. So there must be a file called `/your_local_path/tcga_laml.maf.gz` exists in your computer, it will be treated as `/docker_path/test_maf` in the container.

## Updates

- 2020-09-03: 
  - Use sigminer v1.0.15 and support inputing reference signature index in `fit` and `bt` commands.
  - Allow users to decide if refit the signature exposures after *de novo* extraction with `refit` option.
  - Support matrix as input.
- 2020-08-14: Use sigminer v1.0.11 to use SigProfilerExtractor v1.0.15 to avoid issue raised from SigProfilerExtractor updates.
- 2020-08-05: **Release Sigflow 1.0** and related Docker image. This version is based on Sigminer v1.0.10, R v4.0.2 and SigProfilerExtractor v.1.0.15.
  - Supported SigProfiler.
  - Added `verbose` option.
  - Added `max` option.
  - Added `hyper` option.
  - More flexible and reasonable configuration.
- 2020-07-29: **Release Sigflow 0.1** using Docker. Sigflow 0.1 is based on Sigminer v1.0.9 and R v4.0.2



## Test

There are some example data sets in this repository, you can find how to test different workflows in [test/test.sh](test/test.sh).
It is time consuming to run all tests, just pick an example test similar to your task and see how it works. Before releasing a new version of **Sigflow**, I would run all tests to make sure they work well.

## Citation

The **Sigflow** has not been formly published, if you are using **Sigflow** fow now, please cite:

Sigflow: an automated and comprehensive pipeline for cancer genome mutational signature analysis.
Shixiang Wang, Ziyu Tao, Tao Wu, Xue-Song Liu.
bioRxiv 2020.08.12.247528; doi: https://doi.org/10.1101/2020.08.12.247528

## License

This software is released under [Academic Free License ("AFL") v.3.0](https://opensource.org/licenses/AFL-3.0)

Copyright 2020 © Shixiang Wang, Xue-Song Liu

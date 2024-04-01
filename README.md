# Sigflow: Streamline Analysis Workflows for Mutational Signatures

<img alt="GitHub repo size" src="https://img.shields.io/github/repo-size/shixiangwang/sigflow"> <img alt="Docker Automated build" src="https://img.shields.io/docker/automated/shixiangwang/sigflow"> <img alt="Docker Image Version (latest by date)" src="https://img.shields.io/docker/v/shixiangwang/sigflow?color=blue"> <img alt="Docker Image Size (latest by date)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow"> <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/shixiangwang/sigflow">

```
       .__        _____.__                 
  _____|__| _____/ ____\  |   ______  _  __
 /  ___/  |/ ___\   __\|  |  /  _ \ \/ \/ /
 \___ \|  / /_/  >  |  |  |_(  <_> )     / 
/____  >__\___  /|__|  |____/\____/ \/\_/  
     \/  /_____/  
```

**Sigflow** provides useful mutational signature analysis workflows. It can auto-extract mutational signatures,
fit mutation data to all/specified COSMIC reference signatures (SBS/DBS/INDEL) and run bootstrapping analysis for
studying signature stability. 

> Any bugs or suggestions please report to [GitHub issues](https://github.com/ShixiangWang/sigflow/issues), I will respond as soon as possible.

<details>
<summary>Table of content</summary>

- [Sigflow: Streamline Analysis Workflows for Mutational Signatures](#sigflow-streamline-analysis-workflows-for-mutational-signatures)
  - [Installation](#installation)
  - [Use Sigflow docker image](#use-sigflow-docker-image)
  - [Usage](#usage)
    - [Commands and options](#commands-and-options)
    - [Input](#input)
  - [Use cases](#use-cases)
      - [`extract` command](#extract-command)
      - [`fit` command](#fit-command)
      - [`bt` command](#bt-command)
      - [`show` command](#show-command)
      - [How to use Docker to run Sigflow](#how-to-use-docker-to-run-sigflow)
  - [Updates](#updates)
  - [Test](#test)
  - [Other tools](#other-tools)
  - [Citation](#citation)
  - [License](#license)

</details>

## Installation

If you would like to use [Docker](https://hub.docker.com/r/shixiangwang/sigflow), skip the following installation step and go to PART ['Use Sigflow docker image'](#use-sigflow-docker-image) directly.

> Using Sigflow Docker image is recommended for users without experiences in programming, especially in R.

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
# remotes::install_github("ShixiangWang/sigminer@v1.0.17")
```

3. Clone this repository, run

```bash
$ git clone https://github.com/ShixiangWang/sigflow
```

4. Link the R script as a executable file (command)

```bash
$ cd sigflow
$ ln -s $PWD/sigflow.R /usr/local/bin/sigflow  # You can choose another place instead of /usr/bin/sigflow
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
$ docker pull shixiangwang/sigflow:1.0
```

> NOTE: Sigflow version has no prefix `v`.

Current available tag versions:

- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/0.1?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/0.1">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.0?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.0">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.1?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.1">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.2?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.2">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.3?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.3">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/1.5?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/1.5">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/2.0?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/2.0">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/2.1?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/2.1">
- <img alt="Docker Image Version" src="https://img.shields.io/docker/v/shixiangwang/sigflow/2.2?color=blue"> <img alt="Docker Image Size (tag)" src="https://img.shields.io/docker/image-size/shixiangwang/sigflow/2.2">

Use latest version:

```bash
$ docker pull shixiangwang/sigflow:latest
```

> The latest version uses the latest (successful build) commit from GitHub, so it may have not been
> prepared or fully tested. **So, be careful!**

Run the docker by:

```bash
$ docker run shixiangwang/sigflow
```

See [test/test_docker.sh](test/test_docker.sh) for examples.

If you want to go into the docker container terminal, run

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
  extract - extract signatures by either automatic or semi-automatic way.
            Of note, when you use manual way, you need to run 2 times, 
            firstly you should set --manual to get signature estimation results,
            and secondly you should set --manual --number N to get N signatures.
  ==
  fit     - fit signatures in >=1 samples based on COSMIC reference signatures.
  ==
  bt      - run bootstrap signature fitting analysis in >=1 samples.
  ==
  show    - show some useful information or plots. See README for details.

Usage:
  sigflow extract --input=<file> [--output=<outdir>] [--mode=<class>] [--manual --number <sigs>] [--max <max>] [--genome=<genome>] [--nrun=<runs>] [--cores=<cores>] [--sigprofiler] [--refit] [--hyper] [--verbose]
  sigflow fit --input=<file> [--output=<outdir>] [--index=<index>] [--mode=<class>] [--genome=<genome>] [--verbose]
  sigflow bt  --input=<file> [--output=<outdir>] [--index=<index>] [--mode=<class>] [--genome=<genome>] [--nrun=<runs>] [--verbose]
  sigflow show [--isearch=<keyword>] [--index=<index> --mode=<class>] [--output=<outdir>] [--verbose]
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
  --refit                         refit the de-novo signatures with quadratic programming or nnls (SigProfiler).
  --hyper                         enable hyper mutation handling in COSMIC signatures (not used by SigProfiler approach).
  --sigprofiler                   enable automatic extraction by SigProfiler software.
  --isearch <keyword>             search and how cancer type specific reference signature index by keyword, e.g. breast.
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

Example datasets along with many example code are available in clone repository above (you can read it online at [here](https://github.com/ShixiangWang/sigflow/tree/master/test)).

The following parts give an example for each command.

Result directory of any command has the following structure.

- Files with extension`.RData` and `.rds`  are R related files to reproduce the results, and can be imported into R for further analysis and visualization.
- Files with extension`.pdf` are common visualization results used for communication.
- Files with extension `.csv` are formated data tables used for inspection, communication or further analysis.

#### `extract` command

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigflow
$ cd sigflow/test
$ sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf -m MAF -r 10 -T 4 --max 10
```

This will auto-extract SBS/DBS/INDEL signatures from data `toga_laml.maf.gz` (a gzipped MAF file from [Maftools](https://github.com/PoisonAlien/maftools)) by 10 Bayesian NMF runs with 4 computer threads and output results to directory `test_results/test_maf`. At default, **Bayesian NMF** approach is used, it starts from 10 signatures (set by `--max`) and reduces to a optimal signature number. If `--sigprofiler` is enabled, i.e.

```sh
$ sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf -m MAF -r 10 -T 4 --max 10 --sigprofiler
```

Sigflow will use the [SigProfiler](https://github.com/AlexandrovLab/SigProfilerExtractor) to auto-extract signatures, here it will extract 2 to 10 signatures and determine the optimal solution.

> NOTE, **in practice, set `-r` to a value  `>=10` is recommended** for auto-extraction with Bayesian NMF, `>=100` for semi-automatic extraction with basic NMF and automatic extraction with SigProfiler (enabled by `--sigprofiler`).

Results of `extract` command have the following structure:

![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909141307.png)

> Here, no DBS records found in input data, so no corresponding result files exist.

- Tally: mutation catalogue data and plots of all samples or individual samples are stored in files/directory contains `tally`.

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144059.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909143533.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909143624.png)

- Signature: signature profile (relative contribution in each signature) data and plots of all samples are stored in files contains `signature`.

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144133.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909143815.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909143950.png)

- Exposure: exposure profile (relative and absolute contribution of each signature to each sample) data and plots of all samples are stored in files contains `exposure`.

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144214.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144310.png)

  ![](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144402.png)

- Similarity: to understand the etiologies of extracted signatures, cosine similarity analysis is applied to extracted signatures and reference signatures from COSMIC database. The result files contains `similarity` and `best_match`.

  ![image-20200909144940222](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909144940.png)

  ![image-20200909145001347](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145001.png)

  ![image-20200909145039312](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145039.png)

  ![image-20200909145114261](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145114.png)

  ![image-20200909145214337](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145214.png)

- Clustering: the samples can be clustered based on signature relative exposure. This analysis is done by kmeans and the results are outputed (Note, the cluster number is same as signature number).

  ![image-20200909145522952](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145523.png)

  ![image-20200909145546261](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145546.png)

  ![image-20200909145649068](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909145649.png)

  

#### `fit` command

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigflow
$ cd sigflow/test
$ sigflow fit -i tcga_laml.maf.gz -o test_results/test_fitting -m MAF
```

This will auto-fit input data `tcga_laml.maf.gz` to COSMIC SBS/DBS/INDEL signatures. Signature exposure data tables and plots are outputed.

Results of `fit` command have the following structure:

![image-20200909150037803](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150037.png)

> Here, no DBS records found in input data, so no corresponding result files exist.
>
> `legacy` represents COSMIC v2 SBS signatures and `SBS` represents COSMIC v3 SBS signatures.

- Tally: same as results from `extract` command.

- Fitting: fitted relative/absolute signature exposure, reconstructed error (calculated by Frobenius norm) data and corresponding plots of all samples are stored in files contains `fitting`.

  ![image-20200909150444291](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150444.png)

  ![image-20200909151002018](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909151002.png)

  > Relative signature exposure in each sample.

  ![image-20200909150607896](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150608.png)

  > Absolute signature exposure in each sample.

  ![image-20200909150902466](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150902.png)

  > Visualization of relative signature exposure.

  ![image-20200909150743352](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150743.png)

  > Visualization of absolute signature exposure.

  

  ![image-20200909150646872](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909150646.png)

  > Reconstructed error for each sample.

  

#### `bt` command

Bootstrapping analysis was performed to evaluate the stability of signature exposure. For a tumor, this analysis firstly resamples mutations based on the observed mutation type (component) frequency and then applies signature fitting to the bootstrapping samples. For example, if a tumor harbors 100 mutations and 20 mutations are classified into T[C>T]T, then we resample 100 mutations and the probability to assign these mutation to T[C>T]T is 0.2. If we repeat such process many times, we can estimate the confidence interval of exposure of a signature in this tumor.

> More details please read paper [Detecting presence of mutational signatures in cancer with confidence](https://academic.oup.com/bioinformatics/article/34/2/330/4209996).

```sh
$ # Assume you have done the clone step
$ # git clone https://github.com/ShixiangWang/sigflow
$ cd sigflow/test
$ sigflow bt -i tcga_laml.maf.gz -o test_results/test_bt -m SBS -r 5
```

This will  resample mutation catalogue of each sample based on observed mutation type frequency and run signature fitting using COSMIC SBS/DBS/INDEL signatures. The process is repeated multiple times and controlled by option `-r` (here is 5). This bootstrap analysis is used to estimate the instability of signature exposure. Data tables and plots of bootstrap signature exposures, errors and p values under different exposure cutoff are outputed.

> NOTE, in practice, set `-r` to a value  `>=100` is recommended.

Results of `bt` command have the following structure:

![image-20200909152336875](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909152337.png)

- Tally: same as results from `extract` command.

- Bootstrap fitting: fitted relative/absolute signature exposure, reconstructed error (calculated by Frobenius norm of residue), signature instability data and corresponding plots of all bootstrapping samples and individual samples are stored in files/directory contains `bootstrap`.

  ![image-20200909153710051](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909153710.png)

  ![image-20200909153800225](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909153800.png)

  > Bootstrap signature exposure distribution.

  ![image-20200909154812780](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909154812.png)

  > For each sample, the distribution of bootstrap exposures is plotted as boxplot and the fitting result with original input data is labelled by triangle. 

  ![image-20200909154022111](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909154022.png)

  > Note, the sample data without bootstrapping process are also fitted and labelled as `type = "optimal"`

  ![image-20200909153952230](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909153952.png)

  ![image-20200909154325658](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909154325.png)

  > The p values are calculated as the proportion of how many bootstrapping samples have exposures under specified exposure cutoff.

  ![image-20200909154555142](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200909154555.png)

  > Signature exposure instability is measured as MRSE between exposures in bootstrap samples and exposures in the original samples for each tumor/patient.




#### `show` command

`show` command provides extra information to help user analyze signatures. This includes:

1. Search cancer type specific signature indices, this may help user to set the reference indices in `fit` and `bt` commands. This information could read [online](https://shixiangwang.github.io/sigminer-doc/sigflow.html#cancer-type-specific-signature-index-database).
2. Generate COSMIC reference signature profiles.

For the no.1 task, one could run

```sh
$ sigflow show --isearch breast
```

This will generate the following output:

![image-20200918150257473](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200918150257.png)

For the no.2 task, one could run

```sh
sigflow show --mode SBS --index 1,2,3,7a -o test_show_sig_profile
```

This will generate signature profile for signature 1,2,3,7. For SBS, two versions of plots exist.

![image-20200918150649303](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200918150649.png)

COSMIC v2:

![image-20200918150754770](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200918150754.png)

COSMIC v3:

![image-20200918150825289](https://gitee.com/ShixiangWang/ImageCollection/raw/master/png/20200918150825.png)

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

- 2024-04-01: upgraded sigminer to v2.3.
- 2023-06-08: add thread controls in boostrap fit instead of using all available cores as workers.
- 2022-02-13: add troubleshooting section.
- 2021-06-17: upgraded sigminer to v2.0.2 to fix [issue #31](https://github.com/ShixiangWang/sigflow/issues/31).
- 2021-04-02: upgraded sigminer to v2.0.
- 2021-03-29: 
  - enhanced data check.
  - upgraded sigminer to v1.9 (the alpha version for sigminer v2).
- 2020-12-01: Another typo for hg38 reference genome build, thanks to @PalashPandey.
- 2020-11-30: Fix a typo which incorrectly assign the reference genome, thanks to @PalashPandey.
- 2020-10-19: Update citation paper.
- 2020-09-28: Copy signature estimation output pdf file from SigProfiler to result directory, **release Sigflow 1.3**.
- 2020-09-18: Update doc for `bt` command and add doc for `show` command.
- 2020-09-15: Release **Sigflow 1.2**
  - Upgrade sigminer version (v1.0.18) to fix bug when outputing results for only one signatures ([#17](https://github.com/ShixiangWang/sigflow/issues/17)).
- 2020-09-15: Prepare and try to release **Sigflow 1.1**. This is a version with its co-evolutionary R package Sigminer v1.0.17 and gold-standard de-novo signature extraction tool SigProfilerExtractor v1.0.17 as backends.
- 2020-09-14:
  - Add new command `show` to search cancer type specific reference signature indices and plot COSMIC signatures.
  - Support `--refit` in SigProfiler calling.
  - Fix a bug in identifying COSMIC v2 indices in signature fitting.
  - Upgrade sigminer version to v1.0.17.
- 2020-09-09:
  - Update README and documentation of input and usage.
- 2020-09-03: 
  - Use sigminer v1.0.15 and support inputing reference signature index in `fit` and `bt` commands.
  - Allow users to decide if refit the signature exposures after *de novo* extraction with `refit` option.
  - Support matrix as input.
- 2020-08-14: Use sigminer v1.0.11 to use SigProfilerExtractor v1.0.15 to avoid issue raised from SigProfilerExtractor updates.
- 2020-08-05: **Release Sigflow 1.0** and related Docker image. This version is based on Sigminer v1.0.10, R v4.0.2 and SigProfilerExtractor v.1.0.15.
  - Support SigProfiler.
  - Add `verbose` option.
  - Add `max` option.
  - Add `hyper` option.
  - More flexible and reasonable configuration.
- 2020-07-29: **Release Sigflow 0.1** using Docker. Sigflow 0.1 is based on Sigminer v1.0.9 and R v4.0.2



## Test

There are some example data sets in this repository, you can find how to test different workflows in [test/test.sh](test/test.sh).
It is time consuming to run all tests, just pick an example test similar to your task and see how it works. Before releasing a new version of **Sigflow**, I would run all tests to make sure they work well.

## Troubleshooting

1. Error like `the supplied end is > refwidth`. ([[#32](https://github.com/ShixiangWang/sigflow/issues/32)])

The reference genome for variant calling is not (perfectly) match the specified genome in `sig_tally()`. If you make sure the reference genome is correct, please try finding the variant records with uncompatible position and removing them before rerun.


## Other tools

- For interactive analysis and visualization, please refer to its co-evolutionary R package [sigminer](https://github.com/ShixiangWang/sigminer).
- For mutation analysis, please refer to [Maftools](https://github.com/PoisonAlien/maftools).

## Citation

If you are using **Sigflow** fow now in academic field, please cite:

Shixiang Wang, Ziyu Tao, Tao Wu, Xue-Song Liu, Sigflow: An Automated And Comprehensive Pipeline For Cancer Genome Mutational Signature Analysis, Bioinformatics, btaa895, https://doi.org/10.1093/bioinformatics/btaa895

## License

This software is released under [Academic Free License ("AFL") v.3.0](https://opensource.org/licenses/AFL-3.0)

Copyright 2020 © Shixiang Wang, Xue-Song Liu

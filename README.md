# sigflow: Streamline Analysis Workflow for Mutational Signatures

```
       .__        _____.__                 
  _____|__| _____/ ____\  |   ______  _  __
 /  ___/  |/ ___\   __\|  |  /  _ \ \/ \/ /
 \___ \|  / /_/  >  |  |  |_(  <_> )     / 
/____  >__\___  /|__|  |____/\____/ \/\_/  
     \/  /_____/  
```

**sigflow** provides mutational signature analysis workflows based on R package [sigminer](https://github.com/ShixiangWang/sigminer).


## Installation

Requirements:

- R (>3.6).
- R package [sigminer](https://github.com/ShixiangWang/sigminer).
- R package [docopt]().

Steps:

1. Install R - follow the instructions at <https://cran.r-project.org/>.
2. Install R packages, run

```r
install.packages("docopt")
install.packages("BiocManager")
BiocManager::install("sigminer", dependencies = TRUE)
```

3. Clone this repository, run

```bash
git clone https://github.com/ShixiangWang/sigminer.workflow
```

## Usage

```
sigflow: Streamline Analysis Workflow for Mutational Signatures.

Author: Shixiang Wang (wangshx@shanghaitech.edu.cn)
Copyright: MIT@2020

Desc:
  There are several subcommands.
  ==
  extract - extract signatures by either automated or manual way.
            Of note, when you use manual way, you need to run 2 times, 
            firstly you should set --manual to get signature estimation results,
            and secondly you should set --manual --number N to get N signatures.
  ==
  fit     - fit signatures in >=1 samples.

Usage:
  sigflow.R extract --input=<file> [--output=<outdir>] [--mode=<class>] [--manual --number <sigs>] [--genome=<genome>] [--nrun=<runs>] [--cores=<cores>]
  sigflow.R (-h | --help)
  sigflow.R --version

Options:
  -h --help     Show help message.
  --version     Show version.
  -i <file>, --input <file>       input file/directory path.
  -o <outdir>, --output <outdir>  output directory path [default: ./sigflow_result/].
  -m <class>, --mode <class>      extract mode, can be one of SBS, DBS, ID, MAF (for three types), CN [default: SBS].
  --manual                        enable manual extraction, set -N=0 for outputing signature estimation firstly.
  -N <sigs>, --number <sigs>      extract specified number of signatures [default: 0].
  -g <genome>, --genome <genome>  genome build, can be hg19, hg38 or mm10, [default: hg19].
  -r <runs>, --nrun <runs>        run times of NMF to get results [default: 30].
  -T <cores>, --cores <cores>     cores to run the program, large dataset will benefit from it [default: 1].

```

## License

This software is release under [Academic Free License ("AFL") v. 3.0](https://opensource.org/licenses/AFL-3.0)

Copyright 2020 © Shixiang Wang

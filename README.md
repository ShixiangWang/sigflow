# sigflow: Streamline Analysis Workflow for Mutational Signatures

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

```bash
$ Rscript sigflow.R -h
```

## Copyright

Copyright 2020 © Shixiang Wang. All rights reserved. Without the permission of the owner, do not copy or modify, especially for commercial use.

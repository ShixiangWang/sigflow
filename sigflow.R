#!/usr/bin/env Rscript

"=================================================================
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
  sigflow.R extract --input=<file> [--output=<outdir>] [--mode=<class>] [--manual --number <sigs>][--format=<format>] [--genome=<genome>] [--nrun=<runs>] [--cores=<cores>]
  sigflow.R (-h | --help)
  sigflow.R --version

Options:
  -h --help     Show help message.
  --version     Show version.
  -i <file>, --input <file>       input file path.
  -o <outdir>, --output <outdir>  output directory path [default: ./sigflow_result/].
  -m <class>, --mode <class>      extract mode, can be one of SBS, DBS, ID, MAF (for three types), CN [default: SBS].
  --manual                        enable manual extraction, set -N=0 for outputing signature estimation firstly.
  -N <sigs>, --number <sigs>      extract specified number of signatures [default: 0].
  -g <genome>, --genome <genome>  human genome build, can be hg19 or hg38, [default: hg19].
  -r <runs>, --nrun <runs>        times of NMF to get results [default: 30].
  -T <cores>, --cores <cores>     cores to run the program, large dataset will benefit from it [default: 1].

=================================================================
" -> doc

library(docopt)
arguments <- docopt(doc, version = "sigflow v0.1\n")

message("
=================================================================

       .__        _____.__                 
  _____|__| _____/ ____\\  |   ______  _  __
 /  ___/  |/ ___\\   __\\|  |  /  _ \\ \\/ \\/ /
 \\___ \\|  / /_/  >  |  |  |_(  <_> )     / 
/____  >__\\___  /|__|  |____/\\____/ \\/\\_/  
     \\/  /_____/                           

Name   :       sigflow
Author :       Shixiang Wang
Version:        0.1
License:        MIT
Link   :        https://github.com/ShixiangWang/sigminer.workflow
Doc    :        https://shixiangwang.github.io/sigminer-doc/
============================== START
")


# Parsing input -----------------------------------------------------------

message("Parsing parameters...\n------")
# print(arguments)
ARGS <- arguments[!startsWith(names(arguments), "--")]
for (i in seq_along(ARGS)) {
  message(names(ARGS)[i], "\t\t: ", ARGS[i])
}
message("------\n")

# Function part -----------------------------------------------------------
## Program to go

flow_extraction <- function(obj, genome_build, mode, manual_step, nrun, cores, result_dir) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
    sigminer:::send_success("Result directory ", result_dir, " created.")
  } else {
    sigminer:::send_info("Result directory ", result_dir, " existed.")
  }

  timer <- Sys.time()
  sigminer:::send_info("Started.")
  on.exit(sigminer:::send_elapsed_time(timer))

  if (genome_build == "hg19") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  } else if (genome_build == "hg38") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  }

  if (!require(ref_genome, character.only = TRUE)) {
    sigminer:::send_info("Package ", ref_genome, " not found, try installing.")
    BiocManager::install(ref_genome)
  }

  if (mode != "CN") {
    if (mode == "MAF") {
      mode <- "ALL"
    }
    tally_list <- sig_tally(obj, mode = "ALL", ref_genome = ref_genome)
  } else {
    tally_list <- sig_tally(obj, cores = cores)
  }

  if (manual_step < 0) {
    ## auto-extract
    if (mode == "CN") {
      sigs_CN <- sig_auto_extract(
        nmf_matrix = tally_list$nmf_matrix,
        result_prefix = "BayesianNMF_CN",
        destdir = file.path(result_dir, "BayesianNMF"),
        K0 = ncol(tally_list$nmf_matrix) - 1L,
        nrun = nrun,
        cores = cores,
        skip = TRUE
      )
    } else {
      if (mode == "ALL" | mode == "SBS") {
        mat <- tally_list$SBS_96
        if (!is.null(mat)) {
          sigs_SBS <- sig_auto_extract(
            nmf_matrix = mat,
            result_prefix = "BayesianNMF_SBS",
            destdir = file.path(result_dir, "BayesianNMF"),
            K0 = ncol(mat) - 1L,
            nrun = nrun,
            cores = cores,
            skip = TRUE
          )
        } else {
          sigminer:::send_info("SBS_96 matrix is NULL, skip extracting.")
        }
      } else if (mode == "ALL" | mode == "DBS") {
        mat <- tally_list$DBS_78
        if (!is.null(mat)) {
          sigs_DBS <- sig_auto_extract(
            nmf_matrix = mat,
            result_prefix = "BayesianNMF_DBS",
            destdir = file.path(result_dir, "BayesianNMF"),
            K0 = ncol(mat) - 1L,
            nrun = nrun,
            cores = cores,
            skip = TRUE
          )
        } else {
          sigminer:::send_info("DBS_78 matrix is NULL, skip extracting.")
        }
      } else if (mode == "ALL" | mode == "ID") {
        mat <- tally_list$ID_83
        if (!is.null(mat)) {
          sigs_ID <- sig_auto_extract(
            nmf_matrix = mat,
            result_prefix = "BayesianNMF_ID",
            destdir = file.path(result_dir, "BayesianNMF"),
            K0 = ncol(mat) - 1L,
            nrun = nrun,
            cores = cores,
            skip = TRUE
          )
        } else {
          sigminer:::send_info("ID_83 matrix is NULL, skip extracting.")
        }
      }
    }
  } else {
    ## manual-extract
  }
  
  ## Outputing to result directory
  if (exists("sigs_CN")) {
    
  }
  if (exists("sigs_SBS")) {
    
  }
  if (exists("sigs_DBS")) {
    
  }
  if (exists("sigs_ID")) {
    
  }
  
}

# Action part -------------------------------------------------------------
## Check inputs and call the working functions
suppressPackageStartupMessages(library("sigminer"))

message("Reading file...\n------")

mode <- ARGS$mode
input <- ARGS$input
file_format <- ARGS$format
genome_build <- ARGS$genome

if (any(endsWith(input, c("xls", "xlsx")))) {
  if (!require("readxl")) {
    message("readxl not found, try installing.")
    install.packages("readxl")
    library("readxl")
  }
  input <- readxl::read_excel(input)
}


if (mode == "CN") {
  isCN <- TRUE
} else {
  isCN <- FALSE
}

if (!isCN) {
  obj <- sigminer::read_maf(input, verbose = TRUE)
} else {
  obj <- sigminer::read_copynumber(input, genome_build = genome_build, verbose = TRUE)
}

message("------\n")

if (ARGS$extract) {
  message("Running signature extraction pipeline...\n------")
  if (!ARGS$manual) {
    manual_step <- -1L
  } else {
    manual_step <- as.integer(ARGS$number)
  }
  nrun <- as.integer(ARGS$nrun)
  cores <- min(as.integer(ARGS$cores), parallel::detectCores())
  flow_extraction(
    obj = obj, genome_build = genome_build, mode = ARGS$mode,
    manual_step = manual_step, nrun = nrun, cores = cores,
    result_dir = ARGS$output
  )
}


# End part ----------------------------------------------------------------

message("
============================== END

TODO: Some End messages to be added...

=================================================================
")

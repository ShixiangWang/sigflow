#!/usr/bin/env Rscript

"=================================================================
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
  show    - show some useful information.

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
  --sigprofiler                   enable automatic extraction by SigProfiler software.
  -v, --verbose                   print verbose message.

=================================================================
" -> doc

if (!suppressMessages(require("docopt"))) {
  install.packages("docopt", repos = "https://cloud.r-project.org")
}

library("docopt")
arguments <- docopt(doc, version = "sigflow v1.1\n")

## Stop error parsing
if (!exists("arguments")) {
  quit("no", status = -1)
}

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
Version:       1.1
License:       AFL 3.0
Link   :       https://github.com/ShixiangWang/sigminer.workflow
Doc    :       https://shixiangwang.github.io/sigminer-doc/
============================== START
")


# Parsing input -----------------------------------------------------------

message("Parsing parameters...\n------")
ARGS <- arguments[!startsWith(names(arguments), "--")]
for (i in seq_along(ARGS)) {
  if (names(ARGS)[i] %in% "sigprofiler") {
    message(names(ARGS)[i], "\t: ", ARGS[i])
  } else {
    message(names(ARGS)[i], "\t\t: ", ARGS[i])
  }
}

message("------\n")

# Function part -----------------------------------------------------------
## Program to go

flow_extraction <- function(obj, genome_build, mode, manual_step, nrun, cores, result_dir,
                            max_number = 100, 
                            refit = TRUE,
                            rm_hyper = FALSE, sigprofiler = FALSE) {
  if (!dir.exists(file.path(result_dir, "results"))) {
    dir.create(file.path(result_dir, "results"), recursive = TRUE)
  }

  timer <- Sys.time()
  sigminer:::send_info("Pipeline for extraction started.")
  on.exit(sigminer:::send_elapsed_time(timer))

  if (genome_build == "hg19") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  } else if (genome_build == "hg38") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  } else if (genome_build == "mm10") {
    ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
  }

  if (!suppressMessages(require(ref_genome, character.only = TRUE))) {
    sigminer:::send_info("Package ", ref_genome, " not found, try installing.")
    BiocManager::install(ref_genome)
  }

  if (mode != "CN") {
    if (mode == "MAF") {
      mode <- "ALL"
    }
    if (manual_step <= 0) {
      if (!is.matrix(obj)) {
        tally_list <- sig_tally(obj, mode = "ALL", ref_genome = ref_genome)
      } else {
       ## Construct a fake tally result
       ## In this situation, mode cannot be ALL
       mm <- switch(mode,
                    SBS = "SBS_96",
                    DBS = "DBS_78",
                    ID = "ID_83")
       tally_list <- list()
       tally_list[[mm]] <- obj %>% t()
      }
      save(tally_list, file = file.path(result_dir, "maf_tally.RData"))
      if (!is.null(tally_list$SBS_96) & mode %in% c("SBS", "ALL")) {
        output_tally(tally_list$SBS_96 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "SBS")
      }
      if (!is.null(tally_list$DBS_78) & mode %in% c("DBS", "ALL")) {
        output_tally(tally_list$DBS_78 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "DBS")
      }
      if (!is.null(tally_list$ID_83) & mode %in% c("ID", "ALL")) {
        output_tally(tally_list$ID_83 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "ID")
      }
    }
  } else {
    if (manual_step <= 0) {
      if (!is.matrix(obj)) {
        tally_list <- sig_tally(obj, ignore_chrs = c("chrX", "chrY"), cores = cores)
      } else {
        tally_list <- list()
        tally_list[["nmf_matrix"]] <- obj %>% t()
      }
      save(tally_list, file = file.path(result_dir, "cn_tally.RData"))
      output_tally(tally_list, result_dir = file.path(result_dir, "results"), mut_type = "CN")
    }
  }

  if (manual_step < 0) {
    ## auto-extract
    if (mode == "CN") {
      if (sigprofiler) {
        sigprofiler_extract(
          nmf_matrix = tally_list$nmf_matrix,
          output = file.path(result_dir, "SigProfiler_CN"),
          range = 2:min(30, nrow(tally_list$nmf_matrix) - 1, max_number),
          nrun = nrun,
          refit = refit,
          cores = cores,
          use_conda = TRUE
        )
        sigs_CN <- sigprofiler_import(file.path(result_dir, "SigProfiler_CN"))
        data.table::fwrite(sigs_CN$all_stats, file = file.path(result_dir, "results", "SigProfiler_CN_stats.csv"))
        sigs_CN <- sigs_CN$solution
      } else {
        sigminer:::send_info("Auto extract copy number signatures.")
        sigs_CN <- sig_auto_extract(
          nmf_matrix = tally_list$nmf_matrix,
          result_prefix = "BayesianNMF_CN",
          destdir = file.path(result_dir, "BayesianNMF"),
          strategy = "stable",
          K0 = min(30, nrow(tally_list$nmf_matrix) - 1, max_number),
          nrun = nrun,
          cores = cores,
          optimize = refit,
          skip = TRUE
        )
      }
    } else {
      if (mode == "ALL" | mode == "SBS") {
        mat <- tally_list$SBS_96
        if (!is.null(mat)) {
          sigminer:::send_info("Auto extract SBS signatures.")
          if (sigprofiler) {
            sigprofiler_extract(
              nmf_matrix = mat,
              output = file.path(result_dir, "SigProfiler_SBS"),
              range = 2:min(30, nrow(mat) - 1, max_number),
              nrun = nrun,
              refit = refit,
              cores = cores,
              use_conda = TRUE
            )
            sigs_SBS <- sigprofiler_import(file.path(result_dir, "SigProfiler_SBS"))
            data.table::fwrite(sigs_SBS$all_stats, file = file.path(result_dir, "results", "SigProfiler_SBS_stats.csv"))
            sigs_SBS <- sigs_SBS$solution
          } else {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_SBS <- sig_auto_extract(
              nmf_matrix = mat,
              result_prefix = "BayesianNMF_SBS",
              destdir = file.path(result_dir, "BayesianNMF"),
              strategy = "stable",
              K0 = min(30, nrow(mat) - 1, max_number),
              nrun = nrun,
              cores = cores,
              optimize = refit,
              skip = TRUE
            )
          }
        } else {
          sigminer:::send_info("SBS_96 matrix is NULL, skip extracting.")
        }
      }
      if (mode == "ALL" | mode == "DBS") {
        mat <- tally_list$DBS_78
        if (!is.null(mat)) {
          sigminer:::send_info("Auto extract DBS signatures.")
          if (sigprofiler) {
            sigprofiler_extract(
              nmf_matrix = mat,
              output = file.path(result_dir, "SigProfiler_DBS"),
              range = 2:min(15, nrow(mat) - 1, max_number),
              nrun = nrun,
              refit = refit,
              cores = cores,
              use_conda = TRUE
            )
            sigs_DBS <- sigprofiler_import(file.path(result_dir, "SigProfiler_DBS"))
            data.table::fwrite(sigs_DBS$all_stats, file = file.path(result_dir, "results", "SigProfiler_DBS_stats.csv"))
            sigs_DBS <- sigs_DBS$solution
          } else {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_DBS <- sig_auto_extract(
              nmf_matrix = mat,
              result_prefix = "BayesianNMF_DBS",
              destdir = file.path(result_dir, "BayesianNMF"),
              strategy = "stable",
              K0 = min(15, nrow(mat) - 1, max_number),
              nrun = nrun,
              cores = cores,
              optimize = refit,
              skip = TRUE
            )
          }
        } else {
          sigminer:::send_info("DBS_78 matrix is NULL, skip extracting.")
        }
      }
      if (mode == "ALL" | mode == "ID") {
        mat <- tally_list$ID_83
        if (!is.null(mat)) {
          sigminer:::send_info("Auto extract ID signatures.")
          if (sigprofiler) {
            sigprofiler_extract(
              nmf_matrix = mat,
              output = file.path(result_dir, "SigProfiler_ID"),
              range = 2:min(20, nrow(mat) - 1, max_number),
              nrun = nrun,
              refit = refit,
              cores = cores,
              use_conda = TRUE
            )
            sigs_ID <- sigprofiler_import(file.path(result_dir, "SigProfiler_ID"))
            data.table::fwrite(sigs_ID$all_stats, file = file.path(result_dir, "results", "SigProfiler_ID_stats.csv"))
            sigs_ID <- sigs_ID$solution
          } else {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_ID <- sig_auto_extract(
              nmf_matrix = mat,
              result_prefix = "BayesianNMF_ID",
              destdir = file.path(result_dir, "BayesianNMF"),
              strategy = "stable",
              K0 = min(20, nrow(mat) - 1, max_number),
              nrun = nrun,
              cores = cores,
              optimize = refit,
              skip = TRUE
            )
          }
        } else {
          sigminer:::send_info("ID_83 matrix is NULL, skip extracting.")
        }
      }
    }
  } else {
    ## manual-extract
    if (manual_step == 0) {
      suppressPackageStartupMessages(library("NMF"))
      ## do signature estimation
      if (mode == "CN") {
        est_CN <- sig_estimate(
          nmf_matrix = tally_list$nmf_matrix,
          nrun = nrun,
          range = 2:min(30, max_number),
          cores = cores,
          pConstant = 1e-9,
          save_plots = TRUE,
          plot_basename = file.path(result_dir, "manual_extraction", "CN_sig_number"),
          verbose = TRUE
        )

        data.table::fwrite(est_CN$survey, file = file.path(result_dir, "manual_extraction", "CN_sig_number_survey.csv"))
        p <- show_sig_number_survey(est_CN$survey)
        ggsave(file.path(result_dir, "manual_extraction", "CN_sig_number_survey_simple.pdf"),
          plot = p,
          width = 10, height = 6
        )
      } else {
        if (mode == "ALL" | mode == "SBS") {
          mat <- tally_list$SBS_96

          if (!is.null(mat)) {
            est_SBS <- sig_estimate(
              nmf_matrix = mat,
              nrun = nrun,
              range = 2:min(30, max_number),
              cores = cores,
              pConstant = 1e-9,
              save_plots = TRUE,
              plot_basename = file.path(result_dir, "manual_extraction", "SBS_sig_number"),
              verbose = TRUE
            )

            data.table::fwrite(est_SBS$survey, file = file.path(result_dir, "manual_extraction", "SBS_sig_number_survey.csv"))
            p <- show_sig_number_survey(est_SBS$survey)
            ggsave(file.path(result_dir, "manual_extraction", "SBS_sig_number_survey_simple.pdf"),
              plot = p,
              width = 10, height = 6
            )
          }
        }
        if (mode == "ALL" | mode == "DBS") {
          mat <- tally_list$DBS_78

          if (!is.null(mat)) {
            est_DBS <- sig_estimate(
              nmf_matrix = mat,
              nrun = nrun,
              range = 2:min(15, max_number),
              cores = cores,
              pConstant = 1e-9,
              save_plots = TRUE,
              plot_basename = file.path(result_dir, "manual_extraction", "DBS_sig_number"),
              verbose = TRUE
            )

            data.table::fwrite(est_DBS$survey, file = file.path(result_dir, "manual_extraction", "DBS_sig_number_survey.csv"))
            p <- show_sig_number_survey(est_DBS$survey)
            ggsave(file.path(result_dir, "manual_extraction", "DBS_sig_number_survey_simple.pdf"),
              plot = p,
              width = 10, height = 6
            )
          }
        }
        if (mode == "ALL" | mode == "ID") {
          mat <- tally_list$ID_83

          if (!is.null(mat)) {
            est_ID <- sig_estimate(
              nmf_matrix = mat,
              nrun = nrun,
              range = 2:min(20, max_number),
              cores = cores,
              pConstant = 1e-9,
              save_plots = TRUE,
              plot_basename = file.path(result_dir, "manual_extraction", "ID_sig_number"),
              verbose = TRUE
            )

            data.table::fwrite(est_ID$survey, file = file.path(result_dir, "manual_extraction", "ID_sig_number_survey.csv"))
            p <- show_sig_number_survey(est_ID$survey)
            ggsave(file.path(result_dir, "manual_extraction", "ID_sig_number_survey_simple.pdf"),
              plot = p,
              width = 10, height = 6
            )
          }
        }
      }

      message("==============================")
      message(
        "Checking plots in ",
        file.path(result_dir, "manual_extraction"),
        " to choose a proper signature number for running the next step."
      )
      message("NOTE: if you run all mutation types in this steps, you have to extract them one by one.")
      message("==============================")
    } else {
      suppressPackageStartupMessages(library("NMF"))
      ## extract specified signatures
      if (mode == "CN") {
        load(file = file.path(result_dir, "cn_tally.RData"))

        sigs_CN <- sig_extract(
          nmf_matrix = tally_list$nmf_matrix,
          n_sig = manual_step,
          nrun = nrun,
          cores = cores,
          optimize = refit,
          pConstant = 1e-9
        )
      } else {
        load(file = file.path(result_dir, "maf_tally.RData"))

        if (mode == "ALL") {
          sigminer:::send_warning("In the step 2 of manual extraction, you can only extract one type of signature (SBS/DBS/ID).")
          sigminer:::send_warning("Set mode to 'SBS'.")
          mode <- "SBS"
        }

        if (mode == "SBS") {
          mat <- tally_list$SBS_96

          if (!is.null(mat)) {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_SBS <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = refit,
              pConstant = 1e-9
            )
          }
        }
        if (mode == "DBS") {
          mat <- tally_list$DBS_78

          if (!is.null(mat)) {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_DBS <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = refit,
              pConstant = 1e-9
            )
          }
        }
        if (mode == "ID") {
          mat <- tally_list$ID_83

          if (!is.null(mat)) {
            if (rm_hyper) {
              mat <- sigminer::handle_hyper_mutation(mat)
            }
            sigs_ID <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = refit,
              pConstant = 1e-9
            )
          }
        }
      }
    }
  }

  ## Outputing to result sub-directories: CN, SBS, DBS, ID, etc.
  if (exists("sigs_CN")) {
    output_sig(sigs_CN, result_dir = file.path(result_dir, "results"), mut_type = "CN")
  }
  if (exists("sigs_SBS")) {
    output_sig(sigs_SBS, result_dir = file.path(result_dir, "results"), mut_type = "SBS")
  }
  if (exists("sigs_DBS")) {
    output_sig(sigs_DBS, result_dir = file.path(result_dir, "results"), mut_type = "DBS")
  }
  if (exists("sigs_ID")) {
    output_sig(sigs_ID, result_dir = file.path(result_dir, "results"), mut_type = "ID")
  }

  return(invisible(NULL))
}


flow_fitting <- function(obj, genome_build, mode, result_dir, nrun = NULL,
                         index = "ALL",
                         prog = c("fit", "bootstrap")) {
  prog <- match.arg(prog)

  if (!dir.exists(file.path(result_dir, "results"))) {
    dir.create(file.path(result_dir, "results"), recursive = TRUE)
  }

  timer <- Sys.time()
  sigminer:::send_info("Pipeline for (bootstrap) fitting started.")
  on.exit(sigminer:::send_elapsed_time(timer))

  if (genome_build == "hg19") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  } else if (genome_build == "hg38") {
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  } else if (genome_build == "mm10") {
    ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
  }

  if (!suppressMessages(require(ref_genome, character.only = TRUE))) {
    sigminer:::send_info("Package ", ref_genome, " not found, try installing.")
    BiocManager::install(ref_genome)
  }

  if (mode != "CN") {
    if (mode == "MAF") {
      mode <- "ALL"
    }
    tally_list <- sig_tally(obj, mode = "ALL", ref_genome = ref_genome)
    save(tally_list, file = file.path(result_dir, "maf_tally.RData"))
    if (!is.null(tally_list$SBS_96) & mode %in% c("SBS", "ALL")) {
      output_tally(tally_list$SBS_96 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "SBS")
    }
    if (!is.null(tally_list$DBS_78) & mode %in% c("DBS", "ALL")) {
      output_tally(tally_list$DBS_78 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "DBS")
    }
    if (!is.null(tally_list$ID_83) & mode %in% c("ID", "ALL")) {
      output_tally(tally_list$ID_83 %>% t(), result_dir = file.path(result_dir, "results"), mut_type = "ID")
    }
  } else {
    stop("Fitting reference signatures for CN is not supported for now!")
    # tally_list <- sig_tally(obj, ignore_chrs = c("chrX", "chrY"), cores = parallel::detectCores())
    # save(tally_list, file = file.path(result_dir, "cn_tally.RData"))
    # output_tally(tally_list, result_dir = file.path(result_dir, "results"), mut_type = "CN")
  }

  if (prog == "fit") {
    if (mode == "ALL" | mode == "SBS") {
      mat <- tally_list$SBS_96

      if (!is.null(mat)) {
        mat <- t(mat)
        ## COSMIC V2 SBS
        fit_legacy <- sig_fit(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "legacy",
          mode = "SBS",
          return_class = "data.table",
          return_error = TRUE
        )
        output_fit(fit_legacy, result_dir = file.path(result_dir, "results"), mut_type = "legacy")

        fit_SBS <- sig_fit(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "SBS",
          mode = "SBS",
          return_class = "data.table",
          return_error = TRUE
        )
        output_fit(fit_SBS, result_dir = file.path(result_dir, "results"), mut_type = "SBS")
      }
    }
    if (mode == "ALL" | mode == "DBS") {
      mat <- tally_list$DBS_78

      if (!is.null(mat)) {
        mat <- t(mat)
        fit_DBS <- sig_fit(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "DBS",
          mode = "DBS",
          return_class = "data.table",
          return_error = TRUE
        )
        output_fit(fit_DBS, result_dir = file.path(result_dir, "results"), mut_type = "DBS")
      }
    }
    if (mode == "ALL" | mode == "ID") {
      mat <- tally_list$ID_83

      if (!is.null(mat)) {
        mat <- t(mat)
        fit_ID <- sig_fit(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "ID",
          mode = "ID",
          return_class = "data.table",
          return_error = TRUE
        )
        output_fit(fit_ID, result_dir = file.path(result_dir, "results"), mut_type = "ID")
      }
    }
  } else {
    ## bootstrap fitting
    if (mode == "ALL" | mode == "SBS") {
      mat <- tally_list$SBS_96

      if (!is.null(mat)) {
        mat <- t(mat)
        ## COSMIC V2 SBS
        bt_legacy <- sig_fit_bootstrap_batch(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "legacy",
          mode = "SBS",
          n = nrun,
          p_val_thresholds = c(0.01, 0.05, 0.10, 0.25),
          use_parallel = TRUE,
          job_id = "legacy",
          result_dir = file.path(result_dir, "bootstrap"),
        )

        output_bootstrap(bt_legacy, result_dir = file.path(result_dir, "results"), mut_type = "legacy")

        bt_SBS <- sig_fit_bootstrap_batch(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "SBS",
          mode = "SBS",
          n = nrun,
          p_val_thresholds = c(0.01, 0.05, 0.10, 0.25),
          use_parallel = TRUE,
          job_id = "SBS",
          result_dir = file.path(result_dir, "bootstrap"),
        )

        output_bootstrap(bt_SBS, result_dir = file.path(result_dir, "results"), mut_type = "SBS")
      }
    }
    if (mode == "ALL" | mode == "DBS") {
      mat <- tally_list$DBS_78

      if (!is.null(mat)) {
        mat <- t(mat)
        bt_DBS <- sig_fit_bootstrap_batch(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "DBS",
          mode = "DBS",
          n = nrun,
          p_val_thresholds = c(0.01, 0.05, 0.10, 0.25),
          use_parallel = TRUE,
          job_id = "DBS",
          result_dir = file.path(result_dir, "bootstrap"),
        )

        output_bootstrap(bt_DBS, result_dir = file.path(result_dir, "results"), mut_type = "DBS")
      }
    }
    if (mode == "ALL" | mode == "ID") {
      mat <- tally_list$ID_83

      if (!is.null(mat)) {
        mat <- t(mat)
        bt_ID <- sig_fit_bootstrap_batch(
          catalogue_matrix = mat,
          sig_index = index,
          sig_db = "ID",
          mode = "ID",
          n = nrun,
          p_val_thresholds = c(0.01, 0.05, 0.10, 0.25),
          use_parallel = TRUE,
          job_id = "ID",
          result_dir = file.path(result_dir, "bootstrap"),
        )

        output_bootstrap(bt_ID, result_dir = file.path(result_dir, "results"), mut_type = "ID")
      }
    }
  }
}

# Action part -------------------------------------------------------------
## Check inputs and call the working functions
suppressMessages(library("sigminer"))
suppressMessages(library("ggplot2"))
suppressMessages(library("pheatmap"))

if (ARGS$verbose) {
  message()
  message("Library:\t", paste(.libPaths(), collapse = "\t"))
  print(utils::sessionInfo())
}

message("Reading file...\n------")

mode <- ARGS$mode
input <- ARGS$input
genome_build <- ARGS$genome
result_dir <- path.expand(ARGS$output)

if (any(endsWith(input, c("xls", "xlsx")))) {
  if (!suppressMessages(require("readxl"))) {
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

if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
  message("Result directory ", result_dir, " created.")
} else {
  message("Result directory ", result_dir, " existed.")
}

if (!isCN) {
  if (!file.exists(file.path(result_dir, "maf_obj.RData"))) {
    if (dir.exists(input)) {
      fs <- list.files(input, pattern = "*.vcf", full.names = TRUE)
      if (length(fs) < 1) {
        message("When input is a directory, it should contain VCF files!")
        quit("no", status = -1)
      }
      message("Try parsing VCF files...")
      obj <- sigminer::read_vcf(fs, genome_build = genome_build, keep_only_pass = FALSE, verbose = TRUE)
    } else {
      message("Try reading as a MAF file.")
      if (!is.data.frame(input)) {
        obj <- suppressMessages(data.table::fread(input, header = TRUE, data.table = FALSE))
      } else {
        obj <- input
      }
      if (!"Tumor_Sample_Barcode" %in% colnames(obj)) {
        message("No 'Tumor_Sample_Barcode' column detected, try parsing as a component-by-sample matrix.")
        rownames(obj) <- obj[[1]]
        obj[[1]] <- NULL
        obj <- as.matrix(obj)
        message("Read as a matrix done.")
      } else {
        obj <- tryCatch(
          sigminer::read_maf(obj, verbose = TRUE),
          error = function(e) {
            message("Read input as a MAF file failed, try parsing as a component-by-sample matrix.")
            rownames(obj) <- obj[[1]]
            obj[[1]] <- NULL
            obj <- as.matrix(obj)
            message("Read as a matrix done.")
            return(obj)
          }
        )
      }
    }
    save(obj, file = file.path(result_dir, "maf_obj.RData"))
  } else {
    load(file = file.path(result_dir, "maf_obj.RData"))
  }
} else {
  if (!file.exists(file.path(result_dir, "cn_obj.RData"))) {
    message("Try reading as a Segment file.")
    if (!is.data.frame(input)) {
      obj <- suppressMessages(data.table::fread(input, header = TRUE, data.table = FALSE))
    } else {
      obj <- input
    }
    if (!"sample" %in% colnames(obj)) {
      message("No 'sample' column detected, try parsing as a component-by-sample matrix.")
      rownames(obj) <- obj[[1]]
      obj[[1]] <- NULL
      obj <- as.matrix(obj)
      message("Read as a matrix done.")
    } else {
      obj <- tryCatch(
        sigminer::read_copynumber(obj, genome_build = genome_build, verbose = TRUE),
        error = function(e) {
          message("Read input as a Segment file failed, try parsing as a component-by-sample matrix.")
          rownames(obj) <- obj[[1]]
          obj[[1]] <- NULL
          obj <- as.matrix(obj)
          message("Read as a matrix done.")
          return(obj)
        }
      )
    }
    save(obj, file = file.path(result_dir, "cn_obj.RData"))
  } else {
    load(file = file.path(result_dir, "cn_obj.RData"))
  }
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
  if (as.integer(ARGS$max) < 2) {
    max_number <- 100
  } else {
    max_number <- as.integer(ARGS$max)
  }
  tryCatch(
    {
      if (ARGS$verbose) {
        flow_extraction(
          obj = obj, genome_build = genome_build, mode = ARGS$mode,
          manual_step = manual_step, nrun = nrun, cores = cores,
          result_dir = result_dir,
          max_number = max_number,
          refit = ARGS$refit,
          rm_hyper = ARGS$hyper,
          sigprofiler = ARGS$sigprofiler
        )
      } else {
        suppressMessages(
          flow_extraction(
            obj = obj, genome_build = genome_build, mode = ARGS$mode,
            manual_step = manual_step, nrun = nrun, cores = cores,
            result_dir = result_dir,
            max_number = max_number,
            refit = ARGS$refit,
            rm_hyper = ARGS$hyper,
            sigprofiler = ARGS$sigprofiler
          )
        )
      }
    },
    error = function(e) {
      message("An error is detected in extract process. Quit the program!")
      message(e$message)
      quit("no", -1)
    }
  )
} else if (ARGS$fit) {
  message("Running signature fitting pipeline...\n------")
  tryCatch(
    {
      if (ARGS$verbose) {
        flow_fitting(
          obj = obj, genome_build = genome_build, mode = ARGS$mode,
          result_dir = result_dir, index = ARGS$index, prog = "fit"
        )
      } else {
        suppressMessages(
          flow_fitting(
            obj = obj, genome_build = genome_build, mode = ARGS$mode,
            result_dir = result_dir, index = ARGS$index, prog = "fit"
          )
        )
      }
    },
    error = function(e) {
      message("An error is detected in fitting process. Quit the program!")
      message(e$message)
      quit("no", -1)
    }
  )
} else if (ARGS$bt) {
  message("Running signature bootstrap fitting pipeline...\n------")
  nrun <- as.integer(ARGS$nrun)
  tryCatch(
    {
      if (ARGS$verbose) {
        flow_fitting(
          obj = obj, genome_build = genome_build, mode = ARGS$mode,
          result_dir = result_dir, nrun = nrun, index = ARGS$index, prog = "bootstrap"
        )
      } else {
        suppressMessages(
          flow_fitting(
            obj = obj, genome_build = genome_build, mode = ARGS$mode,
            result_dir = result_dir, nrun = nrun, index = ARGS$index, prog = "bootstrap"
          )
        )
      }
    },
    error = function(e) {
      message("An error is detected in bootstrap fitting process. Quit the program!")
      message(e$message)
      quit("no", -1)
    }
  )
}

# End part ----------------------------------------------------------------

message("
Please check the output directory.
Thanks for using Sigflow!
============================== END
")

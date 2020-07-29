#!/usr/bin/env Rscript

"=================================================================
sigflow: Streamline Analysis Workflow for Mutational Signatures.

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
  fit     - fit signatures in >=1 samples.
  ==
  bt      - run bootstrap signature fitting analysis in >=1 samples.

Usage:
  sigflow extract --input=<file> [--output=<outdir>] [--mode=<class>] [--manual --number <sigs>] [--genome=<genome>] [--nrun=<runs>] [--cores=<cores>]
  sigflow fit --input=<file> [--output=<outdir>] [--mode=<class>] [--genome=<genome>]
  sigflow bt --input=<file> [--output=<outdir>] [--mode=<class>] [--genome=<genome>] [--nrun=<runs>]
  sigflow (-h | --help)
  sigflow --version

Options:
  -h --help     Show help message.
  --version     Show version.
  -i <file>, --input <file>       input file/directory path.
  -o <outdir>, --output <outdir>  output directory path [default: ./sigflow_result/].
  -m <class>, --mode <class>      extract/fit mode, can be one of SBS, DBS, ID, MAF (for three types), CN (not supported in fit subcommand) [default: SBS].
  --manual                        enable manual extraction, set -N=0 for outputing signature estimation firstly.
  -N <sigs>, --number <sigs>      extract specified number of signatures [default: 0].
  -g <genome>, --genome <genome>  genome build, can be hg19, hg38 or mm10, [default: hg19].
  -r <runs>, --nrun <runs>        run times of NMF (extract) or bootstrapping (bt) to get results [default: 30].
  -T <cores>, --cores <cores>     cores to run the program, large dataset will benefit from it [default: 1].

=================================================================
" -> doc

if (!require("docopt")) {
  install.packages('docopt', repos = 'https://cloud.r-project.org')
}

library("docopt")
arguments <- docopt(doc, version = "sigflow v0.1\n")

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
Version:       0.1
License:       AFL 3.0
Link   :       https://github.com/ShixiangWang/sigminer.workflow
Doc    :       https://shixiangwang.github.io/sigminer-doc/
============================== START
")


# Parsing input -----------------------------------------------------------

message("Parsing parameters...\n------")
ARGS <- arguments[!startsWith(names(arguments), "--")]
for (i in seq_along(ARGS)) {
  message(names(ARGS)[i], "\t\t: ", ARGS[i])
}
message("------\n")

# Function part -----------------------------------------------------------
## Program to go

output_tally <- function(x, result_dir, mut_type = "SBS") {
  ## x is a matrix with row representing components (motifs) and column representing samples
  ## output all samples in bar plots
  message("Outputing tally results for ", mut_type)

  if (is.list(x)) {
    x <- x$nmf_matrix %>% t()
  }

  samps <- colnames(x)

  if (mut_type != "CN") {
    p_total <- show_catalogue(x, style = "cosmic", mode = mut_type, x_label_angle = 90, x_label_vjust = 0.5)
  } else {
    p_total <- show_catalogue(x, style = "cosmic", mode = "copynumber")
  }

  ggsave(file.path(result_dir, paste0(mut_type, "_tally_total.pdf")),
    plot = p_total, width = 12, height = 3
  )

  if (mut_type != "CN") {
    p_samps <- show_catalogue(x,
      style = "cosmic", mode = mut_type,
      samples = samps,
      x_label_angle = 90, x_label_vjust = 0.5
    )
  } else {
    p_samps <- show_catalogue(x, 
                              style = "cosmic", mode = "copynumber", 
                              normalize = "feature",
                              samples = samps)
  }

  ggsave(file.path(result_dir, paste0(mut_type, "_tally_samps.pdf")),
    plot = p_samps, width = 12, height = 2 * length(samps), limitsize = FALSE
  )
  
  return(invisible(NULL))
}

output_sig <- function(sig, result_dir, mut_type = "SBS") {
  message("Outputing signature results for ", mut_type, "\n==")
  
  stopifnot(inherits(sig, "Signature"))
  
  message("Outputing signature object, signature and exposure matrix.")
  ## Output Signature object
  saveRDS(sig, file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signatureObj.rds")))

  ## Output data
  if (mut_type != "CN") {
    data.table::fwrite(sig_signature(sig) %>% data.table::as.data.table(keep.rownames = "component"),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature.csv"))
    )
  } else {
    data.table::fwrite(sig_signature(sig, normalize = "feature") %>% data.table::as.data.table(keep.rownames = "component"),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature.csv"))
    )
  }
  data.table::fwrite(get_sig_exposure(sig),
    file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_exposure.csv"))
  )
  
  message("Outputing sample clusters based on signature contribution.")
  if (sig$K > 1) {
    message("=> Running k-means clustering.")
    grp <- tryCatch(
      get_groups(sig, method = "k-means"),
      error = function(e) {
        message("==> Error detected when running k-means, switch to directly assign samples based on exposure")
        get_groups(sig, method = "exposure")
      }
    )
    data.table::fwrite(grp,
                       file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_kmeans_cluster.csv"))
    )
    pdf(file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_kmeans_cluster.pdf")), 
        width = 1.5 * length(unique(grp$group)), height = 4)
    show_groups(grp)
    dev.off()
  } else {
    message("Skip clustering only one signatures.")
  }

  ## Output plots
  message("Outputing signature profile plot.")
  if (mut_type != "CN") {
    p <- show_sig_profile(sig, mode = mut_type, style = "cosmic", x_label_angle = 90, x_label_vjust = 0.5)
  } else {
    p <- show_sig_profile(sig, mode = "copynumber", normalize = "feature", style = "cosmic")
  }
  ggsave(file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature_profile.pdf")),
    plot = p, width = 12, height = 2 * sig$K
  )

  message("Outputing signature exposure plot.")
  if (ncol(sig$Exposure) < 50) {
    p <- show_sig_exposure(sig, style = "cosmic", hide_samps = FALSE)
  } else {
    p <- show_sig_exposure(sig, style = "cosmic", hide_samps = TRUE, rm_space = TRUE)
  }
  cowplot::save_plot(file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_exposure_profile.pdf")),
    plot = p
  )

  ## Similar analysis and output
  if (mut_type != "CN") {
    message("Outputing signature similarity analysis results.")
    sim <- get_sig_similarity(sig, sig_db = mut_type)
    data.table::fwrite(sim$similarity %>% data.table::as.data.table(keep.rownames = "sig"),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_similarity.csv"))
    )
    pheatmap::pheatmap(sim$similarity, cluster_cols = TRUE, cluster_rows = FALSE,
                       filename = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_similarity.pdf")),
                       cellheight = 15, fontsize = 7)
    data.table::fwrite(sim$best_match %>% data.table::as.data.table(),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_COSMIC_best_match.csv"))
    )
  }
  
  return(invisible(NULL))
}

output_fit <- function(x, result_dir, mut_type = "SBS") {
  message("Outputing signature fitting results for ", mut_type)
  
  expo <- x$expo
  error <- x$errors %>% data.table::as.data.table(keep.rownames = TRUE)
  colnames(error) <- c("sample", "error")
  
  rel_expo <- expo %>% 
    dplyr::mutate(s = rowSums(.[, -1])) %>%
    dplyr::mutate_at(dplyr::vars(-c("sample", "s")), ~./.data$s) %>%
    dplyr::select(-"s") %>%
    data.table::as.data.table()
  
  data.table::fwrite(expo, file = file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure.csv")))
  data.table::fwrite(rel_expo, file = file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure.csv")))
  data.table::fwrite(error, file = file.path(result_dir, paste0(mut_type, "_fitting_reconstruction_errors.csv")))
  
  if (mut_type != "SBS") {
    p1 <- show_sig_fit(expo, palette = NULL, plot_fun = "boxplot") + ggpubr::rotate_x_text()
    p2 <- show_sig_fit(expo, palette = NULL, plot_fun = "violin") + ggpubr::rotate_x_text()
    
    p3 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "boxplot") + ggpubr::rotate_x_text()
    p4 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "violin") + ggpubr::rotate_x_text()
  } else {
    z <- get_sig_db("SBS")
    sigs <- rownames(z$aetiology)[!grepl("artefact", z$aetiology$aetiology)]
    message("Removed 'Possible sequencing artefact' signatures in plots.")
    
    p1 <- show_sig_fit(expo, palette = NULL, plot_fun = "boxplot", signatures = sigs) + ggpubr::rotate_x_text()
    p2 <- show_sig_fit(expo, palette = NULL, plot_fun = "violin", signatures = sigs) + ggpubr::rotate_x_text()
    
    p3 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "boxplot", signatures = sigs) + ggpubr::rotate_x_text()
    p4 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "violin", signatures = sigs) + ggpubr::rotate_x_text()
  }
  
  if (mut_type == "legacy") {
    ## 30 signatures
    width <- 10
    height <- 4
  } else if (mut_type == "SBS") {
    ## 49 signatures
    width <- 15
    height <- 4
  } else if (mut_type == "DBS") {
    ## 11 signatures
    width <- 5
    height <- 4
  } else if (mut_type == "ID") {
    ## 17 signatures
    width <- 7
    height <- 4
  }
  
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure_boxplot.pdf")),
         plot = p1, width = width, height = height)
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure_violin.pdf")),
         plot = p2, width = width, height = height)
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure_boxplot.pdf")),
         plot = p3, width = width, height = height)
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure_violin.pdf")),
         plot = p4, width = width, height = height)
}

output_bootstrap <- function(x, result_dir, mut_type = "SBS") {
  message("Outputing signature bootstrap fitting results for ", mut_type)
  
  data.table::fwrite(x$expo, file = file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure.csv")))
  data.table::fwrite(x$error, file = file.path(result_dir, paste0(mut_type, "_bootstrap_reconstruction_errors.csv")))
  data.table::fwrite(x$p_val, file = file.path(result_dir, paste0(mut_type, "_bootstrap_p_values_under_different_exposure_cutoffs.csv")))
  
  if (mut_type != "SBS") {
    p1 <- show_sig_bootstrap_stability(x) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    p2 <- show_sig_bootstrap_exposure(x) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    
  } else {
    z <- get_sig_db("SBS")
    sigs <- rownames(z$aetiology)[!grepl("artefact", z$aetiology$aetiology)]
    message("Removed 'Possible sequencing artefact' signatures in plots.")
    
    p1 <- show_sig_bootstrap_stability(x, signatures = sigs) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    p2 <- show_sig_bootstrap_exposure(x, signatures = sigs) + theme(legend.position = "none") + ggpubr::rotate_x_text()
  }
  
  if (mut_type == "legacy") {
    ## 30 signatures
    width <- 10
    height <- 4
  } else if (mut_type == "SBS") {
    ## 49 signatures
    width <- 15
    height <- 4
  } else if (mut_type == "DBS") {
    ## 11 signatures
    width <- 5
    height <- 4
  } else if (mut_type == "ID") {
    ## 17 signatures
    width <- 7
    height <- 4
  }
  
  ggsave(file.path(result_dir, paste0(mut_type, "_bootstrap_signature_instability_boxplot.pdf")),
         plot = p1, width = width, height = height)
  ggsave(file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure_boxplot.pdf")),
         plot = p2, width = width, height = height)
  
  samps <- unique(x$expo$sample)
  samp_dir <- file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure_per_sample_boxplot"))
  if (!dir.exists(samp_dir)) {
    dir.create(samp_dir, recursive = TRUE)
  }
  for (i in samps) {
    message("Plotting sample: ", i)
    if (mut_type == "SBS") {
      p <- show_sig_bootstrap_exposure(x, signatures = sigs, sample = i) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    } else {
      p <- show_sig_bootstrap_exposure(x, sample = i) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    }
    ggsave(file.path(samp_dir, paste0(i, ".pdf")),
           plot = p, width = width, height = height)
  }
}

flow_extraction <- function(obj, genome_build, mode, manual_step, nrun, cores, result_dir) {
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

  if (!require(ref_genome, character.only = TRUE)) {
    sigminer:::send_info("Package ", ref_genome, " not found, try installing.")
    BiocManager::install(ref_genome)
  }

  if (mode != "CN") {
    if (mode == "MAF") {
      mode <- "ALL"
    }
    if (manual_step <= 0) {
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
    }
  } else {
    if (manual_step <= 0) {
      tally_list <- sig_tally(obj, ignore_chrs = c("chrX", "chrY"), cores = cores)
      save(tally_list, file = file.path(result_dir, "cn_tally.RData"))
      output_tally(tally_list, result_dir = file.path(result_dir, "results"), mut_type = "CN")
    }
  }

  if (manual_step < 0) {
    ## auto-extract
    if (mode == "CN") {
      sigs_CN <- sig_auto_extract(
        nmf_matrix = tally_list$nmf_matrix,
        result_prefix = "BayesianNMF_CN",
        destdir = file.path(result_dir, "BayesianNMF"),
        strategy = "stable",
        K0 = 30L,
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
            strategy = "stable",
            K0 = 30L,
            nrun = nrun,
            cores = cores,
            skip = TRUE
          )
        } else {
          sigminer:::send_info("SBS_96 matrix is NULL, skip extracting.")
        }
      }
      if (mode == "ALL" | mode == "DBS") {
        mat <- tally_list$DBS_78
        if (!is.null(mat)) {
          sigs_DBS <- sig_auto_extract(
            nmf_matrix = mat,
            result_prefix = "BayesianNMF_DBS",
            destdir = file.path(result_dir, "BayesianNMF"),
            strategy = "stable",
            K0 = 30L,
            nrun = nrun,
            cores = cores,
            skip = TRUE
          )
        } else {
          sigminer:::send_info("DBS_78 matrix is NULL, skip extracting.")
        }
      }
      if (mode == "ALL" | mode == "ID") {
        mat <- tally_list$ID_83
        if (!is.null(mat)) {
          sigs_ID <- sig_auto_extract(
            nmf_matrix = mat,
            result_prefix = "BayesianNMF_ID",
            destdir = file.path(result_dir, "BayesianNMF"),
            strategy = "stable",
            K0 = 30L,
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
    if (manual_step == 0) {
      suppressPackageStartupMessages(library("NMF"))
      ## do signature estimation
      if (mode == "CN") {
        est_CN <- sig_estimate(
          nmf_matrix = tally_list$nmf_matrix,
          nrun = nrun,
          range = 2:30,
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
              range = 2:30,
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
              range = 2:30,
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
              range = 2:30,
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
          optimize = TRUE,
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
            sigs_SBS <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = TRUE,
              pConstant = 1e-9
            )
          }
        } 
        if (mode == "DBS") {
          mat <- tally_list$DBS_78

          if (!is.null(mat)) {
            sigs_DBS <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = TRUE,
              pConstant = 1e-9
            )
          }
        }
        if (mode == "ID") {
          mat <- tally_list$ID_83

          if (!is.null(mat)) {
            sigs_ID <- sig_extract(
              nmf_matrix = mat,
              n_sig = manual_step,
              nrun = nrun,
              cores = cores,
              optimize = TRUE,
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



flow_fitting <- function(obj, genome_build, mode, result_dir, nrun = NULL, prog = c("fit", "bootstrap")) {
  
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
  
  if (!require(ref_genome, character.only = TRUE)) {
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
          sig_index = "ALL",
          sig_db = "legacy",
          mode = "SBS",
          return_class = "data.table",
          return_error = TRUE
        )
        output_fit(fit_legacy, result_dir = file.path(result_dir, "results"), mut_type = "legacy")
        
        fit_SBS <- sig_fit(
          catalogue_matrix = mat,
          sig_index = "ALL",
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
          sig_index = "ALL",
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
          sig_index = "ALL",
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
          sig_index = "ALL",
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
          sig_index = "ALL",
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
          sig_index = "ALL",
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
          sig_index = "ALL",
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
suppressPackageStartupMessages(library("sigminer"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("pheatmap"))

message("Reading file...\n------")

mode <- ARGS$mode
input <- ARGS$input
genome_build <- ARGS$genome
result_dir <- path.expand(ARGS$output)

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
      obj <- sigminer::read_maf(input, verbose = TRUE)
    }
    save(obj, file = file.path(result_dir, "maf_obj.RData"))
  } else {
    load(file = file.path(result_dir, "maf_obj.RData"))
  }
} else {
  if (!file.exists(file.path(result_dir, "cn_obj.RData"))) {
    obj <- sigminer::read_copynumber(input, genome_build = genome_build, verbose = TRUE)
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
  tryCatch(
    flow_extraction(
      obj = obj, genome_build = genome_build, mode = ARGS$mode,
      manual_step = manual_step, nrun = nrun, cores = cores,
      result_dir = result_dir
    ),
    error = function(e) {
      message("An error is detected in extract process. Quit the program!")
      message(e$message)
      quit("no", -1)
    }
  )
} else if (ARGS$fit) {
  message("Running signature fitting pipeline...\n------")
  tryCatch(
    flow_fitting(
      obj = obj, genome_build = genome_build, mode = ARGS$mode,
      result_dir = result_dir, prog = "fit"),
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
    flow_fitting(
      obj = obj, genome_build = genome_build, mode = ARGS$mode,
      result_dir = result_dir, nrun = nrun, prog = "bootstrap"),
    error = function(e) {
      message("An error is detected in bootstrap fitting process. Quit the program!")
      message(e$message)
      quit("no", -1)
    }
  )
}

# End part ----------------------------------------------------------------

message("
============================== END
")

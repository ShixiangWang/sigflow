## Modify based on code from https://stackoverflow.com/questions/46902203/verify-r-packages-installed-into-docker-container
c("remotes", "sigminer", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10") -> chk_pkgs

ret <- suppressPackageStartupMessages(
  sapply(chk_pkgs, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
)

missing_pkgs <- names(ret[!ret])

if (length(missing_pkgs) > 0) {
  warning("The following packages are not installed: ",
    paste0(sprintf("  - %s", missing_pkgs), collapse = "\n"),
    immediate. = TRUE
  )
  message("Try installing them...")

  if ("sigminer" %in% missing_pkgs) {
    BiocManager::install("ShixiangWang/sigminer@v2.3.0", dependencies = TRUE)
  }

  rem_pkgs <- setdiff(missing_pkgs, "sigminer")

  if (length(rem_pkgs) > 0) {
    BiocManager::install(rem_pkgs)
  }

  ret <- suppressPackageStartupMessages(
    sapply(chk_pkgs, require, character.only = TRUE, quietly = FALSE, warn.conflicts = FALSE)
  )

  missing_pkgs <- names(ret[!ret])

  if (length(missing_pkgs) > 0) {
    warning("Still cannot fix installing problem, quit with -1.", immediate. = TRUE)
  }
}

data_type = c("chr_size", "centro_loc", "cytobands", "transcript", "gene")

for (i in data_type) {
  sigminer::get_genome_annotation(data_type = i, genome_build = "hg19")
  sigminer::get_genome_annotation(data_type = i, genome_build = "hg38")
  sigminer::get_genome_annotation(data_type = i, genome_build = "mm10")
}

quit(save = "no", status = ifelse(length(missing_pkgs) == 0, 0, -1), runLast = FALSE)

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
    BiocManager::install("ShixiangWang/sigminer@v1.9", dependencies = TRUE)
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

quit(save = "no", status = ifelse(length(missing_pkgs) == 0, 0, -1), runLast = FALSE)

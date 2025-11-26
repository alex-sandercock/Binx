#!/usr/bin/env Rscript

# Compare binx GWAS output to GWASpoly results.
# Default column names follow binx TSV (marker_id, beta, p_value) and GWASpoly tutorial outputs.
# Usage example:
# Rscript scripts/compare_binx_gwaspoly.R \
#   --binx=binx_out.tsv \
#   --gwaspoly=gwaspoly_out.tsv \
#   --trait=TraitName \
#   --out=benchmark/merged.tsv

options(stringsAsFactors = FALSE)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list(
    binx = NULL,
    gwaspoly = NULL,
    trait = NULL,
    out = "merged.tsv",
    binx_marker = "marker_id",
    binx_beta = "beta",
    binx_p = "p_value",
    gwaspoly_marker = "Marker",
    gwaspoly_beta = "EffB",
    gwaspoly_p = "pValue"
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else "TRUE"
    parsed[[key]] <- val
  }
  required <- c("binx", "gwaspoly", "trait")
  missing <- required[sapply(parsed[required], is.null)]
  if (length(missing) > 0) {
    stop(paste("Missing required args:", paste(missing, collapse = ", ")))
  }
  parsed
}

read_table <- function(path) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "", quote = "\"")
}

main <- function() {
  args <- parse_args()
  cat("Reading binx results from", args$binx, "\n")
  b <- read_table(args$binx)
  cat("Reading GWASpoly results from", args$gwaspoly, "\n")
  g <- read_table(args$gwaspoly)

  rename_safe <- function(df, old, new) {
    if (!(old %in% names(df))) {
      stop(paste("Column", old, "not found"))
    }
    names(df)[names(df) == old] <- new
    df
  }

  b <- rename_safe(b, args$binx_marker, "marker")
  b <- rename_safe(b, args$binx_beta, "beta_binx")
  b <- rename_safe(b, args$binx_p, "p_binx")
  g <- rename_safe(g, args$gwaspoly_marker, "marker")
  g <- rename_safe(g, args$gwaspoly_beta, "beta_gwaspoly")
  g <- rename_safe(g, args$gwaspoly_p, "p_gwaspoly")

  b$mlogp_binx <- -log10(b$p_binx)
  g$mlogp_gwaspoly <- -log10(g$p_gwaspoly)

  merged <- merge(b[, c("marker", "beta_binx", "p_binx", "mlogp_binx")],
                  g[, c("marker", "beta_gwaspoly", "p_gwaspoly", "mlogp_gwaspoly")],
                  by = "marker")
  cat("Merged markers:", nrow(merged), "\n")

  corr <- function(x, y) {
    if (length(x) < 2) return(NA_real_)
    stats::cor(x, y, use = "complete.obs")
  }

  summary <- list(
    beta_corr = corr(merged$beta_binx, merged$beta_gwaspoly),
    mlogp_corr = corr(merged$mlogp_binx, merged$mlogp_gwaspoly),
    n = nrow(merged)
  )

  cat("Beta correlation:", summary$beta_corr, "\n")
  cat("-log10(p) correlation:", summary$mlogp_corr, "\n")

  dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)
  write.table(merged, file = args$out, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Wrote merged table to", args$out, "\n")
}

tryCatch(main(), error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

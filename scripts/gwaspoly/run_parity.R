#!/usr/bin/env Rscript

# Generate GWASpoly reference outputs for binx parity tests.
# Inputs are passed as --key=value arguments.
# Example:
# Rscript scripts/gwaspoly/run_parity.R \
#   --geno=tests/parity/data/toy/toy.geno.tsv \
#   --pheno=tests/parity/data/toy/toy.pheno.tsv \
#   --trait=Trait1 \
#   --ploidy=4 \
#   --models=additive \
#   --out=tests/parity/fixtures/toy_additive.tsv
#
# With env fixed effect (repeated IDs):
# Rscript scripts/gwaspoly/run_parity.R \
#   --geno=tests/parity/data/potato/new_potato_geno.csv \
#   --pheno=tests/parity/data/potato/new_potato_pheno.csv \
#   --trait=vine.maturity \
#   --ploidy=4 \
#   --fixed=env \
#   --fixed_type=factor \
#   --models=additive \
#   --out=tests/parity/fixtures/potato_additive_env.tsv

options(stringsAsFactors = FALSE)

parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  args <- list(
    geno = NULL,
    pheno = NULL,
    trait = NULL,
    ploidy = NULL,
    out = NULL,
    models = "additive",
    kinship = NULL,
    kinship_out = NULL,
    fixed = NULL,
    fixed_type = NULL,
    format = "numeric",
    n_traits = 1,
    delim = "",
    min_maf = NA,
    max_missing = NA,
    kinship_out = NULL
  )
  for (a in raw) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else "TRUE"
    if (key %in% names(args)) {
      args[[key]] <- val
    } else {
      stop(paste("Unknown argument:", key))
    }
  }
  required <- c("geno", "pheno", "trait", "ploidy", "out")
  missing <- required[sapply(args[required], function(x) is.null(x) || x == "")]
  if (length(missing) > 0) {
    stop(paste("Missing required args:", paste(missing, collapse = ", ")))
  }
  args$ploidy <- as.integer(args$ploidy)
  args$n_traits <- as.integer(args$n_traits)
  if (!is.na(args$min_maf)) args$min_maf <- as.numeric(args$min_maf)
  if (!is.na(args$max_missing)) args$max_missing <- as.numeric(args$max_missing)
  args
}

detect_delim <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
}

main <- function() {
  args <- parse_args()
  message("Loading GWASpoly...")
  suppressPackageStartupMessages(library(GWASpoly))

  models <- unlist(strsplit(args$models, ","))
  fixed <- if (!is.null(args$fixed) && args$fixed != "") unlist(strsplit(args$fixed, ",")) else NULL
  fixed_type <- if (!is.null(args$fixed_type) && args$fixed_type != "") unlist(strsplit(args$fixed_type, ",")) else NULL

  delim <- if (!is.null(args$delim) && args$delim != "") args$delim else detect_delim(args$geno)

  message("Reading genotype + phenotype (format=", args$format, ", delim='", delim, "')")
  gwas <- read.GWASpoly(
    args$ploidy,
    args$pheno,
    args$geno,
    format = args$format,
    n.traits = args$n_traits,
    delim = delim
  )

  # Compute or load kinship
  if (!is.null(args$kinship) && args$kinship != "") {
    message("Using provided kinship: ", args$kinship)
    gwas <- set.K.from.file(gwas, args$kinship, LOCO = FALSE)
  } else {
    message("Computing kinship via set.K")
    gwas <- set.K(gwas, LOCO = FALSE)
  }

  # Optionally write kinship matrix to file for reuse in parity tests
  if (!is.null(args$kinship_out) && args$kinship_out != "") {
    message("Writing kinship to ", args$kinship_out)
    dir.create(dirname(args$kinship_out), showWarnings = FALSE, recursive = TRUE)
    # gwas@K is a list (length 1 when LOCO=FALSE)
    k_mat <- gwas@K[[1]]
    sid <- rownames(k_mat)
    if (is.null(sid)) {
      sid <- colnames(k_mat)
    }
    if (is.null(sid)) {
      stop("Kinship matrix missing row/column names; cannot write kinship_out")
    }
    kin_df <- data.frame(sample_id = sid, k_mat, check.names = FALSE)
    write.table(kin_df, file = args$kinship_out, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Parameters (mirrors GWASpoly defaults; min_maf and max_missing are hooks if we need to tweak)
  geno_freq <- 1 - 5 / nrow(gwas@geno)
  geno_freq <- max(0.01, min(0.99, geno_freq)) # ensure within (0,1)
  params <- set.params(geno.freq = geno_freq)
  if (!is.null(fixed)) {
    params$fixed <- fixed
    if (!is.null(fixed_type)) {
      params$fixed.type <- fixed_type
    }
  }
  if (!is.na(args$min_maf)) {
    params$min.MAF <- args$min_maf
  }
  if (!is.na(args$max_missing)) {
    params$max.missing <- args$max_missing
  }

  message("Running GWASpoly for models: ", paste(models, collapse = ", "))
  gwas <- GWASpoly(gwas, models = models, traits = args$trait, params = params)

  dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)

  for (model in models) {
    out_path <- args$out
    if (length(models) > 1) {
      pieces <- strsplit(out_path, "\\.", fixed = FALSE)[[1]]
      if (length(pieces) >= 2) {
        out_path <- paste0(paste(head(pieces, -1), collapse = "."), "_", model, ".", tail(pieces, 1))
      } else {
        out_path <- paste0(out_path, "_", model, ".tsv")
      }
    }
    message("Writing results for model '", model, "' to ", out_path)

    scores_list <- gwas@scores[[args$trait]]
    effects_list <- gwas@effects[[args$trait]]
    if (!(model %in% colnames(scores_list))) {
      stop(paste("Model", model, "not found in scores list"))
    }
    score_vec <- scores_list[[model]]
    effect_vec <- if (model %in% colnames(effects_list)) effects_list[[model]] else rep(NA_real_, length(score_vec))

    res <- data.frame(
      marker_id = rownames(scores_list),
      model = model,
      score = score_vec,
      p_value = 10^(-score_vec),
      effect = effect_vec,
      stringsAsFactors = FALSE
    )

    # attach map info if available
    if (!is.null(gwas@map)) {
      map <- gwas@map
      res <- merge(
        map[, c("Marker", "Chrom", "Position")],
        res,
        by.x = "Marker",
        by.y = "marker_id",
        all.y = TRUE
      )
      names(res)[names(res) == "Marker"] <- "marker_id"
      names(res)[names(res) == "Chrom"] <- "chrom"
      names(res)[names(res) == "Position"] <- "pos"
    }

    # n_obs for this trait
    if (args$trait %in% colnames(gwas@pheno)) {
      res$n_obs <- sum(!is.na(gwas@pheno[[args$trait]]))
    }

    write.table(res, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

tryCatch(main(), error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

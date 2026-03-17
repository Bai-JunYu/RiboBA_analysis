# Part 1: similar to the RiboCode paper setup.
# Select 1000 protein-coding genes (RPF count >= 10) as positives.
# Create negatives by randomly shifting RPF positions by +/- 3 nt on selected genes.
# Enzymes/samples: rnasei SRR23242345, mnase SRR7073124, p1 SRR23242346.
## ---- Basics ----
options(max.print = 100)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  if (length(idx) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", args[idx[1]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

SCRIPT_DIR <- get_script_dir()
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
PUBLIC_INPUTS_DIR <- Sys.getenv("PUBLIC_INPUTS_DIR", unset = file.path(PROJECT_DIR, "public_inputs"))
RIBOBC_DIR <- Sys.getenv("RIBOBC_DIR", unset = file.path(PROJECT_DIR, "RiboBC"))
SIM_DIR <- file.path(PUBLIC_INPUTS_DIR, "data", "simulate1", "simu_fastq")

if (!dir.exists(PUBLIC_INPUTS_DIR)) stop("PUBLIC_INPUTS_DIR not found: ", PUBLIC_INPUTS_DIR)
if (!dir.exists(RIBOBC_DIR)) stop("RIBOBC_DIR not found: ", RIBOBC_DIR)

## Load helper utilities
source(file.path(RIBOBC_DIR, "R", "helper.R"))

## Transcript annotation
transcript.file <- file.path(PUBLIC_INPUTS_DIR, "data", "try3", "transcript", "tx_information.RData")
.env <- new.env()
load(transcript.file, envir = .env)
tx_info <- if (exists("tx_info_lst", envir=.env)) .env$tx_info_lst$tx_info else {
  if (exists("tx_info", envir=.env)) .env$tx_info else stop("tx_info not found in RData")
}

## ---- Helpers: find BAM+BAI (supports .bam.bai and .bai) ----
get_bam_bai_by_srr <- function(dir, srr){
  id <- sub("^SRR", "", srr)
  bam <- list.files(dir, pattern = paste0(id, ".*\\.bam$"), full.names = TRUE)
  if (length(bam) < 1) stop("No BAM found for ", srr, " in ", dir)
  bam <- bam[1]  # Use the first match; enforce uniqueness if needed.
  candidates <- c(paste0(bam, ".bai"),
                  sub("\\.bam$", ".bam.bai", bam),
                  sub("\\.bam$", ".bai", bam))
  bai <- candidates[file.exists(candidates)][1]
  if (is.na(bai)) stop("BAI not found for ", bam)
  list(bam = bam, bai = bai)
}

## ---- These 3 samples only ----
samples <- list(
  list(label = "rnasei", srr = "SRR23242345",
       dir = file.path(PUBLIC_INPUTS_DIR, "data", "p1", "rnase", "rnasei", "map_cds")),
  list(label = "mnase",  srr = "SRR7073124",
       dir = file.path(PUBLIC_INPUTS_DIR, "data", "p1", "rnase", "mnase", "map_cds")),
  list(label = "p1",     srr = "SRR23242346",
       dir = file.path(PUBLIC_INPUTS_DIR, "data", "p1", "rnase", "p1", "map_cds"))
)

## ---- Run input_bam() ----
reads_bam3 <- list()
for (s in samples){
  files <- get_bam_bai_by_srr(s$dir, s$srr)
  message("Loading ", s$label, " (", s$srr, ")")
  reads_bam3[[paste0(s$label, "_", s$srr)]] <- input_bam(
    bam_path     = files$bam,
    bai_path     = files$bai,
    tx_info      = tx_info,
    only_uni_map = TRUE,
    add_ratio    = 0.001
  )
}


suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(ShortRead)
  library(dplyr)
  library(tidyr)
  library(readr)
})

## ===================== INPUTS =====================
## Existing objects:
##   - reads_bam3: list( rnasei_* = list(SRR= list(no_add df, add5 df)), mnase_*, p1_* )
##   - tx_info   : list(flt_tx_seqs, tx_lens, mix_tx, mix_tx_pos)

## ---- Build transcript coordinates and attributes on the concatenated reference ----

build_tx_table <- function(tx_info) {
  stopifnot(all(c("flt_tx_seqs","tx_lens","mix_tx","mix_tx_pos") %in% names(tx_info)))
  tl <- tx_info$tx_lens
  mp <- tx_info$mix_tx_pos
  stopifnot(nrow(tl) == nrow(mp))

  tx_tbl <- tl %>%
    mutate(
      tx_start        = mp$utr5_p5,
      tx_end          = mp$utr3_p3,
      cds_gpos_start  = mp$utr5_p3 + 1L,
      cds_gpos_end    = mp$utr3_p5 - 1L
    ) %>%
    select(tx_id, tx_name, gene_id, utr5_len, cds_len,
           tx_start, tx_end, cds_gpos_start, cds_gpos_end)

  gr <- GRanges(seqnames = "mix_tx",
                ranges   = IRanges(start = tx_tbl$tx_start, end = tx_tbl$tx_end),
                tx_id    = tx_tbl$tx_id,
                tx_name  = tx_tbl$tx_name,
                gene_id  = tx_tbl$gene_id,
                tx_start = tx_tbl$tx_start,
                tx_end   = tx_tbl$tx_end)
  list(tx_tbl = tx_tbl, tx_gr = gr)
}

## ---- Flatten reads for one sample ----
flatten_sample_reads <- function(sample_list) {
  stopifnot(is.list(sample_list) && length(sample_list) == 1L)
  leaf <- sample_list[[1]]

  no_add <- leaf$no_add %>%
    mutate(class = "no_add", seq5 = NA_character_)

  add5 <- leaf$add5 %>%
    mutate(class = "add5", seq5 = seq) %>%
    select(-seq)

  bind_rows(no_add, add5)
}

## ---- Map reads to transcript intervals on the concatenated reference (fully contained) ----
map_reads_to_tx <- function(df, tx_gr) {
  read_gr <- GRanges(seqnames = "mix_tx",
                     ranges   = IRanges(start = df$pos, width = df$qwidth))
  hits <- findOverlaps(read_gr, tx_gr, type = "within", select = "first")
  keep <- !is.na(hits)
  df   <- df[keep, , drop = FALSE]
  tx_m <- as.data.frame(mcols(tx_gr))[hits[keep], ]

  df$tx_name  <- tx_m$tx_name
  df$tx_id    <- tx_m$tx_id
  df$gene_id  <- tx_m$gene_id
  df$tx_start <- tx_m$tx_start
  df$tx_end   <- tx_m$tx_end
  df
}

## ---- Pick representative transcripts from CDS-length %% 3 == 0 and sample 1000 genes ----
# Strictly use CDS-overlap read counts to filter genes and choose representative transcripts.
# df_map must include: gene_id, tx_name, pos, qwidth (pos is 5' end on mix_tx coordinates).
# tx_tbl must include: tx_name, gene_id, cds_len.
# mix_tx_pos row order must match tx_tbl row order.

pick_genes_and_tx_cds <- function(df_map, tx_tbl, mix_tx_pos,
                                  n_genes = 1000, min_reads = 10L) {
  stopifnot(all(c("gene_id","tx_name","pos","qwidth") %in% names(df_map)))
  stopifnot(all(c("tx_name","gene_id","cds_len") %in% names(tx_tbl)))
  stopifnot(nrow(tx_tbl) == nrow(mix_tx_pos))

  # Keep only CDS-translatable (>0 and divisible by 3) transcript candidates.
  tx_ok <- tx_tbl[tx_tbl$cds_len > 0L & (tx_tbl$cds_len %% 3L == 0L),
                  c("tx_name","gene_id")]
  tx_ok$tx_idx <- match(tx_ok$tx_name, tx_tbl$tx_name)

  # For each read in df_map, get the corresponding CDS interval on mix_tx coordinates.
  # cds_start = utr5_p3 + 1; cds_end = utr3_p5 - 1
  tx_idx_map <- match(df_map$tx_name, tx_tbl$tx_name)
  cds_start <- mix_tx_pos$utr5_p3[tx_idx_map] + 1L
  cds_end   <- mix_tx_pos$utr3_p5[tx_idx_map] - 1L

  read_start <- df_map$pos
  read_end   <- df_map$pos + df_map$qwidth - 1L

  in_cds <- !is.na(cds_start) & !is.na(cds_end) &
            (read_end >= cds_start) & (read_start <= cds_end)

  df_cds <- df_map[in_cds, c("gene_id","tx_name")]
  if (nrow(df_cds) == 0L) {
    stop("No reads overlap CDS. Check whether df_map coordinates match mix_tx_pos.")
  }

  # Filter: keep only reads from transcript candidates (CDS%%3==0).
  df_cds <- dplyr::inner_join(df_cds, tx_ok[, c("tx_name","gene_id")],
                              by = c("tx_name","gene_id"))

  # Gene-level counts (CDS-overlap based).
  gene_counts <- df_cds |>
    dplyr::count(gene_id, name = "rpf_n") |>
    dplyr::filter(rpf_n >= min_reads)

  if (nrow(gene_counts) < n_genes) {
    stop("Only ", nrow(gene_counts),
         " genes meet: >= ", min_reads, " CDS-overlap reads. Need ", n_genes, ".")
  }

  set.seed(20251110)
  picked_genes <- sample(gene_counts$gene_id, n_genes)

  # Select representative transcript per picked gene by maximal CDS-overlap read count.
  tx_choice <- df_cds |>
    dplyr::filter(gene_id %in% picked_genes) |>
    dplyr::count(gene_id, tx_name, name = "tx_reads") |>
    dplyr::group_by(gene_id) |>
    dplyr::slice_max(tx_reads, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  picked_info <- gene_counts |>
    dplyr::filter(gene_id %in% picked_genes) |>
    dplyr::inner_join(tx_choice, by = "gene_id") |>
    dplyr::select(gene_id, tx_name, rpf_n)

  list(picked_genes = picked_genes, picked_info = picked_info)
}


## ---- Apply +/- 3 nt jitter within the same transcript interval ----
jitter_within_tx <- function(pos, qwidth, tx_start, tx_end) {
  off <- sample(c(-3L,-2L,-1L,0L,1L,2L,3L), length(pos), replace = TRUE)
  p2  <- pos + off
  pmin(pmax(p2, tx_start), tx_end - qwidth + 1L)
}

## ---- Export FASTQ from concatenated reference (chunked writing) ----
write_fastq_from_mix <- function(df, mix_tx, out_path, prefix, chunk_size = 5e6) {
  n <- nrow(df); if (n == 0) return(invisible(0L))
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  for (k in seq_along(idx)) {
    ii <- idx[[k]]
    seqs <- as.character(subseq(DNAStringSet(mix_tx)[rep(1,length(df$pos[ii]))], start = df$pos[ii], width = df$qwidth[ii]))
    is_add <- df$class[ii] == "add5" & !is.na(df$seq5[ii])
    if (any(is_add)) substr(seqs[is_add], 1, 1) <- df$seq5[ii][is_add]

    ids <- sprintf("%s|gene:%s|tx:%s|pos:%d|len:%d|cls:%s",
                   prefix, df$gene_id[ii], df$tx_name[ii],
                   df$pos[ii], df$qwidth[ii], df$class[ii])

    quals <- vapply(df$qwidth[ii], function(L) paste0(rep.int("I", L), collapse = ""), "")
    fq <- ShortReadQ(sread = DNAStringSet(seqs),
                     quality = BStringSet(quals),
                     id = BStringSet(ids))
    writeFastq(fq, out_path, compress = FALSE, mode = if (k == 1) "w" else "a")
  }
  invisible(nrow(df))
}

## ---- Translate AA sequence using representative transcripts (CDS%%3==0) ----
get_tx_aa <- function(tx_name, tx_tbl, tx_seqs) {
  row <- tx_tbl[match(tx_name, tx_tbl$tx_name), ]
  if (is.na(row$utr5_len) || is.na(row$cds_len) || row$cds_len <= 0L) return(NA_character_)
  if (row$cds_len %% 3L != 0L) return(NA_character_)  # Safety check.
  dna <- tx_seqs[[tx_name]]
  if (is.null(dna)) return(NA_character_)
  cds_start_in_tx <- as.integer(row$utr5_len) + 1L
  cds_end_in_tx   <- cds_start_in_tx + as.integer(row$cds_len) - 1L
  cds <- subseq(dna, start = cds_start_in_tx, end = cds_end_in_tx)
  as.character(translate(cds, if.fuzzy.codon = "error"))
}

## ===================== Main Workflow =====================
process_all_samples <- function(reads_bam3, tx_info,
                                out_dir = ".", n_genes = 1000, min_reads = 10L) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tx_idx  <- build_tx_table(tx_info)
  tx_tbl  <- tx_idx$tx_tbl
  tx_gr   <- tx_idx$tx_gr
  tx_seqs <- tx_info$flt_tx_seqs
  mix_tx  <- tx_info$mix_tx

  sample_keys <- names(reads_bam3)

  for (sk in sample_keys) {
    cat("\n========== Processing sample:", sk, "==========\n")
    df_raw <- flatten_sample_reads(reads_bam3[[sk]])
    df_map <- map_reads_to_tx(df_raw, tx_gr)
    cat("Mapped reads:", nrow(df_map), "\n")

    picked <- pick_genes_and_tx_cds(df_map, tx_tbl, mix_tx_pos = tx_info$mix_tx_pos, n_genes = n_genes, min_reads = min_reads)
    picked_genes <- picked$picked_genes
    picked_info  <- picked$picked_info

    # POS = original positions
    pos_df <- df_map

    # NEG = jitter only reads from selected 1000 genes, constrained to same transcript interval
    neg_df <- df_map
    sel_idx <- which(neg_df$gene_id %in% picked_genes)
    if (length(sel_idx) > 0) {
      neg_df$pos[sel_idx] <- jitter_within_tx(
        pos      = neg_df$pos[sel_idx],
        qwidth   = neg_df$qwidth[sel_idx],
        tx_start = neg_df$tx_start[sel_idx],
        tx_end   = neg_df$tx_end[sel_idx]
      )
    }

    # Export FASTQ
    pos_fq <- file.path(out_dir, paste0("POS_", sk, ".fastq"))
    neg_fq <- file.path(out_dir, paste0("NEG_", sk, ".fastq"))
    cat("Writing:", pos_fq, "\n")
    n_pos <- write_fastq_from_mix(pos_df, mix_tx, pos_fq, prefix = paste0("POS|", sk))
    cat("Writing:", neg_fq, "\n")
    n_neg <- write_fastq_from_mix(neg_df, mix_tx, neg_fq, prefix = paste0("NEG|", sk))
    cat("Done FASTQ. POS reads:", n_pos, " NEG reads:", n_neg, "\n")

    # AA sequences for representative transcripts (from cds_len%%3==0 transcripts only)
    picked_info <- picked_info %>%
      mutate(aa_seq = vapply(tx_name, get_tx_aa, character(1),
                             tx_tbl = tx_tbl, tx_seqs = tx_seqs))

    out_tsv <- file.path(out_dir, paste0("picked_1000_genes_", sk, ".tsv"))
    write_tsv(picked_info, out_tsv)
    cat("Wrote:", out_tsv, " (", nrow(picked_info), " rows)\n", sep = "")
  }
  invisible(TRUE)
}

## ===================== Example Call =====================
process_all_samples(reads_bam3, tx_info, 
out_dir = SIM_DIR, 
n_genes = 1000, min_reads = 20L)

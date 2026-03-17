pred_p <- function(par_rpf,
                   tx_info,
                   iter_num) {
  # iterate n times to reach balance
  pred_psite <- function(df_rpf,
                         candi_psite,
                         candi_cut5,
                         candi_cut3,
                         tx_info,
                         ribo_size,
                         par_0,
                         iter_num) {
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
    )
    idx_ij <- match(candi_psite@x, uni_candi_psite)
    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq, bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
    )
    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
    prob_cum <- vector(mode = "numeric", length = iter_num)

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
      tmp_sum <- Matrix::rowSums(candi_p_weight)
      prob_cum[i] <- sum(df_rpf$weight * log(tmp_sum))
      sparse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / tmp_sum)@x
      )
      vec_pnum[] <- Matrix::rowSums(sparse_iter)
    }
    return(
      list(
        vec_pnum = vec_pnum,
        candi_pos = candi_pos,
        prob_cum = prob_cum
      )
    )
  }

  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[
    par_rpf$lst_rpf$rpf_overlap$qwidth %in% unique(convert_idx[4, ]),
  ]
  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = psite_pos
  )
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  ribo_num_0 <- pred_psite(
    df_rpf = par_rpf$lst_rpf$rpf_overlap,
    candi_psite = candi_psite,
    candi_cut5 = candi_cut5,
    candi_cut3 = candi_cut3,
    tx_info = tx_info,
    ribo_size = ribo_size,
    par_0 = par_rpf$par_0,
    iter_num = iter_num
  )
  return(ribo_num_0)
}

get_cleavage_seq <- function(seqs,
                             p_site,
                             ribo_size) {
  seqs <- Biostrings::DNAStringSet(x = seqs)[
    rep(1L, length(p_site))
  ]

  up_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site - sum(ribo_size[1:2]),
    width = ribo_size[1]
  )

  dn_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site + ribo_size[3],
    width = ribo_size[4]
  )

  return(
    list(
      up_seq = up_seq,
      dn_seq = dn_seq
    )
  )
}

maintain_prob <- function(cut_seq,
                          prod_hd,
                          bias) {
  prob_mc <- stringr::str_split(
    string = as.vector(cut_seq),
    pattern = "",
    simplify = TRUE
  )

  prob_mp <- matrix(
    data = bias[prob_mc],
    nrow = nrow(prob_mc)
  )

  prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

  return(prob_mp)
}

get_cleavage_prob <- function(seqs,
                              bias,
                              prob_hd5,
                              prob_hd3) {
  # for 5'
  maintain_prob5 <- maintain_prob(
    cut_seq = seqs$up_seq,
    prod_hd = prob_hd5,
    bias = bias
  )

  cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

  maintain_cumprod5rev <- matrix(
    data = 1,
    nrow = nrow(maintain_prob5),
    ncol = length(prob_hd5)
  )

  maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
    maintain_prob5[, length(prob_hd5):2]
  )

  final_p5 <- maintain_cumprod5rev * cle_p5rev

  final_p5 <- final_p5[, ncol(final_p5):1]

  final_p5 <- final_p5 / Matrix::rowSums(final_p5)

  # for 3'
  maintain_prob3 <- maintain_prob(
    cut_seq = seqs$dn_seq,
    prod_hd = prob_hd3,
    bias = bias
  )

  cle_p3 <- 1 - maintain_prob3

  maintain_cumprod3 <- matrix(
    data = 1,
    nrow = nrow(maintain_prob3),
    ncol = length(prob_hd3)
  )

  maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
    maintain_prob3[, -length(prob_hd3)]
  )

  final_p3 <- maintain_cumprod3 * cle_p3

  final_p3 <- final_p3 / Matrix::rowSums(final_p3)

  return(
    list(
      final_p5 = final_p5,
      final_p3 = final_p3
    )
  )
}

simulate_par0 <- function(
    par_5add,
    par_lig,
    rnase,
    size_pos,
    steel5,
    steel3) {
prob_hinder <- function(pos, steel5, steel3) {
      shiftlen <- c(2, 3)
      pos <- pos[-c(1:shiftlen[1])]
      mpos <- median(pos) + 1
      pos_change <- (pos[1] - mpos):(pos[length(pos)] - mpos + shiftlen[1])
      sigmoid <- 1 / (1 + exp(-steel5 * pos_change))
      p5 <- rev(sigmoid)
      pos <- pos[-c(1:(shiftlen[2] - shiftlen[1]))]
      mpos <- mpos + 1
      pos_change <- (pos[1] - mpos):(pos[length(pos)] - mpos + shiftlen[2])
      sigmoid <- 1 / (1 + exp(-steel3 * pos_change))
      p3 <- sigmoid
      return(list(p5 = p5, p3 = p3))
    }

  if (par_5add == T) {
    par_5add <- c(0.03, 0.04, 0.03, 0.2, 0.7)
  } else {
    par_5add <- c(0.001, 0.002, 0.0003, 0.0004, 1)
  }
  prob_add5 <- par_5add
  names(prob_add5) <- c("A", "C", "G", "T", "")

  if (par_lig == T) {
    eff_intercept <- -0.5
    eff_f5 <- runif(64, 0, 0.25)
    eff_f3 <- runif(64, 0, 0.25)
  } else {
    eff_intercept <- 0
    eff_f5 <- rep(0, 64)
    eff_f3 <- rep(0, 64)
  }
  names(eff_f5) <- names(eff_f3) <- names(Biostrings::GENETIC_CODE)
tmp_p = prob_hinder(1:size_pos, steel5, steel3)
  if (rnase == "rnase-i") {
    cut_bias <- rep(0.1, 4)
    names(cut_bias) <- c("A", "C", "G", "T")
    ribo_size <- c(size_pos, 9L, 15L, size_pos)
  } else if (rnase == "p1") {
    cut_bias <- rep(0.1, 4)
    names(cut_bias) <- c("A", "C", "G", "T")
    ribo_size <- c(size_pos, 10L, 18L, size_pos)
  } else {
    cut_bias <- c(0.003324, 0.72334, 0.6342, 0.1000)
    names(cut_bias) <- c("A", "C", "G", "T")
    ribo_size <- c(size_pos, 9L, 15L, size_pos)
  }
  par_0 <- list(
    prob_add5 = prob_add5,
    ribo_size = ribo_size,
    prob_hd = list(
      p5 = tmp_p$p5,
      p3 = tmp_p$p3
    ),
    base_prob = 0.1,
    cut_bias = list(
      s7 = cut_bias
    ),
    eff_intercept = eff_intercept,
    eff_f5 = eff_f5,
    eff_f3 = eff_f3
  )
  return(par_0)
}



# 已知frame预测p位点
pred_p_frame <- function(par_rpf, tx_info, iter_num) {
  # iterate n times to reach balance
  pred_psite <- function(df_rpf,
                         candi_psite,
                         candi_cut5,
                         candi_cut3,
                         tx_info,
                         ribo_size,
                         par_0,
                         iter_num) {
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
    )
    idx_ij <- match(candi_psite@x, uni_candi_psite)
    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq, bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
    )
    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
    prob_cum <- vector(mode = "numeric", length = iter_num)

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
      tmp_sum <- Matrix::rowSums(candi_p_weight)
      prob_cum[i] <- sum(df_rpf$weight * log(tmp_sum))
      sparse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / tmp_sum)@x
      )
      vec_pnum[] <- Matrix::rowSums(sparse_iter)
    }
    return(
      list(
        vec_pnum = vec_pnum,
        candi_pos = candi_pos,
        prob_cum = prob_cum
      )
    )
  }

  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  # frame
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  rg_psite <- IRanges::IRanges(start = psite_pos, width = 1L)
  hit_psite <- IRanges::findOverlaps(rg_psite, rg_cds)
  idx_p <- hit_psite@from[
    ((psite_pos[hit_psite@from] - rg_cds@start[hit_psite@to]) %% 3L) == 1L
  ]
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx)[idx_p],
    j = sequence(nvec = psite_num_idx, from = 1L)[idx_p],
    x = psite_pos[idx_p]
  )
  judg_v <- vector(length = length(psite_pos))
  judg_v[idx_p] <- TRUE
  candi_psite1 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = judg_v
  )
  idx_psites1 <- Matrix::rowSums(candi_psite1) > 0L
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx)[idx_p],
    j = sequence(nvec = psite_num_idx, from = 1L)[idx_p],
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )[idx_p]
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx)[idx_p],
    j = sequence(nvec = psite_num_idx, from = 1L)[idx_p],
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )[idx_p]
  )
  idx_psites2 <- Matrix::rowSums(candi_psite) > 0
  ribo_num_0 <- pred_psite(
    df_rpf = par_rpf$lst_rpf$rpf_overlap[idx_psites1, ],
    candi_psite = candi_psite[idx_psites2, ],
    candi_cut5 = candi_cut5[idx_psites2, ],
    candi_cut3 = candi_cut3[idx_psites2, ],
    tx_info = tx_info,
    ribo_size = ribo_size,
    par_0 = par_rpf$par_0,
    iter_num = iter_num
  )
  return(ribo_num_0)
}


# 模拟数据准备rpf信号
prep_rpf_simu <- function(
    lst_simu,
    tx_info,
    add5,
    prob_add5) {
  if (!add5) {
    rpf_info <- data.frame(
      pos = lst_simu$reads_add5_pos,
      qwidth = BiocGenerics::width(lst_simu$reads_add5)
    )

    offsets <- estimate_offsets(
      tx_info = tx_info,
      rpf_info = rpf_info,
      number_min = 0.99,
      max_frame_min = 0.3,
      choose_terminal = "start",
      limited_range = 10:18
    )

    rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]

    rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
  } else {
    if (prob_add5[5] == 1) {
      rpf_info <- data.frame(
        pos = lst_simu$reads_add5_pos - BiocGenerics::width(lst_simu$add_base),
        qwidth = width(lst_simu$reads_add5)
      )

      ref_seqs <- Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, nrow(rpf_info))],
        start = rpf_info$pos,
        width = rpf_info$qwidth
      )

      idx_add <- lst_simu$reads_add5 == ref_seqs

      rpf_info_noadd <- rpf_info[idx_add, ]
      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info_noadd,
        number_min = 0.99,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )

      rpf_info <- rpf_info_noadd[rpf_info_noadd$qwidth %in% offsets[1, ], ]

      rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
    } else {
      rpf_info <- data.frame(
        pos = lst_simu$reads_add5_pos - BiocGenerics::width(lst_simu$add_base),
        qwidth = width(lst_simu$reads_add5)
      )

      ref_seqs <- Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, nrow(rpf_info))],
        start = rpf_info$pos,
        width = rpf_info$qwidth
      )

      idx_add <- lst_simu$reads_add5 == ref_seqs

      rpf_info_add <- rpf_info[!idx_add, ]
      rpf_info_noadd <- rpf_info[idx_add, ]

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = 0.99,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )

      rpf_info_add <- rpf_info_add[rpf_info_add$qwidth %in% offsets[1, ], ]
      rpf_info_noadd <- rpf_info_noadd[
        rpf_info_noadd$qwidth %in% offsets[1, ],
      ]

      read_tag_add <- base::table(
        rpf_info_add$pos * 100 + rpf_info_add$qwidth + 99
      )

      read_num_tag_add <- as.numeric(names(read_tag_add))

      read_weight_add <- data.frame(
        pos = as.integer(read_num_tag_add %/% 100),
        qwidth = as.integer(read_num_tag_add %% 100),
        weight = as.vector(read_tag_add),
        tag = read_num_tag_add
      )

      ref_base <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, nrow(read_weight_add))],
          start = read_weight_add$pos - 1,
          width = 1
        )
      )

      tr_base <- sapply(1:4, function(x) {
        1 + prob_add5[x] / sum(prob_add5[1:4][-x])
      })

      read_weight_add$weight <- read_weight_add$weight * tr_base[ref_base]

      read_tag <- base::table(
        rpf_info_noadd$pos * 100 + rpf_info_noadd$qwidth
      )

      read_num_tag <- as.numeric(names(read_tag))

      read_weight <- data.frame(
        pos = as.integer(read_num_tag %/% 100),
        qwidth = as.integer(read_num_tag %% 100),
        weight = as.vector(read_tag),
        tag = read_num_tag
      )

      ref_base1 <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, length(read_num_tag))],
          start = read_num_tag %/% 100,
          width = 1
        )
      )

      tr_base2 <- sapply(1:4, function(x) {
        1 - prob_add5[x] /
          (prob_add5[5] +
            prob_add5[x])
      })

      read_weight$weight <- tr_base2[ref_base1] * read_weight$weight

      rpf_overlap <- rbind(read_weight, read_weight_add)
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]
    }
  }
  return(list(
    rpf_info = rpf_info,
    rpf_overlap = rpf_overlap,
    offsets = offsets
  ))
}

# 给模拟的p位点上核糖体数目产生reads
simulate_reads <- function(
    par_setup,
    tx_info,
    p_vec,
    add_noise,
    noise_ratio) {

  # Separately generate sequences preceding cleavage at both ends of
  # the ribosome
  fun_cut_seq <- function(seqs,
                          p_pos,
                          ribo_size) {
    p_site <- which(p_pos > 0L)
    seqs <- Biostrings::DNAStringSet(seqs)[rep(1L, length(p_site))]
    up_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site - sum(ribo_size[1:2]),
      width = ribo_size[1]
    )
    dn_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site + ribo_size[3],
      width = ribo_size[4]
    )
    return(
      list(
        up_seq = up_seq,
        dn_seq = dn_seq
      )
    )
  }
  # Generating cleavage probabilities for each ribosome terminus
  fun_cut_prob <- function(seqs,
                           bias,
                           prob_hd5,
                           prob_hd3) {
    maintain_prob <- function(cut_seq,
                              prod_hd) {
      prob_mc <- stringr::str_split(
        as.vector(cut_seq),
        pattern = "",
        simplify = TRUE
      )
      prob_mp <- matrix(
        data = bias[prob_mc],
        nrow = nrow(prob_mc)
      )
      prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)
      return(prob_mp)
    }
    # for 5'
    maintain_prob5 <- maintain_prob(
      cut_seq = seqs$up_seq,
      prod_hd = prob_hd5
    )
    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]
    maintain_cumprod5rev <- matrix(
      1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )
    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      maintain_prob5[, length(prob_hd5):2]
    )
    final_p5 <- maintain_cumprod5rev * cle_p5rev
    final_p5 <- cbind(final_p5, 1 - Matrix::rowSums(final_p5))
    final_p5 <- final_p5[, ncol(final_p5):1]

    # for 3'
    maintain_prob3 <- maintain_prob(
      cut_seq = seqs$dn_seq,
      prod_hd = prob_hd3
    )
    cle_p3 <- 1 - maintain_prob3
    maintain_cumprod3 <- matrix(
      1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )
    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      maintain_prob3[, -length(prob_hd3)]
    )
    final_p3 <- maintain_cumprod3 * cle_p3
    final_p3 <- cbind(final_p3, 1 - Matrix::rowSums(final_p3))

    return(list(final_p5 = final_p5, final_p3 = final_p3))
  }
  # Integrate various parameters to generate simulated reads
  fun_create_read <- function(seqs,
                              cut_exp,
                              p_pos,
                              ribo_size) {
    p_site <- which(p_pos > 0L)
    p_pos_expand <- rep(p_site, p_pos[p_site])
    expand_m5 <- matrix(
      rep(
        cut_exp$final_p5,
        rep(
          p_pos[p_site],
          ncol(cut_exp$final_p5)
        )
      ),
      ncol = ncol(cut_exp$final_p5)
    )
    expand_m3 <- matrix(
      rep(
        cut_exp$final_p3,
        rep(
          p_pos[p_site], ncol(cut_exp$final_p3)
        )
      ),
      ncol = ncol(cut_exp$final_p3)
    )
    pos5idx <- as.integer(
      rowSums(
        runif(nrow(expand_m5)) > matrixStats::rowCumsums(expand_m5)
      )
    )
    pos3idx <- as.integer(
      rowSums(
        runif(nrow(expand_m3)) > matrixStats::rowCumsums(expand_m3)
      )
    )
    out_reg <- (pos5idx != 0L) & (pos3idx != ribo_size[4])
    cut_pos5 <- pos5idx + p_pos_expand - sum(ribo_size[1:2]) - 1L
    cut_pos3 <- pos3idx + p_pos_expand + ribo_size[3]

    p_correct <- as.vector(
      base::table(
        c(p_pos_expand[out_reg], unique(p_pos_expand))
      ) - 1L
    )
    p_original <- as.vector(base::table(p_pos_expand))

    reads <- Biostrings::subseq(
      Biostrings::DNAStringSet(seqs)[rep(1L, length(cut_pos5))],
      start = cut_pos5,
      end = cut_pos3 - 1L
    )

    return(list(
      reads = reads[out_reg],
      p_correct = p_correct,
      p_pos = p_site,
      p_original = p_original
    ))
  }

  cut_seq <- fun_cut_seq(
    seqs = tx_info$mix_tx,
    p_pos = p_vec,
    ribo_size = par_setup$ribo_size
  )
  cut_prob <- fun_cut_prob(
    seqs = cut_seq,
    bias = par_setup$cut_bias$s7,
    prob_hd5 = par_setup$prob_hd$p5,
    prob_hd3 = par_setup$prob_hd$p3
  )
  expected_cut_prob = list(
    p5=colMeans(cut_prob$final_p5),
  p3 = colMeans(cut_prob$final_p3))
  # create_read$p_original 与 create_read$p_correct差别就是p_original是实际reads数量，两者差别大只是证明参数选择不好，切割位点超出ribosize设置了
  create_read <- fun_create_read(
    seqs = tx_info$mix_tx,
    cut_exp = cut_prob,
    p_pos = p_vec,
    ribo_size = par_setup$ribo_size
  )
  # Ligase Ligation Probability
  kmer5 <- Biostrings::subseq(
    create_read$reads,
    start = 1,
    width = 3
  )
  kmer3 <- Biostrings::subseq(
    create_read$reads,
    end = BiocGenerics::width(create_read$reads),
    width = 3
  )
  ligation_prob <- exp(
    par_setup$eff_intercept +
      par_setup$eff_f5[as.vector(kmer5)] +
      par_setup$eff_f3[as.vector(kmer3)]
  )
  idx_ligat <- runif(
    length(create_read$reads)
  ) < ligation_prob
  create_read$reads_lig <- create_read$reads[idx_ligat]

# A single base is added to the 5' end， 独立概率
  add_base <- sample(
    names(par_setup$prob_add5),
    length(create_read$reads_lig),
    replace = TRUE,
    prob = par_setup$prob_add5
  )
  create_read$add_base <- Biostrings::DNAStringSet(add_base)
  create_read$reads_add5 <- Biostrings::xscat(
    Biostrings::DNAStringSet(add_base),
    create_read$reads_lig
  )
  create_read$reads_add5_pos <- create_read$reads_lig@ranges@start
print('go')
  # 噪音加入，修改类似于RNA结合蛋白的噪音
  if (add_noise == TRUE) {
    n_pos <- create_read$reads_add5_pos +
      sample(
        -30:30,
        length(create_read$reads_add5_pos),
        replace = TRUE
      )
    n_width <- sample(
      25:34,
      length(create_read$reads_add5_pos),
      replace = TRUE,
      prob = c(1:5, 5:1)
    )
    tmp_pos <- sample(
      n_pos,
      size = noise_ratio * length(n_pos)
    )
    create_read$reads_add5_pos <- c(create_read$reads_add5_pos, tmp_pos)

    tmp_idx <- sample(1:length(n_pos), noise_ratio * length(n_pos))
    create_read$reads_add5 <- c(
      create_read$reads_add5,
      Biostrings::subseq(
        Biostrings::DNAStringSet(tx_info$mix_tx)[
          rep(1L, length(create_read$reads_add5_pos))
        ],
        start = n_pos,
        width = n_width
      )[tmp_idx]
    )
  }
create_read$expected_cut_prob = expected_cut_prob
  return(create_read)
}
# 给样本p位点数量，产生模拟的p位点上核糖体数目并产生reads
create_simu_reads <- function(
    par_0,
    ribo_num,
    tx_info,
    sample_num,
    library_size,
    noise_ratio) {
  psite_input <- rmultinom(
    n = sample_num,
    size = library_size,
    prob = ribo_num
  )

  lst_simulated_reads <- lapply(
    seq.int(sample_num), function(i) {
      simulated_reads <- simulate_reads(
        par_setup = par_0,
        tx_info = tx_info,
        p_vec = psite_input[, i],
        add_noise = TRUE,
        noise_ratio = noise_ratio
      )
      return(simulated_reads)
    }
  )
  return(lst_simulated_reads)
}

# 非压缩、分块写 FASTQ；dna 为 Biostrings::DNAStringSet
write_fastq_chunked <- function(dna,
                                file,
                                chunk_size = 250000L,  # 每块 25 万条，内存友好
                                phred = 40L,           # 统一质量（PHRED33）
                                id_prefix = "read") {
  stopifnot(inherits(dna, "DNAStringSet"))
  n <- length(dna)
  if (n == 0L) stop("No sequences to write.")
  # 首块用覆盖写，其余用追加
  if (file.exists(file)) file.remove(file)

  # 质量字符（PHRED33: 33 + phred）
  phred <- as.integer(phred)
  if (phred < 0L || phred > 60L) stop("`phred` must be in 0..60")
  qchar <- intToUtf8(33L + phred)

  starts <- seq.int(1L, n, by = as.integer(chunk_size))
  for (s in starts) {
    e <- min(s + chunk_size - 1L, n)
    dna_chunk <- dna[s:e]

    # 生成/填充 ID
    nm <- names(dna_chunk)
    if (is.null(nm)) {
      nm <- sprintf("%s_%06d", id_prefix, seq.int(s, e))
    } else {
      empty <- which(is.na(nm) | nm == "")
      if (length(empty)) {
        nm[empty] <- sprintf("%s_%06d", id_prefix, s - 1L + empty)
      }
    }
    ids <- Biostrings::BStringSet(nm)

    # 质量串（每条 read 一整行相同字符）
    wv <- Biostrings::width(dna_chunk)
    quals <- Biostrings::BStringSet(vapply(wv, function(k) strrep(qchar, k), ""))

    fq <- ShortRead::ShortReadQ(
      sread   = dna_chunk,
      quality = quals,
      id      = ids
    )

    # 非压缩 + 追加写
    # ShortRead::writeFastq 支持 mode = "w"/"a"
    ShortRead::writeFastq(fq, file = file, compress = FALSE,
                          mode = if (s == 1L) "w" else "a")
  }
  invisible(file)
}

infer_par <- function(tx_info,
                      lst_bam,
                      rnase_bias,
                      RNase,
                      ligat_par) {
                        #browser()
  correct_kmer <- function(rpf,
                           tx_info,
                           par_0) {
    correct_rpf <- function(kmer5,
                            kmer3,
                            rpf_num,
                            par_0) {
      corrected_rpf <- exp(
        (rpf_num -
          par_0$eff_intercept -
          par_0$eff_f5[kmer5] -
          par_0$eff_f3[kmer3]
        )
      ) - 1
      return(as.vector(corrected_rpf))
    }

    # 3-mer at the terminal of reads
    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf$pos + rpf$qwidth - 1L,
      width = 3L
    ))

    corrected_rpf <- correct_rpf(
      kmer5 = kmer_5,
      kmer3 = kmer_3,
      rpf_num = log(rpf$weight + 1),
      par_0 = par_0
    )

    corrected_rpf[which(corrected_rpf < 0)] <- 0

    if (par_0$eff_intercept != 0) {
      corrected_rpf <- corrected_rpf / (sum(corrected_rpf) / sum(rpf$weight))
    }

    return(corrected_rpf)
  }

  glm2par <- function(fit_coefs,
                      par_0) {
    par_0$eff_intercept[1] <- fit_coefs[1]

    par_0$eff_f5[
      stringr::str_split(
        string = names(fit_coefs)[
          grep(
            pattern = "kmer5",
            x = names(fit_coefs)
          )
        ],
        pattern = "5",
        simplify = TRUE
      )[, 2]
    ] <- fit_coefs[2:64]

    par_0$eff_f3[
      stringr::str_split(
        string = names(fit_coefs)[
          grep(
            pattern = "kmer3",
            x = names(fit_coefs)
          )
        ],
        pattern = "3",
        simplify = TRUE
      )[, 2]
    ] <- fit_coefs[65:127]

    return(
      list(
        par_0 = par_0
      )
    )
  }

  iter_rpf_overlap <- function(potential_rpf_lig2,
                               rpf_overlap) {
    rpf_overlap_iter <- potential_rpf_lig2[potential_rpf_lig2$weight > 0, ]
    rpf_overlap_iter$tag <- rpf_overlap_iter$pos * 100 + rpf_overlap_iter$qwidth
    iter_idx <- match(rpf_overlap_iter$tag, rpf_overlap$tag)
    rpf_overlap_iter$weight <- rpf_overlap$weight[iter_idx]
    return(list(
      iter_idx = iter_idx,
      rpf_overlap_iter = rpf_overlap_iter
    ))
  }

  maintain_prob <- function(cut_seq,
                            prod_hd,
                            bias) {
    prob_mc <- stringr::str_split(
      string = as.vector(cut_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mp <- matrix(
      data = bias[prob_mc],
      nrow = nrow(prob_mc)
    )

    prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

    return(prob_mp)
  }

  get_cleavage_prob <- function(seqs,
                                bias,
                                prob_hd5,
                                prob_hd3) {
    # for 5'
    maintain_prob5 <- maintain_prob(
      cut_seq = seqs$up_seq,
      prod_hd = prob_hd5,
      bias = bias
    )

    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

    maintain_cumprod5rev <- matrix(
      data = 1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )

    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      maintain_prob5[, length(prob_hd5):2]
    )

    final_p5 <- maintain_cumprod5rev * cle_p5rev

    final_p5 <- final_p5[, ncol(final_p5):1]

    final_p5 <- final_p5 / Matrix::rowSums(final_p5)

    # for 3'
    maintain_prob3 <- maintain_prob(
      cut_seq = seqs$dn_seq,
      prod_hd = prob_hd3,
      bias = bias
    )

    cle_p3 <- 1 - maintain_prob3

    maintain_cumprod3 <- matrix(
      data = 1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )

    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      maintain_prob3[, -length(prob_hd3)]
    )

    final_p3 <- maintain_cumprod3 * cle_p3

    final_p3 <- final_p3 / Matrix::rowSums(final_p3)

    return(
      list(
        final_p5 = final_p5,
        final_p3 = final_p3
      )
    )
  }

  get_cleavage_seq <- function(seqs,
                               p_site,
                               ribo_size) {
    seqs <- Biostrings::DNAStringSet(x = seqs)[
      rep(1L, length(p_site))
    ]

    up_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site - sum(ribo_size[1:2]),
      width = ribo_size[1]
    )

    dn_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site + ribo_size[3],
      width = ribo_size[4]
    )

    return(
      list(
        up_seq = up_seq,
        dn_seq = dn_seq
      )
    )
  }

  em_ribo_num <- function(df_rpf,
                          candi_psite,
                          candi_cut5,
                          candi_cut3,
                          tx_info,
                          ribo_size,
                          par_0,
                          pos_psite,
                          iter_num) {
    # filter out read weight 0
    idx_read <- df_rpf$weight > 0L

    df_rpf <- df_rpf[idx_read, ]

    candi_p_weight <- candi_psite <- candi_psite[idx_read, ]

    candi_cut5 <- candi_cut5[idx_read, ]

    candi_cut3 <- candi_cut3[idx_read, ]

    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)

    shrink_pos <- as.integer(shrink_p)

    index_j <- rank(
      x = shrink_pos,
      ties.method = "first"
    ) -
      rank(
        x = shrink_pos,
        ties.method = "min"
      ) +
      1L

    candi_pos <- as.integer(levels(shrink_p))

    vec_pnum <- rep(1, length(candi_pos))

    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = ribo_size
    )

    idx_ij <- match(candi_psite@x, uni_candi_psite)

    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) *
        (candi_cut5@x - 1L) +
        idx_ij
    ] *
      cut_prob$final_p3[
        nrow(cut_prob$final_p3) *
          (candi_cut3@x - 1L) +
          idx_ij
      ]

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]

      rse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / Matrix::rowSums(candi_p_weight))@x
      )

      vec_pnum[] <- Matrix::rowSums(rse_iter)
    }

    idx_psite <- match(pos_psite, candi_pos)

    vec_pnum <- vec_pnum[idx_psite]

    vec_pnum[is.na(vec_pnum)] <- 0

    return(vec_pnum)
  }

  adjust_ribo <- function(ribo_size,
                          prob_hd,
                          limit1 = 0.05,
                          limit2 = 0.9,
                          limit3 = 0.9) {
    ribo_0 <- ribo_size
    prob_hd <- prob_hd
    p5len <- ribo_size[1]
    p3len <- ribo_size[4]
    chg <- 0
    if (chg == 0) {
      if (prob_hd$p3[1] < limit1) {
        if (prob_hd$p3[2] < limit1) {
          ribo_size[4] <- ribo_size[4] - 1L
          ribo_size[3] <- ribo_size[3] + 1L
          chg <- 1
        }
      } else {
        if (min(prob_hd$p3) == prob_hd$p3[1]) {
          ribo_size[4] <- ribo_size[4] + 1L
          ribo_size[3] <- ribo_size[3] - 1L
          chg <- 1
        }
      }
    }
    if (chg == 0) {
      if (prob_hd$p5[p5len] < limit1) {
        if (prob_hd$p5[p5len - 1] < limit1) {
          ribo_size[1] <- ribo_size[1] - 1L
          ribo_size[2] <- ribo_size[2] + 1L
          chg <- 1
        }
      } else {
        if (min(prob_hd$p5) == prob_hd$p5[p5len]) {
          ribo_size[1] <- ribo_size[1] + 1L
          ribo_size[2] <- ribo_size[2] - 1L
          chg <- 1
        }
      }
    }
    if (chg == 0) {
      if (prob_hd$p3[p3len] > limit2) {
        if (prob_hd$p3[p3len - 1] > limit3) {
          ribo_size[4] <- ribo_size[4] - 1L
          chg <- 1
        }
      } else {
        if (max(prob_hd$p3) != prob_hd$p3[p3len]) {
          ribo_size[4] <- ribo_size[4] - 1L
          chg <- 1
        } else {
          ribo_size[4] <- ribo_size[4] + 1L
          chg <- 1
        }
      }
    }

    if (chg == 0) {
      if (prob_hd$p5[1] > limit2) {
        if (prob_hd$p5[2] > limit3) {
          ribo_size[1] <- ribo_size[1] - 1L
          chg <- 1
        }
      } else {
        if (max(prob_hd$p5) != prob_hd$p5[1]) {
          ribo_size[1] <- ribo_size[1] - 1L
          chg <- 1
        } else {
          ribo_size[1] <- ribo_size[1] + 1L
          chg <- 1
        }
      }
    }

    if (identical(ribo_0, ribo_size)) {
      chg_ribo <- FALSE
    } else {
      (
        chg_ribo <- TRUE
      )
    }
    return(list(ribo_size = ribo_size, chg_ribo = chg_ribo))
  }

  assign_par <- function(out_par,
                         par_0,
                         reference_base_value) {
    bias_len <- length(par_0$cut_bias$s7)

    hindrance_len <- length(par_0$prob_hd$p5)

    out_par <- c(
      out_par[1:(bias_len - 1)],
      reference_base_value,
      out_par[-(1:(bias_len - 1))]
    )

    par_0$cut_bias$s7[
      1:bias_len
    ] <- out_par[1:bias_len]

    par_0$prob_hd$p5 <- out_par[
      (bias_len + 1):(bias_len + hindrance_len)
    ]

    par_0$prob_hd$p3 <- out_par[
      (bias_len + hindrance_len + 1):length(out_par)
    ]

    return(par_0)
  }

  maintain_prob_faster <- function(prob_mc,
                                   base_to_pos,
                                   prod_hd,
                                   bias) {
    prob_mp <- matrix(
      data = bias[prob_mc],
      nrow = nrow(prob_mc)
    )

    prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

    return(prob_mp)
  }

  get_cleavage_prob_faster <- function(prob_mc_up,
                                       prob_mc_dn,
                                       base_to_pos,
                                       bias,
                                       prob_hd5,
                                       prob_hd3) {
    # for 5'
    maintain_prob5 <- maintain_prob_faster(
      prob_mc = prob_mc_up,
      base_to_pos = base_to_pos,
      prod_hd = prob_hd5,
      bias = bias
    )

    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

    maintain_cumprod5rev <- matrix(
      data = 1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )

    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      x = maintain_prob5[, length(prob_hd5):2]
    )

    final_p5 <- maintain_cumprod5rev * cle_p5rev

    final_p5 <- final_p5[, length(prob_hd5):1]

    final_p5 <- final_p5 / matrixStats::rowSums2(final_p5)

    # for 3'
    maintain_prob3 <- maintain_prob_faster(
      prob_mc = prob_mc_dn,
      base_to_pos = base_to_pos,
      prod_hd = prob_hd3,
      bias = bias
    )

    cle_p3 <- 1 - maintain_prob3

    maintain_cumprod3 <- matrix(
      data = 1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )

    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      x = maintain_prob3[, -length(prob_hd3)]
    )

    final_p3 <- maintain_cumprod3 * cle_p3

    final_p3 <- final_p3 / matrixStats::rowSums2(final_p3)

    return(
      list(
        final_p5 = final_p5,
        final_p3 = final_p3
      )
    )
  }

  ligation_par <- function(p_num,
                           rpf,
                           candi_psite,
                           candi_cut5,
                           candi_cut3,
                           candi_weight,
                           tx_info,
                           kmer_num,
                           par_0) {
    # 3-mer at the terminal of reads
    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf$pos + rpf$qwidth - 1L,
      width = 3L
    ))

    # choose RPFs, sort expected RPF numbers first
    vec_pos <- sort((c(
      sapply(
        names(par_0$eff_f5), function(kmer_i) {
          which(kmer_5 == kmer_i)[1:kmer_num]
        }
      ),
      sapply(
        names(par_0$eff_f3), function(kmer_i) {
          which(kmer_3 == kmer_i)[1:kmer_num]
        }
      )
    )))

    # prepare to calculate expected RPF numbers
    candi_psite <- candi_psite[vec_pos, ]@x

    candi_weight <- candi_weight[vec_pos, ]

    candi_cut5 <- candi_cut5[vec_pos, ]@x

    candi_cut3 <- candi_cut3[vec_pos, ]@x

    uni_candi_psite <- unique(candi_psite)

    idx_ij <- match(candi_psite, uni_candi_psite)

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    prob_mc_up <- stringr::str_split(
      string = as.character(cut_seq$up_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mc_dn <- stringr::str_split(
      string = as.character(cut_seq$dn_seq),
      pattern = "",
      simplify = TRUE
    )

    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    vec_5 <- nrow(cut_prob$final_p5) * (candi_cut5 - 1L) + idx_ij

    vec_3 <- nrow(cut_prob$final_p3) * (candi_cut3 - 1L) + idx_ij

    obj_num <- rpf$weight[vec_pos]

    vec_weight <- p_num[candi_psite]

    base_to_pos <- 1:4

    bias <- par_0$cut_bias$s7

    hd5_len <- 5:(length(par_0$prob_hd$p5) + 4)

    hd3_len <- seq.int(
      from = length(par_0$prob_hd$p5) + 5,
      to = length(par_0$prob_hd$p5) + length(par_0$prob_hd$p3) + 4
    )

    # calculate expected RPF numbers
    rss_rpf <- function(x) {
      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len],
        prob_hd3 = x[hd3_len]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      return(Matrix::rowSums(candi_weight))
    }

    exp_rpf_num <- rss_rpf(
      x = c(par_0$cut_bias$s7, par_0$prob_hd$p5, par_0$prob_hd$p3)
    )

    # prepare input matrix of model
    mat_train <- data.frame(
      kmer5 = kmer_5[vec_pos],
      kmer3 = kmer_3[vec_pos],
      exp_num = log(exp_rpf_num + 1),
      obj_num = log(obj_num + 1),
      stringsAsFactors = TRUE
    )

    # fit generalized linear model
    glm_fit <- stats::glm(
      formula = obj_num ~ offset(exp_num) + kmer5 + kmer3,
      data = mat_train,
      family = gaussian(link = "identity")
    )

    return(
      list(
        glm_coef = glm_fit,
        mat_train = mat_train,
        kmer_5 = kmer_5,
        kmer_3 = kmer_3
      )
    )
  }

  iter_estimate <- function(par_0,
                            candi_psite,
                            candi_weight,
                            prob_mc_up,
                            prob_mc_dn,
                            vec_5,
                            vec_3,
                            rpf_group_5,
                            rpf_group_3,
                            obj_num,
                            ribo_num,
                            kmer_cof,
                            par0_unknown,
                            init_num,
                            init_tol,
                            sec_tol,
                            sec_diff,
                            iter_times) {
    obj_num <- obj_num + 1

    rss_bias <- function(x) {
      bias[bias_len] <- x[bias_len]

      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len],
        prob_hd3 = x[hd3_len]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      exp_num <- (Matrix::rowSums(candi_weight) + 1) * kmer_cof

      return(
        sum(
          sqrt(
            (exp_num - obj_num)^2
          ) / (exp_num + exp_num^2 / 10)
        )
      )
    }

    rss_hd <- function(x) {
      bias[bias_len] <- par_0$base_prob

      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len - bias_lenght + 1],
        prob_hd3 = x[hd3_len - bias_lenght + 1]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      exp_num <- (Matrix::rowSums(candi_weight) + 1) * kmer_cof

      return(
        sum(
          sqrt(
            (exp_num - obj_num)^2 / (exp_num + exp_num^2 / 10)
          )
        )
      )
    }

    bias_lenght <- length(par_0$cut_bias$s7)
    vec_weight <- ribo_num[candi_psite@x]

    bias_len <- 1:(bias_lenght - 1L)

    base_to_pos <- 1:bias_lenght

    bias <- par_0$cut_bias$s7

    bias[bias_lenght] <- par_0$base_prob

    names(base_to_pos) <- names(par_0$cut_bias$s7)

    hd5_len <- bias_lenght:(length(par_0$prob_hd$p5) + bias_lenght - 1L)

    hd3_len <- seq.int(
      from = length(par_0$prob_hd$p5) + bias_lenght,
      to = length(par_0$prob_hd$p5) +
        length(par_0$prob_hd$p3) +
        bias_lenght - 1L
    )

    if (length(unique(par_0$cut_bias$s7)) == 1) {
      print('no cut bias')
      lst_pars <- list()
      for (round_i in 1:1) {
        if (par0_unknown) {
          lst_par <- lapply(
            1:init_num, function(i) {
              par0 <- c(
                sort(runif(length(hd5_len), 0.01, 0.9), decreasing = TRUE),
                sort(runif(length(hd3_len), 0.01, 0.9))
              )
              res = try(dfoptim::nmkb(
                  par = par0,
                  fn = rss_hd,
                  lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
                  upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
                  control = list(
                    tol = init_tol, # decide search range
                    maxfeval = 5000
                  )
                ),
                silent = TRUE
              )
              if (inherits(res, "try-error")) NULL else res
            }
          )
          lst_par = lst_par[!vapply(lst_par, is.null, logical(1))]
          par_mx = sapply(lst_par, function(opt_out){
            c(opt_out$value, opt_out$par)
          })
          par0_raw <- par_mx[-1, which.min(par_mx[1, ])]
        } else {
          par0_raw <- c(par_0$prob_hd$p5, par_0$prob_hd$p3)
        }

        for (i in seq_along(sec_tol)) {
          diff_value <- 1e4
          value0 <- 0
          round0 <- 1
          par0_raw[which(par0_raw >= (1 - 1e-6))] <- runif(
            length(which(par0_raw >= (1 - 1e-6))),
            (1 - 1e-6 - 2e-7),
            (1 - 1e-6 - 1e-7)
          )
          par0_raw[which(par0_raw <= (1e-8))] <- runif(
            length(which(par0_raw <= (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
          )
          try(
            opt_out <- dfoptim::nmkb(
              par = par0_raw,
              fn = rss_hd,
              lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
              upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
              control = list(
                tol = sec_tol[i], # decide search range
                maxfeval = 5000
              )
            ),
            silent = TRUE
          )

          while ((abs(diff_value) > sec_diff[i]) & (round0 < iter_times)) {
            diff_value <- value0 - opt_out$value
            value0 <- opt_out$value
            par0 <- opt_out$par
            par0[which(par0 > (1 - 1e-6))] <- runif(
              length(which(par0 > (1 - 1e-6))),
              (1 - 1e-6 - 2e-7),
              (1 - 1e-6 - 1e-7)
            )
            par0[which(par0 < (1e-8))] <- runif(
              length(which(par0 < (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
            )
            round0 <- round0 + 1
            lst_opt <- try(
              opt_out <- dfoptim::nmkb(
                par = par0,
                fn = rss_hd,
                lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
                upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
                control = list(
                  tol = sec_tol[i], # decide search range
                  maxfeval = 5000
                )
              ),
              silent = TRUE
            )
            message(!inherits(lst_opt, "try-error"))
            if (!inherits(lst_opt, "try-error")) {
              diff_value <- value0 - opt_out$value
              value0 <- opt_out$value
              par0 <- opt_out$par
            }
          }
        }
        lst_pars[[round_i]] <- par0
      }
      lst_pars <- do.call(rbind, lst_pars)
      par0 <- apply(lst_pars, 2, median)
      par0 <- c(rep(par_0$base_prob, bias_lenght - 1L), par0)
    } else {
      if (par0_unknown) {
        print("cut bias")
        lst_par <- lapply(
          1:init_num, function(i) {
            par0 <- c(
              (par_0$cut_bias$s7 +
                par_0$cut_bias$s7 *
                  runif(length(par_0$cut_bias$s7), -0.5, 0.5))[
                -length(par_0$cut_bias$s7)
              ],
              sort(runif(length(hd5_len), 0.01, 0.9), decreasing = TRUE),
              sort(runif(length(hd3_len), 0.01, 0.9))
            )
            res = try(dfoptim::nmkb(
                par = par0,
                fn = rss_bias,
                lower = c(
                  rep(1e-4, length(bias) - 1),
                  rep(1e-8, length(c(hd3_len, hd5_len)))
                ),
                upper = c(
                  rep(1 - 1e-2, length(bias) - 1),
                  rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
                ),
                control = list(
                  tol = init_tol, # decide search range
                  maxfeval = 5000
                )
              ),
              silent = TRUE
            )
            if (inherits(res, "try-error")) NULL else res
          }
        )
        lst_par = lst_par[!vapply(lst_par, is.null, logical(1))]
        par_mx = sapply(lst_par, function(opt_out){
            c(opt_out$value, opt_out$par)
          })
          par0_raw <- par_mx[-1, which.min(par_mx[1, ])]
      } else {
        par0_raw <- c(
          par_0$cut_bias$s7[-bias_lenght],
          par_0$prob_hd$p5, par_0$prob_hd$p3
        )
      }
      for (i in seq_along(sec_tol)) {
        diff_value <- 1e4
        value0 <- 0
        round0 <- 1
        try(
          opt_out <- dfoptim::nmkb(
            par = par0_raw,
            fn = rss_bias,
            lower = c(
              rep(1e-4, length(bias) - 1),
              rep(1e-8, length(c(hd3_len, hd5_len)))
            ),
            upper = c(
              rep(1 - 1e-2, length(bias) - 1),
              rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
            ),
            control = list(
              tol = sec_tol[i], # decide search range
              maxfeval = 5000
            )
          ),
          silent = TRUE
        )
        while ((abs(diff_value) > sec_diff[i]) & (round0 < iter_times)) {
          diff_value <- value0 - opt_out$value
          value0 <- opt_out$value
          par0 <- opt_out$par

          par0_bias <- par0[1:(length(bias) - 1)]
          par0_hd <- par0[-(1:(length(bias) - 1))]
          par0_bias[which(par0_bias > (1 - 1e-2))] <- runif(
            length(which(par0_bias > (1 - 1e-2))),
            (1 - 1e-2 - 2e-3),
            (1 - 1e-2 - 1e-3)
          )
          par0_bias[1:(length(bias) - 1)][which(par0_bias < (1e-4))] <- runif(
            length(which(par0_bias < (1e-4))), 1e-4 + 1e-5, 1e-4 + 2e-5
          )
          par0_hd[which(par0_hd > (1 - 1e-6))] <- runif(
            length(which(par0_hd > (1 - 1e-6))),
            (1 - 1e-6 - 2e-7),
            (1 - 1e-6 - 1e-7)
          )
          par0_hd[1:(length(bias) - 1)][which(par0_hd < (1e-8))] <- runif(
            length(which(par0_hd < (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
          )
          round0 <- round0 + 1

          lst_opt <- try(
            opt_out <- dfoptim::nmkb(
              par = par0,
              fn = rss_bias,
              lower = c(
                rep(1e-4, length(bias) - 1),
                rep(1e-8, length(c(hd3_len, hd5_len)))
              ),
              upper = c(
                rep(1 - 1e-2, length(bias) - 1),
                rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
              ),
              control = list(
                tol = sec_tol[i], # decide search range
                maxfeval = 5000
              )
            ),
            silent = TRUE
          )
          message(!inherits(lst_opt, "try-error"))
          if (!inherits(lst_opt, "try-error")) {
            diff_value <- value0 - opt_out$value
            value0 <- opt_out$value
            par0 <- opt_out$par
          }
        }
      }
    }

    par0 <- assign_par(
      out_par = par0,
      par_0 = par_0,
      reference_base_value = par_0$base_prob
    )

    return(par0)
  }

  convert_rpf_psite <- function(ribo_size) {
    convert_idx <- matrix(
      data = sapply(
        seq.int(ribo_size[1]), function(cut5_i) {
          sapply(
            seq.int(ribo_size[4]), function(cut3_i) {
              return(c(
                cut5_i,
                cut3_i,
                ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
                ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
              ))
            }
          )
        }
      ),
      nrow = 4
    )
    return(convert_idx)
  }

  filter_prob_matrix <- function(candi_psite,
                                 candi_cut5,
                                 candi_cut3,
                                 candi_weight,
                                 rpf_expected,
                                 ribo_num,
                                 tx_info,
                                 par_0,
                                 num_per_group,
                                 ribo_size) {
    candi_frame <- matrixStats::rowMaxs(as.matrix(candi_cut5 %% 3L))

    len_rg <- seq.int(
      from = (par_0$ribo_size[2] + par_0$ribo_size[3] + 1L),
      to = (sum(par_0$ribo_size) - 1L)
    )

    lst_pos <- lapply(
      len_rg, function(len_i) {
        pos_i <- lapply(
          0:2, function(frame_i) {
            pos_i <- which(
              (rpf_expected$qwidth == len_i) & (candi_frame == frame_i)
            )[1:num_per_group]
            if (is.na(pos_i[1])) {
              return(NULL)
            } else {
              return(pos_i)
            }
          }
        )
        return(pos_i)
      }
    )

    vec_pos <- sort(na.omit(unlist(lst_pos)))

    obj_num <- rpf_expected$weight[vec_pos]

    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf_expected$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf_expected$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf_expected$pos + rpf_expected$qwidth - 1L,
      width = 3L
    ))

    kmer_cof <- exp(
      as.numeric(
        (par_0$eff_intercept +
          par_0$eff_f5[kmer_5] +
          par_0$eff_f3[kmer_3])[vec_pos]
      )
    )

    candi_psite <- candi_psite[vec_pos, ]

    candi_weight <- candi_weight[vec_pos, ]

    candi_cut5 <- candi_cut5[vec_pos, ]

    candi_cut3 <- candi_cut3[vec_pos, ]

    uni_candi_psite <- sort(unique(candi_psite@x))

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    idx_ij <- match(candi_psite@x, uni_candi_psite)

    prob_mc_up <- stringr::str_split(
      string = as.character(cut_seq$up_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mc_dn <- stringr::str_split(
      string = as.character(cut_seq$dn_seq),
      pattern = "",
      simplify = TRUE
    )

    vec_5 <- length(uni_candi_psite) * (candi_cut5@x - 1L) + idx_ij

    vec_3 <- length(uni_candi_psite) * (candi_cut3@x - 1L) + idx_ij

    return(
      list(
        candi_psite = candi_psite,
        candi_weight = candi_weight,
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        vec_5 = vec_5,
        vec_3 = vec_3,
        obj_num = obj_num,
        kmer_cof = kmer_cof
      )
    )
  }

  prepare_prob_matrix <- function(rpf_expected,
                                  convert_idx,
                                  ribo_num,
                                  rg_cds,
                                  mix_tx,
                                  par_0) {
    qwidth <- as.character(rpf_expected$qwidth)

    psite_num_idx <- as.vector(
      sapply(
        split(
          x = convert_idx[1, ], f = convert_idx[4, ]
        ), length
      )[qwidth]
    )

    psite_pos <- rep(x = rpf_expected$pos, times = psite_num_idx) +
      unlist(
        split(
          x = convert_idx[3, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )

    rg_psite <- IRanges::IRanges(start = psite_pos, width = 1L)

    idx_cds <- IRanges::findOverlaps(query = rg_psite, subject = rg_cds)

    idx_inframe <- idx_cds@from[
      (psite_pos[idx_cds@from] - (rg_cds@start + 1L)[idx_cds@to]) %%
        3L == 0L
    ]

    candi_cut3 <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = unlist(
        split(
          x = convert_idx[2, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )[idx_inframe]
    )

    candi_cut5 <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = unlist(
        split(
          x = convert_idx[1, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )[idx_inframe]
    )

    candi_psite <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = psite_pos[idx_inframe]
    )

    candi_p_prob <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = 1.1
    )

    uni_candi_psite <- as.integer(sort(unique(candi_psite@x)))

    cut_seq <- get_cleavage_seq(
      seqs = mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    idx_ij <- match(x = candi_psite@x, table = uni_candi_psite)

    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    candi_p_prob@x <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[
        nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij
      ] *
      ribo_num[candi_psite@x]

    idx_order <- order(Matrix::rowSums(candi_p_prob), decreasing = TRUE)

    candi_psite <- candi_psite[idx_order, ]

    candi_p_prob <- candi_p_prob[idx_order, ]

    candi_cut5 <- candi_cut5[idx_order, ]

    candi_cut3 <- candi_cut3[idx_order, ]

    rpf_expected <- rpf_expected[idx_order, ]

    return(list(
      candi_psite = candi_psite,
      candi_cut5 = candi_cut5,
      candi_cut3 = candi_cut3,
      candi_weight = candi_p_prob,
      rpf_expected = rpf_expected,
      uni_candi_psite = uni_candi_psite
    ))
  }

  p_site2rpf <- function(rpf_integrated,
                         p_site_choosed,
                         ribo_size) {
    convert_idx <- convert_rpf_psite(ribo_size = ribo_size)

    # Calculate the reads distribution of selected ribosomes
    df_read <- data.frame(
      pos = rep(
        x = p_site_choosed,
        each = ncol(convert_idx)
      ) - as.integer(convert_idx[3, ]),
      qwidth = rep(
        x = as.integer(convert_idx[4, ]),
        times = length(p_site_choosed)
      )
    )
    read_num_tag <- as.numeric(
      names(
        base::table(df_read$pos * 100 + df_read$qwidth)
      )
    )

    read_0_tag <- base::setdiff(read_num_tag, rpf_integrated$tag)

    df_exp <- rbind(
      data.frame(
        pos = as.integer(read_0_tag %/% 100),
        qwidth = as.integer(read_0_tag %% 100),
        weight = 0L
      ),
      rpf_integrated[rpf_integrated$tag %in% read_num_tag, -4]
    )

    return(df_exp)
  }

  prepare_par0 <- function(par_5add = c(0, 0, 0, 0, 1),
                           eff_f5 = rep(0, 64),
                           eff_f3 = rep(0, 64),
                           eff_intercept = 0,
                           base_prob = 0.1,
                           cut_bias,
                           mix_tx_len,
                           p_init) {
    prob_add5 <- par_5add
    names(prob_add5) <- c("A", "C", "G", "T")

    eff_intercept <- eff_intercept
    eff_f5 <- eff_f5
    eff_f3 <- eff_f3
    names(eff_f5) <-
      names(eff_f3) <-
      names(Biostrings::GENETIC_CODE)

    cut_bias <- cut_bias
    names(cut_bias) <- c("A", "C", "G", "T")

    par_0 <- list(
      prob_add5 = par_5add,
      ribo_size = p_init$ribo_size,
      prob_hd = list(
        p5 = (seq(810, 10, -800 / p_init$ribo_size[1]) / 1000)[
          1:p_init$ribo_size[1]
        ],
        p3 = (seq(10, 810, 800 / p_init$ribo_size[4]) / 1000)[
          1:p_init$ribo_size[4]
        ]
      ),
      base_prob = base_prob,
      cut_bias = list(
        s7 = cut_bias
      ),
      eff_intercept = eff_intercept,
      eff_f5 = eff_f5,
      eff_f3 = eff_f3
    )

    ribo_num <- vector(mode = "integer", length = mix_tx_len)
    ribo_num[p_init$p_pos] <- p_init$p_num

    return(
      list(
        par_0 = par_0,
        ribo_num = ribo_num
      )
    )
  }

  integrate_rpf <- function(rpf_info) {
    read_tag <- base::table(
      rpf_info$pos * 100 + rpf_info$qwidth
    )

    read_num_tag <- as.numeric(names(read_tag))

    read_weight <- data.frame(
      pos = as.integer(read_num_tag %/% 100),
      qwidth = as.integer(read_num_tag %% 100),
      weight = as.vector(read_tag),
      tag = read_num_tag
    )

    return(read_weight)
  }

  estimate_offsets <- function(tx_info,
                               rpf_info,
                               number_min,
                               max_frame_min = 0.3,
                               choose_terminal = "start",
                               limited_range = 8:18) {
    decide_offset <- function(df_read,
                              rg_annot,
                              choose_terminal,
                              limited_range = limited_range) {
      rg_read <- IRanges::IRanges(
        start = df_read$pos,
        width = df_read$qwidth
      )

      tmp_idx <- IRanges::findOverlaps(
        query = rg_annot,
        subject = rg_read,
        minoverlap = 3L,
        type = "within"
      )

      offsets <- rg_annot@start[tmp_idx@from] - df_read$pos[tmp_idx@to]

      freq_v <- base::table(offsets)[as.character(limited_range)]

      offset_p <- as.integer(as.numeric(names((which.max(freq_v)))))

      rg_len <- df_read$qwidth[1]

      if (choose_terminal == "start") {
        offset_p <- offset_p + 1L
      } else {
        offset_p <- rg_len - offset_p - 2L
      }

      return(c(rg_len, offset_p))
    }
    # Convert coordinate information into IRanges
    rg_reads <- IRanges::IRanges(
      start = rpf_info$pos,
      width = rpf_info$qwidth
    )

    rg_cds <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L,
      end = tx_info$mix_tx_pos$utr3_p5 - 1L
    )

    rg_codon_start <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L,
      width = 3L
    )

    rg_codon_stop <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr3_p5 - 3L,
      width = 3L
    )

    idx_cds <- IRanges::findOverlaps(
      query = rg_reads,
      subject = rg_cds,
      minoverlap = 3L
    )

    rpf_info <- rpf_info[idx_cds@from, ]

    rpf_info$frame_idx <-
      (rg_reads@start[idx_cds@from] - rg_cds@start[idx_cds@to]) %% 3L

    lst_rpf <- base::split(x = rpf_info, f = rpf_info$qwidth)

    # filter groups by number
    gp_num <- sapply(lst_rpf, function(x) {
      nrow(x)
    })

    # choose groups
    max_gp <- which.max(gp_num)
    number_min <- sum(gp_num) * (1 - number_min) * 0.5
    if (number_min < 1000) {
      number_min <- 1000
    }
    tmp1 <- sapply(1:(max_gp - 3), function(i) {
      if ((i == 1 | (i == 2))) {
        FALSE
      } else {
        (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
          (gp_num[i] > number_min) &
          ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
      }
    })
    tmp2 <- sapply((max_gp + 3):length(gp_num), function(i) {
      if ((i == length(gp_num)) | (i == (length(gp_num) - 1))) {
        FALSE
      } else {
        (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
          (gp_num[i] > number_min) &
          ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
      }
    })
    if (sum(tmp1) == 0) {
      rpf_len_5 <- min(which(gp_num > number_min))
    } else {
      rpf_len_5 <- max(which(tmp1))
    }
    if (sum(tmp2) == 0) {
      rpf_len_3 <- max(which(gp_num > number_min))
    } else {
      rpf_len_3 <- min(which(tmp2)) + max_gp
    }

    lst_rpf <- lst_rpf[names(gp_num)[rpf_len_5:rpf_len_3]]

    # filter groups by period
    gp_period <- sapply(
      lst_rpf, function(x) {
        tmp_v <- base::as.vector(base::table(c(0:2, x$frame_idx))) - 1L
        return(tmp_v / sum(tmp_v))
      }
    )

    lst_rpf <- lst_rpf[
      colnames(gp_period)[
        which(apply(gp_period, 2, max) > max_frame_min)
      ]
    ]

    m_offset <- lapply(
      lst_rpf, function(x) {
        if (choose_terminal == "start") {
          offsets <- decide_offset(
            df_read = x,
            rg_annot = rg_codon_start,
            choose_terminal = choose_terminal,
            limited_range = limited_range
          )
        } else {
          offsets <- decide_offset(
            df_read = x,
            rg_annot = rg_codon_stop,
            choose_terminal = choose_terminal,
            limited_range = limited_range
          )
        }
        return(offsets)
      }
    )

    colnames(m_offset) <- NULL

    return(m_offset)
  }

  infer_add5 <- function(df_rpf_add,
                         df_rpf_noadd,
                         tx_info,
                         eff_intercept = 0,
                         eff_f5 = rep(0, 64),
                         eff_f3 = rep(0, 64)) {
    if (unique(eff_f5) == 0) {
      names(eff_f5) <-
        names(eff_f3) <-
        names(Biostrings::GENETIC_CODE)
    }

    # Find the corrected single base added relation matrix
    names_gp <- paste(
      rep(c("A", "C", "G", "T"), each = 4),
      rep(c("A", "C", "G", "T"), 4),
      sep = ""
    )[c(-1, -6, -11, -16)]
    kmer_5 <- Biostrings::xscat(
      Biostrings::DNAStringSet(df_rpf_add$seq),
      Biostrings::subseq(
        Biostrings::DNAStringSet(tx_info$mix_tx)[rep(1L, nrow(df_rpf_add))],
        start = df_rpf_add$pos + 1,
        width = 2L
      )
    )
    kmer_3 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_add))],
      end = df_rpf_add$pos + df_rpf_add$qwidth - 1,
      width = 3
    )
    base_change <- paste(
      as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1L, nrow(df_rpf_add))],
          start = df_rpf_add$pos,
          width = 1L
        )
      ),
      df_rpf_add$seq,
      sep = ""
    )

    correction_num <- exp(
      as.numeric(
        (eff_intercept +
          eff_f5[as.character(kmer_5)] +
          eff_f3[as.character(kmer_3)])
      )
    )

    freq_mat <- matrix(
      sapply(
        names_gp, function(name_gp) {
          sum(correction_num[name_gp == base_change])
        }
      ),
      nrow = 3
    )

    # Ratio of adding different bases, suppose add T probability is 1
    rss <- function(x) {
      freq_mat_exp <- matrix(
        c(x[-1], 1, x[-2], 1, x[-3], 1, x),
        nrow = 3
      )
      obj_m <- t(freq_mat) / Matrix::rowSums(t(freq_mat))
      exp_m <- t(freq_mat_exp) / Matrix::rowSums(t(freq_mat_exp))
      return(
        sum(
          (log(obj_m) - log(exp_m))^2
        )
      )
    }
    par0 <- 3:1
    opt_out <- dfoptim::nmkb(
      par = par0,
      fn = rss,
      lower = rep(.01, 3),
      upper = rep(100, 3)
    )
    add_prob <- c(opt_out$par, 1)

    # Ratio without adding bases
    kmer_5 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_noadd))],
      start = df_rpf_noadd$pos,
      width = 3
    )
    kmer_3 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_noadd))],
      end = df_rpf_noadd$pos + df_rpf_noadd$qwidth - 1,
      width = 3
    )

    correction_num <- exp(
      as.numeric(
        (eff_intercept +
          eff_f5[as.character(kmer_5)] +
          eff_f3[as.character(kmer_3)])
      )
    )

    num_miss_gp <- sum(
      colMeans(
        freq_mat /
          c(add_prob[-1], add_prob[-2], add_prob[-3], add_prob[-4]) *
          rep(add_prob, each = 3)
      )
    )

    add_prob <- c(
      add_prob,
      (sum(correction_num) - num_miss_gp) /
        (sum(freq_mat) + num_miss_gp) *
        sum(add_prob)
    )
    return(add_prob / sum(add_prob))
  }

  prep_rpf <- function(lst_rpf,
                       tx_info,
                       add5,
                       prob_add5,
                       number_min) {
    names(prob_add5) <- c("A", "C", "G", "T", "")
    if (!add5) {
      rpf_info <- data.frame(
        pos = lst_rpf$no_add$pos,
        qwidth = lst_rpf$no_add$qwidth
      )

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = number_min,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 8:18
      )

      tmp_idx <- which(sapply(offsets, length) == 2)
      offsets <- do.call(cbind, offsets[tmp_idx])
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]

      rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
    } else {
      rpf_info_add <- data.frame(
        pos = lst_rpf$add5$pos,
        qwidth = lst_rpf$add5$qwidth
      )

      rpf_info_noadd <- data.frame(
        pos = lst_rpf$no_add$pos,
        qwidth = lst_rpf$no_add$qwidth
      )

      rpf_info <- rbind(rpf_info_noadd, rpf_info_add)

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = number_min,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )
      tmp_idx <- which(sapply(offsets, length) == 2)
      offsets <- do.call(cbind, offsets[tmp_idx])

      rpf_info_add <- rpf_info_add[rpf_info_add$qwidth %in% offsets[1, ], ]
      rpf_info_noadd <- rpf_info_noadd[
        rpf_info_noadd$qwidth %in% offsets[1, ],
      ]

      read_tag_add <- base::table(
        rpf_info_add$pos * 100 + rpf_info_add$qwidth + 99
      )

      read_num_tag_add <- as.numeric(names(read_tag_add))

      read_weight_add <- data.frame(
        pos = as.integer(read_num_tag_add %/% 100),
        qwidth = as.integer(read_num_tag_add %% 100),
        weight = as.vector(read_tag_add),
        tag = read_num_tag_add
      )

      ref_base <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, nrow(read_weight_add))],
          start = read_weight_add$pos - 1,
          width = 1
        )
      )

      tr_base <- sapply(1:4, function(x) {
        1 + prob_add5[x] / sum(prob_add5[1:4][-x])
      })

      read_weight_add$weight <- read_weight_add$weight * tr_base[ref_base]

      read_tag <- base::table(
        rpf_info_noadd$pos * 100 + rpf_info_noadd$qwidth
      )

      read_num_tag <- as.numeric(names(read_tag))

      read_weight <- data.frame(
        pos = as.integer(read_num_tag %/% 100),
        qwidth = as.integer(read_num_tag %% 100),
        weight = as.vector(read_tag),
        tag = read_num_tag
      )

      ref_base1 <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, length(read_num_tag))],
          start = read_num_tag %/% 100,
          width = 1
        )
      )

      tr_base2 <- sapply(1:4, function(x) {
        1 - prob_add5[x] /
          (prob_add5[5] +
            prob_add5[x])
      })

      read_weight$weight <- tr_base2[ref_base1] * read_weight$weight

      rpf_overlap <- rbind(read_weight, read_weight_add)
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]
    }
    return(list(
      rpf_info = rpf_info,
      rpf_overlap = rpf_overlap,
      offsets = offsets
    ))
  }

  initiate_p_site <- function(tx_info,
                              rg_cds,
                              rpf_info,
                              offsets,
                              ribo_size) {
    # initial p-site number
    p_table <- base::table(
      offsets[2, ][match(rpf_info$qwidth, offsets[1, ])] +
        rpf_info$pos
    )

    p_pos <- vector(
      mode = "integer",
      length = length(tx_info$mix_tx)
    )

    p_pos[as.integer(names(p_table))] <- as.vector(p_table)

    idx_p <- sequence(
      nvec = rg_cds@width,
      from = rg_cds@start
    )

    p_cds <- as.integer(matrixStats::colSums2(
      matrix(
        data = p_pos[idx_p],
        nrow = 3
      )
    ))

    idx_p <- matrix(data = idx_p, nrow = 3)[2, ]

    choose_p <- which(p_cds > 0L)

    p_pos <- idx_p[choose_p]

    p_num <- p_cds[choose_p]

    return(
      list(
        rpf_info = rpf_info,
        p_pos = p_pos,
        p_num = p_num,
        ribo_size = ribo_size
      )
    )
  }

  estimate_ribo_size <- function(offsets_mat, rpf_counts, RNase) {
    min_rpf <- as.numeric(names(rpf_counts[1]))
    if (RNase == "rnase-i") {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    } else if (RNase == "mnase") {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    } else {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    }

    S3 <- min_rpf - S5 - 2
    return(as.integer(c(L5, S5, S3, L3)))
  }

  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  trunc_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 16L,
    end = tx_info$mix_tx_pos$utr3_p5 - 16L
  )
  is_add5 <- length(lst_bam) == 2
  if (is_add5) {
    message("Fit 5`add parameters")
    prob_add5 <- infer_add5(
      df_rpf_add = lst_bam$add5,
      df_rpf_noadd = lst_bam$no_add,
      tx_info = tx_info
    )
  } else {
    prob_add5 <- c(0, 0, 0, 0, 1)
  }
  print(prob_add5)
  lst_rpf <- prep_rpf(
    lst_rpf = lst_bam,
    tx_info = tx_info,
    add5 = is_add5,
    prob_add5 = prob_add5,
    number_min = 0.999
  )

  ribo_size <- estimate_ribo_size(
    offsets_mat = lst_rpf$offsets,
    rpf_counts = table(lst_rpf$rpf_info$qwidth),
    RNase = RNase
  )

  init_p_site <- initiate_p_site(
    tx_info = tx_info,
    rg_cds = rg_cds,
    rpf_info = lst_rpf$rpf_info,
    offsets = lst_rpf$offsets,
    ribo_size = ribo_size
  )

  if (!(rnase_bias)) {
    par_0 <- prepare_par0(
      base_prob = 0.1,
      cut_bias = c(0.1, 0.1, 0.1, 0.1),
      mix_tx_len = length(tx_info$mix_tx),
      p_init = init_p_site
    )
  } else {
    par_0 <- prepare_par0(
      base_prob = 0.1,
      cut_bias = c(0.002, 0.7, 0.7, 0.1),
      mix_tx_len = length(tx_info$mix_tx),
      p_init = init_p_site
    )
  }
  convert_idx <- convert_rpf_psite(ribo_size = par_0$par_0$ribo_size)

  choose_p_site <- IRanges::IRanges(
    start = order(par_0$ribo_num, decreasing = TRUE),
    width = 1L
  )

  choosed_p <- choose_p_site@start[
    IRanges::findOverlaps(
      query = choose_p_site,
      subject = trunc_cds
    )@from
  ]

  potential_rpf <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = choosed_p[1:5e2],
    ribo_size = par_0$par_0$ribo_size
  )
  prob_matrix <- prepare_prob_matrix(
    rpf_expected = potential_rpf,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_0$par_0
  )

  prob_matrix_filter <- filter_prob_matrix(
    candi_psite = prob_matrix$candi_psite,
    candi_cut5 = prob_matrix$candi_cut5,
    candi_cut3 = prob_matrix$candi_cut3,
    candi_weight = prob_matrix$candi_weight,
    rpf_expected = prob_matrix$rpf_expected,
    ribo_num = par_0$ribo_num,
    tx_info = tx_info,
    par_0 = par_0$par_0,
    num_per_group = 1e3,
    ribo_size = par_0$par_0$ribo_size
  )

  par_round_1 <- iter_estimate(
    iter_times = 10,
    init_num = 20,
    init_tol = 1e2,
    sec_tol = c(1e2, 1e1, 1e0),
    sec_diff = c(1e0, 1e-1, 1e-2),
    par_0 = par_0$par_0,
    par0_unknown = TRUE,
    candi_psite = prob_matrix_filter$candi_psite,
    candi_weight = prob_matrix_filter$candi_weight,
    prob_mc_up = prob_matrix_filter$prob_mc_up,
    prob_mc_dn = prob_matrix_filter$prob_mc_dn,
    vec_5 = prob_matrix_filter$vec_5,
    vec_3 = prob_matrix_filter$vec_3,
    rpf_group_5 = prob_matrix_filter$rpf_group_5,
    rpf_group_3 = prob_matrix_filter$rpf_group_3,
    obj_num = prob_matrix_filter$obj_num,
    ribo_num = par_0$ribo_num,
    kmer_cof = prob_matrix_filter$kmer_cof
  )
  
  # 大改动 估计——ribonum
  potential_rpf_2 <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = prob_matrix$uni_candi_psite,
    ribo_size = par_round_1$ribo_size
  )

  prob_matrix_2 <- prepare_prob_matrix(
    rpf_expected = potential_rpf_2,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_round_1
  )

  tmp_ribo_num <- em_ribo_num(
    df_rpf = prob_matrix_2$rpf_expected,
    candi_psite = prob_matrix_2$candi_psite,
    candi_cut5 = prob_matrix_2$candi_cut5,
    candi_cut3 = prob_matrix_2$candi_cut3,
    tx_info = tx_info,
    ribo_size = par_round_1$ribo_size,
    par_0 = par_round_1,
    pos_psite = prob_matrix$uni_candi_psite,
    iter_num = 8
  )

  # 第一轮更新ribo_num
  par_0$ribo_num[prob_matrix$uni_candi_psite] <- tmp_ribo_num
  par_0$par_0 <- par_round_1

  # 第二轮估计hd 和 cut_bias
  potential_rpf <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = choosed_p[1:5e2],
    ribo_size = par_0$par_0$ribo_size
  )

  prob_matrix <- prepare_prob_matrix(
    rpf_expected = potential_rpf,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_0$par_0
  )

  prob_matrix_filter <- filter_prob_matrix(
    candi_psite = prob_matrix$candi_psite,
    candi_cut5 = prob_matrix$candi_cut5,
    candi_cut3 = prob_matrix$candi_cut3,
    candi_weight = prob_matrix$candi_weight,
    rpf_expected = prob_matrix$rpf_expected,
    ribo_num = par_0$ribo_num,
    tx_info = tx_info,
    par_0 = par_0$par_0,
    num_per_group = 1e3,
    ribo_size = par_0$par_0$ribo_size
  )

  par_round_2 <- iter_estimate(
    iter_times = 10,
    init_num = 20,
    init_tol = 1e2,
    sec_tol = c(1e2, 1e1, 1e0),
    sec_diff = c(1e0, 1e-1, 1e-2),
    par_0 = par_0$par_0,
    par0_unknown = FALSE,
    candi_psite = prob_matrix_filter$candi_psite,
    candi_weight = prob_matrix_filter$candi_weight,
    prob_mc_up = prob_matrix_filter$prob_mc_up,
    prob_mc_dn = prob_matrix_filter$prob_mc_dn,
    vec_5 = prob_matrix_filter$vec_5,
    vec_3 = prob_matrix_filter$vec_3,
    rpf_group_5 = prob_matrix_filter$rpf_group_5,
    rpf_group_3 = prob_matrix_filter$rpf_group_3,
    obj_num = prob_matrix_filter$obj_num,
    ribo_num = par_0$ribo_num,
    kmer_cof = prob_matrix_filter$kmer_cof
  )
print(par_round_2$ribo_size)
print(par_round_2$cut_bias$s7)
print(par_round_2$prob_hd)
  # 第二轮更新ribo_num
  potential_rpf_2 <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = prob_matrix$uni_candi_psite,
    ribo_size = par_round_2$ribo_size
  )

  prob_matrix_2 <- prepare_prob_matrix(
    rpf_expected = potential_rpf_2,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_round_2
  )

  tmp_ribo_num <- em_ribo_num(
    df_rpf = prob_matrix_2$rpf_expected,
    candi_psite = prob_matrix_2$candi_psite,
    candi_cut5 = prob_matrix_2$candi_cut5,
    candi_cut3 = prob_matrix_2$candi_cut3,
    tx_info = tx_info,
    ribo_size = par_round_2$ribo_size,
    par_0 = par_round_2,
    pos_psite = prob_matrix$uni_candi_psite,
    iter_num = 8
  )
  # 第二轮更新ribo_num后ribo_num最终稳定
  par_0$ribo_num[prob_matrix$uni_candi_psite] <- tmp_ribo_num
  par_0$par_0 <- par_round_2

  if (ligat_par) {
    # 估计ligation参数
    message("Fit ligation parameters")
    potential_rpf_lig_raw <- p_site2rpf(
      rpf_integrated = lst_rpf$rpf_overlap,
      p_site_choosed = choosed_p[1:1e4],
      ribo_size = par_0$par_0$ribo_size
    )

    prob_matrix_lig_raw <- prepare_prob_matrix(
      rpf_expected = potential_rpf_lig_raw,
      convert_idx = convert_idx,
      ribo_num = par_0$ribo_num,
      rg_cds = rg_cds,
      mix_tx = tx_info$mix_tx,
      par_0 = par_0$par_0
    )

    glm_result <- ligation_par(
      p_num = par_0$ribo_num,
      rpf = prob_matrix_lig_raw$rpf_expected,
      candi_psite = prob_matrix_lig_raw$candi_psite,
      candi_cut5 = prob_matrix_lig_raw$candi_cut5,
      candi_cut3 = prob_matrix_lig_raw$candi_cut3,
      candi_weight = prob_matrix_lig_raw$candi_weight,
      tx_info = tx_info,
      kmer_num = 1e4L,
      par_0 = par_0$par_0
    )

    lig_par <- glm2par(
      fit_coefs = glm_result$glm_coef$coefficients,
      par_0 = par_0$par_0
    )

    rpf_overlap_corrected <- lst_rpf$rpf_overlap
  rpf_overlap_corrected$weight <- correct_kmer(
    rpf = lst_rpf$rpf_overlap,
    tx_info = tx_info,
    par_0 = lig_par$par_0
  )

  names(prob_add5) <- c("A", "C", "G", "T", "")
  par_0$par_0$prob_add5 <- prob_add5

  qwidth_range <- unique(convert_rpf_psite(par_0$par_0$ribo_size)[4, ])
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[
    lst_rpf$rpf_overlap$qwidth %in% qwidth_range,
  ]
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[, -4]

  return(list(
    par_0 = par_0$par_0,
    rpf_corrected = rpf_overlap_corrected[, -4],
    lst_rpf = lst_rpf,
    lig_par = lig_par$par_0
  ))
  }else{
    names(prob_add5) <- c("A", "C", "G", "T", "")
  par_0$par_0$prob_add5 <- prob_add5

  qwidth_range <- unique(convert_rpf_psite(par_0$par_0$ribo_size)[4, ])
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[
    lst_rpf$rpf_overlap$qwidth %in% qwidth_range,
  ]
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[, -4]
  return(list(
    par_0 = par_0$par_0,
    rpf_corrected = NULL,
    lst_rpf = lst_rpf,
    lig_par = NULL
  ))
  }
  
}

cut_prob_plot = function(par0){
    base_prob = mean(par0$cut_bias$s7)
    p5 <- par0$prob_hd$p5
    p3 <- par0$prob_hd$p3
    pmf5_raw <- rev(c(1, cumprod(rev(base_prob^p5))[-length(p5)]) * (1 - rev(base_prob^p5)))
    pmf3_raw <-     c(1, cumprod(    base_prob^p3)[-length(p3)]) * (1 -     (base_prob^p3))
    pmf5 <- pmf5_raw / sum(pmf5_raw)
    pmf3 <- pmf3_raw / sum(pmf3_raw)   
  return(cbind(pmf5, pmf3))
}

#' Simulate codon-level P-site occupancies
#'
#' @param candid_orf Candidate ORF structures produced by the preprocessing
#'   helpers.
#' @param par_simu Parameter list returned by [prep_rpf()] or
#'   [prep_rpf_simu()], containing the observed RPF overlap table and offsets.
#' @param tx_info Transcript annotation list with mixed transcript coordinates.
#' @param cds_gp_num Number of quantile groups used to stratify CDS profiles.
#' @param min_sig Minimum number of significant codons required to keep a CDS.
#'
#' @return A list containing simulated positive/negative P-site profiles and the
#'   corresponding genomic ranges for CDS and nORF segments.
#'
#' @keywords internal
simulate_psite <- function(candid_orf,
                           par_simu,
                           tx_info,
                           cds_gp_num = 5,
                           min_sig = 6) {
  # 得到CDS上所有reads数量read_num_orf, CDS有效codon数目sig_num, rpkm_cds, PME
  df_rpf <- par_simu$rpf_info
  offsets <- par_simu$offsets
  idx_qwidth <- df_rpf$qwidth %in% offsets[1, ]
  rg_p <- IRanges::IRanges(
    start = df_rpf$pos[idx_qwidth] +
      offsets[2, ][match(df_rpf$qwidth[idx_qwidth], offsets[1, ])] + 1L,
    width = 1L
  )
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  hit_p <- IRanges::findOverlaps(
    rg_p,
    rg_cds,
    type = "within"
  )
  # split P-site according to CDS
  lst_p_orf <- split(
    (rg_p@start[hit_p@from] - rg_cds@start[hit_p@to]) %/% 3L + 1L,
    hit_p@to
  )
  # 逐组 CDS 长度（与 split 的 names 对齐）
  cds_ids <- as.integer(names(lst_p_orf))
  aa_len_orf <- BiocGenerics::width(rg_cds[cds_ids]) %/% 3L

  # 每组 reads 数
  read_num_orf <- vapply(lst_p_orf, length, integer(1))

  # 单位长度
  unit_len <- as.integer(ceiling(aa_len_orf / pmax(read_num_orf, 1L)))

  rpkm_cds <- read_num_orf / (aa_len_orf * 3L) * 1e3 / sum(read_num_orf) * 1e6
  tx_names <- tx_info$tx_lens$tx_name[cds_ids]
  names(rpkm_cds) <- tx_names

  # PME 计算（逐组对齐）
  pme <- function(p_orf, p_num, len_aa, unit_len) {
    sig_num <- length(unique(p_orf)) # 非零位点数的替代：独立 codon 索引
    freq <- table(p_orf %/% unit_len) / p_num
    vec <- as.numeric(freq)
    pme <- -sum(vec * log2(vec)) / -(log2(1 / ceiling(len_aa / unit_len)) + 1e-10)
    c(sig_num, pme)
  }
  pme_v <- t(mapply(
    function(p, n, L, u) pme(p_orf = p, p_num = n, len_aa = L, unit_len = u),
    lst_p_orf, read_num_orf, aa_len_orf, unit_len,
    SIMPLIFY = TRUE
  ))
  colnames(pme_v) <- c("sig_num", "pme")
  cds_ribo <- as.data.frame(pme_v)
  colnames(cds_ribo) <- c("sig_num", "pme")
  cds_ribo$rpkm_cds <- rpkm_cds
  cds_ribo$reads_num <- read_num_orf
  cds_ribo$aa_length <- aa_len_orf
  cds_ribo$tx_name <- tx_names
  names(lst_p_orf) <- tx_names

  cds_p_num <- lapply(seq_along(aa_len_orf), function(i) {
    as.vector(table(c(1:aa_len_orf[i], lst_p_orf[[i]])) - 1L)
  })
  names(cds_p_num) <- tx_names
  # 函数make_short用CDS的psite分布生成更短的ORF分布
  make_short <- function(sig_num, orf_aa, cds_psite, cds_aa) {
    which_sig <- sample(which(cds_psite > 0), sig_num)
    if (orf_aa <= cds_aa) {
      retain_pos <- c(
        cds_psite[1:3],
        cds_psite[sample(4:(cds_aa - 3), orf_aa - 6)],
        cds_psite[(cds_aa - 2):cds_aa]
      )
    } else {
      retain_pos <- c(
        cds_psite[1:3],
        cds_psite[sample(4:(cds_aa - 3), orf_aa - 6, replace = TRUE)],
        cds_psite[(cds_aa - 2):cds_aa]
      )
    }
    if (sum(retain_pos > 0) < sig_num) {
      retain_pos[sample(1:orf_aa, sig_num)] <- cds_psite[which_sig]
    }
    return(retain_pos)
  }
  # 对CDS分类
  cds_group <- function(data_choose, group_n, min_sig) {
    idx <- data_choose$sig_num >= min_sig
    data_filter <- data_choose[idx, ]
    data_filter <- dplyr::mutate(
      data_filter,
      sig_num_log = log1p(sig_num),
      rpkm_log = log1p(rpkm_cds),
      reads_log = log1p(reads_num),
      aa_len_log = log1p(aa_length)
    )
    df <- dplyr::select(
      data_filter,
      sig_num_log,
      pme,
      rpkm_log,
      reads_log,
      aa_len_log
    )
    df <- dplyr::mutate(
      df,
      dplyr::across(
        dplyr::everything(),
        ~ ifelse(is.finite(.x), .x, NA_real_)
      )
    )
    df <- tidyr::drop_na(df)

    # 2) 标准化
    X <- scale(df)

    # 3) PCA & 方向校正
    pc <- stats::prcomp(X, center = TRUE, scale. = TRUE)
    pc1 <- pc$x[, 1]
    if (stats::cor(pc1, df$rpkm_log) < 0) pc1 <- -pc1 # 让PC1与rpkm正相关

    # 4) 等频5组（1=低，5=高）
    group5 <- dplyr::ntile(pc1, group_n)
    data_choose$class <- 0
    data_choose$class[idx] <- group5
    return(data_choose)
  }
  cds_ribo_filter <- cds_group(
    data_choose = cds_ribo,
    group_n = cds_gp_num,
    min_sig = min_sig
  )
  item_number <- min(c(1000, floor(sum(cds_ribo_filter$class > 0) / cds_gp_num / 2 - 10)))
  # 选择CDS
  cds_choose <- function(data_choose, item_n, lst_psite) {
    lst_df <- split(data_choose$tx_name, data_choose$class)
    lst_df[[1]] <- NULL
    cds_pos_tx <- lapply(lst_df, function(x) {
      matrix(sample(x, 2 * item_n), ncol = 2)
    })
    tx_choose <- do.call(rbind, cds_pos_tx)
    return(list(
      cds_pos = lst_psite[tx_choose[, 1]],
      cds_neg = lst_psite[tx_choose[, 2]]
    ))
  }
  cds_sample <- cds_choose(data_choose = cds_ribo_filter, item_n = item_number, lst_psite = cds_p_num)
  pos_cds_name <- names(cds_sample$cds_pos)
  neg_cds_name <- names(cds_sample$cds_neg)
  rg_tx <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p5,
    end = tx_info$mix_tx_pos$utr3_p3
  )
  rg_tx_pos <- rg_tx[match(pos_cds_name, tx_info$tx_lens$tx_name)]
  norf_choose <- sapply(setdiff(names(candid_orf$orf_relation), "n_extend"), function(class_x) {
    rg_class_x <- candid_orf$orf_relation[[class_x]]
    hit_orf <- IRanges::findOverlaps(rg_class_x, rg_tx_pos, type = "any")
    return(rg_class_x[hit_orf@from])
  }, simplify = FALSE, USE.NAMES = TRUE)
  norf_choose$down_orf <- sample(norf_choose$down_orf, length(norf_choose$up_orf))
  norf_choose$outframe_orf <- sample(norf_choose$outframe_orf, length(norf_choose$up_orf))
  norf_filter <- norf_choose
  names(norf_filter) <- NULL
  norf_filter <- do.call(c, norf_filter)
  # 选择nORF在挑选为阳性的CDS上
  norf_psite <- function(rg_data, min_sig, cds_p_lst, rg_tx_pos, cds_info) {
    idx <- sample(c("y", "n"), length(rg_data), replace = T)
    norf_aa <- (IRanges::width(rg_data) + 2) / 3
    cds_aa <- sapply(cds_p_lst, length)
    cds_i <- IRanges::findOverlaps(rg_data, rg_tx_pos, type = "any")@to
    idx_cds <- match(names(cds_p_lst)[cds_i], cds_info$tx_name)
    prop_i <- ceiling(cds_info$sig_num[idx_cds] / cds_info$aa_length[idx_cds] * norf_aa)
    tp_idx <- prop_i >= cds_info$sig_num[idx_cds]
    prop_i[tp_idx] <- (cds_info$sig_num[idx_cds] / 2)[tp_idx]
    prop_i[prop_i < min_sig] <- min_sig
    orf_psite <- sapply(seq_along(norf_aa), function(i) {
      make_short(
        sig_num = prop_i[i],
        orf_aa = norf_aa[i],
        cds_psite = cds_p_lst[cds_i][[i]],
        cds_aa = cds_info$aa_length[idx_cds][i]
      )
    })
    return(list(
      psite_norf = orf_psite,
      norf_idx = idx
    ))
  }
  norf_psite_lst <- norf_psite(
    rg_data = norf_filter,
    min_sig = min_sig,
    cds_p_lst = cds_sample$cds_pos,
    rg_tx_pos = rg_tx_pos,
    cds_info = cds_ribo_filter
  )
  rg_cds_pos <- rg_cds[match(pos_cds_name, tx_info$tx_lens$tx_name)]
  names(rg_cds_pos) <- pos_cds_name
  rg_cds_neg <- rg_cds[match(neg_cds_name, tx_info$tx_lens$tx_name)]
  names(rg_cds_neg) <- neg_cds_name
  psite_pos <- vector(mode = "integer", length = length(tx_info$mix_tx))
  idx_p1 <- sequence(
    nvec = rg_cds_pos@width %/% 3L,
    from = rg_cds_pos@start + 1L,
    by = 3L
  )
  tp_idx <- norf_psite_lst$norf_idx == "y"
  idx_p2 <- sequence(
    nvec = (norf_filter@width[tp_idx] + 2L) %/% 3L,
    from = norf_filter@start[tp_idx] + 1L,
    by = 3L
  )
  psite_pos[idx_p1] <- unlist(cds_sample$cds_pos, recursive = FALSE, use.names = FALSE)
  psite_pos[idx_p2] <- unlist(norf_psite_lst$psite_norf[tp_idx], recursive = FALSE, use.names = FALSE)

  psite_neg <- vector(mode = "integer", length = length(tx_info$mix_tx))
  idx_p1 <- sequence(
    nvec = rg_cds_neg@width %/% 3L,
    from = rg_cds_neg@start + 1L,
    by = 3L
  )
  idx_p2 <- sequence(
    nvec = (norf_filter@width[!tp_idx] + 2L) %/% 3L,
    from = norf_filter@start[!tp_idx] + 1L,
    by = 3L
  )
  # 用多项式分配对每个psite上核糖体数量拆散成3位平均分配
  psite_bi <- function(psite) {
    ps <- psite
    idx <- which(ps > 0L)
    v <- ps[idx]
    left <- stats::rbinom(length(v), size = v, prob = 1 / 3)
    mid <- stats::rbinom(length(v), size = v - left, prob = 1 / 2)
    right <- v - left - mid
    ans <- integer(length(ps))
    ans[idx - 1L] <- left
    ans[idx] <- ans[idx] + mid
    ans[idx + 1L] <- ans[idx + 1L] + right
    return(ans)
  }
  psite_neg[idx_p1] <- unlist(cds_sample$cds_neg, recursive = FALSE, use.names = FALSE)
  psite_neg[idx_p2] <- unlist(norf_psite_lst$psite_norf[!tp_idx], recursive = FALSE, use.names = FALSE)
  psite_change <- psite_bi(psite = psite_neg)
  simulate_info <- list(
    psite_pos = psite_pos,
    psite_neg = psite_change,
    rg_cds_pos = rg_cds_pos,
    rg_cds_neg = rg_cds_neg,
    rg_norf_pos = norf_filter[tp_idx],
    rg_norf_neg = norf_filter[!tp_idx]
  )
  return(simulate_info)
}

#' Generate simulated ribosome footprint reads
#'
#' Internal helper that reproduces the simulation logic used during parameter
#' learning. The implementation is adapted from the training utilities embedded
#' in the original workflow.
#'
#' @param par_setup Output of [infer_par()]$par_0 containing ligation and
#'   cleavage parameters.
#' @param tx_info Transcript annotation list with mixed transcript coordinates.
#' @param p_vec Integer vector indicating the number of ribosomes at each
#'   codon position.
#' @param add_noise Logical; if `TRUE`, Gaussian noise is added to the sampled
#'   reads.
#' @param noise_ratio Numeric noise level applied when `add_noise` is `TRUE`.
#'
#' @return A list with the simulated reads and P-site correction metadata.
#'
#' @keywords internal
simulate_reads_internal <- function(par_setup,
                                    tx_info,
                                    p_vec,
                                    add_noise,
                                    noise_ratio) {
  # Separately generate sequences preceding cleavage at both ribosome termini
  fun_cut_seq <- function(seqs,
                          p_pos,
                          ribo_size) {
    p_site <- which(p_pos > 0L)
    seqs <- Biostrings::DNAStringSet(seqs)[rep(1L, length(p_site))]
    up_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site - sum(ribo_size[1:2]),
      width = ribo_size[1]
    )
    dn_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site + ribo_size[3],
      width = ribo_size[4]
    )
    return(
      list(
        up_seq = up_seq,
        dn_seq = dn_seq
      )
    )
  }
  # Generating cleavage probabilities for each ribosome terminus
  fun_cut_prob <- function(seqs,
                           bias,
                           prob_hd5,
                           prob_hd3) {
    maintain_prob <- function(cut_seq,
                              prod_hd) {
      prob_mc <- stringr::str_split(
        as.vector(cut_seq),
        pattern = "",
        simplify = TRUE
      )
      prob_mp <- matrix(
        data = bias[prob_mc],
        nrow = nrow(prob_mc)
      )
      prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)
      return(prob_mp)
    }
    # for 5'
    maintain_prob5 <- maintain_prob(
      cut_seq = seqs$up_seq,
      prod_hd = prob_hd5
    )
    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]
    maintain_cumprod5rev <- matrix(
      1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )
    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      maintain_prob5[, length(prob_hd5):2]
    )
    final_p5 <- maintain_cumprod5rev * cle_p5rev
    final_p5 <- cbind(final_p5, 1 - Matrix::rowSums(final_p5))
    final_p5 <- final_p5[, ncol(final_p5):1]

    # for 3'
    maintain_prob3 <- maintain_prob(
      cut_seq = seqs$dn_seq,
      prod_hd = prob_hd3
    )
    cle_p3 <- 1 - maintain_prob3
    maintain_cumprod3 <- matrix(
      1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )
    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      maintain_prob3[, -length(prob_hd3)]
    )
    final_p3 <- maintain_cumprod3 * cle_p3
    final_p3 <- cbind(final_p3, 1 - Matrix::rowSums(final_p3))

    return(list(final_p5 = final_p5, final_p3 = final_p3))
  }
  # Integrate various parameters to generate simulated reads
  fun_create_read <- function(seqs,
                              cut_exp,
                              p_pos,
                              ribo_size) {
    p_site <- which(p_pos > 0L)
    p_pos_expand <- rep(p_site, p_pos[p_site])
    expand_m5 <- matrix(
      rep(
        cut_exp$final_p5,
        rep(
          p_pos[p_site],
          ncol(cut_exp$final_p5)
        )
      ),
      ncol = ncol(cut_exp$final_p5)
    )
    expand_m3 <- matrix(
      rep(
        cut_exp$final_p3,
        rep(
          p_pos[p_site], ncol(cut_exp$final_p3)
        )
      ),
      ncol = ncol(cut_exp$final_p3)
    )
    pos5idx <- as.integer(
      rowSums(
        runif(nrow(expand_m5)) > matrixStats::rowCumsums(expand_m5)
      )
    )
    pos3idx <- as.integer(
      rowSums(
        runif(nrow(expand_m3)) > matrixStats::rowCumsums(expand_m3)
      )
    )
    out_reg <- (pos5idx != 0L) & (pos3idx != ribo_size[4])
    cut_pos5 <- pos5idx + p_pos_expand - sum(ribo_size[1:2]) - 1L
    cut_pos3 <- pos3idx + p_pos_expand + ribo_size[3]

    p_correct <- as.vector(
      base::table(
        c(p_pos_expand[out_reg], unique(p_pos_expand))
      ) - 1L
    )
    p_original <- as.vector(base::table(p_pos_expand))

    reads <- Biostrings::subseq(
      Biostrings::DNAStringSet(seqs)[rep(1L, length(cut_pos5))],
      start = cut_pos5,
      end = cut_pos3 - 1L
    )

    return(list(
      reads = reads[out_reg],
      p_correct = p_correct,
      p_pos = p_site,
      p_original = p_original
    ))
  }

  cut_seq <- fun_cut_seq(
    seqs = tx_info$mix_tx,
    p_pos = p_vec,
    ribo_size = par_setup$ribo_size
  )
  cut_prob <- fun_cut_prob(
    seqs = cut_seq,
    bias = par_setup$cut_bias$s7,
    prob_hd5 = par_setup$prob_hd$p5,
    prob_hd3 = par_setup$prob_hd$p3
  )
  create_read <- fun_create_read(
    seqs = tx_info$mix_tx,
    cut_exp = cut_prob,
    p_pos = p_vec,
    ribo_size = par_setup$ribo_size
  )
  add_base <- sample(
    names(par_setup$prob_add5),
    length(create_read$reads),
    replace = TRUE,
    prob = par_setup$prob_add5
  )
  create_read$add_base <- Biostrings::DNAStringSet(add_base)
  create_read$reads_add5 <- Biostrings::xscat(
    Biostrings::DNAStringSet(add_base),
    create_read$reads
  )
  create_read$reads_add5_pos <- create_read$reads@ranges@start
  kmer5 <- Biostrings::subseq(
    create_read$reads_add5,
    start = 1,
    width = 3
  )
  kmer3 <- Biostrings::subseq(
    create_read$reads_add5,
    end = BiocGenerics::width(create_read$reads_add5),
    width = 3
  )
  ligation_prob <- exp(
    par_setup$eff_intercept +
      par_setup$eff_f5[as.vector(kmer5)] +
      par_setup$eff_f3[as.vector(kmer3)]
  )
  idx_ligat <- runif(
    length(create_read$reads_add5)
  ) < ligation_prob
  create_read$idx_ligat <- idx_ligat
  create_read$add_base <- create_read$add_base[idx_ligat]
  create_read$reads_add5 <- create_read$reads_add5[idx_ligat]
  create_read$reads_add5_pos <- create_read$reads_add5_pos[idx_ligat]
  if (isTRUE(add_noise)) {
    n_pos <- create_read$reads_add5_pos +
      sample(
        -30:30,
        length(create_read$reads_add5_pos),
        replace = TRUE
      )
    n_width <- sample(
      25:34,
      length(create_read$reads_add5_pos),
      replace = TRUE,
      prob = c(1:5, 5:1)
    )
    create_read$reads_add5_pos <- c(
      create_read$reads_add5_pos,
      sample(
        n_pos,
        size = noise_ratio * length(n_pos)
      )
    )
    create_read$reads <- c(
      create_read$reads,
      Biostrings::subseq(
        Biostrings::DNAStringSet(tx_info$mix_tx)[
          rep(
            1L,
            length(create_read$reads_add5_pos)
          )
        ],
        start = n_pos,
        width = n_width
      )[
        sample(
          seq_along(n_pos),
          noise_ratio * length(n_pos)
        )
      ]
    )
  }

  return(create_read)
}

#' Simulate ribosome profiling reads for benchmarking
#'
#' Generates simulated ribosome profiling reads by combining
#' [simulate_psite()] with the internal [simulate_reads()] generator.
#'
#' @param candid_orf Candidate ORFs produced by [find_candid_orf()].
#' @param par_rpf Parameter object returned by [infer_par()], containing the
#'   learned ribosome footprint parameters.
#' @param tx_info Transcript annotation list with mixed transcript coordinates.
#'
#' @return A list with two components: `simulated_reads`, the output of
#'   [simulate_reads()], and `orf_info`, the P-site simulation metadata.
#'
#' @export
simulation_reads <- function(candid_orf,
                             par_rpf,
                             tx_info) {
  psite_simu <- simulate_psite(
    candid_orf = candid_orf,
    par_simu = par_rpf$lst_rpf,
    tx_info = tx_info
  )
  # for positive reads
  psite_input_pos <- stats::rmultinom(
    n = 1,
    size = sum(psite_simu$psite_pos),
    prob = psite_simu$psite_pos
  )[, 1]
  # for negative reads
  psite_input_neg <- stats::rmultinom(
    n = 1,
    size = sum(psite_simu$psite_neg),
    prob = psite_simu$psite_neg
  )[, 1]

  psite_input <- psite_input_pos + psite_input_neg

  par_0 <- par_rpf$par_0
  simulated_reads <- simulate_reads_internal(
    par_setup = par_0,
    tx_info = tx_info,
    p_vec = psite_input,
    add_noise = FALSE,
    noise_ratio = 0.5
  )

  return(list(
    simulated_reads = simulated_reads,
    orf_info = psite_simu
  ))
}

estimate_offsets <- function(tx_info,
                             rpf_info,
                             number_min,
                             max_frame_min = 0.3,
                             choose_terminal = "start",
                             limited_range = 8:15) {
  decide_offset <- function(df_read,
                            rg_annot,
                            choose_terminal,
                            limited_range = limited_range) {
    rg_read <- IRanges::IRanges(
      start = df_read$pos,
      width = df_read$qwidth
    )

    tmp_idx <- IRanges::findOverlaps(
      query = rg_annot,
      subject = rg_read,
      minoverlap = 3L,
      type = "within"
    )

    offsets <- rg_annot@start[tmp_idx@from] - df_read$pos[tmp_idx@to]

    freq_v <- base::table(offsets)[as.character(limited_range)]

    offset_p <- as.integer(as.numeric(names((which.max(freq_v)))))

    rg_len <- df_read$qwidth[1]

    if (choose_terminal == "start") {
      offset_p <- offset_p + 1L
    } else {
      offset_p <- rg_len - offset_p - 2L
    }

    return(c(rg_len, offset_p))
  }
  # Convert coordinate information into IRanges
  rg_reads <- IRanges::IRanges(
    start = rpf_info$pos,
    width = rpf_info$qwidth
  )

  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )

  rg_codon_start <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    width = 3L
  )

  rg_codon_stop <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr3_p5 - 3L,
    width = 3L
  )

  idx_cds <- IRanges::findOverlaps(
    query = rg_reads,
    subject = rg_cds,
    minoverlap = 3L
  )

  rpf_info <- rpf_info[idx_cds@from, ]

  rpf_info$frame_idx <-
    (rg_reads@start[idx_cds@from] - rg_cds@start[idx_cds@to]) %% 3L

  lst_rpf <- base::split(x = rpf_info, f = rpf_info$qwidth)

  # filter groups by number
  gp_num <- sapply(lst_rpf, function(x) {
    nrow(x)
  })

  # choose groups
  max_gp <- which.max(gp_num)
  number_min <- sum(gp_num) * (1 - number_min) * 0.5
  # tmp1 <- sapply(1:(max_gp - 3), function(i) {
  #   if ((i == 1 | (i == 2))) {
  #     FALSE
  #   } else {
  #     (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
  #       (gp_num[i] > number_min) &
  #       ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
  #   }
  # })
  # tmp2 <- sapply((max_gp + 3):length(gp_num), function(i) {
  #   if ((i == length(gp_num)) | (i == (length(gp_num) - 1))) {
  #     FALSE
  #   } else {
  #     (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
  #       (gp_num[i] > number_min) &
  #       ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
  #   }
  # })
  # if (sum(tmp1) == 0) {
  #   rpf_len_5 <- min(which(gp_num > number_min))
  # } else {
  #   rpf_len_5 <- max(which(tmp1))
  # }
  # if (sum(tmp2) == 0) {
  #   rpf_len_3 <- max(which(gp_num > number_min))
  # } else {
  #   rpf_len_3 <- min(which(tmp2)) + max_gp
  # }
  #
  # lst_rpf <- lst_rpf[names(gp_num)[rpf_len_5:rpf_len_3]]
  lst_rpf <- lst_rpf[names(gp_num)[
    min(which(gp_num > number_min)):max(which(gp_num > number_min))
  ]]
  # filter groups by period
  gp_period <- sapply(
    lst_rpf, function(x) {
      tmp_v <- base::as.vector(base::table(c(0:2, x$frame_idx))) - 1L
      return(tmp_v / sum(tmp_v))
    }
  )

  lst_rpf <- lst_rpf[
    colnames(gp_period)[
      which(apply(gp_period, 2, max) > max_frame_min)
    ]
  ]

  m_offset <- sapply(
    lst_rpf, function(x) {
      if (choose_terminal == "start") {
        offsets <- decide_offset(
          df_read = x,
          rg_annot = rg_codon_start,
          choose_terminal = choose_terminal,
          limited_range = limited_range
        )
      } else {
        offsets <- decide_offset(
          df_read = x,
          rg_annot = rg_codon_stop,
          choose_terminal = choose_terminal,
          limited_range = limited_range
        )
      }
      return(offsets)
    }
  )

  colnames(m_offset) <- NULL

  return(m_offset)
}

integrate_rpf <- function(rpf_info) {
  read_tag <- base::table(
    rpf_info$pos * 100 + rpf_info$qwidth
  )

  read_num_tag <- as.numeric(names(read_tag))

  read_weight <- data.frame(
    pos = as.integer(read_num_tag %/% 100),
    qwidth = as.integer(read_num_tag %% 100),
    weight = as.vector(read_tag),
    tag = read_num_tag
  )

  return(read_weight)
}

pred_p <- function(par_rpf,
                   tx_info,
                   iter_num) {
  # iterate n times to reach balance
  pred_psite <- function(df_rpf,
                         candi_psite,
                         candi_cut5,
                         candi_cut3,
                         tx_info,
                         ribo_size,
                         par_0,
                         iter_num) {
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
    )
    idx_ij <- match(candi_psite@x, uni_candi_psite)
    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq, bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
    )
    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
    prob_cum <- vector(mode = "numeric", length = iter_num)

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
      tmp_sum <- Matrix::rowSums(candi_p_weight)
      prob_cum[i] <- sum(df_rpf$weight * log(tmp_sum))
      sparse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / tmp_sum)@x
      )
      vec_pnum[] <- Matrix::rowSums(sparse_iter)
    }
    return(
      list(
        vec_pnum = vec_pnum,
        candi_pos = candi_pos,
        prob_cum = prob_cum
      )
    )
  }

  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[
    par_rpf$lst_rpf$rpf_overlap$qwidth %in% unique(convert_idx[4, ]),
  ]
  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = psite_pos
  )
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  ribo_num_0 <- pred_psite(
    df_rpf = par_rpf$lst_rpf$rpf_overlap,
    candi_psite = candi_psite,
    candi_cut5 = candi_cut5,
    candi_cut3 = candi_cut3,
    tx_info = tx_info,
    ribo_size = ribo_size,
    par_0 = par_rpf$par_0,
    iter_num = iter_num
  )
  return(ribo_num_0)
}


get_cleavage_seq <- function(seqs,
                             p_site,
                             ribo_size) {
  seqs <- Biostrings::DNAStringSet(x = seqs)[
    rep(1L, length(p_site))
  ]

  up_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site - sum(ribo_size[1:2]),
    width = ribo_size[1]
  )

  dn_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site + ribo_size[3],
    width = ribo_size[4]
  )

  return(
    list(
      up_seq = up_seq,
      dn_seq = dn_seq
    )
  )
}

maintain_prob <- function(cut_seq,
                          prod_hd,
                          bias) {
  prob_mc <- stringr::str_split(
    string = as.vector(cut_seq),
    pattern = "",
    simplify = TRUE
  )

  prob_mp <- matrix(
    data = bias[prob_mc],
    nrow = nrow(prob_mc)
  )

  prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

  return(prob_mp)
}

get_cleavage_prob <- function(seqs,
                              bias,
                              prob_hd5,
                              prob_hd3) {
  # for 5'
  maintain_prob5 <- maintain_prob(
    cut_seq = seqs$up_seq,
    prod_hd = prob_hd5,
    bias = bias
  )

  cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

  maintain_cumprod5rev <- matrix(
    data = 1,
    nrow = nrow(maintain_prob5),
    ncol = length(prob_hd5)
  )

  maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
    maintain_prob5[, length(prob_hd5):2]
  )

  final_p5 <- maintain_cumprod5rev * cle_p5rev

  final_p5 <- final_p5[, ncol(final_p5):1]

  final_p5 <- final_p5 / Matrix::rowSums(final_p5)

  # for 3'
  maintain_prob3 <- maintain_prob(
    cut_seq = seqs$dn_seq,
    prod_hd = prob_hd3,
    bias = bias
  )

  cle_p3 <- 1 - maintain_prob3

  maintain_cumprod3 <- matrix(
    data = 1,
    nrow = nrow(maintain_prob3),
    ncol = length(prob_hd3)
  )

  maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
    maintain_prob3[, -length(prob_hd3)]
  )

  final_p3 <- maintain_cumprod3 * cle_p3

  final_p3 <- final_p3 / Matrix::rowSums(final_p3)

  return(
    list(
      final_p5 = final_p5,
      final_p3 = final_p3
    )
  )
}

prep_rpf_simu <- function(lst_simu,
                          tx_info,
                          add5,
                          prob_add5) {
  if (!add5) {
    rpf_info <- data.frame(
      pos = lst_simu$reads_add5_pos,
      qwidth = BiocGenerics::width(lst_simu$reads_add5)
    )

    offsets <- estimate_offsets(
      tx_info = tx_info,
      rpf_info = rpf_info,
      number_min = 0.9,
      max_frame_min = 0.3,
      choose_terminal = "start",
      limited_range = 10:18
    )

    rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]

    rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
  } else {
    if (prob_add5[5] == 1) {
      rpf_info <- data.frame(
        pos = lst_simu$reads_add5_pos - BiocGenerics::width(lst_simu$add_base),
        qwidth = BiocGenerics::width(lst_simu$reads_add5)
      )

      ref_seqs <- Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, nrow(rpf_info))],
        start = rpf_info$pos,
        width = rpf_info$qwidth
      )

      idx_add <- lst_simu$reads_add5 == ref_seqs

      rpf_info_noadd <- rpf_info[idx_add, ]
      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info_noadd,
        number_min = 0.99,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )

      rpf_info <- rpf_info_noadd[rpf_info_noadd$qwidth %in% offsets[1, ], ]

      rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
    } else {
      rpf_info <- data.frame(
        pos = lst_simu$reads_add5_pos - BiocGenerics::width(lst_simu$add_base),
        qwidth = BiocGenerics::width(lst_simu$reads_add5)
      )

      ref_seqs <- Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, nrow(rpf_info))],
        start = rpf_info$pos,
        width = rpf_info$qwidth
      )

      idx_add <- lst_simu$reads_add5 == ref_seqs

      rpf_info_add <- rpf_info[!idx_add, ]
      rpf_info_noadd <- rpf_info[idx_add, ]

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = 0.99,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )

      rpf_info_add <- rpf_info_add[rpf_info_add$qwidth %in% offsets[1, ], ]
      rpf_info_noadd <- rpf_info_noadd[
        rpf_info_noadd$qwidth %in% offsets[1, ],
      ]

      read_tag_add <- base::table(
        rpf_info_add$pos * 100 + rpf_info_add$qwidth + 99
      )

      read_num_tag_add <- as.numeric(names(read_tag_add))

      read_weight_add <- data.frame(
        pos = as.integer(read_num_tag_add %/% 100),
        qwidth = as.integer(read_num_tag_add %% 100),
        weight = as.vector(read_tag_add),
        tag = read_num_tag_add
      )

      ref_base <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, nrow(read_weight_add))],
          start = read_weight_add$pos - 1,
          width = 1
        )
      )

      tr_base <- sapply(1:4, function(x) {
        1 + prob_add5[x] / sum(prob_add5[1:4][-x])
      })

      read_weight_add$weight <- read_weight_add$weight * tr_base[ref_base]

      read_tag <- base::table(
        rpf_info_noadd$pos * 100 + rpf_info_noadd$qwidth
      )

      read_num_tag <- as.numeric(names(read_tag))

      read_weight <- data.frame(
        pos = as.integer(read_num_tag %/% 100),
        qwidth = as.integer(read_num_tag %% 100),
        weight = as.vector(read_tag),
        tag = read_num_tag
      )

      ref_base1 <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, length(read_num_tag))],
          start = read_num_tag %/% 100,
          width = 1
        )
      )

      tr_base2 <- sapply(1:4, function(x) {
        1 - prob_add5[x] /
          (prob_add5[5] +
            prob_add5[x])
      })

      read_weight$weight <- tr_base2[ref_base1] * read_weight$weight

      rpf_overlap <- rbind(read_weight, read_weight_add)
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]
    }
  }
  return(list(
    rpf_info = rpf_info,
    rpf_overlap = rpf_overlap,
    offsets = offsets
  ))
}

# 模拟数据训练模型
train_translation_classifier <- function(tx_info,
                                         ribo_seq_simulate,
                                         candidate_translate_region,
                                         simulate_orf_info,
                                         min_reads_number,
                                         active_codon_number,
                                         n_cds_states,
                                         threshold_method = c("pr_f1", "youden"),
                                         par_rpf_real,
                                         train_iter_num) {
  # Feature extractor for translated candidate regions
  extract_translation_features <- function(ribo_num,
                                           candi_translat,
                                           tx_info,
                                           min_pos_num,
                                           eff_ribo_num,
                                           bin_num) {
    # Project ribosome footprints onto candidate translated segments
    rg_p <- IRanges::IRanges(
      start = ribo_num[[1]]$candi_pos,
      width = 1L
    )
    candi_orf <- do.call(c, candi_translat)
    candi_orf@start <- candi_orf@start + 3L
    candi_orf@width <- candi_orf@width - 1L

    hit_p <- IRanges::findOverlaps(rg_p, candi_orf)
    hit_p <- hit_p[
      hit_p@to %in%
        which(
          IRanges::countOverlaps(candi_orf, rg_p) >= min_pos_num
        )
    ]
    dist_p <- rg_p@start[hit_p@from] - candi_orf@start[hit_p@to]
    max_dist <- max(dist_p) %/% 3L + 3L
    idx_frame <- dist_p %% 3L
    idx_pos <- (dist_p %/% 3L) + 1L
    orf_n <- as.factor(hit_p@to)
    idx_orf <- as.integer(orf_n)
    feat_lst <- lapply(seq_along(ribo_num), function(ribo_num_i) {
      p_num <- ribo_num[[ribo_num_i]]$vec_pnum

      # Build sparse codon-by-codon count matrices, one ORF per row
      spa_frame0 <- Matrix::sparseMatrix(
        i = c(idx_orf[idx_frame == 0L], length(levels(orf_n)) + 1L),
        j = c(idx_pos[idx_frame == 0L], max_dist),
        x = c(p_num[hit_p@from][idx_frame == 0L], 1)
      )
      spa_frame0 <- spa_frame0[-nrow(spa_frame0), ]
      spa_frame1 <- Matrix::sparseMatrix(
        i = c(idx_orf[idx_frame == 1L], length(levels(orf_n)) + 1L),
        j = c(idx_pos[idx_frame == 1L], max_dist),
        x = c(p_num[hit_p@from][idx_frame == 1L], 1)
      )
      spa_frame1 <- spa_frame1[-nrow(spa_frame1), ]
      spa_frame2 <- Matrix::sparseMatrix(
        i = c(idx_orf[idx_frame == 2L], length(levels(orf_n)) + 1L),
        j = c(idx_pos[idx_frame == 2L], max_dist),
        x = c(p_num[hit_p@from][idx_frame == 2L], 1)
      )
      spa_frame2 <- spa_frame2[-nrow(spa_frame2), ]

      # Normalize ribosome counts across frames for each codon
      orf_id <- as.integer(
        levels(as.factor(hit_p@to))
      )
      spa0 <- spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
      spa0@x <- rep(0, length(spa0@x))
      spa_frame0 <- spa_frame0 + spa0
      spa_frame1 <- spa_frame1 + spa0
      spa_frame2 <- spa_frame2 + spa0
      spa_inframe@x <- spa_frame0@x / spa_all@x
      spa_frame0_fq <- spa_inframe
      spa_inframe@x <- spa_frame1@x / spa_all@x
      spa_frame1_fq <- spa_inframe
      spa_inframe@x <- spa_frame2@x / spa_all@x
      spa_frame2_fq <- spa_inframe

      orf_aa <- Biostrings::width(candi_orf[orf_id]) / 3

      orf_bin <- Matrix::summary(spa_all)
      orf_bin$x <- floor((orf_bin$j - 0.5) / orf_aa[orf_bin$i] * bin_num) + 1
      spa_bin <- Matrix::sparseMatrix(
        i = orf_bin$i,
        j = orf_bin$j,
        x = orf_bin$x,
        dims = dim(spa_all)
      )
      bin_lst <- lapply(seq.int(bin_num), function(i) {
        spa_bin == i
      })
      spa_eff <- spa_all > eff_ribo_num
      # Summaries within positional bins
      bin_eff <- function(bin_lst, spa_eff) {
        sapply(bin_lst, function(bin_i) {
          Matrix::rowSums(spa_eff & bin_i)
        })
      }

      bin_prop_sum <- function(spa_prop, bin_lst, spa_eff) {
        sapply(bin_lst, function(bin_i) {
          Matrix::rowSums(spa_prop * (spa_eff & bin_i))
        })
      }
      spa_all_log <- spa_all
      spa_all_log@x <- log(spa_all@x + 1)
      bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst, spa_eff) {
        sapply(bin_lst, function(bin_i) {
          Matrix::rowSums((spa_eff & bin_i) * (spa_ribo * spa_prop))
        })
      }

      # Log-transform ORF length and aggregated read counts
      feature_mat <- cbind(
        log(orf_aa),
        log(1 + Matrix::rowSums(spa_all)),
        log(1 + Matrix::rowSums(spa_frame0)),
        log(1 + Matrix::rowSums(spa_frame1)),
        log(1 + Matrix::rowSums(spa_frame2)),
        bin_eff(bin_lst, spa_eff),
        bin_prop_sum(spa_frame0_fq, bin_lst, spa_eff),
        bin_prop_sum(spa_frame1_fq, bin_lst, spa_eff),
        bin_prop_sum(spa_frame2_fq, bin_lst, spa_eff),
        bin_weight_sum(spa_all_log, spa_frame0_fq, bin_lst, spa_eff),
        bin_weight_sum(spa_all_log, spa_frame1_fq, bin_lst, spa_eff),
        bin_weight_sum(spa_all_log, spa_frame2_fq, bin_lst, spa_eff)
      )
      return(list(
        feature_mat = feature_mat,
        orf_id = orf_id
      ))
    })

    feature_mat_all <- do.call(cbind, lapply(feat_lst, `[[`, "feature_mat"))
    return(
      list(
        feature_mat = feature_mat_all,
        candi_orf = candi_orf[feat_lst[[1]]$orf_id]
      )
    )
  }

  classify_orf_relative_to_cds <- function(candi_ribo_trans, tx_info) {
    # 1) Build transcript intervals
    rg_tx <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p5,
      end   = tx_info$mix_tx_pos$utr3_p3
    )

    # 2) Match candidate ORFs to their transcripts
    hits <- IRanges::findOverlaps(candi_ribo_trans, rg_tx, select = "first")
    n <- length(candi_ribo_trans)

    # Default everything to Variant
    orftype <- rep("Variant", n)

    ok <- !is.na(hits)
    if (any(ok)) {
      # 3) Only compute coordinates for the matched entries
      s <- IRanges::start(candi_ribo_trans)[ok]
      e <- IRanges::end(candi_ribo_trans)[ok]

      # Define UTR/CDS boundaries
      bp5 <- tx_info$mix_tx_pos$utr5_p3[hits[ok]] + 1L # CDS start
      bp3 <- tx_info$mix_tx_pos$utr3_p5[hits[ok]] - 1L # CDS end

      # 4) Four region indicators (-1/0/1)
      idx11 <- ifelse(s < bp5, -1L, ifelse(s == bp5, 0L, 1L))
      idx12 <- ifelse(s < bp3, -1L, 1L) # no 0 branch
      idx21 <- ifelse(e < bp5, -1L, 1L) # no 0 branch
      idx22 <- ifelse(e < bp3, -1L, ifelse(e == bp3, 0L, 1L))

      key <- paste0(idx11, idx12, idx21, idx22)

      # 5) Lookup table keeps the rules centralized and easy to read
      lut <- c(
        "-1-1-1-1" = "uORF",
        "-1-11-1"  = "uoORF",
        "-1-110"   = "Ext",
        "0-110"    = "CDS",
        "1-11-1"   = "intORF",
        "1-110"    = "Trunc",
        "1-111"    = "doORF",
        "1111"     = "dORF"
      )

      mapped <- lut[match(key, names(lut))]
      orftype[ok] <- ifelse(is.na(mapped), "Variant", mapped)
    }

    # Fix the factor levels for downstream consistency
    factor(
      orftype,
      levels = c("CDS", "uORF", "uoORF", "intORF", "Ext", "Trunc", "doORF", "dORF", "Variant")
    )
  }

  map_group <- function(t) {
    ifelse(t %in% c("Variant", "Ext", "CDS", "Trunc"), "G1_main",
      ifelse(t %in% c("uORF", "dORF"), "G2_u_d",
        ifelse(t %in% c("uoORF"), "G3_uo",
          ifelse(t %in% c("doORF"), "G3_do",
            ifelse(t %in% c("intORF"), "G4_int", "G_other")
          )
        )
      )
    )
  }

  train_rf <- function(translation_features, simulate_orf, tx_info,
                       p_train = 0.8,
                       threshold_method = c("pr_f1", "youden"),
                       target_precision = NULL,
                       w_clip = 50,
                       # RF hyperparams (feel free to tune)
                       num.trees = 1000,
                       mtry = NULL, # default: sqrt(p)
                       min.node.size = 10,
                       sample.fraction = 0.8,
                       max.depth = NULL,
                       splitrule = "gini",
                       seed = 1) {
    threshold_method <- match.arg(threshold_method)

    # ---- Build labels (y) from simulation tags ----
    orf_end <- Biostrings::end(translation_features$candi_orf)
    v_type <- integer(length(translation_features$candi_orf))
    v_type[orf_end %in% (Biostrings::end(simulate_orf$rg_norf_pos) + 2L)] <- 11L
    v_type[orf_end %in% (Biostrings::end(simulate_orf$rg_norf_neg) + 2L)] <- 22L
    v_type[orf_end %in% Biostrings::end(simulate_orf$rg_cds_pos)] <- 1L
    v_type[orf_end %in% Biostrings::end(simulate_orf$rg_cds_neg)] <- 2L

    y <- v_type
    y[y %in% c(2L, 22L, 0L)] <- 0L
    y[y %in% c(1L, 11L)] <- 1L
    y <- as.integer(y)

    # ---- ORF type & super-class features ----
    tmporf_class <- classify_orf_relative_to_cds(translation_features$candi_orf, tx_info)

    map_super <- function(t) {
      ifelse(t %in% c("Variant", "Ext", "CDS", "Trunc"), "G1_main",
        ifelse(t %in% c("uORF", "dORF"), "G2_u_d",
          ifelse(t %in% c("uoORF"), "G3_uo",
            ifelse(t %in% c("doORF"), "G3_do",
              ifelse(t %in% c("intORF"), "G4_int", "G_other")
            )
          )
        )
      )
    }
    super_class <- map_super(as.character(tmporf_class))

    # One-hot for super-class; bind numeric features
    X_grp <- Matrix::sparse.model.matrix(~ factor(super_class,
      levels = c("G1_main", "G2_u_d", "G3_uo", "G3_do", "G4_int", "G_other")
    ) - 1)
    colnames(X_grp) <- c("G1_main", "G2_u_d", "G3_uo", "G3_do", "G4_int", "G_other")

    feature_mat <- translation_features$feature_mat
    if (is.null(colnames(feature_mat))) {
      colnames(feature_mat) <- sprintf("feat_%03d", seq_len(ncol(feature_mat)))
    }
    X <- cbind(X_grp, feature_mat)

    # ---- Stratified split by super_class x label ----
    strata <- paste0(super_class, "__", y)
    lv <- unique(strata)
    idx_train <- integer(0)
    idx_val <- integer(0)
    for (s in lv) {
      ix <- which(strata == s)
      m <- length(ix)
      if (m >= 2) {
        n_tr <- max(1L, floor(m * p_train))
        sel <- sample(ix, m)
        idx_train <- c(idx_train, sel[seq_len(n_tr)])
        if (m > n_tr) idx_val <- c(idx_val, sel[(n_tr + 1L):m])
      } else if (m == 1) {
        idx_train <- c(idx_train, ix)
      }
    }
    idx_train <- sample(idx_train)
    idx_val <- sample(idx_val)

    # ---- Per-type positive sample weights on TRAIN only ----
    dt <- data.table::data.table(y = y[idx_train], type = as.character(tmporf_class[idx_train]))
    tab <- dt[, .(pos = sum(y == 1), neg = sum(y == 0)), by = type]
    eps <- 1e-6
    tab[, spw := (neg + eps) / (pos + eps)]
    tab[spw > w_clip, spw := w_clip]
    spw_map <- setNames(tab$spw, tab$type)

    w_train <- ifelse(y[idx_train] == 1, spw_map[as.character(tmporf_class[idx_train])], 1)
    w_val <- rep(1, length(idx_val)) # keep validation unbiased

    # ---- Build data.frames for ranger ----
    # ranger expects a data.frame and factor response for classification
    Xtr <- as.matrix(X[idx_train, , drop = FALSE])
    Xva <- as.matrix(X[idx_val, , drop = FALSE])

    train_df <- data.frame(
      y = factor(y[idx_train], levels = c(0, 1)),
      as.data.frame(Xtr, check.names = FALSE)
    )
    val_df <- data.frame(
      y = factor(y[idx_val], levels = c(0, 1)),
      as.data.frame(Xva, check.names = FALSE)
    )

    # default mtry if not provided
    if (is.null(mtry)) mtry <- max(1L, floor(sqrt(ncol(train_df) - 1L)))

    # ---- Train Random Forest (ranger) ----
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      sample.fraction = sample.fraction,
      max.depth = max.depth,
      splitrule = splitrule,
      probability = TRUE, # get class probabilities
      case.weights = w_train, # per-instance weights (train only)
      respect.unordered.factors = "partition",
      importance = "impurity",
      seed = seed,
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    # ---- Validate & threshold selection ----
    # predicted prob of class "1"
    p_val <- predict(rf, data = val_df)$predictions[, "1"]
    y_val <- as.integer(as.character(val_df$y))

    .bin_metrics <- function(y_true, p, thr) {
      out <- lapply(thr, function(t) {
        pred <- as.integer(p >= t)
        tp <- sum(pred == 1 & y_true == 1)
        fp <- sum(pred == 1 & y_true == 0)
        tn <- sum(pred == 0 & y_true == 0)
        fn <- sum(pred == 0 & y_true == 1)
        prec <- if ((tp + fp) == 0) 1 else tp / (tp + fp)
        rec <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
        tpr <- rec
        fpr <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
        f1 <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
        c(
          tp = tp, fp = fp, tn = tn, fn = fn,
          precision = prec, recall = rec, tpr = tpr, fpr = fpr, f1 = f1
        )
      })
      m <- do.call(rbind, out)
      data.frame(threshold = thr, m, row.names = NULL, check.names = FALSE)
    }

    thr_cand <- sort(unique(c(p_val, seq(0, 1, by = 0.001))), decreasing = TRUE)
    mtab <- .bin_metrics(y_val, p_val, thr_cand)
    # pick threshold
    if (!is.null(target_precision)) {
      ok <- which(mtab$precision >= target_precision)
      if (length(ok) > 0) {
        pick <- ok[length(ok)]
        t_opt <- mtab$threshold[pick]
        crit <- list(method = "precision>=target", target_precision = target_precision)
      } else {
        pick <- which.max(mtab$f1)
        t_opt <- mtab$threshold[pick]
        crit <- list(method = "fallback_maxF1")
      }
    } else if (threshold_method == "youden") {
      j <- mtab$tpr - mtab$fpr
      pick <- which.max(j)
      t_opt <- mtab$threshold[pick]
      crit <- list(method = "youden", J = max(j, na.rm = TRUE))
    } else { # "pr_f1"
      pick <- which.max(mtab$f1)
      t_opt <- mtab$threshold[pick]
      crit <- list(method = "maxF1", F1 = mtab$f1[pick])
    }

    sel <- which(mtab$threshold == t_opt)[1]

    # ---- AUCs on validation (for quick sanity check) ----
    trapz <- function(x, y) {
      o <- order(x)
      x <- x[o]
      y <- y[o]
      sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    }
    ap <- function(recall, precision) {
      o <- order(recall)
      r <- recall[o]
      p <- precision[o]
      sum(p * c(r[1], diff(r)))
    }
    auroc <- trapz(mtab$fpr, mtab$tpr)
    aupr_trapz <- trapz(mtab$recall, mtab$precision)
    aupr_ap <- ap(mtab$recall, mtab$precision)

    val_summary <- list(
      threshold = t_opt,
      precision = mtab$precision[sel],
      recall = mtab$recall[sel],
      f1 = mtab$f1[sel],
      tpr = mtab$tpr[sel],
      fpr = mtab$fpr[sel],
      counts = unclass(mtab[sel, c("tp", "fp", "tn", "fn")]),
      criterion = crit,
      auroc = auroc,
      aupr_trapz = aupr_trapz,
      aupr_ap = aupr_ap,
      best_iteration = NA_integer_,
      algorithm = "ranger_random_forest"
    )

    list(
      rf            = rf,
      t_opt         = t_opt,
      val_metrics   = val_summary,
      val_curve     = mtab,
      idx_train     = idx_train,
      idx_val       = idx_val,
      y             = y,
      v_type        = v_type,
      orf_class     = tmporf_class,
      super_class   = super_class
    )
  }

  message("Start extract translation features for simulated ribo-seq data ...")
  train_translation_features <- extract_translation_features(
    ribo_num = ribo_seq_simulate,
    candi_translat = candidate_translate_region,
    tx_info = tx_info,
    min_pos_num = min_reads_number,
    eff_ribo_num = active_codon_number,
    bin_num = n_cds_states
  )

  message("Start train translation model...")
  train_model <- train_rf(
    translation_features = train_translation_features,
    simulate_orf = simulate_orf_info,
    tx_info = tx_info,
    threshold_method = threshold_method
  )
  message("Start computing codon-level ribosome activity (MLE estimation)...")
  # 真实数据进行特征提取（单独列在外面）
  ribo_num_lst <- lapply(train_iter_num, function(train_i) {
    pred_p(
      par_rpf = par_rpf_real,
      tx_info = tx_info,
      iter_num = train_i
    )
  })
  message("Start extract translation features for real ribo-seq data...")
  translation_features <- extract_translation_features(
    ribo_num = ribo_num_lst,
    candi_translat = candidate_translate_region,
    tx_info = tx_info,
    min_pos_num = min_reads_number,
    eff_ribo_num = active_codon_number,
    bin_num = n_cds_states
  )

  orf_type <- classify_orf_relative_to_cds(
    translation_features$candi_orf,
    tx_info
  )

  group5 <- map_group(as.character(orf_type))
  # One-hot for fine 5-type (drop one level to avoid perfect collinearity)
  X_grp5 <- Matrix::sparse.model.matrix(~ factor(group5,
    levels = c("G1_main", "G2_u_d", "G3_uo", "G3_do", "G4_int", "G_other")
  ) - 1)
  feat <- translation_features$feature_mat
  if (is.null(colnames(feat))) {
    colnames(feat) <- sprintf("feat_%03d", seq_len(ncol(feat)))
  }
  colnames(X_grp5) <- c("G1_main", "G2_u_d", "G3_uo", "G3_do", "G4_int", "G_other")
  X_new <- cbind(X_grp5, feat)

  # 1) 取出训练时用到的特征列名（最稳妥）
  vars <- train_model$rf$forest$independent.variable.names

  # 2) 对齐新矩阵：补缺列为0、去多余列、按vars重排
  X_new_df <- as.data.frame(as.matrix(X_new))
  miss <- setdiff(vars, colnames(X_new_df))
  if (length(miss)) X_new_df[miss] <- 0
  X_new_df <- X_new_df[, vars, drop = FALSE]

  # 3) 预测正类概率（"1"）
  y_pred_prob <- predict(train_model$rf, data = X_new_df)$predictions[, "1"]
  # 4) 用训练阶段挑的阈值二分类
  # y_pred_prob <- predict(train_model$bst, X_new, iterationrange = c(1, train_model$bst$best_iteration))

  return(list(
    train_model = train_model,
    y_pred_prob = y_pred_prob,
    orf_type = orf_type,
    candidate_rg = translation_features$candi_orf,
    ribo_num = ribo_num_lst
  ))
}



cut_prob_plot <- function(par0) {
  base_prob <- mean(par0$cut_bias$s7)
  p5 <- par0$prob_hd$p5
  p3 <- par0$prob_hd$p3
  pmf5_raw <- rev(c(1, cumprod(rev(base_prob^p5))[-length(p5)]) * (1 - rev(base_prob^p5)))
  pmf3_raw <- c(1, cumprod(base_prob^p3)[-length(p3)]) * (1 - (base_prob^p3))
  pmf5 <- pmf5_raw
  pmf3 <- pmf3_raw
  return(cbind(pmf5, pmf3))
}

infer_ribo_num_train <- function(
    par_lst,
    tx_info,
    is_mnase = FALSE) {
  par_rpf <- par_lst
  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[
    par_rpf$lst_rpf$rpf_overlap$qwidth %in% unique(convert_idx[4, ]),
  ]

  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = psite_pos
  )
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  df_rpf <- par_rpf$lst_rpf$rpf_overlap
  par_0 <- par_rpf$par_0
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )

  frame_sum <- lapply(par_rpf$lst_rpf$offsets[1, ], function(len_i) {
    idx_i <- df_rpf$qwidth == len_i
    candi_psite <- candi_psite[idx_i, ]
    candi_cut5 <- candi_cut5[idx_i, ]
    candi_cut3 <- candi_cut3[idx_i, ]
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    if (is_mnase == TRUE) {
      cut_seq <- get_cleavage_seq(
        seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    } else {
      cut_seq <- get_cleavage_seq(
        seqs = tx_info$mix_tx, p_site = uni_candi_psite[1:2], ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      cut_prob$final_p5 <- matrix(
        rep(cut_prob$final_p5[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p5), byrow = T
      )
      cut_prob$final_p3 <- matrix(
        rep(cut_prob$final_p3[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p3), byrow = T
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    }
    M <- candi_p_weight
    sx <- Matrix::summary(M)
    ord <- order(sx$i, -sx$x)
    sx2 <- sx[ord, ]
    sx2 <- sx2[!duplicated(sx2$i), ]
    sx$x <- 0L
    sx$x[as.numeric(rownames(sx2))] <- 1L
    M_argmax <- Matrix::sparseMatrix(
      i    = sx$i,
      j    = sx$j,
      x    = sx$x,
      dims = dim(M)
    )
    sparse_iter <- Matrix::sparseMatrix(
      i = shrink_pos,
      j = index_j,
      x = (df_rpf[idx_i, ]$weight * M_argmax)@x
    )
    vec_pnum[] <- Matrix::rowSums(sparse_iter)
    ribo_num <- list(
      candi_pos = candi_pos,
      vec_pnum = vec_pnum
    )
    tmp_psite_gr <- IRanges::IRanges(start = ribo_num$candi_pos, width = 1L)
    idx_cds_psite <- IRanges::findOverlaps(
      query = tmp_psite_gr,
      subject = rg_cds,
      minoverlap = 1L
    )
    psite_cds <- tmp_psite_gr@start[idx_cds_psite@from] - rg_cds@start[idx_cds_psite@to]
    cds_frame_pnum <- split(ribo_num$vec_pnum[idx_cds_psite@from], psite_cds %% 3L)
    frame_sum <- sapply(cds_frame_pnum, sum)
    print(frame_sum / sum(frame_sum))
    frame_prob = sort(frame_sum / sum(frame_sum), decreasing = T)
len_weight = (log(frame_prob[1] / (1 / 3)) - log(frame_prob[2] / (1 / 3)))
    frame_max <- which.max(frame_sum)
    if (frame_max == 1) {
      ribo_num$candi_pos <- ribo_num$candi_pos + 1L
    } else if (frame_max == 3) {
      ribo_num$candi_pos <- ribo_num$candi_pos - 1L
    }
    if(len_weight < 0.1) {
      ribo_num$vec_pnum <- ribo_num$vec_pnum * 0
    }
    print(sum(ribo_num$vec_pnum))
    #
    tmp_psite_gr <- IRanges::IRanges(start = ribo_num$candi_pos, width = 1L)
    idx_cds_psite <- IRanges::findOverlaps(
      query = tmp_psite_gr,
      subject = rg_cds,
      minoverlap = 1L
    )
    psite_cds <- tmp_psite_gr@start[idx_cds_psite@from] - rg_cds@start[idx_cds_psite@to]
    cds_frame_pnum <- split(ribo_num$vec_pnum[idx_cds_psite@from], psite_cds %% 3L)
    frame_sum <- sapply(cds_frame_pnum, sum)
    print(len_weight)
    return(list(
      ribo_num = ribo_num,
      frame_sum = frame_sum
    ))
  })
  ribo_num_pos <- colSums(do.call(rbind, lapply(frame_sum, function(x) {
    ribo_num_0 <- vector("numeric", length = length(tx_info_lst$tx_info$mix_tx))
    ribo_num_0[x$ribo_num$candi_pos] <- x$ribo_num$vec_pnum
    return(ribo_num_0)
  })))
  return(ribo_num_pos)
}

shift_psite_new <- function(tx_info,
                            rpf_info,
                            offsets,
                            shift_pos,
                            ribo_num = NULL) {
  # create negative RPF
  rpf_info_neg <- rpf_info
  rpf_info_neg$pos <- rpf_info$pos + sample(shift_pos, nrow(rpf_info), replace = TRUE)

  # initial p-site number (POS)
  p_table <- base::table(
    offsets[2, ][match(rpf_info$qwidth, offsets[1, ])] +
      rpf_info$pos
  )
  # NEG
  p_table_neg <- base::table(
    offsets[2, ][match(rpf_info_neg$qwidth, offsets[1, ])] +
      rpf_info_neg$pos
  )

  p_pos <- p_pos_neg <- integer(length(tx_info$mix_tx))
  p_pos[as.integer(names(p_table))] <- as.vector(p_table)
  p_pos_neg[as.integer(names(p_table_neg))] <- as.vector(p_table_neg)

  choose_p <- which(p_pos > 0L)
  choose_p_neg <- which(p_pos_neg > 0L)

  ribo_num_pos <- list(candi_pos = choose_p, vec_pnum = p_pos[choose_p])
  ribo_num_neg <- list(candi_pos = choose_p_neg, vec_pnum = p_pos_neg[choose_p_neg])
  if (is.null(ribo_num)) {
    return(list(
      ribo_num = ribo_num_pos,
      ribo_num_neg = ribo_num_neg
    ))
  } else {
    p_pos <- ribo_num
    choose_p <- which(p_pos > 0L)
    ribo_num_pos <- list(candi_pos = choose_p, vec_pnum = p_pos[choose_p])
    return(list(
      ribo_num = ribo_num_pos,
      ribo_num_neg = ribo_num_neg
    ))
  }
}

run_one_enzyme_new <- function(enzyme,
                               par_pos,
                               par_neg,
                               picked_file,
                               tx_info,
                               candi_translat,
                               ribo_num_pos,
                               ribo_num_neg,
                               is_mnase = FALSE) {
  # ---------- 1) offsets ----------
  off_pos <- par_pos$lst_rpf$offsets
  off_neg <- par_neg$lst_rpf$offsets

  if (is_mnase) {
    # MNase：只保留 23,26,29,32,35 的 RPF 长度
    keep_len <- c("23", "26", "29", "32", "35")
    off_pos <- off_pos[, colnames(off_pos) %in% keep_len, drop = FALSE]
    off_neg <- off_neg[, colnames(off_neg) %in% keep_len, drop = FALSE]
  }
  # RNase I / P1：使用所有 offsets，无需过滤

  # ---------- 2) 生成 P-site 计数（POS/NEG 各一套） ----------
  psite_pos <- shift_psite_new(
    tx_info = tx_info,
    rpf_info = par_pos$lst_rpf$rpf_info,
    offsets = off_pos,
    shift_pos = c(-3L, -2L, -1L, 0L, 1L, 2L, 3L),
    ribo_num = ribo_num_pos
  )

  psite_neg <- shift_psite_new(
    tx_info = tx_info,
    rpf_info = par_neg$lst_rpf$rpf_info,
    offsets = off_neg,
    shift_pos = c(-3L, -2L, -1L, 0L, 1L, 2L, 3L),
    ribo_num = ribo_num_neg
  )

  # ---------- 3) 训练 POS 模型 ----------
  model_pos <- train_features(
    ribo_num_train = psite_pos,
    candi_translat = candi_translat,
    tx_info        = tx_info,
    min_pos_num    = 20,
    eff_ribo_num   = 0.1
  )

  lab_pos <- predict(
    model_pos$train_model$modle_cds,
    data = data.frame(model_pos$feat_m)
  )$predictions[, "1"]

  # ---------- 4) 训练 NEG 模型 ----------
  model_neg <- train_features(
    ribo_num_train = psite_neg,
    candi_translat = candi_translat,
    tx_info        = tx_info,
    min_pos_num    = 20,
    eff_ribo_num   = 0.1
  )

  lab_neg <- predict(
    model_neg$train_model$modle_cds,
    data = data.frame(model_neg$feat_m)
  )$predictions[, "1"]

  # ---------- 5) 只保留 1000 个目标基因 ----------
  picked_genes <- readr::read_tsv(picked_file, show_col_types = FALSE)

  pos_df <- ba_tr(
    train_model = model_pos,
    lab_score = lab_pos,
    tx_info = tx_info,
    picked_genes = picked_genes
  )
  neg_df <- ba_tr(
    train_model = model_neg,
    lab_score = lab_neg,
    tx_info = tx_info,
    picked_genes = picked_genes
  )

  # ---------- 6) 组装成 precrec 输入，1=POS, 0=NEG ----------
  df <- bind_rows(
    pos_df %>% mutate(y = 1L),
    neg_df %>% mutate(y = 0L)
  ) %>%
    filter(is.finite(pred_score), !is.na(y))

  # ---------- 7) evalmod ----------
  res <- precrec::evalmod(
    scores = df$pred_score,
    labels = df$y,
    modnames = enzyme
  )

  list(
    enzyme = enzyme,
    res = res,
    data = df,
    model_pos = model_pos,
    model_neg = model_neg
  )
}


train_features2 <- function(ribo_num_train,
                           candi_translat,
                           tx_info,
                           min_pos_num,
                           eff_ribo_num) {
  # 从spa_frame得到feature_mat
  get_feature <- function(spa_lst) {
    # Summaries within positional bins
    bin_eff <- function(spa_eff) {
      Matrix::rowSums(spa_eff)
    }
    bin_prop_sum <- function(spa_prop, spa_eff) {
      Matrix::rowSums(spa_prop * spa_eff)
    }
    spa_all_log <- spa_lst$spa_all
    spa_all_log@x <- log(spa_lst$spa_all@x + 1)
    bin_weight_sum <- function(spa_ribo, spa_prop, spa_eff) {
      Matrix::rowSums(spa_eff * (spa_ribo * spa_prop))
    }
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      log(spa_lst$orf_aa),
      log(1 + Matrix::rowSums(spa_lst$spa_all)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame0)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame1)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame2)),
      bin_eff(spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, spa_lst$spa_eff)
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    return(feature_mat)
  }
  # Project ribosome footprints onto candidate translated segments
  rg_p <- list(
    rg_p_pos = IRanges::IRanges(
      start = ribo_num_train$ribo_num$candi_pos,
      width = 1L
    ),
    rg_p_neg = IRanges::IRanges(
      start = ribo_num_train$ribo_num_neg$candi_pos,
      width = 1L
    )
  )

  # 取出CDS
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  candi_orf <- do.call(c, candi_translat)
  candi_orf@start <- candi_orf@start + 3L
  candi_orf@width <- candi_orf@width - 1L
  candi_all_orf <- candi_orf
  # 取出CDS所在的stop-stop区域
  cds_idx <- end(candi_orf) %in% end(rg_cds)
  candi_orf <- candi_orf[cds_idx]
  # 得到所有可能翻译区域的特征
  spa_frame_all <- lapply(1, function(i) {
    hit_p <- IRanges::findOverlaps(rg_p[[i]], candi_all_orf)
    hit_p <- hit_p[
      hit_p@to %in%
        which(
          IRanges::countOverlaps(candi_all_orf, rg_p[[i]]) >= min_pos_num
        )
    ]
    dist_p <- rg_p[[i]]@start[hit_p@from] - candi_all_orf@start[hit_p@to]
    max_dist <- max(dist_p) %/% 3L + 3L
    idx_frame <- dist_p %% 3L
    idx_pos <- (dist_p %/% 3L) + 1L
    orf_n <- as.factor(hit_p@to)
    idx_orf <- as.integer(orf_n)

    p_num <- ribo_num_train[[i]]$vec_pnum

    # Build sparse codon-by-codon count matrices, one ORF per row
    spa_frame0 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 0L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 0L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 0L], 1)
    )
    spa_frame0 <- spa_frame0[-nrow(spa_frame0), ]
    spa_frame1 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 1L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 1L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 1L], 1)
    )
    spa_frame1 <- spa_frame1[-nrow(spa_frame1), ]
    spa_frame2 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 2L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 2L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 2L], 1)
    )
    spa_frame2 <- spa_frame2[-nrow(spa_frame2), ]

    # Normalize ribosome counts across frames for each codon
    orf_id <- as.integer(
      levels(as.factor(hit_p@to))
    )
    spa0 <- spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_aa <- Biostrings::width(candi_all_orf[orf_id]) / 3
    spa_eff <- spa_all > eff_ribo_num
    return(list(
      orf_id = orf_id,
      spa_eff = spa_eff,
      orf_aa = orf_aa,
      spa_all = spa_all,
      spa_frame0 = spa_frame0,
      spa_frame1 = spa_frame1,
      spa_frame2 = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    ))
  })
  candi_all_orf = candi_all_orf[spa_frame_all[[1]]$orf_id]
  print(length(candi_all_orf))
  feature_mat_all <- get_feature(spa_lst = spa_frame_all[[1]])

  spa_frame_lst <- lapply(1:2, function(i) {
    hit_p <- IRanges::findOverlaps(rg_p[[i]], candi_orf)
    hit_p <- hit_p[
      hit_p@to %in%
        which(
          IRanges::countOverlaps(candi_orf, rg_p[[i]]) >= min_pos_num
        )
    ]
    dist_p <- rg_p[[i]]@start[hit_p@from] - candi_orf@start[hit_p@to]
    max_dist <- max(dist_p) %/% 3L + 3L
    idx_frame <- dist_p %% 3L
    idx_pos <- (dist_p %/% 3L) + 1L
    orf_n <- as.factor(hit_p@to)
    idx_orf <- as.integer(orf_n)

    p_num <- ribo_num_train[[i]]$vec_pnum

    # Build sparse codon-by-codon count matrices, one ORF per row
    spa_frame0 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 0L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 0L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 0L], 1)
    )
    spa_frame0 <- spa_frame0[-nrow(spa_frame0), ]
    spa_frame1 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 1L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 1L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 1L], 1)
    )
    spa_frame1 <- spa_frame1[-nrow(spa_frame1), ]
    spa_frame2 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 2L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 2L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 2L], 1)
    )
    spa_frame2 <- spa_frame2[-nrow(spa_frame2), ]

    # Normalize ribosome counts across frames for each codon
    orf_id <- as.integer(
      levels(as.factor(hit_p@to))
    )
    spa0 <- spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_aa <- Biostrings::width(candi_orf[orf_id]) / 3
    spa_eff <- spa_all > eff_ribo_num
    return(list(
      orf_id = orf_id,
      spa_eff = spa_eff,
      orf_aa = orf_aa,
      spa_all = spa_all,
      spa_frame0 = spa_frame0,
      spa_frame1 = spa_frame1,
      spa_frame2 = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    ))
  })
  # 首先判断阳性CDS中不是frame1占优势的部分
  spa_frame_p1 <- spa_frame_lst[[1]]
  spa_frame_neg <- spa_frame_lst[[2]]

  feature_mat_p1 <- get_feature(spa_lst = spa_frame_p1)
  feature_mat_neg <- get_feature(spa_lst = spa_frame_neg)

  change_pos_lable <- function(feature_mat_p1, feature_mat_neg) {
    X <- rbind(feature_mat_p1, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p1)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf, data = data.frame(feature_mat_p1))$predictions[, "0"]
    feature_mat_p2 <- feature_mat_p1[p_val < 0.5, ]
    X <- rbind(feature_mat_p2, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p2)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf_2 <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf_2, data = data.frame(feature_mat_p1))$predictions[, "0"]
    return(p_val)
  }

  prob_cds_false <- change_pos_lable(feature_mat_p1 = feature_mat_p1, feature_mat_neg = feature_mat_neg)

  # 得到清洗过后的CDS frame阳性和阴性集
  spa_frame_pos <- lapply(spa_frame_p1[c(1, 3)], function(x) x[prob_cds_false < 0.7])
  spa_frame_pos <- c(spa_frame_pos, lapply(spa_frame_p1[c(-1, -3)], function(x) x[prob_cds_false < 0.7, ]))
  # 得出错位的阴性情况
  spa_frame_ahead <- spa_frame_behind <- spa_frame_pos
  spa_frame_ahead$spa_frame0 <- spa_frame_pos$spa_frame1
  spa_frame_ahead$spa_frame1 <- spa_frame_pos$spa_frame2
  spa_frame_ahead$spa_frame2 <- spa_frame_pos$spa_frame0
  spa_frame_ahead$spa_frame0_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_ahead$spa_frame1_fq <- spa_frame_pos$spa_frame2_fq
  spa_frame_ahead$spa_frame2_fq <- spa_frame_pos$spa_frame0_fq

  spa_frame_behind$spa_frame2 <- spa_frame_pos$spa_frame1
  spa_frame_behind$spa_frame1 <- spa_frame_pos$spa_frame0
  spa_frame_behind$spa_frame0 <- spa_frame_pos$spa_frame2
  spa_frame_behind$spa_frame2_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_behind$spa_frame1_fq <- spa_frame_pos$spa_frame0_fq
  spa_frame_behind$spa_frame0_fq <- spa_frame_pos$spa_frame2_fq
  # 训练CDS翻译模型
  feat_m_pos <- get_feature(spa_lst = spa_frame_pos)
  feat_m_ahead <- get_feature(spa_lst = spa_frame_ahead)
  feat_m_behind <- get_feature(spa_lst = spa_frame_behind)
train_model_cds <- function(m_pos, m_neg,
                            threshold_method = c("maxF1", "youden"),
                            target_precision = NULL,
                            p_train = 0.8,
                            seed = 1L) {

  threshold_method <- match.arg(threshold_method)

  # -------- 组装数据 --------
  X <- rbind(m_pos, m_neg)
  y <- c(rep(1L, nrow(m_pos)), rep(0L, nrow(m_neg)))
  y <- as.integer(y)

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  set.seed(seed)

  # -------- 分层划分 train / val --------
  n_pos <- nrow(m_pos)
  n_neg <- nrow(m_neg)

  idx_pos <- sample(seq_len(n_pos))
  idx_neg <- sample(seq_len(n_neg))

  ntr_pos <- max(1L, floor(n_pos * p_train))
  ntr_neg <- max(1L, floor(n_neg * p_train))

  train_pos_idx <- idx_pos[seq_len(ntr_pos)]
  val_pos_idx   <- idx_pos[(ntr_pos + 1L):n_pos]

  train_neg_idx <- idx_neg[seq_len(ntr_neg)]
  val_neg_idx   <- idx_neg[(ntr_neg + 1L):n_neg]

  # 防一下极端情况：验证集为空
  if (length(val_pos_idx) == 0L || length(val_neg_idx) == 0L) {
    warning("验证集某一类为空，改用非分层随机划分。")
    n <- nrow(X)
    idx_all <- sample(seq_len(n))
    ntr <- max(2L, floor(n * p_train))
    idx_tr <- idx_all[seq_len(ntr)]
    idx_va <- idx_all[(ntr + 1L):n]

    X_train <- X[idx_tr, , drop = FALSE]
    y_train <- y[idx_tr]
    X_val   <- X[idx_va, , drop = FALSE]
    y_val   <- y[idx_va]

    split_info <- list(
      train_idx = idx_tr,
      val_idx   = idx_va
    )
  } else {
    # 正负分别划分后再合并
    X_train <- rbind(
      m_pos[train_pos_idx, , drop = FALSE],
      m_neg[train_neg_idx, , drop = FALSE]
    )
    y_train <- c(rep(1L, length(train_pos_idx)),
                 rep(0L, length(train_neg_idx)))

    X_val <- rbind(
      m_pos[val_pos_idx, , drop = FALSE],
      m_neg[val_neg_idx, , drop = FALSE]
    )
    y_val <- c(rep(1L, length(val_pos_idx)),
               rep(0L, length(val_neg_idx)))

    X_train <- as.matrix(X_train)
    X_val   <- as.matrix(X_val)

    split_info <- list(
      train_pos_idx = train_pos_idx,
      train_neg_idx = train_neg_idx,
      val_pos_idx   = val_pos_idx,
      val_neg_idx   = val_neg_idx
    )
  }

  # -------- 训练集样本权重（简单的正负比） --------
  n_pos_tr <- sum(y_train == 1L)
  n_neg_tr <- sum(y_train == 0L)

  if (n_pos_tr == 0L || n_neg_tr == 0L) {
    stop("训练集没有同时包含正负样本，无法训练 xgboost 模型。")
  }

  pos_w <- n_neg_tr / n_pos_tr
  w_train <- ifelse(y_train == 1L, pos_w, 1)

  dtrain <- xgboost::xgb.DMatrix(
    data   = X_train,
    label  = y_train,
    weight = w_train
  )
  dval <- xgboost::xgb.DMatrix(
    data   = X_val,
    label  = y_val
  )

  # -------- xgboost 参数 --------
  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "aucpr",    # 针对不平衡
    eta              = 0.05,
    max_depth        = 6,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    lambda           = 1
  )

  bst <- xgboost::xgb.train(
    params  = params,
    data    = dtrain,
    nrounds = 5000,
    watchlist = list(train = dtrain, val = dval),
    early_stopping_rounds = 100,
    verbose = 1
  )

  best_iter <- bst$best_iteration
  if (is.null(best_iter) || is.na(best_iter)) {
    best_iter <- bst$niter
  }

  # -------- 在验证集上选最佳阈值 --------
  pred_val <- predict(bst, dval, ntreelimit = best_iter)

  .bin_metrics <- function(y_true, p, thr) {
    out <- lapply(thr, function(t) {
      pred <- as.integer(p >= t)
      tp <- sum(pred == 1 & y_true == 1)
      fp <- sum(pred == 1 & y_true == 0)
      tn <- sum(pred == 0 & y_true == 0)
      fn <- sum(pred == 0 & y_true == 1)
      prec <- if ((tp + fp) == 0) 1 else tp / (tp + fp)
      rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
      tpr  <- rec
      fpr  <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
      f1   <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
      c(
        tp = tp, fp = fp, tn = tn, fn = fn,
        precision = prec, recall = rec,
        tpr = tpr, fpr = fpr, f1 = f1
      )
    })
    m <- do.call(rbind, out)
    data.frame(threshold = thr, m, row.names = NULL, check.names = FALSE)
  }

  thr_cand <- sort(unique(c(pred_val, seq(0, 1, by = 0.001))), decreasing = TRUE)
  mtab <- .bin_metrics(y_val, pred_val, thr_cand)

  if (!is.null(target_precision)) {
    ok <- which(mtab$precision >= target_precision)
    if (length(ok) > 0) {
      pick <- ok[length(ok)]
      t_opt <- mtab$threshold[pick]
      crit  <- list(method = "precision>=target",
                    target_precision = target_precision)
    } else {
      pick <- which.max(mtab$f1)
      t_opt <- mtab$threshold[pick]
      crit  <- list(method = "fallback_maxF1")
    }
  } else if (threshold_method == "youden") {
    j <- mtab$tpr - mtab$fpr
    pick <- which.max(j)
    t_opt <- mtab$threshold[pick]
    crit  <- list(method = "youden", J = max(j, na.rm = TRUE))
  } else {
    pick <- which.max(mtab$f1)
    t_opt <- mtab$threshold[pick]
    crit  <- list(method = "maxF1", F1 = mtab$f1[pick])
  }

  sel <- which(mtab$threshold == t_opt)[1]

  val_summary <- list(
    threshold = t_opt,
    precision = mtab$precision[sel],
    recall    = mtab$recall[sel],
    f1        = mtab$f1[sel],
    tpr       = mtab$tpr[sel],
    fpr       = mtab$fpr[sel],
    counts    = unclass(mtab[sel, c("tp", "fp", "tn", "fn")]),
    criterion = crit,
    best_iteration = best_iter,
    brier = mean((pred_val - y_val)^2)
  )

  list(
    model       = bst,
    t_opt       = t_opt,
    val_metrics = val_summary,
    val_curve   = mtab,
    split_info  = split_info
  )
}

  modle_cds <- train_model_cds(
    m_pos = feat_m_pos,
    m_neg = rbind(feature_mat_neg, feat_m_ahead, feat_m_behind),
    threshold_method = "maxF1"
  )
  # 训练uORF dORF lncORF翻译模型
  pick_cds_spa <- function(spa_lst, aa_min = 6, aa_max = 100, bin_num = 10, rpf_min = 10, rpf_max = 200, eff_pos_min = 3) {
    orf_bin <- Matrix::summary(spa_lst$spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / spa_lst$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_lst$spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_lst$spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- spa_lst$orf_aa / bin_num
    spa_all_log <- spa_lst$spa_all
    spa_all_log@x <- log(spa_lst$spa_all@x + 1)
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(log(norf_aa), bin_num),
      log(1 + rpf_num),
      log(1 + as.vector(bin_f0_num(bin_lst, spa_lst$spa_frame0))),
      log(1 + as.vector(bin_f1_num(bin_lst, spa_lst$spa_frame1))),
      log(1 + as.vector(bin_f2_num(bin_lst, spa_lst$spa_frame2))),
      eff_pos,
      as.vector(bin_prop_sum(spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    # 过滤nORF标准
    shift_k <- function(v, k) {
      len_v <- length(v)
      v_tp <- v
      v_tp[1:(len_v / k)] <- v[(len_v * (k - 1) / k + 1):len_v]
      v_tp[(len_v / k + 1):len_v] <- v[1:(len_v * (k - 1) / k)]
      return(v_tp)
    }
    idx_norf <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= rpf_min) & (rpf_num <= rpf_max) & (eff_pos >= eff_pos_min)
    idx_cds <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= 30) & (eff_pos >= 10)
    idx_norf_change <- shift_k(idx_norf, bin_num)
    idx_cds_change <- shift_k(idx_cds, bin_num)
    # 转化坐标
    change_cord <- function(frame_i, aa_vec, bin_num) {
      orf_bin_f0 <- Matrix::summary(frame_i)
      bin_idx <- floor((orf_bin_f0$j - 0.5) / aa_vec[orf_bin_f0$i] * bin_num) + 1
      idx_tmp <- bin_idx < bin_num
      step_i <- floor(aa_vec[orf_bin_f0$i] / bin_num)
      orf_bin_f0$j[idx_tmp] <- orf_bin_f0$j[idx_tmp] + step_i[idx_tmp]
      orf_bin_f0$j[!idx_tmp] <- orf_bin_f0$j[!idx_tmp] + step_i[!idx_tmp] - aa_vec[orf_bin_f0$i][!idx_tmp] + 1L
      orf_bin_f0 <- orf_bin_f0[orf_bin_f0$j <= dim(frame_i)[2], ]
      spa_bin_f0 <- Matrix::sparseMatrix(
        i = orf_bin_f0$i,
        j = orf_bin_f0$j,
        x = orf_bin_f0$x,
        dims = dim(frame_i)
      )
      return(spa_bin_f0)
    }
    frame_0 <- change_cord(frame_i = spa_lst$spa_frame0, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_1 <- change_cord(frame_i = spa_lst$spa_frame1, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_2 <- change_cord(frame_i = spa_lst$spa_frame2, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    return(list(
      feature_mat = feature_mat[idx_norf, ],
      idx = list(
        idx_norf = idx_norf,
        idx_norf_change = idx_norf_change,
        idx_cds = idx_cds,
        idx_cds_change = idx_cds_change
      ),
      spa_frame = list(
        frame_0 = frame_0,
        frame_1 = frame_1,
        frame_2 = frame_2
      )
    ))
  }
  feat_norf_pos <- pick_cds_spa(spa_lst = spa_frame_pos)
  feat_norf_neg <- pick_cds_spa(spa_lst = spa_frame_neg)
  feat_norf_ahead <- pick_cds_spa(spa_lst = spa_frame_ahead)
  feat_norf_behind <- pick_cds_spa(spa_lst = spa_frame_behind)

  modle_ud_orf <- train_model_cds(
    m_pos = feat_norf_pos$feature_mat,
   m_neg = rbind(feat_norf_neg$feature_mat, feat_norf_ahead$feature_mat, feat_norf_behind$feature_mat)
  )

  # 训练uoORF doORF iORF翻译模型
  #   背景错框翻译，目标分4种情况
  # 阳性
  pick_overlap <- function(back_frame, head_info, bin_num = 10) {
    spa_frame0 <- back_frame$spa_frame0 + head_info$spa_frame$frame_0
    spa_frame1 <- back_frame$spa_frame1 + head_info$spa_frame$frame_1
    spa_frame2 <- back_frame$spa_frame2 + head_info$spa_frame$frame_2

    spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_bin <- Matrix::summary(spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / back_frame$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- back_frame$orf_aa / bin_num
    spa_all_log <- spa_all
    spa_all_log@x <- log(spa_all@x + 1)
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(log(norf_aa), bin_num),
      log(1 + rpf_num),
      log(1 + as.vector(bin_f0_num(bin_lst, spa_frame0))),
      log(1 + as.vector(bin_f1_num(bin_lst, spa_frame1))),
      log(1 + as.vector(bin_f2_num(bin_lst, spa_frame2))),
      eff_pos,
      as.vector(bin_prop_sum(spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    idx_cds <- (norf_aa >= 6) & (rpf_num >= 40) & (eff_pos >= 13) & head_info$idx$idx_norf
    return(feature_mat[idx_cds, ])
  }

  feat_oorf_p1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_pos)
  feat_oorf_p2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_pos)

  # 阴性
  feat_oorf_n1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_behind)
  feat_oorf_n2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_ahead)

  modle_overlap_orf <- train_model_cds(
    m_pos = rbind(feat_oorf_p1, feat_oorf_p2),
    m_neg = rbind(feat_oorf_n1, feat_oorf_n2)
  )
  return(list(
    train_model = list(
      modle_cds = modle_cds,
      modle_ud_orf = modle_ud_orf,
      modle_overlap_orf = modle_overlap_orf
    ),
    candidate_rg = candi_all_orf,
    ribo_num = ribo_num_train$ribo_num,
    feat_m = feature_mat_all
  ))
}


classify_orf_relative_to_cds <- function(candi_ribo_trans, tx_info) {
      # 1) Build transcript intervals
      rg_tx <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$utr5_p5,
        end   = tx_info$mix_tx_pos$utr3_p3
      )

      # 2) Match candidate ORFs to their transcripts
      hits <- IRanges::findOverlaps(candi_ribo_trans, rg_tx, select = "first")
      n <- length(candi_ribo_trans)

      # Default everything to Variant
      orftype <- rep("Variant", n)

      ok <- !is.na(hits)
      if (any(ok)) {
        # 3) Only compute coordinates for the matched entries
        s <- IRanges::start(candi_ribo_trans)[ok]
        e <- IRanges::end(candi_ribo_trans)[ok]

        # Define UTR/CDS boundaries
        bp5 <- tx_info$mix_tx_pos$utr5_p3[hits[ok]] + 1L # CDS start
        bp3 <- tx_info$mix_tx_pos$utr3_p5[hits[ok]] - 1L # CDS end

        # 4) Four region indicators (-1/0/1)
        idx11 <- ifelse(s < bp5, -1L, ifelse(s == bp5, 0L, 1L))
        idx12 <- ifelse(s < bp3, -1L, 1L) # no 0 branch
        idx21 <- ifelse(e < bp5, -1L, 1L) # no 0 branch
        idx22 <- ifelse(e < bp3, -1L, ifelse(e == bp3, 0L, 1L))

        key <- paste0(idx11, idx12, idx21, idx22)

        # 5) Lookup table keeps the rules centralized and easy to read
        lut <- c(
          "-1-1-1-1" = "uORF",
          "-1-11-1"  = "uoORF",
          "-1-110"   = "Ext",
          "0-110"    = "CDS",
          "1-11-1"   = "intORF",
          "1-110"    = "Trunc",
          "1-111"    = "doORF",
          "1111"     = "dORF"
        )

        mapped <- lut[match(key, names(lut))]
        orftype[ok] <- ifelse(is.na(mapped), "Variant", mapped)
      }

      # Fix the factor levels for downstream consistency
      factor(
        orftype,
        levels = c("CDS", "uORF", "uoORF", "intORF", "Ext", "Trunc", "doORF", "dORF", "Variant")
      )
    }

    map_group <- function(t) {
      ifelse(t %in% c("Variant", "Ext", "CDS", "Trunc"), "G1_main",
        ifelse(t %in% c("uORF", "dORF"), "G2_u_d",
          ifelse(t %in% c("uoORF","doORF","intORF"), "G3_uodo", "G_other")
            )
          )
    }



train_features <- function(ribo_num_train,
                           candi_translat,
                           tx_info,
                           min_pos_num,
                           eff_ribo_num) {
  # 从spa_frame得到feature_mat
  get_feature <- function(spa_lst) {
    # Summaries within positional bins
    bin_eff <- function(spa_eff) {
      Matrix::rowSums(spa_eff)
    }
    bin_prop_sum <- function(spa_prop, spa_eff) {
      Matrix::rowSums(spa_prop * spa_eff)
    }
    spa_all_log <- spa_lst$spa_all
    spa_all_log@x <- log(spa_lst$spa_all@x + 1)
    bin_weight_sum <- function(spa_ribo, spa_prop, spa_eff) {
      Matrix::rowSums(spa_eff * (spa_ribo * spa_prop))
    }
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      log(spa_lst$orf_aa),
      log(1 + Matrix::rowSums(spa_lst$spa_all)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame0)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame1)),
      log(1 + Matrix::rowSums(spa_lst$spa_frame2)),
      bin_eff(spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      bin_prop_sum(spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, spa_lst$spa_eff)
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    return(feature_mat)
  }
  # Project ribosome footprints onto candidate translated segments
  rg_p <- list(
    rg_p_pos = IRanges::IRanges(
      start = ribo_num_train$ribo_num$candi_pos,
      width = 1L
    ),
    rg_p_neg = IRanges::IRanges(
      start = ribo_num_train$ribo_num_neg$candi_pos,
      width = 1L
    )
  )

  # 取出CDS
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  candi_orf <- do.call(c, candi_translat)
  candi_orf@start <- candi_orf@start + 3L
  candi_orf@width <- candi_orf@width - 1L
  candi_all_orf <- candi_orf
  # 取出CDS所在的stop-stop区域
  cds_idx <- end(candi_orf) %in% end(rg_cds)
  candi_orf <- candi_orf[cds_idx]
  # 得到所有可能翻译区域的特征
  spa_frame_all <- lapply(1, function(i) {
    hit_p <- IRanges::findOverlaps(rg_p[[i]], candi_all_orf)
    hit_p <- hit_p[
      hit_p@to %in%
        which(
          IRanges::countOverlaps(candi_all_orf, rg_p[[i]]) >= min_pos_num
        )
    ]
    dist_p <- rg_p[[i]]@start[hit_p@from] - candi_all_orf@start[hit_p@to]
    max_dist <- max(dist_p) %/% 3L + 3L
    idx_frame <- dist_p %% 3L
    idx_pos <- (dist_p %/% 3L) + 1L
    orf_n <- as.factor(hit_p@to)
    idx_orf <- as.integer(orf_n)

    p_num <- ribo_num_train[[i]]$vec_pnum

    # Build sparse codon-by-codon count matrices, one ORF per row
    spa_frame0 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 0L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 0L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 0L], 1)
    )
    spa_frame0 <- spa_frame0[-nrow(spa_frame0), ]
    spa_frame1 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 1L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 1L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 1L], 1)
    )
    spa_frame1 <- spa_frame1[-nrow(spa_frame1), ]
    spa_frame2 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 2L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 2L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 2L], 1)
    )
    spa_frame2 <- spa_frame2[-nrow(spa_frame2), ]

    # Normalize ribosome counts across frames for each codon
    orf_id <- as.integer(
      levels(as.factor(hit_p@to))
    )
    spa0 <- spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_aa <- Biostrings::width(candi_all_orf[orf_id]) / 3
    spa_eff <- spa_all > eff_ribo_num
    return(list(
      orf_id = orf_id,
      spa_eff = spa_eff,
      orf_aa = orf_aa,
      spa_all = spa_all,
      spa_frame0 = spa_frame0,
      spa_frame1 = spa_frame1,
      spa_frame2 = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    ))
  })
  candi_all_orf = candi_all_orf[spa_frame_all[[1]]$orf_id]
  print(length(candi_all_orf))
  feature_mat_all <- get_feature(spa_lst = spa_frame_all[[1]])

  spa_frame_lst <- lapply(1:2, function(i) {
    hit_p <- IRanges::findOverlaps(rg_p[[i]], candi_orf)
    hit_p <- hit_p[
      hit_p@to %in%
        which(
          IRanges::countOverlaps(candi_orf, rg_p[[i]]) >= min_pos_num
        )
    ]
    dist_p <- rg_p[[i]]@start[hit_p@from] - candi_orf@start[hit_p@to]
    max_dist <- max(dist_p) %/% 3L + 3L
    idx_frame <- dist_p %% 3L
    idx_pos <- (dist_p %/% 3L) + 1L
    orf_n <- as.factor(hit_p@to)
    idx_orf <- as.integer(orf_n)

    p_num <- ribo_num_train[[i]]$vec_pnum

    # Build sparse codon-by-codon count matrices, one ORF per row
    spa_frame0 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 0L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 0L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 0L], 1)
    )
    spa_frame0 <- spa_frame0[-nrow(spa_frame0), ]
    spa_frame1 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 1L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 1L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 1L], 1)
    )
    spa_frame1 <- spa_frame1[-nrow(spa_frame1), ]
    spa_frame2 <- Matrix::sparseMatrix(
      i = c(idx_orf[idx_frame == 2L], length(levels(orf_n)) + 1L),
      j = c(idx_pos[idx_frame == 2L], max_dist),
      x = c(p_num[hit_p@from][idx_frame == 2L], 1)
    )
    spa_frame2 <- spa_frame2[-nrow(spa_frame2), ]

    # Normalize ribosome counts across frames for each codon
    orf_id <- as.integer(
      levels(as.factor(hit_p@to))
    )
    spa0 <- spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_aa <- Biostrings::width(candi_orf[orf_id]) / 3
    spa_eff <- spa_all > eff_ribo_num
    return(list(
      orf_id = orf_id,
      spa_eff = spa_eff,
      orf_aa = orf_aa,
      spa_all = spa_all,
      spa_frame0 = spa_frame0,
      spa_frame1 = spa_frame1,
      spa_frame2 = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    ))
  })
  # 首先判断阳性CDS中不是frame1占优势的部分
  spa_frame_p1 <- spa_frame_lst[[1]]
  spa_frame_neg <- spa_frame_lst[[2]]

  feature_mat_p1 <- get_feature(spa_lst = spa_frame_p1)
  feature_mat_neg <- get_feature(spa_lst = spa_frame_neg)

  change_pos_lable <- function(feature_mat_p1, feature_mat_neg) {
    X <- rbind(feature_mat_p1, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p1)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf, data = data.frame(feature_mat_p1))$predictions[, "0"]
    feature_mat_p2 <- feature_mat_p1[p_val < 0.5, ]
    X <- rbind(feature_mat_p2, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p2)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf_2 <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf_2, data = data.frame(feature_mat_p1))$predictions[, "0"]
    return(p_val)
  }

  #prob_cds_false <- change_pos_lable(feature_mat_p1 = feature_mat_p1, feature_mat_neg = feature_mat_neg)

  # 得到清洗过后的CDS frame阳性和阴性集
  #spa_frame_pos <- lapply(spa_frame_p1[c(1, 3)], function(x) x[prob_cds_false < 0.7])
  #spa_frame_pos <- c(spa_frame_pos, lapply(spa_frame_p1[c(-1, -3)], function(x) x[prob_cds_false < 0.7, ]))
  # 得出错位的阴性情况
  spa_frame_pos = spa_frame_p1
  spa_frame_ahead <- spa_frame_behind <- spa_frame_pos
  spa_frame_ahead$spa_frame0 <- spa_frame_pos$spa_frame1
  spa_frame_ahead$spa_frame1 <- spa_frame_pos$spa_frame2
  spa_frame_ahead$spa_frame2 <- spa_frame_pos$spa_frame0
  spa_frame_ahead$spa_frame0_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_ahead$spa_frame1_fq <- spa_frame_pos$spa_frame2_fq
  spa_frame_ahead$spa_frame2_fq <- spa_frame_pos$spa_frame0_fq

  spa_frame_behind$spa_frame2 <- spa_frame_pos$spa_frame1
  spa_frame_behind$spa_frame1 <- spa_frame_pos$spa_frame0
  spa_frame_behind$spa_frame0 <- spa_frame_pos$spa_frame2
  spa_frame_behind$spa_frame2_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_behind$spa_frame1_fq <- spa_frame_pos$spa_frame0_fq
  spa_frame_behind$spa_frame0_fq <- spa_frame_pos$spa_frame2_fq
  # 训练CDS翻译模型
  feat_m_pos <- get_feature(spa_lst = spa_frame_pos)
  feat_m_ahead <- get_feature(spa_lst = spa_frame_ahead)
  feat_m_behind <- get_feature(spa_lst = spa_frame_behind)
train_model_cds <- function(m_pos, m_neg,
                            threshold_method = c("maxF1", "youden"),
                            target_precision = NULL) {
  threshold_method <- match.arg(threshold_method)

  # -------- 组装训练数据 --------
  X <- rbind(m_pos, m_neg)
  y <- factor(c(rep(1L, nrow(m_pos)), rep(0L, nrow(m_neg))), levels = c(0L, 1L))
  y_int <- as.integer(as.character(y))  # 0/1 数值，用于计算指标

  train_df <- data.frame(y = y, X)

  # -------- 训练随机森林，开启 OOB 概率 --------
  rf <- ranger::ranger(
    formula   = y ~ .,
    data      = train_df,
    num.trees = 500,
    mtry      = 10,
    min.node.size   = 10,
    sample.fraction = 1,
    max.depth       = NULL,
    splitrule       = "gini",
    probability     = TRUE,  # 输出类别概率
    respect.unordered.factors = "partition",
    importance      = "impurity",
    num.threads     = max(1L, parallel::detectCores() - 1L),
    oob.error       = TRUE
  )

  # ranger 在 classification + probability=TRUE 时，
  # rf$predictions 是 n × 2 矩阵，每列是 P(y = level)
  p_oob <- rf$predictions[, "1"]  # 取阳性类的 OOB 概率

  # -------- 在 OOB 概率上做阈值扫描，选最优 cutoff --------

  .bin_metrics <- function(y_true, p, thr) {
    out <- lapply(thr, function(t) {
      pred <- as.integer(p >= t)
      tp <- sum(pred == 1 & y_true == 1)
      fp <- sum(pred == 1 & y_true == 0)
      tn <- sum(pred == 0 & y_true == 0)
      fn <- sum(pred == 0 & y_true == 1)

      prec <- if ((tp + fp) == 0) 1 else tp / (tp + fp)
      rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
      tpr  <- rec
      fpr  <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
      f1   <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)

      c(tp = tp, fp = fp, tn = tn, fn = fn,
        precision = prec, recall = rec,
        tpr = tpr, fpr = fpr, f1 = f1)
    })
    m <- do.call(rbind, out)
    data.frame(threshold = thr, m, row.names = NULL, check.names = FALSE)
  }

  # 候选阈值：用实际出现过的概率 + 一个细网格
  thr_cand <- sort(unique(c(p_oob, seq(0, 1, by = 0.001))), decreasing = TRUE)
  mtab <- .bin_metrics(y_int, p_oob, thr_cand)

  # -------- 选择最佳阈值 --------
  if (!is.null(target_precision)) {
    # 在 precision >= target_precision 的阈值中，选召回率尽量高的（阈值尽量小）
    ok <- which(mtab$precision >= target_precision)
    if (length(ok) > 0) {
      pick <- ok[length(ok)]    # thr_cand 是降序，这里取“最小的 t”
      t_opt <- mtab$threshold[pick]
      crit <- list(method = "precision>=target", target_precision = target_precision)
    } else {
      # 没有阈值能达到目标精度，就退回到 max F1
      pick <- which.max(mtab$f1)
      t_opt <- mtab$threshold[pick]
      crit <- list(method = "fallback_maxF1")
    }
  } else if (threshold_method == "youden") {
    j <- mtab$tpr - mtab$fpr
    pick <- which.max(j)
    t_opt <- mtab$threshold[pick]
    crit <- list(method = "youden", J = max(j, na.rm = TRUE))
  } else {
    # 默认：最大 F1
    pick <- which.max(mtab$f1)
    t_opt <- mtab$threshold[pick]
    crit <- list(method = "maxF1", F1 = mtab$f1[pick])
  }

  sel <- which(mtab$threshold == t_opt)[1]
  calib_summary <- list(
    threshold = t_opt,
    precision = mtab$precision[sel],
    recall    = mtab$recall[sel],
    f1        = mtab$f1[sel],
    tpr       = mtab$tpr[sel],
    fpr       = mtab$fpr[sel],
    counts    = unclass(mtab[sel, c("tp", "fp", "tn", "fn")]),
    criterion = crit
  )

  list(
    rf           = rf,          # 训练好的随机森林
    t_opt        = t_opt,       # 最佳阈值
    calib_table  = mtab,        # 全部阈值的指标（画 ROC/PR 都可以用）
    calib_summary = calib_summary
  )
}
  modle_cds <- train_model_cds(
    m_pos = feat_m_pos,
    m_neg = rbind(feature_mat_neg, feat_m_ahead, feat_m_behind),
    threshold_method = 'youden'
  )
  # 训练uORF dORF lncORF翻译模型
  pick_cds_spa <- function(spa_lst, aa_min = 6, aa_max = 100, bin_num = 10, rpf_min = 10, rpf_max = 200, eff_pos_min = 3) {
    orf_bin <- Matrix::summary(spa_lst$spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / spa_lst$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_lst$spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_lst$spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- spa_lst$orf_aa / bin_num
    spa_all_log <- spa_lst$spa_all
    spa_all_log@x <- log(spa_lst$spa_all@x + 1)
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(log(norf_aa), bin_num),
      log(1 + rpf_num),
      log(1 + as.vector(bin_f0_num(bin_lst, spa_lst$spa_frame0))),
      log(1 + as.vector(bin_f1_num(bin_lst, spa_lst$spa_frame1))),
      log(1 + as.vector(bin_f2_num(bin_lst, spa_lst$spa_frame2))),
      eff_pos,
      as.vector(bin_prop_sum(spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    # 过滤nORF标准
    shift_k <- function(v, k) {
      len_v <- length(v)
      v_tp <- v
      v_tp[1:(len_v / k)] <- v[(len_v * (k - 1) / k + 1):len_v]
      v_tp[(len_v / k + 1):len_v] <- v[1:(len_v * (k - 1) / k)]
      return(v_tp)
    }
    idx_norf <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= rpf_min) & (rpf_num <= rpf_max) & (eff_pos >= eff_pos_min)
    idx_cds <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= 30) & (eff_pos >= 10)
    idx_norf_change <- shift_k(idx_norf, bin_num)
    idx_cds_change <- shift_k(idx_cds, bin_num)
    # 转化坐标
    change_cord <- function(frame_i, aa_vec, bin_num) {
      orf_bin_f0 <- Matrix::summary(frame_i)
      bin_idx <- floor((orf_bin_f0$j - 0.5) / aa_vec[orf_bin_f0$i] * bin_num) + 1
      idx_tmp <- bin_idx < bin_num
      step_i <- floor(aa_vec[orf_bin_f0$i] / bin_num)
      orf_bin_f0$j[idx_tmp] <- orf_bin_f0$j[idx_tmp] + step_i[idx_tmp]
      orf_bin_f0$j[!idx_tmp] <- orf_bin_f0$j[!idx_tmp] + step_i[!idx_tmp] - aa_vec[orf_bin_f0$i][!idx_tmp] + 1L
      orf_bin_f0 <- orf_bin_f0[orf_bin_f0$j <= dim(frame_i)[2], ]
      spa_bin_f0 <- Matrix::sparseMatrix(
        i = orf_bin_f0$i,
        j = orf_bin_f0$j,
        x = orf_bin_f0$x,
        dims = dim(frame_i)
      )
      return(spa_bin_f0)
    }
    frame_0 <- change_cord(frame_i = spa_lst$spa_frame0, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_1 <- change_cord(frame_i = spa_lst$spa_frame1, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_2 <- change_cord(frame_i = spa_lst$spa_frame2, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    return(list(
      feature_mat = feature_mat[idx_norf, ],
      idx = list(
        idx_norf = idx_norf,
        idx_norf_change = idx_norf_change,
        idx_cds = idx_cds,
        idx_cds_change = idx_cds_change
      ),
      spa_frame = list(
        frame_0 = frame_0,
        frame_1 = frame_1,
        frame_2 = frame_2
      )
    ))
  }
  feat_norf_pos <- pick_cds_spa(spa_lst = spa_frame_pos)
  feat_norf_neg <- pick_cds_spa(spa_lst = spa_frame_neg)
  feat_norf_ahead <- pick_cds_spa(spa_lst = spa_frame_ahead)
  feat_norf_behind <- pick_cds_spa(spa_lst = spa_frame_behind)

  modle_ud_orf <- train_model_cds(
    m_pos = feat_norf_pos$feature_mat,
    m_neg = rbind(feat_norf_neg$feature_mat, feat_norf_ahead$feature_mat, feat_norf_behind$feature_mat)
  )

  # 训练uoORF doORF iORF翻译模型
  #   背景错框翻译，目标分4种情况
  # 阳性
  pick_overlap <- function(back_frame, head_info, bin_num = 10) {
    spa_frame0 <- back_frame$spa_frame0 + head_info$spa_frame$frame_0
    spa_frame1 <- back_frame$spa_frame1 + head_info$spa_frame$frame_1
    spa_frame2 <- back_frame$spa_frame2 + head_info$spa_frame$frame_2

    spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_bin <- Matrix::summary(spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / back_frame$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- back_frame$orf_aa / bin_num
    spa_all_log <- spa_all
    spa_all_log@x <- log(spa_all@x + 1)
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(log(norf_aa), bin_num),
      log(1 + rpf_num),
      log(1 + as.vector(bin_f0_num(bin_lst, spa_frame0))),
      log(1 + as.vector(bin_f1_num(bin_lst, spa_frame1))),
      log(1 + as.vector(bin_f2_num(bin_lst, spa_frame2))),
      eff_pos,
      as.vector(bin_prop_sum(spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    idx_cds <- (norf_aa >= 6) & (rpf_num >= 40) & (eff_pos >= 13) & head_info$idx$idx_norf
    return(feature_mat[idx_cds, ])
  }

  feat_oorf_p1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_pos)
  feat_oorf_p2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_pos)

  # 阴性
  feat_oorf_n1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_behind)
  feat_oorf_n2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_ahead)

  modle_overlap_orf <- train_model_cds(
    m_pos = rbind(feat_oorf_p1, feat_oorf_p2),
    m_neg = rbind(feat_oorf_n1, feat_oorf_n2)
  )
  return(list(
    train_model = list(
      modle_cds = modle_cds,
      modle_ud_orf = modle_ud_orf,
      modle_overlap_orf = modle_overlap_orf
    ),
    candidate_rg = candi_all_orf,
    ribo_num = ribo_num_train$ribo_num,
    feat_m = feature_mat_all
  ))
}


train_features_new <- function(ribo_num_train,
                           candi_translat,
                           tx_info,
                           min_pos_num,
                           eff_ribo_num) {
                            #browser()
                            matrix_sd = function(M){
  M2 <- M
    M2@x <- M2@x^2
    row_n = Matrix::rowSums(M !=0)
  row_mean  <- Matrix::rowSums(M)/row_n
  row_mean2 <- Matrix::rowSums(M2)/row_n
row_var_all <- pmax(row_mean2 - row_mean^2, 0)  # 数值误差纠正
row_sd  <- sqrt(row_var_all)
return(list(row_n=row_n, row_mean = row_mean, row_sd=row_sd))
}
  ## ===================== Helper1: 只构建基础 spa 矩阵 =====================

build_spa_base <- function(rg_p, candi_rg, p_num,
                           eff_ribo_num, min_pos_num) {
                            #browser()
  hits <- IRanges::findOverlaps(rg_p, candi_rg)
  if (length(hits) == 0L) return(NULL)

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  # 每个 ORF 命中的 P 数；只保留 read 数 >= min_pos_num 的 ORF
  orf_counts <- tabulate(sh, nbins = length(candi_rg))
  keep_orf   <- which(orf_counts >= min_pos_num)
  hits       <- hits[sh %in% keep_orf]
  if (length(hits) == 0L) return(NULL)

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  # ORF index → 行号 1..n_orf
  orf_levels <- sort(unique(sh))
  orf_factor <- factor(sh, levels = orf_levels)
  idx_orf    <- as.integer(orf_factor)
  n_orf      <- length(orf_levels)

  # 与 ORF 起点的距离（nt）
  dist_p   <- IRanges::start(rg_p)[qh] - IRanges::start(candi_rg)[sh]
  idx_pos  <- (dist_p %/% 3L) + 1L
  idx_frame <- dist_p %% 3L
  max_pos  <- max(idx_pos)

  w <- p_num[qh]  # 每个 P 的权重

  make_spa_frame <- function(which_frame) {
    sel <- idx_frame == which_frame
    Matrix::sparseMatrix(
      i    = c(idx_orf[sel], n_orf + 1L),
      j    = c(idx_pos[sel], max_pos),
      x    = c(w[sel], 1),
      dims = c(n_orf + 1L, max_pos)
    )[1:n_orf, , drop = FALSE]
  }

  spa_frame0 <- make_spa_frame(0L)
  spa_frame1 <- make_spa_frame(1L)
  spa_frame2 <- make_spa_frame(2L)

  spa_all <- spa_frame0 + spa_frame1 + spa_frame2
  spa_eff <- spa_all > eff_ribo_num   # 有效 codon

  # 对齐模式：三帧都补零到和 spa_all 的非零位一致
  spa0       <- spa_all
  spa0@x     <- rep(0, length(spa0@x))
  spa_frame0 <- spa_frame0 + spa0
  spa_frame1 <- spa_frame1 + spa0
  spa_frame2 <- spa_frame2 + spa0

  eps <- 1e-12
  frac_from <- function(num_mat, den_mat) {
    out    <- den_mat
    out@x  <- (num_mat@x + eps) / (den_mat@x + 3 * eps)
    out
  }
  spa_frame0_fq <- frac_from(spa_frame0, spa_all)
  spa_frame1_fq <- frac_from(spa_frame1, spa_all)
  spa_frame2_fq <- frac_from(spa_frame2, spa_all)

  orf_aa <- Biostrings::width(candi_rg[orf_levels]) / 3L

  list(
    orf_id        = orf_levels,
    orf_aa        = orf_aa,
    spa_all       = spa_all,
    spa_eff       = spa_eff,
    spa_frame0    = spa_frame0,
    spa_frame1    = spa_frame1,
    spa_frame2    = spa_frame2,
    spa_frame0_fq = spa_frame0_fq,
    spa_frame1_fq = spa_frame1_fq,
    spa_frame2_fq = spa_frame2_fq
    # 这里不再放 logodds/entropy/dist；这些交给 Helper2
  )
}

## ===================== Helper2: 在基础 spa 上加帧特征 =====================

add_frame_features <- function(spa_lst, cds_frame = NULL) {
  eps <- 1e-12

  # ---- log-odds ----
  logit_from <- function(fq_mat) {
    out <- fq_mat
    p   <- fq_mat@x
    out@x <- log((p + eps) / (1 - p + eps))
    out
  }
  bias_logodds_0 <- logit_from(spa_lst$spa_frame0_fq)
  bias_logodds_1 <- logit_from(spa_lst$spa_frame1_fq)
  bias_logodds_2 <- logit_from(spa_lst$spa_frame2_fq)

  # ---- 熵 ----
  entropy_from <- function(fq_mat) {
    out <- fq_mat
    p   <- fq_mat@x
    out@x <- -(p + eps) * log(p + eps)
    out
  }
  entropy_0 <- entropy_from(spa_lst$spa_frame0_fq)
  entropy_1 <- entropy_from(spa_lst$spa_frame1_fq)
  entropy_2 <- entropy_from(spa_lst$spa_frame2_fq)

  # ---- 均匀分布 (1/3,1/3,1/3) 的 L2 距离：dist_uniform_1 ----
  dist_uniform_1 <- spa_lst$spa_all
  dist_uniform_1@x <- sqrt(
    (spa_lst$spa_frame0_fq@x - 1/3)^2 +
    (spa_lst$spa_frame1_fq@x - 1/3)^2 +
    (spa_lst$spa_frame2_fq@x - 1/3)^2
  )

  # ---- 如提供 cds_frame，则再算 dist_uniform_2~8 ----
  if (!is.null(cds_frame)) {
    cds_frame       <- cds_frame / sum(cds_frame)
    cds_frame_ahead <- cds_frame[c(2, 3, 1)]
    cds_frame_behind <- cds_frame[c(3, 1, 2)]

    cf1 <- cds_frame + cds_frame_ahead
    cds_frame1 <- cf1 / sum(cf1)

    cf2 <- cds_frame + cds_frame_behind
    cds_frame2 <- cf2 / sum(cf2)

    cf3 <- cds_frame * (1/3) + cds_frame_ahead * (2/3)
    cds_frame3 <- cf3 / sum(cf3)

    cf4 <- cds_frame * (1/3) + cds_frame_behind * (2/3)
    cds_frame4 <- cf4 / sum(cf4)

    mk_dist <- function(target_vec) {
      m <- spa_lst$spa_all
      m@x <- sqrt(
        (spa_lst$spa_frame0_fq@x - target_vec[1])^2 +
        (spa_lst$spa_frame1_fq@x - target_vec[2])^2 +
        (spa_lst$spa_frame2_fq@x - target_vec[3])^2
      )
      m
    }

    dist_uniform_2 <- mk_dist(cds_frame_ahead)
    dist_uniform_3 <- mk_dist(cds_frame_behind)
    dist_uniform_4 <- mk_dist(cds_frame)
    dist_uniform_5 <- mk_dist(cds_frame1)
    dist_uniform_6 <- mk_dist(cds_frame2)
    dist_uniform_7 <- mk_dist(cds_frame3)
    dist_uniform_8 <- mk_dist(cds_frame4)
  } else {
    # 如果没给 cds_frame，就只保留 dist_uniform_1，其余设为 NULL
    dist_uniform_2 <- dist_uniform_3 <- dist_uniform_4 <- NULL
    dist_uniform_5 <- dist_uniform_6 <- dist_uniform_7 <- dist_uniform_8 <- NULL
  }

  # ---- 把这些特征塞回 list ----
  spa_lst$bias_logodds_0  <- bias_logodds_0
  spa_lst$bias_logodds_1  <- bias_logodds_1
  spa_lst$bias_logodds_2  <- bias_logodds_2
  spa_lst$entropy_0       <- entropy_0
  spa_lst$entropy_1       <- entropy_1
  spa_lst$entropy_2       <- entropy_2
  spa_lst$dist_uniform_1  <- dist_uniform_1
  spa_lst$dist_uniform_2  <- dist_uniform_2
  spa_lst$dist_uniform_3  <- dist_uniform_3
  spa_lst$dist_uniform_4  <- dist_uniform_4
  spa_lst$dist_uniform_5  <- dist_uniform_5
  spa_lst$dist_uniform_6  <- dist_uniform_6
  spa_lst$dist_uniform_7  <- dist_uniform_7
  spa_lst$dist_uniform_8  <- dist_uniform_8

  spa_lst
}

## ===================== Helper 3: 从 spa 列表抽 feature 矩阵 =================

  get_feature <- function(spa_lst) {
    bin_eff <- function(spa_eff) {
      Matrix::rowSums(spa_eff)
    }
    bin_prop_sum <- function(spa_prop, spa_eff) {
      Matrix::rowSums(spa_prop * spa_eff)
    }
    bin_weight_sum <- function(spa_ribo, spa_prop, spa_eff) {
      Matrix::rowSums(spa_eff * (spa_ribo * spa_prop))
    }

    spa_all_log <- spa_lst$spa_all

    # 每行ORF生成位置线性权重分别递增和递减
row_linear_fill <- function(M, from = 0, to = 1) {
  stopifnot(inherits(M, "dgCMatrix"))
  # 提取稀疏矩阵的(i, j, x) 三元组
  S <- Matrix::summary(M)  # columns: i, j, x
  
  # 如果全是 0，直接返回同维度的零矩阵
  if (nrow(S) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = dim(M),
      dimnames = dimnames(M)
    ))
  }
  
  # 按行索引 i 分组
  idx_by_row <- split(seq_len(nrow(S)), S$i)
  
  newx <- numeric(nrow(S))  # 存放新的非零值
  
  for (ind in idx_by_row) {
    k <- length(ind)
    if (k == 0L) next
    
    # 按列 j 排序，保证沿着一行从左到右递增
    ord <- order(S$j[ind])
    idx <- ind[ord]
    
    # 在线性区间 [from, to] 上生成 k 个点
    newx[idx] <- seq(from, to, length.out = k)
    # 如果你希望只有一个点时是 1 而不是 0，可以改成：
    # if (k == 1L) newx[idx] <- to else newx[idx] <- seq(from, to, length.out = k)
  }
  
  Matrix::sparseMatrix(
    i = S$i,
    j = S$j,
    x = newx,
    dims = dim(M),
    dimnames = dimnames(M)
  )
}
spa_all_log1 = row_linear_fill(spa_all_log, from = 0, to = 1)
spa_all_log0 = row_linear_fill(spa_all_log, from = 1, to = 0)

    # 每个codon最大frame统计
    max_frame = function(f1, f2, f3){
      return(data.frame(list(
        max1 = Matrix::rowSums((f1<f2)&(f2<f3)),
        max2 = Matrix::rowSums((f2<f1)&(f1<f3)),
        max3 = Matrix::rowSums((f2<f3)&(f3<f1)),
        max4 = Matrix::rowSums((f3<f2)&(f2<f1)),
        max5 = Matrix::rowSums((f1<f3)&(f3<f2)),
        max6 = Matrix::rowSums((f3<f1)&(f1<f2))
      ))
      )
    }

    feature_mat <- cbind(
      orf_aa       = spa_lst$orf_aa,
      rpf_all      = Matrix::rowSums(spa_lst$spa_all),
      rpf_f0       = Matrix::rowSums(spa_lst$spa_frame0),
      rpf_f1       = Matrix::rowSums(spa_lst$spa_frame1),
      rpf_f2       = Matrix::rowSums(spa_lst$spa_frame2),
      eff_codons   = bin_eff(spa_lst$spa_eff),

      f0_prop_sum  = bin_prop_sum(spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_prop_sum  = bin_prop_sum(spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_prop_sum  = bin_prop_sum(spa_lst$spa_frame2_fq, spa_lst$spa_eff),

      f0_weighted  = bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_weighted  = bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_weighted  = bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f00_weighted  = bin_weight_sum(spa_all_log1, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f10_weighted  = bin_weight_sum(spa_all_log1, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f20_weighted  = bin_weight_sum(spa_all_log1, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f10_weighted  = bin_weight_sum(spa_all_log0, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f11_weighted  = bin_weight_sum(spa_all_log0, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f12_weighted  = bin_weight_sum(spa_all_log0, spa_lst$spa_frame2_fq, spa_lst$spa_eff),

      logodds0_sum = bin_prop_sum(spa_lst$bias_logodds_0, spa_lst$spa_eff),
      logodds1_sum = bin_prop_sum(spa_lst$bias_logodds_1, spa_lst$spa_eff),
      logodds2_sum = bin_prop_sum(spa_lst$bias_logodds_2, spa_lst$spa_eff),

      entropy0_sum = bin_prop_sum(spa_lst$entropy_0, spa_lst$spa_eff),
      entropy1_sum = bin_prop_sum(spa_lst$entropy_1, spa_lst$spa_eff),
      entropy2_sum = bin_prop_sum(spa_lst$entropy_2, spa_lst$spa_eff),

      dist_u1_sum  = bin_prop_sum(spa_lst$dist_uniform_1, spa_lst$spa_eff),
      dist_u2_sum  = bin_prop_sum(spa_lst$dist_uniform_2, spa_lst$spa_eff),
      dist_u3_sum  = bin_prop_sum(spa_lst$dist_uniform_3, spa_lst$spa_eff),
      dist_u4_sum  = bin_prop_sum(spa_lst$dist_uniform_4, spa_lst$spa_eff),
      dist_u5_sum  = bin_prop_sum(spa_lst$dist_uniform_5, spa_lst$spa_eff),
      dist_u6_sum  = bin_prop_sum(spa_lst$dist_uniform_6, spa_lst$spa_eff),
      dist_u7_sum  = bin_prop_sum(spa_lst$dist_uniform_7, spa_lst$spa_eff),
      dist_u8_sum  = bin_prop_sum(spa_lst$dist_uniform_8, spa_lst$spa_eff)
    )
feature_mat <- cbind(feature_mat, max_frame(spa_lst$spa_frame0, spa_lst$spa_frame1, spa_lst$spa_frame2))
    colnames(feature_mat) <- paste0("feat_", seq_len(ncol(feature_mat)))
    feature_mat
  }
  ## ===================== 1. P-site & ORF 定义 ================================

  rg_p_pos <- IRanges::IRanges(
    start = ribo_num_train$ribo_num$candi_pos,
    width = 1L
  )
  rg_p_neg <- IRanges::IRanges(
    start = ribo_num_train$ribo_num_neg$candi_pos,
    width = 1L
  )

  # CDS 区域
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end   = tx_info$mix_tx_pos$utr3_p5 - 1L
  )

  # 所有候选翻译区域
  candi_all_orf <- do.call(c, candi_translat)
  candi_all_orf@start <- candi_all_orf@start + 3L
  candi_all_orf@width <- candi_all_orf@width - 1L
  candi_all_orf_back  <- candi_all_orf

  # CDS stop-stop 对应的 ORF（主 CDS）
  cds_idx   <- IRanges::end(candi_all_orf) %in% IRanges::end(rg_cds)
  candi_orf <- candi_all_orf[cds_idx]
  
  ## ===================== 2. 全 ORF（用于最后打分） + cds_frame dist =========

  spa_all_pos <- build_spa_base(
    rg_p        = rg_p_pos,
    candi_rg    = candi_all_orf,
    p_num       = ribo_num_train$ribo_num$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  if (is.null(spa_all_pos)) {
    stop("No ORF passed min_pos_num in positive set (all ORFs).")
  }

# 基于 CDS ORF 的 frame 占比，算全局 cds_frame（和你原来一致）
  cds_rows <- sort(match(candi_orf, candi_all_orf_back[spa_all_pos$orf_id]))
  cds_rows <- cds_rows[!is.na(cds_rows)]

  if (length(cds_rows) > 0L) {
    cds_frame_all <- c(
      sum(spa_all_pos$spa_frame0[cds_rows, ]),
      sum(spa_all_pos$spa_frame1[cds_rows, ]),
      sum(spa_all_pos$spa_frame2[cds_rows, ])
    )
    cds_frame_all <- cds_frame_all / sum(spa_all_pos$spa_all[cds_rows, ])
  } else {
    cds_frame_all <- c(
      sum(spa_all_pos$spa_frame0),
      sum(spa_all_pos$spa_frame1),
      sum(spa_all_pos$spa_frame2)
    )
    cds_frame_all <- cds_frame_all / sum(spa_all_pos$spa_all)
  }
    
spa_all_pos <- add_frame_features(spa_all_pos, cds_frame_all)

  feature_mat_all <- get_feature(spa_all_pos)
  candi_all_orf   <- candi_all_orf_back[spa_all_pos$orf_id]

## ===================== 3. CDS 正/负集 + ahead/behind =======================

  spa_frame_pos <- build_spa_base(
    rg_p        = rg_p_pos,
    candi_rg    = candi_orf,
    p_num       = ribo_num_train$ribo_num$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  spa_frame_neg <- build_spa_base(
    rg_p        = rg_p_neg,
    candi_rg    = candi_orf,
    p_num       = ribo_num_train$ribo_num_neg$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  if (is.null(spa_frame_pos) || is.null(spa_frame_neg)) {
    stop("Not enough positive/negative CDS ORFs for training.")
  }

  # 用正集的全局 frame 占比作为 cds_frame（比原来按 set 各算更直观）
  cds_frame_cds <- c(
    sum(spa_frame_pos$spa_frame0),
    sum(spa_frame_pos$spa_frame1),
    sum(spa_frame_pos$spa_frame2)
  )
  cds_frame_cds <- cds_frame_cds / sum(spa_frame_pos$spa_all)

  spa_frame_pos <- add_frame_features(spa_frame_pos, cds_frame_cds)
  spa_frame_neg <- add_frame_features(spa_frame_neg, cds_frame_cds)

  # 错框版本（ahead / behind）
  spa_frame_ahead  <- spa_frame_pos
  spa_frame_behind <- spa_frame_pos

  spa_frame_ahead$spa_frame0    <- spa_frame_pos$spa_frame1
  spa_frame_ahead$spa_frame1    <- spa_frame_pos$spa_frame2
  spa_frame_ahead$spa_frame2    <- spa_frame_pos$spa_frame0
  spa_frame_ahead$spa_frame0_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_ahead$spa_frame1_fq <- spa_frame_pos$spa_frame2_fq
  spa_frame_ahead$spa_frame2_fq <- spa_frame_pos$spa_frame0_fq

  spa_frame_behind$spa_frame2    <- spa_frame_pos$spa_frame1
  spa_frame_behind$spa_frame1    <- spa_frame_pos$spa_frame0
  spa_frame_behind$spa_frame0    <- spa_frame_pos$spa_frame2
  spa_frame_behind$spa_frame2_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_behind$spa_frame1_fq <- spa_frame_pos$spa_frame0_fq
  spa_frame_behind$spa_frame0_fq <- spa_frame_pos$spa_frame2_fq

  # ahead/behind 也用同一个 cds_frame_cds 生成 dist_uniform_2~8
  spa_frame_ahead  <- add_frame_features(spa_frame_ahead,  cds_frame_cds)
  spa_frame_behind <- add_frame_features(spa_frame_behind, cds_frame_cds)

  # 特征矩阵（CDS 级）
  feat_m_pos    <- get_feature(spa_frame_pos)
  feat_m_ahead  <- get_feature(spa_frame_ahead)
  feat_m_behind <- get_feature(spa_frame_behind)
  feature_mat_neg <- get_feature(spa_frame_neg)

  ## ===================== 4. 主 CDS 模型 ======================================
  modle_cds = train_model_cds_rf(
    m_pos = feat_m_pos,
    m_neg = rbind(feature_mat_neg, feat_m_ahead, feat_m_behind)
  )
## ===================== 5. pick_cds_spa：nORF 分 bin 特征 ===================
feat_norf_pos    <- make_bin_features(spa_lst = spa_frame_pos)
  feat_norf_neg    <- make_bin_features(spa_lst = spa_frame_neg)
  feat_norf_ahead  <- make_bin_features(spa_lst = spa_frame_ahead)
  feat_norf_behind <- make_bin_features(spa_lst = spa_frame_behind)
make_shifted_frames()












  change_pos_lable <- function(feature_mat_p1, feature_mat_neg) {
    X <- rbind(feature_mat_p1, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p1)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf, data = data.frame(feature_mat_p1))$predictions[, "0"]
    feature_mat_p2 <- feature_mat_p1[p_val < 0.5, ]
    X <- rbind(feature_mat_p2, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p2)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf_2 <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 200,
      mtry = 3,
      min.node.size = 10,
      sample.fraction = 0.8,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf_2, data = data.frame(feature_mat_p1))$predictions[, "0"]
    return(p_val)
  }

  prob_cds_false <- change_pos_lable(feature_mat_p1 = feature_mat_p1, feature_mat_neg = feature_mat_neg)

  # 得到清洗过后的CDS frame阳性和阴性集
  spa_frame_pos <- lapply(spa_frame_p1[c(1, 3)], function(x) x[prob_cds_false < 0.7])
  spa_frame_pos <- c(spa_frame_pos, lapply(spa_frame_p1[c(-1, -3)], function(x) x[prob_cds_false < 0.7, ]))
  # 得出错位的阴性情况
  spa_frame_pos = spa_frame_p1
  spa_frame_ahead <- spa_frame_behind <- spa_frame_pos
  spa_frame_ahead$spa_frame0 <- spa_frame_pos$spa_frame1
  spa_frame_ahead$spa_frame1 <- spa_frame_pos$spa_frame2
  spa_frame_ahead$spa_frame2 <- spa_frame_pos$spa_frame0
  spa_frame_ahead$spa_frame0_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_ahead$spa_frame1_fq <- spa_frame_pos$spa_frame2_fq
  spa_frame_ahead$spa_frame2_fq <- spa_frame_pos$spa_frame0_fq

  spa_frame_behind$spa_frame2 <- spa_frame_pos$spa_frame1
  spa_frame_behind$spa_frame1 <- spa_frame_pos$spa_frame0
  spa_frame_behind$spa_frame0 <- spa_frame_pos$spa_frame2
  spa_frame_behind$spa_frame2_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_behind$spa_frame1_fq <- spa_frame_pos$spa_frame0_fq
  spa_frame_behind$spa_frame0_fq <- spa_frame_pos$spa_frame2_fq
  # 训练CDS翻译模型
  feat_m_pos <- get_feature(spa_lst = spa_frame_pos)
  feat_m_ahead <- get_feature(spa_lst = spa_frame_ahead)
  feat_m_behind <- get_feature(spa_lst = spa_frame_behind)

  modle_cds <- train_model_cds(
    m_pos = feat_m_pos,
    m_neg = rbind(feature_mat_neg, feat_m_ahead, feat_m_behind),
  )
library(randomForest)

train_model_cds_rf <- function(m_pos, m_neg,
                               ntree  = 500,
                               mtry   = NULL,
                               classwt = NULL) {
  # -------- 组装训练数据 --------
  X <- rbind(m_pos, m_neg)
  y <- factor(
    c(rep(1L, nrow(m_pos)),  # 正类：翻译 CDS
      rep(0L, nrow(m_neg))), # 负类：非翻译
    levels = c(0L, 1L)
  )

  train_df <- data.frame(y = y, X)

  # -------- 默认 mtry：sqrt(p) --------
  if (is.null(mtry)) {
    mtry <- max(1L, floor(sqrt(ncol(X))))
  }

  # -------- 训练 randomForest --------
  rf <- randomForest::randomForest(
    y ~ .,
    data    = train_df,
    ntree   = 500,
    mtry    = mtry,
    importance = TRUE,    # 方便后面看特征重要性
     # 需要的话可以设置正负类的采样数
    classwt  = c("0"=3, "1"=1),   # 不平衡时可设，例如 c("0"=1, "1"=3)
    keep.forest = TRUE
  )
  list(
    model   = rf,
    n_pos   = nrow(m_pos),
    n_neg   = nrow(m_neg)
  )
}
plot(rf)
rf$importance
predict_ud_orf <- function(model_obj, X_new) {
  rf <- model_obj$model
  new_df <- data.frame(X_new)

  # type = "prob" 返回两列：P(class=0), P(class=1)
  p <- predict(rf, newdata = new_df, type = "prob")
  as.numeric(p[, "1"])  # 取预测为“翻译 ORF”(1) 的概率
}


  # 训练uORF dORF lncORF翻译模型
  pick_cds_spa <- function(spa_lst, aa_min = 6, aa_max = 1e6, bin_num = 10, rpf_min = 1, rpf_max = 1e6, eff_pos_min = 3) {
    orf_bin <- Matrix::summary(spa_lst$spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / spa_lst$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_lst$spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_lst$spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- spa_lst$orf_aa / bin_num
    spa_all_log <- spa_lst$spa_all
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(norf_aa, bin_num),
      rpf_num,
      as.vector(bin_f0_num(bin_lst, spa_lst$spa_frame0)),
      as.vector(bin_f1_num(bin_lst, spa_lst$spa_frame1)),
      as.vector(bin_f2_num(bin_lst, spa_lst$spa_frame2)),
      eff_pos,
      as.vector(bin_prop_sum(spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_lst$spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    # 过滤nORF标准
    shift_k <- function(v, k) {
      len_v <- length(v)
      v_tp <- v
      v_tp[1:(len_v / k)] <- v[(len_v * (k - 1) / k + 1):len_v]
      v_tp[(len_v / k + 1):len_v] <- v[1:(len_v * (k - 1) / k)]
      return(v_tp)
    }
    idx_norf <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= rpf_min) & (rpf_num <= rpf_max) & (eff_pos >= eff_pos_min)
    idx_cds <- (norf_aa >= aa_min) & (norf_aa <= aa_max) & (rpf_num >= 30) & (eff_pos >= 10)
    idx_norf_change <- shift_k(idx_norf, bin_num)
    idx_cds_change <- shift_k(idx_cds, bin_num)
    # 转化坐标
    change_cord <- function(frame_i, aa_vec, bin_num) {
      orf_bin_f0 <- Matrix::summary(frame_i)
      bin_idx <- floor((orf_bin_f0$j - 0.5) / aa_vec[orf_bin_f0$i] * bin_num) + 1
      idx_tmp <- bin_idx < bin_num
      step_i <- floor(aa_vec[orf_bin_f0$i] / bin_num)
      orf_bin_f0$j[idx_tmp] <- orf_bin_f0$j[idx_tmp] + step_i[idx_tmp]
      orf_bin_f0$j[!idx_tmp] <- orf_bin_f0$j[!idx_tmp] + step_i[!idx_tmp] - aa_vec[orf_bin_f0$i][!idx_tmp] + 1L
      orf_bin_f0 <- orf_bin_f0[orf_bin_f0$j <= dim(frame_i)[2], ]
      spa_bin_f0 <- Matrix::sparseMatrix(
        i = orf_bin_f0$i,
        j = orf_bin_f0$j,
        x = orf_bin_f0$x,
        dims = dim(frame_i)
      )
      return(spa_bin_f0)
    }
    frame_0 <- change_cord(frame_i = spa_lst$spa_frame0, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_1 <- change_cord(frame_i = spa_lst$spa_frame1, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    frame_2 <- change_cord(frame_i = spa_lst$spa_frame2, aa_vec = spa_lst$orf_aa, bin_num = bin_num)
    return(list(
      feature_mat = feature_mat[idx_norf, ],
      idx = list(
        idx_norf = idx_norf,
        idx_norf_change = idx_norf_change,
        idx_cds = idx_cds,
        idx_cds_change = idx_cds_change
      ),
      spa_frame = list(
        frame_0 = frame_0,
        frame_1 = frame_1,
        frame_2 = frame_2
      )
    ))
  }
  make_shifted_frames <- function(spa_lst, bin_num = 10) {

  change_cord <- function(frame_i, aa_vec, bin_num) {
    orf_bin <- Matrix::summary(frame_i)
    if (nrow(orf_bin) == 0L) {
      return(Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        x = numeric(0), dims = dim(frame_i)
      ))
    }

    bin_idx <- floor((orf_bin$j - 0.5) / aa_vec[orf_bin$i] * bin_num) + 1
    idx_tmp <- bin_idx < bin_num

    step_i <- floor(aa_vec[orf_bin$i] / bin_num)

    # 前 bin 的 codon 往后平移 step_i
    orf_bin$j[idx_tmp] <- orf_bin$j[idx_tmp] + step_i[idx_tmp]
    # 最后一个 bin 的 codon 往后平移，并保证不超过 ORF 终点
    orf_bin$j[!idx_tmp] <- orf_bin$j[!idx_tmp] +
      step_i[!idx_tmp] - aa_vec[orf_bin$i][!idx_tmp] + 1L

    orf_bin <- orf_bin[orf_bin$j <= dim(frame_i)[2], ]

    Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(frame_i)
    )
  }

  frame_0 <- change_cord(spa_lst$spa_frame0, spa_lst$orf_aa, bin_num)
  frame_1 <- change_cord(spa_lst$spa_frame1, spa_lst$orf_aa, bin_num)
  frame_2 <- change_cord(spa_lst$spa_frame2, spa_lst$orf_aa, bin_num)

  list(
    orf_aa    = spa_lst$orf_aa,
    spa_frame = list(
      frame_0 = frame_0,
      frame_1 = frame_1,
      frame_2 = frame_2
    ),
    bin_num   = bin_num
  )
}
make_bin_features <- function(spa_lst,
                              aa_min      = 6,
                              aa_max      = 1e6,
                              bin_num     = 10,
                              rpf_min     = 1,
                              rpf_max     = 1e6,
                              eff_pos_min = 3) {

  # ---- 1) 把每个 ORF 均分成 bin_num 段，给每个 codon 标 bin 号 ----
  orf_bin <- Matrix::summary(spa_lst$spa_all)       # i=orf, j=codon, x=count
  orf_bin$x <- floor((orf_bin$j - 0.5) / spa_lst$orf_aa[orf_bin$i] * bin_num) + 1

  spa_bin <- Matrix::sparseMatrix(
    i    = orf_bin$i,
    j    = orf_bin$j,
    x    = orf_bin$x,                  # 填 bin index
    dims = dim(spa_lst$spa_all)
  )

  # 每个 bin 一个 mask（同一个 ORF 内的不同区间）
  bin_lst <- lapply(seq_len(bin_num), function(b) spa_bin == b)

  # ---- 2) 按 bin 汇总各种量 ----
  bin_rpf_num <- function(bin_lst, spa_all) {
    sapply(bin_lst, function(mask) Matrix::rowSums(spa_all * mask))
  }

  bin_f_num <- function(bin_lst, spa_frame) {
    sapply(bin_lst, function(mask) Matrix::rowSums(spa_frame * mask))
  }

  bin_eff <- function(bin_lst) {
    sapply(bin_lst, function(mask) Matrix::rowSums(mask))
  }

  bin_prop_sum <- function(spa_prop, bin_lst) {
    sapply(bin_lst, function(mask) Matrix::rowSums(spa_prop * mask))
  }

  bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
    sapply(bin_lst, function(mask) Matrix::rowSums(mask * (spa_ribo * spa_prop)))
  }

  norf_aa    <- spa_lst$orf_aa / bin_num              # 每个 bin 对应的 “局部 aa 长度”
  spa_all    <- spa_lst$spa_all
  spa_all_log <- spa_all   # 如果要 log 可再改

  # 基本 bin 特征
  rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_all))
  eff_pos <- as.vector(bin_eff(bin_lst))

  f0_num <- as.vector(bin_f_num(bin_lst, spa_lst$spa_frame0))
  f1_num <- as.vector(bin_f_num(bin_lst, spa_lst$spa_frame1))
  f2_num <- as.vector(bin_f_num(bin_lst, spa_lst$spa_frame2))

  f0_fq_sum <- as.vector(bin_prop_sum(spa_lst$spa_frame0_fq, bin_lst))
  f1_fq_sum <- as.vector(bin_prop_sum(spa_lst$spa_frame1_fq, bin_lst))
  f2_fq_sum <- as.vector(bin_prop_sum(spa_lst$spa_frame2_fq, bin_lst))

  f0_wsum <- as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame0_fq, bin_lst))
  f1_wsum <- as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame1_fq, bin_lst))
  f2_wsum <- as.vector(bin_weight_sum(spa_all_log, spa_lst$spa_frame2_fq, bin_lst))

  # 与 get_feature 对齐的 “二阶特征”：log-odds, entropy, dist_uniform_1~8
  # log-odds
  logodds0_sum <- as.vector(bin_prop_sum(spa_lst$bias_logodds_0, bin_lst))
  logodds1_sum <- as.vector(bin_prop_sum(spa_lst$bias_logodds_1, bin_lst))
  logodds2_sum <- as.vector(bin_prop_sum(spa_lst$bias_logodds_2, bin_lst))

  # entropy
  ent0_sum <- as.vector(bin_prop_sum(spa_lst$entropy_0, bin_lst))
  ent1_sum <- as.vector(bin_prop_sum(spa_lst$entropy_1, bin_lst))
  ent2_sum <- as.vector(bin_prop_sum(spa_lst$entropy_2, bin_lst))

  # dist_uniform 系列
  du1_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_1, bin_lst))
  du2_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_2, bin_lst))
  du3_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_3, bin_lst))
  du4_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_4, bin_lst))
  du5_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_5, bin_lst))
  du6_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_6, bin_lst))
  du7_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_7, bin_lst))
  du8_sum <- as.vector(bin_prop_sum(spa_lst$dist_uniform_8, bin_lst))

  # 每个 ORF 在每个 bin 的 aa 长度（复制 bin_num 次）
  feat_len <- rep(norf_aa, each = bin_num)

  feature_mat <- cbind(
    feat_len,          # 1
    rpf_num,           # 2
    f0_num, f1_num, f2_num,         # 3-5
    eff_pos,           # 6
    f0_fq_sum, f1_fq_sum, f2_fq_sum,     # 7-9
    f0_wsum,  f1_wsum,  f2_wsum,         # 10-12
    logodds0_sum, logodds1_sum, logodds2_sum,  # 13-15
    ent0_sum, ent1_sum, ent2_sum,            # 16-18
    du1_sum, du2_sum, du3_sum, du4_sum,      # 19-22
    du5_sum, du6_sum, du7_sum, du8_sum       # 23-26
  )

  colnames(feature_mat) <- paste0("feat_", seq_len(ncol(feature_mat)))

  # ---- 3) nORF / CDS 的选择 (按 ORF 维度) ----
  n_orf <- length(spa_lst$orf_aa)

  # 这些条件是按 ORF 级别定义的，因此只在 ORF 长度上定义一次，然后在 bin 维度上重复
  orf_len <- spa_lst$orf_aa
  # ORF 总 RPF、总有效位点，用于 ORF 级别过滤
  rpf_tot <- Matrix::rowSums(spa_all)
  eff_tot <- Matrix::rowSums(spa_all > 0)

  idx_norf_orf <- (orf_len >= aa_min) & (orf_len <= aa_max) &
                  (rpf_tot >= rpf_min) & (rpf_tot <= rpf_max) &
                  (eff_tot >= eff_pos_min)

  idx_cds_orf  <- (orf_len >= aa_min) & (orf_len <= aa_max) &
                  (rpf_tot >= 30) & (eff_tot >= 10)

  # 复制到 bin 级别（每个 ORF 的所有 bin 同进同出）
  idx_norf <- rep(idx_norf_orf, each = bin_num)
  idx_cds  <- rep(idx_cds_orf,  each = bin_num)

  list(
    feature_mat = feature_mat[idx_norf, , drop = FALSE],
    idx = list(
      idx_norf = idx_norf_orf,  # ORF 级
      idx_cds  = idx_cds_orf
    ),
    meta = list(
      orf_aa  = spa_lst$orf_aa,
      bin_num = bin_num
    )
  )
}

  feat_norf_pos <- pick_cds_spa(spa_lst = spa_frame_pos)
  feat_norf_neg <- pick_cds_spa(spa_lst = spa_frame_neg)
  feat_norf_ahead <- pick_cds_spa(spa_lst = spa_frame_ahead)
  feat_norf_behind <- pick_cds_spa(spa_lst = spa_frame_behind)

  modle_ud_orf <- train_model_cds(
    m_pos = feat_norf_pos$feature_mat, 
    #m_neg = feat_norf_neg$feature_mat,
    threshold_method = 'youden',
    m_neg = rbind(feat_norf_neg$feature_mat,feat_norf_ahead$feature_mat, feat_norf_behind$feature_mat)
  )
rf_ud_orf <- train_model_cds_rf(
    m_pos = feat_m_pos[,-c(1:12)],
    m_neg = feature_mat_neg[,-c(1:12)]
  )
colMeans(feat_m_pos[,-c(1:5,7:12)]/feat_m_pos[,6])
colMeans(rbind(feature_mat_neg)[,-c(1:5,7:12)]/rbind(feature_mat_neg)[,6])
(rf_ud_orf$model$importance)


hist(predict(rf_ud_orf$model, newdata = data.frame(feat_m_pos), type = "prob")[, "1"], breaks=50)



cbind(predict(modle_cds$rf, data = data.frame(feature_mat_all))$predictions[, "1"]
  ,predict(modle_ud_orf$rf, data = data.frame(feature_mat_all))$predictions[, "1"], feature_mat_all[,1:12])
cbind(predict(modle_ud_orf$rf, data = data.frame(feat_norf_pos$feature_mat))$predictions[, "1"], feat_norf_pos$feature_mat)
  table(feat_norf_pos$feature_mat[,4] > feat_norf_pos$feature_mat[,3], feat_norf_pos$feature_mat[,4] > feat_norf_pos$feature_mat[,5])
  table(feat_norf_ahead$feature_mat[,4] > feat_norf_ahead$feature_mat[,3], feat_norf_ahead$feature_mat[,4] > feat_norf_ahead$feature_mat[,5])

cbind(predict(modle_cds$rf, data = data.frame(feat_m_pos))$predictions[, "1"],
predict(modle_ud_orf$rf, data = data.frame(feat_m_pos))$predictions[, "1"], feat_m_pos[,1:6])
  options(max.print = 1000)
  hist(predict(modle_ud_orf$rf, data = data.frame(feat_m_pos))$predictions[, "1"])
  
  
  # 训练uoORF doORF iORF翻译模型
  #   背景错框翻译，目标分4种情况
  # 阳性
  pick_overlap <- function(back_frame, head_info, bin_num = 10) {
    spa_frame0 <- back_frame$spa_frame0 + head_info$spa_frame$frame_0
    spa_frame1 <- back_frame$spa_frame1 + head_info$spa_frame$frame_1
    spa_frame2 <- back_frame$spa_frame2 + head_info$spa_frame$frame_2

    spa_inframe <- spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa_inframe@x <- spa_frame0@x / spa_all@x
    spa_frame0_fq <- spa_inframe
    spa_inframe@x <- spa_frame1@x / spa_all@x
    spa_frame1_fq <- spa_inframe
    spa_inframe@x <- spa_frame2@x / spa_all@x
    spa_frame2_fq <- spa_inframe

    orf_bin <- Matrix::summary(spa_all)
    orf_bin$x <- floor((orf_bin$j - 0.5) / back_frame$orf_aa[orf_bin$i] * bin_num) + 1
    spa_bin <- Matrix::sparseMatrix(
      i = orf_bin$i,
      j = orf_bin$j,
      x = orf_bin$x,
      dims = dim(spa_all)
    )
    bin_lst <- lapply(seq.int(bin_num), function(i) {
      spa_bin == i
    })
    # Summaries within positional bins
    bin_rpf_num <- function(bin_lst, spa_all) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_all * bin_i)
      })
    }
    bin_f0_num <- function(bin_lst, spa_frame0) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame0 * bin_i)
      })
    }
    bin_f1_num <- function(bin_lst, spa_frame1) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame1 * bin_i)
      })
    }
    bin_f2_num <- function(bin_lst, spa_frame2) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_frame2 * bin_i)
      })
    }
    rpf_num <- as.vector(bin_rpf_num(bin_lst, spa_all))
    bin_eff <- function(bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i)
      })
    }
    eff_pos <- as.vector(bin_eff(bin_lst))
    bin_prop_sum <- function(spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(spa_prop * bin_i)
      })
    }

    bin_weight_sum <- function(spa_ribo, spa_prop, bin_lst) {
      sapply(bin_lst, function(bin_i) {
        Matrix::rowSums(bin_i * (spa_ribo * spa_prop))
      })
    }
    norf_aa <- back_frame$orf_aa / bin_num
    spa_all_log <- spa_all
    # Log-transform ORF length and aggregated read counts
    feature_mat <- cbind(
      rep(norf_aa, bin_num),
      rpf_num,
      as.vector(bin_f0_num(bin_lst, spa_frame0)),
      as.vector(bin_f1_num(bin_lst, spa_frame1)),
      as.vector(bin_f2_num(bin_lst, spa_frame2)),
      eff_pos,
      as.vector(bin_prop_sum(spa_frame0_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame1_fq, bin_lst)),
      as.vector(bin_prop_sum(spa_frame2_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame0_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame1_fq, bin_lst)),
      as.vector(bin_weight_sum(spa_all_log, spa_frame2_fq, bin_lst))
    )
    colnames(feature_mat) <- paste0('feat_', seq_len(ncol(feature_mat)))
    idx_cds <- (norf_aa >= 6) & (rpf_num >= 40) & (eff_pos >= 13) & head_info$idx$idx_norf
    return(feature_mat[idx_cds, ])
  }

  feat_oorf_p1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_pos)
  feat_oorf_p2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_pos)

  # 阴性
  feat_oorf_n1 <- pick_overlap(back_frame = spa_frame_ahead, head_info = feat_norf_behind)
  feat_oorf_n2 <- pick_overlap(back_frame = spa_frame_behind, head_info = feat_norf_ahead)

  modle_overlap_orf <- train_model_cds(
    m_pos = rbind(feat_oorf_p1, feat_oorf_p2),
    m_neg = rbind(feat_oorf_n1, feat_oorf_n2)
  )
  return(list(
    train_model = list(
      modle_cds = modle_cds,
      modle_ud_orf = modle_ud_orf,
      modle_overlap_orf = modle_overlap_orf
    ),
    candidate_rg = candi_all_orf,
    ribo_num = ribo_num_train$ribo_num,
    feat_m = feature_mat_all
  ))
}


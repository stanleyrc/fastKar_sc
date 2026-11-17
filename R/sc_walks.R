## Single-cell-guided walk sampling.
##
## A balanced integer graph admits many walk decompositions (each node's copies can be wired to its neighbors' copies in any order).
## sample.gwalks() draws these wirings uniformly. When many single cells share the amplicon, the copies of segments that sit on one molecule
## rise and fall together across cells, and the junctions of one molecule are carried by the same cells. sc_sample_gwalks() uses that:
##
##   features   = graph nodes (per-cell copy number of the node's interval) and ALT edges (per-cell junction read rate)
##   R          = Spearman correlation of the features across cells
##   walk score = |F_w| * mean_{i<j in F_w} (R_ij - tau)   (F_w = distinct features on walk w; 0 for walks with < 2 features)
##   wiring     = sum of walk scores over all walk copies of the wiring
##
## and runs simulated annealing (then a Metropolis chain at the final temperature) over wirings with the same move as local.sampling()
## (swap two copies at one node). Walks are scored from co-variation only; their copy numbers come from the graph.
## Cell-level inputs never enter the graph itself, so the graph (CN, junctions) is fitted without them.

#' Node x cell copy-number matrix from a per-cell bin matrix
#'
#' @param gg balanced gGraph
#' @param bins data.table/data.frame chr, start, end of the bins (columns of M)
#' @param M cells x bins copy-number matrix
#' @param min_bins nodes overlapping fewer bins get NA
#' @return cells x nodes matrix, columns named by snode.id
#' @export
sc_node_matrix = function(gg, bins, M, min_bins = 1) {
  nd = gg$nodes$dt
  bgr = GenomicRanges::GRanges(sub("^chr", "", bins$chr), IRanges::IRanges(bins$start, bins$end))
  ngr = GenomicRanges::GRanges(sub("^chr", "", as.character(nd$seqnames)), IRanges::IRanges(nd$start, nd$end))
  ov = GenomicRanges::findOverlaps(ngr, bgr)
  out = matrix(NA_real_, nrow(M), nrow(nd), dimnames = list(rownames(M), as.character(nd$snode.id)))
  qh = S4Vectors::queryHits(ov); sh = S4Vectors::subjectHits(ov)
  for (i in unique(qh)) {
    b = sh[qh == i]
    if (length(b) >= min_bins) out[, i] = rowMeans(M[, b, drop = FALSE], na.rm = TRUE)
  }
  out
}

## edge lookup: signed node pair (exit a, enter b) -> edge row; walks traverse a -> b
.sc_edge_key = function(n1, s1, n2, s2) paste(n1, s1, n2, s2)
.sc_edge_table = function(gg) {
  e = gg$edges$dt[, .(edge.id, n1, n2, n1.side, n2.side, type)]
  ## a walk step from node a (signed) to node b (signed): exit side of a = right if a > 0; entry side of b = left if b > 0
  k1 = .sc_edge_key(e$n1, e$n1.side, e$n2, e$n2.side); k2 = .sc_edge_key(e$n2, e$n2.side, e$n1, e$n1.side)
  list(key = c(k1, k2), edge.id = c(e$edge.id, e$edge.id), alt = c(e$type == "ALT", e$type == "ALT"))
}
.sc_walk_features = function(sn, circular, et) {
  a = sn; b = c(sn[-1], if (circular) sn[1] else integer(0)); a = a[seq_along(b)]
  if (!length(b)) return(paste0("n", unique(abs(sn))))
  k = .sc_edge_key(abs(a), ifelse(a > 0, "right", "left"), abs(b), ifelse(b > 0, "left", "right"))
  m = match(k, et$key)
  alt = unique(et$edge.id[m[!is.na(m) & et$alt[m]]])
  c(paste0("n", unique(abs(sn))), if (length(alt)) paste0("e", alt))
}

#' Feature correlation matrix from per-cell node copy numbers and per-cell junction reads
#'
#' @param node_cn cells x nodes matrix (sc_node_matrix)
#' @param edge_reads optional cells x ALT-edge matrix of junction reads (columns named by edge.id), divided by cell depth if given
#' @param depth optional per-cell read depth (junction reads are divided by it)
#' @param min_cells features measured (finite and, for edges, non-zero) in fewer cells are dropped
#' @return list(R = feature correlation matrix with dimnames n<snode.id> / e<edge.id>)
#' @export
sc_feature_cor = function(node_cn, edge_reads = NULL, depth = NULL, min_cells = 5) {
  X = node_cn; colnames(X) = paste0("n", colnames(node_cn))
  if (!is.null(edge_reads) && ncol(edge_reads)) {
    J = edge_reads; if (!is.null(depth)) J = J / depth
    colnames(J) = paste0("e", colnames(edge_reads))
    J = J[, colSums(edge_reads > 0, na.rm = TRUE) >= min_cells, drop = FALSE]
    X = cbind(X, J[rownames(X), , drop = FALSE])
  }
  ok = colSums(is.finite(X)) >= min_cells & apply(X, 2, function(v) stats::sd(v, na.rm = TRUE)) > 0
  X = X[, ok, drop = FALSE]
  R = suppressWarnings(stats::cor(X, method = "spearman", use = "pairwise.complete.obs"))
  R[!is.finite(R)] = 0
  list(R = R)
}

#' Score walks by single-cell co-variation of their features
#' @param walks list(snode.id = list of signed node vectors, circular = logical)
#' @return numeric vector of walk scores
#' @export
sc_score_walks = function(walks, R, et, tau = 0.5) {
  fn = colnames(R)
  vapply(seq_along(walks$snode.id), function(i) {
    f = intersect(.sc_walk_features(walks$snode.id[[i]], walks$circular[i], et), fn)
    if (length(f) < 2) return(0)
    r = R[f, f, drop = FALSE]
    length(f) * (mean(r[upper.tri(r)]) - tau)
  }, 0)
}

#' Single-cell-guided walk sampling (simulated annealing + Metropolis over wirings)
#'
#' @param gg balanced integer gGraph (as for sample.gwalks)
#' @param R feature correlation matrix (sc_feature_cor()$R)
#' @param tau correlation level that counts as "same molecule" (pairs below it lower the score)
#' @param nsteps annealing steps; nsample further steps at the final temperature, recorded every thin steps
#' @param T0,T1 start and end temperature (score units)
#' @param restarts independent chains (each from a uniform random wiring)
#' @param mc.cores cores over restarts
#' @return list(best = list(score, walks), samples = data.table(chain, step, score, walk, wkey, circular), chains = per-chain traces)
#' @export
sc_sample_gwalks = function(gg, R, tau = 0.5, nsteps = 5000, nsample = 2000, thin = 20, T0 = 5, T1 = 0.05, restarts = 4,
                            mc.cores = 1, seed = 1, frozen.nodes = NULL) {
  wiring = gg.to.wiring(gg); et = .sc_edge_table(wiring$gg)
  chain = function(r) {
    set.seed(seed + r)
    internal.edges = .wiring_shuffle(data.table::copy(wiring$internal.edges))
    st = .wiring_state(internal.edges, wiring$loose.ends)
    pivots = setdiff(st$n[which(st$cn > 1)], frozen.nodes)
    prow = lapply(pivots, function(p) .wiring_rows(st, p))
    right = st$right
    w = .wiring_traverse(st, right); s = sum(sc_score_walks(w, R, et, tau))
    best = list(score = s, walks = w); trace = numeric(nsteps + nsample); rec = list()
    temps = c(exp(seq(log(T0), log(T1), length.out = max(nsteps, 1))), rep(T1, nsample))
    for (j in seq_along(temps)) {
      if (!length(pivots)) break
      k = sample.int(length(pivots), 1); r2 = .wiring_swap(right, prow[[k]])
      w2 = .wiring_traverse(st, r2); s2 = sum(sc_score_walks(w2, R, et, tau))
      if (s2 >= s || stats::runif(1) < exp((s2 - s) / temps[j])) { right = r2; s = s2; w = w2 }
      if (s > best$score) best = list(score = s, walks = w)
      trace[j] = s
      if (j > nsteps && (j - nsteps) %% thin == 0)
        rec[[length(rec) + 1]] = data.table::data.table(chain = r, step = j, score = s, walk = seq_along(w$snode.id),
          wkey = vapply(w$snode.id, function(x) paste(x, collapse = ","), ""), circular = w$circular)
    }
    list(best = best, samples = data.table::rbindlist(rec), trace = trace)
  }
  ch = parallel::mclapply(seq_len(restarts), chain, mc.cores = mc.cores)
  b = ch[[which.max(vapply(ch, function(x) x$best$score, 0))]]$best
  list(best = b, samples = data.table::rbindlist(lapply(ch, `[[`, "samples")), traces = lapply(ch, `[[`, "trace"), edge_table = et)
}

#' Walk frequencies from sc_sample_gwalks samples (or any table of walk instances with chain, step, wkey)
#' @return data.table wkey, circular, n_states (wirings containing the walk), freq, mean_copies
#' @export
sc_walk_freq = function(samples) {
  n = data.table::uniqueN(samples[, .(chain, step)])
  x = samples[, .(copies = .N), by = .(chain, step, wkey, circular)]
  x[, .(n_states = .N, mean_copies = mean(copies)), by = .(wkey, circular)][, freq := n_states / n][order(-freq)]
}

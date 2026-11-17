# Synthetic checks of single-cell-guided walk sampling (sc_walks.R).
#
# Graph: segments A, B (left ends loose), M (2 copies), C, D (right ends loose), with junctions A->M, B->M, M->C, M->D (1 copy each).
# The two wirings at M give {A-M-C, B-M-D} or {A-M-D, B-M-C}. Cells: A and C share one molecule's copy number, B and D another's.
# Uniform sampling picks each wiring half the time; the single-cell sampler should pick A-M-C / B-M-D.

sc_toy_graph <- function() {
  gr <- GenomicRanges::GRanges("1", IRanges::IRanges(c(1, 1001, 2001, 3001, 4001), width = 1000))
  gr$cn <- c(1, 1, 2, 1, 1)                         # A B M C D
  ed <- data.table::data.table(n1 = c(1, 2, 3, 3), n2 = c(3, 3, 4, 5), n1.side = "right", n2.side = "left", cn = 1, type = "ALT")
  gGnome::gG(nodes = gr, edges = ed)
}

sc_toy_cells <- function(n = 60, seed = 1) {
  set.seed(seed)
  m1 <- rlnorm(n, 3, 0.6); m2 <- rlnorm(n, 3, 0.6)   # molecule 1 (A, C) and molecule 2 (B, D), independent
  X <- cbind(`1` = m1, `2` = m2, `3` = m1 + m2, `4` = m1, `5` = m2) * exp(matrix(rnorm(n * 5, 0, 0.05), n))
  rownames(X) <- paste0("c", seq_len(n)); X
}

wiring_of <- function(keys) {
  has <- function(a, b) any(grepl(paste0("(^|,)-?", a, ",-?3,-?", b, "(,|$)"), keys) | grepl(paste0("(^|,)-?", b, ",-?3,-?", a, "(,|$)"), keys))
  if (has(1, 4) && has(2, 5)) "AC/BD" else if (has(1, 5) && has(2, 4)) "AD/BC" else "other"
}

test_that("feature correlation puts A with C and B with D", {
  skip_if_no_loosefix()
  R <- sc_feature_cor(sc_toy_cells())$R
  expect_gt(R["n1", "n4"], 0.9); expect_gt(R["n2", "n5"], 0.9); expect_lt(abs(R["n1", "n5"]), 0.4)
})

test_that("walk score prefers co-varying walks", {
  skip_if_no_loosefix()
  gg <- sc_toy_graph(); R <- sc_feature_cor(sc_toy_cells())$R; et <- fastKar:::.sc_edge_table(gg.to.wiring(gg)$gg)
  good <- list(snode.id = list(c(1L, 3L, 4L), c(2L, 3L, 5L)), circular = c(FALSE, FALSE))
  bad <- list(snode.id = list(c(1L, 3L, 5L), c(2L, 3L, 4L)), circular = c(FALSE, FALSE))
  expect_gt(sum(sc_score_walks(good, R, et)), sum(sc_score_walks(bad, R, et)))
})

test_that("sc_sample_gwalks recovers the co-varying wiring; uniform sampling does not prefer it", {
  skip_if_no_loosefix()
  gg <- sc_toy_graph(); R <- sc_feature_cor(sc_toy_cells())$R
  res <- sc_sample_gwalks(gg, R, nsteps = 200, nsample = 200, thin = 5, restarts = 2, seed = 3)
  expect_identical(wiring_of(vapply(res$best$walks$snode.id, paste, "", collapse = ",")), "AC/BD")
  st <- res$samples[, .(w = wiring_of(wkey)), by = .(chain, step)]
  expect_gt(mean(st$w == "AC/BD"), 0.9)
  set.seed(5)
  u <- sample.gwalks(gg, N = 400, return.gw = FALSE, remove.dups = FALSE, verbose = FALSE)
  uw <- vapply(u, function(w) wiring_of(vapply(w$snode.id, paste, "", collapse = ",")), "")
  expect_lt(abs(mean(uw == "AC/BD") - 0.5), 0.15)
})

test_that("sc_walk_freq counts states and copies", {
  s <- data.table::data.table(chain = 1, step = c(1, 1, 2), wkey = c("1,3,4", "1,3,4", "2,3,5"), circular = FALSE)
  f <- sc_walk_freq(s)
  expect_equal(f[wkey == "1,3,4"]$mean_copies, 2); expect_equal(f[wkey == "2,3,5"]$freq, 0.5)
})

test_that("multimers of one circle score below the monomer copies", {
  skip_if_no_loosefix()
  gg <- sc_toy_graph(); R <- sc_feature_cor(sc_toy_cells())$R; et <- fastKar:::.sc_edge_table(gg.to.wiring(gg)$gg)
  mono <- list(snode.id = list(3L, 3L), circular = c(TRUE, TRUE)); multi <- list(snode.id = list(c(3L, 3L)), circular = TRUE)
  expect_gt(sum(sc_score_walks(mono, R, et)), sum(sc_score_walks(multi, R, et)))
})

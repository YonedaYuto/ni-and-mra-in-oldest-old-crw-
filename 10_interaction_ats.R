R.version.string
library(data.table)
library(here)
library(ggplot2)

install.packages("rankFD")
library(rankFD)

if (!requireNamespace("rankFD", quietly = TRUE))
  stop("Package 'rankFD' is required for the ANOVA-type statistic: ",
       "install.packages('rankFD'). It is the reference implementation of the ",
       "rank-based factorial procedures of Brunner, Bathke and Konietschke ",
       "(2018), and the primary test of this script comes from it.")

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
for (d in c(OUT_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

IN_CSV <- file.path(OUT_DIR, "table_bm_noninferiority_mfim_subgroups.csv")

THR_SEVERE <- 26
THR_ELDER  <- 90

N_EXPECTED <- 2400

CONF_LEVEL <- 0.95

ALPHA_INTERACTION <- 0.05

ATS_EFFECT     <- "unweighted"
ATS_HYP_MAIN   <- "H0p"
ATS_HYP_SENS   <- "H0F"
ATS_CI_METHOD  <- "normal"

BOOT_HL_CI <- TRUE
BOOT_B     <- 2000
BOOT_SEED  <- 20260912

SELECT_BY_HOLM <- FALSE

LEVEL_SELECTION <- "marginal"
if (!LEVEL_SELECTION %in% c("marginal", "noninferior-levels"))
  stop("LEVEL_SELECTION must be 'marginal' or 'noninferior-levels', not '",
       LEVEL_SELECTION, "'.")
GATE_ON_MARGINAL <- identical(LEVEL_SELECTION, "marginal")

REF_WITHIN <- "within-disease"
REF_COHORT <- "cohort-remainder"
SELECT_REFERENCE <- REF_WITHIN

MIN_N <- 5

ALLOW_OVERWRITE <- FALSE

DISEASE_ORDER <- c("MSD", "CVD", "DS")
ALL_LAB       <- "All"
DISEASE_LONG  <- c(MSD = "musculoskeletal disease",
                   CVD = "cerebrovascular disease",
                   DS  = "disuse syndrome")

CARE_INDEP  <- "Independent"
CARE_NEEDED <- "Needed"
CARE_ALL    <- "All"
CARE_ORDER  <- c(CARE_INDEP, CARE_NEEDED)

QUESTION_ORDER <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")

MOD_DISEASE <- "disease"
MOD_CARE    <- "care"

TABLE1_CLASS_PCT <- c("脳血管" = 45.5, "運動器" = 40.3, "廃用" = 14.1)
TABLE1_CARE_N    <- 1051
TABLE1_CARE_PCT  <- 43.8

X_LAB_B <- "Comparator"
X_LAB_A <- "Subgroup of interest"

COL_TEXT  <- "grey15"
COL_LINE  <- "grey15"
LTY_SET   <- c("solid", "22", "12")
SHAPE_SET <- c(16, 17, 15)

Y_LIM    <- c(13, 91)
Y_BREAKS <- seq(20, 90, by = 10)

DODGE      <- 0.16
PT_SIZE    <- 2.0
LW_LINE    <- 0.55
LW_ERR     <- 0.40
ERR_WIDTH  <- 0.07

ANNOTATE_P <- TRUE
FS_BASE    <- 9
FS_ANN     <- 2.5

PANEL_W  <- 2.55
PANEL_H  <- 2.60
NCOL_MAX <- 3
MARGIN_W <- 1.1
MARGIN_H <- 2.5

if (!exists("relabel_levels")) {
  cand <- here::here("01_labels.R")
  hit <- cand[file.exists(cand)]
  if (length(hit)) {
    source(hit[1])
    message("Sourced labels from: ", hit[1])
  } else {
    warning("01_labels.R not found in: ", paste(cand, collapse = ", "),
            " - falling back to the built-in level labels.")
    LEVEL_LABELS_FALLBACK <- c("脳血管" = "CVD",
                               "運動器" = "MSD",
                               "廃用"   = "DS",
                               "なし"   = "Independent",
                               "あり"   = "Needed")
    relabel_levels <- function(var, lv) {
      unname(ifelse(lv %in% names(LEVEL_LABELS_FALLBACK),
                    LEVEL_LABELS_FALLBACK[lv], lv))
    }
  }
}

F_FAM   <- file.path(OUT_DIR, "table_interaction_ats_mfim_noninferior.csv")
F_CON   <- file.path(OUT_DIR, "table_interaction_contrasts_ats_mfim_noninferior.csv")
F_LEV   <- file.path(OUT_DIR, "table_interaction_levels_ats_mfim_noninferior.csv")
F_DUNN  <- file.path(OUT_DIR, "table_interaction_dunn_ats_mfim_noninferior.csv")
F_PTS   <- file.path(OUT_DIR, "fig_interaction_ats_mfim_points.csv")
F_QC    <- file.path(OUT_DIR, "qc_interaction_selection_sec10_ats.csv")
F_FIG_D <- file.path(FIG_DIR, "fig_interaction_ats_mfim_by_disease.tiff")
F_FIG_C <- file.path(FIG_DIR, "fig_interaction_ats_mfim_by_care.tiff")

targets <- c(F_FAM, F_CON, F_LEV, F_DUNN, F_PTS, F_QC, F_FIG_D, F_FIG_C)
exists_already <- targets[file.exists(targets)]
if (length(exists_already) && !ALLOW_OVERWRITE) {
  stop("These outputs already exist and ALLOW_OVERWRITE is FALSE:\n  ",
       paste(exists_already, collapse = "\n  "),
       "\nDelete them, rename them, or set ALLOW_OVERWRITE <- TRUE.")
}

if (!identical(SELECT_REFERENCE, REF_WITHIN))
  warning("SELECT_REFERENCE is '", SELECT_REFERENCE, "'. The interaction cells ",
          "are ALWAYS built within the modifier level, because the ",
          "cohort-remainder comparator is shared between levels and does not ",
          "give a factorial grid. The pairs admitted here would then be gated ",
          "by one comparator and tested against another. See the header.")

RANKFD_VERSION <- as.character(utils::packageVersion("rankFD"))
message("rankFD ", RANKFD_VERSION, "; ANOVA-type statistic with effect = '",
        ATS_EFFECT, "' (", if (identical(ATS_EFFECT, "unweighted"))
          "pseudo-ranks" else "classical ranks", "), primary hypothesis ",
        ATS_HYP_MAIN, ", sensitivity ", ATS_HYP_SENS, ".")

save_tiff <- function(path, plot, width, height) {
  args <- list(filename = path, plot = plot, width = width, height = height,
               units = "in", dpi = 600, compression = "lzw", bg = "white")
  if (!requireNamespace("ragg", quietly = TRUE) && isTRUE(capabilities("cairo")))
    args$type <- "cairo"
  do.call(ggplot2::ggsave, args)
  message(sprintf("Saved %s  (%.1f x %.1f in, 600 dpi, %.1f MB)",
                  path, width, height, file.size(path) / 1024^2))
}

wrap_text <- function(x, width = 150) paste(strwrap(x, width = width), collapse = "\n")

num1 <- function(v) ifelse(is.na(v), "-",
                           ifelse(v == round(v), sprintf("%.0f", v),
                                  sprintf("%.1f", v)))
fmt_p <- function(p) ifelse(is.na(p), "-",
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

if (!exists("BNB_small")) stop("BNB_small not found. Load it first.")
dat <- as.data.table(BNB_small)

need <- c("age", "class", "support_in", "mFIM_in", "mFIM_out")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("BNB_small has no column(s): ", paste(miss, collapse = ", "))

dat[, age        := as.numeric(age)]
dat[, mFIM_in    := as.numeric(mFIM_in)]
dat[, mFIM_out   := as.numeric(mFIM_out)]
dat[, class      := factor(class)]
dat[, support_in := factor(support_in)]

if (nrow(dat) != N_EXPECTED)
  warning(sprintf("BNB_small has %d rows, not the %d analysed admissions.",
                  nrow(dat), N_EXPECTED))

n_before <- nrow(dat)
dat <- dat[is.finite(age) & is.finite(mFIM_in) & is.finite(mFIM_out) &
             !is.na(class) & !is.na(support_in)]
n_dropped <- n_before - nrow(dat)
if (n_dropped > 0)
  warning(sprintf(paste0("%d row(s) dropped for a missing age, class, ",
                         "support_in, mFIM_in or mFIM_out. The plan lists all ",
                         "five as completely observed - check the input."),
                  n_dropped))

dat[, elderly := age     >= THR_ELDER]
dat[, severe  := mFIM_in <= THR_SEVERE]

resolve_levels <- function(d, col) {
  lv_raw <- levels(d[[col]])
  lv_eng <- relabel_levels(col, lv_raw)
  untranslated <- lv_raw[lv_eng == lv_raw]
  if (length(untranslated))
    warning("relabel_levels() did not translate these levels of ", col, ": ",
            paste(untranslated, collapse = ", "),
            " - they are carried through unchanged, which is visible in the ",
            "output rather than silently wrong.")
  list(raw = lv_raw, eng = lv_eng)
}

cls <- resolve_levels(dat, "class")
dat[, class_eng := factor(cls$eng[match(as.character(class), cls$raw)],
                          levels = intersect(DISEASE_ORDER, cls$eng))]

car <- resolve_levels(dat, "support_in")
dat[, care_eng := factor(car$eng[match(as.character(support_in), car$raw)],
                         levels = intersect(CARE_ORDER, car$eng))]
if (anyNA(dat$care_eng))
  warning("Some support_in values did not map onto ",
          paste(CARE_ORDER, collapse = " / "),
          " - those rows will be dropped from the stratified families.")

qc_class <- dat[, .(n = .N), by = .(raw = as.character(class))]
qc_class[, `:=`(variable            = "class",
                english             = cls$eng[match(raw, cls$raw)],
                level_position      = match(raw, cls$raw),
                table1_expected_pct = unname(TABLE1_CLASS_PCT[raw]))]
qc_care <- dat[, .(n = .N), by = .(raw = as.character(support_in))]
qc_care[, `:=`(variable            = "support_in",
               english             = car$eng[match(raw, car$raw)],
               level_position      = match(raw, car$raw))]
qc_care[, table1_expected_pct := fifelse(english == CARE_NEEDED,
                                         TABLE1_CARE_PCT, NA_real_)]
qc_lab <- rbindlist(list(qc_class, qc_care), use.names = TRUE)
qc_lab[, percent_of_cohort := round(100 * n / nrow(dat), 1)]
setorder(qc_lab, variable, level_position)

message("---- 22-1-2.0  Stratifying labels, reconciled against Table 1 ----")
print(qc_lab[, .(variable, level_position, raw, english, n,
                 percent_of_cohort, table1_expected_pct)])
bad <- qc_lab[is.finite(table1_expected_pct) &
                abs(percent_of_cohort - table1_expected_pct) > 0.5]
if (nrow(bad))
  warning("The stratum percentages do not match Table 1. Do not use these ",
          "results until this is resolved: ",
          paste(bad$variable, bad$raw, sep = "/", collapse = ", "))
n_care_needed <- dat[care_eng == CARE_NEEDED, .N]
if (n_care_needed != TABLE1_CARE_N)
  message(sprintf(paste0("Note: %d admissions are certified as needing care ",
                         "here; Table 1 reports %d."),
                  n_care_needed, TABLE1_CARE_N))

if (!file.exists(IN_CSV))
  stop("Section 08's table was not found:\n  ", IN_CSV,
       "\nRun 08_noninferiority_tests.R first.")

ni <- fread(IN_CSV)
need_cols <- c("question", "question_label", "outcome", "outcome_var",
               "disease", "reference", "care_stratum", "group_a", "group_b",
               "n_a", "n_b", "margin", "hl_shift", "hl_ci_lo", "hl_ci_hi",
               "p_noninf", "noninferior", "noninferior_holm", "note")
miss_c <- setdiff(need_cols, names(ni))
if (length(miss_c))
  stop("The input CSV is missing column(s): ", paste(miss_c, collapse = ", "),
       "\nIs ", basename(IN_CSV), " really the output of section 08 ",
       "(Brunner-Munzel version)?")

MARGIN <- unique(ni$margin)
if (length(MARGIN) != 1L || !is.finite(MARGIN))
  stop("Section 08's table does not carry one finite margin: ",
       paste(MARGIN, collapse = ", "))
M_TXT <- format(MARGIN)

bad_lv <- c(setdiff(unique(ni$disease),      c(ALL_LAB, DISEASE_ORDER)),
            setdiff(unique(ni$care_stratum), c(CARE_ALL, CARE_ORDER)),
            setdiff(unique(ni$question),     QUESTION_ORDER),
            setdiff(unique(ni$reference),    c(REF_WITHIN, REF_COHORT)))
if (length(bad_lv))
  stop("Unexpected level(s) in the input: ", paste(bad_lv, collapse = ", "),
       "\nThe family construction in this script is keyed to section 08's labels.")

DECISION_COL <- if (SELECT_BY_HOLM) "noninferior_holm" else "noninferior"
ni[, decision := as.logical(get(DECISION_COL))]

ni[, status := "inconclusive"]
ni[!is.finite(p_noninf), status := "not tested"]
ni[is.finite(p_noninf) & is.finite(hl_ci_hi) & hl_ci_hi < -MARGIN,
   status := "inferior"]
ni[decision %in% TRUE, status := "non-inferior"]

message(sprintf(paste0("Read %d contrasts from %s. Decision column: %s. ",
                       "Status: %s."),
                nrow(ni), basename(IN_CSV), DECISION_COL,
                paste(sprintf("%s %d", names(table(ni$status)),
                              as.integer(table(ni$status))), collapse = ", ")))

median_ci_binom <- function(v, conf = CONF_LEVEL) {
  v <- sort(v[is.finite(v)])
  n <- length(v)
  if (n < 6L) return(c(lo = NA_real_, hi = NA_real_, achieved = NA_real_))
  k <- stats::qbinom((1 - conf) / 2, n, 0.5)
  if (k < 1L) return(c(lo = NA_real_, hi = NA_real_, achieved = NA_real_))
  c(lo       = v[k],
    hi       = v[n - k + 1L],
    achieved = 1 - 2 * stats::pbinom(k - 1L, n, 0.5))
}

bm_rel <- function(a, b, conf = CONF_LEVEL) {
  n1 <- length(b)
  n2 <- length(a)
  ra <- rank(a)
  rb <- rank(b)
  r  <- rank(c(b, a))
  r1 <- r[seq_len(n1)]
  r2 <- r[n1 + seq_len(n2)]
  m1 <- mean(r1)
  m2 <- mean(r2)
  p_hat <- (m2 - (n2 + 1) / 2) / n1
  v1 <- sum((r1 - rb - m1 + (n1 + 1) / 2)^2) / (n1 - 1)
  v2 <- sum((r2 - ra - m2 + (n2 + 1) / 2)^2) / (n2 - 1)
  s2 <- n1 * v1 + n2 * v2
  se <- sqrt(s2) / (n1 * n2)
  if (s2 > 0) {
    df   <- s2^2 / ((n1 * v1)^2 / (n1 - 1) + (n2 * v2)^2 / (n2 - 1))
    crit <- stats::qt(1 - (1 - conf) / 2, df)
    list(p = p_hat, se = se, df = df,
         lo = p_hat - crit * se, hi = p_hat + crit * se, degenerate = FALSE)
  } else {
    list(p = p_hat, se = 0, df = NA_real_,
         lo = p_hat, hi = p_hat, degenerate = TRUE)
  }
}

local({
  ka_Y <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
  ka_N <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)
  k    <- bm_rel(a = ka_N, b = ka_Y, conf = 0.95)
  got  <- c(k$p, k$lo, k$hi, k$df)
  want <- c(0.788961, 0.5952169, 0.9827052, 17.683)
  tol  <- c(5e-7, 5e-8, 5e-8, 5e-4)
  if (any(!is.finite(got)) || any(abs(got - want) > tol))
    stop("bm_rel() does not reproduce the published Brunner-Munzel example ",
         "(estimate, CI, df = ", paste(sprintf("%.7g", got), collapse = ", "),
         "). Do not use this script until that is resolved.")
  message("Known-answer test passed: bm_rel() reproduces the Brunner & Munzel (2000) example.")
})

dunn_all_pairs <- function(y, g) {
  keep <- is.finite(y) & !is.na(g)
  y <- y[keep]; g <- droplevels(factor(g[keep]))
  N  <- length(y)
  r  <- rank(y)
  tt <- as.numeric(table(y))
  tie <- sum(tt^3 - tt)
  sig2_base <- N * (N + 1) / 12 - tie / (12 * (N - 1))
  lv  <- levels(g)
  nk  <- as.integer(table(g))
  rb  <- as.numeric(tapply(r, g, mean))
  names(nk) <- names(rb) <- lv
  cmb <- utils::combn(length(lv), 2)
  out <- data.table(
    cell_i = lv[cmb[1, ]],
    cell_j = lv[cmb[2, ]],
    n_i    = nk[cmb[1, ]],
    n_j    = nk[cmb[2, ]],
    mean_rank_i = rb[cmb[1, ]],
    mean_rank_j = rb[cmb[2, ]])
  out[, sigma := sqrt(sig2_base * (1 / n_i + 1 / n_j))]
  out[, z     := (mean_rank_i - mean_rank_j) / sigma]
  out[, p_raw := 2 * stats::pnorm(-abs(z))]
  out[, p_holm := stats::p.adjust(p_raw, method = "holm")]
  out[]
}

local({
  x <- c(1, 2, 2, 3, 3, 3, 4, 5, 5, 7, 7, 7, 9)
  y <- c(2, 3, 3, 4, 4, 6, 6, 7, 8, 9, 9)
  rs_z <- function(x, y) {
    v  <- c(x, y); r <- rank(v)
    N  <- length(v); n1 <- length(x); n2 <- length(y)
    W  <- sum(r[seq_len(n1)])
    tt <- as.numeric(table(v)); tie <- sum(tt^3 - tt)
    EW <- n1 * (N + 1) / 2
    VW <- n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1)))
    (W - EW) / sqrt(VW)
  }
  g  <- factor(rep(c("x", "y"), c(length(x), length(y))), levels = c("x", "y"))
  dz <- dunn_all_pairs(c(x, y), g)$z
  kw <- stats::kruskal.test(c(x, y), g)$statistic
  if (abs(dz - rs_z(x, y)) > 1e-10)
    stop("dunn_all_pairs() does not reproduce the tie-corrected rank-sum z ",
         "on two groups (", sprintf("%.10g vs %.10g", dz, rs_z(x, y)),
         "). Do not use this script until that is resolved.")
  if (abs(dz^2 - as.numeric(kw)) > 1e-10)
    stop("dunn_all_pairs() does not reproduce the Kruskal-Wallis chi-square ",
         "on two groups (", sprintf("%.10g vs %.10g", dz^2, as.numeric(kw)),
         "). Do not use this script until that is resolved.")
  message("Known-answer test passed: Dunn's z equals the tie-corrected ",
          "rank-sum z, and its square the Kruskal-Wallis chi-square.")
})

ats_one <- function(d, hyp) {
  out <- list(stat = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_,
              w_stat = NA_real_, w_df = NA_real_, w_p = NA_real_,
              term = NA_character_, sizes_ok = NA, msg = "")
  r <- tryCatch(
    rankFD::rankFD(y ~ level * group, data = as.data.frame(d),
                   effect = ATS_EFFECT, hypothesis = hyp,
                   CI.method = ATS_CI_METHOD, info = FALSE),
    error = function(e) {
      out$msg <<- paste0("rankFD failed (", hyp, "): ", conditionMessage(e))
      NULL
    })
  if (is.null(r)) return(out)

  A <- r[["ANOVA.Type.Statistic"]]
  W <- r[["Wald.Type.Statistic"]]
  if (is.null(A) || is.null(dim(A))) {
    out$msg <- paste0("rankFD (", hyp, ") returned no ANOVA.Type.Statistic ",
                      "matrix; components were: ",
                      paste(names(r), collapse = ", "))
    return(out)
  }

  find_row <- function(m) {
    if (is.null(m) || is.null(rownames(m))) return(integer(0))
    rn <- rownames(m)
    which(grepl(":", rn, fixed = TRUE) &
            grepl("level", rn, fixed = TRUE) &
            grepl("group", rn, fixed = TRUE))
  }
  ia <- find_row(A)
  if (length(ia) != 1L) {
    out$msg <- paste0("could not find a unique interaction row in rankFD's ",
                      "ANOVA.Type.Statistic (", hyp, "); its rows were: ",
                      paste(rownames(A), collapse = ", "))
    return(out)
  }
  out$term <- rownames(A)[ia]

  take <- function(m, i, want) {
    cn <- colnames(m)
    hit <- match(want, cn)
    if (anyNA(hit)) return(NULL)
    as.numeric(m[i, hit])
  }
  va <- take(A, ia, c("Statistic", "df1", "df2", "p-Value"))
  if (is.null(va)) {
    out$msg <- paste0("rankFD's ANOVA.Type.Statistic (", hyp, ") does not ",
                      "carry the columns Statistic / df1 / df2 / p-Value; it ",
                      "had: ", paste(colnames(A), collapse = ", "))
    return(out)
  }
  out$stat <- va[1]; out$df1 <- va[2]; out$df2 <- va[3]; out$p <- va[4]

  iw <- find_row(W)
  if (length(iw) == 1L) {
    vw <- take(W, iw, c("Statistic", "df", "p-Value"))
    if (!is.null(vw)) { out$w_stat <- vw[1]; out$w_df <- vw[2]; out$w_p <- vw[3] }
  }

  D <- r[["Descriptive"]]
  if (!is.null(D) && all(c("level", "group") %in% names(D)) &&
      "Size" %in% names(D)) {
    mine <- as.data.table(d)[, .(mine = .N), by = .(level = as.character(level),
                                                    group = as.character(group))]
    theirs <- as.data.table(D)[, .(level = as.character(level),
                                   group = as.character(group),
                                   theirs = as.integer(Size))]
    m <- merge(mine, theirs, by = c("level", "group"), all = TRUE)
    out$sizes_ok <- !anyNA(m$mine) && !anyNA(m$theirs) &&
      isTRUE(all(m$mine == m$theirs))
  }
  out
}

hl_shift <- function(a, b) as.numeric(stats::median(as.vector(outer(a, b, "-"))))

boot_hl_diffs <- function(cells_ab, B = BOOT_B, conf = CONF_LEVEL) {
  levs <- names(cells_ab)
  k    <- length(levs)
  H <- matrix(NA_real_, nrow = B, ncol = k, dimnames = list(NULL, levs))
  for (r in seq_len(B)) {
    for (j in seq_len(k)) {
      a <- cells_ab[[j]]$a
      b <- cells_ab[[j]]$b
      H[r, j] <- hl_shift(sample(a, length(a), replace = TRUE),
                          sample(b, length(b), replace = TRUE))
    }
  }
  qs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  cmb <- utils::combn(k, 2)
  data.table(
    level_i = levs[cmb[1, ]],
    level_j = levs[cmb[2, ]],
    delta_hl_boot_se = apply(cmb, 2, function(ij)
      stats::sd(H[, ij[1]] - H[, ij[2]])),
    delta_hl_ci_lo = apply(cmb, 2, function(ij)
      unname(stats::quantile(H[, ij[1]] - H[, ij[2]], qs[1], names = FALSE))),
    delta_hl_ci_hi = apply(cmb, 2, function(ij)
      unname(stats::quantile(H[, ij[1]] - H[, ij[2]], qs[2], names = FALSE))))
}

QUESTIONS <- list(
  Q1 = list(split = "elderly", outcome = "mFIM_in",  stratified = FALSE),
  Q2 = list(split = "severe",  outcome = "mFIM_out", stratified = FALSE),
  Q3 = list(split = "elderly", outcome = "mFIM_out", stratified = FALSE),
  Q4 = list(split = "elderly", outcome = "mFIM_in",  stratified = TRUE),
  Q5 = list(split = "severe",  outcome = "mFIM_out", stratified = TRUE),
  Q6 = list(split = "elderly", outcome = "mFIM_out", stratified = TRUE))

for (q in names(QUESTIONS)) {
  ov <- unique(ni[question == q, outcome_var])
  if (length(ov) && !identical(sort(ov), QUESTIONS[[q]]$outcome))
    stop("Section 08 says ", q, " has outcome '", paste(ov, collapse = "/"),
         "', this script says '", QUESTIONS[[q]]$outcome,
         "'. The two definitions of the questions have drifted apart.")
}

TWIN_UNSTRATIFIED <- vapply(names(QUESTIONS), function(q) {
  qq <- QUESTIONS[[q]]
  if (!isTRUE(qq$stratified)) return(q)
  cand <- names(QUESTIONS)[vapply(QUESTIONS, function(z)
    !isTRUE(z$stratified) && identical(z$split, qq$split) &&
      identical(z$outcome, qq$outcome), logical(1))]
  if (length(cand) == 1L) cand else NA_character_
}, character(1), USE.NAMES = TRUE)

strat_q <- names(QUESTIONS)[vapply(QUESTIONS, function(z) isTRUE(z$stratified),
                                   logical(1))]
if (GATE_ON_MARGINAL && anyNA(TWIN_UNSTRATIFIED[strat_q]))
  stop("These stratified question(s) have no unique unstratified twin with the ",
       "same split and outcome, so the marginal panel of their care-modifier ",
       "families cannot be located: ",
       paste(strat_q[is.na(TWIN_UNSTRATIFIED[strat_q])], collapse = ", "),
       ". Fix the QUESTIONS list or set LEVEL_SELECTION <- 'noninferior-levels'.")
message("Unstratified twin of each question: ",
        paste(sprintf("%s->%s", names(TWIN_UNSTRATIFIED), TWIN_UNSTRATIFIED),
              collapse = ", "))

fam_defs <- list()
add_fam <- function(...) fam_defs[[length(fam_defs) + 1L]] <<- list(...)

for (q in intersect(QUESTION_ORDER, unique(ni$question))) {
  qq <- QUESTIONS[[q]]
  for (st in intersect(c(CARE_ALL, CARE_ORDER), unique(ni[question == q, care_stratum])))
    add_fam(modifier = MOD_DISEASE, question = q, care_stratum = st,
            disease = NA_character_, levels = DISEASE_ORDER,
            gate_question = q, gate_disease = ALL_LAB, gate_care = st)
  if (qq$stratified)
    for (dz in intersect(c(ALL_LAB, DISEASE_ORDER), unique(ni[question == q, disease])))
      add_fam(modifier = MOD_CARE, question = q, care_stratum = NA_character_,
              disease = dz, levels = CARE_ORDER,
              gate_question = unname(TWIN_UNSTRATIFIED[q]),
              gate_disease = dz, gate_care = CARE_ALL)
}

cells_for <- function(fd, lev) {
  qq  <- QUESTIONS[[fd$question]]
  d   <- dat
  if (identical(fd$modifier, MOD_DISEASE)) {
    if (!identical(fd$care_stratum, CARE_ALL)) d <- d[care_eng == fd$care_stratum]
    d <- d[class_eng == lev]
    key_disease <- lev
    key_care    <- fd$care_stratum
  } else {
    d <- d[care_eng == lev]
    if (!identical(fd$disease, ALL_LAB)) d <- d[class_eng == fd$disease]
    key_disease <- fd$disease
    key_care    <- lev
  }
  flag <- d[[qq$split]]
  flag[is.na(flag)] <- FALSE
  a <- d[[qq$outcome]][flag]
  b <- d[[qq$outcome]][!flag]
  list(a = a[is.finite(a)], b = b[is.finite(b)],
       key_disease = key_disease, key_care = key_care)
}

cells_marginal <- function(fd) {
  if (is.na(fd$gate_question) || is.null(QUESTIONS[[fd$gate_question]]))
    return(list(a = numeric(0), b = numeric(0)))
  qq <- QUESTIONS[[fd$gate_question]]
  d  <- dat
  if (!identical(fd$gate_care,    CARE_ALL)) d <- d[care_eng  == fd$gate_care]
  if (!identical(fd$gate_disease, ALL_LAB))  d <- d[class_eng == fd$gate_disease]
  flag <- d[[qq$split]]
  flag[is.na(flag)] <- FALSE
  a <- d[[qq$outcome]][flag]
  b <- d[[qq$outcome]][!flag]
  list(a = a[is.finite(a)], b = b[is.finite(b)])
}

qc_rows   <- list()
fam_use   <- list()
n_mismatch <- 0L

for (fd in fam_defs) {
  fid <- if (identical(fd$modifier, MOD_DISEASE))
    sprintf("byDisease | %s | care=%s", fd$question, fd$care_stratum)
  else
    sprintf("byCare | %s | disease=%s", fd$question, fd$disease)

  g_id  <- sprintf("%s | disease=%s | care=%s",
                   fd$gate_question, fd$gate_disease, fd$gate_care)
  g_row <- if (!is.na(fd$gate_question))
    ni[reference == SELECT_REFERENCE & question == fd$gate_question &
         disease == fd$gate_disease & care_stratum == fd$gate_care] else ni[0]
  g_status <- if (nrow(g_row)) g_row$status[1] else "absent from 22-1"

  gm   <- cells_marginal(fd)
  g_na <- length(gm$a); g_nb <- length(gm$b)
  g_ok <- NA
  if (nrow(g_row)) {
    g_ok <- (g_na == g_row$n_a[1]) && (g_nb == g_row$n_b[1])
    if (!isTRUE(g_ok)) n_mismatch <- n_mismatch + 1L
  }

  gate_reason <- if (!GATE_ON_MARGINAL) ""
  else if (!nrow(g_row))
    sprintf("the marginal panel %s is not in section 08's table", g_id)
  else if (!identical(g_status, "non-inferior"))
    sprintf("the marginal panel %s has 22-1 status: %s", g_id, g_status)
  else ""

  qc_rows[[length(qc_rows) + 1L]] <- data.table(
    family = fid, role = "marginal panel (gate)", modifier = fd$modifier,
    question = fd$gate_question, level = "(modifier collapsed)",
    disease = fd$gate_disease, care_stratum = fd$gate_care,
    reference = SELECT_REFERENCE,
    status_22_1 = g_status,
    hl_shift_22_1 = if (nrow(g_row)) g_row$hl_shift[1] else NA_real_,
    hl_ci_lo_22_1 = if (nrow(g_row)) g_row$hl_ci_lo[1] else NA_real_,
    hl_ci_hi_22_1 = if (nrow(g_row)) g_row$hl_ci_hi[1] else NA_real_,
    p_noninf_22_1 = if (nrow(g_row)) g_row$p_noninf[1] else NA_real_,
    n_a_22_1 = if (nrow(g_row)) g_row$n_a[1] else NA_integer_,
    n_b_22_1 = if (nrow(g_row)) g_row$n_b[1] else NA_integer_,
    n_a_rebuilt = g_na, n_b_rebuilt = g_nb, counts_match = g_ok,
    included = GATE_ON_MARGINAL && !nzchar(gate_reason),
    exclusion_reason = if (GATE_ON_MARGINAL) gate_reason else
      "not used: LEVEL_SELECTION = 'noninferior-levels' gates on the levels")

  kept <- list()
  for (lev in fd$levels) {
    cc <- cells_for(fd, lev)
    all_row <- ni[reference == SELECT_REFERENCE & question == fd$question &
                    disease == cc$key_disease & care_stratum == cc$key_care]
    st <- if (nrow(all_row)) all_row$status[1] else "absent from 22-1"

    n_a <- length(cc$a)
    n_b <- length(cc$b)

    n_ok <- NA
    if (nrow(all_row)) {
      n_ok <- (n_a == all_row$n_a[1]) && (n_b == all_row$n_b[1])
      if (!isTRUE(n_ok)) n_mismatch <- n_mismatch + 1L
    }

    reason <- if (!nrow(all_row)) "no such contrast in section 08's table"
    else if (min(n_a, n_b) < MIN_N) sprintf("fewer than %d in a cell", MIN_N)
    else if (GATE_ON_MARGINAL) gate_reason
    else if (!identical(st, "non-inferior")) sprintf("22-1 status: %s", st)
    else ""

    qc_rows[[length(qc_rows) + 1L]] <- data.table(
      family = fid, role = "grid level", modifier = fd$modifier,
      question = fd$question,
      level = lev, disease = cc$key_disease, care_stratum = cc$key_care,
      reference = SELECT_REFERENCE,
      status_22_1 = st,
      hl_shift_22_1 = if (nrow(all_row)) all_row$hl_shift[1] else NA_real_,
      hl_ci_lo_22_1 = if (nrow(all_row)) all_row$hl_ci_lo[1] else NA_real_,
      hl_ci_hi_22_1 = if (nrow(all_row)) all_row$hl_ci_hi[1] else NA_real_,
      p_noninf_22_1 = if (nrow(all_row)) all_row$p_noninf[1] else NA_real_,
      n_a_22_1 = if (nrow(all_row)) all_row$n_a[1] else NA_integer_,
      n_b_22_1 = if (nrow(all_row)) all_row$n_b[1] else NA_integer_,
      n_a_rebuilt = n_a, n_b_rebuilt = n_b, counts_match = n_ok,
      included = !nzchar(reason),
      exclusion_reason = reason)

    if (!nzchar(reason))
      kept[[lev]] <- c(cc, list(level = lev,
                                group_a = all_row$group_a[1],
                                group_b = all_row$group_b[1],
                                hl_22_1 = all_row$hl_shift[1],
                                status_22_1 = st))
  }

  if (length(kept) >= 2L)
    fam_use[[fid]] <- list(def = fd, id = fid, cells = kept,
                           gate = list(id = g_id, status = g_status,
                                       n_a = g_na, n_b = g_nb))
}

qc <- rbindlist(qc_rows)
fwrite(qc, F_QC)

if (n_mismatch)
  stop(n_mismatch, " cell(s) rebuilt here do not have the counts section 08 ",
       "wrote for the same contrast (see counts_match in ", basename(F_QC),
       "). The strata are not being built the same way in the two scripts - ",
       "do not use this output.")
message(sprintf(paste0("Cell counts reconciled with section 08 in %d of %d ",
                       "contrasts."),
                qc[counts_match %in% TRUE, .N], qc[!is.na(counts_match), .N]))

message(sprintf(paste0("Families: %d defined, %d analysed, %d set aside. Rule: ",
                       "%s."),
                length(fam_defs), length(fam_use),
                length(fam_defs) - length(fam_use),
                if (GATE_ON_MARGINAL)
                  paste("LEVEL_SELECTION = 'marginal' - a family is admitted by",
                        "its marginal panel and then tested over every level")
                else
                  paste("LEVEL_SELECTION = 'noninferior-levels' - a family needs",
                        ">= 2 levels that are themselves non-inferior")))
if (!length(fam_use))
  message("No family qualifies, so there is nothing to test for interaction. ",
          "The selection table says why, row by row.")

if (GATE_ON_MARGINAL && length(fam_use)) {
  mixed <- vapply(fam_use, function(FM) {
    bad <- vapply(FM$cells, function(z) !identical(z$status_22_1, "non-inferior"),
                  logical(1))
    if (any(bad)) paste0(FM$id, " (", paste(sprintf(
      "%s: %s", names(FM$cells)[bad],
      vapply(FM$cells[bad], function(z) z$status_22_1, character(1))),
      collapse = "; "), ")") else NA_character_
  }, character(1))
  mixed <- mixed[!is.na(mixed)]
  if (length(mixed))
    message("Grids that include a level section 08 did not declare ",
            "non-inferior (reported, not excluded):\n  ",
            paste(mixed, collapse = "\n  "))
}

G_B <- X_LAB_B
G_A <- X_LAB_A

set.seed(BOOT_SEED)

fam_out <- list(); lev_out <- list(); dunn_out <- list()
pts_out <- list();  con_out <- list()
t_start <- Sys.time()

for (FM in fam_use) {
  fd   <- FM$def
  levs <- names(FM$cells)
  k    <- length(levs)

  long <- rbindlist(lapply(levs, function(L) {
    cc <- FM$cells[[L]]
    rbind(data.table(level = L, group = G_A, y = cc$a),
          data.table(level = L, group = G_B, y = cc$b))
  }))
  long[, level := factor(level, levels = levs)]
  long[, group := factor(group, levels = c(G_B, G_A))]
  long[, cell  := factor(paste(level, group, sep = " / "),
                         levels = as.vector(t(outer(levs, c(G_A, G_B),
                                                    paste, sep = " / "))))]

  a_main <- ats_one(long[, .(y, level, group)], ATS_HYP_MAIN)
  a_sens <- ats_one(long[, .(y, level, group)], ATS_HYP_SENS)

  kw <- stats::kruskal.test(y ~ cell, data = long)

  dn <- dunn_all_pairs(long$y, long$cell)
  dn[, family := FM$id]
  dn[, comparison_type := fifelse(
    sub(" / .*$", "", cell_i) == sub(" / .*$", "", cell_j),
    "within level (subgroup vs comparator)",
    fifelse(sub("^.* / ", "", cell_i) == sub("^.* / ", "", cell_j),
            "between levels, same group", "between levels, different group"))]
  setcolorder(dn, c("family", "comparison_type", "cell_i", "cell_j",
                    "n_i", "n_j", "mean_rank_i", "mean_rank_j",
                    "z", "p_raw", "p_holm"))
  dunn_out[[FM$id]] <- dn

  wl <- dn[comparison_type == "within level (subgroup vs comparator)"]
  wl_sig <- wl[p_holm < ALPHA_INTERACTION, sub(" / .*$", "", cell_i)]
  pattern <- if (!nrow(wl)) "-"
  else if (!length(wl_sig)) "none of the levels"
  else if (length(wl_sig) == nrow(wl)) "every level"
  else paste(wl_sig, collapse = ", ")

  re <- lapply(levs, function(L) bm_rel(FM$cells[[L]]$a, FM$cells[[L]]$b))
  names(re) <- levs
  p_i  <- vapply(re, function(z) z$p,  numeric(1))
  se_i <- vapply(re, function(z) z$se, numeric(1))
  df_i <- vapply(re, function(z) z$df, numeric(1))
  hl_i <- vapply(levs, function(L) hl_shift(FM$cells[[L]]$a, FM$cells[[L]]$b),
                 numeric(1))

  q_ok <- all(is.finite(p_i)) && all(is.finite(se_i)) && all(se_i > 0)
  if (q_ok) {
    w      <- 1 / se_i^2
    p_bar  <- sum(w * p_i) / sum(w)
    q_stat <- sum(w * (p_i - p_bar)^2)
    q_df   <- k - 1L
    q_p    <- stats::pchisq(q_stat, q_df, lower.tail = FALSE)
  } else {
    p_bar <- NA_real_; q_stat <- NA_real_; q_df <- NA_integer_; q_p <- NA_real_
  }

  cmb <- utils::combn(k, 2)
  con <- data.table(
    family = FM$id, modifier = fd$modifier, question = fd$question,
    level_i = levs[cmb[1, ]], level_j = levs[cmb[2, ]])
  con[, `:=`(
    n_a_i = vapply(level_i, function(L) length(FM$cells[[L]]$a), integer(1)),
    n_b_i = vapply(level_i, function(L) length(FM$cells[[L]]$b), integer(1)),
    n_a_j = vapply(level_j, function(L) length(FM$cells[[L]]$a), integer(1)),
    n_b_j = vapply(level_j, function(L) length(FM$cells[[L]]$b), integer(1)),
    rel_effect_i = unname(p_i[level_i]), rel_effect_j = unname(p_i[level_j]),
    hl_i = unname(hl_i[level_i]),        hl_j = unname(hl_i[level_j]))]
  con[, delta_rel_effect := rel_effect_i - rel_effect_j]
  con[, se_i_ := unname(se_i[level_i])]
  con[, se_j_ := unname(se_i[level_j])]
  con[, df_i_ := unname(df_i[level_i])]
  con[, df_j_ := unname(df_i[level_j])]
  con[, delta_se := sqrt(se_i_^2 + se_j_^2)]
  con[, delta_df := (se_i_^2 + se_j_^2)^2 /
        (se_i_^4 / df_i_ + se_j_^4 / df_j_)]
  con[, delta_ci_lo := delta_rel_effect -
        stats::qt(1 - (1 - CONF_LEVEL) / 2, delta_df) * delta_se]
  con[, delta_ci_hi := delta_rel_effect +
        stats::qt(1 - (1 - CONF_LEVEL) / 2, delta_df) * delta_se]
  con[, delta_p := 2 * stats::pt(-abs(delta_rel_effect / delta_se), delta_df)]
  con[, delta_p_holm := stats::p.adjust(delta_p, method = "holm")]
  con[, delta_hl := hl_i - hl_j]
  con[, c("se_i_", "se_j_", "df_i_", "df_j_") := NULL]

  if (BOOT_HL_CI) {
    bb <- boot_hl_diffs(FM$cells, B = BOOT_B)
    con <- merge(con, bb, by = c("level_i", "level_j"), all.x = TRUE, sort = FALSE)
    con[, boot_reps := BOOT_B]
  } else {
    con[, `:=`(delta_hl_boot_se = NA_real_, delta_hl_ci_lo = NA_real_,
               delta_hl_ci_hi = NA_real_, boot_reps = NA_integer_)]
  }
  con[, margin := MARGIN]
  con[, delta_hl_within_margin := is.finite(delta_hl_ci_lo) &
        is.finite(delta_hl_ci_hi) &
        delta_hl_ci_lo > -MARGIN & delta_hl_ci_hi < MARGIN]
  con_out[[FM$id]] <- con

  for (L in levs) {
    cc <- FM$cells[[L]]
    for (gg in c(G_B, G_A)) {
      v  <- if (identical(gg, G_A)) cc$a else cc$b
      mc <- median_ci_binom(v)
      pts_out[[length(pts_out) + 1L]] <- data.table(
        family = FM$id, modifier = fd$modifier, question = fd$question,
        outcome = unique(ni[question == fd$question, outcome])[1],
        outcome_var = QUESTIONS[[fd$question]]$outcome,
        care_stratum = cc$key_care, disease = cc$key_disease,
        level = L, group = gg, n = length(v),
        median = as.numeric(stats::median(v)),
        ci_lo = unname(mc["lo"]), ci_hi = unname(mc["hi"]),
        ci_level_achieved = unname(mc["achieved"]),
        q1 = as.numeric(stats::quantile(v, 0.25, type = 7, names = FALSE)),
        q3 = as.numeric(stats::quantile(v, 0.75, type = 7, names = FALSE)),
        mean = mean(v))
    }
    lev_out[[length(lev_out) + 1L]] <- data.table(
      family = FM$id, modifier = fd$modifier, question = fd$question,
      care_stratum = cc$key_care, disease = cc$key_disease, level = L,
      group_a = cc$group_a, group_b = cc$group_b,
      status_22_1 = cc$status_22_1,
      n_a = length(cc$a), n_b = length(cc$b),
      median_a = as.numeric(stats::median(cc$a)),
      median_b = as.numeric(stats::median(cc$b)),
      median_difference = as.numeric(stats::median(cc$a) - stats::median(cc$b)),
      hl_shift = unname(hl_i[L]), hl_shift_22_1 = cc$hl_22_1,
      rel_effect = re[[L]]$p, rel_effect_se = re[[L]]$se,
      rel_effect_ci_lo = re[[L]]$lo, rel_effect_ci_hi = re[[L]]$hi,
      rel_effect_df = re[[L]]$df,
      dunn_z_within = wl[sub(" / .*$", "", cell_i) == L, z][1],
      dunn_p_within = wl[sub(" / .*$", "", cell_i) == L, p_raw][1],
      dunn_p_holm_within = wl[sub(" / .*$", "", cell_i) == L, p_holm][1])
  }

  fam_out[[FM$id]] <- data.table(
    family = FM$id, modifier = fd$modifier, question = fd$question,
    question_label = unique(ni[question == fd$question, question_label])[1],
    outcome = unique(ni[question == fd$question, outcome])[1],
    care_stratum = if (identical(fd$modifier, MOD_DISEASE)) fd$care_stratum else "both levels compared",
    disease = if (identical(fd$modifier, MOD_CARE)) fd$disease else "levels compared",
    reference = SELECT_REFERENCE,
    levels_compared = paste(levs, collapse = " vs "),
    k_levels = k, n_total = nrow(long),
    level_selection = LEVEL_SELECTION,
    gate_panel = FM$gate$id,
    gate_status_22_1 = FM$gate$status,
    gate_n_a = FM$gate$n_a, gate_n_b = FM$gate$n_b,
    levels_all_noninferior_in_22_1 = all(vapply(
      FM$cells, function(z) identical(z$status_22_1, "non-inferior"), logical(1))),
    levels_not_noninferior_in_22_1 = {
      bad <- vapply(FM$cells, function(z)
        !identical(z$status_22_1, "non-inferior"), logical(1))
      if (any(bad)) paste(sprintf("%s (%s)", names(FM$cells)[bad],
                                  vapply(FM$cells[bad], function(z) z$status_22_1,
                                         character(1))), collapse = "; ") else ""
    },
    ats_effect = ATS_EFFECT,
    ats_ranks = if (identical(ATS_EFFECT, "unweighted")) "pseudo-ranks" else "classical ranks",
    ats_term = a_main$term,
    ats_H0p_statistic = a_main$stat, ats_H0p_df1 = a_main$df1,
    ats_H0p_df2 = a_main$df2, ats_H0p_p = a_main$p,
    wts_H0p_statistic = a_main$w_stat, wts_H0p_df = a_main$w_df,
    wts_H0p_p = a_main$w_p,
    ats_H0F_statistic = a_sens$stat, ats_H0F_df1 = a_sens$df1,
    ats_H0F_df2 = a_sens$df2, ats_H0F_p = a_sens$p,
    wts_H0F_statistic = a_sens$w_stat, wts_H0F_df = a_sens$w_df,
    wts_H0F_p = a_sens$w_p,
    rankfd_sizes_match_H0p = a_main$sizes_ok,
    kw_chisq = unname(kw$statistic), kw_df = unname(kw$parameter),
    kw_p = unname(kw$p.value),
    dunn_n_comparisons = nrow(dn),
    dunn_within_level_pattern = pattern,
    q_statistic = q_stat, q_df = q_df, q_p = q_p,
    rel_effect_pooled = p_bar,
    rel_effect_range = if (all(is.finite(p_i))) max(p_i) - min(p_i) else NA_real_,
    hl_range = max(hl_i) - min(hl_i),
    max_abs_delta_rel_effect = max(abs(con$delta_rel_effect)),
    max_abs_delta_hl = max(abs(con$delta_hl)),
    all_delta_hl_within_margin = if (BOOT_HL_CI)
      all(con$delta_hl_within_margin) else NA,
    boot_reps = if (BOOT_HL_CI) BOOT_B else NA_integer_,
    rankfd_version = RANKFD_VERSION,
    note = paste(c(if (nzchar(a_main$msg)) a_main$msg,
                   if (nzchar(a_sens$msg)) a_sens$msg,
                   if (!q_ok) paste0("the homogeneity test of the pairwise ",
                                     "relative effects was not computed: a ",
                                     "relative effect has a zero variance ",
                                     "estimate (complete separation)")),
                 collapse = "; "))
}

fam <- if (length(fam_out)) rbindlist(fam_out) else data.table()
lev <- if (length(lev_out)) rbindlist(lev_out) else data.table()
dun <- if (length(dunn_out)) rbindlist(dunn_out) else data.table()
pts <- if (length(pts_out)) rbindlist(pts_out) else data.table()
con <- if (length(con_out)) rbindlist(con_out, fill = TRUE) else data.table()

if (nrow(fam)) {
  message(sprintf("Tests complete in %.1f s (BOOT_HL_CI = %s, BOOT_B = %s).",
                  as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                  BOOT_HL_CI, format(BOOT_B)))

  fam[, ats_H0p_p_holm := stats::p.adjust(ats_H0p_p, method = "holm"), by = modifier]
  fam[, ats_H0F_p_holm := stats::p.adjust(ats_H0F_p, method = "holm"), by = modifier]
  fam[, kw_p_holm      := stats::p.adjust(kw_p,      method = "holm"), by = modifier]
  fam[, q_p_holm       := stats::p.adjust(q_p,       method = "holm"), by = modifier]
  fam[, interaction_by_ats_H0p      := ats_H0p_p      < ALPHA_INTERACTION]
  fam[, interaction_by_ats_H0p_holm := ats_H0p_p_holm < ALPHA_INTERACTION]
  fam[, interaction_by_ats_H0F      := ats_H0F_p      < ALPHA_INTERACTION]
  fam[, tests := sprintf(paste0(
    "PRIMARY: rank-based ANOVA-type statistic (rankFD %s) for the interaction ",
    "of a 2 x %d design, effect = '%s' (%s), hypothesis H0p (relative effects, ",
    "Behrens-Fisher variance), F approximation, two-sided at %s. SENSITIVITY: ",
    "the same statistic under H0F (distribution functions). The Wald-type ",
    "statistic is reported but not decided on (liberal in small samples). ",
    "ESTIMATE: the difference between two levels' pairwise relative effects ",
    "with a Welch-Satterthwaite interval, and the difference between their ",
    "Hodges-Lehmann shifts in motor FIM points with a percentile bootstrap ",
    "interval. DESCRIPTIVE ONLY: Kruskal-Wallis over the %d cells and Dunn's ",
    "comparisons - neither tests an interaction. ATS, Kruskal-Wallis and ",
    "homogeneity p values Holm-adjusted across the families of one modifier ",
    "type; no interval adjusted. Complete cases; not in the statistical ",
    "analysis plan; observational subgroups, not allocated arms."),
    RANKFD_VERSION, k_levels, ats_effect, ats_ranks,
    format(ALPHA_INTERACTION), 2L * k_levels)]
  fam[, tests := paste0(tests, " SELECTION: ", if (GATE_ON_MARGINAL) paste0(
    "the family was admitted because its marginal panel (", gate_panel,
    ", the same contrast with the modifier collapsed) was declared non-inferior ",
    "by section 08, and every level of the modifier then entered the grid, so ",
    "the contrast is the full interaction over all ", k_levels, " levels",
    fifelse(nzchar(levels_not_noninferior_in_22_1),
            paste0(" - including ", levels_not_noninferior_in_22_1,
                   ", which section 08 did not itself declare non-inferior"),
            ""), ".")
    else paste0("only the levels section 08 declared non-inferior entered the ",
                "grid (LEVEL_SELECTION = 'noninferior-levels'), so the contrast ",
                "is between those levels and not over the whole modifier."))]

  fam[, `:=`(ord_m = match(modifier, c(MOD_DISEASE, MOD_CARE)),
             ord_q = match(question, QUESTION_ORDER),
             ord_c = match(care_stratum, c(CARE_ALL, CARE_ORDER,
                                           "both levels compared")),
             ord_d = match(disease, c(ALL_LAB, DISEASE_ORDER, "levels compared")))]
  setorder(fam, ord_m, ord_q, ord_c, ord_d)
  fam[, c("ord_m", "ord_q", "ord_c", "ord_d") := NULL]
}

if (nrow(fam)) {
  n_ats <- fam[is.finite(ats_H0p_p), .N]
  message(sprintf("Self-check: the primary ATS (H0p) was computed in %d of %d families.",
                  n_ats, nrow(fam)))
  if (n_ats == 0L)
    stop("rankFD returned no usable ANOVA-type statistic for any family. The ",
         "`note` column carries its message. Do not use this output.")
  if (n_ats < nrow(fam))
    warning(nrow(fam) - n_ats, " family/families have no primary ATS; see the ",
            "`note` column of ", basename(F_FAM), ".")

  bad_sz <- fam[rankfd_sizes_match_H0p %in% FALSE, .N]
  if (bad_sz)
    stop(bad_sz, " family/families where rankFD's own Descriptive table does ",
         "not match the cell sizes built here. The data handed to the test is ",
         "not the data intended - do not use this output.")
  message(sprintf("Self-check: rankFD's cell sizes match this script's in %d of %d families.",
                  fam[rankfd_sizes_match_H0p %in% TRUE, .N],
                  fam[!is.na(rankfd_sizes_match_H0p), .N]))

  bad_df <- fam[is.finite(ats_H0p_df1) & (ats_H0p_df1 <= 0 |
                                            ats_H0p_df1 > k_levels - 1 + 1e-6), .N]
  if (bad_df)
    stop(bad_df, " family/families whose ATS numerator df is outside (0, k-1]. ",
         "The row read out of rankFD is probably not the interaction - do not ",
         "use this output.")
  message("Self-check: every ATS numerator df lies in (0, k-1], as an ",
          "interaction contrast requires.")

  bad_sign <- lev[is.finite(hl_shift) & hl_shift != 0 & is.finite(rel_effect) &
                    sign(hl_shift) != sign(rel_effect - 0.5), .N]
  message(sprintf(paste0("Self-check: the HL shift and the relative effect ",
                         "agree in sign in %d of %d levels."),
                  lev[is.finite(hl_shift) & hl_shift != 0, .N] - bad_sign,
                  lev[is.finite(hl_shift) & hl_shift != 0, .N]))
  if (bad_sign)
    stop("The Hodges-Lehmann shift and the relative effect disagree in sign in ",
         bad_sign, " level(s). One of the two is oriented the wrong way round.")

  drift <- lev[is.finite(hl_shift_22_1) &
                 abs(hl_shift - hl_shift_22_1) > 5e-3, .N]
  message(sprintf(paste0("Self-check: the HL shift reproduces section 08 in ",
                         "%d of %d levels."),
                  lev[is.finite(hl_shift_22_1), .N] - drift,
                  lev[is.finite(hl_shift_22_1), .N]))
  if (drift)
    stop(drift, " level(s) do not reproduce section 08's Hodges-Lehmann ",
         "shift. The two scripts are not looking at the same patients.")

  bad_ci <- con[is.finite(delta_ci_lo) &
                  (delta_rel_effect < delta_ci_lo |
                     delta_rel_effect > delta_ci_hi), .N]
  bad_agree <- con[is.finite(delta_p) & abs(delta_rel_effect) > 1e-12 &
                     ((delta_p < 1 - CONF_LEVEL) !=
                        (delta_ci_lo > 0 | delta_ci_hi < 0)), .N]
  if (bad_ci || bad_agree)
    stop("The interaction estimate and its interval are inconsistent in ",
         bad_ci + bad_agree, " contrast(s). Do not use this output.")
  message(sprintf(paste0("Self-check: all %d interaction contrasts have the ",
                         "estimate inside its interval, and interval and p ",
                         "value agree on whether 0 is excluded."), nrow(con)))
  if (BOOT_HL_CI) {
    bad_boot <- con[is.finite(delta_hl_ci_lo) &
                      (delta_hl < delta_hl_ci_lo | delta_hl > delta_hl_ci_hi), .N]
    if (bad_boot)
      warning(bad_boot, " bootstrap interval(s) of delta_hl do not contain the ",
              "point estimate. A percentile interval can do this when the ",
              "resampling distribution is strongly skewed; read those rows ",
              "with care.")
  }

  chk <- merge(dun[, .(n_cmp = .N,
                       n_within = sum(comparison_type ==
                                        "within level (subgroup vs comparator)")),
                   by = family],
               fam[, .(family, k_levels)], by = "family")
  chk <- merge(chk, con[, .(n_con = .N), by = family], by = "family")
  bad_cmp <- chk[n_cmp != choose(2 * k_levels, 2) | n_within != k_levels |
                   n_con != choose(k_levels, 2)]
  if (nrow(bad_cmp))
    stop("The comparison tables do not have the expected number of rows in: ",
         paste(bad_cmp$family, collapse = "; "))
  message(sprintf(paste0("Self-check: %d Dunn comparisons and %d interaction ",
                         "contrasts across %d families, as expected."),
                  nrow(dun), nrow(con), nrow(fam)))

  if (GATE_ON_MARGINAL) {
    want <- qc[role == "grid level" & status_22_1 != "absent from 22-1",
               .(k_available = .N), by = family]
    short <- merge(fam[, .(family, k_levels)], want, by = "family")[
      k_levels < k_available]
    if (nrow(short)) {
      det <- qc[role == "grid level" & included %in% FALSE &
                  family %in% short$family,
                sprintf("%s: %s dropped (%s)", family, level, exclusion_reason)]
      warning("The grid is narrower than the modifier in ", nrow(short),
              " analysed family/families, so their interaction is a comparison ",
              "of the surviving levels only:\n  ",
              paste(det, collapse = "\n  "))
    } else {
      message("Self-check: every analysed family is tested over all the levels ",
              "of its modifier that section 08 carries.")
    }
    message(sprintf(paste0("Self-check: the marginal panel of every analysed ",
                           "family is non-inferior in section 08 (%d of %d)."),
                    fam[gate_status_22_1 == "non-inferior", .N], nrow(fam)))
  }
}

PANEL_LAB <- function(d) {
  contrast <- sub(",\\s*by pre-admission care$", "",
                  sub("^[^,]*,\\s*", "", d$question_label))
  substr(contrast, 1, 1) <- toupper(substr(contrast, 1, 1))
  ifelse(d$modifier == MOD_DISEASE,
         sprintf("%s  %s\n%s\nPre-admission care: %s",
                 d$question, d$outcome, contrast, d$care_stratum),
         sprintf("%s  %s\n%s\nDisease: %s",
                 d$question, d$outcome, contrast, d$disease))
}

draw_modifier <- function(mod, path, title, subtitle) {
  fm <- fam[modifier == mod]
  if (!nrow(fm)) {
    message("No family with modifier '", mod, "'; that figure is skipped.")
    return(invisible(NULL))
  }
  fm <- copy(fm)
  fm[, panel := PANEL_LAB(fm)]
  P <- merge(pts[modifier == mod], fm[, .(family, panel)], by = "family")
  P[, panel := factor(panel, levels = fm$panel)]
  P[, group := factor(group, levels = c(G_B, G_A))]
  lev_all <- unique(P$level)
  lev_all <- c(intersect(DISEASE_ORDER, lev_all), intersect(CARE_ORDER, lev_all))
  P[, level := factor(level, levels = lev_all)]

  pd <- position_dodge(width = DODGE)
  p <- ggplot(P, aes(x = group, y = median, group = level)) +
    geom_line(aes(linetype = level), position = pd,
              linewidth = LW_LINE, colour = COL_LINE) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = pd,
                  width = ERR_WIDTH, linewidth = LW_ERR, colour = COL_LINE) +
    geom_point(aes(shape = level), position = pd,
               size = PT_SIZE, colour = COL_LINE, fill = "white") +
    scale_linetype_manual(values = setNames(LTY_SET[seq_along(lev_all)], lev_all),
                          name = NULL, drop = FALSE) +
    scale_shape_manual(values = setNames(SHAPE_SET[seq_along(lev_all)], lev_all),
                       name = NULL, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.42)) +
    facet_wrap(~ panel, ncol = min(NCOL_MAX, nrow(fm))) +
    labs(x = NULL, y = "Motor FIM (median, 95% CI)",
         title = title, subtitle = wrap_text(subtitle),
         caption = wrap_text(paste0(
           "Points are medians with their distribution-free 95% confidence ",
           "intervals (order statistics of the sign test); the levels are ",
           "offset horizontally so the intervals do not overlap. Lines that ",
           "are not parallel are what an interaction looks like; whether the ",
           "departure is more than sampling variation is what the ATS in each ",
           "panel answers, and how large an interaction the data still allow ",
           "is in ", basename(F_CON), ". ",
           if (GATE_ON_MARGINAL)
             paste0("Every level of the modifier is drawn, because the panel ",
                    "was admitted on its marginal contrast - the same ",
                    "comparison with the modifier collapsed, non-inferior at a ",
                    "margin of ", M_TXT, " motor FIM points; a level that ",
                    "section 08 did not itself declare non-inferior is named ",
                    "in ", basename(F_FAM), ". ")
           else
             paste0("Only the levels section 08 declared non-inferior at a ",
                    "margin of ", M_TXT, " motor FIM points are drawn. "),
           "Complete cases; observational subgroups, not allocated ",
           "arms; no causal reading is supported."))) +
    theme_bw(base_size = FS_BASE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = FS_BASE - 1, colour = COL_TEXT,
                                    lineheight = 1.1, margin = margin(3, 3, 3, 3)),
          axis.text = element_text(colour = COL_TEXT),
          axis.title = element_text(colour = COL_TEXT),
          legend.position = "bottom",
          legend.key.width = unit(0.42, "in"),
          plot.title = element_text(size = FS_BASE + 2, face = "bold",
                                    colour = COL_TEXT),
          plot.subtitle = element_text(size = FS_BASE - 1, colour = COL_TEXT),
          plot.caption = element_text(size = FS_BASE - 2, colour = COL_TEXT,
                                      hjust = 0),
          plot.margin = margin(8, 10, 6, 10))

  if (!is.null(Y_LIM))
    p <- p + scale_y_continuous(limits = Y_LIM, breaks = Y_BREAKS,
                                expand = expansion(mult = c(0.02, 0.30)))
  else
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.32)))

  if (ANNOTATE_P) {
    d2 <- con[, .(nk = .N), by = family]
    one <- merge(con[, .(family, delta_rel_effect, delta_ci_lo, delta_ci_hi)],
                 d2[nk == 1L, .(family)], by = "family")
    fmz <- merge(fm, one, by = "family", all.x = TRUE, sort = FALSE)
    ann <- fmz[, .(panel = factor(panel, levels = fm$panel),
                   lab = paste0(
                     sprintf("ATS (H0p) p = %s\nATS (H0F) p = %s\nKW p = %s",
                             fmt_p(ats_H0p_p), fmt_p(ats_H0F_p), fmt_p(kw_p)),
                     ifelse(is.na(delta_rel_effect), "",
                            sprintf("\nDelta rel. effect = %.3f (%.3f to %.3f)",
                                    delta_rel_effect, delta_ci_lo, delta_ci_hi))))]
    p <- p + geom_text(data = ann, inherit.aes = FALSE,
                       aes(x = -Inf, y = Inf, label = lab),
                       hjust = -0.06, vjust = 1.08, size = FS_ANN,
                       colour = COL_TEXT, lineheight = 1.05)
  }

  ncol <- min(NCOL_MAX, nrow(fm))
  nrow_f <- ceiling(nrow(fm) / ncol)
  save_tiff(path, p,
            width  = MARGIN_W + PANEL_W * ncol,
            height = MARGIN_H + PANEL_H * nrow_f)
  invisible(p)
}

SUB_COMMON <- paste0(
  "Each line joins a comparator group to its subgroup of interest inside one ",
  "level of the modifier, on the motor FIM scale. ",
  if (GATE_ON_MARGINAL)
    paste0("Each panel was admitted because its MARGINAL contrast - the same ",
           "comparison with the modifier collapsed - was declared non-inferior ",
           "by section 08 at a margin of ", M_TXT, " points, so the question ",
           "here is no longer whether a gap exceeds the margin but whether the ",
           "gap is the SAME in every level, asked of every level. ")
  else
    paste0("Every pair drawn was declared non-inferior by section 08 at a ",
           "margin of ", M_TXT, " points, so the question here is no longer ",
           "whether a gap exceeds the margin but whether the gap is the SAME ",
           "in the levels shown. "),
  "MSD, ", DISEASE_LONG[["MSD"]], "; CVD, ",
  DISEASE_LONG[["CVD"]], "; DS, ", DISEASE_LONG[["DS"]], ". ",
  "ATS is the rank-based ANOVA-type statistic for the interaction of the 2 x k ",
  "design, on pseudo-ranks: H0p tests the relative effects with the ",
  "Behrens-Fisher variance and is the primary analysis, H0F tests the ",
  "distribution functions and is the sensitivity analysis. KW p is the ",
  "Kruskal-Wallis test over all cells of the panel, which cannot separate an ",
  "interaction from the group and modifier main effects and is shown for ",
  "description only. None is adjusted for the other panels; the adjusted ",
  "values are in ", basename(F_FAM), ".")

if (nrow(fam)) {
  draw_modifier(MOD_DISEASE, F_FIG_D,
                "Motor FIM by disease class: is the subgroup-comparator gap the same in every class?",
                SUB_COMMON)
  draw_modifier(MOD_CARE, F_FIG_C,
                "Motor FIM by pre-admission care: is the subgroup-comparator gap the same at both levels?",
                SUB_COMMON)
}

round_num <- function(d, digits = 3) {
  if (!nrow(d)) return(d)
  d <- copy(d)
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  for (j in num) set(d, j = j, value = round(d[[j]], digits))
  d[]
}

if (nrow(fam)) {
  out_fam <- round_num(fam)
  for (j in c("ats_H0p_p", "ats_H0F_p", "wts_H0p_p", "wts_H0F_p", "kw_p", "q_p",
              "ats_H0p_p_holm", "ats_H0F_p_holm", "kw_p_holm", "q_p_holm"))
    set(out_fam, j = j, value = signif(fam[[j]], 3))
  front <- c("family", "modifier", "question", "question_label", "outcome",
             "disease", "care_stratum", "reference", "levels_compared",
             "k_levels", "n_total",
             "level_selection", "gate_panel", "gate_status_22_1",
             "gate_n_a", "gate_n_b",
             "levels_all_noninferior_in_22_1", "levels_not_noninferior_in_22_1",
             "ats_effect", "ats_ranks", "ats_term",
             "ats_H0p_statistic", "ats_H0p_df1", "ats_H0p_df2",
             "ats_H0p_p", "ats_H0p_p_holm",
             "interaction_by_ats_H0p", "interaction_by_ats_H0p_holm",
             "ats_H0F_statistic", "ats_H0F_df1", "ats_H0F_df2",
             "ats_H0F_p", "ats_H0F_p_holm", "interaction_by_ats_H0F",
             "max_abs_delta_rel_effect", "max_abs_delta_hl",
             "all_delta_hl_within_margin",
             "wts_H0p_statistic", "wts_H0p_df", "wts_H0p_p",
             "wts_H0F_statistic", "wts_H0F_df", "wts_H0F_p",
             "q_statistic", "q_df", "q_p", "q_p_holm",
             "kw_chisq", "kw_df", "kw_p", "kw_p_holm",
             "dunn_n_comparisons", "dunn_within_level_pattern",
             "rel_effect_pooled", "rel_effect_range", "hl_range")
  setcolorder(out_fam, c(front, setdiff(names(out_fam), front)))
  fwrite(out_fam, F_FAM)

  out_con <- round_num(con, 4)
  for (j in c("delta_p", "delta_p_holm"))
    set(out_con, j = j, value = signif(con[[j]], 3))
  front_c <- c("family", "modifier", "question", "level_i", "level_j",
               "n_a_i", "n_b_i", "n_a_j", "n_b_j",
               "rel_effect_i", "rel_effect_j", "delta_rel_effect",
               "delta_se", "delta_df", "delta_ci_lo", "delta_ci_hi",
               "delta_p", "delta_p_holm",
               "hl_i", "hl_j", "delta_hl", "delta_hl_boot_se",
               "delta_hl_ci_lo", "delta_hl_ci_hi", "boot_reps",
               "margin", "delta_hl_within_margin")
  setcolorder(out_con, c(front_c, setdiff(names(out_con), front_c)))
  fwrite(out_con, F_CON)

  out_lev <- round_num(lev, 4)
  for (j in c("dunn_p_within", "dunn_p_holm_within"))
    set(out_lev, j = j, value = signif(lev[[j]], 3))
  fwrite(out_lev, F_LEV)

  out_dun <- round_num(dun, 4)
  for (j in c("p_raw", "p_holm"))
    set(out_dun, j = j, value = signif(dun[[j]], 3))
  fwrite(out_dun, F_DUNN)

  fwrite(round_num(pts, 3), F_PTS)
} else {
  fwrite(data.table(family = character(0), modifier = character(0),
                    note = character(0)), F_FAM)
  fwrite(data.table(family = character(0), level_i = character(0),
                    level_j = character(0)), F_CON)
  fwrite(data.table(family = character(0), level = character(0)), F_LEV)
  fwrite(data.table(family = character(0), cell_i = character(0),
                    cell_j = character(0)), F_DUNN)
  fwrite(data.table(family = character(0), level = character(0),
                    group = character(0)), F_PTS)
}

if (nrow(fam)) {
  message(sprintf(paste0("\n---- 22-1-2  Interaction across the levels of each ",
                         "family, two-sided alpha %s ----"),
                  format(ALPHA_INTERACTION)))
  print(fam[, .(family, k = k_levels, n = n_total,
                ats_H0p = signif(ats_H0p_p, 3),
                ats_H0F = signif(ats_H0F_p, 3),
                kw = signif(kw_p, 3),
                interaction = interaction_by_ats_H0p,
                max_d_hl = max_abs_delta_hl,
                within_margin = all_delta_hl_within_margin)],
        nrows = 60)

  message("\n---- The interaction, estimated ----")
  print(con[, .(family, i = level_i, j = level_j,
                d_rel = round(delta_rel_effect, 3),
                lo = round(delta_ci_lo, 3), hi = round(delta_ci_hi, 3),
                d_hl = delta_hl,
                hl_lo = round(delta_hl_ci_lo, 1),
                hl_hi = round(delta_hl_ci_hi, 1))],
        nrows = 80)
}

message("\n==== Section 10 complete ====")
for (f in c(F_FAM, F_CON, F_LEV, F_DUNN, F_PTS, F_QC,
            if (file.exists(F_FIG_D)) F_FIG_D, if (file.exists(F_FIG_C)) F_FIG_C))
  message("  ", f)
message("")
message("PRIMARY test: the rank-based ANOVA-type statistic under H0p, on ",
        "pseudo-ranks - the factorial extension of the Brunner-Munzel test ",
        "section 08 decides on. SENSITIVITY: the same statistic under H0F. ",
        "Kruskal-Wallis and Dunn's comparisons are in the output because they ",
        "were asked for and they describe the cells; neither tests an ",
        "interaction, and a within-level comparison significant in one level ",
        "and not in another is not evidence of one.")
message("Read ", basename(F_CON), " before writing that there is no ",
        "interaction. delta_hl is the interaction in motor FIM points and ",
        "delta_hl_ci_lo/hi is what the data still allow; the smallest level ",
        "here can hold a few dozen admissions, and detecting an interaction ",
        "needs roughly four times the sample a main effect needs. ",
        "`delta_hl_within_margin` says whether the interval excludes an ",
        "interaction as large as the ", M_TXT, "-point margin - that, not the ",
        "p value, is the clinically readable statement.")
message("SELECTION (LEVEL_SELECTION = '", LEVEL_SELECTION, "'): ",
        if (GATE_ON_MARGINAL)
          paste0("a family was admitted by its MARGINAL panel - the same ",
                 "contrast with the modifier collapsed - and then tested over ",
                 "EVERY level of that modifier, so each p value is the full ",
                 "2 x k interaction and not a comparison of two named levels. ",
                 "`gate_panel`, `gate_status_22_1` and ",
                 "`levels_not_noninferior_in_22_1` in ", basename(F_FAM),
                 " say what admitted each family and which of its levels ",
                 "section 08 did not itself declare non-inferior; ",
                 basename(F_QC), " carries the marginal panel and every level, ",
                 "used or not.")
        else
          paste0("only the levels section 08 declared non-inferior entered ",
                 "each grid, so a p value here is a contrast between those ",
                 "levels, NOT the interaction over the whole modifier. Set ",
                 "LEVEL_SELECTION <- 'marginal' for that."))
message("rankFD ", RANKFD_VERSION, " is a NEW dependency for this project. Add ",
        "it to 06_session_info_license.R's output before the next code deposit.")
message("These are subgroups of one observational cohort. The x axis is ",
        "labelled Comparator and Subgroup of interest, not control and ",
        "intervention: nobody was allocated to anything, and no row here ",
        "supports a causal reading.")

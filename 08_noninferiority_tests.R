library(data.table)
library(here)

OUT_DIR <- here::here("data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

THR_SEVERE <- 26
THR_ELDER  <- 90

N_EXPECTED <- 2400

MARGIN <- 20

CONF_LEVEL <- 0.95

ALPHA_ONE_SIDED <- (1 - CONF_LEVEL) / 2

MIN_N   <- 5
SMALL_N <- 10

GRID_MAX <- 5000

P_ADJUST_METHOD <- "holm"

ALLOW_OVERWRITE <- FALSE

DISEASE_ORDER <- c("MSD", "CVD", "DS")
ALL_LAB       <- "All"

CARE_INDEP  <- "Independent"
CARE_NEEDED <- "Needed"
CARE_ALL    <- "All"
CARE_ORDER  <- c(CARE_INDEP, CARE_NEEDED)

REF_WITHIN <- "within-disease"
REF_COHORT <- "cohort-remainder"
REF_ORDER  <- c(REF_WITHIN, REF_COHORT)

TABLE1_CLASS_PCT <- c("脳血管" = 45.5, "運動器" = 40.3, "廃用" = 14.1)
TABLE1_CARE_N    <- 1051
TABLE1_CARE_PCT  <- 43.8

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

F_FULL    <- file.path(OUT_DIR, "table_bm_noninferiority_mfim_subgroups.csv")
F_COMPACT <- file.path(OUT_DIR, "table_bm_noninferiority_mfim_subgroups_compact.csv")
F_QC      <- file.path(OUT_DIR, "qc_bm_noninferiority_strata_sec08.csv")

targets <- c(F_FULL, F_COMPACT, F_QC)
exists_already <- targets[file.exists(targets)]
if (length(exists_already) && !ALLOW_OVERWRITE) {
  stop("These outputs already exist and ALLOW_OVERWRITE is FALSE:\n  ",
       paste(exists_already, collapse = "\n  "),
       "\nDelete them, rename them, or set ALLOW_OVERWRITE <- TRUE.")
}

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

n_non_integer <- dat[mFIM_in != round(mFIM_in) | mFIM_out != round(mFIM_out), .N]
if (n_non_integer > 0)
  warning(sprintf(paste0("%d row(s) carry a non-integer motor FIM. The motor ",
                         "FIM is an integer score - check the input."),
                  n_non_integer))

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
          " - those rows will be dropped from the stratified questions.")

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

qc <- rbindlist(list(qc_class, qc_care), use.names = TRUE)
qc[, percent_of_cohort := round(100 * n / nrow(dat), 1)]
setcolorder(qc, c("variable", "level_position", "raw", "english", "n",
                  "percent_of_cohort", "table1_expected_pct"))
setorder(qc, variable, level_position)
fwrite(qc, F_QC)

message("---- 22-1.0  Stratifying labels, reconciled against Table 1 ----")
print(qc)
bad <- qc[is.finite(table1_expected_pct) &
            abs(percent_of_cohort - table1_expected_pct) > 0.5]
if (nrow(bad))
  warning("The stratum percentages do not match Table 1. Do not use these ",
          "tests until this is resolved: ",
          paste(bad$variable, bad$raw, sep = "/", collapse = ", "))
n_care_needed <- dat[care_eng == CARE_NEEDED, .N]
if (n_care_needed != TABLE1_CARE_N)
  message(sprintf(paste0("Note: %d admissions are certified as needing care ",
                         "here; Table 1 reports %d."),
                  n_care_needed, TABLE1_CARE_N))

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

bm_core <- function(a, b, shift = 0, conf = CONF_LEVEL,
                    ra = rank(a), rb = rank(b)) {
  y  <- a + shift
  n1 <- length(b)
  n2 <- length(y)
  N  <- n1 + n2
  r  <- rank(c(b, y))
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
    stat   <- n1 * n2 * (m2 - m1) / (N * sqrt(s2))
    df     <- s2^2 / ((n1 * v1)^2 / (n1 - 1) + (n2 * v2)^2 / (n2 - 1))
    crit   <- stats::qt(1 - (1 - conf) / 2, df)
    lo     <- p_hat - crit * se
    hi     <- p_hat + crit * se
    p_one  <- stats::pt(stat, df, lower.tail = FALSE)
    p_two  <- 2 * stats::pt(-abs(stat), df)
    accept <- abs(stat) <= crit
  } else {
    stat   <- if (p_hat > 0.5) Inf else if (p_hat < 0.5) -Inf else NaN
    df     <- NA_real_
    lo     <- p_hat
    hi     <- p_hat
    p_one  <- if (p_hat > 0.5) 0 else if (p_hat < 0.5) 1 else NA_real_
    p_two  <- if (p_hat == 0.5) NA_real_ else 0
    accept <- p_hat == 0.5
  }
  list(p_hat = p_hat, stat = stat, df = df, se = se, lo = lo, hi = hi,
       p_one = p_one, p_two = p_two, accept = accept,
       degenerate = !(s2 > 0))
}

local({
  ka_Y <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
  ka_N <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)
  k    <- bm_core(a = ka_N, b = ka_Y, shift = 0, conf = 0.95)
  got  <- c(k$stat, k$df, k$p_two, k$p_hat, k$lo, k$hi)
  want <- c(3.1375, 17.683, 0.005786, 0.788961, 0.5952169, 0.9827052)
  tol  <- c(5e-5, 5e-4, 5e-7, 5e-7, 5e-8, 5e-8)
  if (any(!is.finite(got)) || any(abs(got - want) > tol))
    stop("bm_core() does not reproduce the published Brunner-Munzel example ",
         "(statistic, df, p, estimate, CI = ",
         paste(sprintf("%.7g", got), collapse = ", "),
         "). Do not use this script until that is resolved.")
  message("Known-answer test passed: bm_core() reproduces the Brunner & Munzel (2000) example.")
})

hl_ci_bm <- function(a, b, conf = CONF_LEVEL, margin = MARGIN) {
  out <- list(lo = NA_real_, hi = NA_real_, ni_by_interval = NA, msg = "")
  u <- sort(unique(as.vector(outer(unique(a), unique(b), "-"))))
  if (length(u) > GRID_MAX) {
    out$msg <- sprintf(paste0("HL interval not computed: %d distinct pairwise ",
                              "differences exceed GRID_MAX = %d"),
                       length(u), GRID_MAX)
    return(out)
  }
  mids  <- if (length(u) > 1L) (u[-1L] + u[-length(u)]) / 2 else numeric(0)
  theta <- c(u, mids)
  is_u  <- c(rep(TRUE, length(u)), rep(FALSE, length(mids)))
  if (!any(u == -margin)) {
    theta <- c(theta, -margin)
    is_u  <- c(is_u, FALSE)
  }
  o     <- order(theta)
  theta <- theta[o]
  is_u  <- is_u[o]

  ra <- rank(a)
  rb <- rank(b)
  ev <- vapply(theta, function(th) {
    k <- bm_core(a, b, shift = -th, conf = conf, ra = ra, rb = rb)
    c(stat = k$stat, accept = as.numeric(k$accept))
  }, c(stat = 0, accept = 0))
  stat   <- ev["stat", ]
  accept <- ev["accept", ] == 1

  below <- theta <= -margin
  out$ni_by_interval <- all(!accept[below] & stat[below] > 0)

  idx <- which(accept)
  if (!length(idx)) {
    out$msg <- "HL interval not computed: no shift is compatible with the data"
    return(out)
  }
  i_lo <- idx[1L]
  i_hi <- idx[length(idx)]
  out$lo <- if (is_u[i_lo] || !any(u < theta[i_lo])) theta[i_lo] else max(u[u < theta[i_lo]])
  out$hi <- if (is_u[i_hi] || !any(u > theta[i_hi])) theta[i_hi] else min(u[u > theta[i_hi]])
  if (!all(diff(idx) == 1L))
    out$msg <- paste0("the shifts not rejected by the BM test do not form one ",
                      "interval; hl_ci_lo and hl_ci_hi are its outer ends")
  out
}

bm_ni_one <- function(a, b, conf = CONF_LEVEL, margin = MARGIN) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  n1 <- length(a)
  n2 <- length(b)

  mci_a <- median_ci_binom(a, conf)
  mci_b <- median_ci_binom(b, conf)
  qtl <- function(v, p) {
    if (!length(v)) return(NA_real_)
    as.numeric(stats::quantile(v, p, type = 7, names = FALSE))
  }

  skipped <- ""
  flags   <- character(0)
  wmsg    <- character(0)
  comp    <- NULL

  if (n1 < MIN_N || n2 < MIN_N) {
    skipped <- sprintf("not tested: fewer than %d observations in one group", MIN_N)
  } else if (length(unique(c(a, b))) < 2L) {
    skipped <- "not tested: the pooled outcome takes a single value"
  } else {
    comp <- withCallingHandlers(
      tryCatch({
        at_m <- bm_core(a, b, shift = margin, conf = conf)
        no_m <- bm_core(a, b, shift = 0,      conf = conf)
        list(
          hl   = stats::median(as.vector(outer(a, b, "-"))),
          at_m = at_m,
          no_m = no_m,
          ci   = hl_ci_bm(a, b, conf = conf, margin = margin),
          chk  = max(abs(mean(outer(a + margin, b, ">")) +
                           0.5 * mean(outer(a + margin, b, "==")) - at_m$p_hat),
                     abs(mean(outer(a, b, ">")) +
                           0.5 * mean(outer(a, b, "==")) - no_m$p_hat)))
      },
      error = function(e) {
        skipped <<- paste0("Brunner-Munzel computation failed: ", conditionMessage(e))
        NULL
      }),
      warning = function(w) {
        wmsg <<- c(wmsg, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }

  hl <- hl_lo <- hl_hi <- NA_real_
  bm_stat <- bm_df <- p_ni <- NA_real_
  ps_m <- ps_m_lo <- ps_m_hi <- NA_real_
  ps_0 <- ps_0_lo <- ps_0_hi <- NA_real_
  ni_interval <- NA
  chk_pairs   <- NA_real_

  if (!is.null(comp)) {
    hl      <- comp$hl
    hl_lo   <- comp$ci$lo
    hl_hi   <- comp$ci$hi
    bm_stat <- comp$at_m$stat
    bm_df   <- comp$at_m$df
    p_ni    <- comp$at_m$p_one
    ps_m    <- comp$at_m$p_hat
    ps_m_lo <- comp$at_m$lo
    ps_m_hi <- comp$at_m$hi
    ps_0    <- comp$no_m$p_hat
    ps_0_lo <- comp$no_m$lo
    ps_0_hi <- comp$no_m$hi
    ni_interval <- comp$ci$ni_by_interval
    chk_pairs   <- comp$chk

    if (min(n1, n2) < SMALL_N)
      flags <- c(flags, sprintf(paste0("fewer than %d observations in one ",
                                       "group: the BM t approximation may be ",
                                       "liberal"), SMALL_N))
    if (comp$at_m$degenerate)
      flags <- c(flags, paste0("BM variance zero at the margin (complete ",
                               "separation): statistic infinite, df undefined, ",
                               "p is its limit"))
    if (comp$no_m$degenerate)
      flags <- c(flags, paste0("BM variance zero without the margin: the ",
                               "interval of prob_sup_nomargin is degenerate"))
    if (isTRUE(min(ps_m_lo, ps_0_lo) < 0) || isTRUE(max(ps_m_hi, ps_0_hi) > 1))
      flags <- c(flags, paste0("a probability interval extends outside [0, 1] ",
                               "(t approximation near the boundary)"))
    if (nzchar(comp$ci$msg))
      flags <- c(flags, comp$ci$msg)
    if (is.finite(hl_lo) && is.finite(hl_hi) && (hl < hl_lo || hl > hl_hi))
      flags <- c(flags, "the HL estimate lies outside its interval")
  }

  note_parts <- c(skipped, flags, wmsg)
  note_parts <- note_parts[nzchar(note_parts)]

  data.table(
    n_a = n1, n_b = n2,
    median_a = if (n1) as.numeric(stats::median(a)) else NA_real_,
    q1_a = qtl(a, 0.25), q3_a = qtl(a, 0.75),
    median_ci_lo_a = unname(mci_a["lo"]), median_ci_hi_a = unname(mci_a["hi"]),
    median_ci_level_a = unname(mci_a["achieved"]),
    min_a = if (n1) min(a) else NA_real_, max_a = if (n1) max(a) else NA_real_,
    median_b = if (n2) as.numeric(stats::median(b)) else NA_real_,
    q1_b = qtl(b, 0.25), q3_b = qtl(b, 0.75),
    median_ci_lo_b = unname(mci_b["lo"]), median_ci_hi_b = unname(mci_b["hi"]),
    median_ci_level_b = unname(mci_b["achieved"]),
    min_b = if (n2) min(b) else NA_real_, max_b = if (n2) max(b) else NA_real_,
    median_difference = if (n1 && n2)
      as.numeric(stats::median(a) - stats::median(b)) else NA_real_,
    margin = margin,
    hl_shift = hl, hl_ci_lo = hl_lo, hl_ci_hi = hl_hi,
    bm_statistic = bm_stat, bm_df = bm_df, p_noninf = p_ni,
    prob_sup_margin = ps_m,
    prob_sup_margin_ci_lo = ps_m_lo, prob_sup_margin_ci_hi = ps_m_hi,
    prob_sup_nomargin = ps_0,
    prob_sup_nomargin_ci_lo = ps_0_lo, prob_sup_nomargin_ci_hi = ps_0_hi,
    ni_by_interval = ni_interval,
    chk_pairs = chk_pairs,
    note = paste(note_parts, collapse = "; ")
  )
}

SPLITS <- list(
  elderly = list(col   = "elderly",
                 lab_a = sprintf("Oldest-old (age >= %d)", THR_ELDER),
                 lab_b = sprintf("age < %d", THR_ELDER)),
  severe  = list(col   = "severe",
                 lab_a = sprintf("Severe (admission motor FIM <= %d)", THR_SEVERE),
                 lab_b = sprintf("admission motor FIM > %d", THR_SEVERE))
)

OUTCOME_LAB <- c(mFIM_in  = "Motor FIM at admission",
                 mFIM_out = "Motor FIM at discharge")

QUESTIONS <- list(
  list(id = "Q1", split = "elderly", outcome = "mFIM_in",  stratified = FALSE,
       label = "Admission motor FIM, oldest-old vs the rest"),
  list(id = "Q2", split = "severe",  outcome = "mFIM_out", stratified = FALSE,
       label = "Discharge motor FIM, severe vs the rest"),
  list(id = "Q3", split = "elderly", outcome = "mFIM_out", stratified = FALSE,
       label = "Discharge motor FIM, oldest-old vs the rest"),
  list(id = "Q4", split = "elderly", outcome = "mFIM_in",  stratified = TRUE,
       label = "Admission motor FIM, oldest-old vs the rest, by pre-admission care"),
  list(id = "Q5", split = "severe",  outcome = "mFIM_out", stratified = TRUE,
       label = "Discharge motor FIM, severe vs the rest, by pre-admission care"),
  list(id = "Q6", split = "elderly", outcome = "mFIM_out", stratified = TRUE,
       label = "Discharge motor FIM, oldest-old vs the rest, by pre-admission care")
)

facets     <- c(ALL_LAB, intersect(DISEASE_ORDER, levels(dat$class_eng)))
care_lvls  <- intersect(CARE_ORDER, levels(dat$care_eng))

rows <- list()
for (Q in QUESTIONS) {
  sp  <- SPLITS[[Q$split]]
  strata <- if (Q$stratified) care_lvls else CARE_ALL

  for (st in strata) {
    base <- if (identical(st, CARE_ALL)) dat else dat[care_eng == st]
    for (fc in facets) {
      in_class <- if (identical(fc, ALL_LAB)) rep(TRUE, nrow(base)) else base$class_eng == fc
      in_class[is.na(in_class)] <- FALSE
      flag <- base[[sp$col]]
      idx_a <- in_class & flag

      for (rf in REF_ORDER) {
        idx_b <- if (identical(rf, REF_WITHIN)) in_class & !flag else !idx_a

        one <- bm_ni_one(base[[Q$outcome]][idx_a], base[[Q$outcome]][idx_b])

        lab_b <- if (identical(rf, REF_WITHIN)) {
          if (identical(fc, ALL_LAB)) sprintf("Others (%s)", sp$lab_b)
          else sprintf("Others within %s (%s)", fc, sp$lab_b)
        } else {
          if (identical(st, CARE_ALL))
            "Every other admission in the cohort (any disease class)"
          else
            sprintf("Every other admission in the %s stratum (any disease class)", st)
        }

        coincide <- identical(fc, ALL_LAB)
        one[, note := paste(c(if (nzchar(note)) note,
                              if (coincide)
                                "the two reference definitions coincide in this panel"),
                            collapse = "; ")]

        rows[[length(rows) + 1L]] <- cbind(
          data.table(question        = Q$id,
                     question_label  = Q$label,
                     outcome         = unname(OUTCOME_LAB[Q$outcome]),
                     outcome_var     = Q$outcome,
                     disease         = fc,
                     reference       = rf,
                     care_stratum    = st,
                     group_a         = sp$lab_a,
                     group_b         = lab_b),
          one)
      }
    }
  }
}
res <- rbindlist(rows)

res[, holm_family := paste(question, care_stratum, reference, sep = " | ")]
res[, holm_family_size := sum(is.finite(p_noninf)), by = holm_family]
if (identical(P_ADJUST_METHOD, "none")) {
  res[, p_holm := NA_real_]
} else {
  res[, p_holm := stats::p.adjust(p_noninf, method = P_ADJUST_METHOD), by = holm_family]
}
res[, p_adjust_method := if (identical(P_ADJUST_METHOD, "none")) NA_character_ else P_ADJUST_METHOD]
res[, noninferior      := p_noninf < ALPHA_ONE_SIDED]
res[, noninferior_holm := p_holm   < ALPHA_ONE_SIDED]
res[, test := sprintf(paste0(
  "Brunner-Munzel test for non-inferiority, one-sided at %s: ",
  "H0 P(A + %s > B) + 0.5 P(A + %s = B) <= 0.5 (A subgroup, B comparator); ",
  "margin %s motor FIM points (MCID); t approximation with estimated df; ",
  "HL 95%% CI by inverting the two-sided Brunner-Munzel test over shifts; ",
  "complete cases; not in the statistical analysis plan"),
  format(ALPHA_ONE_SIDED), format(MARGIN), format(MARGIN), format(MARGIN))]

setorder(res, question, care_stratum, reference, disease)
front <- c("question", "question_label", "outcome", "outcome_var",
           "disease", "reference", "care_stratum",
           "group_a", "n_a", "median_a", "q1_a", "q3_a",
           "group_b", "n_b", "median_b", "q1_b", "q3_b",
           "median_difference", "margin",
           "hl_shift", "hl_ci_lo", "hl_ci_hi",
           "bm_statistic", "bm_df", "p_noninf", "p_holm",
           "noninferior", "noninferior_holm",
           "holm_family", "holm_family_size", "p_adjust_method",
           "prob_sup_margin", "prob_sup_margin_ci_lo", "prob_sup_margin_ci_hi",
           "prob_sup_nomargin", "prob_sup_nomargin_ci_lo", "prob_sup_nomargin_ci_hi")
setcolorder(res, c(front, setdiff(names(res), front)))

add_note <- function(old, new) ifelse(nzchar(old), paste(old, new, sep = "; "), new)

if (any(is.finite(res$chk_pairs))) {
  worst <- max(res$chk_pairs, na.rm = TRUE)
  message(sprintf(paste0("Self-check: rank-based and pair-count stochastic ",
                         "superiority agree to %.3g (%d contrasts)."),
                  worst, sum(is.finite(res$chk_pairs))))
  if (worst > 1e-10)
    stop("The Brunner-Munzel estimate does not reproduce the count over all ",
         "pairs. The orientation of the test is wrong - do not use this output.")
}

tst <- res[is.finite(p_noninf) & is.finite(prob_sup_margin_ci_lo) &
             abs(prob_sup_margin_ci_lo - 0.5) > 1e-12]
if (nrow(tst)) {
  n_mism <- tst[(p_noninf < ALPHA_ONE_SIDED) != (prob_sup_margin_ci_lo > 0.5), .N]
  message(sprintf(paste0("Self-check: p value and interval of prob_sup_margin ",
                         "agree on non-inferiority in %d of %d contrasts."),
                  nrow(tst) - n_mism, nrow(tst)))
  if (n_mism)
    stop("The p value and the confidence interval of prob_sup_margin disagree ",
         "on non-inferiority in ", n_mism, " row(s). Do not use this output.")
}

disc <- res[!is.na(ni_by_interval) & !is.na(noninferior) &
              ni_by_interval != noninferior, which = TRUE]
message(sprintf(paste0("Self-check: the HL interval and the p value agree on ",
                       "non-inferiority in %d of %d contrasts."),
                res[!is.na(ni_by_interval) & !is.na(noninferior), .N] - length(disc),
                res[!is.na(ni_by_interval) & !is.na(noninferior), .N]))
if (length(disc)) {
  res[disc, note := add_note(note, paste0(
    "the HL interval and p_noninf disagree on non-inferiority; ",
    "the decision is `noninferior`"))]
  warning(length(disc), " contrast(s) where the HL interval and the p value ",
          "disagree on non-inferiority (named in `note`).")
}

edge <- res[is.finite(hl_ci_lo) & hl_ci_lo == -MARGIN, which = TRUE]
if (length(edge)) {
  res[edge, note := add_note(note, paste0(
    "hl_ci_lo equals -MARGIN exactly: read `noninferior`, not the interval"))]
  message(sprintf("Note: %d contrast(s) have an HL lower limit of exactly -%s.",
                  length(edge), format(MARGIN)))
}

res[, c("chk_pairs", "ni_by_interval") := NULL]

dup <- merge(res[disease == ALL_LAB & reference == REF_WITHIN,
                 .(question, care_stratum, p_w = p_noninf, h_w = hl_shift)],
             res[disease == ALL_LAB & reference == REF_COHORT,
                 .(question, care_stratum, p_c = p_noninf, h_c = hl_shift)],
             by = c("question", "care_stratum"))
if (nrow(dup) && (!isTRUE(all.equal(dup$p_w, dup$p_c)) ||
                  !isTRUE(all.equal(dup$h_w, dup$h_c))))
  warning("The two reference definitions disagree in the All panel, where they ",
          "should coincide by construction. Check the group construction.")

n_expected_rows <- 3L * length(facets) * length(REF_ORDER) *
  (1L + length(care_lvls))
if (nrow(res) != n_expected_rows)
  warning(sprintf(paste0("Expected %d rows, produced %d. A disease panel or a ",
                         "care stratum is missing from the data."),
                  n_expected_rows, nrow(res)))
message(sprintf("Contrasts computed: %d rows (%d tested, %d not tested); %d distinct tests.",
                nrow(res), sum(is.finite(res$p_noninf)), sum(!is.finite(res$p_noninf)),
                nrow(res) - nrow(res[disease == ALL_LAB & reference == REF_COHORT])))

round_num <- function(d, digits = 3) {
  d <- copy(d)
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  for (j in num) set(d, j = j, value = round(d[[j]], digits))
  d[]
}
out_full <- round_num(res)
out_full[, p_noninf := signif(res$p_noninf, 3)]
out_full[, p_holm   := signif(res$p_holm, 3)]
PROB_COLS <- c("prob_sup_margin", "prob_sup_margin_ci_lo", "prob_sup_margin_ci_hi",
               "prob_sup_nomargin", "prob_sup_nomargin_ci_lo", "prob_sup_nomargin_ci_hi")
for (j in PROB_COLS) set(out_full, j = j, value = round(res[[j]], 4))
fwrite(out_full, F_FULL)

num1 <- function(v) ifelse(is.na(v), "-",
                           ifelse(v == round(v), sprintf("%.0f", v),
                                  sprintf("%.1f", v)))
fmt_iqr  <- function(m, lo, hi) ifelse(is.na(m), "-",
                                       paste0(num1(m), " (", num1(lo), "-", num1(hi), ")"))
fmt_ci   <- function(e, lo, hi) ifelse(is.na(e), "-",
                                       ifelse(is.na(lo) | is.na(hi),
                                              sprintf("%.1f (interval not computed)", e),
                                              sprintf("%.1f (%.1f to %.1f)", e, lo, hi)))
fmt_prob <- function(e, lo, hi) ifelse(is.na(e), "-",
                                       sprintf("%.3f (%.3f to %.3f)", e, lo, hi))
fmt_stat <- function(s) ifelse(is.na(s), "-", sprintf("%.2f", s))
fmt_p    <- function(p) ifelse(is.na(p), "-",
                               ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_yn   <- function(v) ifelse(is.na(v), "-", ifelse(v, "yes", "no"))

compact <- res[, .(question, question_label, outcome, disease, reference,
                   care_stratum,
                   group_a, n_a,
                   median_iqr_a  = fmt_iqr(median_a, q1_a, q3_a),
                   group_b, n_b,
                   median_iqr_b  = fmt_iqr(median_b, q1_b, q3_b),
                   median_difference,
                   hl_fmt        = fmt_ci(hl_shift, hl_ci_lo, hl_ci_hi),
                   bm_stat_fmt   = fmt_stat(bm_statistic),
                   bm_df_fmt     = ifelse(is.na(bm_df), "-", sprintf("%.1f", bm_df)),
                   psup_m_fmt    = fmt_prob(prob_sup_margin,
                                            prob_sup_margin_ci_lo, prob_sup_margin_ci_hi),
                   psup_0_fmt    = fmt_prob(prob_sup_nomargin,
                                            prob_sup_nomargin_ci_lo, prob_sup_nomargin_ci_hi),
                   p_ni_fmt      = fmt_p(p_noninf),
                   p_holm_fmt    = fmt_p(p_holm),
                   ni_fmt        = fmt_yn(noninferior),
                   ni_holm_fmt   = fmt_yn(noninferior_holm),
                   note)]
setnames(compact,
         c("median_iqr_a", "median_iqr_b", "hl_fmt", "bm_stat_fmt", "bm_df_fmt",
           "psup_m_fmt", "psup_0_fmt", "p_ni_fmt", "p_holm_fmt",
           "ni_fmt", "ni_holm_fmt"),
         c("median_a (IQR)", "median_b (IQR)",
           "HL shift A-B (95% CI)",
           "BM statistic", "BM df",
           sprintf("P(A+%s > B) (95%% CI)", format(MARGIN)),
           "P(A > B) (95% CI)",
           "p non-inferiority (one-sided)", "p Holm",
           sprintf("non-inferior (margin %s)", format(MARGIN)),
           "non-inferior (Holm)"))
fwrite(compact, F_COMPACT)

message(sprintf(paste0("\n---- 22-1  Brunner-Munzel non-inferiority tests, margin %s, ",
                       "one-sided alpha %s; unadjusted and Holm-adjusted ----"),
                format(MARGIN), format(ALPHA_ONE_SIDED)))
print(res[, .(question, disease, reference, care = care_stratum, n_a, n_b,
              hl = hl_shift, hl_lo = hl_ci_lo, hl_hi = hl_ci_hi,
              bm = round(bm_statistic, 2),
              psup_m = round(prob_sup_margin, 3),
              lo = round(prob_sup_margin_ci_lo, 3),
              hi = round(prob_sup_margin_ci_hi, 3),
              p = signif(p_noninf, 3), p_holm = signif(p_holm, 3),
              NI = noninferior)],
      nrows = 100)

message("\n==== Section 08 complete ====")
message("  ", F_FULL)
message("  ", F_COMPACT)
message("  ", F_QC)
message("")
message("Check before use: qc_bm_noninferiority_strata_sec08.csv must reproduce ",
        "Table 1 (CVD 45.5, MSD 40.3, DS 14.1; care needed 43.8). If it does ",
        "not, the strata are mislabelled and nothing computed here should be used.")
message("Read `reference` before quoting a by-disease row: within-disease and ",
        "cohort-remainder answer different questions, and Figure 3 uses the latter.")
message("Section 09 still reads the Mann-Whitney table of the earlier version ",
        "and draws its null line at 0; it does not show these results.")

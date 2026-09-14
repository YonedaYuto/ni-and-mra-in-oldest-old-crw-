library(data.table)
library(here)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE))
  stop("Package 'patchwork' is required: install.packages('patchwork')")

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
for (d in c(OUT_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

IN_CSV     <- file.path(OUT_DIR, "table_bm_noninferiority_mfim_subgroups.csv")
IN_COMPACT <- file.path(OUT_DIR, "table_bm_noninferiority_mfim_subgroups_compact.csv")

ALLOW_OVERWRITE <- TRUE

ALPHA_ONE_SIDED <- 0.025

DISEASE_ORDER  <- c("All", "MSD", "CVD", "DS")
CARE_ORDER     <- c("All", "Independent", "Needed")
QUESTION_ORDER <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")
REF_WITHIN <- "within-disease"
REF_COHORT <- "cohort-remainder"

DISEASE_LONG <- c(MSD = "musculoskeletal disease",
                  CVD = "cerebrovascular disease",
                  DS  = "disuse syndrome")

COL_TEXT   <- "grey15"
COL_RULE   <- "grey35"
COL_NULL   <- "grey45"
COL_MARK   <- "black"
COL_MUTED  <- "grey45"

SHAPE_SUB <- 15
SHAPE_ALL <- 18
SIZE_SUB  <- 1.9
SIZE_ALL  <- 3.0
LW_CI     <- 0.45

FIG_WIDTH      <- 13.0
HEIGHT_PER_ROW <- 0.17
HEIGHT_BASE    <- 2.3
PANEL_WIDTHS   <- c(4.6, 3.4, 2.0)

FS_BASE <- 9
FS_ROW  <- 2.85
FS_HEAD <- 2.85

X_LABEL <- 0.000
X_NA    <- 0.520
X_MEDA  <- 0.715
X_NB    <- 0.805
X_MEDB  <- 1.000
X_NWID  <- 0.055

X_R_HL <- 0.250
X_R_P  <- 0.643

NULL_X   <- 0.5
X_LIM    <- c(-0.02, 1.02)
X_BREAKS <- seq(0, 1, by = 0.25)

set.seed(20260909)

wrap_text <- function(x, width = 150) paste(strwrap(x, width = width), collapse = "\n")

num1 <- function(v) ifelse(is.na(v), "-",
                           ifelse(v == round(v), sprintf("%.0f", v),
                                  sprintf("%.1f", v)))
fmt_n   <- function(n) ifelse(is.na(n), "-", format(n, big.mark = ",", trim = TRUE))
fmt_iqr <- function(m, lo, hi) ifelse(is.na(m), "-",
                                      paste0(num1(m), " (", num1(lo), "-", num1(hi), ")"))
fmt_hl  <- function(e) ifelse(is.na(e), "-", sprintf("%.1f", e))

F_WITHIN <- file.path(FIG_DIR, "fig_forest_mfim_within_disease.tiff")
F_COHORT <- file.path(FIG_DIR, "fig_forest_mfim_cohort_remainder.tiff")
F_ROWS   <- file.path(OUT_DIR, "fig_forest_mfim_rows.csv")

targets <- c(F_WITHIN, F_COHORT, F_ROWS)
exists_already <- targets[file.exists(targets)]
if (length(exists_already) && !ALLOW_OVERWRITE) {
  stop("These outputs already exist and ALLOW_OVERWRITE is FALSE:\n  ",
       paste(exists_already, collapse = "\n  "),
       "\nDelete them, rename them, or set ALLOW_OVERWRITE <- TRUE.")
}

save_tiff <- function(path, plot, width, height) {
  args <- list(filename = path, plot = plot, width = width, height = height,
               units = "in", dpi = 600, compression = "lzw", bg = "white")
  if (!requireNamespace("ragg", quietly = TRUE) && isTRUE(capabilities("cairo")))
    args$type <- "cairo"
  do.call(ggplot2::ggsave, args)
  message(sprintf("Saved %s  (%.1f x %.1f in, 600 dpi, %.1f MB)",
                  path, width, height, file.size(path) / 1024^2))
}

for (f in c(IN_CSV, IN_COMPACT))
  if (!file.exists(f))
    stop("Section 08's table was not found:\n  ", f,
         "\nRun 08_noninferiority_tests.R first.")

res <- fread(IN_CSV)

need <- c("question", "question_label", "outcome", "disease", "reference",
          "care_stratum", "group_a", "group_b", "n_a", "n_b",
          "median_a", "q1_a", "q3_a", "median_b", "q1_b", "q3_b",
          "margin", "hl_shift",
          "prob_sup_margin", "prob_sup_margin_ci_lo", "prob_sup_margin_ci_hi",
          "test", "note")
miss <- setdiff(need, names(res))
if (length(miss))
  stop("The input CSV is missing column(s): ", paste(miss, collapse = ", "),
       "\nIs ", basename(IN_CSV), " really the output of section 08 (Brunner-Munzel version)?")

MARGIN <- unique(res$margin)
if (length(MARGIN) != 1L || !is.finite(MARGIN))
  stop("Section 08's table does not carry one finite margin: ",
       paste(MARGIN, collapse = ", "))
M_TXT   <- format(MARGIN)
P_LABEL <- sprintf("P(A+%s > B)", M_TXT)

if (!all(grepl(sprintf("one-sided at %s", format(ALPHA_ONE_SIDED)), res$test, fixed = TRUE)))
  stop("Section 08's `test` column does not say one-sided at ",
       format(ALPHA_ONE_SIDED), ". The subtitle of this figure would be wrong; ",
       "set ALPHA_ONE_SIDED to the level section 08 used.")

bad_d <- setdiff(unique(res$disease),      DISEASE_ORDER)
bad_c <- setdiff(unique(res$care_stratum), CARE_ORDER)
bad_q <- setdiff(unique(res$question),     QUESTION_ORDER)
bad_r <- setdiff(unique(res$reference),    c(REF_WITHIN, REF_COHORT))
if (length(c(bad_d, bad_c, bad_q, bad_r)))
  stop("Unexpected level(s) in the input: ",
       paste(c(bad_d, bad_c, bad_q, bad_r), collapse = ", "),
       "\nThe row order in this script is keyed to section 08's labels.")

cmp   <- fread(IN_COMPACT)
P_COL <- sprintf("P(A+%s > B) (95%% CI)", M_TXT)
KEYS  <- c("question", "disease", "reference", "care_stratum")
miss_c <- setdiff(c(KEYS, "n_a", "n_b", P_COL), names(cmp))
if (length(miss_c))
  stop("The compact CSV is missing column(s): ", paste(miss_c, collapse = ", "),
       "\nIs ", basename(IN_COMPACT), " really the output of section 08?")

chk <- merge(res[, c(KEYS, "n_a", "n_b"), with = FALSE],
             cmp[, c(KEYS, "n_a", "n_b"), with = FALSE],
             by = KEYS, all = TRUE, suffixes = c("_full", "_compact"))
same_run <- nrow(res) == nrow(cmp) && nrow(chk) == nrow(res) &&
  !anyNA(chk$n_a_full) && !anyNA(chk$n_a_compact) &&
  isTRUE(all(chk$n_a_full == chk$n_a_compact)) &&
  isTRUE(all(chk$n_b_full == chk$n_b_compact))
if (!same_run)
  stop("The full and compact tables of section 08 do not describe the same ",
       "contrasts. Re-run section 08 so that both come from one run.")

res <- merge(res, cmp[, c(KEYS, P_COL), with = FALSE], by = KEYS, all.x = TRUE)
setnames(res, P_COL, "p_ci_text")

res[, disease      := factor(disease,      levels = DISEASE_ORDER)]
res[, care_stratum := factor(care_stratum, levels = CARE_ORDER)]
res[, question     := factor(question,     levels = QUESTION_ORDER)]
setorder(res, question, care_stratum, disease)

message(sprintf("Read %d contrasts from %s (%d within-disease, %d cohort-remainder).",
                nrow(res), basename(IN_CSV),
                res[reference == REF_WITHIN, .N], res[reference == REF_COHORT, .N]))
n_untested <- res[!(is.finite(prob_sup_margin) & is.finite(prob_sup_margin_ci_lo) &
                      is.finite(prob_sup_margin_ci_hi)), .N]
if (n_untested)
  message(sprintf(paste0("%d contrast(s) have no interval and will be printed ",
                         "without a marker (see the `note` column of the input)."),
                  n_untested))

build_lines <- function(d) {
  out <- list()
  add <- function(...) out[[length(out) + 1L]] <<- data.table(...)

  for (q in levels(droplevels(d$question))) {
    dq <- d[question == q]
    if (!nrow(dq)) next

    qlab <- sub(",\\s*by pre-admission care$", "", unique(dq$question_label)[1])
    add(kind = "qhead",
        label = sprintf("%s.  %s", q, qlab),
        n_a = NA_character_, med_a = NA_character_,
        n_b = NA_character_, med_b = NA_character_,
        hl = NA_character_, p_ci = NA_character_,
        est = NA_real_, lo = NA_real_, hi = NA_real_,
        is_all = FALSE, tested = FALSE)

    strata <- levels(droplevels(dq$care_stratum))
    stratified <- !identical(strata, "All")

    for (st in strata) {
      ds <- dq[care_stratum == st]
      if (!nrow(ds)) next
      if (stratified)
        add(kind = "chead",
            label = sprintf("   Pre-admission care: %s", st),
            n_a = NA_character_, med_a = NA_character_,
            n_b = NA_character_, med_b = NA_character_,
            hl = NA_character_, p_ci = NA_character_,
            est = NA_real_, lo = NA_real_, hi = NA_real_,
            is_all = FALSE, tested = FALSE)

      setorder(ds, disease)
      indent <- if (stratified) "      " else "   "
      for (i in seq_len(nrow(ds))) {
        r <- ds[i]
        est_i <- as.numeric(r$prob_sup_margin)
        lo_i  <- as.numeric(r$prob_sup_margin_ci_lo)
        hi_i  <- as.numeric(r$prob_sup_margin_ci_hi)
        add(kind = "data",
            label = paste0(indent, as.character(r$disease)),
            n_a   = fmt_n(r$n_a),
            med_a = fmt_iqr(r$median_a, r$q1_a, r$q3_a),
            n_b   = fmt_n(r$n_b),
            med_b = fmt_iqr(r$median_b, r$q1_b, r$q3_b),
            hl    = fmt_hl(r$hl_shift),
            p_ci  = ifelse(is.na(r$p_ci_text) | !nzchar(r$p_ci_text), "-",
                           as.character(r$p_ci_text)),
            est   = est_i,
            lo    = lo_i,
            hi    = hi_i,
            is_all = identical(as.character(r$disease), "All"),
            tested = is.finite(est_i) && is.finite(lo_i) && is.finite(hi_i))
      }
    }

    add(kind = "spacer", label = "",
        n_a = NA_character_, med_a = NA_character_,
        n_b = NA_character_, med_b = NA_character_,
        hl = NA_character_, p_ci = NA_character_,
        est = NA_real_, lo = NA_real_, hi = NA_real_,
        is_all = FALSE, tested = FALSE)
  }

  L <- rbindlist(out)
  if (nrow(L) && L[nrow(L), kind] == "spacer") L <- L[-nrow(L)]
  L[, ord := .I]
  L[, y := -ord]
  L[]
}

Y_HDR1 <-  0.75
Y_HDR2 <- -0.05
Y_RULE <- -0.55

theme_panel <- function() {
  theme_void(base_size = FS_BASE) +
    theme(plot.margin = margin(2, 2, 2, 2),
          plot.title  = element_blank())
}

y_scale <- function(L) {
  scale_y_continuous(limits = c(min(L$y) - 0.9, Y_HDR1 + 0.9),
                     expand = c(0, 0))
}

panel_left <- function(L) {
  span_a <- c(X_NA - X_NWID, X_MEDA)
  span_b <- c(X_NB - X_NWID, X_MEDB)
  hdr1 <- data.table(x0  = c(span_a[1], span_b[1]),
                     x1  = c(span_a[2], span_b[2]),
                     lab = c("Subgroup of interest", "Comparator"))
  hdr1[, x := (x0 + x1) / 2]
  hdr2 <- data.table(x    = c(X_LABEL, X_NA, X_MEDA, X_NB, X_MEDB),
                     hj   = c(0, 1, 1, 1, 1),
                     lab  = c("Subgroup", "No.", "Median (IQR)",
                              "No.", "Median (IQR)"))
  ggplot() +
    geom_text(data = hdr1, aes(x = x, y = Y_HDR1, label = lab),
              hjust = 0.5, size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_segment(data = hdr1, aes(x = x0, xend = x1,
                                  y = Y_HDR1 - 0.45, yend = Y_HDR1 - 0.45),
                 linewidth = 0.25, colour = COL_RULE) +
    geom_text(data = hdr2, aes(x = x, y = Y_HDR2, label = lab, hjust = hj),
              size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_segment(aes(x = 0, xend = 1, y = Y_RULE, yend = Y_RULE),
                 linewidth = 0.35, colour = COL_RULE) +
    geom_text(data = L[kind == "qhead"], aes(x = X_LABEL, y = y, label = label),
              hjust = 0, size = FS_ROW, fontface = "bold", colour = COL_TEXT) +
    geom_text(data = L[kind == "chead"], aes(x = X_LABEL, y = y, label = label),
              hjust = 0, size = FS_ROW, fontface = "italic", colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_LABEL, y = y, label = label),
              hjust = 0, size = FS_ROW, colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_NA,   y = y, label = n_a),
              hjust = 1, size = FS_ROW, colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_MEDA, y = y, label = med_a),
              hjust = 1, size = FS_ROW, colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_NB,   y = y, label = n_b),
              hjust = 1, size = FS_ROW, colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_MEDB, y = y, label = med_b),
              hjust = 1, size = FS_ROW, colour = COL_TEXT) +
    scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    y_scale(L) + coord_cartesian(clip = "off") + theme_panel()
}

panel_forest <- function(L, xlim, xbreaks) {
  D <- L[kind == "data" & tested == TRUE]
  D[, lo_d := pmax(lo, xlim[1])]
  D[, hi_d := pmin(hi, xlim[2])]
  D[, cut_lo := lo < xlim[1]]
  D[, cut_hi := hi > xlim[2]]
  D[, est_d := pmin(pmax(est, xlim[1]), xlim[2])]

  p <- ggplot() +
    geom_segment(aes(x = xlim[1], xend = xlim[2], y = Y_RULE, yend = Y_RULE),
                 linewidth = 0.35, colour = COL_RULE) +
    geom_segment(aes(x = NULL_X, xend = NULL_X, y = Y_RULE, yend = min(L$y) - 0.9),
                 linetype = "dashed", linewidth = 0.4, colour = COL_NULL) +
    geom_segment(data = D[cut_lo == FALSE & cut_hi == FALSE],
                 aes(x = lo_d, xend = hi_d, y = y, yend = y),
                 linewidth = LW_CI, colour = COL_MARK)
  if (nrow(D[cut_lo == TRUE]))
    p <- p + geom_segment(data = D[cut_lo == TRUE],
                          aes(x = hi_d, xend = lo_d, y = y, yend = y),
                          linewidth = LW_CI, colour = COL_MARK,
                          arrow = arrow(length = unit(0.055, "in"), type = "closed"))
  if (nrow(D[cut_hi == TRUE]))
    p <- p + geom_segment(data = D[cut_hi == TRUE],
                          aes(x = lo_d, xend = hi_d, y = y, yend = y),
                          linewidth = LW_CI, colour = COL_MARK,
                          arrow = arrow(length = unit(0.055, "in"), type = "closed"))

  p +
    geom_point(data = D[is_all == FALSE], aes(x = est_d, y = y),
               shape = SHAPE_SUB, size = SIZE_SUB, colour = COL_MARK) +
    geom_point(data = D[is_all == TRUE], aes(x = est_d, y = y),
               shape = SHAPE_ALL, size = SIZE_ALL, colour = COL_MARK) +
    geom_text(data = L[kind == "data" & tested == FALSE],
              aes(x = xlim[1] + 0.25 * diff(xlim), y = y), label = "not tested",
              size = FS_ROW - 0.35, colour = COL_MUTED, hjust = 0.5) +
    scale_x_continuous(limits = xlim, breaks = xbreaks, expand = c(0, 0)) +
    y_scale(L) +
    labs(x = sprintf(paste0("<-- Shortfall of %s or more          %s          ",
                            "Shortfall of less than %s -->"),
                     M_TXT, P_LABEL, M_TXT)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = FS_BASE) +
    theme(axis.line.x  = element_line(linewidth = 0.35, colour = COL_RULE),
          axis.ticks.x = element_line(linewidth = 0.35, colour = COL_RULE),
          axis.ticks.length.x = unit(0.06, "cm"),
          axis.text.x  = element_text(size = FS_BASE - 1, colour = COL_TEXT,
                                      margin = margin(t = 2)),
          axis.title.x = element_text(size = FS_BASE - 1, colour = COL_TEXT,
                                      margin = margin(t = 4)),
          plot.margin  = margin(2, 6, 2, 6))
}

panel_right <- function(L) {
  ggplot() +
    geom_text(aes(x = X_R_HL, y = Y_HDR1, label = "HL shift"),
              hjust = 1, size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_text(aes(x = X_R_HL, y = Y_HDR2, label = "(points)"),
              hjust = 1, size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_text(aes(x = X_R_P, y = Y_HDR1, label = P_LABEL),
              hjust = 0.5, size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_text(aes(x = X_R_P, y = Y_HDR2, label = "(95% CI)"),
              hjust = 0.5, size = FS_HEAD, fontface = "bold", colour = COL_TEXT) +
    geom_segment(aes(x = 0, xend = 1, y = Y_RULE, yend = Y_RULE),
                 linewidth = 0.35, colour = COL_RULE) +
    geom_text(data = L[kind == "data"], aes(x = X_R_HL, y = y, label = hl),
              hjust = 1, size = FS_ROW, colour = COL_TEXT) +
    geom_text(data = L[kind == "data"], aes(x = X_R_P, y = y, label = p_ci),
              hjust = 0.5, size = FS_ROW, colour = COL_TEXT) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale(L) + coord_cartesian(clip = "off") + theme_panel()
}

TITLE_TAIL <- c("each subgroup against the rest of its OWN disease class",
                "each subgroup against the rest of the COHORT")
names(TITLE_TAIL) <- c(REF_WITHIN, REF_COHORT)

COMPARATOR_NOTE <- c(
  "The comparator of every row is the remainder of the SAME disease class.",
  paste0("The comparator of every row is every OTHER admission in the cohort ",
         "- other disease classes included - and, in the stratified questions, ",
         "every other admission in the same pre-admission care stratum."))
names(COMPARATOR_NOTE) <- c(REF_WITHIN, REF_COHORT)

FILE_FOR <- c(F_WITHIN, F_COHORT)
names(FILE_FOR) <- c(REF_WITHIN, REF_COHORT)

drawn <- list()
for (rf in c(REF_WITHIN, REF_COHORT)) {
  d <- res[reference == rf]
  if (!nrow(d)) {
    warning("No rows for reference '", rf, "'; that figure is skipped.")
    next
  }
  L <- build_lines(d)

  if (is.null(X_LIM)) {
    rng <- range(c(L$lo, L$hi, NULL_X), na.rm = TRUE)
    pad <- 0.04 * diff(rng)
    xlim <- c(rng[1] - pad, rng[2] + pad)
    wide <- diff(rng)
    span <- diff(range(L$est, na.rm = TRUE))
    if (is.finite(span) && span > 0 && wide > 3 * span)
      message("Note: the widest interval is more than three times the spread ",
              "of the point estimates, so the axis is dominated by one row. ",
              "Set X_LIM to clamp it (clamped intervals get an arrowhead).")
  } else {
    xlim <- X_LIM
  }
  xbreaks <- if (is.null(X_BREAKS)) pretty(xlim, n = 6) else X_BREAKS
  xbreaks <- xbreaks[xbreaks > xlim[1] & xbreaks < xlim[2]]

  n_lines <- nrow(L)
  height  <- HEIGHT_BASE + HEIGHT_PER_ROW * n_lines

  combined <- patchwork::wrap_plots(
    panel_left(L), panel_forest(L, xlim, xbreaks), panel_right(L),
    nrow = 1, widths = PANEL_WIDTHS) +
    patchwork::plot_annotation(
      title = paste0("Motor FIM non-inferiority, margin ", M_TXT,
                     " points, in the subset-defining strata: ",
                     TITLE_TAIL[[rf]]),
      subtitle = wrap_text(paste0(
        P_LABEL, ", shown with its 95% confidence interval, is P(A + ", M_TXT,
        " > B) + 0.5 P(A + ", M_TXT, " = B) for a patient A from the subgroup ",
        "of interest and a patient B from its comparator. ",
        COMPARATOR_NOTE[[rf]], " The dashed line at 0.5 is the non-inferiority ",
        "boundary, where the subgroup falls short by the whole margin of ",
        M_TXT, " points; an interval lying entirely to its right shows ",
        "non-inferiority at one-sided alpha ", format(ALPHA_ONE_SIDED), ". ",
        "HL shift is the Hodges-Lehmann estimate of the subgroup-minus-comparator ",
        "difference in motor FIM points - the median of all pairwise ",
        "differences, not the difference between the two medians on the left. ",
        "Diamonds are the whole-cohort row of each block, squares the disease ",
        "classes. ",
        "MSD, ", DISEASE_LONG[["MSD"]], "; CVD, ", DISEASE_LONG[["CVD"]],
        "; DS, ", DISEASE_LONG[["DS"]], ".")),
      caption = wrap_text(paste0(
        "Brunner-Munzel test for non-inferiority, one-sided, t approximation; ",
        "confidence intervals not adjusted for multiplicity; complete cases. ")),
      theme = theme(
        plot.title    = element_text(size = FS_BASE + 2, face = "bold",
                                     colour = COL_TEXT),
        plot.subtitle = element_text(size = FS_BASE - 1, colour = COL_TEXT),
        plot.caption  = element_text(size = FS_BASE - 2, colour = COL_TEXT,
                                     hjust = 0),
        plot.margin   = margin(8, 10, 6, 10)))

  save_tiff(FILE_FOR[[rf]], combined, width = FIG_WIDTH, height = height)

  L[, reference := rf]
  L[, x_lim_lo := xlim[1]]
  L[, x_lim_hi := xlim[2]]
  drawn[[rf]] <- L
}

rows_out <- rbindlist(drawn, fill = TRUE)
setcolorder(rows_out, c("reference", "ord", "kind", "label",
                        "n_a", "med_a", "n_b", "med_b", "hl", "p_ci",
                        "est", "lo", "hi", "is_all", "tested"))
fwrite(rows_out, F_ROWS)

message("\n---- 22-1-1  Rows drawn ----")
print(rows_out[kind == "data",
               .(reference, label = trimws(label), n_a, med_a, n_b, med_b, hl, p_ci)],
      nrows = 80)

message("\n==== Section 09 complete ====")
message("  ", F_WITHIN)
message("  ", F_COHORT)
message("  ", F_ROWS)
message("")
message("The x axis is ", P_LABEL, ", a probability, so the reference line is at ",
        "0.5: an interval lying entirely to its right shows non-inferiority. ",
        "HL shift is printed without its interval, and it is not the difference ",
        "between the two medians in the left-hand columns.")
message("Read `reference` before quoting a row: the two figures answer ",
        "different questions (EuGM draft of 2026-09-10: within-disease = Fig. 3, ",
        "cohort-remainder = Online Resource 9).")

library(data.table)
library(here)
library(ggplot2)

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
for (d in c(OUT_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

THR_GOOD   <- 65
THR_SEVERE <- 26
THR_ELDER  <- 90

N_EXPECTED <- 2400

ALLOW_OVERWRITE <- FALSE

SEED    <- 20260909
JIT_A_X <- 0.16
JIT_A_Y <- 0.35
JIT_B_X <- 0.32
JIT_B_Y <- 0.32

VIOLIN_SCALE <- "width"
SHOW_BOX     <- TRUE
SHOW_N       <- TRUE
STRIP_N      <- TRUE

GREY_OLD   <- "grey76"
GREY_YOUNG <- "grey94"
COL_LINE   <- "black"
COL_PT     <- "black"
COL_REF    <- "black"
COL_TEXT   <- "grey15"
RAMP_LO    <- "grey80"
RAMP_HI    <- "grey5"

DISEASE_ORDER <- c("MSD", "CVD", "DS")
DISEASE_LONG  <- c(MSD = "musculoskeletal disease",
                   CVD = "cerebrovascular disease",
                   DS  = "disuse syndrome")
ALL_LAB       <- "All"

AGE_LAB_OLD   <- sprintf("Oldest-old\n(age >= %d)", THR_ELDER)
AGE_LAB_YOUNG <- sprintf("Others\n(age < %d)",      THR_ELDER)
AGE_ORDER     <- c(AGE_LAB_OLD, AGE_LAB_YOUNG)

FIM_BREAKS <- c(13, THR_SEVERE, 50, 75, 91)

set.seed(SEED)

wrap_text <- function(x, width = 108) paste(strwrap(x, width = width), collapse = "\n")

if (!exists("relabel_levels")) {
  cand <- here::here("01_labels.R")
  hit <- cand[file.exists(cand)]
  if (length(hit)) {
    source(hit[1])
    message("Sourced labels from: ", hit[1])
  } else {
    warning("01_labels.R not found in: ", paste(cand, collapse = ", "),
            " - falling back to the built-in disease labels.")
    LEVEL_LABELS_FALLBACK <- c("脳血管" = "CVD",
                               "運動器" = "MSD",
                               "廃用"   = "DS")
    relabel_levels <- function(var, lv) {
      unname(ifelse(lv %in% names(LEVEL_LABELS_FALLBACK),
                    LEVEL_LABELS_FALLBACK[lv], lv))
    }
  }
}

F_VIOLIN  <- file.path(FIG_DIR, "fig_elderly_mfim_violin.tiff")
F_SCATTER <- file.path(FIG_DIR, "fig_age_mfim_discharge_scatter.tiff")
F_SUM_A   <- file.path(OUT_DIR, "fig_elderly_mfim_violin_summary.csv")
F_SUM_B   <- file.path(OUT_DIR, "fig_age_mfim_discharge_quadrants.csv")
F_QC      <- file.path(OUT_DIR, "qc_disease_class_labelling_sec07.csv")

targets <- c(F_VIOLIN, F_SCATTER, F_SUM_A, F_SUM_B, F_QC)
exists_already <- targets[file.exists(targets)]
if (length(exists_already) && !ALLOW_OVERWRITE) {
  stop("These outputs already exist and ALLOW_OVERWRITE is FALSE:\n  ",
       paste(exists_already, collapse = "\n  "),
       "\nDelete them, rename them, or set ALLOW_OVERWRITE <- TRUE.")
}

theme_grey_paper <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor  = element_blank(),
          panel.grid.major  = element_line(linewidth = 0.25, colour = "grey90"),
          panel.border      = element_rect(colour = "grey35", fill = NA, linewidth = 0.4),
          strip.background  = element_rect(fill = "grey92", colour = NA),
          strip.text        = element_text(face = "bold", size = base_size - 0.5,
                                           colour = COL_TEXT),
          axis.text         = element_text(colour = COL_TEXT),
          axis.title        = element_text(colour = COL_TEXT),
          plot.title        = element_text(face = "bold", size = base_size + 1,
                                           colour = COL_TEXT),
          plot.subtitle     = element_text(size = base_size - 1.5, colour = COL_TEXT),
          plot.caption      = element_text(size = base_size - 2, colour = COL_TEXT,
                                           hjust = 0),
          plot.title.position = "plot")
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

if (!exists("BNB_small")) stop("BNB_small not found. Load it first.")
dat <- as.data.table(BNB_small)

need <- c("age", "class", "mFIM_in", "mFIM_out")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("BNB_small has no column(s): ", paste(miss, collapse = ", "))

dat[, age      := as.numeric(age)]
dat[, mFIM_in  := as.numeric(mFIM_in)]
dat[, mFIM_out := as.numeric(mFIM_out)]
dat[, class    := factor(class)]
dat[, elderly  := age      >= THR_ELDER]
dat[, good     := mFIM_out >= THR_GOOD]

if (nrow(dat) != N_EXPECTED)
  warning(sprintf("BNB_small has %d rows, not the %d analysed admissions.",
                  nrow(dat), N_EXPECTED))

n_before <- nrow(dat)
dat <- dat[is.finite(age) & is.finite(mFIM_in) & is.finite(mFIM_out) & !is.na(class)]
if (nrow(dat) != n_before)
  warning(sprintf("%d row(s) dropped for a missing age, class, mFIM_in or mFIM_out.",
                  n_before - nrow(dat)))

lv_raw <- levels(dat$class)
lv_eng <- relabel_levels("class", lv_raw)
if (anyNA(lv_eng) || any(lv_eng == lv_raw))
  warning("Some disease levels were not translated by relabel_levels(): ",
          paste(lv_raw[lv_eng == lv_raw], collapse = ", "))
dat[, class_eng := factor(lv_eng[match(as.character(class), lv_raw)],
                          levels = intersect(DISEASE_ORDER, lv_eng))]

qc <- dat[, .(n = .N), by = .(raw = as.character(class))]
qc[, english          := lv_eng[match(raw, lv_raw)]]
qc[, percent_of_cohort := round(100 * n / nrow(dat), 1)]
qc[, level_position   := match(raw, lv_raw)]
qc[, table1_expected  := c("脳血管" = 45.5,
                           "運動器" = 40.3,
                           "廃用"   = 14.1)[raw]]
setcolorder(qc, c("level_position", "raw", "english", "n",
                  "percent_of_cohort", "table1_expected"))
setorder(qc, level_position)
fwrite(qc, F_QC)

message("---- 22.0 Disease-class labelling, reconciled against Table 1 ----")
print(qc)
bad <- qc[is.finite(table1_expected) & abs(percent_of_cohort - table1_expected) > 0.5]
if (nrow(bad))
  warning("The class percentages do not match Table 1. Do not use these figures ",
          "until this is resolved: ", paste(bad$raw, collapse = ", "))

to_long <- function(d) {
  out <- rbindlist(list(copy(d)[, disease_facet := ALL_LAB],
                        copy(d)[, disease_facet := as.character(class_eng)]),
                   fill = TRUE)
  out <- out[!is.na(disease_facet)]
  out[, disease_facet := factor(disease_facet, levels = c(ALL_LAB, DISEASE_ORDER))]
  out[]
}

panel_n <- function(d) d[, .(n = .N), by = disease_facet]

strip_labeller <- function(d) {
  if (!STRIP_N) return(ggplot2::label_value)
  pn  <- panel_n(d)
  map <- setNames(sprintf("%s (n = %s)", as.character(pn$disease_facet),
                          format(pn$n, big.mark = ",", trim = TRUE)),
                  as.character(pn$disease_facet))
  ggplot2::as_labeller(map)
}

figA <- to_long(dat)
figA[, age_group := factor(fifelse(elderly, AGE_LAB_OLD, AGE_LAB_YOUNG),
                           levels = AGE_ORDER)]

sumA <- figA[, .(n        = .N,
                 median   = as.numeric(median(mFIM_in)),
                 q1       = as.numeric(quantile(mFIM_in, 0.25, type = 7)),
                 q3       = as.numeric(quantile(mFIM_in, 0.75, type = 7)),
                 mean     = mean(mFIM_in),
                 sd       = stats::sd(mFIM_in),
                 min      = min(mFIM_in),
                 max      = max(mFIM_in),
                 n_severe = sum(mFIM_in <= THR_SEVERE),
                 pct_severe = 100 * mean(mFIM_in <= THR_SEVERE)),
             by = .(disease_facet, age_group)]
setorder(sumA, disease_facet, age_group)
fwrite(sumA[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], F_SUM_A)
message("\n---- 22.1 Admission motor FIM by panel and age stratum ----")
print(sumA[, .(disease_facet, age_group = sub("\n", " ", age_group), n,
               median, q1, q3, pct_severe = round(pct_severe, 1))])

y_lo_A <- min(figA$mFIM_in); y_hi_A <- max(figA$mFIM_in)
pad_A  <- 0.10 * (y_hi_A - y_lo_A)
n_lab_A <- figA[, .(n = .N), by = .(disease_facet, age_group)]
n_lab_A[, lab := sprintf("n = %s", format(n, big.mark = ",", trim = TRUE))]

pA <- ggplot(figA, aes(x = age_group, y = mFIM_in)) +
  geom_violin(aes(fill = age_group), trim = TRUE, scale = VIOLIN_SCALE,
              width = 0.88, colour = COL_LINE, linewidth = 0.35) +
  geom_point(position = position_jitter(width = JIT_A_X, height = JIT_A_Y,
                                        seed = SEED),
             shape = 16, size = 0.42, alpha = 0.32, colour = COL_PT) +
  scale_fill_manual(values = setNames(c(GREY_OLD, GREY_YOUNG), AGE_ORDER),
                    guide = "none") +
  scale_y_continuous(breaks = FIM_BREAKS,
                     limits = c(y_lo_A - pad_A, y_hi_A + 0.02 * (y_hi_A - y_lo_A)),
                     expand = c(0, 0)) +
  facet_grid(. ~ disease_facet, labeller = strip_labeller(figA)) +
  labs(x = NULL, y = "Motor FIM at admission",
       title = "Admission motor FIM in the oldest-old and in the rest of the cohort",
       subtitle = wrap_text(paste0(
         "Each panel shows every admission in that disease class as one point, ",
         "overlaid on the violin of its distribution. The white box is the ",
         "median and the interquartile range. MSD, musculoskeletal disease; ",
         "CVD, cerebrovascular disease; DS, disuse syndrome.")),
       caption = wrap_text(paste0(
         "Violins are scaled to a common maximum width, so the two strata are ",
         "comparable in shape and not in area. Points are jittered by up to ",
         JIT_A_Y, " motor FIM points to show the local density."))) +
  theme_grey_paper() +
  theme(axis.text.x   = element_text(size = 9, lineheight = 0.95),
        panel.spacing = grid::unit(0.6, "lines"))

if (SHOW_BOX)
  pA <- pA + geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white",
                          colour = COL_LINE, linewidth = 0.32, coef = 0)

if (SHOW_N)
  pA <- pA + geom_text(data = n_lab_A, aes(x = age_group, y = y_lo_A - pad_A,
                                           label = lab),
                       inherit.aes = FALSE, vjust = -0.25, size = 2.6,
                       colour = COL_TEXT)

save_tiff(F_VIOLIN, pA, width = 10.0, height = 5.0)

figB <- to_long(dat)
figB[, good_f := factor(fifelse(good, "TRUE", "FALSE"), levels = c("FALSE", "TRUE"))]

setorder(figB, disease_facet, good_f)

quad <- figB[, {
  n <- .N
  .(n = n,
    pct_young_mild   = 100 * mean(age <  THR_ELDER & mFIM_in >  THR_SEVERE),
    pct_young_severe = 100 * mean(age <  THR_ELDER & mFIM_in <= THR_SEVERE),
    pct_old_mild     = 100 * mean(age >= THR_ELDER & mFIM_in >  THR_SEVERE),
    pct_old_severe   = 100 * mean(age >= THR_ELDER & mFIM_in <= THR_SEVERE),
    pct_good         = 100 * mean(good))
}, by = disease_facet]
fwrite(quad[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)], F_SUM_B)
message("\n---- 22.2 Share of each panel in the four quadrants of Figure B ----")
print(quad[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])

x_lo <- floor(min(figB$age)   / 10) * 10
x_hi <- ceiling(max(figB$age) / 10) * 10
y_lo <- min(figB$mFIM_in); y_hi <- max(figB$mFIM_in)

x_breaks <- seq(x_lo, x_hi, by = 10)
x_breaks <- x_breaks[x_breaks > x_lo & x_breaks < x_hi]

ramp_breaks <- pretty(range(figB$mFIM_out), n = 4)
ramp_breaks <- ramp_breaks[abs(ramp_breaks - THR_GOOD) > 8]
ramp_breaks <- sort(unique(c(ramp_breaks, THR_GOOD)))
ramp_breaks <- ramp_breaks[ramp_breaks >= min(figB$mFIM_out) &
                           ramp_breaks <= max(figB$mFIM_out)]

pB <- ggplot(figB, aes(x = age, y = mFIM_in)) +
  geom_point(aes(colour = mFIM_out, shape = good_f, size = good_f),
             position = position_jitter(width = JIT_B_X, height = JIT_B_Y,
                                        seed = SEED),
             alpha = 0.80) +
  geom_hline(yintercept = THR_SEVERE, linetype = "dashed", linewidth = 0.45,
             colour = COL_REF) +
  geom_vline(xintercept = THR_ELDER,  linetype = "dashed", linewidth = 0.45,
             colour = COL_REF) +
  scale_colour_gradient(low = RAMP_LO, high = RAMP_HI,
                        breaks = ramp_breaks, limits = range(figB$mFIM_out),
                        name = "Motor FIM at discharge (shading)",
                        guide = guide_colourbar(barheight = grid::unit(0.32, "cm"),
                                                barwidth  = grid::unit(3.6, "cm"),
                                                title.position = "top",
                                                order = 1)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 15),
                     labels = c("FALSE" = sprintf("< %d", THR_GOOD),
                                "TRUE"  = sprintf(">= %d", THR_GOOD)),
                     name = "Motor FIM at discharge (symbol)",
                     guide = guide_legend(title.position = "top", order = 2,
                                          override.aes = list(size = 2.2,
                                                              alpha = 1,
                                                              colour = RAMP_HI))) +
  scale_size_manual(values = c("FALSE" = 0.85, "TRUE" = 0.70), guide = "none") +
  scale_x_continuous(breaks = x_breaks,
                     limits = c(x_lo, x_hi), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = FIM_BREAKS,
                     limits = c(y_lo - 1, y_hi + 1), expand = c(0, 0)) +
  facet_grid(. ~ disease_facet, labeller = strip_labeller(figB)) +
  labs(x = "Age at admission (years)", y = "Motor FIM at admission",
       title = "Admission motor FIM against age, shaded by the discharge motor FIM",
       subtitle = wrap_text(paste0(
         "Dashed lines: the severity threshold (admission motor FIM ", THR_SEVERE,
         ") and age ", THR_ELDER, ". MSD, musculoskeletal disease; ",
         "CVD, cerebrovascular disease; DS, disuse syndrome.")),
       caption = wrap_text(paste0(
         "One point per admission, jittered by up to ", JIT_B_X, " years and ",
         JIT_B_Y, " motor FIM points. Points are shaded from light grey (lowest ",
         "discharge motor FIM) to black (highest); admissions reaching a ",
         "discharge motor FIM of ", THR_GOOD, " or more are drawn as squares."))) +
  theme_grey_paper() +
  theme(panel.spacing     = grid::unit(0.6, "lines"),
        legend.position   = "bottom",
        legend.box        = "horizontal",
        legend.title      = element_text(size = 9, colour = COL_TEXT),
        legend.text       = element_text(size = 8.5, colour = COL_TEXT),
        legend.key        = element_blank(),
        legend.background = element_blank())

save_tiff(F_SCATTER, pB, width = 11.0, height = 5.2)

message("\n==== Section 07 complete ====")
message("  ", F_VIOLIN)
message("  ", F_SCATTER)
message("  ", F_SUM_A)
message("  ", F_SUM_B)
message("  ", F_QC)
message("")
message("Check before use: qc_disease_class_labelling_sec07.csv must reproduce ",
        "Table 1 (CVD 45.5, MSD 40.3, DS 14.1). If it does not, the panels are ",
        "mislabelled and nothing drawn here should be used.")

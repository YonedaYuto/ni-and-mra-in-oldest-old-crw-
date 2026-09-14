library(data.table)
library(here)

OUT_DIR <- here::here("data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

THR_SEVERE <- 26
THR_ELDER  <- 90

REPORTED_MEDIAN <- c(Overall = 81, Severe = 99.5, Elderly = 81)
REPORTED_N      <- c(Overall = 2400, Severe = 576, Elderly = 337)

rds <- file.path(OUT_DIR, "BNB_smallplus.rds")
csv <- file.path(OUT_DIR, "BNB_smallplus.csv")

if (exists("BNB_smallplus")) {
  dat <- as.data.table(BNB_smallplus)
  src <- "the object BNB_smallplus in the session"
} else if (file.exists(rds)) {
  dat <- as.data.table(readRDS(rds))
  src <- rds
} else if (file.exists(csv)) {
  dat <- fread(csv)
  src <- csv
} else {
  stop("BNB_smallplus not found. Run script 03 first, or load it into the session.")
}
message("Read from: ", src)

need <- c("term", "age", "mFIM_in")
if (length(setdiff(need, names(dat))) > 0) {
  stop("BNB_smallplus is missing: ", paste(setdiff(need, names(dat)), collapse = ", "))
}

dat[, `:=`(term    = as.numeric(term),
           age     = as.numeric(age),
           mFIM_in = as.numeric(mFIM_in))]

if (nrow(dat) != REPORTED_N[["Overall"]]) {
  warning(sprintf("Expected %d analysed admissions but read %d rows. Check the input.",
                  REPORTED_N[["Overall"]], nrow(dat)))
}

strata <- list(
  "Whole analysis cohort"        = rep(TRUE, nrow(dat)),
  "Severe subset"                = dat$mFIM_in <= THR_SEVERE,
  "Oldest-old subset"            = dat$age     >= THR_ELDER,
  "Severe and oldest-old"        = dat$mFIM_in <= THR_SEVERE & dat$age >= THR_ELDER
)

person_time <- function(label, keep) {
  keep <- keep & !is.na(keep)
  d    <- dat[keep]
  n    <- nrow(d)
  n_na <- sum(is.na(d$term))
  v    <- d$term[!is.na(d$term)]
  q    <- if (length(v)) as.numeric(quantile(v, c(0.25, 0.5, 0.75), type = 7)) else rep(NA_real_, 3)
  data.table(
    stratum                 = label,
    n                       = n,
    n_missing_length_of_stay= n_na,
    n_contributing          = length(v),
    total_person_days       = sum(v),
    total_person_years      = round(sum(v) / 365.25, 1),
    mean_days               = round(mean(v), 1),
    sd_days                 = round(stats::sd(v), 1),
    median_days             = q[2],
    q1_days                 = q[1],
    q3_days                 = q[3],
    min_days                = min(v),
    max_days                = max(v)
  )
}

res <- rbindlist(Map(person_time, names(strata), strata))
print(res)
fwrite(res, file.path(OUT_DIR, "table_person_time.csv"))

check <- data.table(
  stratum        = c("Overall", "Severe", "Elderly"),
  n_recomputed   = res[c(1, 2, 3), n],
  n_reported     = as.numeric(REPORTED_N[c("Overall", "Severe", "Elderly")]),
  med_recomputed = res[c(1, 2, 3), median_days],
  med_reported   = as.numeric(REPORTED_MEDIAN[c("Overall", "Severe", "Elderly")])
)
check[, ok := n_recomputed == n_reported & abs(med_recomputed - med_reported) < 0.05]
print(check)
message(sprintf("Reconciliation with Table 1: %s",
                if (all(check$ok)) "OK" else "*** MISMATCH - resolve before reporting ***"))

o <- res[stratum == "Whole analysis cohort"]
if (o$n_missing_length_of_stay > 0) {
  message(sprintf("NOTE: %d patient(s) have a missing length of stay and are excluded from the total.",
                  o$n_missing_length_of_stay))
}

message("\n---- paste into Results, Participant flow and baseline characteristics ----\n")
message(sprintf(
  paste0("The observation period corresponds to the length of stay; the %s patients of the ",
         "analysis cohort contributed a total of %s patient-days of follow-up (%s patient-years; ",
         "mean %s days, standard deviation %s; median %s days, first to third quartile %s to %s)."),
  format(o$n, big.mark = ","),
  format(o$total_person_days, big.mark = ","),
  format(o$total_person_years, big.mark = ","),
  o$mean_days, o$sd_days, o$median_days, o$q1_days, o$q3_days))

message("\n---- optional, for the two subsets ----\n")
for (s in c("Severe subset", "Oldest-old subset")) {
  r <- res[stratum == s]
  message(sprintf("The %s (n = %d) contributed %s patient-days (median %s days, %s to %s).",
                  tolower(s), r$n, format(r$total_person_days, big.mark = ","),
                  r$median_days, r$q1_days, r$q3_days))
}

message("\nWritten: ", file.path(OUT_DIR, "table_person_time.csv"))

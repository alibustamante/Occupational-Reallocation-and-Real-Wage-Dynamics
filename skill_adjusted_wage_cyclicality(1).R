# ==============================================================================
# Multidimensional Skill-Adjusted Real Wage Cyclicality in the United States
# GitHub replication version | Final paper specification | 1982-2025
# ==============================================================================
#
# Replication code for the paper analysis. The script uses a disk-cached, low-memory
# workflow: it processes one CPS extract at a time and immediately collapses person
# records into compact occupation-month and demographic-cell-month panels. The full
# CPS microdata are never held in memory simultaneously.
#
# PRIMARY DESIGN
#   * Primary outcome: real weekly earnings, 1982-2025.
#   * Secondary outcome: real earnings per actual hour, beginning when valid
#     AHRSWORKT observations are available (normally 1989).
#   * Preferred adjustment: entropy calibration on eight O*NET latent factors.
#   * Transparent benchmark: fixed deciles of the equal-domain scalar skill score.
#   * Conventional worker-composition benchmark: fixed demographic cells based on
#     age group, education, sex, and race/ethnicity; industry-cell robustness is
#     optional and uses IND1990.
#   * Monthly cyclicality: 12-month changes on the 12-month unemployment change,
#     month fixed effects, and Newey-West standard errors.
#   * Quarterly corroboration: four-quarter changes on unemployment and the
#     four-quarter output-gap change.
#   * Episode analysis: clean pre-COVID, pandemic/recovery transition,
#     current-date and strict post-transition samples, corrected recession
#     influence-window exclusions, and factor-employment responses.
#
# INPUTS
#   ipums_cps_1981_1985.csv ... ipums_cps_2021_2025.csv
#   SOCXOCC2010.csv
#   onet_skill_measure_results_v2/data/onet_skill_measures_by_soc.csv
#   Unemployment Rate.csv, Output Gap.csv, Recessions.csv, CPIAUCSL.csv
#
# IMPORTANT
#   Run with source(), not by pasting expressions into the console:
#     source("skill_adjusted_wage_cyclicality.R", echo = FALSE)
#
#   DATA_DIR may be changed below to point to a private/local data directory. Raw
#   IPUMS CPS microdata should not be committed to a public GitHub repository.
# ==============================================================================

# ------------------------------- 0. CONFIG --------------------------------------

DATA_DIR <- Sys.getenv("SKILL_WAGE_DATA_DIR", unset = ".")
OUTPUT_DIR <- Sys.getenv(
  "SKILL_WAGE_OUTPUT_DIR",
  unset = file.path(DATA_DIR, "skill_wage_results")
)
# The public replication version maintains its own cache. With REBUILD_CACHE <- FALSE,
# existing cache files are reused and missing cache files are created automatically.
CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
CHUNK_DIR <- file.path(CACHE_DIR, "extract_panels")

CROSSWALK_FILE <- file.path(DATA_DIR, "SOCXOCC2010.csv")
SKILL_MEASURES_FILE <- file.path(
  DATA_DIR, "onet_skill_measure_results_v2", "data", "onet_skill_measures_by_soc.csv"
)
ALT_FACTOR_DIR <- file.path(DATA_DIR, "onet_skill_measure_results_v2", "data")
UNEMP_FILE <- file.path(DATA_DIR, "Unemployment Rate.csv")
OUTPUT_GAP_FILE <- file.path(DATA_DIR, "Output Gap.csv")
RECESSION_FILE <- file.path(DATA_DIR, "Recessions.csv")
CPI_FILE <- file.path(DATA_DIR, "CPIAUCSL.csv")

START_YEAR <- 1982L
END_YEAR <- 2025L
BASE_DOLLAR_YEAR <- 2025L
PRIMARY_BASELINE <- c(1982L, 1986L)
ALT_BASELINE <- c(1991L, 1995L)
MAIN_AGE_MIN <- 16L
MAIN_AGE_MAX <- 64L
ACTUAL_HOURS_START_YEAR <- 1989L
USUAL_ORG_HOURS_START_YEAR <- 1994L

PRIMARY_FACTOR_COLS <- sprintf("FACTOR_%02d", 1:8)
PRIMARY_FACTOR_LABELS <- c(
  "Physical and manual capability",
  "Management and interpersonal coordination",
  "Education and professional knowledge",
  "Engineering and technical systems",
  "Perceptual and cognitive vigilance",
  "Quantitative and scientific capability",
  "Information analysis and administrative execution",
  "Transportation and public-safety operations"
)
PRIMARY_SCALAR_SKILL <- "SCALAR_EQUAL_DOMAIN_Z"

# Pipeline controls.
REBUILD_CACHE <- FALSE         # Reuse this replication run's cache when available.
RUN_SECOND_PASS <- TRUE        # Needed for common-cap outcomes.
RUN_BLOCK_BOOTSTRAP <- TRUE
BOOTSTRAP_REPS <- 999L
BOOTSTRAP_BLOCK_MONTHS <- 24L
BOOTSTRAP_SEED <- 20260803L
RUN_ALT_FACTOR_COUNTS <- TRUE
RUN_SCALAR_GRID <- TRUE
RUN_BASELINE_SENSITIVITY <- TRUE
RUN_EPISODE_ANALYSIS <- TRUE
RUN_RECESSION_ANALYSIS <- TRUE
RUN_DEMOGRAPHIC_CELL_ADJUSTMENT <- TRUE
RUN_INDUSTRY_CELL_ROBUSTNESS <- FALSE  # Can create many cells; off by default.

# Wage support and upper-tail handling.
MIN_REAL_WEEKLY <- 20
MAX_REAL_WEEKLY <- 15000
MIN_REAL_HOURLY <- 2
MAX_REAL_HOURLY <- 500
MIN_REAL_ANNUAL <- 500
MAX_REAL_ANNUAL <- 750000
COMMON_CAP_QUANTILE <- 0.995

# Crosswalk coverage.
MIN_MATCH_RATE_HARD <- 0.80
MIN_MATCH_RATE_TARGET <- 0.95

# Entropy calibration.
CALIBRATION_RIDGE <- 1e-8
CALIBRATION_MAXIT <- 2000L
CALIBRATION_TOL <- 0.005
CALIBRATION_HARD_TOL <- 0.05
CALIBRATION_MIN_EFFECTIVE_OCCUPATIONS <- 25
CALIBRATION_HARD_MIN_EFFECTIVE_OCCUPATIONS <- 10
CALIBRATION_MAX_WEIGHT_RATIO <- 100
CALIBRATION_HARD_MAX_WEIGHT_RATIO <- 1000
MIN_VALID_CALIBRATION_SHARE <- 0.95

UNEMPLOYMENT_UNITS <- "auto"
PANDEMIC_TRANSITION_START <- as.Date("2020-03-01")
PANDEMIC_TRANSITION_END <- as.Date("2022-12-01")
POST_TRANSITION_START <- as.Date("2023-01-01")
# A 12-month change dated January 2024 is the first comparison for which both
# the current and lagged month fall strictly after December 2022.
STRICT_POST_TRANSITION_START <- as.Date("2024-01-01")
MIN_EPISODE_OBSERVATIONS <- 18L
RECESSION_INFLUENCE_PRE_MONTHS <- 12L
RECESSION_INFLUENCE_POST_MONTHS <- 12L
COVID_RECESSION_START <- as.Date("2020-03-01")
COVID_RECESSION_INFLUENCE_END <- PANDEMIC_TRANSITION_END
DATA_TABLE_THREADS <- 2L

# ----------------------------- 1. PACKAGES --------------------------------------

required_packages <- c("data.table", "fixest", "ggplot2", "lubridate", "zoo", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(data.table)
library(fixest)
library(ggplot2)
library(lubridate)
library(zoo)
library(scales)
setDTthreads(DATA_TABLE_THREADS)
setFixest_nthreads(DATA_TABLE_THREADS)
options(datatable.print.nrows = 50L, scipen = 999, datatable.alloccol = 256L)

for (d in c("data", "tables", "figures", "models", "logs")) {
  dir.create(file.path(OUTPUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(OUTPUT_DIR, "logs", "run_status.txt")
writeLines("INCOMPLETE", status_file)
log_file <- file.path(OUTPUT_DIR, "logs", "analysis_log.txt")
if (file.exists(log_file)) file.remove(log_file)
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_message("Paper replication analysis started.")

# ------------------------------- 2. HELPERS -------------------------------------

clean_names <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  toupper(gsub("^_|_$", "", x))
}

assert_columns <- function(dt, required, object_name) {
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop(object_name, " missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_sd_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(NA_real_)
  mu <- weighted_mean_safe(x[ok], w[ok])
  sqrt(sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok]))
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}

normalize_soc <- function(x) {
  digits <- gsub("[^0-9]", "", trimws(as.character(x)))
  digits[nchar(digits) == 6L] <- paste0(digits[nchar(digits) == 6L], "00")
  digits[nchar(digits) == 7L] <- paste0(digits[nchar(digits) == 7L], "0")
  out <- rep(NA_character_, length(digits))
  valid <- nchar(digits) == 8L
  out[valid] <- paste0(substr(digits[valid], 1, 6), ".", substr(digits[valid], 7, 8))
  out
}

parse_flexible_date <- function(raw, object_name = "series", allow_all_na = FALSE) {
  if (inherits(raw, "Date")) return(as.Date(raw))
  x <- trimws(as.character(raw))
  x[x %in% c("", "NA", "N/A", ".", "NULL")] <- NA_character_
  out <- as.Date(rep(NA_character_, length(x)))
  for (fmt in c("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d", "%Y/%m/%d", "%Y-%m", "%m/%Y")) {
    idx <- is.na(out) & !is.na(x)
    if (!any(idx)) break
    out[idx] <- suppressWarnings(as.Date(x[idx], format = fmt))
  }
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) out[idx] <- as.Date(suppressWarnings(parse_date_time(x[idx], orders = c("mdy", "ymd"), quiet = TRUE)))
  if (all(is.na(out)) && !allow_all_na) stop("Could not parse dates in ", object_name, call. = FALSE)
  out
}

parse_monthly_date <- function(dt, object_name = "series") {
  dc <- intersect(c("DATE", "OBSERVATION_DATE", "MONTH_DATE", "TIME", "PERIOD"), names(dt))
  if (length(dc)) return(as.Date(floor_date(parse_flexible_date(dt[[dc[1]]], object_name), "month")))
  if (all(c("YEAR", "MONTH") %in% names(dt))) {
    return(as.Date(sprintf("%04d-%02d-01", as.integer(dt$YEAR), as.integer(dt$MONTH))))
  }
  stop(object_name, " lacks date or YEAR/MONTH.", call. = FALSE)
}

parse_quarterly_date <- function(dt, object_name = "series") {
  dc <- intersect(c("DATE", "OBSERVATION_DATE", "QUARTER_DATE", "TIME", "PERIOD"), names(dt))
  if (length(dc)) {
    raw <- as.character(dt[[dc[1]]])
    parsed <- parse_flexible_date(raw, object_name, TRUE)
    idx <- is.na(parsed) & !is.na(raw)
    if (any(idx)) {
      q <- suppressWarnings(as.yearqtr(raw[idx], "%Y Q%q"))
      q2 <- suppressWarnings(as.yearqtr(raw[idx], "%YQ%q"))
      q[is.na(q)] <- q2[is.na(q)]
      parsed[idx] <- as.Date(q)
    }
    if (all(is.na(parsed))) stop("Could not parse quarterly dates in ", object_name, call. = FALSE)
    return(as.Date(as.yearqtr(parsed)))
  }
  if (all(c("YEAR", "QUARTER") %in% names(dt))) return(as.Date(as.yearqtr(paste(dt$YEAR, dt$QUARTER), "%Y %q")))
  stop(object_name, " lacks quarterly date.", call. = FALSE)
}

select_value_column <- function(dt, candidates, exclude = character()) {
  nm <- setdiff(names(dt), exclude)
  hit <- intersect(clean_names(candidates), nm)
  if (length(hit)) return(hit[1])
  nums <- nm[vapply(dt[, ..nm], is.numeric, logical(1))]
  nums <- setdiff(nums, c("YEAR", "MONTH", "QUARTER"))
  if (length(nums) == 1L) return(nums)
  stop("Could not identify value column in series.", call. = FALSE)
}

read_monthly_series <- function(path, candidates, output_name) {
  dt <- fread(path); setnames(dt, clean_names(names(dt)))
  dt[, DATE := parse_monthly_date(dt, basename(path))]
  vc <- select_value_column(dt, candidates, "DATE")
  out <- dt[, .(DATE, VALUE = suppressWarnings(as.numeric(get(vc))))][!is.na(DATE)]
  if (anyDuplicated(out$DATE)) stop("Duplicate monthly dates in ", basename(path), call. = FALSE)
  setnames(out, "VALUE", output_name); setorder(out, DATE); out
}

read_quarterly_series <- function(path, candidates, output_name) {
  dt <- fread(path); setnames(dt, clean_names(names(dt)))
  dt[, QUARTER_DATE := parse_quarterly_date(dt, basename(path))]
  vc <- select_value_column(dt, candidates, "QUARTER_DATE")
  out <- dt[, .(QUARTER_DATE, VALUE = suppressWarnings(as.numeric(get(vc))))][!is.na(QUARTER_DATE)]
  if (anyDuplicated(out$QUARTER_DATE)) stop("Duplicate quarters in ", basename(path), call. = FALSE)
  setnames(out, "VALUE", output_name); setorder(out, QUARTER_DATE); out
}

coef_row <- function(model, term, model_name, outcome, cycle_variable, sample = "main", type = "") {
  b <- coef(model); s <- fixest::se(model); p <- fixest::pvalue(model)
  data.table(
    model = model_name, specification_type = type, outcome = outcome,
    cycle_variable = cycle_variable, sample = sample, term = term,
    estimate = if (term %in% names(b)) unname(b[term]) else NA_real_,
    std_error = if (term %in% names(s)) unname(s[term]) else NA_real_,
    p_value = if (term %in% names(p)) unname(p[term]) else NA_real_,
    nobs = nobs(model), r2 = fitstat(model, "r2")[[1L]]
  )
}

save_plot <- function(p, name, width = 10, height = 6) {
  ggsave(file.path(OUTPUT_DIR, "figures", name), p, width = width, height = height, dpi = 320, bg = "white")
}

# -------------------------- 3. MACRO AND SKILLS --------------------------------

required_files <- c(CROSSWALK_FILE, SKILL_MEASURES_FILE, UNEMP_FILE, OUTPUT_GAP_FILE, RECESSION_FILE, CPI_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

unemployment <- read_monthly_series(UNEMP_FILE, c("UNEMPLOYMENT_RATE", "UNRATE", "RATE", "VALUE"), "UNEMP_RATE")
unemployment <- unemployment[is.finite(UNEMP_RATE)]
if (UNEMPLOYMENT_UNITS == "auto" && max(unemployment$UNEMP_RATE, na.rm = TRUE) <= 1.5) unemployment[, UNEMP_RATE := 100 * UNEMP_RATE]
if (UNEMPLOYMENT_UNITS == "proportion") unemployment[, UNEMP_RATE := 100 * UNEMP_RATE]

recessions <- read_monthly_series(RECESSION_FILE, c("RECESSION", "USREC", "FLAG", "VALUE"), "RECESSION")
recessions[, RECESSION := as.integer(RECESSION > 0)]
output_gap <- read_quarterly_series(OUTPUT_GAP_FILE, c("OUTPUT_GAP", "GAP", "VALUE"), "OUTPUT_GAP")
cpi <- read_monthly_series(CPI_FILE, c("CPIAUCSL", "CPI", "INDEX", "VALUE"), "CPI")[is.finite(CPI) & CPI > 0]
base_cpi <- mean(cpi[year(DATE) == BASE_DOLLAR_YEAR, CPI], na.rm = TRUE)
if (!is.finite(base_cpi)) stop("No CPI base-year observations.", call. = FALSE)
cpi[, DEFLATOR := base_cpi / CPI]
macro_monthly <- Reduce(function(x, y) merge(x, y, by = "DATE", all = TRUE), list(unemployment, recessions, cpi))
setorder(macro_monthly, DATE)
fwrite(macro_monthly, file.path(OUTPUT_DIR, "data", "macro_monthly_clean.csv"))
fwrite(output_gap, file.path(OUTPUT_DIR, "data", "output_gap_clean.csv"))

skill_soc <- fread(SKILL_MEASURES_FILE); setnames(skill_soc, clean_names(names(skill_soc)))
assert_columns(skill_soc, c("SOC_CODE", PRIMARY_SCALAR_SKILL, PRIMARY_FACTOR_COLS), "skill measures")
skill_soc[, SOC_CODE := normalize_soc(SOC_CODE)]
skill_soc <- unique(skill_soc[!is.na(SOC_CODE)], by = "SOC_CODE")

# Add candidate factor files when available.
if (RUN_ALT_FACTOR_COUNTS) {
  for (k in c(5L, 10L, 13L)) {
    f <- file.path(ALT_FACTOR_DIR, sprintf("candidate_factor_scores_f%d.csv", k))
    if (file.exists(f)) {
      a <- fread(f); setnames(a, clean_names(names(a))); a[, SOC_CODE := normalize_soc(SOC_CODE)]
      prefix <- sprintf("F%02d_", k)
      cols <- grep(paste0("^", prefix, "[0-9]{2}$"), names(a), value = TRUE)
      if (length(cols) == k) skill_soc <- merge(skill_soc, a[, c("SOC_CODE", cols), with = FALSE], by = "SOC_CODE", all.x = TRUE)
    }
  }
}

crosswalk <- fread(CROSSWALK_FILE); setnames(crosswalk, clean_names(names(crosswalk)))
assert_columns(crosswalk, "SOC_CODE", "crosswalk")
occ_col <- if ("OCC2010" %in% names(crosswalk)) "OCC2010" else if ("OCC" %in% names(crosswalk)) "OCC" else stop("Crosswalk needs OCC2010 or OCC.")
crosswalk[, `:=`(SOC_CODE = normalize_soc(SOC_CODE), OCC = as.integer(get(occ_col)))]
wc <- intersect(c("CROSSWALK_WEIGHT", "EMPLOYMENT_WEIGHT", "CONVERSION_FACTOR", "SHARE", "WEIGHT"), names(crosswalk))
if (length(wc)) {
  crosswalk[, MAP_WEIGHT_RAW := suppressWarnings(as.numeric(get(wc[1])))]
  weight_source <- wc[1]
} else {
  crosswalk[, MAP_WEIGHT_RAW := 1]
  weight_source <- "equal SOC weights"
}
crosswalk <- crosswalk[is.finite(OCC) & OCC > 0 & !is.na(SOC_CODE)]
if (weight_source == "equal SOC weights") {
  crosswalk <- unique(crosswalk[, .(SOC_CODE, OCC)])
  crosswalk[, MAP_WEIGHT_RAW := 1]
} else {
  crosswalk <- crosswalk[, .(MAP_WEIGHT_RAW = sum(MAP_WEIGHT_RAW[is.finite(MAP_WEIGHT_RAW) & MAP_WEIGHT_RAW > 0], na.rm = TRUE)), by = .(SOC_CODE, OCC)]
}
crosswalk[, MAP_WEIGHT := {
  w <- MAP_WEIGHT_RAW; w[!is.finite(w) | w <= 0] <- 0
  if (sum(w) <= 0) rep(1 / .N, .N) else w / sum(w)
}, by = OCC]

crosswalk_skill <- merge(crosswalk, skill_soc, by = "SOC_CODE", all.x = TRUE)
skill_numeric_cols <- names(crosswalk_skill)[vapply(crosswalk_skill, is.numeric, logical(1))]
skill_numeric_cols <- setdiff(skill_numeric_cols, c("OCC", "MAP_WEIGHT_RAW", "MAP_WEIGHT"))
weighted_map_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
occ_skill <- crosswalk_skill[, c(
  list(N_SOC_MAPPED = uniqueN(SOC_CODE)),
  lapply(.SD, weighted_map_mean, w = MAP_WEIGHT)
), by = OCC, .SDcols = skill_numeric_cols]
assert_columns(occ_skill, c("OCC", PRIMARY_SCALAR_SKILL, PRIMARY_FACTOR_COLS), "OCC skill table")
fwrite(occ_skill, file.path(OUTPUT_DIR, "data", "occ2010_skill_measures.csv"))
log_message("Prepared OCC2010 skills using ", weight_source, ".")

# Columns that are actually copied into compact panels.
alt_factor_cols <- grep("^(F05|F10|F13)_[0-9]{2}$", names(occ_skill), value = TRUE)
panel_skill_cols <- unique(c(PRIMARY_FACTOR_COLS, PRIMARY_SCALAR_SKILL, alt_factor_cols))

# ------------------------- 4. CPS FILE DISCOVERY --------------------------------

cps_files <- list.files(DATA_DIR, "^ipums_cps_[0-9]{4}_[0-9]{4}\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (!length(cps_files)) stop("No CPS files found.", call. = FALSE)
file_year_start <- as.integer(sub(".*_([0-9]{4})_[0-9]{4}\\.csv$", "\\1", cps_files))
cps_files <- cps_files[order(file_year_start)]
log_message("Found ", length(cps_files), " CPS extracts.")

core_vars <- c(
  "YEAR", "MONTH", "EARNWT", "AGE", "SEX", "RACE", "HISPAN", "EDUC",
  "STATEFIP", "OCC", "OCC2010", "IND", "IND1990", "EMPSTAT",
  "EARNWEEK2", "HOURWAGE2", "PAIDHOUR"
)
optional_vars <- c(
  "AHRSWORKT", "UHRSWORKORG", "WKSWORKORG", "UNION", "METFIPS", "METRO",
  "CLASSWKR", "MISH", "ELIGORG"
)

read_cps_extract <- function(path) {
  h <- fread(path, nrows = 0L, showProgress = FALSE); setnames(h, clean_names(names(h)))
  missing <- setdiff(core_vars, names(h))
  if (length(missing)) stop(basename(path), " missing core variables: ", paste(missing, collapse = ", "), call. = FALSE)
  keep <- intersect(c(core_vars, optional_vars), names(h))
  d <- fread(path, select = keep, showProgress = TRUE)
  setnames(d, clean_names(names(d)))
  for (v in setdiff(optional_vars, names(d))) d[, (v) := NA_real_]
  d[, SOURCE_FILE := basename(path)]
  d
}

# Recode conventional worker-composition cells. EDUC codes are grouped broadly so
# the 1992 educational-attainment redesign does not create a detailed-code break.
recode_worker_cells <- function(d) {
  d[, AGE_GROUP := cut(
    AGE, breaks = c(15, 24, 34, 44, 54, 64),
    labels = c("16-24", "25-34", "35-44", "45-54", "55-64"), right = TRUE
  )]
  d[, EDU4 := fifelse(
    EDUC < 73, "Less than high school",
    fifelse(EDUC < 81, "High school",
      fifelse(EDUC < 111, "Some college/associate", "Bachelor or more")
    )
  )]
  d[, RACE_ETH4 := fifelse(
    HISPAN > 0 & HISPAN < 900, "Hispanic",
    fifelse(RACE == 100, "Non-Hispanic White",
      fifelse(RACE == 200, "Non-Hispanic Black", "Non-Hispanic Other")
    )
  )]
  d[, SEX2 := fifelse(SEX == 1, "Male", "Female")]
  # Coarse IND1990 major groups by hundreds; retained only for optional robustness.
  d[, IND1990_MAJOR := fifelse(is.finite(IND1990) & IND1990 > 0, sprintf("%02d", as.integer(IND1990) %/% 100), NA_character_)]
  d
}

construct_wages <- function(d) {
  d[EARNWEEK2 >= 999999 | EARNWEEK2 <= 0, EARNWEEK2 := NA_real_]
  d[HOURWAGE2 >= 999 | HOURWAGE2 <= 0, HOURWAGE2 := NA_real_]
  d[AHRSWORKT <= 0 | AHRSWORKT > 168, AHRSWORKT := NA_real_]
  d[UHRSWORKORG <= 0 | UHRSWORKORG >= 997, UHRSWORKORG := NA_real_]
  d[WKSWORKORG <= 0 | WKSWORKORG >= 98, WKSWORKORG := NA_real_]
  d[YEAR < ACTUAL_HOURS_START_YEAR, AHRSWORKT := NA_real_]
  d[YEAR < USUAL_ORG_HOURS_START_YEAR, UHRSWORKORG := NA_real_]

  d[, WEEKLY_EARNINGS_REAL := EARNWEEK2 * DEFLATOR]
  d[, HOURLY_ACTUAL_REAL := EARNWEEK2 * DEFLATOR / AHRSWORKT]
  d[, HOURLY_USUAL_REAL := EARNWEEK2 * DEFLATOR / UHRSWORKORG]
  d[, HOURLY_REPORTED_REAL := fifelse(PAIDHOUR == 2, HOURWAGE2 * DEFLATOR, NA_real_)]
  d[, ANNUAL_EARNINGS_REAL := EARNWEEK2 * WKSWORKORG * DEFLATOR]

  d[WEEKLY_EARNINGS_REAL < MIN_REAL_WEEKLY | WEEKLY_EARNINGS_REAL > MAX_REAL_WEEKLY, WEEKLY_EARNINGS_REAL := NA_real_]
  for (v in c("HOURLY_ACTUAL_REAL", "HOURLY_USUAL_REAL", "HOURLY_REPORTED_REAL")) {
    d[get(v) < MIN_REAL_HOURLY | get(v) > MAX_REAL_HOURLY, (v) := NA_real_]
  }
  d[ANNUAL_EARNINGS_REAL < MIN_REAL_ANNUAL | ANNUAL_EARNINGS_REAL > MAX_REAL_ANNUAL, ANNUAL_EARNINGS_REAL := NA_real_]
  d
}

outcome_vars <- c(
  "WEEKLY_EARNINGS_REAL", "HOURLY_ACTUAL_REAL", "HOURLY_USUAL_REAL",
  "HOURLY_REPORTED_REAL", "ANNUAL_EARNINGS_REAL"
)

# -------------------- 5. PASS 1: CAPS AND AVAILABILITY --------------------------

pass1_cache <- file.path(CACHE_DIR, "pass1_annual_diagnostics.rds")
if (REBUILD_CACHE || !file.exists(pass1_cache)) {
  pass1_list <- vector("list", length(cps_files))
  availability_list <- vector("list", length(cps_files))
  file_availability <- vector("list", length(cps_files))

  for (i in seq_along(cps_files)) {
    path <- cps_files[i]
    log_message("Pass 1 reading ", basename(path), ".")
    d <- read_cps_extract(path)
    file_availability[[i]] <- data.table(
      SOURCE_FILE = basename(path), VARIABLE = c(core_vars, optional_vars),
      PRESENT = c(core_vars, optional_vars) %in% names(d)
    )
    d <- d[
      YEAR >= START_YEAR & YEAR <= END_YEAR & AGE >= MAIN_AGE_MIN & AGE <= MAIN_AGE_MAX &
        EMPSTAT %in% c(10, 12) & is.finite(EARNWT) & EARNWT > 0
    ]
    d[, `:=`(
      ANALYSIS_WEIGHT = as.numeric(EARNWT),
      DATE = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))
    )]
    d <- merge(d, cpi[, .(DATE, DEFLATOR)], by = "DATE", all.x = TRUE)
    d <- construct_wages(d)

    availability_list[[i]] <- d[, {
      valid_n <- lapply(.SD, function(x) sum(is.finite(x)))
      valid_w <- lapply(.SD, function(x) sum(ANALYSIS_WEIGHT[is.finite(x)]))
      names(valid_n) <- paste0("N_VALID_", names(valid_n))
      names(valid_w) <- paste0("WT_VALID_", names(valid_w))
      c(list(N = .N, WEIGHT = sum(ANALYSIS_WEIGHT)), valid_n, valid_w)
    }, by = YEAR, .SDcols = outcome_vars]

    pass1_list[[i]] <- rbindlist(lapply(outcome_vars, function(v) {
      d[is.finite(get(v)), .(
        ANNUAL_CAP = weighted_quantile(get(v), ANALYSIS_WEIGHT, COMMON_CAP_QUANTILE),
        P50 = weighted_quantile(get(v), ANALYSIS_WEIGHT, 0.50),
        P90 = weighted_quantile(get(v), ANALYSIS_WEIGHT, 0.90),
        P99 = weighted_quantile(get(v), ANALYSIS_WEIGHT, 0.99),
        MAX_VALUE = max(get(v)),
        N = .N,
        WEIGHT = sum(ANALYSIS_WEIGHT)
      ), by = YEAR][, OUTCOME := v]
    }), fill = TRUE)
    rm(d); gc(full = TRUE)
  }
  pass1 <- rbindlist(pass1_list, fill = TRUE)
  availability <- rbindlist(availability_list, fill = TRUE)
  file_availability_dt <- rbindlist(file_availability, fill = TRUE)
  saveRDS(list(pass1 = pass1, availability = availability, file_availability = file_availability_dt), pass1_cache, compress = FALSE)
} else {
  p1 <- readRDS(pass1_cache); pass1 <- p1$pass1; availability <- p1$availability; file_availability_dt <- p1$file_availability
}

common_caps <- pass1[is.finite(ANNUAL_CAP), .(COMMON_CAP = min(ANNUAL_CAP)), by = OUTCOME]
fwrite(pass1, file.path(OUTPUT_DIR, "tables", "annual_wage_distribution_diagnostics.csv"))
fwrite(common_caps, file.path(OUTPUT_DIR, "tables", "common_cap_diagnostics.csv"))
fwrite(availability, file.path(OUTPUT_DIR, "tables", "outcome_availability_by_year.csv"))
fwrite(file_availability_dt, file.path(OUTPUT_DIR, "tables", "file_variable_availability.csv"))

# Flag large year-to-year changes in observed upper support. These diagnostics do
# not alter the primary estimates, but they make historical topcode transitions
# visible alongside the common-cap robustness specification.
weekly_upper_support <- copy(pass1[OUTCOME == "WEEKLY_EARNINGS_REAL"])
setorder(weekly_upper_support, YEAR)
weekly_upper_support[, `:=`(
  MAX_VALUE_L1 = shift(MAX_VALUE),
  P99_L1 = shift(P99)
)]
weekly_upper_support[, `:=`(
  MAX_VALUE_CHANGE_PERCENT = 100 * (MAX_VALUE / MAX_VALUE_L1 - 1),
  P99_CHANGE_PERCENT = 100 * (P99 / P99_L1 - 1),
  LARGE_UPPER_SUPPORT_BREAK = (
    abs(MAX_VALUE_CHANGE_PERCENT) >= 25 |
      abs(P99_CHANGE_PERCENT) >= 25
  )
)]
fwrite(
  weekly_upper_support,
  file.path(OUTPUT_DIR, "tables", "weekly_earnings_upper_support_breaks.csv")
)

# ------------------- 6. PASS 2: EXTRACT-LEVEL PANELS ----------------------------

aggregate_occ_panel_wide <- function(d, sample_name) {
  sample_ok <- switch(sample_name,
    MAIN = rep(TRUE, nrow(d)),
    PRIME_AGE = d$AGE >= 25 & d$AGE <= 54,
    FULL_TIME = is.finite(d$AHRSWORKT) & d$AHRSWORKT >= 35,
    stop("Unknown sample.")
  )
  valid <- sample_ok & d$SKILL_MATCH & d$ANALYSIS_WEIGHT > 0
  if (!any(valid)) return(data.table())
  x <- d[valid]
  out <- x[, {
    ans <- list()
    for (v in outcome_vars) {
      ok <- is.finite(get(v))
      ans[[paste0("EW__", v)]] <- sum(ANALYSIS_WEIGHT[ok])
      ans[[paste0("WS__", v)]] <- sum(ANALYSIS_WEIGHT[ok] * get(v)[ok])
      ans[[paste0("N__", v)]] <- sum(ok)
      cap <- common_caps[OUTCOME == v, COMMON_CAP][1]
      if (RUN_SECOND_PASS && is.finite(cap)) {
        capped <- pmin(get(v), cap)
        okc <- is.finite(capped)
        vc <- paste0(v, "_COMMON_CAP")
        ans[[paste0("EW__", vc)]] <- sum(ANALYSIS_WEIGHT[okc])
        ans[[paste0("WS__", vc)]] <- sum(ANALYSIS_WEIGHT[okc] * capped[okc])
        ans[[paste0("N__", vc)]] <- sum(okc)
      }
    }
    ans <- c(ans, lapply(.SD, function(z) z[1L]))
    ans
  }, by = .(DATE, YEAR, MONTH, OCC), .SDcols = panel_skill_cols]
  out[, `:=`(QUARTER_DATE = as.Date(as.yearqtr(DATE)), SAMPLE = sample_name)]
  out
}

extract_outcome_panel <- function(wide_panel, outcome, sample) {
  ew <- paste0("EW__", outcome)
  ws <- paste0("WS__", outcome)
  nn <- paste0("N__", outcome)
  if (!all(c(ew, ws, nn) %in% names(wide_panel))) return(data.table())
  d <- wide_panel[SAMPLE == sample & is.finite(get(ew)) & get(ew) > 0]
  if (!nrow(d)) return(data.table())
  keep <- c("DATE", "YEAR", "MONTH", "QUARTER_DATE", "OCC", "SAMPLE", panel_skill_cols)
  out <- d[, ..keep]
  out[, `:=`(
    EMPLOYMENT_WEIGHT = d[[ew]],
    WAGE_SUM = d[[ws]],
    N = as.integer(d[[nn]]),
    MEAN_WAGE = d[[ws]] / d[[ew]],
    OUTCOME = outcome
  )]
  out
}

aggregate_demographic_panel <- function(d, outcome, include_industry = FALSE) {
  valid <- is.finite(d[[outcome]]) & d$SKILL_MATCH & d$ANALYSIS_WEIGHT > 0 &
    !is.na(d$AGE_GROUP) & !is.na(d$EDU4) & !is.na(d$SEX2) & !is.na(d$RACE_ETH4)
  if (include_industry) valid <- valid & !is.na(d$IND1990_MAJOR)
  if (!any(valid)) return(data.table())
  keys <- c("DATE", "YEAR", "MONTH", "AGE_GROUP", "EDU4", "SEX2", "RACE_ETH4")
  if (include_industry) keys <- c(keys, "IND1990_MAJOR")
  x <- d[valid]
  out <- x[, .(
    EMPLOYMENT_WEIGHT = sum(ANALYSIS_WEIGHT),
    WAGE_SUM = sum(ANALYSIS_WEIGHT * get(outcome)),
    N = .N
  ), by = keys]
  out[, `:=`(
    QUARTER_DATE = as.Date(as.yearqtr(DATE)),
    MEAN_WAGE = WAGE_SUM / EMPLOYMENT_WEIGHT,
    OUTCOME = outcome,
    CELL_TYPE = if (include_industry) "demographic_industry" else "demographic_core"
  )]
  out
}

chunk_manifest <- vector("list", length(cps_files))
if (REBUILD_CACHE) unlink(CHUNK_DIR, recursive = TRUE, force = TRUE)
dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(cps_files)) {
  path <- cps_files[i]
  cache_path <- file.path(CHUNK_DIR, paste0(tools::file_path_sans_ext(basename(path)), ".rds"))
  if (!REBUILD_CACHE && file.exists(cache_path)) {
    chunk_manifest[[i]] <- data.table(SOURCE_FILE = basename(path), CACHE_FILE = cache_path, REUSED = TRUE)
    next
  }
  log_message("Pass 2 building panels from ", basename(path), ".")
  d <- read_cps_extract(path)
  d <- d[
    YEAR >= START_YEAR & YEAR <= END_YEAR & AGE >= MAIN_AGE_MIN & AGE <= MAIN_AGE_MAX &
      EMPSTAT %in% c(10, 12) & is.finite(EARNWT) & EARNWT > 0
  ]
  d[, `:=`(
    ANALYSIS_WEIGHT = as.numeric(EARNWT),
    OCC = as.integer(OCC2010),
    DATE = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))
  )]
  d <- merge(d, cpi[, .(DATE, DEFLATOR)], by = "DATE", all.x = TRUE)
  d <- construct_wages(d)
  d <- recode_worker_cells(d)
  d <- merge(d, occ_skill[, c("OCC", panel_skill_cols), with = FALSE], by = "OCC", all.x = TRUE)
  d[, SKILL_MATCH := complete.cases(.SD), .SDcols = PRIMARY_FACTOR_COLS]

  coverage <- d[, .(
    N = .N, WEIGHT = sum(ANALYSIS_WEIGHT),
    MATCHED_N = sum(SKILL_MATCH), MATCHED_WEIGHT = sum(ANALYSIS_WEIGHT[SKILL_MATCH])
  ), by = YEAR]
  coverage[, `:=`(
    MATCH_RATE_N = MATCHED_N / N,
    MATCH_RATE_WEIGHT = MATCHED_WEIGHT / WEIGHT,
    SOURCE_FILE = basename(path)
  )]

  panels <- lapply(c("MAIN", "PRIME_AGE", "FULL_TIME"), function(sample_name) {
    aggregate_occ_panel_wide(d, sample_name)
  })
  occ_panel_chunk <- rbindlist(panels, fill = TRUE)

  demo_core <- if (RUN_DEMOGRAPHIC_CELL_ADJUSTMENT) aggregate_demographic_panel(d, "WEEKLY_EARNINGS_REAL", FALSE) else data.table()
  demo_ind <- if (RUN_INDUSTRY_CELL_ROBUSTNESS) aggregate_demographic_panel(d, "WEEKLY_EARNINGS_REAL", TRUE) else data.table()

  saveRDS(
    list(occ_panel = occ_panel_chunk, demographic_core = demo_core, demographic_industry = demo_ind, coverage = coverage),
    cache_path, compress = FALSE
  )
  chunk_manifest[[i]] <- data.table(SOURCE_FILE = basename(path), CACHE_FILE = cache_path, REUSED = FALSE)
  rm(d, occ_panel_chunk, demo_core, demo_ind, panels); gc(full = TRUE)
}

manifest <- rbindlist(chunk_manifest)
fwrite(manifest, file.path(OUTPUT_DIR, "tables", "cache_manifest.csv"))

# Load only compact panels.
chunk_objects <- lapply(manifest$CACHE_FILE, readRDS)
occ_panel_wide <- rbindlist(
  lapply(chunk_objects, function(x) x[["occ_panel"]]),
  use.names = TRUE,
  fill = TRUE
)
demo_core_panel <- rbindlist(
  lapply(chunk_objects, function(x) x[["demographic_core"]]),
  use.names = TRUE,
  fill = TRUE
)
demo_ind_panel <- rbindlist(
  lapply(chunk_objects, function(x) x[["demographic_industry"]]),
  use.names = TRUE,
  fill = TRUE
)
coverage_all <- rbindlist(
  lapply(chunk_objects, function(x) x[["coverage"]]),
  use.names = TRUE,
  fill = TRUE
)
coverage_by_year <- coverage_all[, .(
  N = sum(N), WEIGHT = sum(WEIGHT), MATCHED_N = sum(MATCHED_N), MATCHED_WEIGHT = sum(MATCHED_WEIGHT)
), by = YEAR]
coverage_by_year[, `:=`(MATCH_RATE_N = MATCHED_N / N, MATCH_RATE_WEIGHT = MATCHED_WEIGHT / WEIGHT)]
rm(chunk_objects); gc(full = TRUE)

fwrite(coverage_by_year, file.path(OUTPUT_DIR, "tables", "skill_crosswalk_coverage_by_year.csv"))
overall_match <- sum(coverage_by_year$MATCHED_WEIGHT) / sum(coverage_by_year$WEIGHT)
if (!is.finite(overall_match) || overall_match < MIN_MATCH_RATE_HARD) stop("Weighted skill coverage below hard threshold.", call. = FALSE)
if (min(coverage_by_year$MATCH_RATE_WEIGHT) < MIN_MATCH_RATE_TARGET) warning("At least one year has skill coverage below target.")

setorder(occ_panel_wide, SAMPLE, DATE, OCC)
saveRDS(occ_panel_wide, file.path(CACHE_DIR, "combined_occupation_month_panel_wide.rds"), compress = FALSE)
fwrite(occ_panel_wide[, .(DATE, YEAR, MONTH, QUARTER_DATE, OCC, SAMPLE)],
       file.path(OUTPUT_DIR, "data", "occupation_month_panel_keys.csv"))
log_message("Combined compact wide panel has ", format(nrow(occ_panel_wide), big.mark = ","), " rows.")

# ---------------------- 7. SCALAR DECOMPOSITION --------------------------------

assign_fixed_groups_panel <- function(panel, skill_var, n_groups, baseline, outcome, sample) {
  d <- panel[OUTCOME == outcome & SAMPLE == sample & is.finite(get(skill_var))]
  base <- d[YEAR >= baseline[1] & YEAR <= baseline[2], .(
    BASE_WEIGHT = sum(EMPLOYMENT_WEIGHT), SKILL_SCORE = weighted_mean_safe(get(skill_var), EMPLOYMENT_WEIGHT)
  ), by = OCC]
  if (nrow(base) < n_groups) stop("Too few occupations for scalar groups.")
  setorder(base, SKILL_SCORE, OCC)
  base[, CUM_SHARE := cumsum(BASE_WEIGHT) / sum(BASE_WEIGHT)]
  base[, GROUP := pmin(n_groups, pmax(1L, ceiling(CUM_SHARE * n_groups)))]
  cuts <- base[, .(UPPER = max(SKILL_SCORE)), by = GROUP][GROUP < n_groups, UPPER]
  all_occ <- unique(d[, .(OCC, SKILL_SCORE = get(skill_var))], by = "OCC")
  all_occ[, GROUP := pmin(n_groups, pmax(1L, 1L + findInterval(SKILL_SCORE, cuts, left.open = TRUE)))]
  all_occ
}

build_scalar_group_panel <- function(panel, skill_var, n_groups, baseline, outcome, sample, specification) {
  gm <- assign_fixed_groups_panel(panel, skill_var, n_groups, baseline, outcome, sample)
  d <- merge(panel[OUTCOME == outcome & SAMPLE == sample], gm[, .(OCC, GROUP)], by = "OCC", all.x = TRUE)
  out <- d[!is.na(GROUP), .(
    EMPLOYMENT_WEIGHT = sum(EMPLOYMENT_WEIGHT),
    MEAN_WAGE = sum(WAGE_SUM) / sum(EMPLOYMENT_WEIGHT),
    N = sum(N)
  ), by = .(DATE, YEAR, MONTH, QUARTER_DATE, GROUP)]
  out[, EMPLOYMENT_SHARE := EMPLOYMENT_WEIGHT / sum(EMPLOYMENT_WEIGHT), by = DATE]
  out[, `:=`(SPECIFICATION = specification, OUTCOME = outcome, SAMPLE = sample, N_GROUPS = n_groups)]
  out
}

exact_group_decomposition <- function(gp, time_var = "DATE", lag_periods = 12L) {
  ids <- c("GROUP", "SPECIFICATION", "OUTCOME", "SAMPLE", "N_GROUPS")
  lag <- gp[, c(ids, time_var, "EMPLOYMENT_SHARE", "MEAN_WAGE"), with = FALSE]
  setnames(lag, c("EMPLOYMENT_SHARE", "MEAN_WAGE"), c("SHARE_LAG", "WAGE_LAG"))
  if (time_var == "DATE") lag[, (time_var) := get(time_var) %m+% months(lag_periods)]
  else lag[, (time_var) := as.Date(as.yearqtr(get(time_var)) + lag_periods / 4)]
  x <- merge(gp, lag, by = c(ids, time_var), all.x = TRUE)
  x[, `:=`(
    WITHIN = 0.5 * (EMPLOYMENT_SHARE + SHARE_LAG) * (MEAN_WAGE - WAGE_LAG),
    COMPOSITION = 0.5 * (MEAN_WAGE + WAGE_LAG) * (EMPLOYMENT_SHARE - SHARE_LAG),
    CURRENT = EMPLOYMENT_SHARE * MEAN_WAGE,
    PREVIOUS = SHARE_LAG * WAGE_LAG
  )]
  bycols <- c(time_var, "SPECIFICATION", "OUTCOME", "SAMPLE", "N_GROUPS")
  out <- x[, .(
    OBSERVED_WAGE = sum(CURRENT, na.rm = TRUE),
    OBSERVED_WAGE_LAG = sum(PREVIOUS, na.rm = TRUE),
    WITHIN_CHANGE = sum(WITHIN, na.rm = TRUE),
    COMPOSITION_CHANGE = sum(COMPOSITION, na.rm = TRUE),
    VALID_GROUPS = sum(is.finite(WITHIN) & is.finite(COMPOSITION))
  ), by = bycols]
  out[, `:=`(
    TOTAL_CHANGE = OBSERVED_WAGE - OBSERVED_WAGE_LAG,
    TOTAL_PCT = 100 * (OBSERVED_WAGE / OBSERVED_WAGE_LAG - 1),
    WITHIN_PCT = 100 * WITHIN_CHANGE / OBSERVED_WAGE_LAG,
    COMPOSITION_PCT = 100 * COMPOSITION_CHANGE / OBSERVED_WAGE_LAG,
    DECOMP_ERROR = (OBSERVED_WAGE - OBSERVED_WAGE_LAG) - WITHIN_CHANGE - COMPOSITION_CHANGE,
    METHOD = "fixed scalar groups"
  )]
  out[VALID_GROUPS < N_GROUPS, c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT", "DECOMP_ERROR") := NA_real_]
  setorderv(out, time_var); out
}

quarterly_group_panel <- function(gp) {
  q <- gp[, .(
    EMPLOYMENT_WEIGHT = sum(EMPLOYMENT_WEIGHT),
    MEAN_WAGE = weighted_mean_safe(MEAN_WAGE, EMPLOYMENT_WEIGHT),
    N = sum(N)
  ), by = .(QUARTER_DATE, GROUP, SPECIFICATION, OUTCOME, SAMPLE, N_GROUPS)]
  q[, EMPLOYMENT_SHARE := EMPLOYMENT_WEIGHT / sum(EMPLOYMENT_WEIGHT), by = .(QUARTER_DATE, SPECIFICATION, OUTCOME, SAMPLE)]
  q
}

# -------------------- 8. CONTINUOUS CALIBRATION --------------------------------

build_moment_matrix <- function(x, order = 2L, cross_products = TRUE) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  parts <- list(x); names_out <- paste0("MEAN_", colnames(x))
  if (order >= 2L) {
    sq <- x^2; colnames(sq) <- paste0("SQ_", colnames(x)); parts[[2]] <- sq; names_out <- c(names_out, colnames(sq))
    if (cross_products && ncol(x) > 1L) {
      pairs <- combn(seq_len(ncol(x)), 2L)
      cr <- vapply(seq_len(ncol(pairs)), function(k) x[, pairs[1, k]] * x[, pairs[2, k]], numeric(nrow(x)))
      if (is.vector(cr)) cr <- matrix(cr, ncol = 1L)
      colnames(cr) <- vapply(seq_len(ncol(pairs)), function(k) paste0("CROSS_", colnames(x)[pairs[1, k]], "_", colnames(x)[pairs[2, k]]), character(1))
      parts[[length(parts) + 1L]] <- cr; names_out <- c(names_out, colnames(cr))
    }
  }
  out <- do.call(cbind, parts); colnames(out) <- names_out; out
}

entropy_balance_dual <- function(q, moments, target) {
  q <- as.numeric(q); x <- as.matrix(moments); target <- as.numeric(target)
  ok <- is.finite(q) & q > 0 & complete.cases(x)
  full <- rep(NA_real_, length(q)); q2 <- q[ok]; x2 <- x[ok, , drop = FALSE]; q2 <- q2 / sum(q2)
  if (nrow(x2) <= ncol(x2) + 2L) return(list(weights = full, convergence = 99L, max_error = Inf, effective_n = NA, ratio = NA))
  fn <- function(lambda) {
    eta <- as.vector(x2 %*% lambda); m <- max(eta)
    m + log(sum(q2 * exp(eta - m))) - sum(lambda * target) + 0.5 * CALIBRATION_RIDGE * sum(lambda^2)
  }
  gr <- function(lambda) {
    eta <- as.vector(x2 %*% lambda); eta <- eta - max(eta)
    p <- q2 * exp(eta); p <- p / sum(p)
    as.vector(crossprod(p, x2)) - target + CALIBRATION_RIDGE * lambda
  }
  fit <- tryCatch(optim(rep(0, ncol(x2)), fn, gr, method = "BFGS", control = list(maxit = CALIBRATION_MAXIT, reltol = 1e-10)), error = function(e) NULL)
  if (is.null(fit)) return(list(weights = full, convergence = 98L, max_error = Inf, effective_n = NA, ratio = NA))
  eta <- as.vector(x2 %*% fit$par); eta <- eta - max(eta)
  p <- q2 * exp(eta); p <- p / sum(p); full[ok] <- p
  achieved <- as.vector(crossprod(p, x2))
  list(weights = full, convergence = fit$convergence, max_error = max(abs(achieved - target)), effective_n = 1 / sum(p^2), ratio = max(p / q2))
}

calibrate_month <- function(q, factor_matrix, targets) {
  attempts <- list(
    full = list(x = build_moment_matrix(factor_matrix, 2, TRUE), target = targets$full),
    means_squares = list(x = build_moment_matrix(factor_matrix, 2, FALSE), target = targets$means_squares),
    means = list(x = build_moment_matrix(factor_matrix, 1, FALSE), target = targets$means)
  )
  last <- NULL
  for (nm in names(attempts)) {
    a <- attempts[[nm]]; fit <- entropy_balance_dual(q, a$x, a$target); last <- fit
    accepted <- fit$convergence == 0L && is.finite(fit$max_error) && fit$max_error <= CALIBRATION_TOL &&
      is.finite(fit$effective_n) && fit$effective_n >= CALIBRATION_MIN_EFFECTIVE_OCCUPATIONS &&
      is.finite(fit$ratio) && fit$ratio <= CALIBRATION_MAX_WEIGHT_RATIO
    if (accepted) { fit$moment_set <- nm; fit$accepted <- TRUE; return(fit) }
  }
  last$moment_set <- "means_failed_thresholds"; last$accepted <- FALSE; last
}

choose_baseline <- function(panel, requested) {
  years <- sort(unique(panel[is.finite(MEAN_WAGE), YEAR]))
  use <- years[years >= requested[1] & years <= requested[2]]
  if (length(use) >= 3L) return(c(min(use), max(use)))
  if (length(years) < 3L) stop("Insufficient years for outcome baseline.")
  c(years[1], years[min(5L, length(years))])
}

build_continuous_index_panel <- function(panel, outcome, sample, factor_cols, baseline, specification) {
  d <- copy(panel[OUTCOME == outcome & SAMPLE == sample])
  d <- d[complete.cases(d[, ..factor_cols])]
  if (!nrow(d)) stop("No compact panel rows for ", outcome, " / ", sample)
  baseline <- choose_baseline(d, baseline)
  base <- d[YEAR >= baseline[1] & YEAR <= baseline[2], c(
    list(BASE_WEIGHT = sum(EMPLOYMENT_WEIGHT)), lapply(.SD, function(z) z[1L])
  ), by = OCC, .SDcols = factor_cols]
  center <- vapply(factor_cols, function(v) weighted_mean_safe(base[[v]], base$BASE_WEIGHT), numeric(1))
  scale <- vapply(factor_cols, function(v) weighted_sd_safe(base[[v]], base$BASE_WEIGHT), numeric(1)); scale[!is.finite(scale) | scale <= 0] <- 1
  cal_cols <- paste0(factor_cols, "_CAL")
  for (j in seq_along(factor_cols)) {
    d[, (cal_cols[j]) := (get(factor_cols[j]) - center[j]) / scale[j]]
    base[, (cal_cols[j]) := (get(factor_cols[j]) - center[j]) / scale[j]]
  }
  bq <- base$BASE_WEIGHT / sum(base$BASE_WEIGHT); bx <- as.matrix(base[, ..cal_cols])
  targets <- list(
    full = as.vector(crossprod(bq, build_moment_matrix(bx, 2, TRUE))),
    means_squares = as.vector(crossprod(bq, build_moment_matrix(bx, 2, FALSE))),
    means = as.vector(crossprod(bq, build_moment_matrix(bx, 1, FALSE)))
  )
  dates <- sort(unique(d$DATE)); indices <- vector("list", length(dates)); diagnostics <- vector("list", length(dates))
  for (i in seq_along(dates)) {
    dd <- d[DATE == dates[i]]; q <- dd$EMPLOYMENT_WEIGHT / sum(dd$EMPLOYMENT_WEIGHT)
    fit <- calibrate_month(q, as.matrix(dd[, ..cal_cols]), targets)
    hard <- all(is.finite(fit$weights)) && is.finite(fit$max_error) && fit$max_error <= CALIBRATION_HARD_TOL &&
      is.finite(fit$effective_n) && fit$effective_n >= CALIBRATION_HARD_MIN_EFFECTIVE_OCCUPATIONS &&
      is.finite(fit$ratio) && fit$ratio <= CALIBRATION_HARD_MAX_WEIGHT_RATIO
    observed <- sum(q * dd$MEAN_WAGE)
    standardized <- if (hard) sum(fit$weights * dd$MEAN_WAGE) else NA_real_
    indices[[i]] <- data.table(
      DATE = dates[i], YEAR = year(dates[i]), MONTH = month(dates[i]), QUARTER_DATE = as.Date(as.yearqtr(dates[i])),
      OBSERVED_WAGE = observed, SKILL_STANDARDIZED_WAGE = standardized,
      COMPOSITION_GAP = observed - standardized, EMPLOYMENT_WEIGHT = sum(dd$EMPLOYMENT_WEIGHT),
      N_OCCUPATIONS = nrow(dd), SPECIFICATION = specification, OUTCOME = outcome, SAMPLE = sample,
      N_FACTORS = length(factor_cols), BASELINE_START = baseline[1], BASELINE_END = baseline[2]
    )
    diagnostics[[i]] <- data.table(
      DATE = dates[i], MOMENT_SET = fit$moment_set, ACCEPTED = fit$accepted, HARD_VALID = hard,
      CONVERGENCE = fit$convergence, MAX_MOMENT_ERROR = fit$max_error,
      EFFECTIVE_OCCUPATIONS = fit$effective_n, MAX_WEIGHT_RATIO = fit$ratio,
      N_OCCUPATIONS = nrow(dd), SPECIFICATION = specification, OUTCOME = outcome, SAMPLE = sample
    )
  }
  index <- rbindlist(indices); diag <- rbindlist(diagnostics)
  valid_share <- mean(diag$HARD_VALID)
  if (valid_share < MIN_VALID_CALIBRATION_SHARE) stop("Too few valid calibrations for ", specification)
  list(index = index, diagnostics = diag, occupation_panel = d)
}

continuous_decomposition <- function(index, time_var = "DATE", lag_periods = 12L) {
  value_cols <- c("OBSERVED_WAGE", "SKILL_STANDARDIZED_WAGE", "COMPOSITION_GAP")
  lag <- index[, c(time_var, value_cols), with = FALSE]
  setnames(lag, value_cols, paste0(value_cols, "_LAG"))
  if (time_var == "DATE") lag[, (time_var) := get(time_var) %m+% months(lag_periods)]
  else lag[, (time_var) := as.Date(as.yearqtr(get(time_var)) + lag_periods / 4)]
  out <- merge(index, lag, by = time_var, all.x = TRUE)
  out[, `:=`(
    TOTAL_PCT = 100 * (OBSERVED_WAGE / OBSERVED_WAGE_LAG - 1),
    WITHIN_PCT = 100 * (SKILL_STANDARDIZED_WAGE / SKILL_STANDARDIZED_WAGE_LAG - 1),
    COMPOSITION_PCT = 100 * ((COMPOSITION_GAP - COMPOSITION_GAP_LAG) / OBSERVED_WAGE_LAG)
  )]
  # Force exact additive accounting in percentage-point units.
  out[, COMPOSITION_PCT := TOTAL_PCT - WITHIN_PCT]
  out[, DECOMP_ERROR := TOTAL_PCT - WITHIN_PCT - COMPOSITION_PCT]
  out
}

quarterly_continuous_index <- function(index) {
  index[, .(
    OBSERVED_WAGE = weighted_mean_safe(OBSERVED_WAGE, EMPLOYMENT_WEIGHT),
    SKILL_STANDARDIZED_WAGE = weighted_mean_safe(SKILL_STANDARDIZED_WAGE, EMPLOYMENT_WEIGHT),
    COMPOSITION_GAP = weighted_mean_safe(COMPOSITION_GAP, EMPLOYMENT_WEIGHT),
    EMPLOYMENT_WEIGHT = sum(EMPLOYMENT_WEIGHT), N_OCCUPATIONS = mean(N_OCCUPATIONS),
    SPECIFICATION = SPECIFICATION[1], OUTCOME = OUTCOME[1], SAMPLE = SAMPLE[1], N_FACTORS = N_FACTORS[1]
  ), by = QUARTER_DATE]
}

# ------------------------ 9. MAIN ESTIMATION ------------------------------------

main_weekly_panel <- extract_outcome_panel(occ_panel_wide, "WEEKLY_EARNINGS_REAL", "MAIN")
main_weekly <- build_continuous_index_panel(
  main_weekly_panel, "WEEKLY_EARNINGS_REAL", "MAIN", PRIMARY_FACTOR_COLS, PRIMARY_BASELINE, "continuous_8factor_weekly"
)
main_weekly_monthly <- continuous_decomposition(main_weekly$index, "DATE", 12L)
main_weekly_quarterly <- continuous_decomposition(quarterly_continuous_index(main_weekly$index), "QUARTER_DATE", 4L)

secondary_actual_panel <- extract_outcome_panel(occ_panel_wide, "HOURLY_ACTUAL_REAL", "MAIN")
secondary_actual <- build_continuous_index_panel(
  secondary_actual_panel, "HOURLY_ACTUAL_REAL", "MAIN", PRIMARY_FACTOR_COLS, c(1989L, 1993L), "continuous_8factor_actual_hourly"
)
secondary_actual_monthly <- continuous_decomposition(secondary_actual$index, "DATE", 12L)
secondary_actual_quarterly <- continuous_decomposition(quarterly_continuous_index(secondary_actual$index), "QUARTER_DATE", 4L)

scalar_weekly_gp <- build_scalar_group_panel(
  main_weekly_panel, PRIMARY_SCALAR_SKILL, 10L, PRIMARY_BASELINE, "WEEKLY_EARNINGS_REAL", "MAIN", "scalar_deciles_weekly"
)
scalar_weekly_monthly <- exact_group_decomposition(scalar_weekly_gp, "DATE", 12L)
scalar_weekly_quarterly <- exact_group_decomposition(quarterly_group_panel(scalar_weekly_gp), "QUARTER_DATE", 4L)

fwrite(main_weekly$index, file.path(OUTPUT_DIR, "data", "main_continuous_weekly_monthly_index.csv"))
fwrite(main_weekly_monthly, file.path(OUTPUT_DIR, "data", "main_continuous_weekly_monthly_decomposition.csv"))
fwrite(main_weekly_quarterly, file.path(OUTPUT_DIR, "data", "main_continuous_weekly_quarterly_decomposition.csv"))
fwrite(main_weekly$diagnostics, file.path(OUTPUT_DIR, "tables", "main_continuous_calibration_diagnostics.csv"))
fwrite(secondary_actual_monthly, file.path(OUTPUT_DIR, "data", "secondary_actual_hourly_monthly_decomposition.csv"))
fwrite(scalar_weekly_monthly, file.path(OUTPUT_DIR, "data", "benchmark_scalar_weekly_monthly_decomposition.csv"))

# --------------------- 10. CONVENTIONAL CELL ADJUSTMENT -------------------------

cell_standardized_index <- function(panel, baseline, specification) {
  if (!nrow(panel)) return(NULL)
  id_cols <- setdiff(names(panel), c("DATE", "YEAR", "MONTH", "QUARTER_DATE", "EMPLOYMENT_WEIGHT", "WAGE_SUM", "N", "MEAN_WAGE", "OUTCOME", "CELL_TYPE"))
  panel[, CELL_ID := do.call(paste, c(.SD, sep = "|")), .SDcols = id_cols]
  baseline <- choose_baseline(panel, baseline)
  target <- panel[YEAR >= baseline[1] & YEAR <= baseline[2], .(TARGET_WEIGHT = sum(EMPLOYMENT_WEIGHT)), by = CELL_ID]
  target[, TARGET_SHARE := TARGET_WEIGHT / sum(TARGET_WEIGHT)]
  dates <- sort(unique(panel$DATE)); out <- vector("list", length(dates))
  for (i in seq_along(dates)) {
    current_cells <- panel[DATE == dates[i]]
    current <- sum(current_cells$EMPLOYMENT_WEIGHT * current_cells$MEAN_WAGE) / sum(current_cells$EMPLOYMENT_WEIGHT)
    d <- merge(current_cells, target[, .(CELL_ID, TARGET_SHARE)], by = "CELL_ID")
    coverage <- sum(d$TARGET_SHARE)
    adjusted <- if (coverage >= 0.95) sum((d$TARGET_SHARE / coverage) * d$MEAN_WAGE) else NA_real_
    out[[i]] <- data.table(
      DATE = dates[i], YEAR = year(dates[i]), MONTH = month(dates[i]), QUARTER_DATE = as.Date(as.yearqtr(dates[i])),
      OBSERVED_WAGE = current, SKILL_STANDARDIZED_WAGE = adjusted, COMPOSITION_GAP = current - adjusted,
      EMPLOYMENT_WEIGHT = sum(d$EMPLOYMENT_WEIGHT), N_OCCUPATIONS = uniqueN(d$CELL_ID),
      SPECIFICATION = specification, OUTCOME = "WEEKLY_EARNINGS_REAL", SAMPLE = "MAIN", N_FACTORS = 0L,
      TARGET_CELL_COVERAGE = coverage
    )
  }
  idx <- rbindlist(out)
  list(index = idx, decomposition = continuous_decomposition(idx, "DATE", 12L))
}

demo_core_result <- if (RUN_DEMOGRAPHIC_CELL_ADJUSTMENT) cell_standardized_index(demo_core_panel, PRIMARY_BASELINE, "conventional_demographic_cells") else NULL
demo_ind_result <- if (RUN_INDUSTRY_CELL_ROBUSTNESS) cell_standardized_index(demo_ind_panel, PRIMARY_BASELINE, "conventional_demographic_industry_cells") else NULL
if (!is.null(demo_core_result)) fwrite(demo_core_result$decomposition, file.path(OUTPUT_DIR, "data", "conventional_demographic_cell_decomposition.csv"))
if (!is.null(demo_ind_result)) fwrite(demo_ind_result$decomposition, file.path(OUTPUT_DIR, "data", "conventional_demographic_industry_cell_decomposition.csv"))

# ----------------------- 11. REGRESSION FUNCTIONS -------------------------------

prepare_monthly_regression <- function(decomp) {
  d <- merge(decomp, unemployment, by = "DATE", all.x = TRUE)
  d <- merge(d, recessions, by = "DATE", all.x = TRUE)
  ul <- unemployment[, .(DATE = DATE %m+% months(12), UNEMP_L12 = UNEMP_RATE)]
  d <- merge(d, ul, by = "DATE", all.x = TRUE)
  d[, `:=`(
    D12_UNEMP = UNEMP_RATE - UNEMP_L12,
    D12_REFERENCE_DATE = DATE %m-% months(12L),
    TIME_ID = 12L * (year(DATE) - START_YEAR) + month(DATE),
    MONTH_FACTOR = factor(month(DATE))
  )]
  # The current-date episode is useful descriptively. STRICT_POST_TRANSITION
  # separately identifies observations for which both t and t-12 are after the
  # pandemic/recovery transition.
  d[, EPISODE := fcase(
    DATE < PANDEMIC_TRANSITION_START, "Clean pre-COVID",
    DATE <= PANDEMIC_TRANSITION_END, "Pandemic/recovery transition",
    DATE >= POST_TRANSITION_START, "Post-transition current date",
    default = NA_character_
  )]
  d[, `:=`(
    STRICT_POST_TRANSITION = (
      DATE >= STRICT_POST_TRANSITION_START &
        D12_REFERENCE_DATE > PANDEMIC_TRANSITION_END
    ),
    STRICT_NONTRANSITION_SAMPLE = (
      DATE < PANDEMIC_TRANSITION_START |
        D12_REFERENCE_DATE > PANDEMIC_TRANSITION_END
    )
  )]
  d
}

prepare_quarterly_regression <- function(decomp) {
  d <- merge(decomp, output_gap, by = "QUARTER_DATE", all.x = TRUE)
  ql <- output_gap[, .(QUARTER_DATE = as.Date(as.yearqtr(QUARTER_DATE) + 1), OUTPUT_GAP_L4 = OUTPUT_GAP)]
  d <- merge(d, ql, by = "QUARTER_DATE", all.x = TRUE)
  uq <- unemployment[, .(UNEMP_RATE = mean(UNEMP_RATE, na.rm = TRUE)), by = .(QUARTER_DATE = as.Date(as.yearqtr(DATE)))]
  uql <- uq[, .(QUARTER_DATE = as.Date(as.yearqtr(QUARTER_DATE) + 1), UNEMP_L4 = UNEMP_RATE)]
  d <- merge(d, uq, by = "QUARTER_DATE", all.x = TRUE); d <- merge(d, uql, by = "QUARTER_DATE", all.x = TRUE)
  d[, `:=`(
    D4_OUTPUT_GAP = OUTPUT_GAP - OUTPUT_GAP_L4,
    D4_UNEMP = UNEMP_RATE - UNEMP_L4,
    QUARTER_FACTOR = factor(quarter(QUARTER_DATE)),
    QTIME_ID = 4L * (year(QUARTER_DATE) - START_YEAR) + quarter(QUARTER_DATE)
  )]
  d
}

fit_monthly_components <- function(d, prefix, type = "continuous") {
  outcomes <- c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT")
  names_out <- c("aggregate", "standardized_or_within", "composition")
  rbindlist(lapply(seq_along(outcomes), function(i) {
    m <- feols(as.formula(paste0(outcomes[i], " ~ D12_UNEMP | MONTH_FACTOR")), data = d, vcov = NW(12) ~ TIME_ID)
    coef_row(m, "D12_UNEMP", paste0(prefix, "_", names_out[i]), outcomes[i], "12-month unemployment change", type = type)
  }))
}

main_reg <- prepare_monthly_regression(main_weekly_monthly)
scalar_reg <- prepare_monthly_regression(scalar_weekly_monthly)
actual_reg <- prepare_monthly_regression(secondary_actual_monthly)
main_coefficients <- rbindlist(list(
  fit_monthly_components(main_reg, "main_continuous_weekly", "continuous 8-factor"),
  fit_monthly_components(scalar_reg, "benchmark_scalar_weekly", "scalar deciles"),
  fit_monthly_components(actual_reg, "secondary_actual_hourly", "continuous 8-factor")
), fill = TRUE)
if (!is.null(demo_core_result)) main_coefficients <- rbind(main_coefficients, fit_monthly_components(prepare_monthly_regression(demo_core_result$decomposition), "conventional_demographic_cells", "conventional cells"), fill = TRUE)
fwrite(main_coefficients, file.path(OUTPUT_DIR, "tables", "main_cyclicality_coefficients.csv"))

# Attribution ratios. Build one row per estimation method rather than grouping
# by the component outcome. The Version 3.1 grouping created three partial rows
# per method and therefore left COMPOSITION_SHARE blank.
attribution_long <- copy(main_coefficients)
attribution_long[, COMPONENT := fcase(
  grepl("_aggregate$", model), "aggregate",
  grepl("_standardized_or_within$", model), "standardized",
  grepl("_composition$", model), "composition",
  default = NA_character_
)]
attribution_long[, METHOD_ID := sub(
  "_(aggregate|standardized_or_within|composition)$", "", model
)]
attribution_long[, WAGE_OUTCOME := fcase(
  grepl("actual_hourly", METHOD_ID), "HOURLY_ACTUAL_REAL",
  default = "WEEKLY_EARNINGS_REAL"
)]
attribution_long[, METHOD_LABEL := fcase(
  METHOD_ID == "main_continuous_weekly", "Eight-factor continuous: weekly earnings",
  METHOD_ID == "benchmark_scalar_weekly", "Scalar deciles: weekly earnings",
  METHOD_ID == "secondary_actual_hourly", "Eight-factor continuous: actual-hour earnings",
  METHOD_ID == "conventional_demographic_cells", "Demographic cells: weekly earnings",
  default = METHOD_ID
)]

attribution_estimates <- dcast(
  attribution_long[!is.na(COMPONENT)],
  METHOD_ID + METHOD_LABEL + specification_type + WAGE_OUTCOME + cycle_variable + sample ~ COMPONENT,
  value.var = "estimate"
)
attribution_ses <- dcast(
  attribution_long[!is.na(COMPONENT)],
  METHOD_ID + METHOD_LABEL + specification_type + WAGE_OUTCOME + cycle_variable + sample ~ COMPONENT,
  value.var = "std_error"
)
se_component_cols <- intersect(
  c("aggregate", "standardized", "composition"),
  names(attribution_ses)
)
setnames(
  attribution_ses,
  se_component_cols,
  paste0(se_component_cols, "_SE")
)
attribution_pvalues <- dcast(
  attribution_long[!is.na(COMPONENT)],
  METHOD_ID + METHOD_LABEL + specification_type + WAGE_OUTCOME + cycle_variable + sample ~ COMPONENT,
  value.var = "p_value"
)
pvalue_component_cols <- intersect(
  c("aggregate", "standardized", "composition"),
  names(attribution_pvalues)
)
setnames(
  attribution_pvalues,
  pvalue_component_cols,
  paste0(pvalue_component_cols, "_P_VALUE")
)

attribution_statistics <- merge(
  attribution_estimates,
  attribution_ses,
  by = c("METHOD_ID", "METHOD_LABEL", "specification_type", "WAGE_OUTCOME", "cycle_variable", "sample"),
  all = TRUE
)
attribution_statistics <- merge(
  attribution_statistics,
  attribution_pvalues,
  by = c("METHOD_ID", "METHOD_LABEL", "specification_type", "WAGE_OUTCOME", "cycle_variable", "sample"),
  all = TRUE
)
attribution_statistics[, `:=`(
  STANDARDIZED_SHARE = fifelse(is.finite(aggregate) & aggregate != 0, standardized / aggregate, NA_real_),
  COMPOSITION_SHARE = fifelse(is.finite(aggregate) & aggregate != 0, composition / aggregate, NA_real_),
  COMPOSITION_SHARE_PERCENT = fifelse(is.finite(aggregate) & aggregate != 0, 100 * composition / aggregate, NA_real_),
  ACCOUNTING_ERROR = aggregate - standardized - composition
)]
setorder(attribution_statistics, WAGE_OUTCOME, METHOD_ID)

if (any(!is.finite(attribution_statistics$COMPOSITION_SHARE))) {
  stop("Attribution shares could not be computed for one or more main methods.", call. = FALSE)
}
if (max(abs(attribution_statistics$ACCOUNTING_ERROR), na.rm = TRUE) > 1e-8) {
  warning("Attribution coefficient accounting error exceeds tolerance.", call. = FALSE)
}
fwrite(attribution_statistics, file.path(OUTPUT_DIR, "tables", "main_attribution_statistics.csv"))
fwrite(attribution_long, file.path(OUTPUT_DIR, "tables", "main_attribution_components_long.csv"))

# Quarterly corroboration.
qmain <- prepare_quarterly_regression(main_weekly_quarterly)
quarterly_coefficients <- rbindlist(lapply(c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT"), function(y) {
  m1 <- feols(as.formula(paste0(y, " ~ D4_OUTPUT_GAP | QUARTER_FACTOR")), qmain, vcov = NW(4) ~ QTIME_ID)
  m2 <- feols(as.formula(paste0(y, " ~ D4_UNEMP | QUARTER_FACTOR")), qmain, vcov = NW(4) ~ QTIME_ID)
  rbind(
    coef_row(m1, "D4_OUTPUT_GAP", paste0("quarterly_output_gap_", y), y, "four-quarter output-gap change"),
    coef_row(m2, "D4_UNEMP", paste0("quarterly_unemployment_", y), y, "four-quarter unemployment change")
  )
}))
fwrite(quarterly_coefficients, file.path(OUTPUT_DIR, "tables", "quarterly_cyclicality_coefficients.csv"))

# ----------------------- 12. EPISODE ANALYSIS -----------------------------------

fit_episode_components <- function(d, label, model_slug) {
  valid <- d[complete.cases(TOTAL_PCT, WITHIN_PCT, COMPOSITION_PCT, D12_UNEMP)]
  n_valid <- nrow(valid)
  if (n_valid < MIN_EPISODE_OBSERVATIONS || uniqueN(valid$D12_UNEMP) < 3L) {
    return(data.table())
  }
  out <- fit_monthly_components(
    valid,
    paste0("episode_", model_slug),
    "episode"
  )
  out[, `:=`(
    EPISODE = label,
    PERIOD_START = min(valid$DATE),
    PERIOD_END = max(valid$DATE),
    N_COMPLETE_MONTHS = n_valid
  )]
  out
}

episode_definitions <- data.table(
  EPISODE = c(
    "Clean pre-COVID",
    "Pandemic/recovery transition",
    "Post-transition current date",
    "Strict clean post-transition",
    "Strict non-transition sample"
  ),
  DEFINITION = c(
    "Current month is before March 2020.",
    "Current month is March 2020 through December 2022.",
    "Current month is January 2023 or later; some 12-month references may still fall in 2022.",
    "Both the current month and its 12-month reference are after December 2022; begins January 2024.",
    "Clean pre-COVID observations combined with strict clean post-transition observations."
  )
)
fwrite(episode_definitions, file.path(OUTPUT_DIR, "tables", "episode_definitions.csv"))

episode_coefficients <- data.table()
if (RUN_EPISODE_ANALYSIS) {
  episode_parts <- list(
    fit_episode_components(
      main_reg[EPISODE == "Clean pre-COVID"],
      "Clean pre-COVID", "clean_pre_covid"
    ),
    fit_episode_components(
      main_reg[EPISODE == "Pandemic/recovery transition"],
      "Pandemic/recovery transition", "pandemic_recovery_transition"
    ),
    fit_episode_components(
      main_reg[EPISODE == "Post-transition current date"],
      "Post-transition current date", "post_transition_current_date"
    ),
    fit_episode_components(
      main_reg[STRICT_POST_TRANSITION == TRUE],
      "Strict clean post-transition", "strict_clean_post_transition"
    ),
    fit_episode_components(
      main_reg[STRICT_NONTRANSITION_SAMPLE == TRUE],
      "Strict non-transition sample", "strict_nontransition_sample"
    )
  )
  episode_coefficients <- rbindlist(episode_parts, fill = TRUE)
  fwrite(
    episode_coefficients,
    file.path(OUTPUT_DIR, "tables", "episode_specific_cyclicality_coefficients.csv")
  )
}

# Leave-one-recession-influence-window-out. For 12-month changes, omitting only
# official recession months leaves comparisons whose current or reference month is
# still affected by that recession. Each window therefore begins 12 months before
# the NBER start and ends 12 months after the NBER end. The COVID window is
# extended through December 2022 to cover the pandemic/reopening transition.
leave_one_recession <- data.table()
recession_influence_windows <- data.table()
if (RUN_RECESSION_ANALYSIS) {
  rec <- recessions[RECESSION == 1][order(DATE)]
  rec[, SPELL := cumsum(c(1L, as.integer(diff(DATE) > 35)))]
  spells <- rec[, .(RECESSION_START = min(DATE), RECESSION_END = max(DATE)), by = SPELL]
  spells[, `:=`(
    INFLUENCE_START = RECESSION_START %m-% months(RECESSION_INFLUENCE_PRE_MONTHS),
    INFLUENCE_END = RECESSION_END %m+% months(RECESSION_INFLUENCE_POST_MONTHS)
  )]
  spells[
    RECESSION_START == COVID_RECESSION_START & INFLUENCE_END < COVID_RECESSION_INFLUENCE_END,
    INFLUENCE_END := COVID_RECESSION_INFLUENCE_END
  ]
  spells[, `:=`(
    REGRESSION_SAMPLE_START = min(main_reg$DATE, na.rm = TRUE),
    REGRESSION_SAMPLE_END = max(main_reg$DATE, na.rm = TRUE)
  )]
  spells[, N_OMITTED_COMPLETE_MONTHS := vapply(seq_len(nrow(spells)), function(i) {
    main_reg[
      DATE >= spells$INFLUENCE_START[i] & DATE <= spells$INFLUENCE_END[i] &
        complete.cases(TOTAL_PCT, WITHIN_PCT, COMPOSITION_PCT, D12_UNEMP),
      .N
    ]
  }, integer(1))]
  recession_influence_windows <- copy(spells)
  fwrite(
    recession_influence_windows,
    file.path(OUTPUT_DIR, "tables", "recession_influence_windows.csv")
  )

  leave_one_recession <- rbindlist(lapply(seq_len(nrow(spells)), function(i) {
    if (spells$N_OMITTED_COMPLETE_MONTHS[i] == 0L) return(data.table())
    d <- main_reg[
      DATE < spells$INFLUENCE_START[i] |
        DATE > spells$INFLUENCE_END[i]
    ]
    valid_remaining <- d[complete.cases(TOTAL_PCT, WITHIN_PCT, COMPOSITION_PCT, D12_UNEMP), .N]
    if (valid_remaining < 60L) return(data.table())
    out <- fit_monthly_components(
      d,
      paste0("omit_recession_influence_", format(spells$RECESSION_START[i], "%Y%m")),
      "leave-one-recession-influence-window-out"
    )
    out[, `:=`(
      OMITTED_RECESSION_START = spells$RECESSION_START[i],
      OMITTED_RECESSION_END = spells$RECESSION_END[i],
      OMITTED_INFLUENCE_START = spells$INFLUENCE_START[i],
      OMITTED_INFLUENCE_END = spells$INFLUENCE_END[i],
      N_OMITTED_COMPLETE_MONTHS = spells$N_OMITTED_COMPLETE_MONTHS[i]
    )]
    out
  }), fill = TRUE)
  fwrite(
    leave_one_recession,
    file.path(OUTPUT_DIR, "tables", "leave_one_recession_influence_window_out_coefficients.csv")
  )
  # Retain the prior filename for downstream workflows, but its contents now use
  # the corrected influence-window definition.
  fwrite(
    leave_one_recession,
    file.path(OUTPUT_DIR, "tables", "leave_one_recession_out_coefficients.csv")
  )
}

# ------------------------ 13. FACTOR EMPLOYMENT ---------------------------------

factor_employment <- main_weekly_panel[, c(
  list(EMPLOYMENT_WEIGHT = sum(EMPLOYMENT_WEIGHT)),
  lapply(.SD, function(z) weighted_mean_safe(z, EMPLOYMENT_WEIGHT))
), by = DATE, .SDcols = PRIMARY_FACTOR_COLS]
factor_employment <- merge(factor_employment, unemployment, by = "DATE", all.x = TRUE)
fl <- factor_employment[, c("DATE", PRIMARY_FACTOR_COLS), with = FALSE]
fl[, DATE := DATE %m+% months(12)]; setnames(fl, PRIMARY_FACTOR_COLS, paste0(PRIMARY_FACTOR_COLS, "_L12"))
factor_employment <- merge(factor_employment, fl, by = "DATE", all.x = TRUE)
ul <- unemployment[, .(DATE = DATE %m+% months(12), UNEMP_L12 = UNEMP_RATE)]
factor_employment <- merge(factor_employment, ul, by = "DATE", all.x = TRUE)
factor_employment[, `:=`(D12_UNEMP = UNEMP_RATE - UNEMP_L12, TIME_ID = 12L * (year(DATE) - START_YEAR) + month(DATE), MONTH_FACTOR = factor(month(DATE)))]
factor_employment_coefficients <- rbindlist(lapply(seq_along(PRIMARY_FACTOR_COLS), function(j) {
  v <- PRIMARY_FACTOR_COLS[j]; dv <- paste0("D12_", v)
  factor_employment[, (dv) := get(v) - get(paste0(v, "_L12"))]
  m <- feols(as.formula(paste0(dv, " ~ D12_UNEMP | MONTH_FACTOR")), factor_employment, vcov = NW(12) ~ TIME_ID)
  out <- coef_row(m, "D12_UNEMP", paste0("factor_employment_", v), dv, "12-month unemployment change", type = "factor employment")
  out[, FACTOR_LABEL := PRIMARY_FACTOR_LABELS[j]]; out
}))
fwrite(factor_employment_coefficients, file.path(OUTPUT_DIR, "tables", "factor_employment_cyclicality.csv"))

# -------------------------- 14. ROBUSTNESS --------------------------------------

robustness_coefficients <- data.table()

# Alternative samples and outcomes.
robust_specs <- list(
  prime_age_weekly = list(outcome = "WEEKLY_EARNINGS_REAL", sample = "PRIME_AGE", factors = PRIMARY_FACTOR_COLS, baseline = PRIMARY_BASELINE),
  full_time_weekly = list(outcome = "WEEKLY_EARNINGS_REAL", sample = "FULL_TIME", factors = PRIMARY_FACTOR_COLS, baseline = c(1989L, 1993L)),
  reported_hourly = list(outcome = "HOURLY_REPORTED_REAL", sample = "MAIN", factors = PRIMARY_FACTOR_COLS, baseline = PRIMARY_BASELINE),
  usual_hourly = list(outcome = "HOURLY_USUAL_REAL", sample = "MAIN", factors = PRIMARY_FACTOR_COLS, baseline = c(1994L, 1998L)),
  weekly_common_cap = list(outcome = "WEEKLY_EARNINGS_REAL_COMMON_CAP", sample = "MAIN", factors = PRIMARY_FACTOR_COLS, baseline = PRIMARY_BASELINE),
  actual_common_cap = list(outcome = "HOURLY_ACTUAL_REAL_COMMON_CAP", sample = "MAIN", factors = PRIMARY_FACTOR_COLS, baseline = c(1989L, 1993L))
)
for (nm in names(robust_specs)) {
  s <- robust_specs[[nm]]
  rp <- extract_outcome_panel(occ_panel_wide, s$outcome, s$sample)
  if (nrow(rp)) {
    r <- build_continuous_index_panel(rp, s$outcome, s$sample, s$factors, s$baseline, nm)
    rd <- prepare_monthly_regression(continuous_decomposition(r$index, "DATE", 12L))
    robustness_coefficients <- rbind(robustness_coefficients, fit_monthly_components(rd, nm, "outcome/sample robustness"), fill = TRUE)
  }
}

# Alternative factor counts.
if (RUN_ALT_FACTOR_COUNTS) {
  for (k in c(5L, 10L, 13L)) {
    cols <- sprintf("F%02d_%02d", k, seq_len(k)); cols <- intersect(cols, names(main_weekly_panel))
    if (length(cols) == k) {
      r <- build_continuous_index_panel(main_weekly_panel, "WEEKLY_EARNINGS_REAL", "MAIN", cols, PRIMARY_BASELINE, paste0("continuous_", k, "factor"))
      rd <- prepare_monthly_regression(continuous_decomposition(r$index, "DATE", 12L))
      robustness_coefficients <- rbind(robustness_coefficients, fit_monthly_components(rd, paste0("continuous_", k, "factor"), "factor-count robustness"), fill = TRUE)
    }
  }
}

# Alternative baseline.
if (RUN_BASELINE_SENSITIVITY) {
  r <- build_continuous_index_panel(main_weekly_panel, "WEEKLY_EARNINGS_REAL", "MAIN", PRIMARY_FACTOR_COLS, ALT_BASELINE, "baseline_1991_1995")
  rd <- prepare_monthly_regression(continuous_decomposition(r$index, "DATE", 12L))
  robustness_coefficients <- rbind(robustness_coefficients, fit_monthly_components(rd, "baseline_1991_1995", "baseline robustness"), fill = TRUE)
}

# Scalar group grid.
if (RUN_SCALAR_GRID) {
  for (g in c(5L, 10L, 20L)) {
    gp <- build_scalar_group_panel(main_weekly_panel, PRIMARY_SCALAR_SKILL, g, PRIMARY_BASELINE, "WEEKLY_EARNINGS_REAL", "MAIN", paste0("scalar_", g, "groups"))
    d <- prepare_monthly_regression(exact_group_decomposition(gp, "DATE", 12L))
    robustness_coefficients <- rbind(robustness_coefficients, fit_monthly_components(d, paste0("scalar_", g, "groups"), "scalar group-count robustness"), fill = TRUE)
  }
}
fwrite(robustness_coefficients, file.path(OUTPUT_DIR, "tables", "robustness_cyclicality_coefficients.csv"))

# ---------------------- 15. MOVING-BLOCK BOOTSTRAP ------------------------------

block_bootstrap <- data.table()
if (RUN_BLOCK_BOOTSTRAP) {
  boot_data <- main_reg[complete.cases(TOTAL_PCT, WITHIN_PCT, COMPOSITION_PCT, D12_UNEMP)]
  n <- nrow(boot_data); set.seed(BOOTSTRAP_SEED)
  draw_blocks <- function() {
    starts <- sample(seq_len(max(1L, n - BOOTSTRAP_BLOCK_MONTHS + 1L)), ceiling(n / BOOTSTRAP_BLOCK_MONTHS), replace = TRUE)
    unlist(lapply(starts, function(s) s:min(n, s + BOOTSTRAP_BLOCK_MONTHS - 1L)))[seq_len(n)]
  }
  block_bootstrap <- rbindlist(lapply(seq_len(BOOTSTRAP_REPS), function(b) {
    idx <- draw_blocks(); d <- copy(boot_data[idx]); d[, `:=`(TIME_ID_BOOT = .I, MONTH_FACTOR_BOOT = factor(month(DATE)))]
    coefs <- vapply(c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT"), function(y) {
      m <- tryCatch(feols(as.formula(paste0(y, " ~ D12_UNEMP | MONTH_FACTOR_BOOT")), d, vcov = NW(12) ~ TIME_ID_BOOT), error = function(e) NULL)
      if (is.null(m)) NA_real_ else unname(coef(m)["D12_UNEMP"])
    }, numeric(1))
    data.table(REPLICATION = b, AGGREGATE = coefs[1], STANDARDIZED = coefs[2], COMPOSITION = coefs[3])
  }))
  bootstrap_summary <- melt(block_bootstrap, id.vars = "REPLICATION", variable.name = "COMPONENT", value.name = "ESTIMATE")[, .(
    MEAN = mean(ESTIMATE, na.rm = TRUE), SD = sd(ESTIMATE, na.rm = TRUE),
    P025 = quantile(ESTIMATE, 0.025, na.rm = TRUE), MEDIAN = median(ESTIMATE, na.rm = TRUE), P975 = quantile(ESTIMATE, 0.975, na.rm = TRUE),
    VALID_REPS = sum(is.finite(ESTIMATE))
  ), by = COMPONENT]
  fwrite(block_bootstrap, file.path(OUTPUT_DIR, "tables", "main_block_bootstrap_draws.csv"))
  fwrite(bootstrap_summary, file.path(OUTPUT_DIR, "tables", "main_block_bootstrap_summary.csv"))
}

# ------------------------------ 16. FIGURES -------------------------------------

plot_main <- melt(
  main_reg, id.vars = "DATE", measure.vars = c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT"),
  variable.name = "COMPONENT", value.name = "PERCENT_CHANGE"
)
plot_main[, COMPONENT := factor(COMPONENT, levels = c("TOTAL_PCT", "WITHIN_PCT", "COMPOSITION_PCT"), labels = c("Observed", "Skill-standardized", "Skill composition"))]
p1 <- ggplot(plot_main, aes(DATE, PERCENT_CHANGE, linetype = COMPONENT)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_line(linewidth = 0.65, na.rm = TRUE) +
  labs(title = "Real weekly earnings growth and multidimensional skill composition", x = NULL, y = "12-month change (percent)", linetype = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
save_plot(p1, "figure_1_weekly_earnings_decomposition.png", 11, 6)

coef_plot_data <- copy(attribution_statistics)
p2 <- ggplot(coef_plot_data, aes(composition, reorder(METHOD_LABEL, composition))) +
  geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(size = 2.5) +
  geom_errorbarh(
    aes(
      xmin = composition - 1.96 * composition_SE,
      xmax = composition + 1.96 * composition_SE
    ),
    height = 0.15
  ) +
  labs(
    title = "Composition contribution across adjustment methods",
    x = "Coefficient on 12-month unemployment change",
    y = NULL
  ) + theme_minimal(base_size = 12)
save_plot(p2, "figure_2_method_comparison.png", 9, 5)

if (nrow(episode_coefficients) > 0L) {
  episode_plot_data <- episode_coefficients[grepl("_composition$", model)]
  episode_plot_data[, EPISODE := factor(
    EPISODE,
    levels = c(
      "Clean pre-COVID",
      "Pandemic/recovery transition",
      "Post-transition current date",
      "Strict clean post-transition",
      "Strict non-transition sample"
    )
  )]
  p3 <- ggplot(episode_plot_data, aes(estimate, EPISODE)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(size = 2.5) +
    geom_errorbarh(
      aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error),
      height = 0.15
    ) +
    labs(
      title = "Occupational skill-composition cyclicality by historical episode",
      x = "Composition coefficient on 12-month unemployment change",
      y = NULL
    ) + theme_minimal(base_size = 12)
  save_plot(p3, "figure_3_episode_specific_composition.png", 10, 6)
}

# -------------------------- 17. FINAL SUMMARY -----------------------------------

analysis_summary <- data.table(
  ITEM = c(
    "Start year", "End year", "CPS extracts", "Compact wide occupation-month rows",
    "Overall weighted skill match", "Primary baseline", "Primary calibration valid share",
    "Maximum continuous decomposition error", "Maximum scalar decomposition error",
    "Pandemic transition window", "Strict post-transition start",
    "Attribution rows with finite composition shares"
  ),
  VALUE = c(
    START_YEAR, END_YEAR, length(cps_files), nrow(occ_panel_wide),
    overall_match, paste(PRIMARY_BASELINE, collapse = "-"), mean(main_weekly$diagnostics$HARD_VALID),
    max(abs(main_weekly_monthly$DECOMP_ERROR), na.rm = TRUE),
    max(abs(scalar_weekly_monthly$DECOMP_ERROR), na.rm = TRUE),
    paste(PANDEMIC_TRANSITION_START, PANDEMIC_TRANSITION_END, sep = " to "),
    as.character(STRICT_POST_TRANSITION_START),
    sum(is.finite(attribution_statistics$COMPOSITION_SHARE))
  )
)
fwrite(analysis_summary, file.path(OUTPUT_DIR, "tables", "analysis_summary.csv"))

# Critical completion checks.
critical_files <- c(
  file.path(OUTPUT_DIR, "tables", "main_cyclicality_coefficients.csv"),
  file.path(OUTPUT_DIR, "tables", "main_attribution_statistics.csv"),
  file.path(OUTPUT_DIR, "tables", "robustness_cyclicality_coefficients.csv"),
  file.path(OUTPUT_DIR, "data", "main_continuous_weekly_monthly_decomposition.csv")
)
if (RUN_EPISODE_ANALYSIS) {
  critical_files <- c(
    critical_files,
    file.path(OUTPUT_DIR, "tables", "episode_specific_cyclicality_coefficients.csv")
  )
}
if (RUN_RECESSION_ANALYSIS) {
  critical_files <- c(
    critical_files,
    file.path(OUTPUT_DIR, "tables", "recession_influence_windows.csv"),
    file.path(OUTPUT_DIR, "tables", "leave_one_recession_influence_window_out_coefficients.csv")
  )
}
if (!all(file.exists(critical_files))) {
  missing_critical <- critical_files[!file.exists(critical_files)]
  stop(
    "Critical output files were not created:\n",
    paste(missing_critical, collapse = "\n"),
    call. = FALSE
  )
}
if (!exists("main_coefficients") || nrow(main_coefficients) == 0L) {
  stop("main_coefficients was not created or is empty.", call. = FALSE)
}
if (!all(is.finite(main_coefficients$estimate))) {
  warning("Some main coefficient estimates are nonfinite; inspect tables.")
}
if (!exists("attribution_statistics") || nrow(attribution_statistics) == 0L ||
    any(!is.finite(attribution_statistics$COMPOSITION_SHARE))) {
  stop("Attribution statistics are missing or contain nonfinite composition shares.", call. = FALSE)
}
if (RUN_EPISODE_ANALYSIS && (!exists("episode_coefficients") || nrow(episode_coefficients) == 0L)) {
  stop("Episode analysis was requested but produced no coefficient rows.", call. = FALSE)
}
if (RUN_RECESSION_ANALYSIS && (!exists("leave_one_recession") || nrow(leave_one_recession) == 0L)) {
  stop("Recession influence-window analysis was requested but produced no coefficient rows.", call. = FALSE)
}

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "logs", "session_info.txt"))
writeLines("SUCCESS", status_file)
log_message("Analysis completed successfully. Results: ", OUTPUT_DIR)

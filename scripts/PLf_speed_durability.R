#####################################################################################################
#                                                                                                   #
#                    Fatigue-integrated Power Law (PLf) and Speed Durability                        #
#                    -------------------------------------------------------                        #
#                                                                                                   #
#   0. Load the libraries and datasets                                                              # 
#   1. Prepare the 800-1500m datasets                                                               #
#   2. Estimation of the PL(f) parameters: goodness-of-fit and uncertainty analysis                 #
#   3. Data description                                                                             #
#   4. Fatigue-intergrated Power Law (PLf) model computation                                        #
#   5. S durability / S decline analysis                                                            #
#   6. Risk ratio analysis                                                                          #
#   7. Performance relation to Fresh S, S durability and position                                   #
#   8. ECSS YIA talk - slide figures                                                                #
#                                                                                                   #
#                                                                        Maxime Walt, 2025-2026     #
#                                                                                                   #
#####################################################################################################

# -----------------------------------------------------------------------------
# READ ME:
#  In this script, we are applying the Fatigue-integrated PL (PLf) model 
#  (based on Drake eq.) DL and MC 100-m splits.
#  We assume downshifting S (based on Drake eq.), downshifting also depends 
#  on the positions changes (tau as a tactical multiplier).
#  Purpose: Does the model align with pacing literature? 
#           Does it link with performance?
#  -----------------------------------------------------------------------------

# =============================================================================
# 0. Load the libraries and datasets
# =============================================================================

# 0.1. Load the libraries

library(readxl); library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(patchwork)
library(ggpubr); library(purrr); library(minpack.lm); library(broom); library(scales); library(tibble)
library(pracma); library(splines); library(car); library(corrplot); library(effects); library(ggeffects)
library(sjPlot); library(extrafont); library(ggiraphExtra); library(randomForest); library(lme4)
library(rstatix); library(ggtext); library(stringr); library(stringi); library(ggh4x); library(logistf)
library(geomtextpath); library(influence.ME); library(performance); library(emmeans); library(broom.mixed)
library(pROC); library(DHARMa)  

# 0.2. Load the datasets

# SBs datasets

load("C:/Users/maxim/Documents/UNIL - Doctorat/1_Research_data/3_Pacing_profiles/Pacing_profiles/WA_SB_best_PL")
load("C:/Users/maxim/Documents/UNIL - Doctorat/1_Research_data/3_Pacing_profiles/Pacing_profiles/results_PL_WA_SB")
load("C:/Users/maxim/Documents/UNIL - Doctorat/1_Research_data/3_Pacing_profiles/Pacing_profiles/results_PL_WA_SB_unnested")

# Splits datasets

Pacing_splits_PHD_800 <- read_excel("~/UNIL - Doctorat/1_Research_data/3_Pacing_profiles/Splits data/Pacing_splits_PHD_800m.xlsx")

Pacing_splits_PHD_800$Competition_date <- as.Date(Pacing_splits_PHD_800$Competition_date)

Pacing_splits_PHD_1500 <- read_excel(
  "~/UNIL - Doctorat/1_Research_data/3_Pacing_profiles/Splits data/Pacing_splits_PHD_1500m.xlsx",
  col_types = "text"
)

# Convert split columns to numeric (adapt column names to your sheet)
num_cols <- c("0_100m", "100_200m", "200_300m", "300_400m", "400_500m", "500_600m", "600_700m", "700_800m", "800_900m", "900_1000m",
              "1000_1100m", "1100_1200m", "1200_1300m", "1300_1400m", "1400_1500m", "Final_time", "Sum_time", "Diff", "Competition_year",
              "Temperature_C", "Humidity")
Pacing_splits_PHD_1500[num_cols] <- lapply(Pacing_splits_PHD_1500[num_cols], as.numeric)

Pacing_splits_PHD_1500$Competition_date <- as.Date(
  as.numeric(Pacing_splits_PHD_1500$Competition_date), origin = "1899-12-30"
)

str(Pacing_splits_PHD_800)
str(Pacing_splits_PHD_1500)

# =============================================================================
# 1. Prepare the 800-1500m datasets
# =============================================================================

# 1.0. Build explicit keys once and keep them

add_keys <- function(df) {
  df %>%
    mutate(
      # Normalize Athletes and Nationality to ASCII for safe joins
      Athletes_clean   = stri_trans_general(Athletes, "Latin-ASCII"),
      Nationality_clean= stri_trans_general(Nationality, "Latin-ASCII"),
      
      # Build merge keys on clean versions
      athlete_key = paste(Athletes_clean, Nationality_clean, sep = "_"),
      season_id   = paste(athlete_key, Competition_year, sep = "_"),
      race_id     = Comp_ID,
      entry_id    = paste(athlete_key, race_id, sep = "_")
    )
}

Pacing_splits_PHD_800  <- add_keys(Pacing_splits_PHD_800)
Pacing_splits_PHD_1500 <- add_keys(Pacing_splits_PHD_1500)

Pacing_splits_PHD_800 <- subset(Pacing_splits_PHD_800, Pacing_splits_PHD_800$Competition_level=="World")
Pacing_splits_PHD_1500 <- subset(Pacing_splits_PHD_1500, Pacing_splits_PHD_1500$Competition_level=="World")

# 1.1. Compute speed in m/s

# General function: computes speed for any df with split columns like "0_400m"
compute_split_speeds <- function(df) {
  df %>%
    mutate(across(
      matches("^\\d+_\\d+m$"),
      .fns = ~ {
        # Extract distance from column name
        dist <- as.numeric(str_extract(cur_column(), "(?<=_)\\d+(?=m$)")) -
          as.numeric(str_extract(cur_column(), "^\\d+"))
        dist / .
      },
      .names = "speed_{.col}"
    ))
}

# Apply to all relevant datasets
Pacing_splits_PHD_800   <- compute_split_speeds(Pacing_splits_PHD_800)
Pacing_splits_PHD_1500  <- compute_split_speeds(Pacing_splits_PHD_1500)

# 1.2. Select variables of interest 

prepare_split_data <- function(df) {
  
  df <- df %>%
    group_by(race_id, Gender) %>%
    arrange(Final_time, .by_group = TRUE) %>%
    ungroup()
  
  return(df)
}

splits_800  <- prepare_split_data(Pacing_splits_PHD_800)
splits_1500 <- prepare_split_data(Pacing_splits_PHD_1500)

#splits_800 <- splits_800[!is.na(splits_800$`0_100m`),]

# 1.3. Long format transform function

transform_splits <- function(data) {
  meta_keep <- c("athlete_key","season_id","race_id","entry_id",
                 "Athletes","Gender","Nationality","Competition",
                 "Competition_type","Championship_step","Competition_level",
                 "Competition_year","Q","Final_time","Temperature_C",
                 "Humidity")
  
  data %>%
    select(any_of(meta_keep), matches("^\\d+_\\d+m$")) %>%
    tidyr::pivot_longer(
      cols = matches("^\\d+_\\d+m$"),
      names_to   = c("start_m","end_m"),
      names_pattern = "^(\\d+)_(\\d+)m$",
      values_to  = "split_time_s",
      values_drop_na = TRUE
    ) %>%
    mutate(
      start_m = as.numeric(start_m),
      end_m   = as.numeric(end_m),
      distance_m = end_m - start_m,
      split_speed_m_s = distance_m / split_time_s
    ) %>%
    group_by(race_id, entry_id) %>%
    arrange(start_m, .by_group = TRUE) %>%
    mutate(
      split_index     = dplyr::row_number(),
      cum_distance_m  = end_m,
      race_prop_raw   = cum_distance_m / max(cum_distance_m),
      race_prop       = pmin(pmax(race_prop_raw, 0), 1),
      cum_time_s      = cumsum(split_time_s)
    ) %>% 
    select(-race_prop_raw) %>%
    ungroup()
}

# Use
long_800  <- transform_splits(splits_800)
long_1500 <- transform_splits(splits_1500)

# 1.4. Compute a rank, position change function and position range

# Position change
compute_position_and_change <- function(df) {
  # Prefer race_id; fallback to Competition
  race_key <- if ("race_id" %in% names(df)) "race_id" else "Competition"
  
  df1 <- df %>%
    arrange(.data[[race_key]], entry_id, split_index) %>%
    group_by(.data[[race_key]], entry_id) %>%
    # deterministic cumulative time within entry
    mutate(cumulative_time = cumsum(split_time_s)) %>%
    ungroup() %>%
    # position at each split within race (lower cum time = better)
    group_by(.data[[race_key]], split_index) %>%
    mutate(position = rank(cumulative_time, ties.method = "min")) %>%
    ungroup() %>%
    # position change within entry across splits
    group_by(.data[[race_key]], entry_id) %>%
    arrange(split_index, .by_group = TRUE) %>%
    mutate(position_change = position - dplyr::lag(position, default = dplyr::first(position))) %>%
    ungroup()
  
  # final rank per entry = position at the race's last split
  last_idx <- df1 %>%
    group_by(.data[[race_key]]) %>%
    summarise(last_split_index = max(split_index, na.rm = TRUE), .groups = "drop")
  
  finals <- df1 %>%
    inner_join(last_idx, by = race_key) %>%
    filter(split_index == last_split_index) %>%
    select(!!race_key, entry_id, rank = position) %>%
    mutate(rank = as.integer(rank))
  
  df1 %>%
    left_join(finals, by = c(setNames(nm = race_key), "entry_id"))
}

long_800  <- compute_position_and_change(long_800)
long_1500 <- compute_position_and_change(long_1500)

# Position range
summary(long_800$rank) # max = 12
summary(long_1500$rank) # max = 16

get_range_rank <- function(rank, event) {
  cut_points <- switch(as.character(event),
                       "800"  = c(0, 3, 6, 9, 12),
                       "1500" = c(0, 3, 6, 9, 12, 16),
                       c(0, 8)) # fallback
  labels <- switch(as.character(event),
                   "800"  = c("1-3rd", "4-6th", "7-9th", "10-12th"),
                   "1500" = c("1-3rd", "4-6th", "7-9th", "10-12th", "13-16th"),
                   "Other")
  cut(rank, breaks = cut_points, labels = labels, include.lowest = TRUE, right = TRUE)
}

long_800$range_rank   <- get_range_rank(long_800$rank, "800")
long_1500$range_rank  <- get_range_rank(long_1500$rank, "1500")

# 1.5. Compute pacing behavior: CV and Speed 0-200m

long_800 <- long_800 %>%
  arrange(race_id, entry_id, split_index) %>%
  group_by(race_id, entry_id) %>%
  mutate(CV = (sd(split_speed_m_s[1:8]) / mean(split_speed_m_s[1:8]))*100,
         Speed_200 = mean(split_speed_m_s[1:2]))

long_1500 <- long_1500 %>%
  arrange(race_id, entry_id, split_index) %>%
  group_by(race_id, entry_id) %>%
  mutate(CV = (sd(split_speed_m_s[1:15]) / mean(split_speed_m_s[1:15]))*100,
         Speed_200 = mean(split_speed_m_s[1:2]))

# 1.6. Final merge and filtering

# Prepare SB table once
results_PL_WA_SB2 <- results_PL_WA_SB %>%
  mutate(
    # Normalize full ID to ASCII
    ID_clean   = stri_trans_general(ID, "Latin-ASCII"),
    athlete_key = sub("_(\\d{4})$", "", ID_clean),
    Competition_year = as.integer(sub("^.*_(\\d{4})$", "\\1", ID_clean)),
    season_id   = paste(athlete_key, Competition_year, sep = "_")
  ) %>%
  # Drop duplicates by season_id (keep the first occurrence)
  distinct(season_id, .keep_all = TRUE) %>%
  select(season_id, E, S, a, b, r_squared, residual_se, n_performances)

merge_params <- function(df, sb2 = results_PL_WA_SB2) {
  df %>%
    left_join(sb2, by = "season_id") %>%
    filter(!is.na(E))
}

long_800  <- merge_params(long_800, results_PL_WA_SB2)
long_1500 <- merge_params(long_1500, results_PL_WA_SB2)

# 1.7. Compute performance predictions of PL model

predict_PL <- function(df, dist) {
  df %>%
    mutate(!!paste0("pred_PL_", dist) := (a / dist)^(1 / (b - 1)))
}

long_800  <- predict_PL(long_800, 800)
long_1500 <- predict_PL(long_1500, 1500)

# 1.8. Flag entries with any split v > S

flag_overS_entries <- function(df, tol = 0.00) {
  df %>%
    group_by(race_id, entry_id) %>%
    summarise(
      S_entry   = first(S),
      any_overS = any(split_speed_m_s > S_entry * (1 + tol), na.rm = TRUE),
      n_overS   = sum(split_speed_m_s > S_entry * (1 + tol), na.rm = TRUE),
      max_vS    = max(split_speed_m_s / S_entry, na.rm = TRUE),
      .groups = "drop"
    )
}

# Apply to 800 & 1500 and remove bad entries

overS_800  <- flag_overS_entries(long_800,  tol = 0.00)
overS_1500 <- flag_overS_entries(long_1500, tol = 0.00)

bad_800  <- overS_800  %>% filter(any_overS) %>% select(race_id, entry_id)
bad_1500 <- overS_1500 %>% filter(any_overS) %>% select(race_id, entry_id)

# Keep a log for Methods/Limitations
excluded_overS_800  <- inner_join(long_800  %>% distinct(race_id, entry_id, Athletes, Gender, season_id),
                                  bad_800,  by = c("race_id","entry_id"))
excluded_overS_1500 <- inner_join(long_1500 %>% distinct(race_id, entry_id, Athletes, Gender, season_id),
                                  bad_1500, by = c("race_id","entry_id"))

# quick counts
nrow(excluded_overS_800)
nrow(excluded_overS_1500)

# Delete from database
long_800  <- anti_join(long_800,  bad_800,  by = c("race_id","entry_id"))
long_1500 <- anti_join(long_1500, bad_1500, by = c("race_id","entry_id"))

# =============================================================================
# 2. Estimation of the PL(f) parameters: goodness-of-fit and uncertainty analysis
# =============================================================================

# 2.1. Create quality fit datasets

pl_quality_fit_800 <- long_800 %>%
  select(season_id, E, S, a, b, r_squared, residual_se, n_performances)

pl_quality_fit_1500 <- long_1500 %>%
  select(season_id, E, S, a, b, r_squared, residual_se, n_performances)

pl_quality_fit_963 <- rbind(pl_quality_fit_800, pl_quality_fit_1500)  
pl_quality_fit_963$race_id <- NULL
pl_quality_fit_963$entry_id <- NULL
pl_quality_fit_963 <- pl_quality_fit_963 %>% distinct(season_id, .keep_all = TRUE)

# 2.2. Identify and remove implausible PL fits
# -----------------------------------------------------------------------------
# From the LMMs diagnostics, the conclusion was that some PL fits (estimated by S)
# were implausible from a physiological perspective. All those implausible fits
# come from a 2-point PL fitting.
# Instead of removing all the 2-point fits (which would weaken the upcoming analysis,
# mostly MC analysis by cutting the sample too hard), we decided to delete the 
# season_id that have a fitted S above the maximum of the 3-point PL fits. This 
# lands to ~25, which is in line with our first paper (Walt et al., 2025, IJSPP)'s 
# benchmark tables.
# ------------------------------------------------------------------------------

# Identify the implausible S
reliable_S <- pl_quality_fit_963 %>% filter(n_performances >= 3) %>% pull(S)
max(reliable_S) # 25.05
implausible_seasons <- pl_quality_fit_963 %>% filter(S > max(reliable_S)) %>% pull(season_id)
length(implausible_seasons)  # 12

# Remove the implausible S, sample = 951
pl_quality_fit_951 <- pl_quality_fit_963 %>% filter(!(season_id %in% implausible_seasons))
long_800_951  <- long_800  %>% filter(!(season_id %in% implausible_seasons))
long_1500_951 <- long_1500 %>% filter(!(season_id %in% implausible_seasons))

# 2.3. Goodness-of-fit quality: log-log regressions

# Quality fit summary (all n)
fit_quality_summary <- pl_quality_fit_951 %>%
  summarise(
    n_athletes = n(),
    mean_r2 = mean(r_squared, na.rm = TRUE), sd_r2 = sd(r_squared, na.rm = TRUE),
    mean_residual_se = mean(residual_se, na.rm = TRUE), sd_residual_se = sd(residual_se, na.rm = TRUE),
    mean_n = mean(n_performances), sd_n = sd(n_performances),
    min_n = min(n_performances), max_n = max(n_performances),
    pct_two_point = mean(n_performances == 2) * 100
  )

# Quality fit summary by n of SBs
fit_quality_by_n <- pl_quality_fit_951 %>%
  mutate(fit_group = if_else(n_performances == 2,
                             "n = 2 (perfect fit by construction)",
                             "n \u2265 3")) %>%
  group_by(fit_group) %>%
  summarise(
    n_athletes = n(),
    mean_r2 = mean(r_squared, na.rm = TRUE), sd_r2 = sd(r_squared, na.rm = TRUE),
    mean_residual_se = mean(residual_se, na.rm = TRUE), sd_residual_se = sd(residual_se, na.rm = TRUE),
    mean_n = mean(n_performances), min_n = min(n_performances), max_n = max(n_performances)
  )

# 2.4. PL parameters uncertainty

# Supplementary characterization of PL parameter uncertainty

# Create the dataset
pl_param_uncertain_800 <- long_800_951[!duplicated(long_800_951$season_id),]
pl_param_uncertain_800 <- pl_param_uncertain_800 %>% select(season_id, race_id, Athletes, Gender)

pl_param_uncertain_1500 <- long_1500_951[!duplicated(long_1500_951$season_id),]
pl_param_uncertain_1500 <- pl_param_uncertain_1500 %>% select(season_id, race_id, Athletes, Gender)

pl_param_uncertain <- rbind(pl_param_uncertain_800, pl_param_uncertain_1500)
pl_param_uncertain <- pl_param_uncertain[!duplicated(pl_param_uncertain$season_id),]

WA_SB_best_PL2 <- WA_SB_best_PL %>%
  mutate(
    # Normalize full ID to ASCII
    ID_clean   = stri_trans_general(ID, "Latin-ASCII"),
    athlete_key = sub("_(\\d{4})$", "", ID_clean),
    Competition_year = as.integer(sub("^.*_(\\d{4})$", "\\1", ID_clean)),
    season_id   = paste(athlete_key, Competition_year, sep = "_")) %>%
  select(-Athletes)

pl_param_uncertain <- left_join(pl_param_uncertain, WA_SB_best_PL2, by = c("season_id"))

# Analytic parameter precision per season (n>=3)
# SE(intercept) -> precision of S; SE(slope) -> precision of E
param_precision <- pl_param_uncertain %>%
  inner_join(pl_quality_fit_951 %>% filter(n_performances >= 3) %>% select(season_id), by = "season_id") %>%
  group_by(season_id) %>%
  group_modify(~ {
    fit <- lm(log_speed ~ log_time, data = .x)
    co  <- summary(fit)$coefficients
    tibble(
      n_performances = nrow(.x),
      S_full  = exp(co["(Intercept)", "Estimate"]),
      SE_logS = co["(Intercept)", "Std. Error"],
      E_full  = co["log_time", "Estimate"] + 1,
      SE_E    = co["log_time", "Std. Error"]
    )
  }) %>%
  ungroup()

param_precision %>%
  group_by(n_performances) %>%
  summarise(n = n(), median_SE_logS = median(SE_logS), median_SE_E = median(SE_E))

# Correlations with time

sb_separation <- pl_param_uncertain %>%  
  inner_join(pl_quality_fit_951 %>% filter(n_performances >= 3) %>% select(season_id), by = "season_id") %>%
  group_by(season_id) %>%
  summarise(
    n_performances = n(),
    log_time_range = max(log_time) - min(log_time),
    log_time_min_gap = min(diff(sort(log_time))),   # closest pair, in log-time - the key variable
    dist_ratio = max(distance) / min(distance)
  ) %>%
  ungroup()

param_precision_vs_separation <- param_precision %>% left_join(sb_separation[, c("season_id", "log_time_range", "log_time_min_gap", "dist_ratio")],
                                                      by = "season_id")

param_precision_vs_separation %>%
  filter(n_performances == 3) %>%
  summarise(rho_cv = cor(log_time_min_gap, SE_logS, method = "spearman"),
            rho_bias = cor(log_time_min_gap, SE_E, method = "spearman"))

param_precision_vs_separation %>%
  filter(n_performances == 3) %>%
  summarise(rho_cv = cor(log_time_range, SE_logS, method = "spearman"),
            rho_bias = cor(log_time_range, SE_E, method = "spearman"))

# Criterion (full n=3) model vs. a 2-point model using only the two most-separated
# distances (i.e., drop the single in-between distance) - same logic Kordi et al.
# (2018) used to validate their 2-point vs 3-point CS/D' models. n=3 only: this is
# where the separation analysis showed the real instability lives; n>=4 is already stable.
extremes_vs_full <- pl_param_uncertain %>%
  inner_join(pl_quality_fit_951 %>% filter(n_performances == 3) %>% select(season_id), by = "season_id") %>%
  group_by(season_id) %>%
  group_modify(~ {
    df <- .x %>% arrange(log_time)
    fit_full <- lm(log_speed ~ log_time, data = df)
    fit_2pt  <- lm(log_speed ~ log_time, data = df[c(1, nrow(df)), ])  # drop the middle distance
    tibble(
      S_full = exp(coef(fit_full)[1]), E_full = coef(fit_full)[2] + 1,
      S_2pt  = exp(coef(fit_2pt)[1]),  E_2pt  = coef(fit_2pt)[2] + 1
    )
  }) %>%
  ungroup() %>%
  mutate(S_bias_pct = 100 * (S_2pt - S_full) / S_full,
         E_bias_pct = 100 * (E_2pt - E_full) / E_full)

extremes_vs_full %>%
  summarise(n = n(),
            median_S_bias_pct = median(S_bias_pct), mean_S_bias_pct = mean(S_bias_pct),
            median_E_bias_pct = median(E_bias_pct), mean_E_bias_pct = mean(E_bias_pct))

# =============================================================================
# 3. Data description
# =============================================================================

# 3.0. Final datasets renaming

long_800_963  <- long_800 # to be used in Appendix
long_1500_963 <- long_1500 # to be used in Appendix

long_800  <- long_800  %>% filter(!(season_id %in% implausible_seasons))
long_1500 <- long_1500 %>% filter(!(season_id %in% implausible_seasons))

# 3.1. Number of athletes

unique_800 <- long_800[!duplicated(long_800$season_id),]
table(unique_800$Gender)
unique_800 <- unique_800 %>% select(season_id, race_id, Athletes, Gender)

unique_1500 <- long_1500[!duplicated(long_1500$season_id),]
table(unique_1500$Gender)
unique_1500 <- unique_1500 %>% select(season_id, race_id, Athletes, Gender)

# Total race performance: 1012 (532 F; 480 M)
unique_pacing <- rbind(unique_800, unique_1500)
table(unique_pacing$Gender)

# Total unique season-athlete: 951 unique athlete-season PL fits (488 F, 463 M season_ids) (possible multiple appearance within a year)
unique_pacing_2 <- unique_pacing[!duplicated(unique_pacing$season_id),]
table(unique_pacing_2$Gender)

races_per_season_800 <- long_800 %>%
  group_by(season_id) %>%
  summarise(n_races = n_distinct(race_id))

races_per_season_1500 <- long_1500 %>%
  group_by(season_id) %>%
  summarise(n_races = n_distinct(race_id))

# Total unique athlete: 495 (233 F, 262 M)  (possible multiple season_id)
unique_pacing_3 <- unique_pacing[!duplicated(unique_pacing$Athletes),]
table(unique_pacing_3$Gender)

# Delete athlete v > S
which(unique_pacing_2$season_id == "Rose Mary ALMANZA_CUB_2021") # Deleted for one race only
table(long_800$season_id == "Rose Mary ALMANZA_CUB_2021") # 24 rows (8 rows [splits] by race_id) = 3 races
long_800 %>%
  filter(season_id == "Rose Mary ALMANZA_CUB_2021") %>%
  summarise(n_races = n_distinct(race_id))

which(unique_pacing_2$season_id == "Cory Ann MCGEE_USA_2023") # Whole season-id deleted

# 3.2. Number of competition

races_800 <- long_800[!duplicated(long_800$race_id),]
table(races_800$Gender)

races_1500 <- long_1500[!duplicated(long_1500$race_id),]
table(races_1500$Gender)

# 3.3. Table 1: SB performances mean and SD 

WA_SB_best_PL2 <- WA_SB_best_PL %>%
  mutate(
    # Normalize full ID to ASCII
    ID_clean   = stri_trans_general(ID, "Latin-ASCII"),
    athlete_key = sub("_(\\d{4})$", "", ID_clean),
    Competition_year = as.integer(sub("^.*_(\\d{4})$", "\\1", ID_clean)),
    season_id   = paste(athlete_key, Competition_year, sep = "_")) %>%
  select(-Athletes)

unique_pacing_2 <- left_join(unique_pacing_2, WA_SB_best_PL2, by = c("season_id"))

Performance_distr <- unique_pacing_2 %>%
  group_by(Gender, distance, field) %>%
  summarise(
    n = n(),
    mean_sec = mean(time, na.rm = TRUE),
    sd = sd(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Break mean_sec into h:m:s
    hours   = floor(mean_sec / 3600),
    minutes = floor((mean_sec %% 3600) / 60),
    seconds = mean_sec %% 60,
    
    # Conditional formatting
    mean_formatted = case_when(
      distance %in% c(42195, 21095) ~ sprintf("%d:%02d:%02.0f", hours, minutes, seconds),
      TRUE ~ sprintf("%02d:%04.1f", minutes, seconds)
    ),
    
    # Decimal comma for non-marathon formats
    mean_formatted = ifelse(distance %in% c(42195, 21095),
                            mean_formatted,
                            str_replace(mean_formatted, "\\.", ",")),
    
    sd = round(sd, 1)
  ) %>%
  select(Gender, distance, field, n, mean_formatted, sd)

sum(Performance_distr$n)

# =============================================================================
# 4. Fatigue-intergrated Power Law (PLf) model computation
# =============================================================================

# 4.1. Raw fatigue computation

compute_basic_fatigue <- function(df, eps = 1e-8) {
  df %>%
    arrange(race_id, entry_id, split_index) %>%
    group_by(race_id, entry_id) %>%
    mutate(
      alpha  = 1 / (pmax(1 - E[1], eps)),
      v_i    = pmax(split_speed_m_s, 0),
      dt_i   = pmax(split_time_s,    0),
      fatigue_split     = (v_i / S[1])^alpha * dt_i,
      acc_fatigue_basic = cumsum(fatigue_split)
    ) %>%
    ungroup() %>%
    select(-alpha, -v_i, -dt_i)
}

long_800_fatigue <- compute_basic_fatigue(long_800)
long_1500_fatigue <- compute_basic_fatigue(long_1500)

# 4.2. Compute Drake (2024) Recursive PL Downshift
# Based on: https://github.com/jonahd99/PL_vs_CP/blob/main/rmd/figure_12.Rmd

compute_dynamic_fatigue_drake <- function(df, tactical_effect = 0.05, eps = 1e-8) {
  df <- df %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(race_id, entry_id, split_index)
  
  df %>%
    dplyr::group_by(race_id, entry_id) %>%
    dplyr::group_split() %>%
    lapply(function(g) {
      n <- nrow(g)
      S0 <- g$S[1]
      E <- g$E[1]
      alpha <- 1 / max(1 - E, eps)
      
      v  <- pmax(g$split_speed_m_s, 0)
      dt <- pmax(g$split_time_s,    0)
      v_alpha <- v^alpha
      
      S0 <- max(S0, v[1])
      Sdyn <- numeric(n)
      U <- numeric(n)
      Fd <- numeric(n)
      Facc <- numeric(n)
      
      Sdyn[1] <- S0
      U[1] <- S0^alpha
      
      for (i in 1:n) {
        Fd[i]  <- (if (Sdyn[i] > 0) (v[i] / Sdyn[i])^alpha else Inf) * dt[i]
        if (!is.finite(Fd[i])) Fd[i] <- NA_real_
        Facc[i] <- if (i == 1) Fd[i] else Facc[i-1] + Fd[i]
        
        if (i < n) {
          phi <- if (g$position_change[i] >= 2) { # lost ???2 positions 
            1 + tactical_effect
          } else if (g$position_change[i] <= -2) { # gained ???2 positions 
            1 - tactical_effect
          } else {
            1
          }
          
          tryU  <- U[i] - dt[i] * v_alpha[i] * phi
          nextU <- if (is.finite(tryU) && tryU > 0) min(tryU, U[i]) else v[i+1]^alpha
          
          S_next <- max(nextU^(1/alpha), v[i+1])
          if (S_next > Sdyn[i] && S_next - Sdyn[i] < 1e-8) S_next <- Sdyn[i]
          
          Sdyn[i+1] <- S_next
          U[i+1]    <- Sdyn[i+1]^alpha
        }
      }
      
      g$S_dyn_drake             <- Sdyn
      g$Fatigue_dyn_drake       <- Fd
      g$Fatigue_dyn_accum_drake <- Facc
      g
    }) %>% dplyr::bind_rows()
}

long_800_fatigue <- compute_dynamic_fatigue_drake(long_800_fatigue)
long_1500_fatigue <- compute_dynamic_fatigue_drake(long_1500_fatigue)

# Real-time predictions

predict_remaining_time <- function(S_dyn_drake, b, remaining_distance) {
  (S_dyn_drake / remaining_distance)^(1 / (b - 1))
}

long_800_fatigue <- long_800_fatigue %>%
  group_by(entry_id, race_id) %>%
  mutate(
    remaining_distance = 800 - split_index * 100,
    predicted_remaining_time = predict_remaining_time(S_dyn_drake, b, remaining_distance),
    predicted_total_time = cumsum(split_time_s) + predicted_remaining_time
  )

long_1500_fatigue <- long_1500_fatigue %>%
  group_by(entry_id, race_id) %>%
  mutate(
    remaining_distance = 1500 - split_index * 100,
    predicted_remaining_time = predict_remaining_time(S_dyn_drake, b, remaining_distance),
    predicted_total_time = cumsum(split_time_s) + predicted_remaining_time
  )

# Compute MARE for real-time predictions

# Compute performance error metrics between row 2 and n-1
compute_mare_full <- function(df, true_col, pred_cols) {
  df %>%
    group_by(entry_id, race_id) %>%
    group_modify(~ {
      .x_trimmed <- .x[2:(nrow(.x)-1), ]
      results <- lapply(pred_cols, function(pred_col) {
        actual <- .x_trimmed[[true_col]]
        predicted <- .x_trimmed[[pred_col]]
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        mare <- mean(abs(1 - predicted / actual), na.rm = TRUE)
        tibble(MAE = mae, MARE = mare)
      })
      bind_rows(results)
    }) %>%
    ungroup()
}

# Usage

MARE_800_real_time_pred <- compute_mare_full(long_800_fatigue, true_col = "Final_time", pred_cols = "predicted_total_time")
mean(MARE_800_real_time_pred$MARE)

MARE_1500_real_time_pred <- compute_mare_full(long_1500_fatigue, true_col = "Final_time", pred_cols = "predicted_total_time")
mean(MARE_1500_real_time_pred$MARE)

# 4.3. Compute curve shape metrics

compute_fatigue_features <- function(df, distance = 800) {
  # Add starting line row
  df_start <- df %>%
    group_by(race_id, entry_id) %>%
    slice(1) %>%
    mutate(
      race_prop = 0,
      split_index = 0,
      split_speed_m_s = NA,
      split_time_s = 0,
      Fatigue_dyn_drake = 0,
      Fatigue_dyn_accum_drake = 0,
      S_dyn_drake = first(S),
      risk_ratio = NA,
      risk_ratio_hazard = NA
    )
  
  df_full <- bind_rows(df_start, df) %>%
    arrange(race_id, entry_id, race_prop) %>%
    group_by(race_id, entry_id) %>%
    mutate(
      risk_ratio = split_speed_m_s / S_dyn_drake,
      risk_ratio_hazard = risk_ratio * (1 - race_prop)
    ) %>%
    ungroup()
  
  # Compute features per race_id x entry_id
  features <- df_full %>%
    group_by(race_id, entry_id, season_id, Athletes) %>%   # keep season_id & Athletes
    do({
      d <- .
      race_prop <- d$race_prop
      fatigue   <- d$Fatigue_dyn_accum_drake
      S_dyn     <- d$S_dyn_drake
      
      # Fatigue proportions
      max_fatigue <- max(fatigue, na.rm = TRUE)
      prop_low  <- mean(fatigue / max_fatigue <= 0.33, na.rm = TRUE)
      prop_mod  <- mean(fatigue / max_fatigue > 0.33 & fatigue / max_fatigue <= 0.66, na.rm = TRUE)
      prop_high <- mean(fatigue / max_fatigue > 0.66, na.rm = TRUE)
      
      # T50
      T50 <- tryCatch(
        approx(fatigue, race_prop, xout = 0.5 * max_fatigue)$y,
        error = function(e) NA_real_
      )
      
      # Slopes and curvature
      if (length(fatigue) > 2) {
        dy_dx <- diff(fatigue) / diff(race_prop)
        slope_max  <- max(dy_dx, na.rm = TRUE)
        slope_mean <- mean(dy_dx, na.rm = TRUE)
        d2y_dx2    <- diff(dy_dx) / diff(race_prop[-1])
        curvature  <- mean(d2y_dx2, na.rm = TRUE)
      } else {
        slope_max <- slope_mean <- curvature <- NA
      }
      
      # AUC
      auc_total <- pracma::trapz(race_prop, fatigue)
      auc_cumsum <- sapply(seq_along(fatigue),
                           function(i) pracma::trapz(race_prop[1:i], fatigue[1:i]))
      auc_cumsum_pct <- 100 * auc_cumsum / auc_total
      
      # S dynamics
      if (length(S_dyn) > 1 && !all(is.na(S_dyn))) {
        S_dyn_slope <- coef(lm(S_dyn ~ race_prop))[2]
        S_drop <- 100 * (first(S_dyn) - last(S_dyn)) / first(S_dyn)
      } else {
        S_dyn_slope <- S_drop <- NA
      }
      
      tibble(
        season_id = d$season_id[1],
        Athletes  = d$Athletes[1],
        T50 = T50,
        prop_low_fatigue = prop_low,
        prop_moderate_fatigue = prop_mod,
        prop_high_fatigue = prop_high,
        slope_max = slope_max,
        slope_mean = slope_mean,
        curvature = curvature,
        auc = auc_total,
        S_dyn_slope = S_dyn_slope,
        S_drop = S_drop
      )
    }) %>%
    ungroup()
  
  # Join back to full df
  df_full <- df_full %>%
    left_join(features, by = c("race_id", "entry_id", "season_id", "Athletes"))
  
  return(list(results = features, df_with_preds = df_full))
}

# Apply
res_800  <- compute_fatigue_features(long_800_fatigue, distance = 800)
df_features_800  <- res_800$results
long_800_fatigue <- res_800$df_with_preds

res_1500 <- compute_fatigue_features(long_1500_fatigue, distance = 1500)
df_features_1500  <- res_1500$results
long_1500_fatigue <- res_1500$df_with_preds

# Compute Mid-race and Last-lap Fatigue & S

event_mid_last <- list(
  '800'  = list(mid = 6, last = 5),
  '1500' = list(mid = 11, last = 12)
)

add_mid_last_metrics <- function(df, event_name, fatigue_var = "Fatigue_dyn_accum_drake", S_var = "S_dyn_drake", risk_var = "risk_ratio", pos_var = "position") {
  mid_idx  <- event_mid_last[[as.character(event_name)]]$mid
  last_idx <- event_mid_last[[as.character(event_name)]]$last
  df %>%
    group_by(race_id, entry_id) %>%
    mutate(
      Fatigue_mid = .data[[fatigue_var]][mid_idx],
      auc_mid = trapz(race_prop[1:mid_idx], .data[[fatigue_var]][1:mid_idx]),
      S_mid = .data[[S_var]][mid_idx],
      risk_ratio_mid = .data[[risk_var]][[mid_idx]],
      Fatigue_last_lap = .data[[fatigue_var]][last_idx],
      auc_last_lap = trapz(race_prop[1:last_idx], .data[[fatigue_var]][1:last_idx]),
      S_last_lap = .data[[S_var]][last_idx],
      risk_ratio_last_lap = .data[[risk_var]][[last_idx]],
      S_drop_last_lap = 100 * (first(.data[[S_var]]) - .data[[S_var]][last_idx]) / first(.data[[S_var]]),
      position_last_lap = .data[[pos_var]][[last_idx]]
    ) %>% ungroup()
}

# Apply
long_800_fatigue   <- add_mid_last_metrics(long_800_fatigue, "800")
long_1500_fatigue  <- add_mid_last_metrics(long_1500_fatigue, "1500")

# 4.4. Final preparations on datasets

# Joining with Features (Generalized)

str(long_800)

cols_to_join_static <- c("race_id", "Competition_type", "entry_id", "Gender", "E",
                         "Final_time", "rank", "range_rank", "CV", "Speed_200", 
                         "Fatigue_mid", "auc_mid", "S_mid", "risk_ratio_mid",
                         "Fatigue_last_lap", "auc_last_lap", "S_last_lap", "risk_ratio_last_lap",
                         "S_drop_last_lap", "position_last_lap")

events <- c("800", "1500")
df_features_list <- list(df_features_800, df_features_1500)
long_list <- list(long_800_fatigue, long_1500_fatigue)

for (i in seq_along(df_features_list)) {
  event <- events[i]
  # Build event-specific column name
  pred_col <- paste0("pred_PL_", event)
  cols_to_join <- c(cols_to_join_static, pred_col)
  
  # Do the join
  df_features_list[[i]] <- df_features_list[[i]] %>%
    left_join(long_list[[i]][, cols_to_join], by = c("race_id", "entry_id")) %>%
    mutate(
      event = event,
      Pred_PL = !!sym(pred_col) # Make a generic column for all events
    ) %>%
    distinct()
}

# Assign back
df_features_800   <- df_features_list[[1]]
df_features_1500  <- df_features_list[[2]]

# Bind togheter

df_features_analysis <- bind_rows(df_features_800, df_features_1500) %>%
  mutate(across(c(Gender, range_rank, event), as.factor))

# Gender df
df_features_800_male <- subset(df_features_800, df_features_800$Gender =="Male")
df_features_800_female <- subset(df_features_800, df_features_800$Gender =="Female")

df_features_1500_male <- subset(df_features_1500, df_features_1500$Gender =="Male")
df_features_1500_female <- subset(df_features_1500, df_features_1500$Gender =="Female")

# Combine long datasets

long_800_fatigue$event <- "800"
long_1500_fatigue$event <- "1500"

# Stack all events, keeping only final split (where S_drop computed)
all_fatigue <- bind_rows(
  long_800_fatigue, long_1500_fatigue) %>%
  group_by(race_id, entry_id) %>%
  ungroup()

all_fatigue$event <- factor(all_fatigue$event, levels = c("800", "1500"))
all_fatigue$Gender <- factor(all_fatigue$Gender, levels = c("Male", "Female"))

all_fatigue$Competition_type2 <- as.factor(ifelse(all_fatigue$Competition_type == "Meeting", "Meeting", "Championship"))

all_fatigue <- subset(all_fatigue, !is.na(all_fatigue$risk_ratio)) # only keep the "real" splits

# Re-check the fitting and sample size (to match section 2.)
pl_quality_fit2 <- all_fatigue %>% distinct(season_id, .keep_all = TRUE)

fit_quality_summary2 <- pl_quality_fit2 %>%
  summarise(
    n_athletes = n(),
    mean_r2 = mean(r_squared, na.rm = TRUE), sd_r2 = sd(r_squared, na.rm = TRUE),
    mean_residual_se = mean(residual_se, na.rm = TRUE), sd_residual_se = sd(residual_se, na.rm = TRUE),
    mean_n = mean(n_performances), sd_n = sd(n_performances),
    min_n = min(n_performances), max_n = max(n_performances),
    pct_two_point = mean(n_performances == 2) * 100
  )

fit_quality_by_2 <- pl_quality_fit2 %>%
  mutate(fit_group = if_else(n_performances == 2,
                             "n = 2 (perfect fit by construction)",
                             "n \u2265 3")) %>%
  group_by(fit_group) %>%
  summarise(
    n_athletes = n(),
    mean_r2 = mean(r_squared, na.rm = TRUE), sd_r2 = sd(r_squared, na.rm = TRUE),
    mean_residual_se = mean(residual_se, na.rm = TRUE), sd_residual_se = sd(residual_se, na.rm = TRUE),
    mean_n = mean(n_performances), min_n = min(n_performances), max_n = max(n_performances)
  )

# =============================================================================
# 5. S durability / S decline analysis  
# =============================================================================
# -----------------------------------------------------------------------------
# Goals: fatigue causes the decline of S drop, where does that come from?
# Variables: CV (= pacing), Gender (= sex differences), 
#            range_rank (= too aggressive to follow), risk_ratio (= too close to max),
#            Athletes (= ID variability), E (= endurance parameter as an interaction factor), 
#            Events, etc.
# Methods: multilevel models using times series, all event and gender together in the same model
# -----------------------------------------------------------------------------

# 5.0. Splitwise % drop

all_fatigue <- all_fatigue %>%
  group_by(race_id, entry_id) %>%
  mutate(
    S_fresh = first(S_dyn_drake),
    S_drop_splitwise = 100 * (S_fresh - S_dyn_drake) / S_fresh
  ) %>%
  ungroup()

# 5.1. Descriptive stats for all events

summary(long_800_fatigue$S_drop)
summary(long_1500_fatigue$S_drop)
sd(long_800_fatigue$S_drop)
sd(long_1500_fatigue$S_drop)
quantile(long_800_fatigue$S_drop, 0.05)
quantile(long_800_fatigue$S_drop, 0.95)
quantile(long_1500_fatigue$S_drop, 0.05)
quantile(long_1500_fatigue$S_drop, 0.95)

# Retrieve the minimum S_drop row(s) for 1500 m
min_Sdrop_1500 <- long_1500_fatigue %>%
  filter(S_drop == min(S_drop, na.rm = TRUE)) %>%
  select(season_id, race_id, Athletes, Gender, E, S_drop)

min_Sdrop_1500

max_Sdrop_1500 <- long_1500_fatigue %>%
  filter(S_drop == max(S_drop, na.rm = TRUE)) %>%
  select(season_id, race_id, Athletes, Gender, E, S_drop)

max_Sdrop_1500

# Summarize S drop by athlete/race
S_summary <- all_fatigue %>%
  group_by(race_id, entry_id, season_id, Competition_type, event, Gender, range_rank, Q, E, S_last_lap, position_last_lap, Final_time) %>%
  summarize(
    fresh_S = first(S),
    final_S = last(S_dyn_drake),
    S_drop_raw = first(S) - last(S_dyn_drake),
    S_drop_pct = 100 * (first(S) - last(S_dyn_drake)) / first(S)
  ) %>%
  ungroup()

# Total distinct athlete-seasons in S_summary
dplyr::n_distinct(S_summary$season_id)

desc_stats <- S_summary %>%
  group_by(event, Gender) %>%
  summarize(
    N = n(),
    mean_fresh_S = mean(fresh_S, na.rm = TRUE),
    mean_final_S = mean(final_S, na.rm = TRUE),
    mean_S_drop_raw = mean(S_drop_raw, na.rm = TRUE),
    sd_S_drop_raw = sd(S_drop_raw, na.rm = TRUE),
    mean_S_drop_pct = mean(S_drop_pct, na.rm = TRUE),
    sd_S_drop_pct = sd(S_drop_pct, na.rm = TRUE),
    median_S_drop_pct = median(S_drop_pct, na.rm = TRUE),
    min_S_drop_pct = min(S_drop_pct, na.rm = TRUE),
    max_S_drop_pct = max(S_drop_pct, na.rm = TRUE)
  ) %>%
  ungroup()

desc_stats_comp <- S_summary %>%
  group_by(event, Gender, Competition_type) %>%
  summarize(
    N = n(),
    mean_fresh_S = mean(fresh_S, na.rm = TRUE),
    mean_final_S = mean(final_S, na.rm = TRUE),
    mean_S_drop_raw = mean(S_drop_raw, na.rm = TRUE),
    sd_S_drop_raw = sd(S_drop_raw, na.rm = TRUE),
    mean_S_drop_pct = mean(S_drop_pct, na.rm = TRUE),
    sd_S_drop_pct = sd(S_drop_pct, na.rm = TRUE),
    median_S_drop_pct = median(S_drop_pct, na.rm = TRUE),
    min_S_drop_pct = min(S_drop_pct, na.rm = TRUE),
    max_S_drop_pct = max(S_drop_pct, na.rm = TRUE)
  ) %>%
  ungroup()

# Individual variability of S_drop (on season_id level)
vc_model_pct <- lmer(S_drop_pct ~ 1 + (1 | season_id), data = S_summary_sensitivity)
VarCorr(vc_model_pct)

vc <- as.data.frame(VarCorr(vc_model_pct))
between_season_var <- vc$vcov[vc$grp == "season_id"]
residual_var       <- vc$vcov[vc$grp == "Residual"]
icc               <- between_season_var / (between_season_var + residual_var)
between_season_sd <- sqrt(between_season_var)

between_season_sd   # the actual "inter-individual variability" number, in percentage points
icc                 # fraction of total variance that's between-athlete vs. race-to-race noise

# 5.2. Plots: Highlight variability between Athletes

fresh_vs_last_S <- all_fatigue %>%
  group_by(race_id, entry_id, event, Competition_type) %>%
  summarize(
    fresh_S = first(S),
    last_S  = last(S_dyn_drake),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(fresh_S, last_S),
    names_to = "S_type",
    values_to = "S_value"
  ) %>%
  mutate(
    S_type = factor(S_type, levels = c("fresh_S", "last_S"),
                    labels = c("Fresh S", "Final S")),
    event2 = factor(event, labels = c("800 m", "1500 m"))
  )

fresh_vs_last_S$Competition_type <- factor(
  fresh_vs_last_S$Competition_type,
  levels = c("Meeting", "Championship - Heat", "Championship - Semi-final", "Championship - Final"),
  labels = c("DL meeting", "MC heat", "MC semi-final", "MC final")
)

# 800m vs. 1500m plot
ggplot(fresh_vs_last_S, aes(x = S_type, y = S_value, group = entry_id)) +
  # Bar: mean per group (event2, S_type)
  stat_summary(
    aes(group = event2),
    fun = mean,
    geom = "bar",
    fill = "grey90", color = "grey80", width = 0.5, alpha = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  # Trajectory lines per athlete
  geom_line(aes(color = entry_id), size = 0.6, alpha = 0.5, show.legend = FALSE) +
  # Individual points
  geom_jitter(
    aes(color = entry_id), 
    size = 2, width = 0.08, alpha = 0.9, show.legend = FALSE
  ) +
  # Optionally: mean labels on bars
  stat_summary(
    aes(group = event2, label = round(..y.., 1)),
    fun = mean,
    geom = "text",
    vjust = -15, size = 4.5, family = "Times New Roman", fontface = "bold",
    color = "black"
  ) +
  facet_wrap(~ event2, nrow = 1, strip.position = "top") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 32)) +
  scale_x_discrete(labels = c("Fresh S", "Final S")) +
  labs(
    x = NULL,
    y = "Speed Parameter (S, m/s)",
    title = "Decrease in S: Fresh vs. End of Race by Event",
    subtitle = "Individual variability (points/lines) and group means (bars)"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 18, family = "Times New Roman"),
    plot.subtitle = element_text(size = 14, family = "Times New Roman"),
    strip.text = element_text(size = 16, family = "Times New Roman", face = "bold"),
    axis.title.x = element_text(size = 14, family = "Times New Roman", face = "bold"),
    axis.title.y = element_text(size = 14, family = "Times New Roman", face = "bold"),
    text = element_text(family = "Times New Roman", size = 14)
  )

# 800m vs. 1500m vs. competition type plot (clean greys + tomato mean)
ggplot(fresh_vs_last_S, aes(x = S_type, y = S_value, group = entry_id)) +
  stat_summary(
    aes(group = interaction(event2, Competition_type)),
    fun = mean,
    geom = "bar",
    fill = "grey90", color = "grey90", width = 0.5, alpha = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  # Individual athlete-season trajectories (grey)
  geom_line(color = "grey55", linewidth = 0.5, alpha = 0.25, show.legend = FALSE) +
  geom_point(color = "grey55", size = 1.3, alpha = 0.25, show.legend = FALSE) +
  
  # Mean trajectory per event x competition type (tomato)
  stat_summary(
    aes(group = interaction(event2, Competition_type)),
    fun = mean,
    geom = "line",
    color = "tomato",
    linewidth = 1.0,
    alpha = 0.95
  ) +
  stat_summary(
    aes(group = interaction(event2, Competition_type)),
    fun = mean,
    geom = "point",
    color = "tomato",
    size = 2.2,
    alpha = 0.95
  ) +
  
  # Mean labels 
  stat_summary(
    aes(group = interaction(event2, Competition_type), label = round(after_stat(y), 1)),
    fun = mean,
    geom = "text",
    vjust = -10, size = 4.2, family = "Times New Roman", fontface = "bold",
    color = "black"
  ) +
  
  # Forces x-axis on both 800 and 1500m
  ggh4x::facet_grid2(
    event2 ~ Competition_type,
    scales = "free_x", space = "free_x",
    labeller = labeller(Competition_type = label_value, event2 = label_value),
    switch = "y",
    axes = "x"   # <- forces x-axis on both rows
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) +
  labs(
    x = NULL,
    y = "Speed Parameter (S, m/s)",
    title = "Decrease in S across competition contexts",
    subtitle = "Fresh vs. final S by event and competition type"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none",
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.line.x  = element_line()
  )

# Compute mean S trajectory per group
S_traj_summary <- all_fatigue %>%
  filter(!is.na(S_dyn_drake), split_index > 0) %>%
  group_by(event, Competition_type2, Gender, race_prop) %>%
  summarise(S_mean = mean(S_dyn_drake, na.rm = TRUE), .groups = "drop")

# Plot
ggplot() +
  # Individual trajectories (grey)
  geom_line(
    data = all_fatigue %>% filter(!is.na(S_dyn_drake), split_index > 0),
    aes(x = race_prop, y = S_dyn_drake, group = entry_id),
    color = "grey60", linewidth = 0.3, alpha = 0.15
  ) +
  # Mean trajectory (tomato)
  geom_line(
    data = S_traj_summary,
    aes(x = race_prop, y = S_mean),
    color = "tomato", linewidth = 1.2
  ) +
  geom_point(
    data = S_traj_summary,
    aes(x = race_prop, y = S_mean),
    color = "tomato", size = 2.0
  ) +
  # Fresh S reference line per entry (horizontal dashed)
  geom_hline(
    data = all_fatigue %>%
      filter(split_index == 1) %>%
      group_by(event, Competition_type2, Gender) %>%
      summarise(S_fresh = mean(S, na.rm = TRUE), .groups = "drop"),
    aes(yintercept = S_fresh),
    linetype = "dashed", color = "grey30", linewidth = 0.6
  ) +
  facet_grid(Gender ~ event + Competition_type2) +
  scale_x_continuous(
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_y_continuous(limits = c(5, 30)) +
  labs(
    x = "Race progression (%)",
    y = "Speed parameter S (m/s)",
    title = "Real-time decline of the speed parameter (S) during racing",
    subtitle = "Individual trajectories (grey) and group mean (red) by event and competition context"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    panel.spacing = unit(1, "lines")
  )

# 5.3. Compute LMM models

all_fatigue$Competition_type <- as.factor(all_fatigue$Competition_type)

# Step 1: Deciding the fixed effect structure

# Models 1 and 2: with and without weathers metrics
common_data <- all_fatigue %>% filter(!is.na(Temperature_C), !is.na(Humidity))  # n=894 for both

model_S_dyn_drake_1A <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event + Temperature_C + Humidity +
                              (1 + race_prop * risk_ratio | season_id), data = common_data, control = lmerControl(optimizer = "bobyqa"))

model_S_dyn_drake_2A <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop * risk_ratio | season_id), data = common_data, control = lmerControl(optimizer = "bobyqa"))

anova(model_S_dyn_drake_1A, model_S_dyn_drake_2A) 
# Weather variables: inclusion did not improve fit (??²(2)=2.36, p=.31) and both AIC and BIC favored their exclusion
# Weather variables were therefore excluded.

# Fit model 2 on full dataset (published model, Table 2)
model_S_dyn_drake <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event + 
                              (1 + race_prop * risk_ratio | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

summary(model_S_dyn_drake)
isSingular(model_S_dyn_drake) # No singularity issue

# Step 2: Deciding the random effects structure 
# Modelling alternating RE structure
model_S_dyn_drake_3 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop * risk_ratio * race_prop:risk_ratio | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_S_dyn_drake_4 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop * risk_ratio * race_prop:risk_ratio || season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_S_dyn_drake_5 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop + risk_ratio | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_S_dyn_drake_6 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop + risk_ratio || season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_S_dyn_drake_7 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_S_dyn_drake_8 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
                              (1 | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

anova(model_S_dyn_drake, model_S_dyn_drake_3)  # => the exact same
anova(model_S_dyn_drake, model_S_dyn_drake_4)  # => model_S_dyn_drake statistically better (p < 0.001)
anova(model_S_dyn_drake, model_S_dyn_drake_5)  # => model_S_dyn_drake statistically better (p < 0.001)
anova(model_S_dyn_drake, model_S_dyn_drake_6)  # => model_S_dyn_drake statistically better (p < 0.001)
anova(model_S_dyn_drake, model_S_dyn_drake_7)  # => model_S_dyn_drake statistically better (p < 0.001)
anova(model_S_dyn_drake, model_S_dyn_drake_8)  # => model_S_dyn_drake statistically better (p < 0.001)

# Conclusion: model_S_dyn_drake is statistically the best model.

# 5.4. LMMs checks: singularity, convergence, residual diagnostics and influence

# A) Singularity check
isSingular(model_S_dyn_drake)

# B) Convergence issues
model_S_dyn_drake@optinfo$conv$lme4$messages

# C) Residual diagnostics
res <- residuals(model_S_dyn_drake, type = "pearson")     
resid_std <- resid(model_S_dyn_drake) / sigma(model_S_dyn_drake)
summary(res); summary(resid_std)

plot(model_S_dyn_drake) # residuals vs fitted: check homoscedasticity
qqnorm(resid(model_S_dyn_drake)); qqline(resid(model_S_dyn_drake))  # normality
lattice::qqmath(ranef(model_S_dyn_drake)) # normality of the random effects themselves

length(which(abs(resid_std) > 2)); mean(abs(resid_std) > 2) * 100   # observed vs. nominal 4.55%
length(which(abs(resid_std) > 3)); mean(abs(resid_std) > 3) * 100   # observed vs. nominal 0.27%
car::vif(model_S_dyn_drake)

all_fatigue$resid_S <- resid(model_S_dyn_drake)
all_fatigue$resid_S2 <- all_fatigue$resid_S^2

ggplot(all_fatigue, aes(race_prop, resid_S2)) + geom_smooth(method = "loess") + labs(y = "Squared residual")
ggplot(all_fatigue, aes(risk_ratio, resid_S2)) + geom_smooth(method = "loess") + labs(y = "Squared residual")
ggplot(all_fatigue, aes(Competition_type, resid_S2)) + stat_summary() + labs(y = "Mean squared residual")
ggplot(all_fatigue, aes(event, resid_S2)) + stat_summary() + labs(y = "Mean squared residual")

# D) Influence: screen for candidate groups
all_fatigue$resid_S <- resid(model_S_dyn_drake)
flagged_seasons <- all_fatigue %>% filter(abs(resid_S) > 1) %>% distinct(season_id) %>% pull(season_id)
length(flagged_seasons)

all_fatigue %>% filter(abs(resid_S) > 1) %>%
  select(season_id, race_prop, risk_ratio, Competition_type, event) %>%
  print(n = 30)

# Leave-one-group-out refit on flagged candidates
targeted_influence <- purrr::map_dfr(flagged_seasons, function(sid) {
  df_sub <- all_fatigue %>% filter(season_id != sid)
  m <- lmer(S_dyn_drake ~ race_prop * risk_ratio * C(Competition_type, base = 4) + CV + Speed_200 + Gender + event +
              (1 + race_prop * risk_ratio | season_id), data = df_sub,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  corr_mat <- attr(VarCorr(m)$season_id, "correlation")
  tibble(season_id = sid,
         coef_race_prop_risk = fixef(m)["race_prop:risk_ratio"],
         corr_int_raceprop   = corr_mat["(Intercept)", "race_prop"],
         corr_int_risk       = corr_mat["(Intercept)", "risk_ratio"],
         corr_risk_inter     = corr_mat["risk_ratio", "race_prop:risk_ratio"],
         singular = isSingular(m))
})

# Interpretation: does removal move the FE
full_coef <- fixef(model_S_dyn_drake)["race_prop:risk_ratio"]
full_se   <- summary(model_S_dyn_drake)$coefficients["race_prop:risk_ratio", "Std. Error"]
targeted_influence <- targeted_influence %>%
  mutate(dfbeta_like = (coef_race_prop_risk - full_coef) / full_se)  # shift in SE units, DFBETAS-style
range(targeted_influence$coef_race_prop_risk); full_coef
summary(targeted_influence$dfbeta_like)   # values near 0 = negligible influence on this coefficient

# Interpretation: does removal move the RE
table(targeted_influence$singular)
summary(targeted_influence$corr_int_raceprop)
summary(targeted_influence$corr_int_risk)
summary(targeted_influence$corr_risk_inter)
targeted_influence %>% arrange(corr_int_risk) %>% print(n = nrow(targeted_influence))

# 5.5. Fixed effect analysis

# Releveling
all_fatigue$Competition_type3 <- factor(all_fatigue$Competition_type, ordered = FALSE)
all_fatigue$Competition_type3 <- relevel(all_fatigue$Competition_type3, ref = "Meeting")

model_S_dyn_drake_relevel <- lmer(S_dyn_drake ~ race_prop * risk_ratio * Competition_type3 + CV + Speed_200 + Gender + event +
                                    (1 + race_prop * risk_ratio | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

summary(model_S_dyn_drake_relevel)   # sanity check against model_S_dyn_drake before trusting emtrends

# Question 1: As the race progresses, how fast does S decline, in each competition type?
# Simple slope of race_prop by Competition_type, at representative risk_ratio values
# => Q1 averaged over three risk levels (0.4, 0.7, 1.0)
emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, var = "race_prop", pbkrtest.limit = 28647,
         at = list(risk_ratio = c(0.4, 0.7, 1.0)))

emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, var = "race_prop", pbkrtest.limit = 28647,
         at = list(risk_ratio = c(0.7)))

# => Question 1 separately at risk=0.5, risk=0.7, and risk=0.9
# With formal pairwise test
emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, pbkrtest.limit = 28647, var = "race_prop",
         at = list(risk_ratio = 0.5))
emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, pbkrtest.limit = 28647, var = "race_prop",
         at = list(risk_ratio = 0.7))
emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, pbkrtest.limit = 28647, var = "race_prop",
         at = list(risk_ratio = 0.9))

# Question 2: How much does taking on one more unit of risk cost you, in extra S decline, in each competition type?
# Simple slope of risk_ratio by Competition_type, at representative race_prop values
emtrends(model_S_dyn_drake_relevel, pairwise ~ Competition_type3, var = "risk_ratio", pbkrtest.limit = 28647, 
         at = list(race_prop = c(0.3, 0.6, 0.9)))

# Question 3: the risk_ratio slope at four race-progression points, separately per competition type
emtrends(model_S_dyn_drake_relevel, ~ race_prop | Competition_type3, var = "risk_ratio",
         at = list(race_prop = c(0.25, 0.5, 0.75, 0.9)), pbkrtest.limit = 28647)

emtrends(model_S_dyn_drake_relevel, pairwise ~ race_prop | Competition_type3, var = "risk_ratio",
         at = list(race_prop = c(0.25, 0.5, 0.75, 0.9)), pbkrtest.limit = 28647)

# Question 4: Is pacing the common root of the changes in S?
emmeans(model_S_dyn_drake_relevel, ~ race_prop, at = list(race_prop = c(0.25, 0.5, 0.75, 0.9)))
emmeans(model_S_dyn_drake_relevel, pairwise ~ race_prop, at = list(race_prop = c(0, 0.25, 0.5, 0.75, 0.9)))
emmeans(model_S_dyn_drake_relevel, pairwise ~ race_prop | Competition_type3, at = list(race_prop = c(0, 0.25, 0.5, 0.75, 0.9)))

# 5.6. Random effect analysis

# VarCorr(): this function calculates the estimated variances, standard deviations, and correlations between the random-effects
VarCorr(model_S_dyn_drake)
corr_mat_full <- attr(VarCorr(model_S_dyn_drake)$season_id, "correlation")
corr_mat_ful
sigma(model_S_dyn_drake)

# Signs checking
re <- as.data.frame(ranef(model_S_dyn_drake)$season_id)
re$season_id <- rownames(ranef(model_S_dyn_drake)$season_id)
fe <- fixef(model_S_dyn_drake)
r_ref <- median(all_fatigue$risk_ratio, na.rm = TRUE)   # or your chosen representative value

re <- re %>%
  mutate(
    total_intercept  = fe["(Intercept)"] + `(Intercept)`,
    realistic_slope   = (fe["race_prop"] + race_prop) +
      (fe["race_prop:risk_ratio"] + `race_prop:risk_ratio`) * r_ref
  )
cor.test(re$total_intercept, re$realistic_slope)
plot(re$total_intercept, re$realistic_slope)

# Extract random effects with season_id as explicit col
RF_S_dyn <- ranef(model_S_dyn_drake)$season_id %>%
  tibble::rownames_to_column("season_id")

summary(RF_S_dyn)

# Or tidy version already does it:
re_tidy <- broom.mixed::tidy(
  model_S_dyn_drake,
  effects = "ran_vals",
  conf.int = TRUE
) %>%
  filter(group == "season_id")  # athlete-season grouping

# Event lookup
event_lookup <- all_fatigue %>%
  group_by(season_id, event) %>%
  summarise(
    E = mean(E, na.rm = TRUE),
    S = mean(S, na.rm = TRUE),
    S_drop = mean(S_drop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(event2 = factor(event,
                         levels = c("800","1500"),
                         labels = c("800 m","1500 m")))

# Join event info to random effects
re_tidy <- re_tidy %>%
  left_join(event_lookup, by = c("level" = "season_id"))

# Correlation between RE (random effects) and E, by term
cor_results <- re_tidy %>%
  group_by(event2, term) %>%
  summarise(
    cor = cor(estimate, E, use = "complete.obs"),
    p   = cor.test(estimate, E)$p.value,
    n   = n(),
    .groups = "drop"
  )

cor_results

# Relabel terms and events with case_when()
re_plot <- re_tidy %>%
  mutate(
    term_lab = case_when(
      term == "(Intercept)"          ~ "Intercept (baseline S)",
      term == "race_prop"            ~ "Slope: race_prop",
      term == "risk_ratio"           ~ "Slope: risk ratio",
      term == "race_prop:risk_ratio" ~ "Interaction: race_prop × risk",
      TRUE                           ~ term
    ),
    event_lab = case_when(
      event == "800"  ~ "800 m",
      event == "1500" ~ "1500 m",
      TRUE            ~ as.character(event)
    )
  )

# Compute r and p per panel for the annotations
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001",
                            paste0("= ", formatC(p, format = "f", digits = 3)))

cor_annot <- re_plot %>%
  group_by(event_lab, term_lab) %>%
  summarise(
    r  = cor(estimate, E, use = "complete.obs"),
    p  = cor.test(estimate, E)$p.value,
    # where to place the label inside each panel
    x  = quantile(E, 0.15, na.rm = TRUE),
    y  = quantile(estimate, 0.85, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", sprintf("%.2f", r), "\n", "p ", fmt_p(p)))

# Scatter + lm line with r/p annotation, faceted by Event × RE term
ggplot(re_plot, aes(x = E, y = estimate)) +
  geom_point(alpha = 0.45, size = 1.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black", fill = "grey70") +
  geom_text(data = cor_annot,
            aes(x = x, y = y, label = label),
            family = "Times New Roman", size = 4.2, hjust = 1.5) +
  facet_grid(event_lab ~ term_lab, scales = "free_y") +
  labs(
    x = "Endurance parameter (E)",
    y = "Random-effect value (BLUP)",
    title = "Association between endurance (E) and random effects of the S model",
    subtitle = "Points = athlete-seasons; line = OLS fit; shaded = 95% CI; r and p per panel"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

# =============================================================================
# 6. Risk ratio analysis
# =============================================================================
# -----------------------------------------------------------------------------
# Goals: What, how, when is the risk ratio driven by?
# A question: Should the faster athletes (higher S) start faster? How much?
# Variables: CV (= pacing), Gender (= sex differences), 
#            position, proportion of the race,
#            Athletes (= ID variability), E (= endurance parameter as an interaction factor), 
#            S (= speed parameter), events, etc.
# Methods: multilevel models using times series, all event and gender together in the same model
# -----------------------------------------------------------------------------

# 6.1. Compute LMM models

# Step 1: Decide the fixed effect structure

# As for S_dyn model, weather variables are removed.

# Models 1 and 2: with and without PL parameters

model_risk_1 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E + S + event +
                     (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

model_risk_1a <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E + event +
                       (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

model_risk_1b <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + S + event +
                       (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

model_risk_2 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + event +
                       (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

anova(model_risk_1, model_risk_2) # model_risk_1 is statistically better => PL parameters should be included in the risk model
anova(model_risk_1, model_risk_1a) 
anova(model_risk_1, model_risk_1b)
anova(model_risk_1a, model_risk_1b)

# Compute the residuals of E/S

cor.test(all_fatigue$S, all_fatigue$E)
all_fatigue$E_resid <- resid(lm(E ~ S, data = all_fatigue))
cor.test(all_fatigue$E_resid, all_fatigue$E)

model_risk_1c <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + event +
                        (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

model_risk_1d <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + S + E_resid + event +
                        (1 + race_prop | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

anova(model_risk_1, model_risk_1c)
anova(model_risk_1, model_risk_1d) # Identical from a AIC/BIC point of view -> VIF to decide

car::vif(model_risk_1)
car::vif(model_risk_1d) # VIF is better with model_risk_1d

# Step 2: Deciding the random effects structure 
# Modelling alternating RE structure

model_risk_3 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
                     (1 + race_prop || season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_risk_4 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
                       (1 + race_prop + position | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_risk_5 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
                        (1 + race_prop + position || season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))
model_risk_6 <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
                       (1 | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

anova(model_risk_1, model_risk_3)  # => model_risk_1 is statistically better (p < 0.001)
anova(model_risk_1, model_risk_4)  # => model_risk_4 is statistically better (p < 0.001) 
anova(model_risk_4, model_risk_5)  # => model_risk_4 is statistically better (p < 0.001)
anova(model_risk_4, model_risk_6)  # => model_risk_4 is statistically better (p < 0.001)

# Conclusion: model_risk_4 RE is statistically the best.
# Fit model_risk_4 as "model_risk" (published model, Table 3)
model_risk <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
                       (1 + race_prop + position | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

summary(model_risk)

# 6.2. LMMs checks: singularity, convergence, residual diagnostics and influence

# A) Singularity check
isSingular(model_risk)

# B) Convergence issues
model_risk@optinfo$conv$lme4$messages

# C) Residual diagnostics
res_risk <- residuals(model_risk, type = "pearson")     
resid_std_risk <- resid(model_risk) / sigma(model_risk)
summary(res_risk); summary(resid_std_risk)

plot(model_risk) # residuals vs fitted: check homoscedasticity
qqnorm(resid(model_risk)); qqline(resid(model_risk))  # normality
lattice::qqmath(ranef(model_risk)) # normality of the random effects themselves

length(which(abs(resid_std_risk) > 2)); mean(abs(resid_std_risk) > 2) * 100 # observed vs. nominal 4.55%
length(which(abs(resid_std_risk) > 3)); mean(abs(resid_std_risk) > 3) * 100 # observed vs. nominal 0.27%

car::vif(model_risk)

all_fatigue$resid_risk <- resid(model_risk)
all_fatigue$resid_risk2 <- all_fatigue$resid_risk^2

ggplot(all_fatigue, aes(race_prop, resid_risk2)) + geom_smooth(method = "loess") + labs(y = "Squared residual")
ggplot(all_fatigue, aes(Competition_type, resid_risk2)) + stat_summary() + labs(y = "Mean squared residual")
ggplot(all_fatigue, aes(event, resid_risk2)) + stat_summary() + labs(y = "Mean squared residual")
ggplot(all_fatigue, aes(risk_ratio, resid_risk2)) + geom_smooth(method = "loess") + labs(y = "Squared residual")

# D) Influence: screen for candidate groups
flagged_seasons_risk <- all_fatigue %>%
  filter(risk_ratio < 0.999, abs(resid_std_risk) > 5) %>%   # excludes floor-bound splits before flagging
  distinct(season_id) %>% pull(season_id)
length(flagged_seasons_risk)

# Leave-one-group-out refit on flagged candidates
targeted_influence_risk <- purrr::map_dfr(flagged_seasons_risk, function(sid) {
  df_sub <- all_fatigue %>% filter(season_id != sid)
  m <- lmer(risk_ratio ~ race_prop * C(Competition_type, base = 4) * position + CV + Speed_200 + Gender + E_resid + S + event +
              (1 + race_prop + position | season_id), data = df_sub,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  corr_mat <- attr(VarCorr(m)$season_id, "correlation")
  tibble(season_id = sid,
         coef_race_prop_position = fixef(m)["race_prop:position"],
         corr_int_raceprop       = corr_mat["(Intercept)", "race_prop"],
         corr_int_position       = corr_mat["(Intercept)", "position"],
         corr_raceprop_position  = corr_mat["race_prop", "position"],
         singular = isSingular(m))
})

# Interpretation: does removal move the FE
full_coef <- fixef(model_risk)["race_prop:position"]
full_se   <- summary(model_risk)$coefficients["race_prop:position", "Std. Error"]
targeted_influence_risk <- targeted_influence_risk %>%
  mutate(dfbeta_like = (coef_race_prop_position - full_coef) / full_se)  # shift in SE units, DFBETAS-style
range(targeted_influence_risk$coef_race_prop_position); full_coef
summary(targeted_influence_risk$dfbeta_like)   # values near 0 = negligible influence on this coefficient

# Interpretation: does removal move the RE
table(targeted_influence_risk$singular)
summary(targeted_influence_risk$corr_int_raceprop)
summary(targeted_influence_risk$corr_int_position)
summary(targeted_influence_risk$corr_raceprop_position)
targeted_influence_risk %>% arrange(corr_int_raceprop) %>% print(n = nrow(targeted_influence_risk))

# 6.3. Fixed effect analysis
# -----------------------------------------------------------------------------
# Conditional slopes for the headline term (race_prop:position)
# -----------------------------------------------------------------------------

# Releveling
model_risk_relevel <- lmer(risk_ratio ~ race_prop * Competition_type3 * position + CV + Speed_200 + Gender + E_resid + S + event +
       (1 + race_prop + position | season_id), data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

# Sanity check: this refit should be numerically identical to model_risk
summary(model_risk_relevel) 
anova(model_risk, model_risk_relevel)

# Check the actual range of position before picking "at" values
summary(all_fatigue$position)

# How much risk ratio differ by type
emmeans(model_risk_relevel, pairwise ~ Competition_type3)

# Simple slopes of position on risk_ratio, at representative race_prop values
# how much does predicted risk_ratio change for a one-unit change in position, holding everything else fixed
position_trends <- emtrends(model_risk_relevel, ~ race_prop, var = "position",
                            at = list(race_prop = c(0.25, 0.50, 0.75, 0.90)))
position_trends

# Simple slopes of race_prop on risk_ratio, at representative position values
# how much does risk-taking ramp up as the race progresses, and does that ramp-up depend on position
raceprop_trends <- emtrends(model_risk_relevel, ~ position, var = "race_prop",
                            at = list(position = c(1, 4, 8)))
raceprop_trends

# Slopes split by Competition_type: three-way race_prop:Competition_type:position interaction
position_trends_by_comp <- emtrends(model_risk_relevel, ~ Competition_type3 | race_prop,
                                    var = "position",
                                    at = list(race_prop = c(0.25, 0.50, 0.75, 0.90)))
position_trends_by_comp

# Formal test of whether the position-slope differs across competition types,
# at each race_prop value (pairwise contrasts, not just eyeballing separate CIs)
pairs(position_trends_by_comp)

# Simple slopes of race_prop on risk_ratio, at specific competition context values
race_prop_by_comp <- emtrends(model_risk_relevel, ~ Competition_type3, var = "race_prop")
race_prop_by_comp
pairs(race_prop_by_comp)

# Actual marginal position-trend across race_prop
emtrends(model_risk_relevel, ~ 1, var = "position", at = list(race_prop = seq(0, 1, 0.05)))
emtrends(model_risk_relevel, ~ race_prop, var = "position", at = list(race_prop = seq(0, 1, 0.05)))

# 6.4. Random effect analysis

# VarCorr(): this function calculates the estimated variances, standard deviations, and correlations between the random-effects
VarCorr(model_risk)

vc_risk <- as.data.frame(VarCorr(model_risk))
print(vc_risk)

corr_mat_risk_full <- attr(VarCorr(model_risk)$season_id, "correlation")
corr_mat_risk_full

sigma(model_risk)

# Signs checking
# Extract season-level random effects for model_risk
re_risk <- ranef(model_risk)$season_id

# Realized (BLUP-based) correlation between intercept and position slope
cor.test(re_risk$`(Intercept)`, re_risk$position)

# 6.5. Plots

long_800_fatigue %>%
  filter(entry_id == "Emmanuel WANYONYI_KEN_Diamond League Lausanne 2024__45526_Male") %>%
  ggplot(aes(x = race_prop, y = risk_ratio)) +
  geom_line() +
  labs(title = "Risk Ratio by Split", x = "Race Proportion", y = "Risk Ratio")

# =============================================================================
# 7. Performance relation to Fresh S, S durability and position
# =============================================================================
# -----------------------------------------------------------------------------
# Goals: Do top finishers start with a higher "fresh S" (speed parameter) than 
#        others, or do they just "manage" their S better?
# Methods: Logisitic regressions, bc MC vs. DL.
# -----------------------------------------------------------------------------

# 7.0. New variables computation

# New variables computations
S_summary$Top3 <- as.factor(ifelse(S_summary$range_rank == "1-3rd", "Yes", "No"))
S_summary$S_drop_neg <- -S_summary$S_drop_raw
S_summary$position_neg <- -S_summary$position_last_lap

# 7.1. Fitting the regressions
# -----------------------------------------------------------------------------
# DV: Top3 finish (yes/no). Predictors: fresh_S, S_drop_neg, position_last_lap.
# glmer() (not glm()): season_id repeats across race entries within each
# event x gender stratum.
#
# "Remaining S" (S at the last lap) is deliberately excluded: it is
# arithmetically determined by fresh_S and S_drop_neg (remaining_S ~=
# fresh_S + S_drop_neg), so including all three is structural collinearity,
# not a manageable correlation -- confirmed by needing separate models by
# predictor pair to get sane VIFs.
# -----------------------------------------------------------------------------

# Data, stratified by event x gender

data_list <- list(
  "800_Male"    = S_summary %>% filter(event == "800",  Gender == "Male"),
  "800_Female"  = S_summary %>% filter(event == "800",  Gender == "Female") %>%
    mutate(S_drop_neg_c = S_drop_neg - mean(S_drop_neg)),
  "1500_Male"   = S_summary %>% filter(event == "1500", Gender == "Male")  %>%
    mutate(fresh_S_c = fresh_S - mean(fresh_S)),
  "1500_Female" = S_summary %>% filter(event == "1500", Gender == "Female")
)

# Base models (all three predictors, linear)

base_models <- purrr::map(data_list, ~
                            glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                                  family = "binomial", data = .x,
                                  control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))
)

purrr::map(base_models, car::vif)

# Linearity of the logit
# -----------------------------------------------------------------------------
# Established previously (2-predictor models, before position_neg was
# added):
#   - 800_Female: S_drop_neg nonlinear (smooth, accelerating near 0)
#   - 1500_Male:  fresh_S nonlinear (peak near ~11.5; survives exclusion of
#                 the sparse fresh_S>18 tail: p=0.0006 trimmed vs p=0.00016 full)
# Everything else showed no defensible nonlinearity.
# -----------------------------------------------------------------------------

test_predictor_linearity <- function(linear_model) {
  spl_freshS <- update(linear_model, . ~ . - fresh_S + ns(fresh_S, df = 3))
  spl_Sdrop  <- update(linear_model, . ~ . - S_drop_neg + ns(S_drop_neg, df = 3))
  lrt <- function(full, reduced) {
    stat <- as.numeric(2 * (logLik(full) - logLik(reduced)))
    dfd  <- attr(logLik(full), "df") - attr(logLik(reduced), "df")
    pchisq(stat, df = dfd, lower.tail = FALSE)
  }
  tibble::tibble(fresh_S_p = lrt(spl_freshS, linear_model),
                 S_drop_neg_p = lrt(spl_Sdrop, linear_model))
}

test_predictor_linearity(base_models$`800_Female`)   # Becomes linear once position is accounted
test_predictor_linearity(base_models$`1500_Male`)    # fresh_S_p still sign => quadratic model

# Final models
d_800_male    <- data_list$`800_Male`
d_800_female  <- data_list$`800_Female`
d_1500_male   <- data_list$`1500_Male`
d_1500_female <- data_list$`1500_Female`

S_perf_model_800_male <- glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                               family = "binomial", data = d_800_male,
                               control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

S_perf_model_800_female <- glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                                 family = "binomial", data = d_800_female,
                                 control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

S_perf_model_1500_male <- glmer(Top3 ~ fresh_S_c + I(fresh_S_c^2) + S_drop_neg + position_neg + (1 | season_id),
                                family = "binomial", data = d_1500_male,
                                control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

S_perf_model_1500_female <- glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                                  family = "binomial", data = d_1500_female,
                                  control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

top3_models <- list(
  "800_Male"    = S_perf_model_800_male,
  "800_Female"  = S_perf_model_800_female,
  "1500_Male"   = S_perf_model_1500_male,
  "1500_Female" = S_perf_model_1500_female
)

summary(S_perf_model_1500_male)

purrr::map(top3_models, ~ .x@optinfo$conv$lme4$messages)  

# compute the LRT on top3_models$`1500_Male` 
lin_1500M_full <- glmer(Top3 ~ fresh_S_c + S_drop_neg + position_neg + (1 | season_id),
                        family = "binomial", data = data_list$`1500_Male`,
                        control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

lrt_stat <- as.numeric(2 * (logLik(top3_models$`1500_Male`) - logLik(lin_1500M_full)))
lrt_df   <- attr(logLik(top3_models$`1500_Male`), "df") - attr(logLik(lin_1500M_full), "df")
lrt_p    <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

lrt_stat; lrt_df; lrt_p

# Odds ratios 
or_tables <- purrr::map(top3_models, ~
                          broom.mixed::tidy(.x, conf.int = TRUE, conf.method = "Wald", exponentiate = TRUE))

plot(ggpredict(top3_models$`1500_Male`,  terms = "fresh_S_c [all]"))

# 7.2. Diagnostics on the final models

# A) Singularity / convergence
purrr::map(top3_models, isSingular)

# B) Multicollinearity
purrr::map(top3_models, car::vif)

# C) Residual diagnostics
sim_resids <- purrr::map(top3_models, ~ simulateResiduals(.x, n = 750))
purrr::iwalk(sim_resids, ~ { plot(.x); title(main = .y, line = -1) })

DHARMa::testQuantiles(sim_resids$`800_Male`)$p.value
DHARMa::testQuantiles(sim_resids$`800_Female`)$p.value
DHARMa::testQuantiles(sim_resids$`1500_Male`)$p.value
DHARMa::testQuantiles(sim_resids$`1500_Female`)$p.value

# D) Separation -- eyeball for absurdly large estimates/SEs
purrr::map(top3_models, ~ summary(.x)$coefficients)

# E) Influential athlete-seasons
infl <- purrr::map(top3_models, ~ influence.ME::influence(.x, group = "season_id"))
purrr::map(infl, cooks.distance)

# Check without Anass ESSAYI_MAR_2025 (8x the next-highest point in that stratum)

data_list$`1500_Male` %>% filter(season_id == "Anass ESSAYI_MAR_2025")

d_1500_male_sens <- data_list$`1500_Male` %>% filter(season_id != "Anass ESSAYI_MAR_2025")

model_1500_male_sens <- glmer(Top3 ~ fresh_S_c + I(fresh_S_c^2) + S_drop_neg + position_neg + (1 | season_id),
                              family = "binomial", data = d_1500_male_sens,
                              control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

summary(model_1500_male_sens)$coefficients
summary(S_perf_model_1500_male)$coefficients

model_1500_male_linear_sens <- glmer(Top3 ~ fresh_S_c + S_drop_neg + position_neg + (1 | season_id),
                                     family = "binomial", data = d_1500_male_sens,
                                     control = glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL)))

lrt_stat <- as.numeric(2 * (logLik(model_1500_male_sens) - logLik(model_1500_male_linear_sens)))
lrt_df   <- attr(logLik(model_1500_male_sens), "df") - attr(logLik(model_1500_male_linear_sens), "df")
pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

summary(model_1500_male_linear_sens)$coefficients

# => Reinforces the quadratic version 

# F) Discrimination -- report the marginal (fixed-effects-only) AUC as the
#    honest number; conditional AUC mostly reflects athlete identity, not
#    the durability/position metrics themselves.
purrr::map2(top3_models, data_list, ~ pROC::auc(pROC::roc(.y$Top3, fitted(.x))))
purrr::map2(top3_models, data_list, ~ pROC::auc(pROC::roc(.y$Top3, predict(.x, re.form = NA, type = "response"))))

# 7.3. Plotting the results (Figure 2)
# -----------------------------------------------------------------------------
# Figure 2 will plot the OR for all model and all covariates, except Fresh_S
# in the male 1500m, bc of the quadratic shape. 
# OR change across the spectrum of Fresh S and must be draw in a table (see below)
# -----------------------------------------------------------------------------

# Extract OR/CI/p from the final glmer models ---
extract_or <- function(model, model_name) {
  est <- coef(summary(model))
  tibble(
    Model = model_name,
    term = rownames(est),
    estimate = exp(est[, "Estimate"]),
    conf.low = exp(est[, "Estimate"] - 1.96 * est[, "Std. Error"]),
    conf.high = exp(est[, "Estimate"] + 1.96 * est[, "Std. Error"]),
    p_value = est[, "Pr(>|z|)"]
  )
}

OR_df <- purrr::imap_dfr(top3_models, extract_or) %>%
  filter(term %in% c("fresh_S", "S_drop_neg", "position_neg")) %>%
  mutate(
    Predictor = case_when(
      term == "fresh_S"           ~ "Fresh S (pure speed)",
      term == "S_drop_neg"        ~ "-\u0394S (S durability)",
      term == "position_neg"      ~ "Position at 400 m to go"
    ),
    Event  = ifelse(grepl("800", Model), "800 m", "1500 m"),
    Gender = ifelse(grepl("Male", Model), "Male", "Female")
  )

# 1500_Male's fresh_S is quadratic: no single valid OR (see table 4 instead). 
# Add a placeholder row so it still appears in the layout.
OR_df <- OR_df %>%
  bind_rows(tibble(
    Model = "1500_Male", term = "fresh_S", estimate = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, p_value = NA_real_,
    Predictor = "Fresh S (pure speed)", Event = "1500 m", Gender = "Male"
  ))

OR_df$Predictor <- factor(OR_df$Predictor,
                          levels = c("Fresh S (pure speed)", "-\u0394S (S durability)", "Position at 400 m to go"))

# Text labels
OR_df_forest <- OR_df %>%
  mutate(
    group_label = paste0(Event, " - ", Gender),
    or_ci_text = ifelse(is.na(estimate), "Nonlinear (see Figure 3)",
                        sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high)),
    p_text = case_when(
      is.na(p_value)  ~ "",
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01  ~ sprintf("p = %.3f", p_value),
      TRUE            ~ sprintf("p = %.2f", p_value)
    )
  )

# Y positions: 3 predictors x 4 strata
OR_df_forest <- OR_df_forest %>%
  mutate(
    y_pos = case_when(
      group_label == "800 m - Male"    & Predictor == "Fresh S (pure speed)"    ~ 11.0,
      group_label == "800 m - Male"    & Predictor == "-\u0394S (S durability)" ~ 10.5,
      group_label == "800 m - Male"    & Predictor == "Position at 400 m to go"     ~ 10.0,
      group_label == "800 m - Female"  & Predictor == "Fresh S (pure speed)"    ~ 9.0,
      group_label == "800 m - Female"  & Predictor == "-\u0394S (S durability)" ~ 8.5,
      group_label == "800 m - Female"  & Predictor == "Position at 400 m to go"     ~ 8.0,
      group_label == "1500 m - Male"   & Predictor == "Fresh S (pure speed)"    ~ 7.0,
      group_label == "1500 m - Male"   & Predictor == "-\u0394S (S durability)" ~ 6.5,
      group_label == "1500 m - Male"   & Predictor == "Position at 400 m to go"     ~ 6.0,
      group_label == "1500 m - Female" & Predictor == "Fresh S (pure speed)"    ~ 5.0,
      group_label == "1500 m - Female" & Predictor == "-\u0394S (S durability)" ~ 4.5,
      group_label == "1500 m - Female" & Predictor == "Position at 400 m to go"     ~ 4.0
    )
  )

group_headers <- tibble(
  group_label = c("800 m - Male", "800 m - Female", "1500 m - Male", "1500 m - Female"),
  y_pos = c(11.5, 9.5, 7.5, 5.5)
)

# X range (exclude the NA placeholder row)
x_min <- min(OR_df_forest$conf.low, na.rm = TRUE)
x_max <- max(OR_df_forest$conf.high, na.rm = TRUE)
x_range <- x_max - x_min

# Plot -- same styling as your original
p <- ggplot(OR_df_forest, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.7, color = "gray20", na.rm = TRUE) +
  geom_point(size = 3.5, shape = 18, color = "black", na.rm = TRUE) +
  geom_text(aes(x = x_max + x_range * 0.05, label = or_ci_text), hjust = 0, size = 3.5, family = "Times New Roman") +
  geom_text(aes(x = x_max + x_range * 0.35, label = p_text), hjust = 0, size = 3.5, family = "Times New Roman") +
  geom_text(aes(x = x_min - x_range * 0.05, label = as.character(Predictor)), hjust = 1, size = 3.5, family = "Times New Roman") +
  geom_text(data = group_headers, aes(x = x_min - x_range * 0.05, y = y_pos, label = group_label),
            hjust = 1, size = 4, fontface = "bold", family = "Times New Roman") +
  annotate("text", x = x_max + x_range * 0.35, y = 12.2, label = "p-value", hjust = 0, size = 3.5, fontface = "bold", family = "Times New Roman") +
  annotate("text", x = x_max + x_range * 0.05, y = 12.2, label = "OR (95% CI)", hjust = 0, size = 3.5, fontface = "bold", family = "Times New Roman") +
  scale_x_continuous(name = "Odds Ratio", limits = c(x_min - x_range * 0.30, x_max + x_range * 0.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(3.3, 12.9), expand = c(0, 0)) +
  labs(title = "Determinants of Top-3 finish",
       subtitle = "Odds ratios, 95% confidence intervals and p-values from logistic regression") +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = unit(c(5, 5, 5, 5), "mm")
  )

print(p)
ggsave("Figure_2_forest_plot_OR_top3.png", p, width = 15, height = 9, dpi = 600, bg = "white")

# 7.4. Computing OR from quadratic model (Table and Figure)

# Function to compute the conditional OR
compute_conditional_OR <- function(model, x_values, delta = 1,
                                   linear_term = "fresh_S_c", quad_term = "I(fresh_S_c^2)") {
  b <- fixef(model)
  V <- vcov(model)
  b1 <- b[linear_term]; b2 <- b[quad_term]
  Vsub <- V[c(linear_term, quad_term), c(linear_term, quad_term)]
  
  purrr::map_dfr(x_values, function(x) {
    log_or <- b1 * delta + b2 * ((x + delta)^2 - x^2)
    grad <- c(delta, (x + delta)^2 - x^2)
    se_log_or <- sqrt(as.numeric(t(grad) %*% Vsub %*% grad))
    z <- log_or / se_log_or
    tibble::tibble(
      x = x, OR = exp(log_or),
      conf.low = exp(log_or - 1.96 * se_log_or),
      conf.high = exp(log_or + 1.96 * se_log_or),
      p_value = 2 * pnorm(-abs(z))
    )
  })
}

# Data description
mean_freshS_1500M <- mean(data_list$`1500_Male`$fresh_S)
summary(data_list$`1500_Male`$fresh_S)
raw_x <- 8:17
centered_x <- raw_x - mean_freshS_1500M

# OR computation
or_results <- compute_conditional_OR(top3_models$`1500_Male`, centered_x)
or_results$raw_x <- raw_x
or_results$raw_x_to <- raw_x + 1
print(or_results)

# Table
table4 <- compute_conditional_OR(top3_models$`1500_Male`, 8:17 - mean(data_list$`1500_Male`$fresh_S)) %>%
  mutate(
    `Fresh S range` = paste0(8:17, " \u2192 ", 9:18),
    `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", OR, conf.low, conf.high),
    `p-value` = case_when(p_value < 0.001 ~ "p < 0.001", p_value < 0.01 ~ sprintf("p = %.3f", p_value), TRUE ~ sprintf("p = %.2f", p_value))
  ) %>%
  select(`Fresh S range`, `OR (95% CI)`, `p-value`)

print(table4)

# Plotting conditional OR curve

mean_freshS_1500M <- mean(
  data_list$`1500_Male`$fresh_S,
  na.rm = TRUE
)

raw_grid <- seq(from = 8, to = 17, length.out = 300)
centered_grid <- raw_grid - mean_freshS_1500M

or_curve <- compute_conditional_OR(top3_models$`1500_Male`, centered_grid) %>%
  mutate(raw_x = raw_grid)

# find where the CI stops excluding 1 (precise, not eyeballed)
first_ns_x <- or_curve %>% filter(conf.low <= 1 & conf.high >= 1) %>% slice(1) %>% pull(raw_x)

p2 <- ggplot(or_curve, aes(x = raw_x, y = OR)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.20, fill = "gray50") +
  geom_line(linewidth = 0.9, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, color = "gray40") +
  geom_vline(xintercept = first_ns_x, linetype = "dotted", linewidth = 0.5, color = "gray40") +
  annotate("text", x = first_ns_x, y = 6, label = paste0("CI includes 1\nfrom ", round(first_ns_x, 1)),
           hjust = -0.1, size = 3, family = "Times New Roman", color = "gray30") +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  scale_x_continuous(breaks = 8:17) +
  labs(
    title = "Conditional odds ratio for fresh S, 1500 m Male",
    x = expression(paste("Fresh ", italic(S), " (m" %.% s^{-1}, ")")),
    y = expression("OR for +1 m" %.% s^{-1} * " in Fresh " * italic(S))
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))

ggsave("Figure_3_conditional_OR_1500MAN_top3.png", p2, width = 15, height = 9, dpi = 600, bg = "white")

# =============================================================================
# 8. ECSS YIA talk - slide figures
# =============================================================================
# -----------------------------------------------------------------------------
# In this speech, we choose to work with a specific race as a case study.
# Race: Tokyo 2025 World Championships, women's 800 m final
# -----------------------------------------------------------------------------

# 8.0. Theme, subset and athletes

# Reusable slide theme (legible from the back of a hall)
#   Print figures use Times + size 14; projection needs sans + larger.
#   base_family = "" keeps it portable (no font registration needed).
theme_slide <- function(base_size = 20, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size * 1.20),
      plot.subtitle   = element_text(size = base_size * 0.75, colour = "grey30"),
      axis.title      = element_text(face = "bold"),
      axis.text       = element_text(colour = "grey20"),
      axis.line       = element_line(linewidth = 0.6),
      legend.position = "right",
      legend.title    = element_blank(),
      legend.text     = element_text(size = base_size * 0.70)
    )
}

# Subset the final
final_id <- "World Championship Tokyo 2025_Final_45921_Female"

final_800 <- long_800_fatigue %>%
  filter(race_id == final_id, !is.na(split_speed_m_s)) %>%  # drops the index-0 init row
  arrange(rank, split_index)

# Who is actually in the modelled data? (check for Hodgkinson, the two Moraas)
final_800 %>% distinct(Athletes, rank) %>% arrange(rank) %>% print(n = Inf)

# 8.1. The hook: "what the field did"

# Fig 1: Pacing profiles, coloured by athlete
# Exploration version: every athlete a distinct colour so you can SEE the
# Shapes and then choose the cold-open cast.
fig1 <- ggplot(final_800,
               aes(x = cum_distance_m, y = split_speed_m_s,
                   colour = reorder(Athletes, rank), group = Athletes)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(100, 800, 100)) +
  labs(x = "Distance (m)",
       y = expression(Section~speed~(m%.%s^-1)),
       title = "What the field did?",
       subtitle = "Women's 800 m World Championships 2025 final - 100 m section speeds") +
  theme_slide()

print(fig1)
#ggsave("fig1A_ECSS_pacing.png", fig1, width = 12, height = 7, dpi = 300)

# Figure 1b: highlight version (use AFTER you pick the cast)
# Edit `focus` to the exact athlete strings printed above.
focus <- c("Mary MORAA", "Lilian ODIRA") 

plot_dat <- final_800 %>% mutate(hl = ifelse(Athletes %in% focus, Athletes, "Other"))

fig1b <- ggplot() +
  geom_line(data = filter(plot_dat, hl == "Other"),
            aes(cum_distance_m, split_speed_m_s, group = Athletes),
            colour = "grey82", linewidth = 0.8) +
  geom_line(data = filter(plot_dat, hl != "Other"),
            aes(cum_distance_m, split_speed_m_s, colour = hl, group = Athletes),
            linewidth = 1.5) +
  geom_point(data = filter(plot_dat, hl != "Other"),
             aes(cum_distance_m, split_speed_m_s, colour = hl), size = 2.6) +
  scale_x_continuous(breaks = seq(100, 800, 100)) +
  scale_colour_manual(values = c("Mary MORAA" = "#DF2935", 
                                 "Lilian ODIRA" = "#2A8459")) +
  labs(x = "Distance (m)", y = expression(Section~speed~(m%.%s^-1)),
       title = "What the field did?",
       subtitle = "Women's 800 m World Championships 2025 final - 100 m section speeds") +
  annotate("text", x = 215, y = 7.62, label = "led at 400m and 600m", colour = "#DF2935",
           fontface = "italic", size = 3.8, hjust = 0) +
  annotate("text", x = 760, y = 6.3, label = "collapsed", colour = "#DF2935",
           fontface = "italic", size = 3.8, hjust = 0) +
  annotate("text", x = 740, y = 7, label = "won the race", colour = "#2A8459",
           fontface = "italic", size = 3.8, hjust = 0) +
  theme_slide()

print(fig1b)
#ggsave("fig1B_ECSS_pacing_highlight.png", fig1b, width = 12, height = 7, dpi = 300)

# 8.2. The model

# Fig. 2: Showing the downshift of S in image rather than in math

E <- 0.90
t <- seq(5, 300, 1)
pl <- function(S, E, t) S * t^(E - 1)

mod <- dplyr::bind_rows(
  tibble::tibble(t, v = pl(12.6, E, t), state = "Fresh"),
  tibble::tibble(t, v = pl(9.7,  E, t), state = "Fatigued"))

mod1 <- dplyr::bind_rows(
  tibble::tibble(t, v = pl(12.6, E, t), state = "Fresh"))

fig_model <- ggplot(mod, aes(t, v, colour = state)) +
  geom_line(linewidth = 1.3) +
  scale_colour_manual(values = c(Fresh = "#2A8459", Fatigued = "#D7191C"), name = NULL,
                      breaks = c("Fresh", "Fatigued")) +
  annotate("segment", x = 110, xend = 110, y = pl(12.6,E,110), yend = pl(9.7,E,110),
           arrow = arrow(length = unit(.5,"cm")), colour = "grey45") +
  annotate("text", x = 116, y = mean(c(pl(9.7,E,110), pl(12.6,E,110))),
           label = "S downshifts with fatigue\n(E unchanged)", hjust = 0, colour = "grey40", size = 4) +
  labs(x = "Race duration (s)", y = expression(Speed~(m%.%s^-1)),
       title = "Power-Law speed-duration model: fresh vs fatigued") +
  theme_slide()

fig_fresh <- ggplot(mod1, aes(t, v, colour = state)) +
  geom_line(linewidth = 1.3) +
  ylim(c(5,11)) +
  scale_colour_manual(values = c(Fresh = "#2A8459"), name = NULL,
                      breaks = c("Fresh")) +
  labs(x = "Race duration (s)", y = expression(Speed~(m%.%s^-1)),
       title = "Power-Law speed-duration model: fresh vs fatigued") +
  theme_slide()
# add scale_x_log10() + scale_y_log10() if you want the two curves exactly parallel (proves E constant)

# Fig. 2B: showing S decrease and risk ration increase in image rather than in maths
# -----------------------------------------------------------------------------
# Two-panel feature figure for ONE athlete (kept generic / unnamed).
#   Panel A: S_dyn_drake (speed ceiling) falls
#   Panel B: risk_ratio = split_speed_m_s / S_dyn_drake rises toward 1.00
# Both coloured by risk_ratio, using the deck palette (risk_cols).
#
# Uses objects you already have: the per-split data frame produced upstream
# in Dynamic_PL_global.R (one row per split_index per entry), with columns
#   split_index, split_speed_m_s, S_dyn_drake, risk_ratio
# Distance is split_index * 100 for the 800 m.
#
# No athlete name is drawn, by design (keeps slide 8 as the reveal).
# -----------------------------------------------------------------------------

risk_cols <- c("#2A8459", "#F3CA40", "#D7191C")   # low -> high risk (deck palette)

plot_two_panel_one_athlete <- function(df,
                                       v_col     = "split_speed_m_s",
                                       S_col     = "S_dyn_drake",
                                       risk_col  = "risk_ratio",
                                       idx_col   = "split_index",
                                       section_m = 100,
                                       n_grid    = 400,
                                       smooth    = TRUE) {
  
  d <- df %>%
    filter(!is.na(.data[[risk_col]]), .data[[idx_col]] > 0) %>%
    arrange(.data[[idx_col]]) %>%
    transmute(dist = .data[[idx_col]] * section_m,
              v    = .data[[v_col]],
              S    = .data[[S_col]],
              risk = .data[[risk_col]])
  
  stopifnot(nrow(d) >= 3)
  
  # --- -dS from RAW endpoints (matches the table, not the smoothed grid) ----
  S_first <- d$S[1]
  S_last  <- d$S[nrow(d)]
  dS_pct  <- (S_first - S_last) / S_first * 100      # loss as POSITIVE %
  x_end   <- max(d$dist)
  x_brk   <- x_end + section_m * 0.55                # bracket x, just right of curve
  
  gx <- seq(min(d$dist), max(d$dist), length.out = n_grid)
  f  <- if (smooth && nrow(d) >= 4) splinefun else function(x, y) approxfun(x, y, rule = 2)
  vg <- f(d$dist, d$v)(gx)
  Sg <- f(d$dist, d$S)(gx)
  rg <- vg / Sg
  g  <- data.frame(dist = gx, v = vg, S = Sg, risk = rg)
  
  rng <- range(c(d$risk, rg), na.rm = TRUE)
  
  base <- function() {
    list(
      scale_colour_gradientn(colours = risk_cols, limits = c(min(0.3, rng[1]), 1),
                             oob = scales::squish, guide = "none"),
      scale_x_continuous(breaks = seq(200, max(gx), by = 200)),
      labs(x = "Distance (m)"),
      theme_minimal(base_size = 15),
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0),
            plot.subtitle = element_text(colour = "grey30"))
    )
  }
  
  # bracket geometry (a squared-off right bracket spanning S_last..S_first)
  tick <- section_m * 0.15
  brk <- data.frame(
    x    = c(x_brk - tick, x_brk, x_brk, x_brk - tick),
    y    = c(S_first, S_first, S_last, S_last)
  )
  
  pA <- ggplot(g, aes(dist, S, colour = risk)) +
    geom_line(linewidth = 2) +
    geom_path(data = brk, aes(x, y), inherit.aes = FALSE,
              colour = "grey30", linewidth = 0.6) +
    annotate("text", x = x_brk + section_m * 0.25, y = (S_first + S_last) / 2,
             hjust = 0, colour = "grey20", size = 4.4, fontface = "bold",
             label = sprintf("-\u0394S = %.0f%%", dS_pct)) +
    base() +
    coord_cartesian(xlim = c(min(gx), x_brk + section_m * 1.4), clip = "off") +
    labs(y = expression("S  (m" %.% s^-1 * ")"),
         title = "Speed ceiling S falls",
         subtitle = expression("Durability of S  (\u2212\u0394S)"))
  
  pB <- ggplot(g, aes(dist, risk, colour = risk)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "#C0392B", linewidth = 0.7) +
    annotate("text", x = min(gx), y = 1.02, hjust = 0, colour = "#C0392B",
             size = 4.3, label = "risk = 1.00 (ceiling)") +
    geom_line(linewidth = 2) +
    base() +
    coord_cartesian(ylim = c(min(0.3, rng[1]), 1.08)) +
    labs(y = "risk = v / S", title = "Risk ratio rises")
  
  pA | pB
}

# ---- usage -----------------------------------------------------------------
# Replace `per_split_df` with your upstream per-split frame (e.g. all_fatigue
# before the na-filter, or the same object feeding Fig 6A / slide 8), and
# `moraa_key` with her entry_id / Athletes value for the Tokyo 2025 final.

werro <- final_800 %>%
  filter(Athletes == "Audrey WERRO")

fig <- plot_two_panel_one_athlete(werro)
#ggsave("fig2B_ECSS_features_def.png", fig, width = 10.6, height = 4.5, dpi = 300)

# 8.3. Chapter 1: "what it cost them" - risk drives the collapse

# Fig. 3B: The mean split by risk so the driver shows

fa <- all_fatigue %>%                       # all races, long format
  group_by(Athletes, race_id) %>% mutate(race_risk = mean(risk_ratio, na.rm = TRUE)) %>% ungroup() %>%
  mutate(risk_grp = ifelse(race_risk >= median(race_risk, na.rm = TRUE), "High risk", "Low risk"))

fa <- fa %>%
  filter(S <= 22)

means <- fa %>%
  group_by(event, Gender, risk_grp, split_index) %>%
  summarise(S = mean(S_dyn_drake, na.rm = TRUE),
            x = mean(race_prop,   na.rm = TRUE), .groups = "drop")

ev_lab <- ggplot2::as_labeller(c("800" = "800 m", "1500" = "1500 m",
                                 "Male" = "Male", "Female" = "Female"))

fig3 <- ggplot() +
  geom_line(data = fa, aes(race_prop, S_dyn_drake, group = interaction(Athletes, race_id)),
            colour = "grey88", linewidth = .25, alpha = .4) +                 # field recedes
  geom_line(data = means, aes(x, S, colour = risk_grp, group = risk_grp),
            linewidth = 1.7) +                                                # means dominate
  scale_colour_manual(values = c("Low risk" = "#2A8459", "High risk" = "#D7191C"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +        # drop this line if race_prop is already 0-100
  facet_grid(Gender ~ event, labeller = ev_lab) +
  labs(x = "Race progression", y = expression(Speed~parameter~S~(m%.%s^-1)),
       title = "Risk drives the real-time decline of S",
       subtitle = "Grey: individual races. Coloured: mean trajectory by race-level risk\nHigh risk ??? median split on race-level risk ratio; Low risk < median.") +
  theme_slide() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "grey20"),
        legend.position = "top")

#ggsave("fig3_ECSS_S_traj_risk.png", fig3, width = 12, height = 7, dpi = 300)

# 8.4. Fresh S x depletion, the speed/endurance trade-off
# Colored by E = speed-endurance trade-off

# Fig. 4A: population trade-off (n = 896)

# one row per 800 athlete-race: fresh_S, dS (= ?????S, %), E
sum_fat800 <- long_800_fatigue %>%
  filter(!is.na(split_speed_m_s)) %>%
  arrange(race_id, Athletes, split_index) %>%
  group_by(race_id, Athletes, Gender, rank, E) %>%               # race_id => one row per athlete-race
  summarise(fresh_S = first(S_dyn_drake),
            final_S = last(S_dyn_drake),
            dS      = 100 * (fresh_S - final_S) / fresh_S,
            .groups = "drop") %>%
  group_by(last_name = str_remove(Athletes, "^\\S+\\s+")) %>%
  mutate(short = if (n_distinct(Athletes) > 1)
    paste0(str_sub(Athletes, 1, 1), ". ", last_name) else last_name) %>%
  ungroup() %>% select(-last_name)

# Plots
fig4A_all <- ggplot(sum_fat800, aes(fresh_S, dS)) +
  geom_point(aes(colour = E), size = 2.2, alpha = .5) +
  geomtextpath::geom_textsmooth(aes(label = "population trend"),
                                method = "lm", formula = y ~ x, se = FALSE,
                                colour = "grey40", linetype = "dashed", linewidth = .9,
                                fontface = "italic", size = 4, hjust = .85, vjust = -.2) +
  scale_colour_gradientn(colours = c("#D7191C","#F3CA40","#2A8459"), name = "Endurance E") +
  labs(x = expression(Fresh~S~(m%.%s^-1)~"[maximal theoretical speed]"), 
       y = expression(-Delta*S~"(%) [speed durability]"),
       title = "Speed buys risk; endurance decides if it pays",
       subtitle = "All 800 m entries (n = 896); colour = endurance parameter E") +
  theme_slide()

fig_4A_gender <- ggplot(sum_fat800, aes(fresh_S, dS)) +
  geom_point(aes(colour = E), size = 2.2, alpha = .5) +
  geomtextpath::geom_textsmooth(aes(label = Gender, group = Gender),
                                method = "lm", formula = y ~ x, se = FALSE,
                                colour = "grey25", linewidth = 1,
                                fontface = "bold", size = 4.4, hjust = .9, vjust = -.3) +
  scale_colour_gradientn(colours = c("#D7191C","#F3CA40","#2A8459"), name = "Endurance E") +
  labs(x = expression(Fresh~S~(m%.%s^-1)~"[maximal theoretical speed]"), 
       y = expression(-Delta*S~"(%) [speed durability]"),
       title = "Speed buys risk; endurance decides if it pays",
       subtitle = "All 800 m entries (n = 896); colour = endurance parameter E") +
  theme_slide()

#ggsave("fig4_ECSS_tradeoff.png", fig_4A_gender, width = 11, height = 7, dpi = 300)

# Fig. 4B: race specific trade-off

# one row per finalist
summ <- summ %>%
  mutate(last_name  = str_remove(Athletes, "^\\S+\\s+"),   # drop first token only (non-greedy)
         first_init = str_sub(Athletes, 1, 1)) %>%
  group_by(last_name) %>%
  mutate(short = if (n_distinct(Athletes) > 1)
    paste0(first_init, ". ", last_name)      # "M. MORAA", "S. MORAA"
    else last_name) %>%                        # "ODIRA", "HUNTER BELL", ...
  ungroup() %>%
  select(-last_name, -first_init)

# population relationship, female 800 m, all entries (not just this final)
pop <- long_800_fatigue %>%
  filter(Gender == "Female", !is.na(split_speed_m_s)) %>%
  group_by(entry_id) %>%
  arrange(split_index, .by_group = TRUE) %>%
  summarise(fresh_S = first(S_dyn_drake),
            dS      = 100 * (first(S_dyn_drake) - last(S_dyn_drake)) / first(S_dyn_drake),
            .groups = "drop")

pop_fit <- lm(dS ~ fresh_S, data = pop)     # you already have S_drop_pct ~line 834 if you prefer
coef(pop_fit)                                # intercept, slope

# Plot fig. 4B
fig4B <- ggplot(summ, aes(fresh_S, dS)) +
  # population trend goes here instead of a 7-point lm. Replace the slope/intercept
  # with the values implied by the random-effects (baseline S vs rate-of-decline,
  # r = -0.92) so the line is "Table 2, with these athletes placed on it".
  geom_textabline(intercept = coef(pop_fit)[1], slope = coef(pop_fit)[2],
                  label = "population trend",
                  linetype = "dashed", colour = "grey50",
                  size = 4.2, fontface = "italic",
                  hjust = 0.85, vjust = -0.3) +  # hjust slides along the line; vjust lifts it off 
  geom_point(aes(colour = E), size = 7, alpha = 0.95) +
  geom_text(aes(label = short), nudge_y = 2.2, size = 4, fontface = "bold",
            colour = "grey20") +
  scale_colour_gradientn(colours = c("#D7191C","#F3CA40","#2A8459"),
                         name = "Endurance E") +
  labs(x = expression(Fresh~S~(m%.%s^-1)~"[maximal theoretical speed]"),
       y = expression(-Delta*S~"(%) [speed durability]"),
       title = "Speed buys risk, endurance decides if it pays",
       subtitle = "Each point a finalist; dashed line = all female entries (n = 564); colour = endurance parameter E") +
  theme_slide()

# 8.5. Forest plot: relation to performance

# Fig. 5: Speed vs. durability for a Top-3 finish

# Create dataset
OR_df_forest_2 <- OR_df_forest %>% 
  mutate(sig   = p_value < .05)

x_max <- max(OR_df_forest_2$conf.high, na.rm = TRUE)
x_min <- min(OR_df_forest_2$conf.low, na.rm = TRUE)
x_range <- x_max - x_min

OR_df_forest_2 <- OR_df_forest_2 %>% 
  mutate(
    Predictor = ifelse(Predictor == "Fresh S (pure speed)", "Fresh S [maximal theoretical speed]", 
                       ifelse(Predictor == "Position at 400 m to go", "Position at 400 m to go", "-\u0394S [speed durability]")))

# Center of each event-sex group, computed from the y_pos actually in use
y_breaks <- OR_df_forest_2 %>%
  distinct(group_label, y_pos) %>%
  group_by(group_label) %>%
  summarise(y_mid = mean(range(y_pos)), .groups = "drop") %>%
  arrange(desc(y_mid))

# Header sits just above the topmost row, whatever that row's y_pos is
header_y <- max(OR_df_forest_2$y_pos, na.rm = TRUE) + 0.7

forest_fig <- ggplot(OR_df_forest_2, aes(x = estimate, y = y_pos, colour = Predictor)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, linewidth = 1) +
  geom_point(size = 4.4) +
  geom_text(aes(x = x_max + x_range * 0.35, label = p_text), hjust = 0, size = 5, colour = "black") +
  annotate("text", x = x_max + x_range * 0.35, y = header_y, label = "p-value",
           hjust = 0, size = 5, fontface = "bold") +
  scale_colour_manual(values = c("Fresh S [maximal theoretical speed]" = "#2A8459",
                                 "-\u0394S [speed durability]" = "#D7191C",
                                 "Position at 400 m to go" = "#F3CA40"), name = NULL) +
  scale_x_continuous(trans = "log", breaks = c(.8, 1, 1.25, 1.6)) +
  coord_cartesian(xlim = c(0.75, x_max + x_range * 0.75), clip = "off") +
  scale_y_continuous(breaks = y_breaks$y_mid, labels = y_breaks$group_label) +
  labs(x = "Odds ratio for a Top-3 finish",
       y = NULL,
       title = "Raw speed and durability predict success; their balance shifts by event",
       subtitle = "OR (95% CI); dashed line: OR = 1") +
  theme_slide()

print(forest_fig)
#ggsave("fig5_ECSS_forest_slide.png", forest_fig, width = 15, height = 6, dpi = 300)

# 8.6. The closing: what does it mean for our race?

# Fig. 6A: Speed parameter trajectory coloured by risk

focus     <- c("Mary MORAA", "Lilian ODIRA")
risk_cols <- c("#2A8459", "#F3CA40", "#D7191C")   # low -> high risk (blue-yellow-red)

dat6 <- final_800 %>%
  filter(!is.na(risk_ratio), split_index > 0) %>%
  group_by(Athletes) %>%
  mutate(S_fresh = first(S_dyn_drake),
         S_pct   = 100 * S_dyn_drake / S_fresh) %>%
  ungroup() %>%
  mutate(hl    = ifelse(Athletes %in% focus, Athletes, "Other"),
         short = paste0(str_sub(Athletes, 1, 1), ". ",
                        str_remove(Athletes, "^\\S+\\s+")))   # "G. HUNTER BELL", "M. MORAA"

build_fig6 <- function(yvar, ylab) {
  labd <- dat6 %>% filter(hl != "Other") %>%
    group_by(Athletes) %>% filter(split_index == max(split_index)) %>% ungroup()
  
  ggplot() +
    # background field
    geom_line(data = filter(dat6, hl == "Other"),
              aes(cum_distance_m, .data[[yvar]], group = Athletes),
              colour = "grey82", linewidth = 0.8) +
    # focal athletes: line coloured by risk gradient (turns "hot" before the plunge)
    geom_path(data = filter(dat6, hl != "Other"),
              aes(cum_distance_m, .data[[yvar]], colour = risk_ratio, group = Athletes),
              linewidth = 1.9, lineend = "round") +
    geom_text(data = labd,
              aes(cum_distance_m, .data[[yvar]], label = short),
              hjust = -0.12, size = 4.2, fontface = "bold", colour = "grey20") +
    geom_vline(xintercept = 400, linetype = "dashed", colour = "grey55") +
    annotate("text", x = 400, y = Inf, label = "400 m to go",
             vjust = 1.4, hjust = 1.08, size = 4, colour = "grey45") +
    scale_colour_gradientn(colours = risk_cols, limits = c(0.3, 1.0),
                           name = "Risk ratio\n(v / S)") +
    scale_x_continuous(breaks = seq(100, 800, 100),
                       expand = expansion(mult = c(0.02, 0.26))) +  # room for labels
    labs(x = "Distance (m)", y = ylab) +
    theme_slide()
}

# NORMALISED
# survives (Moraa -> ~38 %), and Odira's inflated fresh S (~21.6) does not blow
# up the axis. Fresh S itself is carried by Figs 3-4, so Fig 2 need not show it.
fig6_norm <- build_fig6("S_pct", "Speed ceiling (% of fresh S)") +
  labs(title = "What it cost them",
       subtitle = "Speed-ceiling decay, coloured by proximity to the ceiling")

# ABSOLUTE
# the rest of the field compresses.
fig6_abs <- build_fig6("S_dyn_drake", expression(Speed~ceiling~S~(m%.%s^-1))) +
  labs(title = "What it cost them",
       subtitle = "Speed parameter decay, coloured by proximity to the ceiling (risk-ratio)")

print(fig6_abs)

ggsave("fig6_ECSS_S_risk.png", fig6_abs, width = 11, height = 7, dpi = 300)


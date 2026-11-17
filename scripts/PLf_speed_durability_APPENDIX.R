#####################################################################################################
#                                                                                                   #
#                 APPENDIX: Fatigue-integrated Power Law (PLf) and Speed Durability                 #
#                 -----------------------------------------------------------------                 #
#                                                                                                   #
#   0. Load the libraries and datasets                                                              #
#   1. Setting the tactical multiplier (phi) value                                                  #
#   2. Sensitivity of Delta S to the tactical multiplier (adjustment 5% vs 0%)                      #
#   3. Feasibility floor hitting                                                                    #
#   4. Risk ratio and S: fundamental mathematical issue                                             #
#   5. Relation to performance: sensitivity to the definition of S durability                       #
#   6. Sensitivity of the Top-3 models to the feasibility floor                                     #
#                                                                                                   #
#                                                                        Maxime Walt, 2025-2026     #
#                                                                                                   #
#####################################################################################################

# -----------------------------------------------------------------------------
# READ ME:
#  This script rely entirely on the "PLf_speed_durability.R" script.
#  It has been created to answer "appendix" questions and avoid the analysis 
#  process to be overcomplex.
#  -----------------------------------------------------------------------------

# =============================================================================
# 0. Load the libraries and datasets
# =============================================================================

# 0.1. Load the libraries

library(dplyr); library(purrr); library(tibble); library(future.apply); library(emmeans)

# =============================================================================
# 1. Setting the tactical multiplier (phi) value
# =============================================================================
# -----------------------------------------------------------------------------
# Requires: long_800, long_1500, compute_dynamic_fatigue_drake(), 
# predict_remaining_time()
# all sourced from PLf_speed_durability.R (sections 1 and 4).
# -----------------------------------------------------------------------------

set.seed(2026)  # reproducibility for fold assignment

tactical_grid <- seq(0.00, 0.10, by = 0.01)  # check all possibilities bet 0 and 10%

# 1.1. Distribution of |position_change| >= 2

mean(abs(long_800$position_change) >= 2, na.rm = TRUE)
mean(abs(long_1500$position_change) >= 2, na.rm = TRUE)

# 1.2. MARE for a single (dataset, tactical_effect) pair
compute_mare <- function(df, tactical_effect, total_distance) {
  tmp <- compute_dynamic_fatigue_drake(df, tactical_effect = tactical_effect) %>%
    group_by(entry_id, race_id) %>%
    arrange(split_index, .by_group = TRUE) %>%
    mutate(remaining_distance       = total_distance - split_index * 100,
           predicted_remaining_time = predict_remaining_time(S_post, b, remaining_distance),
           predicted_total_time     = cumsum(split_time_s) + predicted_remaining_time) %>%
    ungroup()
  mean(compute_mare_full(tmp, true_col = "Final_time", pred_cols = "predicted_total_time")$MARE)
}

# Consistency check with the main script: must reproduce exactly
compute_mare(long_800,  0.05, 800)    # 0.01505652
compute_mare(long_1500, 0.05, 1500)   # 0.01382418

# 1.3. Full-sample (in-sample) curve 
insample_curve <- function(df, total_distance, grid = tactical_grid) {
  map_dfr(grid, function(te) {
    tibble(tactical_effect = te, MARE = compute_mare(df, te, total_distance))
  })
}

insample_800  <- insample_curve(long_800,  800)
insample_1500 <- insample_curve(long_1500, 1500) 

baseline_800  <- compute_mare(long_800,  0.00, 800)
baseline_1500 <- compute_mare(long_1500, 0.00, 1500)

insample_best_800  <- insample_800$tactical_effect[which.min(insample_800$MARE)]
insample_best_1500 <- insample_1500$tactical_effect[which.min(insample_1500$MARE)]

# 1.4. Grouped k-fold CV (grouped by athlete_key)

run_cv <- function(df, total_distance, athlete_id_col, k = 5, grid = tactical_grid) {
  set.seed(2026)
  athletes <- unique(df[[athlete_id_col]])
  folds <- sample(rep(1:k, length.out = length(athletes)))
  fold_map <- setNames(folds, athletes)
  df$.fold <- fold_map[as.character(df[[athlete_id_col]])]
  
  results <- map_dfr(1:k, function(test_fold) {
    train_df <- df %>% filter(.fold != test_fold)
    test_df  <- df %>% filter(.fold == test_fold)
    
    train_mare <- map_dfr(grid, function(te) {
      tibble(tactical_effect = te, MARE = compute_mare(train_df, te, total_distance))
    })
    best_te <- train_mare$tactical_effect[which.min(train_mare$MARE)]
    
    test_mare <- compute_mare(test_df, best_te, total_distance)
    
    tibble(fold = test_fold, selected_tactical_effect = best_te, held_out_MARE = test_mare)
  })
  
  results
}

cv_800  <- run_cv(long_800,  800,  athlete_id_col = "athlete_key")
cv_1500 <- run_cv(long_1500, 1500, athlete_id_col = "athlete_key")

# 1.5. Manual 5% k-fold CV (grouped by athlete_key)

manual_cv_fixed <- function(df, total_distance, athlete_id_col, k = 5, fixed_te = 0.05, seed = 2026) {
  set.seed(seed)
  athletes <- unique(df[[athlete_id_col]])
  folds <- sample(rep(1:k, length.out = length(athletes)))
  fold_map <- setNames(folds, athletes)
  df$.fold <- fold_map[as.character(df[[athlete_id_col]])]
  
  map_dfr(1:k, function(test_fold) {
    test_df <- df %>% filter(.fold == test_fold)
    tibble(fold = test_fold, MARE_at_5pct = compute_mare(test_df, fixed_te, total_distance))
  })
}

fixed_cv_800  <- manual_cv_fixed(long_800,  800,  athlete_id_col = "athlete_key")
fixed_cv_1500 <- manual_cv_fixed(long_1500, 1500, athlete_id_col = "athlete_key")

# 1.6. Compare
print(insample_800);  print(cv_800)
print(insample_1500); print(cv_1500)

print(fixed_cv_800);  mean(fixed_cv_800$MARE_at_5pct);  sd(fixed_cv_800$MARE_at_5pct)
print(fixed_cv_1500); mean(fixed_cv_1500$MARE_at_5pct); sd(fixed_cv_1500$MARE_at_5pct)

# =============================================================================
# 2. Sensitivity of Delta S to the tactical multiplier (adjustment 5% vs 0%, 
#    i.e. phi = 1 in all sections)
# =============================================================================

dS_by_phi <- function(df, te) {
  compute_dynamic_fatigue_drake(df, tactical_effect = te) %>%
    group_by(race_id, entry_id) %>%
    arrange(split_index, .by_group = TRUE) %>%
    summarise(dS_raw = first(S) - last(S_dyn_drake),                          
              dS_pct = 100 * (first(S) - last(S_dyn_drake)) / first(S),       
              .groups = "drop")
}

compare_phi <- function(df, label) {
  j <- inner_join(dS_by_phi(df, 0.05), dS_by_phi(df, 0),
                  by = c("race_id", "entry_id"), suffix = c("_phi5", "_phi0"))
  tibble(event        = label,
         n_entries    = nrow(j),
         r_raw        = cor(j$dS_raw_phi5, j$dS_raw_phi0),
         r_pct        = cor(j$dS_pct_phi5, j$dS_pct_phi0),
         mean_abs_diff_pct = mean(abs(j$dS_pct_phi5 - j$dS_pct_phi0)),   # percentage points
         max_abs_diff_pct  = max(abs(j$dS_pct_phi5 - j$dS_pct_phi0)),
         mean_abs_diff_ms = mean(abs(j$dS_raw_phi5 - j$dS_raw_phi0)))  # m/s
}

bind_rows(compare_phi(long_800, "800"), compare_phi(long_1500, "1500"))

# =============================================================================
# 3. Feasibility floor hitting
# =============================================================================
# -----------------------------------------------------------------------------
# Does the tactical multiplier affect the proportion of feasibility floor hit?
# -----------------------------------------------------------------------------

floor_rate_curve <- function(df, grid = tactical_grid) {
  map_dfr(grid, function(te) {
    out <- compute_dynamic_fatigue_drake(df, tactical_effect = te)
    tibble(tactical_effect = te,
           floor_rate = mean(out$floor_bound[out$split_index >= 2]))  # transitions 2..n
  })
}

floor_curve_800  <- floor_rate_curve(long_800)
floor_curve_1500 <- floor_rate_curve(long_1500)

print(floor_curve_800)
print(floor_curve_1500)

# =============================================================================
# 4. Risk ratio and S: fundamental mathematical issue
# =============================================================================
# -----------------------------------------------------------------------------
# Risk ratio is defined as RRi = vi/Si, and the recursive PLf equation updates 
# S using vi and Si. Because vi = RRi*Si, the decline in S is mathematically 
# linked to RRi by construction.
# This part of the script aim to analyse to which extend do this dependency 
# impact the story of the paper. Concretly: will we remove risk ratio from the
# S model?
# A two-step procedure will be computed
#  1. Compute the S model (Table 3, model_S_dyn_drake) without the RR
#  2. Create a null-effect simulation 
# -----------------------------------------------------------------------------

# 4.1. Computing Table 3 (S model) without risk ratio

# Fit the model
model_S_noRR <- lmer(S_dyn_drake ~ race_prop * Competition_type3 +
                       CV + Speed_200 + Gender + event +
                       (1 + race_prop | season_id),
                     data = all_fatigue, control = lmerControl(optimizer = "bobyqa"))

# Adjust p-val using BH
table_3_noRR_adj_pval <- bh_table(model_S_noRR)

table_3_noRR_out <- table_3_noRR_adj_pval %>%
  dplyr::transmute(term,
                   Estimate = round(estimate, 3),
                   SE       = round(std.error, 3),
                   p        = fmt_p(p.value),
                   p_BH     = fmt_p(p_BH))
print(table_3_noRR_out, n = Inf)

# Context slopes of S over race progression (compare the PATTERN with Table 3:
# {MC final, DL meeting} steeper than {MC heat, MC semi-final})
emtrends(model_S_noRR, pairwise ~ Competition_type3, var = "race_prop")
# -4.42 (final), -3.96 (meeting), -3.46 (semi-final), -3.20 (heat) => story holds

# 4.2. Null simulation (R2 point 6.5)
# -----------------------------------------------------------------------------
# Pacing shapes (section speed / entry mean speed) are permuted among entries
# of the same event x sex. Each entry keeps its own mean speed, athlete
# (S, E, season_id), context, sex and event. phi = 1 (tactical_effect = 0):
# positions are meaningless after permutation. PLf, floor and both LMMs are
# identical to the observed analysis. Output: observed coefficient / context
# slope vs null median and 2.5-97.5% band.
# -----------------------------------------------------------------------------

# 4.2.1. Base data: pre-PLf section data, same factor coding as all_fatigue

null_base <- bind_rows(long_800  %>% mutate(event = "800"),
                       long_1500 %>% mutate(event = "1500")) %>%
  mutate(event             = factor(event, levels = c("800", "1500")),
         Gender            = factor(Gender, levels = c("male", "female")),
         Competition_type3 = relevel(factor(Competition_type), ref = "Meeting"))

# 4.2.2. Simulated dataset for one replicate

simulate_null_data <- function(base_df, b, permute = TRUE, tactical_effect = 0) {
  set.seed(2026 + b)
  
  shp <- base_df %>%
    arrange(race_id, entry_id, split_index) %>%
    group_by(event, Gender, race_id, entry_id) %>%
    summarise(vbar  = mean(split_speed_m_s),
              shape = list(split_speed_m_s / mean(split_speed_m_s)), .groups = "drop") %>%
    group_by(event, Gender) %>%
    mutate(shape = if (permute) sample(shape) else shape) %>%   # permutation, not bootstrap
    ungroup() %>%
    mutate(v_sim = purrr::map2(vbar, shape, ~ .x * .y)) %>%
    select(race_id, entry_id, v_sim) %>%
    tidyr::unnest_longer(v_sim) %>%
    group_by(race_id, entry_id) %>%
    mutate(split_index = row_number()) %>%
    ungroup()
  
  sim <- base_df %>%
    select(-any_of(c("CV", "Speed_200"))) %>%
    left_join(shp, by = c("race_id", "entry_id", "split_index")) %>%
    mutate(split_speed_m_s = v_sim,
           split_time_s    = 100 / v_sim,
           position_change = if (permute) 0L else position_change) %>%
    select(-v_sim) %>%
    group_by(race_id, entry_id) %>%
    arrange(split_index, .by_group = TRUE) %>%
    mutate(CV        = 100 * sd(split_speed_m_s) / mean(split_speed_m_s),
           Speed_200 = mean(split_speed_m_s[1:2])) %>%
    ungroup()
  
  compute_dynamic_fatigue_drake(sim, tactical_effect = tactical_effect) %>%
    mutate(risk_ratio = split_speed_m_s / S_dyn_drake)
}

# 4.2.3. Validation: without permutation and with phi = 5%, the simulator must
# reproduce the observed S_dyn_drake exactly

chk <- simulate_null_data(null_base, b = 0, permute = FALSE, tactical_effect = 0.05) %>%
  select(race_id, entry_id, split_index, S_sim = S_dyn_drake) %>%
  inner_join(all_fatigue %>% select(race_id, entry_id, split_index, S_dyn_drake),
             by = c("race_id", "entry_id", "split_index"))
nrow(chk) == nrow(all_fatigue)                                            # must be TRUE
isTRUE(all.equal(chk$S_sim, chk$S_dyn_drake, check.attributes = FALSE))   # must be TRUE

# 4.2.4. Statistics extracted identically from observed and simulated fits

extract_null_stats <- function(m_rr, m_no) {
  sl_rr <- summary(emmeans::emtrends(m_rr, ~ Competition_type3, var = "race_prop",
                                     at = list(risk_ratio = c(0.4, 0.7, 1.0)),
                                     lmer.df = "asymptotic"))
  sl_no <- summary(emmeans::emtrends(m_no, ~ Competition_type3, var = "race_prop",
                                     lmer.df = "asymptotic"))
  bind_rows(
    tibble(model = "RR",   type = "coef",  term = names(lme4::fixef(m_rr)), value = unname(lme4::fixef(m_rr))),
    tibble(model = "noRR", type = "coef",  term = names(lme4::fixef(m_no)), value = unname(lme4::fixef(m_no))),
    tibble(model = "RR",   type = "slope", term = as.character(sl_rr$Competition_type3), value = sl_rr$race_prop.trend),
    tibble(model = "noRR", type = "slope", term = as.character(sl_no$Competition_type3), value = sl_no$race_prop.trend)
  )
}

# 4.2.5. One replicate: simulate, fit both models, extract, diagnostics

fit_null_once <- function(b, base_df) {
  d    <- simulate_null_data(base_df, b)
  ctrl <- lme4::lmerControl(optimizer = "bobyqa")
  
  m_rr <- lme4::lmer(S_dyn_drake ~ race_prop * risk_ratio * Competition_type3 + CV + Speed_200 + Gender + event +
                       (1 + race_prop * risk_ratio | season_id), data = d, control = ctrl)
  m_no <- lme4::lmer(S_dyn_drake ~ race_prop * Competition_type3 + CV + Speed_200 + Gender + event +
                       (1 + race_prop | season_id), data = d, control = ctrl)
  
  diag <- tibble(model = "diag", type = "diag",
                 term  = c("floor_rate", "singular_RR", "singular_noRR", "conv_warn_RR", "conv_warn_noRR"),
                 value = c(mean(d$floor_bound[d$split_index >= 2]),
                           lme4::isSingular(m_rr), lme4::isSingular(m_no),
                           length(m_rr@optinfo$conv$lme4$messages) > 0,
                           length(m_no@optinfo$conv$lme4$messages) > 0))
  
  bind_rows(extract_null_stats(m_rr, m_no), diag) %>% mutate(b = b)
}

# 4.2.6. Run: time 3 replicates first, then B = 100 (results saved to disk)

plan(multisession, workers = max(1, parallel::detectCores() - 1))
pkgs <- c("dplyr", "tidyr", "purrr", "tibble", "lme4", "emmeans")

t0 <- Sys.time()
null_test <- future_lapply(1:3, fit_null_once, base_df = null_base,
                           future.seed = TRUE, future.packages = pkgs)

null_test_df <- bind_rows(null_test)

null_test_df %>% filter(type == "diag") %>% tidyr::pivot_wider(names_from = term, values_from = value)
null_test_df %>% filter(model == "noRR", type == "slope")    # context slopes under the null

Sys.time() - t0   # x (100 / number of workers) ~ total run time

null_res <- future_lapply(1:100, fit_null_once, base_df = null_base,
                          future.seed = TRUE, future.packages = pkgs) %>% bind_rows()
#saveRDS(null_res, "null_simulation_B100.rds")
plan(sequential)

# null_res <- readRDS("null_simulation_B100.rds")   # reuse without rerunning

# 4.2.7. Observed vs null

obs_stats <- extract_null_stats(model_S_dyn_drake_relevel, model_S_noRR)

null_summary <- null_res %>%
  filter(type != "diag") %>%
  group_by(model, type, term) %>%
  summarise(null_median = median(value),
            null_lo     = quantile(value, 0.025),
            null_hi     = quantile(value, 0.975), .groups = "drop") %>%
  left_join(obs_stats %>% rename(observed = value), by = c("model", "type", "term")) %>%
  mutate(outside_null = observed < null_lo | observed > null_hi)

print(null_summary %>% filter(type == "coef",  model == "RR"),   n = Inf)   # Table 3 terms
print(null_summary %>% filter(type == "slope"),                  n = Inf)   # context slopes, both models

# Diagnostics across replicates (observed floor rate 2.6% / 4%)
null_res %>% filter(type == "diag") %>%
  group_by(term) %>%
  summarise(mean = mean(value), min = min(value), max = max(value))

# 4.2.8. Keeping tactical multiplier to 1 in both
obs_phi1 <- simulate_null_data(null_base, b = 0, permute = FALSE, tactical_effect = 0)
ctrl <- lmerControl(optimizer = "bobyqa")
m_rr_phi1 <- lmer(S_dyn_drake ~ race_prop * risk_ratio * Competition_type3 + CV + Speed_200 + Gender + event +
                    (1 + race_prop * risk_ratio | season_id), data = obs_phi1, control = ctrl)
m_no_phi1 <- lmer(S_dyn_drake ~ race_prop * Competition_type3 + CV + Speed_200 + Gender + event +
                    (1 + race_prop | season_id), data = obs_phi1, control = ctrl)
null_summary %>% filter(type == "slope") %>% select(model, term, null_lo, null_hi) %>%
  left_join(extract_null_stats(m_rr_phi1, m_no_phi1) %>% rename(obs_phi1 = value),
            by = c("model", "term")) %>%
  mutate(outside_null = obs_phi1 < null_lo | obs_phi1 > null_hi)

# =============================================================================
# 5. Relation to performance: sensitivity to the definition of S durability
# =============================================================================
# -----------------------------------------------------------------------------
# Goal: check the sensitivity of the Top-3 results to the definition of S durability,
#       using the S at the end of the race (last S_post, including the final section,
#       floor-bound in 59% of entries) instead of the start of the final section
#       (last S_dyn_drake, main analysis).
# Result: the end-of-race version reduced fit (dAIC +4 to +14.2 in three of four
#       strata) and attenuated Delta S associations; Fresh S and position ORs were
#       virtually unchanged. The main definition was retained.
# Same model specifications as the main analysis (incl. quadratic Fresh S in the
# 1500 m male model). All objects carry the suffix "_post" so that the main
# analysis objects are not overwritten. Requires: all_fatigue, top3_models,
# OR_df and extract_or() from the main script.
# -----------------------------------------------------------------------------

# 5.0. Dataset and variables

S_summary_post <- all_fatigue %>%
  group_by(race_id, entry_id, season_id, Competition_type, event, Gender, range_rank, position_last_lap) %>%
  summarize(
    fresh_S    = first(S, order_by = split_index),
    end_race_S = last(S_post, order_by = split_index),
    .groups = "drop"
  ) %>%
  mutate(
    Top3         = as.factor(ifelse(range_rank == "1-3rd", "Yes", "No")),
    S_drop_neg   = end_race_S - fresh_S,          # Delta S, end-of-race version (<= 0)
    position_neg = -position_last_lap
  )

# Share of entries whose last S_post is floor-bound (= speed of the final 100 m)
all_fatigue %>%
  group_by(race_id, entry_id) %>%
  filter(split_index == max(split_index)) %>%
  ungroup() %>%
  summarise(exhausted_at_finish = mean(exhausted_finish))

# 5.1. Data, stratified by event x sex

data_list_post <- list(
  "800_Male"    = S_summary_post %>% filter(event == "800",  Gender == "male"),
  "800_Female"  = S_summary_post %>% filter(event == "800",  Gender == "female"),
  "1500_Male"   = S_summary_post %>% filter(event == "1500", Gender == "male") %>%
    mutate(fresh_S_c = fresh_S - mean(fresh_S)),
  "1500_Female" = S_summary_post %>% filter(event == "1500", Gender == "female")
)

# Same samples as the main analysis
stopifnot(all(purrr::map2_lgl(data_list_post, top3_models, ~ nrow(.x) == nobs(.y))))

# 5.2. Models (same specifications as the main analysis)

ctrl_glmer <- glmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3, relTol = NULL))

top3_models_post <- list(
  "800_Male"    = glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_post$`800_Male`,    control = ctrl_glmer),
  "800_Female"  = glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_post$`800_Female`,  control = ctrl_glmer),
  "1500_Male"   = glmer(Top3 ~ fresh_S_c + I(fresh_S_c^2) + S_drop_neg + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_post$`1500_Male`,   control = ctrl_glmer),
  "1500_Female" = glmer(Top3 ~ fresh_S + S_drop_neg + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_post$`1500_Female`, control = ctrl_glmer)
)

purrr::map(top3_models_post, ~ .x@optinfo$conv$lme4$messages)
purrr::map(top3_models_post, isSingular)

# 5.3. Odds ratios, end-of-race version

OR_df_post <- purrr::imap_dfr(top3_models_post, extract_or) %>%
  filter(term %in% c("fresh_S", "S_drop_neg", "position_neg"))

# 5.4. Comparison with the main analysis

# (a) Delta S: OR per 1 m/s and per 1 SD (the two versions have different spreads)
sd_dS <- function(dl) purrr::map_dbl(dl, ~ sd(.x$S_drop_neg))

cmp_dS <- tibble(
  Model      = names(top3_models),
  OR_main    = OR_df$estimate[match(paste(names(top3_models), "S_drop_neg"), paste(OR_df$Model, OR_df$term))],
  OR_post    = OR_df_post$estimate[match(paste(names(top3_models), "S_drop_neg"), paste(OR_df_post$Model, OR_df_post$term))],
  SD_main    = sd_dS(data_list),
  SD_post    = sd_dS(data_list_post)
) %>%
  mutate(OR_perSD_main = OR_main^SD_main,
         OR_perSD_post = OR_post^SD_post)
print(cmp_dS)

# (b) Model fit on the same data (non-nested models, same sample: AIC comparable)
tibble(Model    = names(top3_models),
       AIC_main = purrr::map_dbl(top3_models, AIC),
       AIC_post = purrr::map_dbl(top3_models_post, AIC)) %>%
  mutate(dAIC_post_minus_main = AIC_post - AIC_main)

# (c) Fresh S and position: do they change when Delta S changes definition?
OR_df %>% filter(term %in% c("fresh_S", "position_neg"), !is.na(estimate)) %>%
  select(Model, term, OR_main = estimate) %>%
  left_join(OR_df_post %>% select(Model, term, OR_post = estimate), by = c("Model", "term"))

# 5.5. Sensitivity: Delta S expressed as % of Fresh S (same specifications, suffix "_pct")

data_list_pct <- purrr::map(data_list, ~ .x %>% mutate(dS_pct = 100 * S_drop_neg / fresh_S))   # <= 0, % of Fresh S

top3_models_pct <- list(
  "800_Male"    = glmer(Top3 ~ fresh_S + dS_pct + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_pct$`800_Male`,    control = ctrl_glmer),
  "800_Female"  = glmer(Top3 ~ fresh_S + dS_pct + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_pct$`800_Female`,  control = ctrl_glmer),
  "1500_Male"   = glmer(Top3 ~ fresh_S_c + I(fresh_S_c^2) + dS_pct + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_pct$`1500_Male`,   control = ctrl_glmer),
  "1500_Female" = glmer(Top3 ~ fresh_S + dS_pct + position_neg + (1 | season_id),
                        family = "binomial", data = data_list_pct$`1500_Female`, control = ctrl_glmer)
)

OR_df_pct <- purrr::imap_dfr(top3_models_pct, extract_or) %>%
  filter(term %in% c("fresh_S", "dS_pct", "position_neg"))

# Delta S per 1 SD (comparable across units) and AIC, main (m/s) vs %
tibble(
  Model        = names(top3_models),
  OR_perSD_ms  = OR_df$estimate[match(paste(names(top3_models), "S_drop_neg"), paste(OR_df$Model, OR_df$term))]^
    purrr::map_dbl(data_list, ~ sd(.x$S_drop_neg)),
  OR_perSD_pct = OR_df_pct$estimate[match(paste(names(top3_models), "dS_pct"), paste(OR_df_pct$Model, OR_df_pct$term))]^
    purrr::map_dbl(data_list_pct, ~ sd(.x$dS_pct)),
  p_pct        = OR_df_pct$p_value[match(paste(names(top3_models), "dS_pct"), paste(OR_df_pct$Model, OR_df_pct$term))],
  AIC_ms       = purrr::map_dbl(top3_models, AIC),
  AIC_pct      = purrr::map_dbl(top3_models_pct, AIC)
)

# =============================================================================
# 6. Sensitivity of the Top-3 models to the feasibility floor
# =============================================================================

# 6.1. Flooring hit descriptive stats (Run after section 7.0 [S_summary$Top3 exists]
last_floor <- all_fatigue %>%
  group_by(race_id, entry_id) %>%
  filter(split_index == max(split_index)) %>%
  ungroup() %>%
  select(race_id, entry_id, event, Gender, floor_last = floor_bound) %>%
  left_join(S_summary %>% select(race_id, entry_id, Top3), by = c("race_id", "entry_id"))

last_floor %>%
  group_by(event, Gender, Top3) %>%
  summarise(n = n(), pct_floor_last = 100 * mean(floor_last), .groups = "drop")

last_floor %>%
  group_by(event) %>%
  summarise(n = n(), pct_floor_last = 100 * mean(floor_last), .groups = "drop")

last_floor %>%
  group_by(Top3) %>%
  summarise(n = n(), pct_floor_last = 100 * mean(floor_last), .groups = "drop")

# 6.3. Fitting without floor hits
data_list_nofloor <- purrr::map(data_list, ~ .x %>%
                                  anti_join(last_floor %>% filter(floor_last), by = c("race_id", "entry_id")))
top3_nofloor <- purrr::map2(top3_models, data_list_nofloor, ~ update(.x, data = .y))
purrr::imap_dfr(top3_nofloor, extract_or) %>%
  filter(term %in% c("fresh_S", "S_drop_neg", "position_neg"))

# The 1500 m male Fresh S terms were filtered out of the exclusion output (they are fresh_S_c terms)
summary(top3_nofloor$`1500_Male`)$coefficients
#####################################################################################################
#                                                                                                   #
#                 APPENDIX: Fatigue-integrated Power Law (PLf) and Speed Durability                 #
#                 -----------------------------------------------------------------                 #
#                                                                                                   #
#   0. Load the libraries and datasets                                                              #
#   1. Setting the tactical multiplier (phi) value                                                  #
#   2. Feasibility floor hitting                                                                    #
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

library(dplyr); library(purrr); library(tibble)

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
  tmp <- compute_dynamic_fatigue_drake(df, tactical_effect = tactical_effect)
  tmp <- tmp %>%
    group_by(race_id, entry_id) %>%
    mutate(
      cum_distance = cumsum(split_time_s * 0 + distance_m),  # verify distance_m exists; adjust if named differently
      remaining_distance = pmax(total_distance - cum_distance, 0.0001),
      predicted_remaining_time = predict_remaining_time(S_dyn_drake, b, remaining_distance),
      predicted_total_time = cumsum(split_time_s) + predicted_remaining_time
    ) %>%
    ungroup() %>%
    mutate(abs_rel_error = abs(predicted_total_time - Final_time) / Final_time)
  
  mean(tmp$abs_rel_error, na.rm = TRUE)
}

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

# 1.4. Grouped k-fold CV (grouped by entry_id)

run_cv <- function(df, total_distance, athlete_id_col, k = 5, grid = tactical_grid) {
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

cv_800  <- run_cv(long_800,  800,  athlete_id_col = "entry_id") 
cv_1500 <- run_cv(long_1500, 1500, athlete_id_col = "entry_id")

# 1.5. Manual 5% k-fold CV (grouped by entry_id)

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

fixed_cv_800  <- manual_cv_fixed(long_800,  800,  athlete_id_col = "entry_id")
fixed_cv_1500 <- manual_cv_fixed(long_1500, 1500, athlete_id_col = "entry_id")

# 1.6. Compare
print(insample_800);  print(cv_800)
print(insample_1500); print(cv_1500)

print(fixed_cv_800);  mean(fixed_cv_800$MARE_at_5pct);  sd(fixed_cv_800$MARE_at_5pct)
print(fixed_cv_1500); mean(fixed_cv_1500$MARE_at_5pct); sd(fixed_cv_1500$MARE_at_5pct)

# =============================================================================
# 2. Feasibility floor hitting
# =============================================================================

# 2.1. Compute the diagnostic function

compute_dynamic_fatigue_drake_diag <- function(df, tactical_effect = 0.00, eps = 1e-8) {
  df <- df %>% dplyr::ungroup() %>% dplyr::arrange(race_id, entry_id, split_index)
  
  df %>%
    dplyr::group_by(race_id, entry_id) %>%
    dplyr::group_split() %>%
    lapply(function(g) {
      n <- nrow(g)
      S0 <- g$S[1]; E <- g$E[1]
      alpha <- 1 / max(1 - E, eps)
      v  <- pmax(g$split_speed_m_s, 0)
      dt <- pmax(g$split_time_s, 0)
      v_alpha <- v^alpha
      S0 <- max(S0, v[1])
      Sdyn <- numeric(n); U <- numeric(n)
      floor_bound <- logical(n)  # NEW: TRUE if the U-recursion was overridden
      Sdyn[1] <- S0; U[1] <- S0^alpha; floor_bound[1] <- NA
      
      for (i in 1:n) {
        if (i < n) {
          phi <- if (g$position_change[i] >= 2) 1 + tactical_effect else
            if (g$position_change[i] <= -2) 1 - tactical_effect else 1
          tryU <- U[i] - dt[i] * v_alpha[i] * phi
          fallback_triggered <- !(is.finite(tryU) && tryU > 0)
          nextU <- if (!fallback_triggered) min(tryU, U[i]) else v[i+1]^alpha
          S_next <- max(nextU^(1/alpha), v[i+1])
          floor_bound[i+1] <- fallback_triggered || isTRUE(all.equal(S_next, v[i+1]))
          if (S_next > Sdyn[i] && S_next - Sdyn[i] < 1e-8) S_next <- Sdyn[i]
          Sdyn[i+1] <- S_next
          U[i+1] <- Sdyn[i+1]^alpha
        }
      }
      g$S_dyn_drake <- Sdyn
      g$floor_bound <- floor_bound
      g
    }) %>% dplyr::bind_rows()
}

# 2.2. Floor-binding rate across the full tactical_grid
# -----------------------------------------------------------------------------
# Does the tactical multiplier affect the proportion of feasibility floor hit?
# -----------------------------------------------------------------------------

floor_rate_curve <- function(df, grid = tactical_grid) {
  map_dfr(grid, function(te) {
    diag <- compute_dynamic_fatigue_drake_diag(df, tactical_effect = te)
    tibble(tactical_effect = te, floor_rate = mean(diag$floor_bound, na.rm = TRUE))
  })
}

floor_curve_800  <- floor_rate_curve(long_800)
floor_curve_1500 <- floor_rate_curve(long_1500)

print(floor_curve_800)
print(floor_curve_1500)

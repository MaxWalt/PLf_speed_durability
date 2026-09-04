############################################################
#                                                          #
#          Pacing datasets: Power law computation          #
#                                                          #
#     1. Prepare dataset                                   #
#     2. Compute PL                                        #
#     3. Save                                              #
#                                                          #
#                                 Maxime Walt, 2025        #
#                                                          #
############################################################

#=== 0. Load libraries and dataset

# 0.1. Load libraries

library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(table1)
library(metrica)
library(performance)
library(ggpubr)
library(digest)
library(ggtext)
library(stringr)
library(stringi)

# 0.2. Load datasets

# WA scraped database
load("C:/Users/maxim/Documents/UNIL - Doctorat/1_Research_data/1_Data_collection/1_Scraping/Script R_Extractor 5000/WA_21_25_all")

# NA's manual check
PL_NA <- read_excel("C:/Users/maxim/Documents/UNIL - Doctorat/1_Research_data/3_Study_B_Pacing_profiles/Pacing_profiles/NA_PL.xlsx")

#=== 1. Prepare dataset 

# 1.1. Work on WA database 

# Normalize spaces

normalize_spaces <- function(x) {
  x %>%
    stri_replace_all_charclass("\\p{Zs}", " ") %>%  # all Unicode spaces ??? space
    str_replace_all("[\\t\\r\\n]+", " ") %>%        # tabs/newlines ??? space
    str_squish()                                     # trim & collapse
}

WA_21_25_all$Competitor <- normalize_spaces(WA_21_25_all$Competitor)

# 1.2. Compute new variables and filter useless ones

WA_21_25_all <- WA_21_25_all %>%
  # drop columns only if present
  select(-any_of(c("WIND", "V8"))) %>%
  # normalize a few names (only if they exist)
  rename(
    Athletes = any_of("Competitor"),
    Nat      = any_of("V6"),
    time     = any_of("Mark_seconds")
  ) %>%
  # compute field, ID, distance and speed
  mutate(
    field = case_when(
      str_ends(event_name, "_R") ~ "road",
      str_ends(event_name, "_i") ~ "indoor",
      TRUE                       ~ "track"
    ),
    ID       = paste(Athletes, Nat, season, sep = "_"),
    # distance is the numeric prefix of event_name (e.g., "400", "1609", "42195_R")
    distance = as.numeric(str_extract(event_name, "^[0-9]+")),
    # speed in m/s (requires numeric time in seconds)
    speed    = distance / as.numeric(time),
    log_speed = log(speed),
    log_time = log(time)
  )

WA_SB <- WA_21_25_all %>%
  select(ID, Athletes, gender, distance, field, time, speed, log_speed, log_time)

# 1.2. Work on manual check

PL_NA <- subset(PL_NA, PL_NA$N_event >= 2)

# Select variables of interest

PL_NA <- PL_NA[, c(1:3, 7:24)]

# Compute speed

PL_NA <- PL_NA %>%
  mutate(
    speed_400_track = 400 / time_400,
    speed_600_track = 600 / time_600,
    speed_800_track = 800 / time_800,
    speed_1000_track = 1000 / time_1000,
    speed_1500_track = 1500 / time_1500,
    speed_1609_track = 1609 / time_1609,
    speed_1609_road = 1609 / time_1609_R,
    speed_2000_track = 2000 / time_2000,
    speed_3000_track = 3000 / time_3000, 
    speed_3218_track = 3218 / time_3218,
    speed_5000_track = 5000 / time_5000,
    speed_5000_road = 5000 / time_5000_R,
    speed_10000_track = 10000 / time_10000,
    speed_10000_road = 10000 / time_10000_R,
    speed_15000_road = 15000 / time_15000_R,
    speed_16090_road = 16090 / time_16090_R,
    speed_21095_road = 21095 / time_21095_R,
    speed_42195_road = 42195 / time_42195_R
  )

# Delete all NA speed

PL_NA <- PL_NA %>%
  filter(if_any(starts_with("speed_"), ~ !is.na(.)))

# Pivot in long format

# PIVOT TIME
time_long <- PL_NA %>%
  pivot_longer(
    cols = starts_with("time_"),
    names_to = "raw",
    values_to = "time"
  ) %>%
  mutate(
    field = if_else(str_detect(raw, "_R$"), "road", "track"),
    distance = str_remove(raw, "_R$"),
    distance = str_remove(distance, "time_")
  ) %>%
  select(ID, Athletes, Gender, distance, field, time) %>%
  distinct()

# PIVOT SPEED
speed_long <- PL_NA %>%
  pivot_longer(
    cols = starts_with("speed_"),
    names_to = "raw",
    values_to = "speed"
  ) %>%
  separate(raw, into = c("discard", "distance", "field"), sep = "_", fill = "right", remove = TRUE) %>%
  mutate(field = if_else(is.na(field), "track", field)) %>%
  select(ID, Athletes, Gender, distance, field, speed) %>%
  distinct()

# JOIN and CLEAN
PL_NA_long <- inner_join(time_long, speed_long,
                             by = c("ID", "Athletes", "Gender", "distance", "field")) %>%
  mutate(distance = as.numeric(distance)) %>%
  arrange(ID, distance)

PL_NA_long <- subset(PL_NA_long, !is.na(PL_NA_long$speed))

# Each ID must have at leat 2 rows

PL_NA_long <- PL_NA_long %>%
  group_by(ID) %>%
  filter(n() >= 2) %>%
  ungroup()

# Compute log-transformed speed and duration for the Power Law model and scale the Best_perf
PL_NA_long <- PL_NA_long %>%
  mutate(
    log_speed = log(speed),
    log_time = log(time)
  )

# Rename gender to match str(WA_SB)

names(PL_NA_long)[3] <- c("gender")

# 1.3. Combine, delete duplicates and save

# Combine

WA_SB <- rbind(WA_SB, PL_NA_long)

# Delete duplicates

WA_SB <- WA_SB %>%
  arrange(ID, distance) %>%
  mutate(
    ID2 = paste(ID, distance, field, time)
  )

WA_SB <- WA_SB[!duplicated(WA_SB$ID2),]

WA_SB$ID2 <- NULL

# Save

#save(WA_SB, file = "WA_SB")

# 1.4. Choose your field!

# If an athlete has both indoor and track performances for the same distance, we keep only the fastest one and set the other one to NA
# If an athlete has both road and track performances for the same distance, we keep only the track one and set the road one to NA
# If an athlete has both road and indoor performances for the same distance, we keep only the indoor one and set the road one to NA

na_out_by_field <- function(df) {
  df %>%
    mutate(
      # Priority values (lower = better)
      priority = case_when(
        field == "indoor" ~ 0L,
        field == "track"  ~ 1L,
        TRUE              ~ 2L
      ),
      time_rank = ifelse(is.na(time), Inf, time)
    ) %>%
    group_by(ID, distance) %>%
    mutate(
      # Find the "best" row index per group
      best_row = {
        # If track exists, best_row is the fastest track
        if (any(field == "track")) {
          which.min(ifelse(field == "track", time_rank, Inf))
          # Else if indoor exists, fastest indoor wins
        } else if (any(field == "indoor")) {
          which.min(ifelse(field == "indoor", time_rank, Inf))
          # Else pick fastest overall
        } else {
          which.min(time_rank)
        }
      },
      # Blank everything except the best row
      time  = if_else(row_number() != best_row, NA_real_, time),
      speed = if_else(row_number() != best_row, NA_real_, speed)
    ) %>%
    ungroup() %>%
    select(-priority, -time_rank, -best_row)
}

# use on your WA_SB
WA_SB_best <- na_out_by_field(WA_SB)

WA_SB_best <- WA_SB_best[!is.na(WA_SB_best$speed), ]

# 1.4. Check data logic

# For each ID, if this speed is smaller than the speed for ANY bigger distance, set to NA

WA_SB_best <- WA_SB_best %>%
  group_by(ID) %>%
  group_modify(~ {
    df <- .x
    
    # For each row i in this athlete's subset:
    df$speed <- sapply(seq_len(nrow(df)), function(i) {
      this_dist  <- df$distance[i]
      this_speed <- df$speed[i]
      
      # Look at speeds in bigger distances for this same athlete
      bigger_speeds <- df$speed[df$distance > this_dist]
      
      # If no bigger distance exists, keep current speed
      if (length(bigger_speeds) == 0) {
        return(this_speed)
      } 
      
      # Otherwise, if the current speed is less than 
      # the max speed of those bigger distances, set to NA
      bigger_max <- max(bigger_speeds, na.rm = TRUE)
      if (this_speed < bigger_max) NA_real_ else this_speed
    })
    
    df
  }) %>%
  ungroup()

WA_SB_best <- WA_SB_best[!is.na(WA_SB_best$speed), ]

# Each ID must have at leat 2 rows

WA_SB_best <- WA_SB_best %>%
  group_by(ID) %>%
  filter(n() >= 2) %>%
  ungroup()

# Ensure positive & finite
WA_SB_best <- WA_SB_best %>%
  filter(is.finite(time), time > 0, is.finite(speed), speed > 0)

#=== 2. Compute PL

results_PL_WA_SB <- WA_SB_best %>%
  group_by(ID) %>%
  do({
    pl_lr <- lm(log_speed ~ log_time, data = .)
    fit_summary <- summary(pl_lr)
    
    intercept <- coef(pl_lr)[1]
    slope <- coef(pl_lr)[2]
    
    S_lm_ALL <- exp(intercept)
    E_lm_ALL <- slope + 1
    
    a <- exp(intercept)
    b <- -slope
    
    # Goodness-of-fit metrics
    r_squared     <- fit_summary$r.squared
    adj_r_squared <- fit_summary$adj.r.squared
    residual_se   <- fit_summary$sigma
    n_performances <- nrow(.)
    df_resid       <- df.residual(pl_lr)
    
    data_with_predictions <- mutate(., pl_speed = a / (time^b))
    predicted_performance <- mutate(.,
                                    predicted_time = (a / distance)^(1 / (b - 1)),
                                    predicted_speed = distance / predicted_time
    )
    
    tibble(
      ID = unique(.$ID),
      intercept = intercept, slope = slope,
      S = S_lm_ALL, E = E_lm_ALL, a = a, b = b,
      r_squared = r_squared, adj_r_squared = adj_r_squared,
      residual_se = residual_se,
      n_performances = n_performances, df_resid = df_resid,
      data_with_predictions = list(data_with_predictions),
      predicted_performance = list(predicted_performance)
    )
  }) %>%
  ungroup()

results_PL_WA_SB_unnested <- results_PL_WA_SB %>%
  select(ID, S, E, a, b, predicted_performance) %>%
  unnest(c(predicted_performance), .id = "ID") %>% 
  select(ID, Athletes, gender, S, E, a, b, distance, field, time, speed, log_speed, log_time, predicted_time, predicted_speed)

results_PL_WA_SB_stats <- results_PL_WA_SB %>%
  dplyr::select(ID, S, E, a, b, predicted_performance, r_squared, residual_se, n_performances) %>%
  unnest(c(predicted_performance), .id = "ID") %>% 
  dplyr::select(ID, Athletes, gender, S, E, a, b, distance, field, time, speed, log_speed, log_time, predicted_time, predicted_speed, r_squared, residual_se, n_performances)

# Merge the results

WA_SB_best_PL <- left_join(WA_SB_best, results_PL_WA_SB_unnested[, c(1:9, 14)], by = c("ID", "Athletes", "gender", "distance", "field"))

#=== 3. Save

save(WA_SB_best, file = "WA_SB_best")
save(WA_SB_best_PL, file = "WA_SB_best_PL")
save(results_PL_WA_SB, file = "results_PL_WA_SB")
save(results_PL_WA_SB_stats, file = "results_PL_WA_SB_stats")
save(results_PL_WA_SB_unnested, file = "results_PL_WA_SB_unnested")



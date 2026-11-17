## PLf_speed_durability
Quantify real-time speed durability and examine how fresh speed, speed durability and positioning correlate with performance.

_The full code and dataset are part of a Submission in International Journal of Sports Physiology and Performance (IJSPP.2026-0263), for the PhD thesis of Maxime Walt._

### Data
#### Pacing datasets
Pacing_splits_PHD_800m and Pacing_splits_PHD_1500m contain the 100m section times of Diamond League and Major Championship races for the 800m and 1500m respectively.
#### PL parameters datasets
WA_SB_best_PL and results_PL_WA_SB_unnested are big athletes databases that retrieve: SBs, PL parameters and PL predictions respectively.

### Scripts
`PL_computation_WA` : used to compute the PL parameters databased (WA_SB_best_PL, results_PL_WA_SB, results_PL_WA_SB_unnested and results_PL_WA_SB_stats). Final databases join, raw datasets not, script attached as informative.

`PLf_speed_durability` : main analyses script

`PLf_speed_durability_APPENDIX` : appendix analyses

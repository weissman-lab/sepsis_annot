# This file create a visualization summary of a particular patient 
# during  a period around the time when they might have sepsis

# -------------------------
SEPSIS_DATA_PATH = "/data/uphs_sepsis/single_file_Sep_2019/single_file_v10-21-2019.csv"
SEPSIS3_EVENTS = "sepsis_cases_sepsis3.csv"
# -------------------------

set.seed(24601)

library(data.table)
library(ggplot2)

sep_all <- fread(SEPSIS_DATA_PATH)
sep_events <- fread(SEPSIS3_EVENTS)


trial_mrn <- sample(sep_events$mrn, size = 1)
trial_ts <- sep_events[mrn == trial_mrn]

make_viz <- function(MRN, center_time) {



}

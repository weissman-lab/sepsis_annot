# This file create a visualization summary of a particular patient 
# during  a period around the time when they might have sepsis

# -------------------------
SEPSIS_DATA_PATH = "/data/uphs_sepsis/single_file_Sep_2019/single_file_v10-21-2019.csv"
SEPSIS3_EVENTS = "sepsis_cases_sepsis3.csv"
# -------------------------

set.seed(24601)

library(data.table)
library(ggplot2)
library(patchwork)

sep_all <- fread(SEPSIS_DATA_PATH)
sep_events <- fread(SEPSIS3_EVENTS)

# Get a random MRN and onset time
trial_mrn <- sample(sep_events$mrn, size = 1)
trial_dt <- sep_events[mrn == trial_mrn]
onset_char_ts <- NULL
if (! any(is.na(trial_dt$sofa_bcx)) & all(trial_dt$sofa_bcx == 1)) onset_char_ts <- trial_dt$bcx_datetime[1]
if (! any(is.na(trial_dt$sofa_abx)) & all(trial_dt$sofa_abx == 1)) onset_char_ts <- trial_dt$abx_datetime[1]
if (! any(is.na(trial_dt$sofa_sofa)) & all(trial_dt$sofa_sofa == 1)) onset_char_ts <- trial_dt$sofa_datetime[1]

if (is.null(onset_char_ts)) print('Warning: No sepsis onset identified') 

onset_ts <- as.POSIXct(strptime(onset_char_ts, format = '%Y-%m-%dT%H:%M:%S'))


# Now make a plot
make_viz <- function(pat_mrn, center_time) {

        temp_dt <- sep_all[MRN == pat_mrn]
        time_idx <- which(temp_dt$PD_BEG_TIMESTAMP == center_time)
       print(paste('Patient has', nrow(temp_dt), 'rows'))
       print(paste('Sepsis-3 onset is at hour', time_idx ))
       # Randomize endpoints so Sep-3 onset isn't exactly in the middle
       lb <- time_idx - 24 - round(runif(1, 2, 8))
       ub <- time_idx + 24 + round(runif(1, 2, 8))
       # Get range
       temp_dt <- temp_dt[lb:ub, sub_hour := lb:ub]
       # Plot
       p_hr <- ggplot(temp_dt, aes(sub_hour, HEART_RATE)) +
        geom_point() + geom_line() + theme_bw() + xlab('hour')
       p_sbp <- ggplot(temp_dt, aes(sub_hour, SYSTOLIC_BP)) +
        geom_point() + geom_line() + theme_bw() + xlab('hour')
       p_t <- ggplot(temp_dt, aes(sub_hour, TEMPERATURE)) +
        geom_point() + geom_line() + theme_bw() + xlab('hour')
       p_rr <- ggplot(temp_dt, aes(sub_hour, RESPIRATORY_RATE)) +
        geom_point() + geom_line() + theme_bw() + xlab('hour')
       p_hr + p_sbp + p_t + p_rr + 
            plot_layout(ncol = 1, guides = 'collect')
        
}

# Test it out
make_viz(trial_mrn, onset_ts)


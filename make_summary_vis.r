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
       
       plot_title <- with(temp_dt, paste(AGE[1], 'year-old',
                     RACE_1[1], GENDER[1], '\nadmitted for',
                     ADMIT_DX_DESC[1]))

       # Randomize offsets so Sep-3 onset is at the RIGHT side
       # This is biologically plausible based on current definition
       # i.e. unlikely that relevant antibiotic starttsime would be 
       # AFTER Sepsis-3 onset given organ dysfunction has already occurred
       lb <- max(1, time_idx - 48 - round(runif(1, 2, 8)))
       ub <- min(nrow(temp_dt), time_idx + 4 + round(runif(1, 2, 8)))
       
       # Get range
       temp_dt <- temp_dt[lb:ub, sub_hour := lb:ub]

       # Plot function
       make_series <- function(dt, val) {
       		   ggplot(dt[! is.na(dt[[val]])], aes_string('sub_hour', val)) +
		   	      geom_point() + geom_line(na.rm = TRUE) +
			      theme_bw() + xlab('hour')
       }

       # List of variable names to plot in order
       feat_vec <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
       		'RESPIRATORY_RATE')		
       
       # Generate and combine plots
       pl_lst <- lapply(feat_vec, function(val) make_series(temp_dt, val)) 
       wrap_elements(Reduce(`+`, pl_lst) + 
        plot_layout(ncol = 1, guides = 'collect')) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(hjust = 0.5))
        
}

# Test it out
make_viz(trial_mrn, onset_ts)


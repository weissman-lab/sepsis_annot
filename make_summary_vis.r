# This file create a visualization summary of a sample of patients 
# during  a period around the time when they might have sepsis
# It will save the de-identified output to an image file
# And generate a crosswalk csv file with PHI

# --------------------------------------------------------------------
SEPSIS_DATA_PATH = "/data/uphs_sepsis/single_file_Sep_2019/single_file_v10-21-2019.csv"
SEPSIS3_EVENTS = "sepsis_cases_sepsis3.csv"
NUM_CASES = 10
# --------------------------------------------------------------------

set.seed(24601)

library(data.table)
library(ggplot2)
library(plotly)
library(dplyr)

sep_all <- fread(SEPSIS_DATA_PATH)
sep_events <- fread(SEPSIS3_EVENTS)

# Identify sepsis onset timestamp
sep_onset_dt <- sep_events[, .(onset_char_ts = dplyr::case_when(all(sofa_bcx == 1) ~ bcx_datetime[1],
                                                               all(sofa_abx == 1) ~ abx_datetime[1],
                                                               all(sofa_sofa == 1) ~ sofa_datetime[1]),
                               mrn = mrn[1],
                               antbx_start_char_ts = abx_datetime[1],
                               mortality = mortality[1],
                               adm_datetime = adm_datetime[1]), 
                           by = visit]

# Identify sepsis events that occur at least 72 hours into the hospitalization
sep_onset_dt[, sepsis_onset_dt := as.POSIXct((strptime(onset_char_ts, format = '%Y-%m-%dT%H:%M:%S')))]
sep_onset_dt[, hosp_adm_dt := as.POSIXct((strptime(adm_datetime, format = '%Y-%m-%dT%H:%M:%S')))]
sep_onset_dt[, hosp_adm_to_sepsis_onset_hrs := as.numeric(sepsis_onset_dt - hosp_adm_dt, units = 'hours')]
sep_onset_dt <- sep_onset_dt[hosp_adm_to_sepsis_onset_hrs >= 72]

# Get sampling of cases
cases_dt <- sep_onset_dt[sample(.N, NUM_CASES)]

# Some data clean up
sep_all[RESPIRATORY_RATE == 'WDL', RESP_RATE := 18]
sep_all[RESPIRATORY_RATE == '', RESP_RATE := NA]
sep_all[! is.na(as.numeric(RESPIRATORY_RATE)),
        RESP_RATE := as.numeric(RESPIRATORY_RATE)]
                                                                       
# Now make a function to generate a plot for each user
make_viz <- function(pat_mrn, onset_time) {
  
  temp_dt <- sep_all[MRN == pat_mrn]
  time_idx <- which(temp_dt$PD_BEG_TIMESTAMP == onset_time)
  
  plot_title <- with(temp_dt, paste(AGE[1], 'year-old',
                                    # RACE_1[1], GENDER[1], # TOO MUCH BIAS RISK
                                    'admitted for',
                                    ADMIT_DX_DESC[1]))
  
  # Randomize offsets so Sep-3 onset is at the RIGHT side
  # This is biologically plausible based on current definition
  # i.e. unlikely that relevant antibiotic starttsime would be 
  # AFTER Sepsis-3 onset given organ dysfunction has already occurred
  lb <- max(1, time_idx - 48 - round(runif(1, 2, 8)))
  ub <- min(nrow(temp_dt), time_idx + 4 + round(runif(1, 2, 8)))
  
  # Get range
  temp_dt <- temp_dt[lb:ub][, sub_hour := lb:ub]
  
  # List of variable names to plot in order
  feat_vec <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
                'RESP_RATE', 'LACTATE_RESULT')
  normal_ranges <- data.table(measure = feat_vec,
                              norm_lo = c(60,  90,  97, 10, 0),
                              norm_hi = c(100, 120, 99, 20, 2))
  
  # Prepare data
  temp_dt <- temp_dt[, c('sub_hour', feat_vec), with = FALSE]
  temp_dt[, hour := sub_hour][, sub_hour := NULL]
  temp_dt_m <- melt(temp_dt, id.vars = 'hour')
  temp_dt_m <- merge(temp_dt_m, normal_ranges, 
                     by.x = 'variable', by.y = 'measure',
                     all.x = TRUE)
  
  # Make plot
  pp <- ggplot(temp_dt_m, aes(hour, value)) +
    geom_line(data = temp_dt_m[! is.na(value)]) +
    theme_bw() + 
    scale_x_continuous('hour', 
                       limits = c(lb, ub),
                       expand = c(.03,.03),
                       breaks = scales::pretty_breaks(n = round((ub - lb) / 10))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    geom_rect(aes(xmin = lb, xmax = ub, ymin = norm_lo, ymax = norm_hi),
             color = 'gray', alpha = 0.2) +
    geom_point() + 
    # ^^ NB putting this point layer at the end 
    # makes the points appear on "top" so the spike lines are more easily activated
    facet_wrap(~ variable, ncol = 1, scales = 'free_y') +
    ggtitle(plot_title) 
  
  return(ggplotly(pp) %>% 
           layout(xaxis = list(showspikes = TRUE)) %>% 
           config(displayModeBar = FALSE))
}

# Test it out
images_list <- lapply(1:NUM_CASES, function(i) make_viz(cases_dt[i]$mrn, cases_dt[i]$sepsis_onset_dt))

# Save images to file
img_path <- 'images/'
invisible(lapply(1:NUM_CASES, function(i) ggsave(filename = paste0(img_path, 'image', i, '.png'), 
                                       plot = images_list[[i]],
                                       width = 6,
                                       height = 8)))


# Save crosswalk data to file - NB. CONTAINS PHI!
cases_dt[, filename := paste0(img_path, 'image', .I, '.png')]
fwrite(cases_dt, 'crosswalk_sepsis.csv')


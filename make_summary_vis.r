# This file create a visualization summary of a sample of patients 
# during  a period around the time when they might have sepsis
# It will save the de-identified output to an image file
# And generate a crosswalk csv file with PHI

# --------------------------------------------------------------------
SEPSIS_DATA_PATH = "/data/uphs_sepsis/updated_combo_Jan_2020/single_file_augmented_data_Jan_2020.csv"
SEPSIS3_EVENTS = "sepsis_cases_sepsis3.csv"
NUM_CASES = 8
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

# Identify sepsis events that occur at least 48 hours into the hospitalization
sep_onset_dt[, sepsis_onset_dt := as.POSIXct(onset_char_ts, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, hosp_adm_dt := as.POSIXct(adm_datetime, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, antbx_start_dt := as.POSIXct(antbx_start_char_ts, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, hosp_adm_to_sepsis_onset_hrs := as.numeric(sepsis_onset_dt - hosp_adm_dt, 
                                                          units = 'hours')]
sep_onset_dt[, antbx_delay_hrs := as.numeric(antbx_start_dt - sepsis_onset_dt,
                                             units = 'hours')]
sep_onset_dt <- sep_onset_dt[antbx_delay_hrs <= 48]
sep_onset_dt <- sep_onset_dt[hosp_adm_to_sepsis_onset_hrs >= 48]

# Generate some features
sep_all[, adm_dt := as.POSIXct(ADMISSION_DTTM, format = '%Y-%m-%d %H:%M:%S')]
sep_all[, disch_dt := as.POSIXct(DISCHARGE_DTTM, format = '%Y-%m-%d %H:%M:%S')]
sep_all[, hosp_los_hrs := as.numeric(disch_dt - adm_dt, units = 'hours')]
  
# Some data clean up
sep_all[RESPIRATORY_RATE == 'WDL', RESP_RATE := 18]
sep_all[RESPIRATORY_RATE == '', RESP_RATE := NA]
sep_all[! is.na(as.numeric(RESPIRATORY_RATE)),
        RESP_RATE := as.numeric(RESPIRATORY_RATE)]

# Get only the updated creatinine
# helper function
get_first_in_seq <- function(x) {
  results <- NULL
  rr <- rle(x)
  gr <- length(rr$lengths)
  res_list <- lapply(1:gr, function(i) { 
   if(is.na(rr$values[i])) {
     return(NA)
   } else {
     return(c(rr$values[i], rep(NA, rr$lengths[i] - 1)))
   }
  })
  return(unlist(res_list))
}

sep_all <- sep_all[order(PAT_ENC_CSN, PD_BEG_TIMESTAMP)][, creat_first_obs := get_first_in_seq(SOFA_CREATININE_RESULT), 
                                                   by = PAT_ENC_CSN]
sep_all <- merge(sep_all, 
                 sep_onset_dt[, .(PAT_ENC_CSN = visit, sepsis_onset_dt)],
                 by = 'PAT_ENC_CSN', 
                 all.x = TRUE)
# Identify which cases have at least some data:
required_fields <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
                     'RESP_RATE', 'LACTATE_RESULT', 'SOFA_RESP_SPO2', 
                     'creat_first_obs', 'max_WBC', 'SOFA_CNS_GLASGOW_RESULT')
sep_all[, has_all_fields := all(unlist(lapply(.SD[seq(max(1, which(PD_BEG_TIMESTAMP == sepsis_onset_dt) - 48),
                                                      min(.N,which(PD_BEG_TIMESTAMP == sepsis_onset_dt) + 48)),
                                                  required_fields, with = FALSE], 
                                              function(x) any(! is.na(x))))), by = PAT_ENC_CSN]

# Get sampling of cases
# Per Rebecca's suggestion, sample in a stratified way across antibiotic delays
sep_onset_dt <- sep_onset_dt[, antbx_delay_cat := cut(antbx_delay_hrs, breaks = c(-1, 6, 12, 24, 49))]
N_PER_GROUP <- round(NUM_CASES / length(unique(sep_onset_dt$antbx_delay_cat)))
cases_dt <- sep_onset_dt[visit %in% sep_all[has_all_fields == TRUE]$PAT_ENC_CSN, 
                         .SD[sample(.N, N_PER_GROUP)], by = antbx_delay_cat]

# Initialize bounds
cases_dt[, `:=`(window_lower = -1, window_upper = -1)]

# Now make a function to generate a plot for each user
make_viz <- function(this_visit, onset_time, img_idx) {
  
  temp_dt <- sep_all[PAT_ENC_CSN == this_visit][order(PD_BEG_TIMESTAMP)]
  time_idx <- which(temp_dt$PD_BEG_TIMESTAMP == onset_time)
  
  plot_title <- with(temp_dt, paste0(AGE[1], '-year-old',
                                    # RACE_1[1], GENDER[1], # TOO MUCH BIAS RISK
                                    ' admitted for\n',
                                    ADMIT_DX_DESC[1]))
  
  # If admission diagnosis is very long, truncate it after the first comma
  truncate_title_if_needed <- function(x) {
    if (nchar(x) > 65) {
      first_comma <- regexpr(',', x)
      new_string <- substr(x, 1, first_comma - 1)
      return(new_string)
    } else {
      return(x)
    }
  }
  
  plot_title <- truncate_title_if_needed(plot_title)
  
  # Randomize offsets so Sep-3 onset is at the RIGHT side
  # This is biologically plausible based on current definition
  # i.e. unlikely that relevant antibiotic starttsime would be 
  # AFTER Sepsis-3 onset given organ dysfunction has already occurred
 # lb <- max(1, time_idx - 48 - round(runif(1, 2, 8)))
 # ub <- min(nrow(temp_dt), time_idx + 4 + round(runif(1, 2, 8)))
 
  # Based on pilot 01 results, let's keep onset right around the middle
  lb <- max(1, time_idx - 48 - round(runif(1, 1, 6)))
  ub <- min(nrow(temp_dt), time_idx + 48 + round(runif(1,1,6)))
  
  # This is hackish, but update the cases file to capture the window for each case
  cases_dt[visit == this_visit, window_lower := lb]
  cases_dt[visit == this_visit, window_upper := ub]
  
   
  # Get range
  temp_dt <- temp_dt[lb:ub][, sub_hour := lb:ub]
  
  # List of variable names to plot in order
  feat_vec <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
                'RESP_RATE', 'LACTATE_RESULT', 'SOFA_RESP_SPO2', 
                'creat_first_obs', 'max_WBC', 'SOFA_CNS_GLASGOW_RESULT')
  normal_ranges <- data.table(measure = feat_vec,
                              norm_lo = c(60,  90,  97, 10, 0, 90, 0.55, 4, 14.8),
                              norm_hi = c(100, 120, 99, 20, 2, 100, 1.15, 11, 15.2))
  
  # Prepare data
  temp_dt <- temp_dt[, c('sub_hour', feat_vec), with = FALSE]
  temp_dt[, hour := sub_hour][, sub_hour := NULL]
  temp_dt_m <- melt(temp_dt, id.vars = 'hour')
  temp_dt_m <- merge(temp_dt_m, normal_ranges, 
                     by.x = 'variable', by.y = 'measure',
                     all.x = TRUE)
  
  # Fix order and labels
  temp_dt_m[, var_fixed := factor(variable, 
                                     levels = feat_vec,
                                     labels = c('Heart rate', 'Systolic blood pressure',
                                     'Temperature (F)', 'Respiratory rate',
                                     'Lactate', 'SpO2(%)',
                                     'Creatinine', 'White blood cell count',
                                     'Glasgow Coma Scale'))]
  
  # Hlper function
  "%nin%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)
  # NB SOFA_RESP_SPO2 is carried forward for each hour
  # Probably only real (and new) when other vital signs are also being observed
  # This is probably true MOST of the time
  # Fix this to avoid signaling too much informative presence:
  temp_dt_m <- temp_dt_m[! (variable == 'SOFA_RESP_SPO2' & 
                              (hour %nin% temp_dt_m[variable %in% 
                                                       c('RESP_RATE', 'SYSTOLIC_BP', 
                                                         'TEMPERATURE')][!is.na(value)]$hour))]
  # This is probably also true for the GCS:
  temp_dt_m <- temp_dt_m[! (variable == 'SOFA_CNS_GLASGOW_RESULT' & 
                              (hour %nin% temp_dt_m[variable %in% 
                                                      c('RESP_RATE', 'SYSTOLIC_BP', 
                                                        'TEMPERATURE')][!is.na(value)]$hour))]
  
  # Determine breaks
  inc <- round((ub - lb) / 8)
  
  # Make plot
  pp <- ggplot(temp_dt_m, aes(hour, value)) +
    theme_bw() + 
    scale_x_continuous('hour', 
                       limits = c(lb, ub),
                       expand = c(.03, .03),
                       breaks = round(c(lb + 0:8 * inc)),
                       labels = round(c(lb + 0:8 * inc - lb + 1))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0.03, 0.03)) +
    geom_rect(aes(xmin = lb, xmax = ub, ymin = norm_lo, ymax = norm_hi),
             fill = 'lightgray', alpha = 0.2) +
    geom_line(data = temp_dt_m[! is.na(value)]) +
    geom_point() + 
    # ^^ NB putting this point layer at the end 
    # makes the points appear on "top" so the spike lines are more easily activated
    lemon::facet_rep_wrap(~ var_fixed, ncol = 1, 
                          scales = 'free_y', repeat.tick.labels = TRUE) +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 12),
          axis.text=element_text(size = 8),
          strip.background =element_rect(fill="white")) 

  ggsave(pp, filename = paste0('images/image', img_idx, '.png'), width = 5, height = 9)

  
  return(ggplotly(pp) %>% 
           layout(xaxis = list(showspikes = TRUE)) %>% 
           config(displayModeBar = FALSE))
} 

# Test it out
images_list <- lapply(1:NUM_CASES, function(i) make_viz(cases_dt[i]$visit, cases_dt[i]$sepsis_onset_dt, i))

# Save images to file
img_path <- 'images/'
setwd(img_path)
invisible(lapply(1:NUM_CASES, function(i) {
  htmlwidgets::saveWidget(images_list[[i]], file = paste0('image', i, '.html'))
}))
setwd('../')


# Save crosswalk data to file - NB. CONTAINS PHI!
cases_dt[, png_filename := paste0(img_path, 'image', .I, '.png')]
cases_dt[, html_filename := paste0(img_path, 'image', .I, '.html')]
fwrite(cases_dt, 'crosswalk_sepsis.csv')

# Also generate 1 control vignette, i.e. those that definitely don't have sepsis
control_visits <- sep_all[! PAT_ENC_CSN %in% sep_events$visits
                          ][hosp_los_hrs > 72
                            ][has_all_fields == TRUE
                              ][, low_severity := max(SOFA_TOTAL_SCORE) < 1, by = PAT_ENC_CSN
                                ][
                                  ][low_severity == TRUE
                                    ][PAT_ENC_CSN == sample(unique(PAT_ENC_CSN), 1)]

make_control_viz <- function(this_visit, onset_time, img_idx) {
  
  temp_dt <- sep_all[PAT_ENC_CSN == this_visit][order(PD_BEG_TIMESTAMP)]
  time_idx <- which(temp_dt$PD_BEG_TIMESTAMP == onset_time)
  
  plot_title <- with(temp_dt, paste0(AGE[1], '-year-old',
                                     # RACE_1[1], GENDER[1], # TOO MUCH BIAS RISK
                                     ' admitted for\n',
                                     ADMIT_DX_DESC[1]))
  
  # If admission diagnosis is very long, truncate it after the first comma
  truncate_title_if_needed <- function(x) {
    if (nchar(x) > 65) {
      first_comma <- regexpr(',', x)
      new_string <- substr(x, 1, first_comma - 1)
      return(new_string)
    } else {
      return(x)
    }
  }
  
  #plot_title <- truncate_title_if_needed(plot_title)
  
  # Randomize offsets so Sep-3 onset is at the RIGHT side
  # This is biologically plausible based on current definition
  # i.e. unlikely that relevant antibiotic starttsime would be 
  # AFTER Sepsis-3 onset given organ dysfunction has already occurred
  # lb <- max(1, time_idx - 48 - round(runif(1, 2, 8)))
  # ub <- min(nrow(temp_dt), time_idx + 4 + round(runif(1, 2, 8)))
  
  # Based on pilot 01 results, let's keep onset right around the middle
  lb <- max(1, time_idx - 48 - round(runif(1, 1, 6)))
  ub <- min(nrow(temp_dt), time_idx + 48 + round(runif(1,1,6)))
  
  # This is hackish, but update the cases file to capture the window for each case
  cases_dt[visit == this_visit, window_lower := lb]
  cases_dt[visit == this_visit, window_upper := ub]
  
  # Get range
  temp_dt <- temp_dt[lb:ub][, sub_hour := lb:ub]
  
  # List of variable names to plot in order
  feat_vec <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
                'RESP_RATE', 'LACTATE_RESULT', 'SOFA_RESP_SPO2', 
                'creat_first_obs', 'max_WBC', 'SOFA_CNS_GLASGOW_RESULT')
  normal_ranges <- data.table(measure = feat_vec,
                              norm_lo = c(60,  90,  97, 10, 0, 90, 0.55, 4, 14.8),
                              norm_hi = c(100, 120, 99, 20, 2, 100, 1.15, 11, 15.2))
  
  # Prepare data
  temp_dt <- temp_dt[, c('sub_hour', feat_vec), with = FALSE]
  temp_dt[, hour := sub_hour][, sub_hour := NULL]
  temp_dt_m <- melt(temp_dt, id.vars = 'hour')
  temp_dt_m <- merge(temp_dt_m, normal_ranges, 
                     by.x = 'variable', by.y = 'measure',
                     all.x = TRUE)
  
  # Fix order and labels
  temp_dt_m[, var_fixed := factor(variable, 
                                  levels = feat_vec,
                                  labels = c('Heart rate', 'Systolic blood pressure',
                                             'Temperature (F)', 'Respiratory rate',
                                             'Lactate', 'SpO2(%)',
                                             'Creatinine', 'White blood cell count',
                                             'Glasgow Coma Scale'))]
  
  # Hlper function
  "%nin%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)
  # NB SOFA_RESP_SPO2 is carried forward for each hour
  # Probably only real (and new) when other vital signs are also being observed
  # This is probably true MOST of the time
  # Fix this to avoid signaling too much informative presence:
  temp_dt_m <- temp_dt_m[! (variable == 'SOFA_RESP_SPO2' & 
                              (hour %nin% temp_dt_m[variable %in% 
                                                      c('RESP_RATE', 'SYSTOLIC_BP', 
                                                        'TEMPERATURE')][!is.na(value)]$hour))]
  # This is probably also true for the GCS:
  temp_dt_m <- temp_dt_m[! (variable == 'SOFA_CNS_GLASGOW_RESULT' & 
                              (hour %nin% temp_dt_m[variable %in% 
                                                      c('RESP_RATE', 'SYSTOLIC_BP', 
                                                        'TEMPERATURE')][!is.na(value)]$hour))]
  
  # Determine breaks
  inc <- round((ub - lb) / 8)
  
  # Make plot
  pp <- ggplot(temp_dt_m, aes(hour, value)) +
    theme_bw() + 
    scale_x_continuous('hour', 
                       limits = c(lb, ub),
                       expand = c(.03, .03),
                       breaks = round(c(lb + 0:8 * inc)),
                       labels = round(c(lb + 0:8 * inc - lb + 1))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0.03, 0.03)) +
    geom_rect(aes(xmin = lb, xmax = ub, ymin = norm_lo, ymax = norm_hi),
              fill = 'lightgray', alpha = 0.2) +
    geom_line(data = temp_dt_m[! is.na(value)]) +
    geom_point() + 
    # ^^ NB putting this point layer at the end 
    # makes the points appear on "top" so the spike lines are more easily activated
    lemon::facet_rep_wrap(~ var_fixed, ncol = 1, 
                          scales = 'free_y', repeat.tick.labels = TRUE) + 
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 12),
          axis.text=element_text(size = 8),
          strip.background =element_rect(fill="white")) 
  
  ggsave(pp, filename = paste0('images/image_control_', img_idx, '.png'), width = 5, height = 9)
  
  
  return(ggplotly(pp) %>% 
           layout(xaxis = list(showspikes = TRUE)) %>% 
           config(displayModeBar = FALSE))
}

make_control_viz(control_visits$PAT_ENC_CSN[1], 
                 control_visits$PD_BEG_TIMESTAMP[round(nrow(control_visits) / 2)],
                 1)


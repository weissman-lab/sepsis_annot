# This file create a visualization summary of a sample of patients 
# during  a period around the time when they might have sepsis
# It will save the de-identified output to an image file
# And generate a crosswalk csv file with PHI

# --------------------------------------------------------------------
SEPSIS_DATA_PATH = "/data/uphs_sepsis/updated_combo_Jan_2020/single_file_augmented_data_Jan_2020.csv"
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
sep_onset_dt[, sepsis_onset_dt := as.POSIXct(onset_char_ts, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, hosp_adm_dt := as.POSIXct(adm_datetime, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, antbx_start_dt := as.POSIXct(antbx_start_char_ts, format = '%Y-%m-%dT%H:%M:%S')]
sep_onset_dt[, hosp_adm_to_sepsis_onset_hrs := as.numeric(sepsis_onset_dt - hosp_adm_dt, 
                                                          units = 'hours')]
sep_onset_dt[, antbx_delay_hrs := as.numeric(antbx_start_dt - sepsis_onset_dt,
                                             units = 'hours')]
sep_onset_dt <- sep_onset_dt[antbx_delay_hrs <= 48]
sep_onset_dt <- sep_onset_dt[hosp_adm_to_sepsis_onset_hrs >= 48]

# Get sampling of cases
# Per Rebecca's suggestion, sample in a stratified way across antibiotic delays
sep_onset_dt <- sep_onset_dt[, antbx_delay_cat := cut(antbx_delay_hrs, breaks = c(-1, 0, 6, 12, 24, 49))]
N_PER_GROUP <- round(NUM_CASES / length(unique(sep_onset_dt$antbx_delay_cat)))
cases_dt <- sep_onset_dt[, .SD[sample(.N, N_PER_GROUP)], by = antbx_delay_cat]
  
# Some data clean up
sep_all[RESPIRATORY_RATE == 'WDL', RESP_RATE := 18]
sep_all[RESPIRATORY_RATE == '', RESP_RATE := NA]
sep_all[! is.na(as.numeric(RESPIRATORY_RATE)),
        RESP_RATE := as.numeric(RESPIRATORY_RATE)]
                                                                       
# Now make a function to generate a plot for each user
make_viz <- function(visit, onset_time, img_idx) {
  
  temp_dt <- sep_all[PAT_ENC_CSN == visit][order(PD_BEG_TIMESTAMP)]
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
                'RESP_RATE', 'LACTATE_RESULT', 'SOFA_RESP_SPO2')
  normal_ranges <- data.table(measure = feat_vec,
                              norm_lo = c(60,  90,  97, 10, 0, 90),
                              norm_hi = c(100, 120, 99, 20, 2, 100))
  
  # Prepare data
  temp_dt <- temp_dt[, c('sub_hour', feat_vec), with = FALSE]
  temp_dt[, hour := sub_hour][, sub_hour := NULL]
  temp_dt_m <- melt(temp_dt, id.vars = 'hour')
  temp_dt_m <- merge(temp_dt_m, normal_ranges, 
                     by.x = 'variable', by.y = 'measure',
                     all.x = TRUE)
  # NB SOFA_RESP_SPO2 is carried forward for each hour
  # Probably only real (and new) when ! is.na(RESP_RATE)
  # Fix this to avoid signaling too much informative presence:
  temp_dt_m <- temp_dt_m[variable == 'SOFA_RESP_SPO2']
  # TODO: FIX THIS!!!!!
  
  # Make plot
  pp <- ggplot(temp_dt_m, aes(hour, value)) +
    theme_bw() + 
    scale_x_continuous('hour', 
                       limits = c(lb, ub),
                       expand = c(.03,.03),
                       breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.02, 0.02)) +
    geom_rect(aes(xmin = lb, xmax = ub, ymin = norm_lo, ymax = norm_hi),
             fill = 'gray', alpha = 0.3) +
    geom_line(data = temp_dt_m[! is.na(value)]) +
    geom_point() + 
    # ^^ NB putting this point layer at the end 
    # makes the points appear on "top" so the spike lines are more easily activated
    facet_wrap(~ variable, ncol = 1, scales = 'free_y') +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 9))

  ggsave(pp, filename = paste0('images/image', img_idx, '.png'), width = 6, height = 8)

  
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


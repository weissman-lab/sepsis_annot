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
library(patchwork)

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
                     'admitted for\n',
                     ADMIT_DX_DESC[1]))

       # Randomize offsets so Sep-3 onset is at the RIGHT side
       # This is biologically plausible based on current definition
       # i.e. unlikely that relevant antibiotic starttsime would be 
       # AFTER Sepsis-3 onset given organ dysfunction has already occurred
       lb <- max(1, time_idx - 48 - round(runif(1, 2, 8)))
       ub <- min(nrow(temp_dt), time_idx + 4 + round(runif(1, 2, 8)))

       # Get range
       temp_dt <- temp_dt[lb:ub][, sub_hour := lb:ub]

       # Plot function
       make_series <- function(dd, val) {
                
                # Make plot
       	        pp <- ggplot(dd[! is.na(dd[[val]])], aes_string('sub_hour', val)) +
	   	        geom_point() + 
       	                geom_line(na.rm = TRUE) +
		        theme_bw() + 
       	                scale_x_continuous('hour', 
       	                                   limits = c(lb, ub),
       	                                   expand = c(.03,.03),
       	                                   breaks = scales::pretty_breaks(n = round((ub - lb) / 10))) +
                        scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
                        annotate('rect',
                                 xmin = lb,
                                     xmax = ub,
                                     ymin = normal_ranges[[val]][1] , 
                                        ymax = normal_ranges[[val]][2], 
                                 color = 'gray', alpha = 0.2)
       	        return(pp)
        }

       # List of variable names to plot in order
       feat_vec <- c('HEART_RATE', 'SYSTOLIC_BP', 'TEMPERATURE',
       		'RESP_RATE', 'LACTATE_RESULT')
       
       normal_ranges <- list(HEART_RATE = c(60, 100),
                             SYSTOLIC_BP = c(100, 140),
                             TEMPERATURE = c(97, 99),
                             RESP_RATE = c(10,20),
                             LACTATE_RESULT = c(0, 2))
       
       # Generate and combine plots
       pl_lst <- lapply(feat_vec, function(val) make_series(temp_dt, val)) 
  
       final_plot <- wrap_elements(Reduce(`+`, pl_lst) + 
        plot_layout(ncol = 1, guides = 'collect')) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(hjust = 0.5))
        
        #return(subplot(pl_lst, nrows = length(feat_vec)))
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


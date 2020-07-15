## Sepsis Annotation Project

This directory contains files related to the sepsis annotation project. The goal is to gather human annotations of "When would have been reasonable to start antibiotics?" and to compare those with other established definitions of sepsis "time zero".

## New version

After multiple revisions to the implementation of the Sepsis-3 definition, sepsis cases are now identified in the file `sepsis_dataset_sepsis3_cohort.dta` that was developed by Jen. Important fields include: 

- sofasepsis_case = sepsis case (community and hospital acquired)
- hosp_sofasepsis_case = hospital acquired sepsis case
- sofa_cirt_first = first criteria to occur out of bcx draw, organ dysfunction, and abx initiation in a sepsis case

- sepsis_abx_datetime = date & time of abx start in a sepsis episode
- sepsis_bcx_datetime = date & time of bcx draw in a sepsis episode
- sepsis_sofa_datetime = date & time of onset of organ dysfunction in a sepsis episode (first delta sofa >=2 within 96h bcx window)

- sepsis_abx_time = time of abx start in a sepsis episode
- sepsis_bcx_time = time of bcx draw in a sepsis episode
- sepsis_sofa_time = time of onset of organ dysfunction in a sepsis episode

- datetime_zero = date and time of first sepsis criteria met in a sepsis case
- time_zero = time of first sepsis criteria met in a sepsis case
- time_to_abx = time btw time zero and abx start

## Old version (NB this approach is now deprecated as of the final version)

To acquire annotations we'll have to pull out sample cases. We'll do this by sampling around sepsis cases identified by the Sepsis-3 criteria. The `sepsis_cases_sepsis3.csv` is a csv version of Rachel's `2019.12.04_sofa18` file.

The definition for the Sepsis-3 onset of sepsis is given as:

-if sepsis was met by BCx first (`sofa_bcx=1`), then `bcx_datetime` is the onset of sepsis

-if sepsis was met by Abx first (`sofa_abx=1`), then `abx_datetime` is the onset of sepsis

-if sepsis was met by SOFA score first (`sofa_sofa=1`), then `sofa_datetime` is the onset of sepsis



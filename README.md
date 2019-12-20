## Sepsis Annotation Project

This directory contains files related to the sepsis annotation project. The goal is to gather human annotations of "When would have been reasonable to start antibiotics?" and to compare those with other established definitions of sepsis "time zero".

To acquire annotations we'll have to pull out sample cases. We'll do this by sampling around sepsis cases identified by the Sepsis-3 criteria. The `sepsis_cases_sepsis3.csv` is a csv version of Rachel's `2019.12.04_sofa18` file.

The definition for the Sepsis-3 onset of sepsis is given as:

-if sepsis was met by BCx first (`sofa_bcx=1`), then `bcx_datetime` is the onset of sepsis

-if sepsis was met by Abx first (`sofa_abx=1`), then `abx_datetime` is the onset of sepsis

-if sepsis was met by SOFA score first (`sofa_sofa=1`), then `sofa_datetime` is the onset of sepsis



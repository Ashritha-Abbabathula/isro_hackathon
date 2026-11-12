# Quantifying the Predictive Contribution of Hard X-ray, Soft X-ray and EUV Observations for Early Solar Flare Forecasting

Evidence from Aditya-L1 (SoLEXS + HEL1OS), GOES-16 EUVS and SDO/HMI, 1 July to 31 December 2024 (183 days, 1,598 flares of class C or above). Author: Ashritha Abbabathula, VIT Chennai. Analysis in R, preprocessing in Python. This grew out of Problem Statement 15 of the Bharatiya Antariksh Hackathon 2026.

## Main findings
- **Nowcasting:** HEL1OS 20-40 keV follows the derivative of the SoLEXS soft X-ray signal in 95% of 372 M/X flares (median rank correlation 0.41; Neupert effect). When both instruments detect a flare, HEL1OS registers it a median of 1 minute earlier.
- **Forecasting:** skill is small (TSS 0.05-0.10 for C+ targets, 0.15-0.25 for M+). Beyond flare history and magnetic state, only soft X-ray adds a little (about 0.04-0.05 TSS), and a single SoLEXS-level feature matches it, so it is largely persistence. Hard X-ray, EUV and full-disk magnetic summaries add nothing consistent. Forward-in-time validation lowers TSS by about 0.02 on average.

## Contents
- `paper/IEEE_paper_Aditya-L1.pdf`: the paper (IEEE conference format). `paper/latex/`: LaTeX source and auto-generated tables.
- `R/`: analysis pipeline (`20_features_v2.R` -> `21_models_v2.R` -> `12`/`13`/`15`/`16` -> `14`/`17`; `lib_models.R` holds the models and metrics).
- `scripts_py/`: preprocessing of the PRADAN light curves and the SDO/HMI magnetic data.
- `notebook/`: `Research_Textbook.html`, a step-by-step walkthrough.
- `outputs_v2/` (leave-one-block-out), `outputs_v2_forward/` (forward-in-time), `outputs_v2_sepincl/` (SEP periods included): result tables and figures.
- `data/`: flare list and SEP periods. The one-minute instrument tables and the 385 MB feature table are not included; rebuild them with the scripts from the raw data (ISRO PRADAN, NOAA NCEI, JSOC).

## Limits
One phase of one solar cycle, 21 X-class flares, full-disk (not per-region) magnetic data only, exploratory statistics. Not an operational forecasting system.

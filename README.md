# Optimising Fast Radio Burst Search Pipelines: A Comparative Study of Digital Image Processing Filters for Structural Morphology Preservation

An offline data processing framework designed to optimize Fast Radio Burst (FRB) detection using public CHIME telescope data. This project combines 2D Digital Image Processing (DIP) with astrophysics to preserve signal morphology, mitigate RFI, and enhance Signal-to-Noise Ratios (SNR) over native 1D search methods.

---

## Project Overview

Fast Radio Bursts (FRBs) are millisecond-duration extragalactic radio transients recorded as 2D time-frequency dynamic spectra (waterfalls). Traditional pipelines often collapse these observations into 1D time series prematurely, risking the loss of fine structural details due to instrumental noise and Radio Frequency Interference (RFI).

This framework introduces a tiered offline strategy:
1. Adaptive Preprocessing: Uses parametric Median Absolute Deviation (MAD) and Interquartile Range (IQR) methods for RFI masking and local-window background detrending.
2. Asymmetric Matched-Filtering: Replaces standard boxcar kernels with a Fast Rise, Exponential Decay (FRED) template and non-linear Exponentially Modified Gaussian (EMG) fitting to eliminate peak-alignment errors and timing slips caused by scattering tails.

---

## Key Results & Findings

* SNR Improvements: Achieved a 1.3x to 1.5x boost in mean SNR (Pipeline A) over native baseline search pipelines.
* Filter Rankings: The Wiener filter provided the best overall compromise between noise reduction and structural morphology preservation.
* Timing Reassignment: Zero-bias timing reassignment using FRED templates effectively eliminated temporal arrival-time bias on scattered bursts.
* Statistical Validation: Validated on 535 authentic FRB events from the CHIME/FRB Open Data Catalog. Kolmogorov-Smirnov (KS) tests (p = 0.35) confirmed no significant Dispersion Measure (DM) distribution bias between the processing subset (200 events) and the full catalog.

---

## Repository Structure

```text
.
├── FRBs ABSTRACT.pdf             # Extended abstract detailing methodology and findings
├── FRB FINAL POSTER.pdf          # Research poster presented at DU/DSKC
├── pipeline_detailed_flowchart.pdf # Full system architecture and pipeline flow diagram
│
├── 200_shortscript.ipynb         # Batch processing script for 200 FRB sample subset
├── ChimeSNRcal.ipynb             # Calculation routines for Peak & Integrated SNR metrics
├── code1.ipynb                   # Native chime offline pipeline (USED from the CHIME website)
├── biastester.py                 # Script executing Kolmogorov-Smirnov (KS) bias testing
│
├── Results(A).png                # Generated plot: Bland Altman plot for the pipelines 
├── Result(B).png                 # Generated plot: RFI leakage comparisons across 2 methods used 
├── Result(C).png                 # Generated plot: SNR gain comparisons across pipelines
├── Biastest_result.png           # Generated plot: DM distribution KS-test validation
│
├── pipeline_comparison_report.csv# Quantitative benchmarking summary across all filters
├── dm_and_snr_audit_report.csv   # Audited DM and SNR metadata for sample catalog
│
├── LICENSE                       # MIT License
└── README.md                     # Project documentation

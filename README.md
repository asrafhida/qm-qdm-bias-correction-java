
/


























Readme · MD
# Bias Correction Scripts for CMIP6 Temperature Data over Java, Indonesia
 
MATLAB implementation of Quantile Mapping (QM) and Quantile Delta Mapping (QDM) used to bias-correct MPI-ESM1-2-HR temperature output against ERA5-Land reanalysis over Java, Indonesia. The scripts support the analysis reported in the accompanying manuscript comparing trend signal preservation between the two methods across lowland and highland zones.
 
## Files
 
**fitter.m**
Fits ten candidate parametric distributions to a data vector and selects the best fit by a chosen criterion (Anderson-Darling statistic, AIC, or Kolmogorov-Smirnov p-value). Domain (real or positive) is inferred from the sign of the input data. Called internally by QM.m and QDM.m to estimate the empirical distributions of observed and modeled data.
 
**QM.m**
Applies static Quantile Mapping. Fits a distribution to the observed series and to the modeled series over the calibration period, then maps modeled values through the modeled CDF and the observed inverse CDF. The transfer function is fixed after calibration and applied unchanged to validation or projection data.
 
**QDM.m**
Applies Quantile Delta Mapping following Cannon et al. (2015). Distributions are refit within a moving 240-month window (±10 years) centered on each time step, preserving the modeled trend by correcting the quantile-mapped value with the model's own change signal (delta) at that time step, rather than relying on a single static transfer function.
 
**ValidasiKoreksiBias.m**
Out-of-sample validation script. Loads the calibration-period transfer functions (fitted on ERA5-Land and MPI-ESM1-2-HR, 1970-1999) from the calibration output file, applies them unchanged to CMIP6 data over 2000-2014, and compares the corrected output against ERA5-Land observations for the same period. This is a differential split-sample test of the stationarity assumption (Teutschbein and Seibert 2012; Maraun 2016). Falls back to an empirical QM or QDM implementation, defined as local functions within the script, whenever the parametric fit from QM.m or QDM.m fails or returns near-constant output for a grid cell.
 
**SSP_emp_par.m**
Projection script. Applies the same 1970-1999 calibration transfer functions to CMIP6 output under SSP2-4.5 and SSP5-8.5, 2015-2100. No observational data are used in this period, since none exist for future scenarios. Computes regional warming trends for the raw model, QM-corrected, and QDM-corrected series, and reports the percentage of the raw trend each method retains. Includes a switchable empirical QDM mode (`cfg.metode_qdm = 'empiris'`) that avoids per-time-step distribution refitting, since the parametric mode can take 9 to 12 hours per scenario at this data volume.
 
## Execution order
 
1. `fitter.m` is not run standalone; it is invoked by `QM.m` and `QDM.m`.
2. A calibration step (not included in this repository) fits `QM.m` and `QDM.m` on ERA5-Land and MPI-ESM1-2-HR data for 1970-1999 and saves the result, including the fitted transfer functions and grid metadata, to a `.mat` file.
3. `ValidasiKoreksiBias.m` reads that calibration file, applies the fitted transfer functions to 2000-2014 CMIP6 data, and evaluates against ERA5-Land observations for the same period.
4. `SSP_emp_par.m` reads the same calibration file and applies the transfer functions to CMIP6 output under SSP2-4.5 and SSP5-8.5 through 2100.
5. Both scripts expect the calibration `.mat` file path set via `cfg.f_cal_mat`, and NetCDF input directories set via `BASE` near the top of each script.
## Requirements
 
- MATLAB (developed and tested on R2023a or later)
- Statistics and Machine Learning Toolbox (required for `fitdist`, `adtest`, `kstest`)
- MATLAB NetCDF support (`ncread`, `ncinfo`, `ncreadatt`; part of base MATLAB)
## Citation
 
If you use these scripts, please cite the associated manuscript:
 
Permadi, F. D., et al. Bias correction of CMIP6 temperature projections over Java, Indonesia: a comparison of Quantile Mapping and Quantile Delta Mapping. *Theoretical and Applied Climatology* (submitted).
 
## License
 
MIT
 










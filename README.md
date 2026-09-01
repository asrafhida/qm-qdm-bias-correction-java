# Bias Correction Scripts for CMIP6 Temperature Data over Java, Indonesia
 
MATLAB implementation of Quantile Mapping (QM) and Quantile Delta Mapping (QDM) used to bias-correct MPI-ESM1-2-HR temperature output against ERA5-Land reanalysis over Java, Indonesia. The scripts support the analysis reported in the accompanying manuscript comparing trend signal preservation between the two methods across lowland and highland zones.
 
## Files
 
**fitter.m**
Fits ten candidate parametric distributions to a data vector and selects the best fit by a chosen criterion (Anderson-Darling statistic, AIC, or Kolmogorov-Smirnov p-value). Domain (real or positive) is inferred from the sign of the input data. Called internally by both QM.m and QDM.m to estimate the empirical distributions of observed and modeled data.
 
**QM.m**
Applies static Quantile Mapping. Fits a distribution to the observed series and to the modeled series over the calibration period, then maps modeled values through the modeled CDF and the observed inverse CDF. The transfer function is fixed after calibration and applied unchanged to projection-period data.
 
**QDM.m**
Applies Quantile Delta Mapping following Cannon et al. (2015). Distributions are refit within a moving 240-month window (±10 years) centered on each time step, preserving the modeled trend by correcting the quantile-mapped value with the model's own change signal (delta) at that time step, rather than relying on a single static transfer function.
 
**SKRIPSIFIORENZADP.m**
Full analysis pipeline. Loads ERA5-Land and MPI-ESM1-2-HR NetCDF data, applies land masking, runs QM and QDM grid cell by grid cell over the historical, validation, and SSP2-4.5/SSP5-8.5 projection periods, and produces the trend, spatial noise, and topographic stratification results reported in the manuscript.
 
## Execution order
 
1. `fitter.m` and the distribution-fitting logic are not run standalone; they are invoked by `QM.m` and `QDM.m`.
2. `QM.m` and `QDM.m` are called cell by cell inside `SKRIPSIFIORENZADP.m` for each time period (historical, validation, SSP2-4.5, SSP5-8.5).
3. `SKRIPSIFIORENZADP.m` is the entry point. It expects ERA5-Land and MPI-ESM1-2-HR NetCDF files and GADM v4.0 province boundaries as input; paths are set near the top of the script via the `BASE` variable.
## Requirements
 
- MATLAB (developed and tested on R2023a or later)
- Statistics and Machine Learning Toolbox (required for `fitdist`, `adtest`, `kstest`)
## Citation
 
If you use these scripts, please cite the associated manuscript:
 
Permadi, F. D., et al. Bias correction of CMIP6 temperature projections over Java, Indonesia: a comparison of Quantile Mapping and Quantile Delta Mapping. *Theoretical and Applied Climatology* (submitted).
 
## License
 
MIT

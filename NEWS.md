# marpoly 0.1.1

- Fixed package installation on R 4.x by calling `base::polyroot()` and removing the invalid `stats::polyroot` import.

# marpoly 0.1.0

- Initial GitHub-ready package.
- Added MAR(r,s) simulation with stationarity checks.
- Added periodograms and spectral densities of orders 2, 3, and 4.
- Added pseudo-causal AR order selection and regression-style summaries.
- Added root allocation enumeration, polyspectral estimation, allocation ranking, error extraction, and residual diagnostics.
- Added examples for MAR(1,1), MAR(1,2), skew Student t and Student t innovations, and order combinations 2/3/4.

# marpoly 0.2.0

- Replaced the previous estimation engine with the MAR spectral/bispectral/trispectral minimum-distance logic from the working research code.
- Added canonical frequency-domain selection for the spectrum, bispectrum, and trispectrum; preliminary pseudo-causal AR normalization; root-allocation starts; random starts; and finite-difference sandwich standard errors.
- Added exported workflow functions `MAR_spectral_estimate()`, `MAR_select_all_by_second_order()`, `MAR_select_by_order_sets()`, `simulate_MAR()`, and `rskt()`.
- Updated `marpoly_fit()`, `marpoly_rank()`, `marpoly_errors()`, and plotting helpers to use the revised estimation objects.
- Added p-values to the parameter table and aligned the DFT convention used by periodograms with the paper convention.

# marpoly 0.1.1

- Fixed package installation on R 4.x by calling `base::polyroot()` and removing the invalid `stats::polyroot` import.

# marpoly 0.1.0

- Initial GitHub-ready package.
- Added MAR(r,s) simulation with stationarity checks.
- Added periodograms and spectral densities of orders 2, 3, and 4.
- Added pseudo-causal AR order selection and regression-style summaries.
- Added root allocation enumeration, polyspectral estimation, allocation ranking, error extraction, and residual diagnostics.
- Added examples for MAR(1,1), MAR(1,2), skew Student t and Student t innovations, and order combinations 2/3/4.

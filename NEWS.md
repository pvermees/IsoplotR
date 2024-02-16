# IsoplotR 6.1

Adds four new U-Pb formats for decoupled 206Pb/238U and 207Pb/235U dating. Makes the Th/U ratio optional for U-Pb formats 8, 11 and 12. No longer applies a common Pb correction to discordia fits.

# IsoplotR 6.0

Adds anchored regression to all isochron functions, as well as model-3 fits whereby overdispersion can be attributed to either diachronous isotopic closure or the composition of the inherited component. Also adds an option to anchor U-Pb isochrons to the Stacey-Kramers mantle evolution model. Adds four new U-Pb input formats. Finally, flips the Pb208/Pb207- and Pb208/Pb206- ratios around in the settings, which is a breaking change.

# IsoplotR 5.6

November 10, 2023: Improved initialisation of Ludwig regression for U-Pb formats 1, 2 and 3.

# IsoplotR 5.5

September 29, 2023: Redefines the Aitchison and concordia distances in terms of isometric logratios. Modifies the default cutoffs accordingly.

# IsoplotR 5.4

August 28, 2023: Renamed U-Th-He 'central' ages in helioplots as 'barycentric' ages to avoid confusion with the central ages reported in radial plots. Ads OGLS regression.

# IsoplotR 5.3

May 19, 2023: Improved model-3 regression for York and Ludwig algorithms.

# IsoplotR 5.2 

Feb 21, 2023: Adds Wasserstein distance to MDS plots; propagates disequilibrium activity ratio errors for U-Pb isochrons, as well as common-Pb and concordia age anchor errors. Replaces ordinary with total least squares in model-2 regression. Recasts model-3 regression terms of ages rather than y-intercepts. Changes default output errors back from 2-sigma to Studentised confidence intervals.

# IsoplotR 5.1 

Oct 18, 2022: Fixes (mostly minor) issues with the major upgrade to 5.0. 

# IsoplotR 5.0 

Sept 13, 2022: Errors are now reported as absolute and relative uncertainties at 1σ, 2σ and arbitrary confidence levels. Improved security by the use of shinylight under the hood of IsoplotRgui. Replaced .png with .svg graphics. Fixes small but nasty bug in the minimum age algorithm.

# IsoplotR 4.4 

May 21, 2022: Additional options for U-series evolution diagrams and isochron plots, including initial 230Th-corrections based on the magmatic Th/U activity ratio.

# IsoplotR 4.3 

Mar 3, 2022: Fixed issue with detrital Stacey-Kramers correction, removed gradient function from Ludwig regression.

# IsoplotR 4.2 

Oct 20, 2021: Added 2D U-Pb isochron regression.

# IsoplotR 4.1 

Jul 9, 2021: Removed student t-multiplier from confidence intervals, improved corrections for secular disequilibrium and common Pb.

# IsoplotR 4.0 

Jun 6, 2021: Expanded API; more robust handling of measured initial U-Pb disequilibrium; lots of minor bug fixes.

# IsoplotR 3.8 

Apr 5, 2021: Calculates intercepts of Pb-Pb isochrons with Stacey-Kramers mantle evolution curve; adds Th-Pb dataset to examples; fixes colour scale of radial plots.

# IsoplotR 3.7 

Mar 15, 2021: Changed overdispersion reporting threshold for model-1 regression from MSWD-based to p-value based. Simplified API for anchored discordia regression.

# IsoplotR 3.6 

Dec 17, 2020: Formally added rug plot option to KDEs. Upgraded mirrors to shinyless version, increased robustness of disequilibrium correction, updated whitelist of allowable functions.

# IsoplotR 3.5 

Nov 10, 2020: Minor bug fixes. Compatible with `shinyless' IsoplotRgui

# IsoplotR 3.4 

Jul 8, 2020: Adds options to change the stroke colour of error ellipses (API change); introduces a log-transform to the random effects model for weighted means; and replaces the three parameter minimum age algorithm with a four parameter alternative.

# IsoplotR 3.3 

Mar 18, 2020: Adds Th-Pb geochronology to the toolbox.

# IsoplotR 3.2 

Jan 22, 2020: Fixes several minor bugs and increases stability.

# IsoplotR 3.1 

Oct 8, 2019: Unified Ludwig regression, including a new model-3 option that redefines the overdispersion parameter in terms of time rather than common Pb composition. Various minor bug fixes.

# IsoplotR 3.0 

Aug 5, 2019: Adds isochrons for U-(Th-)Pb data, implements U-Pb initial disequilibrium corrections using matrix exponentials, completely rewritten Ludwig regression functions.

# IsoplotR 2.7

Jun 6, 2019: Adds new input format for U-Pb data that include 208Pb and 232Th; option to perform discordance filter before common Pb correction. This update changes the .json format, which may cause problems for data saved in previous IsoplotR versions.

# IsoplotR 2.6

May 4, 2019: Adds inverse isochron option to K-Ca, Rb-Sr, Sm-Nd, Lu-Hf and Re-Os data. Improved dispersion estimates in model-3 regression.

# IsoplotR 2.5 

Apr 3, 2019: Several improvements, especially for the common-Pb corrections, which improve IsoplotR's robustness to very discordant data. Prepared in the run-up to the detrital zircon U-Pb dating workshop at the Beijing SHRIMP centre.

# IsoplotR 2.4 

Mar 4, 2019: Implements initial disequilibrium corrections for U-Pb data, improved common Pb-corrections and initial ('common') non-radiogenic isotope compositions for other chronometers; increases robustness of the mixture modelling functions.

# IsoplotR 2.3 

Dec 8, 2018: Allows input uncertainties to be specified as absolute and relative errors at 1 or 2&sigma;.

# IsoplotR 2.2 

Nov 9, 2018: Adds ability to hide or omit samples from plots and calculations.

# IsoplotR 2.1 

Oct 20, 2018: Rewrote discordia regression function, which now uses the same maximum likelihood approach for 2D and 3D regression. Added option to anchor the discordia line to the concordia line, or to a specificy common Pb composition.

# IsoplotR 2.0 

Sep 10, 2018: Adds K-Ca geochronology to the toolbox.

# IsoplotR 1.4 

Aug 13, 2018: Adds options for detrital 230Th-corrections from assumed or measured detritus; more flexible weighted mean plots.

# IsoplotR 1.4 

Jun 4, 2018: Includes bug fixes relating to the 'York fit' algorithm and error formatting, and added functionality for Stacey-Kramers evolution curves and ordinary weighted mean calculation.

# IsoplotR 1.2 

May 1, 2018: Includes bug fixes relating to the MSWD of isochrons, and discordia lines in Tera-Wasserburg space. Improved flexibility for loading data from R matrices and data frames.

# IsoplotR 1.1 

Mar 28, 2018: Adds profile log-likelihood confidence intervals for dispersion parameters, confidence envelopes for regression lines. 

# IsoplotR 1.0 

Nov 11, 2017: Adds model-2 and model-3 regression for all data, common-Pb correction for Pb-Pb data.

# IsoplotR 0.18 

Sep 25, 2017: Adds confidence intervals.

# IsoplotR 0.17 

Sep 1, 2017: Added helium isochrons, volcanic Th-U isochrons, better and more flexible scaling of concordia axes.

# IsoplotR 0.16 

Jul 21, 2017: Adds total Pb-U isochrons and common-Pb corrections.

# IsoplotR 0.15 

Jul 11, 2017: Added U-series disequilibrium dating.

# IsoplotR 0.14 

Jun 18, 2017: Adds Pb-Pb dating functionality.

# IsoplotR 0.13 

Jun 25, 2017: Added Lu-Hf functionality.

# IsoplotR 0.12 

May 6, 2017: Various bug fixes, improved documentation.

# IsoplotR 0.11 

Apr 17, 2017: Additional input formatting options added for Ar-Ar, Rb-Sr, Sm-Nd and Re-Os chronometers.

# IsoplotR 0.10 

Jan 8, 2017: Added Rb-Sr and Sm-Nd functionality.

# IsoplotR 0.9 

Dec 20, 2016: Adds Re-Os isochrons. Improved error catching.

# IsoplotR 0.8 

Nov 16, 2016: Adds functionality for fission track data using the external detector method and LA-ICP-MS. Plots these data and all other input formats on radial plots. Deconvolves age distributions into finite mixtures and constraints minimum ages from heteroscedastic data.

# IsoplotR 0.7 

Aug 13, 2016: Added functionality for U-Th-He data and MDS analysis of detrital age distributions.

# IsoplotR 0.6 

Jul 27, 2016: Adds functions to plot age spectra and calculate weighted means.

# IsoplotR 0.5 

Jul 3, 2016: Calculates 40Ar/39Ar ages and isochrons.

# IsoplotR 0.4 

Jun 14, 2016: Added functions to plot KDEs and CADs functions. Renamed and simplified existing functions.

# IsoplotR 0.3 

May 23, 2016: Adds functionality to compute concordia and discordia ages, and perform York-style linear regression of data with correlated errors.

# IsoplotR 0.2 

May 5, 2016: Adds uncertainty associated with decay constants and the 238U/235U ratio to the concordia diagrams.

# IsoplotR 0.1 

Apr 21, 2016: Plots U-Pb data on Wetherill and Tera-Wasserburg concordia diagrams, taking into account error correlations.
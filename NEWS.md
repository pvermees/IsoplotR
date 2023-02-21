# IsoplotR

* Feb 21, 2023 (5.2): Adds Wasserstein distance to MDS plots; propagates disequilibrium activity ratio errors for U-Pb isochrons, as well as common-Pb and concordia age anchor errors. Replaces ordinary with total least squares in model-2 regression. Recasts model-3 regression terms of ages rather than y-intercepts. Changes default output errors back from 2-sigma to Studentised confidence intervals.

* Oct 18, 2022 (5.1): Fixes (mostly minor) issues with the major upgrade to 5.0. 

* Sept 13, 2022 (5.0): Errors are now reported as absolute and relative uncertainties at 1σ, 2σ and arbitrary confidence levels. Improved security by the use of shinylight under the hood of IsoplotRgui. Replaced .png with .svg graphics. Fixes small but nasty bug in the minimum age algorithm.

* May 21, 2022 (4.4): Additional options for U-series evolution diagrams and isochron plots, including initial 230Th-corrections based on the magmatic Th/U activity ratio.

* Mar 3, 2022 (4.3): Fixed issue with detrital Stacey-Kramers correction, removed gradient function from Ludwig regression.

* Oct 20, 2021 (4.2): Added 2D U-Pb isochron regression.

* Jul 9, 2021 (4.1): Removed student t-multiplier from confidence intervals, improved corrections for secular disequilibrium and common Pb.

* Jun 6, 2021 (4.0): Expanded API; more robust handling of measured initial U-Pb disequilibrium; lots of minor bug fixes.

* Apr 5, 2021 (3.8): Calculates intercepts of Pb-Pb isochrons with Stacey-Kramers mantle evolution curve; adds Th-Pb dataset to examples; fixes colour scale of radial plots.

* Mar 15, 2021 (3.7): Changed overdispersion reporting threshold for model-1 regression from MSWD-based to p-value based. Simplified API for anchored discordia regression.

* Dec 17, 2020 (3.6): Formally added rug plot option to KDEs. Upgraded mirrors to shinyless version, increased robustness of disequilibrium correction, updated whitelist of allowable functions.

* Nov 10, 2020 (3.5): Minor bug fixes. Compatible with `shinyless' IsoplotRgui

* Jul 8, 2020 (3.4): Adds options to change the stroke colour of error ellipses (API change); introduces a log-transform to the random effects model for weighted means; and replaces the three parameter minimum age algorithm with a four parameter alternative.

* Mar 18, 2020 (3.3): Adds Th-Pb geochronology to the toolbox.

* Jan 25, 2020: Release of a Chinese language version of the graphical user interface. Thanks to Shuan Yan (Guangzhou), Pei Zhang and Hao Fang (Beijing), Qiuye Yu (Changsha), Hongcheng Guo (Lehigh), and Shihu Li (Lancaster) for the translation!

* Jan 22, 2020 (3.2): Fixes several minor bugs and increases stability.

* Oct 8, 2019 (3.1): Unified Ludwig regression, including a new model-3 option that redefines the overdispersion parameter in terms of time rather than common Pb composition. Various minor bug fixes.

* Aug 5, 2019 (3.0): Adds isochrons for U-(Th-)Pb data, implements U-Pb initial disequilibrium corrections using matrix exponentials, completely rewritten Ludwig regression functions.

* Jun 6, 2019 (2.7): Adds new input format for U-Pb data that include 208Pb and 232Th; option to perform discordance filter before common Pb correction. This update changes the .json format, which may cause problems for data saved in previous IsoplotR versions.

* May 4, 2019 (2.6): Adds inverse isochron option to K-Ca, Rb-Sr, Sm-Nd, Lu-Hf and Re-Os data. Improved dispersion estimates in model-3 regression.

* Apr 3, 2019 (2.5): Several improvements, especially for the common-Pb corrections, which improve IsoplotR's robustness to very discordant data. Prepared in the run-up to the detrital zircon U-Pb dating workshop at the Beijing SHRIMP centre.

* Mar 4, 2019 (2.4): Implements initial disequilibrium corrections for U-Pb data, improved common Pb-corrections and initial ('common') non-radiogenic isotope compositions for other chronometers; increases robustness of the mixture modelling functions.

* Dec 8, 2018 (2.3): Allows input uncertainties to be specified as absolute and relative errors at 1 or 2&sigma;.

* Nov 9, 2018 (2.2): Adds ability to hide or omit samples from plots and calculations.

* Oct 20, 2018 (2.1): Rewrote discordia regression function, which now uses the same maximum likelihood approach for 2D and 3D regression. Added option to anchor the discordia line to the concordia line, or to a specificy common Pb composition.

* Sep 10, 2018 (2.0): Adds K-Ca geochronology to the toolbox.

* Aug 13, 2018 (1.4): Adds options for detrital 230Th-corrections from assumed or measured detritus; more flexible weighted mean plots.

* Jun 4, 2018 (1.3): Includes bug fixes relating to the 'York fit' algorithm and error formatting, and added functionality for Stacey-Kramers evolution curves and ordinary weighted mean calculation.

* May 1, 2018 (1.2): Includes bug fixes relating to the MSWD of isochrons, and discordia lines in Tera-Wasserburg space. Improved flexibility for loading data from R matrices and data frames.

* Mar 28, 2018 (1.1): Adds profile log-likelihood confidence intervals for dispersion parameters, confidence envelopes for regression lines. 

* Nov 11, 2017 (1.0): Adds model-2 and model-3 regression for all data, common-Pb correction for Pb-Pb data.

* Sep 25, 2017 (0.18): Adds confidence intervals.

* Sep 1, 2017 (0.17): Added helium isochrons, volcanic Th-U isochrons, better and more flexible scaling of concordia axes.

* Jul 21, 2017 (0.16): Adds total Pb-U isochrons and common-Pb corrections.

* Jul 11, 2017 (0.15): Added U-series disequilibrium dating.

* Jun 18, 2017 (0.14): Adds Pb-Pb dating functionality.

* Jun 25, 2017 (0.13): Added Lu-Hf functionality.

* May 6, 2017 (0.12): Various bug fixes, improved documentation.

* Apr 17, 2017 (0.11): Additional input formatting options added for Ar-Ar, Rb-Sr, Sm-Nd and Re-Os chronometers.

* Jan 8, 2017 (0.10): Added Rb-Sr and Sm-Nd functionality.

* Dec 20, 2016 (0.9): Adds Re-Os isochrons. Improved error catching.

* Nov 16, 2016 (0.8): Adds functionality for fission track data using the external detector method and LA-ICP-MS. Plots these data and all other input formats on radial plots. Deconvolves age distributions into finite mixtures and constraints minimum ages from heteroscedastic data.

* Aug 13, 2016 (0.7): Added functionality for U-Th-He data and MDS analysis of detrital age distributions.

* Jul 27, 2016 (0.6): Adds functions to plot age spectra and calculate weighted means.

* Jul 3, 2016 (0.5): Calculates 40Ar/39Ar ages and isochrons.

* Jun 14, 2016 (0.4): Added functions to plot KDEs and CADs functions. Renamed and simplified existing functions.

* May 23, 2016 (0.3): Adds functionality to compute concordia and discordia ages, and perform York-style linear regression of data with correlated errors.

* May 5, 2016 (0.2): Adds uncertainty associated with decay constants and the 238U/235U ratio to the concordia diagrams.

* Apr 21, 2016 (0.1): Plots U-Pb data on Wetherill and Tera-Wasserburg concordia diagrams, taking into account error correlations.
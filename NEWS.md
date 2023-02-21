# IsoplotR

Feb 21, 2023: (v.5.2) Adds Wasserstein distance to MDS plots; propagates disequilibrium activity ratio errors for U-Pb isochrons, as well as common-Pb and concordia age anchor errors. Replaces ordinary with total least squares in model-2 regression. Recasts model-3 regression terms of ages rather than y-intercepts. Changes default output errors back from 2-sigma to Studentised confidence intervals.

Sept 13, 2022: (v.5.0) Errors are now reported as absolute and relative uncertainties at 1σ, 2σ and arbitrary confidence levels. Improved security by the use of shinylight under the hood of IsoplotRgui. Replaced .png with .svg graphics. Fixes small but nasty bug in the minimum age algorithm.

May 21, 2022: (v.4.4) Additional options for U-series evolution diagrams and isochron plots, including initial 230Th-corrections based on the magmatic Th/U activity ratio.

Mar 3, 2022: (v.4.3) Fixed issue with detrital Stacey-Kramers correction, removed gradient function from Ludwig regression.

Oct 20, 2021: (v.4.2) Added 2D U-Pb isochron regression.

Jul 9, 2021: (v.4.1) Removed student t-multiplier from confidence intervals, improved corrections for secular disequilibrium and common Pb.

Jun 6, 2021: (v.4.0) Expanded API; more robust handling of measured initial U-Pb disequilibrium; lots of minor bug fixes.

Apr 5, 2021: (v.3.8) Calculates intercepts of Pb-Pb isochrons with Stacey-Kramers mantle evolution curve; adds Th-Pb dataset to examples; fixes colour scale of radial plots.

Apr 1, 2021: Dengfeng Li hosts a new mirror at Sun Yat-Sen University in Guangzhou (China): `https://isoplotr.sysu.edu.cn/isoplotr/`

Mar 15, 2021: (v.3.7) Changed overdispersion reporting threshold for model-1 regression from MSWD-based to p-value based. Simplified API for anchored discordia regression.

Dec 17, 2020: (v.3.6) Formally added rug plot option to KDEs. Upgraded mirrors to shinyless version, increased robustness of disequilibrium correction, updated whitelist of allowable functions.

Nov 10, 2020: (v.3.5) Minor bug fixes. Compatible with `shinyless' IsoplotRgui

Sep 10, 2020: Danny Stockli hosts a new mirror at UT Austin: `http://isoplotr.geo.utexas.edu`

Jul 8, 2020 (v.3.4): Adds options to change the stroke colour of error ellipses (API change); introduces a log-transform to the random effects model for weighted means; and replaces the three parameter minimum age algorithm with a four parameter alternative.

May 4, 2020: LA.TE.ANDES host a new IsoplotR mirror in Salta, Argentina: `http://isoplotr.lateandes.com/`

Mar 18, 2020 (v.3.3): Adds Th-Pb geochronology to the toolbox.

Jan 25, 2020: Release of a Chinese language version of the graphical user interface. Thanks to Shuan Yan (Guangzhou), Pei Zhang and Hao Fang (Beijing), Qiuye Yu (Changsha), Hongcheng Guo (Lehigh), and Shihu Li (Lancaster) for the translation!

Jan 22, 2020 (v.3.2): Fixes several minor bugs and increases stability.

Oct 30, 2019: Yang Li and Liguang Wu host a new IsoplotR mirror at the Chinese Academy of Sciences in Beijing: `http://www.isoplotr.com/isoplotr/`

Oct 8, 2019 (v.3.1): Unified Ludwig regression, including a new model-3 option that redefines the overdispersion parameter in terms of time rather than common Pb composition. Various minor bug fixes.

Aug 5, 2019 (v.3.0): Adds isochrons for U-(Th-)Pb data, implements U-Pb initial disequilibrium corrections using matrix exponentials, completely rewritten Ludwig regression functions.

Jun 6, 2019 (v2.7): Adds new input format for U-Pb data that include 208Pb and 232Th; option to perform discordance filter before common Pb correction. This update changes the .json format, which may cause problems for data saved in previous IsoplotR versions.

May 24, 2019: Addition of an American IsoplotR mirror at UC Santa Barbara:
`http://isoplotr.geol.ucsb.edu/isoplotr/`

May 4, 2019 (v2.6): Adds inverse isochron option to K-Ca, Rb-Sr, Sm-Nd, Lu-Hf and Re-Os data. Improved dispersion estimates in model-3 regression.

Apr 11, 2019: Two new mirror IsoplotR sites come online:
UK (BGS): `https://shiny.bgs.ac.uk/IsoplotR/`
China (Qiuye Yu): `http://chinageology.org:8080/`

Apr 3, 2019 (v2.5): Several improvements, especially for the common-Pb corrections, which improve IsoplotR's robustness to very discordant data. Prepared in the run-up to the detrital zircon U-Pb dating workshop at the Beijing SHRIMP centre.

Mar 4, 2019 (v2.4): Implements initial disequilibrium corrections for U-Pb data, improved common Pb-corrections and initial ('common') non-radiogenic isotope compositions for other chronometers; increases robustness of the mixture modelling functions.

Dec 8, 2018 (v2.3): Allows input uncertainties to be specified as absolute and relative errors at 1 or 2σ.

Nov 9, 2018 (v2.2): Adds ability to hide or omit samples from plots and calculations.

Oct 20, 2018 (v2.1): Rewrote discordia regression function, which now uses the same maximum likelihood approach for 2D and 3D regression. Added option to anchor the discordia line to the concordia line, or to a specificy common Pb composition.

Sep 10, 2018 (v2.0): Adds K-Ca geochronology to the toolbox.

Aug 13, 2018 (v1.4): Adds options for detrital 230Th-corrections from assumed or measured detritus; more flexible weighted mean plots.

Jun 4, 2018 (v1.3): Includes bug fixes relating to the 'York fit' algorithm and error formatting, and added functionality for Stacey-Kramers evolution curves and ordinary weighted mean calculation.

May 1, 2018 (v1.2): Includes bug fixes relating to the MSWD of isochrons, and discordia lines in Tera-Wasserburg space. Improved flexibility for loading data from R matrices and data frames.

Apr 11, 2018: IsoplotR formally introduced to the Geoscience community at EGU (Vienna).

Apr 3, 2018: IsoplotR paper accepted for publication in Geoscience Frontiers (a free and open journal for free and open software).
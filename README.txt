

scripts to use:


script_displayPhaseLockPercVsdbsTrialCond_barplot.m
Shows relative population of neuron TRIAL-OBSERVATIONS % CLASSIFIED as phase-locking vs non.


script_displayPopRates_PREandDBS_barplot.m
Shows Comparison of spike rates pre-DBS vs during-DBS, arranged by dbsFrequency condition.


script_displayRateChangeCategory_barplotStack.m
Shows Categorizations of DBS-induced rate-changes according to "excite", "inhib", or "no change", and the % of neuron population in each category, according to DBS condition


script_displayPhaseLockEntropyChangeVsdbsCond_boxplots.m
Shows DirectISI entropy deltaH's grouped by DBS condition in boxplots, and can filter the results based on whether a neuron phase-locked or not.
Also performs KruskalWallis test, and post-hoc t-test/ranksum on LFS, HFS


script_CompareRatePsthIsi.m
Shows example activity around DBS, in terms of Rate, PSTH, and log-ISI. Cycles thru each record from a table one at a time. Cancel/Exit with "ctrl + c" in Matlab Command Window.
Also lists unique file's ID (objectID)


script_CompareRatePsthIsi_specificIteration.m
Does exactly what "script_CompareRatePsthIsi.m" does, but at the user-specified file in <specID>.


script_displayPhaseLockNeuronPerc_pieChart.m
Shows pie-chart of number of NEURONS classified as phase-locked vs NOT phase locked
(if NEURON in question has at least one TRIAL-OBSERVATION that showed phase-locking, it gets categorized as phaselocked)


script_displayPercPop_isiAndPhaseLckCategory_pieChart
Shows pie chart of relative proportions of Entropy responses across neuron TRIAL-OBSERVATIONS classifed into 6 categories: nochange, PhsLck_isiHincr, PhsLck_isiHdecr, PhsLckONLY, noPhsLck_isiHincr, noPhsLck_isiHdecr
Also outputs a table showing the relative counts and %'s of TRIAL OBSERVATIONS according to category.


script_displayRateVsEntropy_scatter.m
Shows Rate Vs entropy (and various versions of this), as well as giving Pearson Corr. coefficient and p value of hypothesis test (matlab function: corr)


script_plot3uniqueNeuronLogISIs_freqSwp_PREandDBSon.m
Cycles thru showing log-ISI histograms across frequency conditions

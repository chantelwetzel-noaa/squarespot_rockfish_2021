
# Squarespot Rockfish California Model Runs
library(r4ss)

wd = "N:/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/models"
wd = "C:/Assessments/2021/squarespot_rockfish_2021/models"

# This run incorporates - updated catches with constant rec hist. dicard rate,
# hkl index, and recreational lengths based on LNGTH in the MRFSS period.
model = "1.0_data_catches_index_rec_len"
base = SS_output(file.path(wd, model))
SS_plots(base)

# Add new length estimated from Jason's IGOR analysis for males and females:
# start with the fixed t0
# Male k set equal to the females given the poor sex specific estimate
# Change max age in the model - 34 based on 95% quantile of aged fish
model = "1.1_data_bio_new_len_ests"
base = SS_output(file.path(wd, model))
SS_plots(base) 

# Update the M value based on new max age of 34 (5.4/34) = 0.159
model = "1.2_data_bio_len_max_age"
base = SS_output(file.path(wd, model))
SS_plots(base) 

# Apply data weights with rec devs turned off
model = "1.3_dw_francis"
base = SS_output(file.path(wd, model))
SS_plots(base) 

# Turn on rec dev estimates - removed early devs due to wild behavior
# Suggested sigmaR = 1.6 (although after increasing to this value the new value 
# recommended is now 1.98)
# Really poor gradient on R0 and selectivity parameters
# Population shoots up with large recent recruitments
model = "2.1_rec_devs"
base = SS_output(file.path(wd, model))
SS_plots(base) 

# Apply data weighting to try to tamp down some of the wild behavior of model 2.1
# Also, lambda = 0 for the NWFSC WCGBTS for now
model = "2.2_dw_francis"
base = SS_output(file.path(wd, model))
SS_plots(base) 
# After dw the recommended sigma R decreased to 1.38

# Got the hessian to run on this limited model
model = "2.3_hessian"
base = SS_output(file.path(wd, model))
SS_plots(base) 
# The stock goes low for much to the timeseries with it only shooting up in recent yrs.

# Only allow rec devs starting in 2010
model = "2.4_rec_devs_2010"
base = SS_output(file.path(wd, model))
SS_plots(base) 


#============================================================
# Sensitivities
#============================================================
base = SS_output(file.path(wd, "2.5_rec_devs_sigmaR60"))
catch = SS_output(file.path(wd, "_sensitivities", "2.5_recdevs_sigmaR60_catcheshalf"))
index = SS_output(file.path(wd, "_sensitivities", "2.5_recdevs_sigmaR60_minus2index"))
comps = SS_output(file.path(wd, "_sensitivities", "2.5_recdevs_sigmaR60_minus2comps"))
hkl = SS_output(file.path(wd, "_sensitivities", "2.5_recdevs_sigmaR60_rm_hkl"))

modelnames <- c("Base", "Reduce Recent Catch", "Remove 2yrs Index", 
				"Remove 2yrs Comps", "No HKL")
mysummary <- SSsummarize(list(base, catch, index, comps, hkl))
SSplotComparisons(mysummary, 
			      filenameprefix = "2.5_rec_devs_sigmaR60_",
				  legendlabels = modelnames, 
				  plotdir = file.path(wd, "_sensitivities"),
				  pdf = TRUE)

# HKL Sensitivities
# The single or multiple year runs only removed the lengths but kept index
base = SS_output(file.path(wd, "2.5_rec_devs_sigmaR60"))
hkl1819 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus20182019HKLcomps"))
hkl1719 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus20172019HKLcomps"))
hkl1718 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus20172018HKLcomps"))
hkl2017 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus2017HKLcomps"))
hkl2018 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus2018HKLcomps"))
hkl2019 = SS_output(file.path(wd, "_sensitivities", "2.5_rec_devs_sigR0.6_minus2019HKLcomps"))
hkl = SS_output(file.path(wd, "_sensitivities", "2.5_recdevs_sigmaR60_rm_hkl"))

modelnames <- c("Base",  "-2017 HKL Lengths", "-2018 HKL Lengths", "-2019 HKL Lengths",
				 "-2017 & 2019 HKL Lengths", "-2017 & 2018 HKL Lengths", "-2018 & 2019 HKL Lengths", "No HKL Lengths & Index")
mysummary <- SSsummarize(list(base, hkl2017, hkl2018, hkl2019, hkl1719, hkl1718, hkl1819,   hkl))
SSplotComparisons(mysummary, 
				  legendloc = 'bottomleft',
			      filenameprefix = "2.5_rec_devs_sigmaR60_hkl_",
				  legendlabels = modelnames, 
				  plotdir = file.path(wd, "_sensitivities"),
				  pdf = TRUE)


#============================================================
# Retrospecitves
#============================================================
base = SS_output(file.path(wd, "2.5_rec_devs_sigmaR60"))
retro1 = SS_output(file.path(wd, "_retro", "2.5_rec_devs_sigmaR60", "retro-1"))
retro2 = SS_output(file.path(wd, "_retro", "2.5_rec_devs_sigmaR60", "retro-2"))
retro3 = SS_output(file.path(wd, "_retro", "2.5_rec_devs_sigmaR60", "retro-3"))
retro4 = SS_output(file.path(wd, "_retro", "2.5_rec_devs_sigmaR60", "retro-4"))
retro5 = SS_output(file.path(wd, "_retro", "2.5_rec_devs_sigmaR60", "retro-5"))
retro10 = SS_output(file.path(wd, "_retro","2.5_rec_devs_sigmaR60", "retro-10"))

modelnames <- c("Base", "-1", "-2", "-3", "4", "-5", "-10")
mysummary <- SSsummarize(list(base, retro1, retro2, retro3, retro4, retro5, retro10))
SSplotComparisons(mysummary, 
			      filenameprefix = "2.5_rec_devs_sigmaR60_",
				  endyrvec = c(rev(2016:2021), 2011),
				  legendlabels = modelnames, 
				  plotdir = file.path(wd, "_retro"),
				  pdf = TRUE)
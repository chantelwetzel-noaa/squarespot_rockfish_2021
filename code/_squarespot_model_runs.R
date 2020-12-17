
# Squarespot Rockfish California Model Runs

library(r4ss)

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


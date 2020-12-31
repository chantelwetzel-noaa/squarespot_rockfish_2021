# Squarespot Rockfish  Model Runs

library(r4ss)

wd = "N:/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/models/"

model = "0.01_init_model_updated_catches"

base = SS_output(file.path(wd, model))
SS_plots(base)

# Change age-based selectivity to be age0+ rather than 1+
model = "0.021_selex_age_0"
base = SS_output(file.path(wd, model))
SS_plots(base)

# Update the low bin for para 1 to match data bins
# Allow the model to estimate the init selectivity p5 
model = "0.022_selex_est_p5_survey"
base = SS_output(file.path(wd, model))
SS_plots(base)


# Try to fix the kink in growth at the largest lengths
model = "0.031_L2"
base = SS_output(file.path(wd, model))
SS_plots(base)

# Increase cv young to see if that improves the fits
model = "0.032_cv_young"
base = SS_output(file.path(wd, model))
SS_plots(base)
# this was a mixed bag, improved the trawl survey fits for small fish but decreased the rec fits

# Add growth patterns
model = "0.033_add_gpatterns"
base = SS_output(file.path(wd, model))
SS_plots(base)



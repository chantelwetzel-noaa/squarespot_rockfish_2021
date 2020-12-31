###################################################################################
#
#	Squarespot rockfish 
# Trawl survey exploration and processing
# Look at trawl survey data and the hook & line
#
#############################################################################################

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/nwfscSurvey")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/HandyCode")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/dataModerate_2021")

source("C:/Assessments/2020/survey_summary/code_package/functions/plot_cpue.R")

dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data/"


############################################################################################
#	Load the Survey Data
############################################################################################

hkl = read.csv(file.path(dir, "Hook Line Data", "qryGrandUnifiedThru2019_06182020.csv"))
sub_hkl = hkl[hkl$COMNAME == 'Squarespot Rockfish', ]
sub_hkl = rename_hook_and_line(data = sub_hkl, survey_name = "nwfsc_hkl")
sub_hkl$Length_cm = sub_hkl$Length

load(file.path(dir, "Trawl Survey Catch", "Catch__NWFSC.Combo_2020-07-30.rda"))
load(file.path(dir, "Trawl Survey Bio", "Bio_All_NWFSC.Combo_2020-07-30.rda"))
catch_wcgbt = Out
bio_wcgbt = Data

load(file.path(dir, "Trawl Survey Catch", "Catch__Triennial_2020-07-30.rda"))
load(file.path(dir, "Trawl Survey Bio", "Bio_All_Triennial_2020-07-30.rda"))
catch_tri = Out
bio_tri = Data$Lengths

############################################################################################
# Define the length bins
############################################################################################

len_bins = seq(8, 30, 1)

############################################################################################
# Create visual plots of the survey data
############################################################################################

plot_cpue_fn(dir = file.path(dir, "Trawl Survey Catch", "plots"), 
			 name = "Squarespot rockfish", 
			 catch = catch_wcgbt, bio = bio_wcgbt, 
			 plot = 1:3, n = 20000)

# Does not work for the triennial - likely due to unsexed fish
# plot_cpue_fn(dir = file.path(dir, "Trawl Survey Catch", "plots"), 
# 			 name = "Squarespot rockfish", 
# 			 catch = catch_tri, bio = bio_tri, 
# 			 plot = 1:3, n = 20000)

############################################################################################
# NWFSC WCGBTS data
############################################################################################

# Start with the NWFSC WCGBTS data
strata = CreateStrataDF.fn(names=c("shallow"), 
                           depths.shallow = c( 55),
                           depths.deep    = c(183),
                           lats.south     = c(32.0),
                           lats.north     = c(38.0))

write.csv(strata, file = file.path(dir, "Trawl Survey Catch", "nwfcs_wgbts_strata.csv"))

num.strata = CheckStrata.fn(dir = file.path(dir, "Trawl Survey Catch"), 
							dat = catch_wcgbt, 
							strat.vars = c("Depth_m","Latitude_dd"), 
							strat.df = strata, printfolder = "forSS",  
							verbose = TRUE)

# Calculate the design based index
biomass.nwfsc = Biomass.fn(dir = file.path(dir, "Trawl Survey Catch"), 
						   dat = catch_wcgbt,  
						   strat.df = strata, 
						   printfolder = "forSS", 
						   outputMedian = T) 

# Plot the biomass index
PlotBio.fn(dir = file.path(dir, "Trawl Survey Catch"), dat = biomass.nwfsc, main = "NWFSC WCGBTS", dopng = TRUE)

# Calculate the expanded length composition data 
len = bio_wcgbt

# Calculate the effN
n = GetN.fn(dir = file.path(dir, "Trawl Survey Bio"), dat = len, type = "length", species = "others", printfolder = "forSS")


# This version assigns unsexed fish to a sex - but may not want to do it that way
# See below for alternative versions
#LFs <- SurveyLFs.fn(dir = file.path(dir, "Trawl Survey Bio"), 
#					datL = len, datTows = catch_wcgbt,  
#          strat.df = strata, lgthBins = len_bins, sex = 3, remove999 = FALSE,
#          sexRatioStage = 1, sexRatioUnsexed = 0.5, maxSizeUnsexed = 10, 
#          nSamps = n)

# This version does not do sex assignment for unsexed fish and produces both sexed and unsexed fish comps where:
# sexed fish in the SS file represent only fish that were sexed at sea
# unsexed fish in the SS file represent only unsexed fish 
LFs <- SurveyLFs.fn(dir = file.path(dir, "Trawl Survey Bio"), 
                    datL = len, datTows = catch_wcgbt,  
                    strat.df = strata, lgthBins = len_bins, sex = 3, remove999 = FALSE,
                    sexRatioStage = 999, nSamps = n)

# May want to look at unexpanded lengths since sampling can be low
lfs = UnexpandedLFs.fn(dir = file.path(dir,  "Trawl Survey Bio"), 
                       datL = len, lgthBins = len_bins,
                       sex = 3,  partition = 0, fleet = 2, month = 7)

PlotFreqData.fn(dir = file.path(dir, "Trawl Survey Bio"), 
	  dat = LFs, ylim=c(0, max(len.bins)), inch = 0.10,
      main = "NWFSC WCGBTS", yaxs="i", ylab="Length (cm)", dopng = TRUE)

PlotSexRatio.fn(dir = file.path(dir, "Trawl Survey Bio"), 
	  dat = len, data.type = "length", dopng = TRUE, 
      main = "NWFSC WCGBTS")

############################################################################################
# Triennial Data
############################################################################################
#catch_tri = Out
#bio_tri = Data$Lengths
pos = catch_tri[catch_tri$total_catch_wt_kg >0, ]
par(mfrow = c(3,1))
plot(pos$Latitude_dd, pos$Depth_m)
plot(pos$Latitude_dd, log(pos$cpue_kg_km2))
plot(pos$Depth_m, log(pos$cpue_kg_km2))

table(pos$Year)
#1989 1992 1995 2001 2004 
#   2    4    4    4    1 


#####################################################################################
# Unexpanded Hook & Line Composition Data
#####################################################################################
hk_len = UnexpandedLFs.fn(dir = file.path(dir, "Hook Line Data"), 
                          datL = sub_hkl, lgthBins = len_bins, sex = 3,  
                          partition = 0, fleet = 3,  month = 7, printfolder = "forSS")

PlotFreqData.fn(dir = file.path(dir, "Hook Line Data"), 
      dat = hk_len$comps, ylim=c(0, max(len_bins)), 
      main = "NWFSC HKL", ylab="Length (cm)", dopng = TRUE)


#####################################################################################################
# Run VAST
#####################################################################################################

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/nwfscSurvey")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/VASTWestCoast")
library(VASTWestCoast)
dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data/"
wd = file.path(dir, "Trawl Survey Catch", "vast")

sci_name = "Sebastes_hopkinsi"

# This is probably not necessarily needed
strata.limits = data.frame('STRATA' = c("ca"), 
                          'south_border'   = c(32.0),
                          'north_border'   = c(42.0), 
                          'shallow_border' = c(  55), 
                          'deep_border'    = c(1280) )


obs_model = c(2,0) # do the gamma only for now
obs_model = c(1,0) # lognormal
#cheater way is to do two calls:
#debugonce(VAST_spp)
#VAST_spp(species = "square ..")
#Q

Sim_Settings <- list( "Species" = paste0("WCGBTS_", sci_name), 
                      "ObsModelcondition" = obs_model,
                      #"nknots" = 250,
                      #"depth" = c("no", "linear", "squared")[1],
                      "strata" = strata.limits,
                      "overdispersion" = c("eta1" = 0, "eta2" = 0),
                      "FieldConfig" = c("Omega1"=0, "Epsilon1"=0, "Omega2"="IID", "Epsilon2"=0)
                      #"version" = "VAST_v5_4_0",
                      #"Passcondition" = TRUE,
                      #"overdisperion" = NULL
                      )

Settings_all = list("Species" = paste0("WCGBTS_", sci_name), 
					"overdispersion" = NULL)

Settings_all[["overdispersion"]] <- switch(
      get_spp(Settings_all$Species)["survey"],
      WCGBTS = c("eta1" = 1, "eta2" = 1),
      AFSC.Slope = c("Delta1" = 0, "Delta2" = 0),
      NWFSC.Slope = c("Delta1" = 0, "Delta2" = 0),
      # todo: think about how to write an if statement for this in VAST_do
      Triennial = c("Delta1" = 0, "Delta2" = 0),
      # provide a default to not estimate a vessel-year catchability
      c(eta1 = 0, eta2 = 0))
#Please turn off factor-model variance parameters `L_` that are approaching zero and re-run the model


# VAST_condition calls the VAST_do function  
test <- VAST_condition(conditiondir = wd,
                       settings = Sim_Settings, 
                       spp = Sim_Settings$Species)
                       #datadir = wd,
                       #overdispersion = NULL)

VAST_diagnostics(wd)

# An alternative approach is to use the VAST_spp function which is bigger wrapper around
# VAST_condition and VAST_do.  The VAST_spp works with the nwfscSurvey package to create stratas 
# using the GetSpp.fn function included there.  This function has a suite of pre-defined strata
# selections for common trawl survey species. The GetSpp.fn strata definition is then used by the 
# GetStrata.fn to create the strata matrix and the convert_stata4vast function translate that into
# VAST terms.  

# Example call:
wd = "C:/Assessments/vast/squarespot"
VAST_spp(dir = wd, species = "squarespot_rockfish")

# New function to evaluate and understand if you get an error when running VAST:
get_error(directoryofshittymodel)
# This function will  load every object you need and open a web browser with the R script for the 
# function that created those objects so you can change the settings and rerun the model.
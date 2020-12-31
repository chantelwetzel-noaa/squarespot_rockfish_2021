##############################################################################################################
#
# 	Purpose: Evaluate squarespot rockfish discarding
# 		by source, fishery, and across time.
#
#			  by Chantel Wetzel 
#
##############################################################################################################

library(HandyCode)
library(dplyr)
options(stringsAsFactors = FALSE)


#-----------------------------------------------------------------------------------
# Load the GEMM - the GEMM includes information for commercial 
#-----------------------------------------------------------------------------------

# The gemm splits data north and south of 40 10
dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/copper_rockfish/data/catches"
gemm_all = read.csv(file.path(dir, "GEMM_2019_8_15_2020.csv"))


#---------------------------------------------------------------------------------------------
# Squarespt rockfish in the GEMM data
#---------------------------------------------------------------------------------------------
dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data/PacFIN Catch"
gemm = gemm_all[gemm_all$Species == "Squarespot Rockfish", ]

# Remove the research removals -- 
# Research removals are generally not included with commercial landings (although this does not need to be the case)
# however, removing them here allows you to correctly calculate the discard rate based on commercial data only
# this was an issue because commercial landings for squarespot are fairly small and retaining the research removals
# impacted the rate calculation
#r = gemm[gemm$Sector %in% "Research", ]
gemm = gemm[!gemm$Sector %in% "Research", ] 

aggregate(Landings~Sector, data = gemm, FUN = sum)
aggregate(Discards~Sector, data = gemm, FUN = sum)
aggregate(Discard.Mortality~Sector, data = gemm, FUN = sum)

gemm$grouped_sector = NA
gemm$grouped_sector[gemm$Sector == "California Recreational"] = "ca_rec"
gemm$grouped_sector[is.na(gemm$grouped_sector)] = "commercial"

landings  = aggregate(Landings ~ Year + grouped_sector, data = gemm, drop = FALSE, FUN = sum)
discards  = aggregate(Discards ~ Year + grouped_sector, data = gemm, drop = FALSE, FUN = sum)
disc_mort = aggregate(Discard.Mortality ~ Year + grouped_sector, data = gemm, drop = FALSE, FUN = sum)
all_dead  = aggregate(Mortality..Landings.and.Discard.Mortality. ~ Year + grouped_sector, data = gemm, drop = FALSE, FUN = sum)

all = data.frame(Year = landings$Year,
				 Area = landings$grouped_sector,
				 Landings = landings$Landings,
				 Discard = discards$Discards,
				 Dead_Discard = disc_mort$Discard.Mortality,
				 Tot_Dead = all_dead$Mortality..Landings.and.Discard.Mortality.)

all[is.na(all)] = 0

all$Discard_Mort_Rate = round(all[,"Dead_Discard"] / all[,"Tot_Dead"], 3)
all[is.na(all)] = 0

write.csv(all, file = file.path(dir, "squarespot_gemm_mortality_and_discard.csv"))

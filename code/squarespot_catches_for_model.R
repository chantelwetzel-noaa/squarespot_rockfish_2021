##############################################################################################################
#
# 	Purpose: Summarize Squarespot Rockfish Landing
# 			formate for use by SS
#
#			  by Chantel Wetzel 
#
##############################################################################################################

library(HandyCode)
library(dplyr)
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/dataModerate_2021")


dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data"
dir.create(file.path(dir,"catch_all", "plots"))

# Load the Recreational Data
#-----------------------------------------------------------------------------------
# Miller and Gotshall estimated rec discard in 1965 to be 26.7%
# Ally et al. 1991 estimated rec discard betwen 1985-87 to be 3%
# Apply a ramp down approach in the discard rate between these years 
ca_rec_discard_rate_1965 = 1.267 
# After discussion with Budrick and based on explortation of potential discard rates post-2005
# the estimate from Ally et. al. seems implausibly low -- apply a constant discard rate of 26.7% for 
# all recreational landings prior to the start of MRFSS/RecFIN
#r = (0.267 - 0.03) / (1985 - 1965)
ca_rec_discard_rate = data.frame(year = 1928:1980, rate  = ca_rec_discard_rate_1965)

#step = 1
#for (y in 1966:1985){
#	ca_rec_discard_rate[ca_rec_discard_rate$year == y, "rate"] = ca_rec_discard_rate_1965 - step * r
#	step = step + 1
#}


# California - has both mt and fish numbers (but the fish numbers have decimals)
ca_rec = read.csv(paste0("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_rec_catch_all_2005-2019_crfs.csv"))
ca_rec = ca_rec[ca_rec$Species.Name == "SQUARESPOT ROCKFISH", ]

# mrfss catches
ca_rec_early = read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/mrfss_catch_estimates_1980_2004_final.csv")
ca_rec_early = ca_rec_early[ca_rec_early$AGENCY_CODE == 6 & ca_rec_early$SPECIES_NAME == "Squarespot Rockfish", ]

ca_rec_hist = read.csv(file.path(dir, "RecFIN Catch", "ca_hist_recreational_1928_1980_ej.csv"))
ca_rec_hist$Total = ca_rec_hist[,"SQRSmt_North"] + ca_rec_hist[,"SQRSmt_South"]

# Load the Commercial Data
#-----------------------------------------------------------------------------------

# calculated by averaging the discard ratio from WCGOP data
# this was originally calculated as 20% in error because research removals were not removed when calculating the rate
com_discard_rate = 1.28
# PacFIN Commercial
load(file.path(dir, "PacFIN Catch", "SQRS.CompFT.17.Aug.2020.RData"))
com = SQRS.CompFT.17.Aug.2020

# California Historical Catches 1969 - 1980
com_hist  = read.csv(file.path(dir, "PacFIN Catch", "squarespot_ca_historical_1969_80_ej.csv"))

# California Historical Catches 1916- 1968
com_hist_early  = read.csv(file.path(dir, "PacFIN Catch", "ca_hist_commercial_1916_1968_ej.csv"))

# Grab the WGCOP GEMM data
gemm = read.csv(file.path(dir, "PacFIN Catch", "squarespot_gemm_mortality_and_discard.csv"))


#################################################################################################################
# Evaluate the recreational data first
#################################################################################################################

tmp = aggregate(WGT_AB1..mt. ~ YEAR_, data = ca_rec_early, drop = FALSE, FUN = sum)
colnames(tmp) = c("year", "catch_mt")

# Need to linear interpolate the missing years 1990-1992
x = (tmp[which(tmp$year %in% c(1993)), "catch_mt"] - tmp[which(tmp$year %in% c(1989)), "catch_mt"]) / 4

add = NULL
z = tmp[which(tmp$year == 1989), "catch_mt"]
step = 1
for (y in 1990:1992){
	new = c(y, z + step*x)
	add = rbind(add, new)
	step = step + 1
}
add = as.data.frame(add)
colnames(add) = c('year', 'catch_mt'); rownames(add) = NULL
tmp_early = rbind(tmp, as.data.frame(add))
tmp_early = tmp_early[order(tmp_early$year),]

tmp_late = aggregate(Total.Mortality.MT~RecFIN.Year, data = ca_rec, drop = FALSE, FUN = sum)
colnames(tmp_late) = c('year', 'catch_mt')


# Create a single data frame with all the california rec catches
ca_rec_all = merge(tmp_late, tmp_early, by = c("year",   "catch_mt"), all = TRUE)

# Apply the discard reate to years 1980 and before only
all_rec_ca = data.frame(year = 1928:2020,						
						catch = c(ca_rec_discard_rate[ca_rec_discard_rate$year < 1981, "rate"] * ca_rec_hist[,"Total"], 
								  tmp_early[which(tmp_early$year > 1980 ), "catch_mt"], 
								  tmp_late[, "catch_mt"], 0))


#################################################################################################################
# Evaluate the commercial landings
#################################################################################################################

com$round_mt = com$CATCH.LBS / 2204.62262

tmp_com = aggregate(round_mt ~ YEAR, data = com, FUN = sum, drop = FALSE)
tmp_com$total_mt = NA
colnames(tmp_com) = c('year', 'round_mt', 'total_mt')
# Missing any removals for 2001, 2002, 2006, 2007 - according to the GEMM the landings in 2002 was 0.
tmp_com = rbind(tmp_com, c(2001, 0.1, 0), c(2002, 0, 0), c(2006, 0, 0), c(2007, 0, 0))
tmp_com = tmp_com[order(tmp_com$year),]

#add discards by year
tmp_com[which(tmp_com$year %in% gemm$Year),'total_mt'] = tmp_com[which(tmp_com$year %in% gemm$Year),'round_mt'] +
	gemm[gemm$Area == "commercial","Dead_Discard"]

# Need to get the amount for 2020 from the GMT
tmp_com[which(tmp_com$year == 2020), "total_mt"] =  tmp_com[which(tmp_com$year == 2020), "round_mt"]

# Apply the default discard rate to PacFIN years pre-GEMM
tmp_com[which(tmp_com$year %in% 1981:2001), 'total_mt'] = 
	com_discard_rate * tmp_com[which(tmp_com$year %in% 1981:2001), 'round_mt']

# Need to adjust the historical catches for north and south
# 1916 - 1968 Ralston et al. removals
com_hist_early$total = com_discard_rate * com_hist_early[,"SQRSmt"]

# 1969 - 1980 CalCOM removals
com_hist$total = com_discard_rate * com_hist[,"catch.mt"]

# Knit all the data sources together in a ugly fashion
all_com = matrix(0, length(1916:2020), 2)
colnames(all_com) = c('year', 'catch')

all_com[,"year"] = 1916:2020
# California 1916 - 1968
all_com[which(all_com[,"year"] %in% com_hist_early$Year), "catch"] = com_hist_early[, "total"]
# California 1969-1980
all_com[which(all_com[,"year"] %in% com_hist$Year), "catch"] = com_hist[, "total"]
# California 1981 - 2020
all_com[which(all_com[,"year"] %in% tmp_com$year), "catch"] = tmp_com[, "total_mt"]

#################################################################################################################
# For the catches for SS by each model area:
# Year, Season, Fleet, Catch, Catch SE
#################################################################################################################

catch_ss = rbind( cbind(as.numeric(all_com[,"year"]), 1, 1, round(as.numeric(all_com[,"catch"]), 2), 0.01),
						cbind(as.numeric(all_rec_ca$year), 1, 2, round(as.numeric(all_rec_ca$catch), 2), 0.01))
colnames(catch_ss) = c("Year", "Season", "Fleet", "Catch", "SE")

write.csv(catch_ss, file = file.path(dir, "catch_all", "forSS", "catches_for_ss.csv"), row.names = FALSE)











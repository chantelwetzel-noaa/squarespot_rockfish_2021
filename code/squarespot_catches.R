##############################################################################################################
#
# 	Purpose: Summarize Squarespot Rockfish Landing
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
# California - has both mt and fish numbers (but the fish numbers have decimals)
ca_rec = read.csv(paste0("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_rec_catch_all_2005-2019_crfs.csv"))
ca_rec = ca_rec[ca_rec$Species.Name == "SQUARESPOT ROCKFISH", ]

# mrfss catches
ca_rec_early = read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/mrfss_catch_estimates_1980_2004_final.csv")
ca_rec_early = ca_rec_early[ca_rec_early$AGENCY_CODE == 6 & ca_rec_early$SPECIES_NAME == "Squarespot Rockfish", ]



# Load the Commercial Data
#-----------------------------------------------------------------------------------
# PacFIN Commercial
load(file.path(dir, "PacFIN Catch", "SQRS.CompFT.17.Aug.2020.RData"))
com = SQRS.CompFT.17.Aug.2020

# California Historical Catches
com_hist  = read.csv(file.path(dir, "PacFIN Catch", "squarespot_ca_historical_1969_80_ej.csv"))

# Grab the WGCOP GEMM data
gemm = read.csv(file.path(dir, "PacFIN Catch", "squarespot_gemm_mortality_and_discard.csv"))



#################################################################################################################
# Evaluate the recreational data first
#################################################################################################################

# California - RecFIN years
ca_rec$areas = recfin_areas(data = ca_rec, 
			  area_grouping = list(c("CHANNEL", "SOUTH"), c("BAY AREA", "WINE", "CENTRAL", "REDWOOD")), 
			  area_names  = c("south", "north"), 
			  column_name = "RecFIN.Port.Name")

mode_grouping = list(c("BEACH/BANK", "MAN-MADE/JETTY"), c("PARTY/CHARTER BOATS", "PRIVATE/RENTAL BOATS"))
mode_names = c("rec_shore", "rec_boat")
mode_column_name = "RecFIN.Mode.Name"
ca_rec$mode <- NA
if(!is.null(mode_grouping)){
	for (a in 1:length(mode_grouping)){
		get <- paste(mode_grouping[[a]], collapse = "|")
		find = grep(get, ca_rec[, mode_column_name], ignore.case = TRUE)
		ca_rec$mode[find] = mode_names[a]
	}
}
# all the pacfin data are from either charter boats or private vessels

# Look at the number of records by area
table(ca_rec$RecFIN.Year, ca_rec$areas)
# Almost of the recreeational observations are coming from the southern area

# California - MRFSS years
ca_rec_early$areas = recfin_areas(data = ca_rec_early, 
			  area_grouping = list("Northern California", "Southern California"), 
			  area_names  = c( "north", "south"), 
			  column_name = "RECFIN_SUB_REGION_NAME")

mode_grouping = list(c("Beach/Bank", "Man-Made/Jetty", "Shore"), c("Party/Charter Boats", "Private/Rental Boats"))
mode_names = c("rec_shore", "rec_boat")
mode_column_name = "RECFIN_MODE_NAME"
ca_rec_early$mode <- NA
if(!is.null(mode_grouping)){
	for (a in 1:length(mode_grouping)){
		get <- paste(mode_grouping[[a]], collapse = "|")
		find = grep(get, ca_rec_early[, mode_column_name], ignore.case = TRUE)
		ca_rec_early$mode[find] = mode_names[a]
	}
}

# Look at the number of records by area
table(ca_rec_early$YEAR, ca_rec_early$areas)

check = ca_rec[,"Total.Mortality.MT"]/ (ca_rec[, "Retained.MT"] + ca_rec[, "Released.Dead.MT"])
plot(check)
abline(h = 1, col = 'red')
ca_rec$diff = ca_rec[,"Total.Mortality.MT"] / (ca_rec[, "Retained.MT"] + ca_rec[, "Released.Dead.MT"])
# It looks like these columns may not match due to rounding in the preliminary calcs.

# Cut down the data frame to a manageble size of key items
tmp = data.frame(year = ca_rec_early$YEAR, 
			     areas = ca_rec_early$areas, 
			     mode = ca_rec_early$mode, 
			     catch_mt = ca_rec_early$WGT_AB1..mt.)

tmp_late = data.frame(year = ca_rec$RecFIN.Year, 
			     areas = ca_rec$areas, 
			     mode = ca_rec$mode, 
			     catch_mt = ca_rec$Total.Mortality.MT)

# Create a single data frame with all the california rec catches
ca_rec_all = merge(tmp_late, tmp, by = c("year", "mode", "areas", "catch_mt"), all = TRUE)
a = aggregate(catch_mt ~ year + areas , data = ca_rec_all, FUN = sum, drop = FALSE)
ca_rec_df = data.frame(year = a$year,
					   area = a$areas,
					   catch_mt = a$catch_mt)
ca_rec_df[is.na(ca_rec_df)] = 0


# Look at the mortality in MT
tot = aggregate(Total.Mortality.MT ~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
ret = aggregate(Retained.MT~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
#rel_alive = aggregate(Released.Alive.MT~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
rel_dead  = aggregate(Released.Dead.MT~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)

rec_mt = cbind(tot, ret$Retained.MT, rel_dead$Released.Dead.MT)
colnames(rec_mt) = c("year", "areas", "tot_mort", "retsained_mt", "rel_dead_mt")
rec_mt[is.na(rec_mt)] = 0

# Now lets compare this information by the Num category
tot_num = aggregate(Total.Mortality.Num ~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
ret_num = aggregate(Retained.Num~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
rel_dead_num  = aggregate(Released.Dead.Num~ RecFIN.Year + areas, data = ca_rec, FUN = sum, drop = FALSE)
rec_num = cbind(tot_num, ret_num$Retained.Num, rel_dead_num$Released.Dead.Num)
colnames(rec_num) = c("year", "areas", "tot_mort", "retained_num", "rel_dead_num")
rec_num[is.na(rec_num)] = 0

pngfun(wd = file.path(dir, "catch_all", "plots"), file = "CA_rec_mt_vs_num_2005_2019.png")
par(mfrow = c(2,2))
plot(rec_mt[rec_mt$areas == "south", "year"], rec_mt[rec_mt$areas == "south", "tot_mort"], type = 'l', lwd = 2, xlab = "Year", 
	ylab = "Total Mortality (mt)", ylim = c(0, 50), main = "South of Pt. Conception")
plot(rec_mt[rec_mt$areas == "north", "year"], rec_mt[rec_mt$areas == "north", "tot_mort"], type = 'l', lwd = 2, xlab = "Year", 
	 ylab = "Total Mortality (mt)", ylim = c(0, 50), main = "North of Pt. Conception")

plot(rec_num[rec_num$areas == "south", "year"], rec_num[rec_num$areas == "south", "tot_mort"], type = 'l', lwd = 2, xlab = "Year", 
	 ylab = "Total Mortality (numbers)", ylim = c(0, 225000), main = "South of Pt. Conception")
plot(rec_num[rec_num$areas == "north", "year"], rec_num[rec_num$areas == "north", "tot_mort"], type = 'l', lwd = 2, xlab = "Year", 
	ylab = "Total Mortality (numbers)", ylim = c(0, 225000), main = "North of Pt. Conception")
dev.off()

# Budrick suggested just looking at metric ton values
pngfun(wd = file.path(dir, "catch_all", "plots"), file = "Rec_Catch_Area_2005_2019.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1,1))
barplot(rec_mt[rec_mt$areas == "south", "tot_mort"] + rec_mt[rec_mt$areas == "north", "tot_mort"], main = "Recreational",
	col = "red", ylim = c(0, 25), names.arg = rec_mt[rec_mt$areas == "north", "year"])
barplot(rec_mt[rec_mt$areas == "south", "tot_mort"] , col = "purple", add = TRUE)
legend("topleft", bty = 'n', legend = c("North of Pt. Conception", "South of Pt. Conception"),
	lwd = 2, lty = 1, col = c('red', 'purple'))
mtext(side = 1, "Year", line = 2.5)
mtext(side = 2, "Total Mortality (mt)", line = 2.5)
dev.off()

#--------------------------------------------------------------------------------------------------------------
# Plot all the recreational catches for California
#--------------------------------------------------------------------------------------------------------------
pngfun(wd = file.path(dir, "catch_all", "plots"), file = "CA_Rec_Catch.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1,1))
barplot(ca_rec_df[ca_rec_df$area == 'north', "catch_mt"] + ca_rec_df[ca_rec_df$area == 'south', "catch_mt"], 
	main = "Recreational Catch", col = "blue", ylim = c(0, 25), 
	names.arg = ca_rec_df[ca_rec_df$area == 'north', "year"])
barplot(ca_rec_df[ca_rec_df$area == 'south', "catch_mt"],  col = "red", add = TRUE)
mtext(side = 1, "Year", line = 2.5)
mtext(side = 2, "Total Mortality (mt)", line = 2.5)
legend("topleft", bty = 'n', legend = c("North of Pt. Conception", "South of Pt. Conception"),
	lwd = 6, lty = 1, col = c('blue', 'red'))
dev.off()

write.csv(ca_rec_df, file = file.path(dir, "catch_all", "squarespot_recreational_total_mortality_by_area.csv"))

tmp = aggregate(catch_mt ~ year, data = ca_rec_df, FUN = sum)
write.csv(tmp, file = file.path(dir, "catch_all", "squarespot_recreational_total_mortality.csv"))

#################################################################################################################
# Evaluate the commercial landings
#################################################################################################################

# Look at landed alive vs. dead
table(com$DISP)
#   F     H     I     P 
#   1 10491     7    28 
# Only 1 record of an alive fish landed in the pacfin data.

com$round_mt = com$CATCH.LBS / 2204.62262

# Evaluate the amount of catch by areas south and north of pt conception
# The model will most likely be one area but just checking the patterns
# I think these are the areas that are south of point conception
pcid_south = c("DNA", "HNM", "LGB", "NWB", "OBV", "OLA", "OSD", "OXN", "SB", "SD", "SP", "TRM", "VEN", "WLN")
pcid = sort(unique(com$PCID))
pcid_north = pcid[!(pcid %in% pcid_south)]
com$area = NA
com$area[com$PCID %in% pcid_south] = "south"
com$area[com$PCID %in% pcid_north] = "north"

catch_area = aggregate(round_mt ~ YEAR + area, data = com, FUN = sum, drop = FALSE)
catch_area[is.na(catch_area)] = 0
catch_area_df = data.frame(year = sort(unique(catch_area$YEAR)),
						   north = catch_area[catch_area$area == "north", "round_mt"],
						   south = catch_area[catch_area$area == "south", "round_mt"],
						   percent_south = catch_area[catch_area$area == "south", "round_mt"] / 
						   				   (catch_area[catch_area$area == "south", "round_mt"] + catch_area[catch_area$area == "north", "round_mt"]))


# Evaluate the amount of catch by gear
# I suspect there will be one lumped "commercail" fleet but lets see what the patterns are
table(com$GRGROUP)
#  HKL  NET  POT  TWL 
# 5252 4706    3  576 

tmp = aggregate(com$round_mt, list(year = com$YEAR, gear = com$GRGROUP), FUN = sum, drop = FALSE)
catch_gear_df = data.frame(year = sort(unique(tmp$year)),
						   hkl = tmp[tmp$gear == "HKL", "x"],
						   net = tmp[tmp$gear == "NET", "x"], 
						   pot = tmp[tmp$gear == "POT", "x"], 
						   twl = tmp[tmp$gear == "TWL", "x"])
catch_gear_df[is.na(catch_gear_df)] = 0
percent_by_gear = round(catch_gear_df[ , 2:dim(catch_gear_df)[2]] / 
					apply(catch_gear_df[ , 2:dim(catch_gear_df)[2]], 1, sum), 2)

pngfun(wd = file.path(dir, "catch_all", "plots"), file = "Commercial_Landings_Percent_Gear.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1,1))
barplot(apply(percent_by_gear,1,sum), col = "red", ylim = c(0,1.3), 
	names.arg = catch_area_gear_df[,1])
barplot(apply(percent_by_gear[,2:4],1,sum), col = "purple", add = TRUE)
barplot(apply(percent_by_gear[,3:4],1,sum), col = "green", add = TRUE)
barplot(percent_by_gear[,4],, col = "orange", add = TRUE)
legend("topleft", bty = 'n', legend = c("HKL", "NET", "POT", "TWL"),
	lwd = 2, lty = 1, col = c('red', 'purple', 'green', 'orange'))
mtext(side = 2, "Percent of Total by Gear", line = 2.5)
dev.off()

# apply(catch_gear_df[,-1], 2, sum, na.rm = TRUE)
#        hkl        net        pot        twl 
# 2.19914679 1.19479904 0.01678292 2.13667457 
# By percent:
# apply(percent_by_gear, 2, mean)
#       hkl       net       pot       twl 
# 0.6006667 0.1096667 0.0110000 0.2786667

tmp = aggregate(round_mt ~ YEAR + area + GRGROUP, data = com, FUN = sum, drop = FALSE)
tmp[is.na(tmp)] = 0
catch_area_gear_df = data.frame(year = sort(unique(tmp$YEAR)),
						        n_hkl = tmp[tmp$area == "north" & tmp$GRGROUP == "HKL", "round_mt"],
						        n_net = tmp[tmp$area == "north" & tmp$GRGROUP == "NET", "round_mt"],
						        n_pot = tmp[tmp$area == "north" & tmp$GRGROUP == "POT", "round_mt"],
						        n_twl = tmp[tmp$area == "north" & tmp$GRGROUP == "TWL", "round_mt"],
						        s_hkl = tmp[tmp$area == "south" & tmp$GRGROUP == "HKL", "round_mt"],
						        s_net = tmp[tmp$area == "south" & tmp$GRGROUP == "NET", "round_mt"],
						        s_pot = tmp[tmp$area == "south" & tmp$GRGROUP == "POT", "round_mt"],
						        s_twl = tmp[tmp$area == "south" & tmp$GRGROUP == "TWL", "round_mt"])

#apply(catch_area_gear_df[,-1], 2, sum)
# n_hkl n_net n_pot n_twl s_hkl s_net s_pot s_twl 
# 0.151 0.011 0.000 2.137 2.048 1.184 0.017 0.000 
pngfun(wd = file.path(dir, "catch_all", "plots"), file = "Commercial_Landings_Gear_Area.png", w = 12, h = 7, pt = 12)
par(mfrow = c(2, 1))
barplot(apply(catch_area_gear_df[,2:5],1,sum), col = "red", ylim = c(0, 0.8), main = "North of Pt. Conception",
	names.arg = catch_area_gear_df[,1])
barplot(apply(catch_area_gear_df[,3:5],1,sum), col = "purple", add = TRUE)
barplot(apply(catch_area_gear_df[,4:5],1,sum), col = "green", add = TRUE)
barplot(catch_area_gear_df[,5], col = "orange", add = TRUE)
mtext(side = 2, "Landings (mt)", line = 2.5)
legend("topright", bty = 'n', legend = c("HKL", "NET", "POT", "TWL"),
	lwd = 2, lty = 1, col = c('red', 'purple', 'green', 'orange'))

barplot(apply(catch_area_gear_df[,6:9],1,sum), col = "red", ylim = c(0, 0.8), main = "South of Pt. Conception",
	names.arg = catch_area_gear_df[,1])
barplot(apply(catch_area_gear_df[,7:9],1,sum), col = "purple", add = TRUE)
barplot(apply(catch_area_gear_df[,8:9],1,sum), col = "green", add = TRUE)
barplot(catch_area_gear_df[,9], col = "orange", add = TRUE)
mtext(side = 2, "Landings (mt)", line = 2.5)
dev.off()


tmp = aggregate(round_mt ~ YEAR + area , data = com, FUN = sum, drop = FALSE)
tmp[is.na(tmp)] = 0
catch_area_df = data.frame(year = sort(unique(tmp$YEAR)),
						        north = tmp[tmp$area == "north", "round_mt"],
						        south = tmp[tmp$area == "north", "round_mt"])




#---------------------------------------------------------------------------
# Combine the historical catches to the PacFIN years
#---------------------------------------------------------------------------

all_com = matrix(NA, length(1969:2020), 1)
colnames(all_com) = c("catch_mt")
rownames(all_com) = 1969:2020

all_com[rownames(all_com) %in% catch_area_df$year, "catch_mt"] = catch_area_df[,"south"] + catch_area_df[,"north"]
all_com[rownames(all_com) %in% com_hist$Year, "catch_mt"] = com_hist$catch.mt
all_com[is.na(all_com)] = 0

write.csv(all_com, file = file.path(dir, "catch_all", "squarespot_commercial_landings.csv"))


pngfun(wd = file.path(dir, "catch_all", "plots"), file = "Commercial_Landings.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1, 1))
barplot(all_com[,"catch_mt"], col = "red", ylim = c(0, 1.5), las = 1,
	names.arg = rownames(all_com))
mtext(side = 2, "Landings (mt)", line = 2.5, cex = 1.5)
mtext(side = 1, "Year", line = 2.5, cex = 1.5)
dev.off()


#################################################################################################################
# Create plot to visualize the GEMM data
#################################################################################################################
pngfun(wd = file.path(dir, "catch_all", "plots"), file = "Discard_Data_from_GEMM.png", w = 14, h = 7, pt = 12)
par(mfrow = c(1,2), mar = c(4,5,4,1), oma = c(1,1,1,1))
barplot(gemm[gemm$Area == "commercial", "Dead_Discard"], las = 1, cex.axis = 1.5,
	col = "red", ylim = c(0, 1.1), 
	names.arg = gemm[gemm$Area == "commercial" ,"Year"])
mtext(side = 3, "Coastwide Commercial Discard: WCGOP", line = -2, outer = TRUE, cex = 1.5)
mtext(side = 1, "Year", line = 3, cex = 1.5)
mtext(side = 2, "Total Dead Discard (mt)", line = 3, cex = 1.5)
plot(gemm[gemm$Area == "commercial", "Year"], gemm[gemm$Area == "commercial", "Discard_Mort_Rate"], las = 1,
	 type = 'p', col = 'darkgrey', cex = 2, cex.axis = 1.5, cex.lab = 1.5, ylab = "", xlab = "Year", ylim = c(0, 1))
points(gemm[gemm$Area == "commercial", "Year"], gemm[gemm$Area == "commercial", "Discard_Mort_Rate"], pch = 16, cex = 2, col = 'darkgrey')
mtext(side = 2, "Discard Rate", line = 3.5, cex = 1.5)
dev.off()

#---------------------------------------------------------------------------
# Historical commercial catches
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
# Combine the historical catches to the PacFIN years
#---------------------------------------------------------------------------





###################################################################################
#
#			Squarespot Rockfish
#
#############################################################################################

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/nwfscSurvey")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/PacFIN.Utilities")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/HandyCode")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/dataModerate_2021")

library(ggplot2)

dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish"

pacfin_abbr = "SQRS"
hkl_name = "Squarespot Rockfish"
recfin_name = "SQUARESPOT ROCKFISH"
or_rec_name = "Squarespot"
ca_mrfss_code = 8826010148


############################################################################################
#	Load Data
############################################################################################
load(file.path(dir, "data", "PacFIN BDS", "PacFIN.SQRS.bds.13.Aug.2020.RData"))
pacfin 	 = PacFIN.SQRS.bds.13.Aug.2020

south_ca = c("DNA","HNM","LGB","NWB","OBV","OLA","OSD","OXN","SB","SD","SP","TRM","VEN","WLM")
north_ca = c("ALB","ALM","ARE","AVL","BDG","BKL","BOL","BRG","CRS","CRZ","ERK",
			 "FLN","MNT","MOS","MRO","OAK","OCA","OCM","OCN","ODN","OHB","OMD","OSF","OSL","OSM","OWC","PRN","RCH","RYS","SF","SLT","TML","TRN")
#GFS, GFT = trawl 14 records 
#GLN = gill net 11 records
#HKL = hook & line 85 records
#LGL = longline 31 records

pacfin = rename_pacfin(data = pacfin,
					   area_grouping = list(south_ca, north_ca),
					   area_names = c("south_pt_concep", "north_pt_concep"),
					   fleet_grouping = list(c("GFS", "GFT"), c("GLN"), c("LGL", "HKL")), 
					   fleet_names = c("trawl", "gill_net", "hkl_longline"), 
					   fleet_column_name = "GEAR")

# There are some unexpected large fish (>) in the pacfin data sampled in 1985 only from the gill net gear.
# Need to decide if these lengths should be tossed given that they are well outside expectations.

#recfin = read.csv(file.path(dir, "data", "RecFIN Sample Data", "squarespot_SD051--1980-2019.csv"))
#newrec = read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/data_from_Budrick/ca_rec_retained_lengths_all_2004-2019.csv")
#newrec = newrec[newrec$Species.Name == "SQUARESPOT ROCKFISH",]

##RecFIN
#California
ca_recfin = rename_budrick_recfin(read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_rec_lengths_2004_2020_updated.csv", header=T, na.strings = "-"))
ca_recfin = ca_recfin[ca_recfin$SPECIES_NAME == recfin_name, ]
ca_recfin_data = rename_recfin(data = ca_recfin,
                               area_grouping = list(c("CHANNEL", "SOUTH"), c("BAY AREA", "WINE", "CENTRAL", "REDWOOD", "NOT KNOWN")),
                               area_names = c("south_pt_concep", "north_pt_concep"),
                               area_column_name = "RECFIN_PORT_NAME",
                               mode_grouping = list(c("BEACH/BANK", "MAN-MADE/JETTY"), c("PARTY/CHARTER BOATS", "PRIVATE/RENTAL BOATS"), "NOT KNOWN"),
                               mode_names = c("shore", "boat", "unknown"),
                               mode_column_name = "RecFIN.Mode.Name" )
#table(ca_recfin_data$RecFIN.Mode.Name,ca_recfin_data$Fleet)
#table(ca_recfin_data$State,ca_recfin_data$State_Areas)
#table(ca_recfin_data$SPECIES_NAME,useNA="always")

##MRFSS
#California
ca_mrfss_full = read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_mrfss_bio_1980_2003.csv")
ca_mrfss = ca_mrfss_full[ca_mrfss_full$ST == 6 & ca_mrfss_full$SP_CODE == ca_mrfss_code, ]
ca_mrfss = ca_mrfss[!is.na(ca_mrfss$CNTY), ] # remove records without a county
spc = c(59, 73, 37, 111, 83)
npc = unique(ca_mrfss[!ca_mrfss$CNTY %in% spc, "CNTY"]) 
ca_mrfss$STATE_NAME = "CA"
ca_mrfss_data = rename_mrfss(data = ca_mrfss,
                             len_col = "LNGTH", #"T_LEN",
                             area_grouping = list(spc, npc), 
                             area_names = c("south_pt_concep", "north_pt_concep"), 
                             area_column_name = "CNTY", 
                             mode_grouping = list(c(1,2), c(6, 7)),
                             mode_names = c("shore", "boat"),
                             mode_column_name = "MODE_FX" )
#table(ca_mrfss_data$DIST,ca_mrfss_data$State_Areas)
#table(ca_mrfss_data$MODE_FX,ca_mrfss_data$Fleet)


hkl = read.csv(file.path(dir, "data", "Hook Line Data", 
      "qryGrandUnifiedThru2019_For2021Assessments_DWarehouse version_12142020.csv"))
sub_hkl = hkl[hkl$common_name == 'Squarespot Rockfish', ]
sub_hkl = rename_hook_and_line(data = sub_hkl, survey_name = "nwfsc_hkl")

##Combo Survey
load(file.path(dir, "data", "Trawl Survey Bio", "Bio_All_NWFSC.Combo_2020-12-31.rda"))
combo = Data
rm(Data)
combo_data <- rename_survey_data(data = combo,
                   area_split = c(32.5, 42, 46), 
                   area_names = c("south_pt_concep", "north_pt_concep", "OR", "WA") )
combo_data = combo_data[combo_data$State_Areas != "OR", ]

##Triennial Survey
load(file.path(dir, "data", "Trawl Survey Bio", "Bio_All_Triennial_2020-07-30.rda"))
trien = Data
rm(Data)
trien_data <- rename_survey_data(data = trien[[1]],
                                 area_split = c(32.5, 42, 46), 
                                 area_names = c("south_pt_concep", "north_pt_concep", "OR", "WA") )
table(trien_data$State,trien_data$State_Areas)


input = list()
input[[1]] = combo_data
input[[2]] = sub_hkl
input[[3]] = trien_data
input[[4]] = pacfin
input[[5]] = ca_recfin_data
input[[6]] = ca_mrfss_data

############################################################################################
#	Create data frame with all the input data
############################################################################################
out = create_data_frame(data_list = input)

max = as.numeric(quantile(out$Length, na.rm = TRUE, 0.999))
# Remove the commercial lengths from 1985 that do not appear to be right
# Received confirmation during pre-assessment data meeting
rm = which(out$Length > 40)
out = out[-rm, ] 

############################################################################################
#	Look at the released vs. retained fish prior to removal 
############################################################################################

# table(out[out$Data_Type == "RELEASED", "State_Areas"])
# north_pt_concep south_pt_concep 
#              34             369

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Retained_vs_Released_Lengths.png", w = 7, h = 7, pt = 12)
tmp = out[out$Source == "RecFIN", ]
ggplot(tmp, aes(Length, fill = Data_Type, color = Data_Type)) +
	geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5)
dev.off()

# Remove the released for the rest of the summaries for now:
out = out[out$Data_Type %in% c("RETAINED", NA), ]

############################################################################################
#	Create data frame with all the input data
############################################################################################

sum_data = summarize_data(dir = file.path(dir, "data", "biology"), data = out)

############################################################################################
#	Plot length-at-weight data by source and year
############################################################################################

length_weight_plot(dir = file.path(dir, "data", "biology"), data = out)


############################################################################################
# Estimate Growth Using only Survey data
############################################################################################

survey_dat <- out[out$Source %in% c("NWFSC_WCGBTS", "NWFSC_HKL", "Triennial"),]

est_growth <- estimate_length_weight(data = survey_dat)

save(est_growth, file = file.path(dir, "data", "biology", "growth_estimates_survey.Rdat"))

length_weight_plot(dir = file.path(dir, "data", "biology"),
				   nm_append = "Survey", data = survey_dat, ests = est_growth)

############################################################################################
#	Sex Ratio by Age or Length
############################################################################################
out$round_length = round(out$Length,0)
temp = table(out$round_length, out$Sex)
ratioF = temp[,"F"] / (temp[,"M"] + temp[,"F"])
nobs = temp[,"F"] + temp[,"M"]

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Length_fraction_female.png", w = 7, h = 7, pt = 12)
par(mfrow = c(1,1))
plot(x = names(ratioF), y = ratioF, type="l", col="red", 
	xlab = "Length (cm)", ylab="Fraction female")
abline(h = 0.50, col = "grey", lty = 2, lwd = 2)
symbols(x = names(ratioF), y = ratioF, circles = nobs, 
	inches = 0.1, fg="red", bg = rgb(1,0,0, alpha=0.5), add=T)
dev.off()


age_tmp = out[which(!is.na(out$Age)), ]
age_tmp$source_state = paste0(age_tmp$Source, "_", age_tmp$Sex)
tbl = table(age_tmp$Year, age_tmp$source_state)
write.csv(tbl, file = file.path(dir, "data", "biology", "age_samples_by_sex_source.csv"),
		  row.names = TRUE)

############################################################################################
# Estimate Length-at-Age: These estimates will be based on OR & WA rec boat samples 
############################################################################################

len_age <- estimate_length_age(data = out)

save(len_age, file = file.path(dir, "data", "biology", "length_age_estimates.Rdat"))

line_col = c("red", 'blue', 'grey')
sex_col = alpha(line_col, 0.5)
keep = which(!is.na(out$Age))
tmp = out[keep, ]
lens = 1:max(tmp$Length, na.rm = TRUE)
xmax = max(tmp$Age,    na.rm = TRUE)
ymax = max(tmp$Length + 2, na.rm = TRUE)

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "doc_Length_Age_by_Sex_RE.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1, 1))
plot(tmp[tmp$Sex == "F", "Age"], tmp[tmp$Sex == "F", "Length"], xlab = "Age", ylab = "Length (cm)",
	xaxs = "i", yaxs = "i",ylim = c(0, ymax), xlim = c(0, xmax), pch = 16, col = sex_col[1]) 
points(tmp[tmp$Sex == "M", "Age"], tmp[tmp$Sex == "M", "Length"], pch = 16, col = sex_col[2])
lines(0:xmax, 26.74 * (1 - exp(-0.1239 * (0:xmax + 2.845))), 
	col = line_col[1], lty = 1, lwd = 2)
lines(0:xmax, 20.82 * (1 - exp(-0.2459 * (0:xmax + 1.663))), 
	col = line_col[2], lty = 1, lwd = 2)	
leg = c(paste0("F : Linf = 26.7, k = 0.124"), 
		paste0("M : Linf = 20.8, k = 0.246")) 
legend("bottomright", bty = 'n', legend = leg, lty = c(1,1), col = c(rep(line_col,2)), lwd = 2, cex = 1.25)
dev.off()

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "doc_Length_Age_by_Sex.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1, 1))
plot(tmp[tmp$Sex == "F", "Age"], tmp[tmp$Sex == "F", "Length"], xlab = "Age", ylab = "Length (cm)",
	xaxs = "i", yaxs = "i",ylim = c(0, ymax), xlim = c(0, xmax), pch = 16, col = sex_col[1]) 
points(tmp[tmp$Sex == "M", "Age"], tmp[tmp$Sex == "M", "Length"], pch = 16, col = sex_col[2])
lines(0:xmax, vb_fn(age = 0:xmax, Linf = len_age$all_F[1], L0 = len_age$all_F[2], k = len_age$all_F[3]), 
	col = line_col[1], lty = 1, lwd = 2)
lines(0:xmax, vb_fn(age = 0:xmax, Linf = len_age$all_M[1], L0 = len_age$all_M[2], k = len_age$all_M[3]), 
	col = line_col[2], lty = 1, lwd = 2)	
leg = c(paste0("F : Linf = ", round(len_age$all_F[1], 1),  " L1 = ", round(len_age$all_F[2], 1)," k = ", round(len_age$all_F[3], 3)),
		paste0("M : Linf = ", round(len_age$all_M[1], 1),  " L1 = ", round(len_age$all_M[2], 1)," k = ", round(len_age$all_M[3], 3)))
legend("bottomright", bty = 'n', legend = leg, lty = c(1,1), col = c(rep(line_col,2)), lwd = 2, cex = 1.25)
dev.off()

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "doc_Length_Age_Method_Comparison.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1, 1))
plot(tmp[tmp$Sex == "F", "Age"], tmp[tmp$Sex == "F", "Length"], xlab = "Age", ylab = "Length (cm)",
	xaxs = "i", yaxs = "i",ylim = c(0, ymax), xlim = c(0, xmax), pch = 16, col = sex_col[1]) 
points(tmp[tmp$Sex == "M", "Age"], tmp[tmp$Sex == "M", "Length"], pch = 16, col = sex_col[2])
lines(0:xmax, vb_fn(age = 0:xmax, Linf = len_age$all_F[1], L0 = len_age$all_F[2], k = len_age$all_F[3]), 
	col = line_col[1], lty = 1, lwd = 2)
lines(0:xmax, vb_fn(age = 0:xmax, Linf = len_age$all_M[1], L0 = len_age$all_M[2], k = len_age$all_M[3]), 
	col = line_col[2], lty = 1, lwd = 2)	
lines(0:xmax, 26.74 * (1 - exp(-0.1239 * (0:xmax + 2.845))), 
	col = line_col[1], lty = 2, lwd = 2)
lines(0:xmax, 20.82 * (1 - exp(-0.2459 * (0:xmax + 1.663))), 
	col = line_col[2], lty = 2, lwd = 2)
legend("bottomright", bty = 'n', legend = c("Sum of Squares", "Random Effects"), lty = c(1,2), col = c(rep(line_col,2)), lwd = 2, cex = 1.25)
dev.off()

source_col = c(rgb(red = 1, green = 160/255, blue=0, alpha = 0.80),
			   rgb(red = 160/255, green = 32/255, blue=240/255, alpha = 0.80))
pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "doc_Data_Length_Age_by_Source.png", w = 12, h = 7, pt = 12)
par(mfrow = c(1, 1))
ind = tmp$Source == "NWFSC_WCGBTS"
plot(tmp[ind, "Age"], tmp[ind, "Length"], xlab = "Age", ylab = "Length (cm)",
	xaxs = "i", yaxs = "i",ylim = c(0, ymax), xlim = c(0, xmax+1), pch = 16, col = source_col[1]) 
ind = which(tmp$Source == "NWFSC_HKL")
points(tmp[ind, "Age"], tmp[ind, "Length"], pch = 16, col = source_col[2])
legend("bottomright", bty = 'n', 
		legend = c('NWFSC WCGBT', 'NWFSC HKL'), 
		col = source_col, pch = c(16, 16), cex = 1.25)
dev.off()



############################################################################################
# Create a combined length-weight-plot with the Love estimates added
############################################################################################
data = survey_dat
ests = est_growth
lens = 1:max(data$Length, na.rm = TRUE)
ymax = max(data$Weight, na.rm = TRUE)
xmax = max(data$Length, na.rm = TRUE)

line_col = c("red", 'blue', "grey")
sex_col = c(alpha(line_col[1:2], 0.6), alpha(line_col[3], 0.20))

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Length_Weight_All_w_Love_Ests.png", w = 7, h = 7, pt = 12)
par(mfrow = c(1, 1) )
plot(out[out$Sex == 'U', "Length"], out[out$Sex == 'U', "Weight"], las = 1,
	cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
     xlab = "Length (cm)", ylab = "Weight (kg)", main = "", 
     ylim = c(0, ymax), xlim = c(0, xmax), pch = 16, col = sex_col[3]) 
points(out[out$Sex == 'F', "Length"], out[out$Sex == 'F', "Weight"], pch = 16, col = sex_col[1])
points(out[out$Sex == 'M', "Length"], out[out$Sex == 'M', "Weight"], pch = 16, col = sex_col[2])
lines(lens, ests[paste0("all_F")][[1]][1] * lens ^ ests[paste0("all_F")][[1]][2], col = line_col[1], lwd = 3, lty = 1)
lines(lens, ests[paste0("all_M")][[1]][1] * lens ^ ests[paste0("all_M")][[1]][2], col = line_col[2], lwd = 3, lty = 1)
lines(lens, ests[paste0("all")][[1]][1] * lens ^ ests[paste0("all")][[1]][2], col = 1, lwd = 3, lty = 2)
lines(lens, 1.46e-5 * lens ^ 2.984, col = 'green', lty = 2, lwd = 3)
leg = c(paste0("F: a = ", signif(ests[paste0("all_F")][[1]][1], digits = 3),  
								" b = ", round(ests[paste0("all_F")][[1]][2], 2) ), 
		paste0("M: a = ", signif(ests[paste0("all_M")][[1]][1], digits = 3),  
								" b = ", round(ests[paste0("all_M")][[1]][2], 2) ),
	    paste0("Unsexed: a = ", signif(ests[paste0("all")][[1]][1], digits = 3),  
								" b = ", round(ests[paste0("all")][[1]][2], 2) ),
		"Unsexed: Love 1990, a = 1.45e-05, b = 2.984")
legend("topleft", bty = 'n', legend = leg, lty = c(1, 1, 2, 2), col = c(line_col[1], line_col[2], 1, 'green'), lwd = 3, cex = 1.25)
dev.off()

############################################################################################
#	Comparison between lengths inside vs. outside of the CCA for the HKL
############################################################################################

compare_length_cca (dir = file.path(dir, "data", "biology"), 
					data = survey_dat, file = "hkl_cca_comparison")

############################################################################################
#	Comparison between lengths by sex and across data sources
############################################################################################

length_freq_plot(dir = file.path(dir, "data", "biology"), data = data)

length_by_depth_plot(dir = file.path(dir, "data", "biology"), data = data)

############################################################################################
#	Length or Age frequency plots By Sex, Area, and Fleet. Minimum sampling size for plotting by sex is 10
############################################################################################
data_hist(dir = file.path(dir, "data", "biology", "plots"), 
          data = out, 
          data_type = "Length", 
          group_column = "State", 
          fleet_column = "Source", 
          ymax = c(0.22), 
          do_abline = TRUE)

############################################################################################
# Look at length patterns by latitude - this will be based solely on survey data because
# other data does not have latitude (could add something for "area")
############################################################################################

library(plyr)
out$lat_bin = plyr::round_any(out$Lat, 0.5)
pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Len_by_Latitude.png", w = 7, h = 7, pt = 12)
par(mfrow = c(1, 1), mar = c(4,4,2,2), oma = c(1,1,1,1))
find = which(out$lat_bin < 42)
boxplot(out$Length[find] ~ out$lat_bin[find], ylim = c(0, 40), ylab = "Length (cm)", xlab = "Latitude")
dev.off()


############################################################################################
# Create the length-at-age figure with lit growth estimates added
############################################################################################

line_col = c("red", 'blue')
sex_col = alpha(line_col, 0.5)
xmax = 20
ymax = 30

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Length_Age_by_Sex_Love_Ests.png", w = 7, h = 7, pt = 12)
par(mfrow = c(1, 1))
plot(0:xmax, 25.25 * (1 - exp(-0.18 * (0:xmax + 3.38))), xlab = "Age", ylab = "Length (cm)", las = 1, cex = 1.5,
	cex.lab = 1.5, cex.axis = 1.5, xaxs = "i", yaxs = "i",ylim = c(0, ymax), xlim = c(0, xmax), 
	type = 'l', col = 'red', lwd = 2) 
lines(0:xmax, 24.71 * (1 - exp(-0.06 * (0:xmax + 10.33))), col = 'blue',  lty = 1, lwd = 2)
leg = c("F: Love 1990", "M: Love 1990")
legend("bottomright", bty = 'n', legend = leg, lty = c(1,1), col = c('red', 'blue'), lwd = 2, cex = 1.5)
dev.off()


############################################################################################
#	Quickly look at the commercial samples by gear to see if the amount of data for each
#   and if there looks to be different selectivity by gear type
############################################################################################

# GFS GFT GLN HKL LGL 
#   1  13  11  85  31
#twl = c("GFS", "GFT")
hkl = c("HKL")
lgl = c("LGL")
other = unique(pacfin$GRID)[which(!unique(pacfin$GRID) %in% c(hkl, lgl))]
pacfin$gear[pacfin$GRID %in% hkl] = "hkl"
pacfin$gear[pacfin$GRID %in% lgl] = "lgl"
pacfin$gear[pacfin$GRID %in% other] = "other"

# gear  mean(Length)   N
#  hkl  23.76706      85
#  lgl  22.53548      31
#other  30.87200      25

library(ggplot2)
pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Gear_by_Area.png", w = 7, h = 7, pt = 12)
ggplot(pacfin, aes(Length, fill = gear, color = gear)) +
	geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5) + 
	scale_x_continuous(limits = c(0, 50))
dev.off()


pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Length_Dists_by_Source.png", w = 7, h = 7, pt = 12)
tmp = out[out$Source %in% c("NWFSC_HKL", "NWFSC_WCGBTS", "Triennial"), ]
ggplot(tmp, aes(Length, fill = Source, color = Source)) +
	geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5, bw = 1) + 
	scale_x_continuous(limits = c(0, 50))
dev.off()

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "RecFIN_vs_MRFSS_Length_Dists_by_Area.png", w = 7, h = 7, pt = 12)
tmp = out[out$Source %in% c("RecFIN", "RecFIN_MRFSS", "PacFIN"), ]
tmp$Source[which(tmp$Source %in% c("RecFIN", "RecFIN_MRFSS"))] = "RecFIN"
#tmp = tmp[which(tmp$Length < 37), ]
ggplot(tmp, aes(Length, fill = Source, color = Source)) +
	geom_density(alpha = 0.4, lwd = 0.8, adjust = 0.5, bw = 1)
dev.off()


##################################################################################
# Plot fecundity at length from Dick et al 2017
#################################################################################

len = 0:40
fecund = 3.93184e-7 * len ^ 3.702

pngfun(wd = file.path(dir, "data", "biology", "plots"), file = "Fecundity.png", w = 7, h = 7, pt = 12)
par(mfrow = c(1, 1), oma = c(0, 3, 0, 0))
plot(len, fecund,  xlab = "Length (cm)",
	ylab = "", type = 'l', col = 'red', lwd = 2, las = 1,
	cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, xaxs = 'i', yaxs = 'i')
mtext(side = 2, "Fecundity (millions of eggs)", line = 4, cex = 1.5)
legend('topleft', legend = "Dick et al. 2017", bty = 'n', cex = 1.5)
dev.off()
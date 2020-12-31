################################################################################################
#			PacFIN Comps for the Squarespot Rockfish Assessment 2021	
#		
#					Written by: Chantel Wetzel (11/12/20)
#
#		This code is only used to generate the length comps for squarespot
################################################################################################			

#devtools::install_github("nwfsc-assess/PacFIN.Utilities")
#library(PacFIN.Utilities)
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/PacFIN.Utilities")

setwd("N:/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data")

# Load in the current weight-at-length estimates by sex
# load("C:/Assessments/2019/petrale_2019/Data/Biology/alpha_beta_ests.rda") 
femalea = 1.08e-5; femaleb = 3.09 
malea   = 1.17e-5; maleb = 3.04    
unsexa  = (femalea + malea)/2;  unsexb = (femaleb + maleb)/2        

load(file.path(getwd(), "PacFIN BDS", "PacFIN.SQRS.bds.13.Aug.2020.RData"))
pacfin = PacFIN.SQRS.bds.13.Aug.2020

catch.file = read.csv(file.path(getwd(), "catch_all", "commercial_catches_for_expansion.csv"))
colnames(catch.file) = c("Year", "CA.ALL")

#plotRawData(pacfin)

Pdata = cleanPacFIN(Pdata = pacfin, 
					keep_length_type = c("", "A", "F", "U", "T", NA),
					keep_missing_lengths = FALSE,
					keep_INPFC = c("VUS","CL","VN","COL","NC","SC","EU","CP","EK","MT","PS"))
# Removal Report
# 
# Records in input:                  141 
# Records not in USINPFC             0 
# Records not in INPFC_AREA:         0 
# Records in bad INPFC_AREA:         0 
# Records in badRecords list:        0 
# Records with bad SAMPLE_TYPE       0 
# Records with bad SAMPLE_METHOD     0 
# Records with no SAMPLE_NO          0 
# Records with no usable length      0 
# Records remaining:                 141 

# remove the wild large fish from 1985
remove = which(Pdata$lengthcm > 35)
Pdata = Pdata[-remove,]

MasterPdata = Pdata

#Pdata = getGearGroup(Pdata)

#Pdata$usegear = 'ALL'
#Pdata$mygear = "ALL"
Pdata$fleet = "ALL"
Pdata$stratification = paste(Pdata$state, Pdata$fleet, sep=".")

#################################################################################
# Length comp expansions
#################################################################################

Pdata =  getExpansion_1(Pdata = Pdata, 
						maxExp = 0.95,
						Exp_WA = TRUE, 
						Indiv_Wgts = TRUE,
						plot = FALSE,
						fa = femalea, fb = femaleb, ma = malea, mb = maleb, ua = unsexa, ub = unsexb)

# The convert input will change the catch from external file into pounds
Pdata = getExpansion_2(Pdata = Pdata, 
					   Catch = catch.file, 
					   Units = "MT",
					   maxExp = 0.80)

Pdata$Final_Sample_Size <- capValues(Pdata$Expansion_Factor_1_L * Pdata$Expansion_Factor_2, maxVal = 0.80)

# Look for consistency between lengths and ages of sampled fish
myLbins = c(seq(8, 30, 1))

Lcomps = getComps(Pdata, defaults = c("fishyr", "fleet"), Comps = "LEN")

masterLcomps = Lcomps

Lcomps = doSexRatio(Lcomps, findRatio = TRUE)

writeComps(inComps = Lcomps, 
		   fname = file.path(getwd(), "PacFIN BDS", "forSS", "Lcomps.SQSP.Nov.2020.csv"), 
		   lbins = myLbins, 
		   partition = 2, 
		   sum1 = TRUE,
		   digits = 4)

##############################################################################################################
# Format and rewrite
##############################################################################################################
out = read.csv(file.path(getwd(), "PacFIN BDS", "forSS", "Lcomps.SQSP.Nov.2020.csv"), skip = 3, header = TRUE)
start = which(as.character(out[,1]) %in% c(" Sexes combined ")) + 2
end   = which(as.character(out[,1]) %in% c(" Females then males ")) -1 #nrow(out)
cut_out = out[start:end,]

# format the length comps
cut_out$fleetnum = 1
cut_out$month = 1

ind = which(colnames(cut_out) %in% "L8"):which(colnames(cut_out) %in% "L30.1")

format = cbind(cut_out$fishyr, cut_out$month, cut_out$fleetnum, cut_out$sex, cut_out$partition, 
			   cut_out$Ntows, cut_out$Nsamps, cut_out$InputN, cut_out[,ind])
colnames(format) = c("fishyr", "month", "fleet", "sex", "part", "Nsamps", "Ntows", "InputN", colnames(cut_out[ind]))
format = format[format$fishyr != 2021, ]
write.csv(format, file = paste0(getwd(), "/PacFIN BDS/forSS/Lcomps_unsexed_8_30_formatted.csv"), row.names = FALSE)


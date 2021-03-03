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
fa = 1.08e-5; fb = 3.09 
ma   = 1.17e-5; mb = 3.04    
ua  = (fa + ma)/2;  ub = (fb + mb)/2        

#load(file.path(getwd(), "PacFIN BDS", "PacFIN.SQRS.bds.13.Aug.2020.RData"))
load(file.path(getwd(), "PacFIN BDS", "PacFIN.SQRS.bds.21.Feb.2021.RData"))
bds.file = "SQRS.bds.21.Feb.2021"
pacfin = bds.pacfin

catch.file = read.csv(file.path(getwd(), "catch_all", "commercial_catches_for_expansion.csv"))
colnames(catch.file) = c("Year", "CA.ALL")

Pdata = cleanPacFIN(Pdata = pacfin, 
					CLEAN = TRUE,
					verbose = TRUE)

# remove the wild large fish from 1985
remove = which(Pdata$lengthcm > 35)
Pdata[remove, 'lengthcm'] <- NA

MasterPdata = Pdata

Pdata$fleet = "ALL"
Pdata$stratification = paste(Pdata$state, Pdata$fleet, sep=".")

#################################################################################
# Length comp expansions
#################################################################################

Pdata_exp <- getExpansion_1(Pdata = Pdata,
					   fa = fa, fb = fb, ma = ma, mb = mb, ua = ua, ub = ub)

Pdata_exp <- getExpansion_2(Pdata = Pdata_exp, 
					   Catch = catch.file, 
					   Units = "MT",
					   maxExp = 0.80)

Pdata_exp$Final_Sample_Size <- capValues(Pdata_exp$Expansion_Factor_1_L * Pdata_exp$Expansion_Factor_2, maxVal = 0.80)

# There are only 7 sexed fish - set them to U for simplicity
Pdata_exp[Pdata_exp$SEX != "U", "SEX"] <- "U"

# Look for consistency between lengths and ages of sampled fish
myLbins = c(seq(8, 30, 1))

Lcomps = getComps(Pdata = Pdata_exp, 
				  Comps = "LEN")

writeComps(inComps = Lcomps, 
		   fname = file.path(getwd(), "PacFIN BDS", "forSS", "Lcomps.SQSP.Feb.21.2021.csv"), 
		   lbins = myLbins, 
		   partition = 0, 
		   sum1 = TRUE,
		   digits = 4)

##############################################################################################################
# Format and rewrite
##############################################################################################################
out = read.csv(file.path(getwd(), "PacFIN BDS", "forSS", "Lcomps.SQSP.Feb.21.2021.csv"), skip = 3, header = TRUE)
start = which(as.character(out[,1]) %in% c(" Usexed only ")) + 2
end   = nrow(out) #which(as.character(out[,1]) %in% c(" Females then males ")) -1 #nrow(out)
cut_out = out[start:end,]

cut_out$fleet = 1
ind = which(colnames(cut_out) %in% "L8"):which(colnames(cut_out) %in% "L30.1")

format = cbind(cut_out$fishyr, cut_out$month, cut_out$fleet, cut_out$sex, cut_out$partition, 
			   cut_out$Ntows, cut_out$Nsamps, cut_out$InputN, cut_out[,ind])
colnames(format) = c("fishyr", "month", "fleet", "sex", "part", "Nsamps", "Ntows", "InputN", colnames(cut_out[ind]))
format = format[format$fishyr != 2021, ]
write.csv(format, file = paste0(getwd(), "/PacFIN BDS/forSS/Lcomps_unsexed_8_30_formatted_Feb_2021.csv"), row.names = FALSE)

#########################################################################################
# Calculate the number of trips and fish
#########################################################################################
Pdata = Pdata_exp
temp = Pdata[!is.na(Pdata$lengthcm) & Pdata$SAMPLE_YEAR < 2021,]

Nfish = table(temp$SAMPLE_YEAR, temp$SEX, temp$fleet)

aa = unique(temp$fleet)
yy = sort(unique(temp$SAMPLE_YEAR))
Ntows = matrix(0, length(yy), length(aa))
for(y in 1:length(yy)){
	for(a in 1:length(aa)){
		ind = which(temp$SAMPLE_YEAR == yy[y] & temp$fleet == aa[a])
		if(length(ind) > 0) {Ntows[y, a] = length(unique(temp$SAMPLE_NO[ind])) }
	}
}
colnames(Ntows) = aa
rownames(Ntows) = yy


samps = cbind(rownames(Ntows), Ntows[,"ALL"], Nfish[,,"ALL"])
colnames(samps) = c("Year", "N_Trips", "N_Fish_Unsexed")

write.csv(samps, file = paste0(getwd(), "/PacFIN BDS/forSS/Com_Samples.csv"), row.names = FALSE)

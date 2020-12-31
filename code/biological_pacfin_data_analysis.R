library(ggplot2)
options(stringsAsFactors = TRUE)
library(HandyCode)

dir = "N://Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish/data/PacFIN BDS"
setwd(dir)
load(paste0(dir, "/PacFIN.SQRS.bds.27.Jul.2020.RData"))

data = PacFIN.SQRS.bds.27.Jul.2020
data$length_cm = data$FISH_LENGTH / 10

# Evaluate the available length samples by state
pngfun(wd = getwd(), file = paste0('squarespot_samples_ca_pacfin.png'), h = 12, w = 12)
find = data$SOURCE_AGID == "C"
ggplot(data[find,], aes(x = length_cm)) + 
		geom_histogram() + #facet_grid(~ RECFIN_YEAR) + 
		facet_wrap(facets = c("SAMPLE_YEAR", "SOURCE_AGID"), nrow = 12, ncol = 5) +
		theme_bw() + stat_bin(bins = 60, binwidth = 2)
dev.off()

###########################################################################################################
# Create table of samples by area and year
###########################################################################################################

temp = data[!is.na(data$length_cm) & data$SAMPLE_YEAR < 2021,]

# Once I figure out how to parse by finer area this next line should be replaced
temp$strat = temp$SOURCE_AGID

Nfish = table(temp$SAMPLE_YEAR, temp$strat)

aa = unique(temp$strat)
yy = sort(unique(temp$SAMPLE_YEAR))
Ntows = matrix(0, length(yy), length(aa))
for(y in 1:length(yy)){
	for(a in 1:length(aa)){
		ind = which(temp$SAMPLE_YEAR == yy[y] & temp$strat == aa[a])
		if(length(ind) > 0) {Ntows[y, a] = length(unique(temp$SAMPLE_NO[ind])) }
	}
}
colnames(Ntows) = aa
rownames(Ntows) = yy

samples = NULL
for (a in 1:length(aa)){
	samples = cbind(samples, Ntows[,aa[a]], Nfish[,aa[a]])
}

samples = cbind(Ntows[,"C"], Nfish[,"C"])

colnames(samples) = c("C_NTows", "C_Nfish")

write.csv(samples, file = file.path(getwd(),"PacFIN_Length_Samples.csv"), row.names = TRUE)
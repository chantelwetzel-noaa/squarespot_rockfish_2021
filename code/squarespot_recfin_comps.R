###################################################################################
#
#			Squarespot Rockfish
#   Format recreation sample for use in the model
#
#############################################################################################

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/nwfscSurvey")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/PacFIN.Utilities")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/HandyCode")
devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/dataModerate_2021")

devtools::load_all("C:/Users/Jason.Cope/Documents/Github/nwfscSurvey")
devtools::load_all("C:/Users/Jason.Cope/Documents/Github/PacFIN.Utilities")
devtools::load_all("C:/Users/Jason.Cope/Documents/Github/HandyCode")
devtools::load_all("C:/Users/Jason.Cope/Documents/Github/dataModerate_2021")


library(ggplot2)

dir = "//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Squarespot_Rockfish"

recfin_name = "SQUARESPOT ROCKFISH"
ca_mrfss_code = 8826010148


############################################################################################
#	Load Data
############################################################################################
##RecFIN
#California
ca_recfin = rename_budrick_recfin(read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_rec_lengths_2004_2020_updated.csv", header=T, na.strings = "-"))
ca_recfin = rename_budrick_recfin(read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/Squarespot Revised CRFS Lengths No Region SD501-CALIFORNIA-1980-2020.csv", header=T, na.strings = "-"))
ca_recfin = ca_recfin[ca_recfin$SPECIES_NAME == recfin_name, ]
ca_recfin_data = rename_recfin(data = ca_recfin,
                               area_grouping = list(c("CHANNEL", "SOUTH"), c("BAY AREA", "WINE", "CENTRAL", "REDWOOD", "NOT KNOWN")),
                               area_names = c("south_pt_concep", "north_pt_concep"),
                               #area_column_name = "RECFIN_PORT_NAME",
                               mode_grouping = list(c("BEACH/BANK", "MAN-MADE/JETTY"), c("PARTY/CHARTER BOATS", "PRIVATE/RENTAL BOATS"), "NOT KNOWN"),
                               mode_names = c("shore", "boat", "unknown"),
                               mode_column_name = "RecFIN.Mode.Name" )
ca_recfin_data$Source = "RecFIN_MRFSS" 

##MRFSS
#California
ca_mrfss_full = read.csv("//nwcfile/FRAM/Assessments/CurrentAssessments/DataModerate_2021/Data_From_States/ca/ca_mrfss_bio_1980_2003.csv")
ca_mrfss = ca_mrfss_full[ca_mrfss_full$ST == 6 & ca_mrfss_full$SP_CODE == ca_mrfss_code, ]
ca_mrfss = ca_mrfss[!is.na(ca_mrfss$CNTY), ] # remove records without a county
spc = c(59, 73, 37, 111, 83)
npc = unique(ca_mrfss[!ca_mrfss$CNTY %in% spc, "CNTY"]) 
ca_mrfss$STATE_NAME = "CA"
diff = ca_mrfss$T_LEN / 10 - ca_mrfss$LNGTH / 10
hist(diff)
plot(diff); abline(h = 0, col = 'red', lty = 2)
plot(ca_mrfss$LNGTH, ca_mrfss$T_LEN)
check_1 = ca_mrfss$LNGTH / floor(ca_mrfss$LNGTH) 
check_2 = ca_mrfss$T_LEN / floor(ca_mrfss$T_LEN) 
ca_mrfss_data = rename_mrfss(data = ca_mrfss,
                             len_col = "LNGTH",
                             area_grouping = list(spc, npc), 
                             area_names = c("south_pt_concep", "north_pt_concep"), 
                             area_column_name = "CNTY", 
                             mode_grouping = list(c(1,2), c(6, 7)),
                             mode_names = c("shore", "boat"),
                             mode_column_name = "MODE_FX" )


input = list()
input[[1]] = ca_recfin_data
input[[2]] = ca_mrfss_data

############################################################################################
#	Create data frame with all the input data
############################################################################################
out = create_data_frame(data_list = input)

max = as.numeric(quantile(out$Length, na.rm = TRUE, 0.999))
remove = which(out$Length > max)
out = out[-remove, ]

# Remove the released for the rest of the summaries for now:
out = out[out$Data_Type %in% c("RETAINED", NA), ]

out$Length_cm = out$Length

out$Trawl_id = 1:nrow(out)
n = GetN.fn(dir = file.path(dir, "data", "RecFIN Sample Data"), 
        dat = out, type = "length", species = "others", printfolder = "forSS")
samps = read.csv(file.path(dir, "data", "RecFIN Sample Data", "forSS", "length_SampleSize.csv"))
samps = samps[,-2]
write.csv(samps, file = file.path(dir, "data", "RecFIN Sample Data", "forSS", "length_SampleSize.csv"), row.names = FALSE)

len_bin = seq(8, 30, 1)

lfs = UnexpandedLFs.fn(dir = file.path(dir, "data", "RecFIN Sample Data"), 
                       datL = out, lgthBins = len_bin,
                       sex = 3,  partition = 2, fleet = 2, month = 1)

PlotFreqData.fn(dir = file.path(dir, "data", "RecFIN Sample Data"), 
    dat = lfs$comps_u, ylim=c(0, max(len_bin)), 
    main = "Recreational", yaxs="i", ylab="Length (cm)", dopng = TRUE)

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/nwfscDiag")
#library(nwfscDiag)

#######################################################################################################
# California
#######################################################################################################

mydir = "C:/Assessments/2021/squarespot_rockfish_2021/models"
base_name = "2.5_rec_devs_sigmaR60"

get = get_settings_profile( parameters =  c("NatM_p_1_Fem_GP_1", "SR_BH_steep", "SR_LN(R0)", 
								"L_at_Amax_Fem_GP_1", "VonBert_K_Fem_GP_1"),
							low =  c(0.12, 0.30, -1.2, 20, 0.13),
							high = c(0.30, 1.0,  1.2, 30, 0.23),
							step_size = c(0.005, 0.10, 0.20, 1, 0.01),
							param_space = c('real', 'real', 'relative', 'real', 'real'))

model_settings = get_settings(settings = list(base_name = base_name,
							  run = c("jitter", "profile", "retro"),
							  profile_details = get ))

get = get_settings_profile( parameters =  c("SR_LN(R0)", 
								"L_at_Amax_Fem_GP_1", "VonBert_K_Fem_GP_1"),
							low =  c(-1.2, 20, 0.13),
							high = c(1.2, 30, 0.23),
							step_size = c(0.20, 1, 0.01),
							param_space = c('relative', 'real', 'real'))

model_settings = get_settings(settings = list(base_name = base_name,
							  run = c("jitter", "profile"),
							  profile_details = get ))

run_diagnostics(mydir = mydir, model_settings = model_settings)

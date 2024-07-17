
################################################################################
## TITLE: 2_Imputation.R                                                       #
## AUTHOR: Robin Reichmann                                                     #
## R-Version: version 4.1.1                                                    #
#                                                                              #
## PURPOSE: Imputation of missing values for anthropometry, lifestyle factors, #
#           and biomarkers other than lipocalin 2                              #
#                                                                              #
## OUTPUTS: the original dataset with missing values replaced by random forest #
#           imputed values                                                     #
################################################################################


library("data.table")
library("missForest")


lipo_data <- fread(input = "Data/lipo_data.csv",
                   sep = ",", header = TRUE, data.table = FALSE)

imputation_vars <- c("CACO", "Country", "Center", "Fasting_C",
                     "Phrt_Bld", "Menop_Bld", "Phase_Mnscycle", 
                     "Sex", "Age_Blood",
                     "Height_Adj", "Weight_Adj", "Bmi_Adj",
                     "Waist_Adj", "L_School", "Smoke_Stat",
                     # "Dur_Smok",
                     "Pa_Score", "Pa_Index", "Alc_Drinker", "Alc_Re",
                     "Energy", "Vegs", "Fruits", "Dairy", "Redmeat",
                     "Procmeat", "Fish", "Fibre",
                     "lipocalin2_ngml",
                     "Crp_CLRT_05",
                     # "Igf1_CLRT_02", "Igfbp1_CLRT_02",
                     # "Igfbp2_CLRT_02", "Igfbp3_CLRT_02", "Igfbp3i_CLRT_02",
                     # "Adipo_CLRT_04", "Hmw_Adipo_CLRT_04",
                     "nonHmw_Adipo_CLRT_04",
                     # "Leptin_CLRT_04",
                     # "Sleptin_R_CLRT_04",
                     "Tnfa_CLRT_03",
                     # "Chol_CLRT_08",
                     "Chol_Hdl_CLRT_08", 
                     # "Chol_Ldl_CLRT_08", "Tg_CLRT_08",
                     # "Cpeptide_CLRT_02", "Hba1c_Ngsp_CLRT_01", "Frap_CLRT_10",
                     "Rom_CLRT_10", "Neopterin_CLRT_14"
                     # "Vitd_CLRT_09"
)

sapply(lipo_data[,imputation_vars], function(x) length(which(is.na(x))))


imputed_vars <- c("Height_Adj", "Weight_Adj", "Bmi_Adj", "Waist_Adj",
                  "Alc_Re", "Energy", "Vegs", "Fruits", "Dairy", "Redmeat",
                  "Procmeat", "Fish", "Fibre",
                  "Crp_CLRT_05", "nonHmw_Adipo_CLRT_04", "Tnfa_CLRT_03",
                  "Chol_Hdl_CLRT_08", "Rom_CLRT_10", "Neopterin_CLRT_14")

cat_names <- c("CACO", "Sex", "Country", "Center", "Fasting_C",
               "Phrt_Bld", "Menop_Bld", "Phase_Mnscycle", 
               "L_School", "Smoke_Stat", "Alc_Drinker", "Pa_Index")
lipo_data[, cat_names] <- lapply(lipo_data[, cat_names], function(x) as.factor(x))




# ------------------------------------------------------------------------------


set.seed(101)

# Missing data imputation
pred_RF <- missForest(lipo_data[,imputation_vars],
                      variablewise = TRUE, ntree = 1000, verbose = TRUE, maxiter = 100)
saveRDS(object = pred_RF, file = "Data/missForest_obj.csv")
pred_RF <- readRDS("Data/missForest_obj.csv")


# Calculate imputation error statistics
pred_error <- pred_RF$OOBerror
names(pred_error) <- names(pred_RF$ximp)
pred_error 

lipo_data_imputed <- lipo_data
for(marker_i in imputation_vars) {
  lipo_data_imputed[,marker_i] <- pred_RF$ximp[,marker_i]
}

options(scipen = 99999)
error_table <- cbind("mean" = sapply(lipo_data_imputed[,imputed_vars],
                                     mean),
                     "mean_orig" = sapply(lipo_data[,imputed_vars],
                                          function(x) mean(x, na.rm = T)),
                     "sd" = sapply(lipo_data_imputed[,imputed_vars],
                                   sd),
                     "sd_orig" = sapply(lipo_data[,imputed_vars],
                                        function(x) sd(x, na.rm = T)),
                     "var" = sapply(lipo_data_imputed[,imputed_vars],
                                    var),
                     "var_orig" = sapply(lipo_data[,imputed_vars],
                                         function(x) var(x, na.rm = T)),
                     "IQR" = sapply(lipo_data_imputed[,imputed_vars],
                                    IQR),
                     "IQR_orig" = sapply(lipo_data_imputed[,imputed_vars],
                                         function(x) IQR(x, na.rm = T)),
                     "MSE" = pred_error[imputed_vars],
                     "RMSE" = sqrt(pred_error[imputed_vars]),
                     "NRMSE" = sqrt(pred_error[imputed_vars] / 
                                      sapply(lipo_data_imputed[,imputed_vars],
                                             var)),
                     "NRMSE_m" = sqrt(pred_error[imputed_vars]) / 
                       sapply(lipo_data_imputed[,imputed_vars],
                              mean),
                     "NRMSE_iqr" = sqrt(pred_error[imputed_vars]) / 
                       sapply(lipo_data_imputed[,imputed_vars],
                              IQR),
                     "NRMSE_orig" = sqrt(pred_error[imputed_vars] / 
                                           sapply(lipo_data[,imputed_vars],
                                                  function(x) var(x, na.rm = T))),
                     "NRMSE_m_orig" = sqrt(pred_error[imputed_vars]) / 
                       sapply(lipo_data[,imputed_vars],
                              function(x) mean(x, na.rm = T)),
                     "NRMSE_iqr_orig" = sqrt(pred_error[imputed_vars]) / 
                       sapply(lipo_data[,imputed_vars],
                              function(x) IQR(x, na.rm = T)))
error_table <- as.data.frame(error_table)
error_table
round(error_table, 2)


write.csv(lipo_data_imputed, "Data/lipo_data_RFimputed.csv",
          row.names = FALSE, quote = FALSE, na="")



################################################################################
## TITLE: 4c_IPW proximal colon cancer.R                                       #
## AUTHOR: Robin Reichmann                                                     #
## R-Version: version 4.1.1                                                    #
#                                                                              #
## PURPOSE: Calculation sampling weights using inverse probability weighting   #
#           for the nested case-control dataset for proximal colon cancer      #
#           -> calculate the inverse of the probability that an individual in  #
#              the nested case-control data was selected from the original     #
#              full study data                                                 #
#              - for cases the probability is the proportion of the number     #
#                of selected cases and the number of cases in the full study   #
#                data                                                          #
#              - for controls the probability is based on the number of other  #
#                non-case individuals in the original full study data that     #
#                share the same matching criteria as the case the control was  #
#                matched with                                                  #
#                                                                              #
## OUTPUTS: the original dataset with a column with sampling weights           #
################################################################################


library("data.table")

lipo_data <- read.csv("Data/lipo_CC_prox_data.csv")

lipo_data$D_Recrui <- as.Date(lipo_data$D_Recrui, format = "%d/%m/%Y")
lipo_data$D_Endfup <- as.Date(lipo_data$D_Endfup, format = "%d/%m/%Y")
lipo_data$D_Bld_Coll <- as.Date(lipo_data$D_Bld_Coll, format = "%d/%m/%Y")
lipo_data$T_Bld_Coll <- as.POSIXct(lipo_data$T_Bld_Coll,format="%H:%M:%S")
# lipo_data$D_Recrui <- as.Date(lipo_data$D_Recrui, format = "%Y-%m-%d")
# lipo_data$D_Endfup <- as.Date(lipo_data$D_Endfup, format = "%Y-%m-%d")
# lipo_data$D_Bld_Coll <- as.Date(lipo_data$D_Bld_Coll, format = "%Y-%m-%d")
# lipo_data$T_Bld_Coll <- substr(lipo_data$T_Bld_Coll, 9, 16)
# lipo_data$T_Bld_Coll <- as.POSIXct(lipo_data$T_Bld_Coll,format="%H:%M:%S")


# ------------------------------------------------------------------------------

cohort <- fread(input = "Data/cohort_IPWbase.csv",
                sep = ",", header = TRUE, data.table = FALSE)
# cohort <- cohort[which(cohort$Age_Recr < 70),]

cohort$D_Recrui <- as.Date(cohort$D_Recrui, format = "%d/%m/%Y")
cohort$D_Endfup <- as.Date(cohort$D_Endfup, format = "%d/%m/%Y")
cohort$D_Bld_Coll <- as.Date(cohort$D_Bld_Coll, format = "%d/%m/%Y")
cohort$T_Bld_Coll <- as.POSIXct(cohort$T_Bld_Coll,format="%H:%M:%S")


## generate dataset of unique individuals
#-------------------------------------------------------------------------------
# remove every duplicate entry except the "last" (will be the case, if control turned into case)
# NOT NEEDED IN THIS DATASET
lipo_data_uniq <- lipo_data
# Idepic_Crypt_dupl <- lipo_data_uniq$Idepic_Crypt[which(duplicated(lipo_data_uniq$Idepic_Crypt))]
# for (Idepic_Crypt_i in Idepic_Crypt_dupl) {
#   dupl_idx <- which(lipo_data_uniq$Idepic_Crypt == Idepic_Crypt_i)
#   num_remove <- length(dupl_idx) - 1
#   lipo_data_uniq <- lipo_data_uniq[-dupl_idx[1:num_remove],]
# }


## adjust the time variable for the weighted method:
#-------------------------------------------------------------------------------
case_list <- which(lipo_data_uniq$CC_prox_outcome == 1)
control_list <- which(lipo_data_uniq$CC_prox_outcome == 0)
# cases:    individual date of end of follow up
# controls: latest date of end of follow up of all cases
max_case_D_Endfup <- max(lipo_data_uniq$D_Endfup[case_list])
lipo_data_uniq$D_Endfup_w <- lipo_data_uniq$D_Endfup
lipo_data_uniq$D_Endfup_w[control_list] <- max_case_D_Endfup
lipo_data_uniq$time_w <- (lipo_data_uniq$D_Endfup_w - lipo_data_uniq$D_Bld_Coll) / 365.25


# sampling weights
#-------------------------------------------------------------------------------
# sampling_weights <- c(rep(1, length(lipo_data_uniq$Idepic_Crypt))) # weight 1 for every case
case_prob <- length(which(lipo_data$CC_prox_outcome == 1)) / length(which(cohort$CC_prox_outcome == 1))
case_weight <- 1 / case_prob
sampling_weights <- c(rep(case_weight, length(lipo_data_uniq$Idepic))) # weight for every case

# compute sampling weight for every control
for (control_i in control_list) {
  # show progress (takes a while...)
  cat(paste("\r", which(control_i == control_list)), "of", length(control_list))
  
  # compute weight for each eligible matching
  prob_nincl_list_i <- c()
  for (case_i in case_list) {
    
    # check for eligible matching (fitting matching variables)
    if(lipo_data_uniq$D_Endfup[case_i] <= lipo_data_uniq$D_Endfup[control_i] & # still healthy when matched
       lipo_data_uniq$D_Endfup[case_i] >= lipo_data_uniq$D_Recrui[control_i] & # already in study when case got cancer
       # lipo_data_uniq$FollowUp_Time[case_i] <= lipo_data_uniq$FollowUp_Time[control_i] &
       lipo_data_uniq$Sex[case_i] == lipo_data_uniq$Sex[control_i] &
       lipo_data_uniq$Center[case_i] == lipo_data_uniq$Center[control_i] &
       abs(lipo_data_uniq$Age_Blood[case_i] - lipo_data_uniq$Age_Blood[control_i]) <= 2 & # max 2y difference
       abs(lipo_data_uniq$T_Bld_Coll[case_i] - lipo_data_uniq$T_Bld_Coll[control_i]) <= 14400 & # max 4h difference
       # lipo_data_uniq$Phrt_Bld[case_i] == lipo_data_uniq$Phrt_Bld[control_i] &
       # lipo_data_uniq$Menop_Bld[case_i] == lipo_data_uniq$Menop_Bld[control_i] &
       lipo_data_uniq$Fasting_C[case_i] == lipo_data_uniq$Fasting_C[control_i]
    ) {
      # count number of entries in the full cohort that have similar/same matching
      # characteristics as the selected control (number of "competing" entries)
      M <- length(which(cohort$D_Endfup >= lipo_data_uniq$D_Endfup[case_i] & # still healthy when matched
                          cohort$D_Recrui <= lipo_data_uniq$D_Endfup[case_i] & # already in study when case got cancer
                          # cohort$FollowUp_Time >= lipo_data_uniq$FollowUp_Time[case_i] &
                          cohort$Sex == lipo_data_uniq$Sex[case_i] &
                          cohort$Center == lipo_data_uniq$Center[case_i] &
                          abs(cohort$Age_Blood - lipo_data_uniq$Age_Blood[case_i]) <= 2 & # max 2y difference
                          abs(cohort$T_Bld_Coll - lipo_data_uniq$T_Bld_Coll[case_i]) <= 14400 & # max 4h difference
                          # cohort$Phrt_Bld == lipo_data_uniq$Phrt_Bld[case_i] &
                          # cohort$Menop_Bld == lipo_data_uniq$Menop_Bld[case_i] &
                          cohort$Fasting_C == lipo_data_uniq$Fasting_C[case_i]
      ))
      # timepoint specific prob_nincl (~ chance of not being selected for inclusion)
      prob_nincl <- 1 - (1 / (M - 1))
      prob_nincl_list_i <- c(prob_nincl_list_i, prob_nincl)
    } else { # non eligible matching (different matching variables)
      prob_nincl <- 1
      prob_nincl_list_i <- c(prob_nincl_list_i, prob_nincl)
    }
  }
  # compute sampling prob_incl of the control (~ probability of being selected) over all matchings
  prob_incl_i <- 1 - prod(prob_nincl_list_i)
  sampling_weights[control_i] <-  ifelse(is.infinite(1 / prob_incl_i), 1, 1 / prob_incl_i)
}

# Export weighted dataset for SAS
lipo_data_uniq$sampling_weights <- sampling_weights
write.csv(lipo_data_uniq, "Data/lipo_CC_prox_data_weighted.csv", row.names = FALSE, na="")


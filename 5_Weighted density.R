
################################################################################
## TITLE: 5_Weighted density.R                                                 #
## AUTHOR: Robin Reichmann                                                     #
## R-Version: version 4.1.1                                                    #
#                                                                              #
## PURPOSE: Calculation of density function for the lipocalin 2 distribution   #
#           taking into account the IPW sampling weights for CRC and subsite   #
#           datasets, overall and by sex                                       #
#           -> these may be plotted in the spline panel figure                 #
#                                                                              #
## OUTPUTS: datasets with coordinates for the calculated weighted              #
#           density functions                                                  #
################################################################################


lipo_data <- read.csv("Data/lipo_data_weighted.csv")
lipo_men <- lipo_data[which(lipo_data$Sex == 1),]
lipo_women <- lipo_data[which(lipo_data$Sex == 2),]

lipo_data$w_cut <- ifelse(lipo_data$Sex == 1 & lipo_data$Waist_Adj < 94,
                          0,
                          ifelse(lipo_data$Sex == 2 & lipo_data$Waist_Adj < 80,
                                 0,
                                 1))
lipo_w0 <- lipo_data[which(lipo_data$w_cut == 0),]
lipo_w1 <- lipo_data[which(lipo_data$w_cut == 1),]


density_fun <- function(var, weight, prob.min = 0.03, prob.max = 0.97) {
  q_min <- unname(quantile(var, weights = weight, probs = prob.min))
  q_max <- unname(quantile(var, weights = weight, probs = prob.max))
  idx_within <- which(var >= q_min & var <= q_max)
  dens_obj <- density(var[idx_within],
                      weights = weight[idx_within] / sum(weight[idx_within]))
  return(dens_obj)
}


dens_both <- density_fun(var = lipo_data$lipocalin2_ngml,
                         weight = lipo_data$sampling_weights)
plot(dens_both)
dens_table_both <- cbind(x = dens_both$x, y = dens_both$y)
write.csv(dens_table_both, "Data/Densities/lipo_density_both.csv", row.names = FALSE)


dens_men <- density_fun(var = lipo_men$lipocalin2_ngml,
                        weight = lipo_men$sampling_weights)
plot(dens_men)
dens_table_men <- cbind(x = dens_men$x, y = dens_men$y)
write.csv(dens_table_men, "Data/Densities/lipo_density_men.csv", row.names = FALSE)


dens_women <- density_fun(var = lipo_women$lipocalin2_ngml,
                          weight = lipo_women$sampling_weights)
plot(dens_women)
dens_table_women <- cbind(x = dens_women$x, y = dens_women$y)
write.csv(dens_table_women, "Data/Densities/lipo_density_women.csv", row.names = FALSE)


dens_w0 <- density_fun(var = lipo_w0$lipocalin2_ngml,
                          weight = lipo_w0$sampling_weights)
plot(dens_w0)
dens_table_w0 <- cbind(x = dens_w0$x, y = dens_w0$y)
write.csv(dens_table_w0, "Data/Densities/lipo_density_w0.csv", row.names = FALSE)

dens_w1 <- density_fun(var = lipo_w1$lipocalin2_ngml,
                       weight = lipo_w1$sampling_weights)
plot(dens_w1)
dens_table_w1 <- cbind(x = dens_w1$x, y = dens_w1$y)
write.csv(dens_table_w1, "Data/Densities/lipo_density_w1.csv", row.names = FALSE)


#===============================================================================


lipo_data <- read.csv("Data/lipo_CC_data_weighted.csv")
lipo_men <- lipo_data[which(lipo_data$Sex == 1),]
lipo_women <- lipo_data[which(lipo_data$Sex == 2),]

lipo_data$w_cut <- ifelse(lipo_data$Sex == 1 & lipo_data$Waist_Adj < 94,
                          0,
                          ifelse(lipo_data$Sex == 2 & lipo_data$Waist_Adj < 80,
                                 0,
                                 1))
lipo_w0 <- lipo_data[which(lipo_data$w_cut == 0),]
lipo_w1 <- lipo_data[which(lipo_data$w_cut == 1),]


dens_both <- density_fun(var = lipo_data$lipocalin2_ngml,
                         weight = lipo_data$sampling_weights)
dens_table_both <- cbind(x = dens_both$x, y = dens_both$y)
write.csv(dens_table_both, "Data/Densities/lipo_CC_density_both.csv", row.names = FALSE)


dens_men <- density_fun(var = lipo_men$lipocalin2_ngml,
                        weight = lipo_men$sampling_weights)
dens_table_men <- cbind(x = dens_men$x, y = dens_men$y)
write.csv(dens_table_men, "Data/Densities/lipo_CC_density_men.csv", row.names = FALSE)


dens_women <- density_fun(var = lipo_women$lipocalin2_ngml,
                          weight = lipo_women$sampling_weights)
dens_table_women <- cbind(x = dens_women$x, y = dens_women$y)
write.csv(dens_table_women, "Data/Densities/lipo_CC_density_women.csv", row.names = FALSE)



dens_w0 <- density_fun(var = lipo_w0$lipocalin2_ngml,
                       weight = lipo_w0$sampling_weights)
plot(dens_w0)
dens_table_w0 <- cbind(x = dens_w0$x, y = dens_w0$y)
write.csv(dens_table_w0, "Data/Densities/lipo_CC_density_w0.csv", row.names = FALSE)

dens_w1 <- density_fun(var = lipo_w1$lipocalin2_ngml,
                       weight = lipo_w1$sampling_weights)
plot(dens_w1)
dens_table_w1 <- cbind(x = dens_w1$x, y = dens_w1$y)
write.csv(dens_table_w1, "Data/Densities/lipo_CC_density_w1.csv", row.names = FALSE)


#===============================================================================


lipo_data <- read.csv("Data/lipo_CC_prox_data_weighted.csv")
lipo_men <- lipo_data[which(lipo_data$Sex == 1),]
lipo_women <- lipo_data[which(lipo_data$Sex == 2),]

lipo_data$w_cut <- ifelse(lipo_data$Sex == 1 & lipo_data$Waist_Adj < 94,
                          0,
                          ifelse(lipo_data$Sex == 2 & lipo_data$Waist_Adj < 80,
                                 0,
                                 1))
lipo_w0 <- lipo_data[which(lipo_data$w_cut == 0),]
lipo_w1 <- lipo_data[which(lipo_data$w_cut == 1),]


dens_both <- density_fun(var = lipo_data$lipocalin2_ngml,
                         weight = lipo_data$sampling_weights)
dens_table_both <- cbind(x = dens_both$x, y = dens_both$y)
write.csv(dens_table_both, "Data/Densities/lipo_CC_prox_density_both.csv", row.names = FALSE)


dens_men <- density_fun(var = lipo_men$lipocalin2_ngml,
                        weight = lipo_men$sampling_weights)
dens_table_men <- cbind(x = dens_men$x, y = dens_men$y)
write.csv(dens_table_men, "Data/Densities/lipo_CC_prox_density_men.csv", row.names = FALSE)


dens_women <- density_fun(var = lipo_women$lipocalin2_ngml,
                          weight = lipo_women$sampling_weights)
dens_table_women <- cbind(x = dens_women$x, y = dens_women$y)
write.csv(dens_table_women, "Data/Densities/lipo_CC_prox_density_women.csv", row.names = FALSE)



dens_w0 <- density_fun(var = lipo_w0$lipocalin2_ngml,
                       weight = lipo_w0$sampling_weights)
plot(dens_w0)
dens_table_w0 <- cbind(x = dens_w0$x, y = dens_w0$y)
write.csv(dens_table_w0, "Data/Densities/lipo_CC_prox_density_w0.csv", row.names = FALSE)

dens_w1 <- density_fun(var = lipo_w1$lipocalin2_ngml,
                       weight = lipo_w1$sampling_weights)
plot(dens_w1)
dens_table_w1 <- cbind(x = dens_w1$x, y = dens_w1$y)
write.csv(dens_table_w1, "Data/Densities/lipo_CC_prox_density_w1.csv", row.names = FALSE)


#===============================================================================


lipo_data <- read.csv("Data/lipo_CC_dist_data_weighted.csv")
lipo_men <- lipo_data[which(lipo_data$Sex == 1),]
lipo_women <- lipo_data[which(lipo_data$Sex == 2),]

lipo_data$w_cut <- ifelse(lipo_data$Sex == 1 & lipo_data$Waist_Adj < 94,
                          0,
                          ifelse(lipo_data$Sex == 2 & lipo_data$Waist_Adj < 80,
                                 0,
                                 1))
lipo_w0 <- lipo_data[which(lipo_data$w_cut == 0),]
lipo_w1 <- lipo_data[which(lipo_data$w_cut == 1),]


dens_both <- density_fun(var = lipo_data$lipocalin2_ngml,
                         weight = lipo_data$sampling_weights)
dens_table_both <- cbind(x = dens_both$x, y = dens_both$y)
write.csv(dens_table_both, "Data/Densities/lipo_CC_dist_density_both.csv", row.names = FALSE)


dens_men <- density_fun(var = lipo_men$lipocalin2_ngml,
                        weight = lipo_men$sampling_weights)
dens_table_men <- cbind(x = dens_men$x, y = dens_men$y)
write.csv(dens_table_men, "Data/Densities/lipo_CC_dist_density_men.csv", row.names = FALSE)


dens_women <- density_fun(var = lipo_women$lipocalin2_ngml,
                          weight = lipo_women$sampling_weights)
dens_table_women <- cbind(x = dens_women$x, y = dens_women$y)
write.csv(dens_table_women, "Data/Densities/lipo_CC_dist_density_women.csv", row.names = FALSE)



dens_w0 <- density_fun(var = lipo_w0$lipocalin2_ngml,
                       weight = lipo_w0$sampling_weights)
plot(dens_w0)
dens_table_w0 <- cbind(x = dens_w0$x, y = dens_w0$y)
write.csv(dens_table_w0, "Data/Densities/lipo_CC_dist_density_w0.csv", row.names = FALSE)

dens_w1 <- density_fun(var = lipo_w1$lipocalin2_ngml,
                       weight = lipo_w1$sampling_weights)
plot(dens_w1)
dens_table_w1 <- cbind(x = dens_w1$x, y = dens_w1$y)
write.csv(dens_table_w1, "Data/Densities/lipo_CC_dist_density_w1.csv", row.names = FALSE)


#===============================================================================


lipo_data <- read.csv("Data/lipo_RC_data_weighted.csv")
lipo_men <- lipo_data[which(lipo_data$Sex == 1),]
lipo_women <- lipo_data[which(lipo_data$Sex == 2),]

lipo_data$w_cut <- ifelse(lipo_data$Sex == 1 & lipo_data$Waist_Adj < 94,
                          0,
                          ifelse(lipo_data$Sex == 2 & lipo_data$Waist_Adj < 80,
                                 0,
                                 1))
lipo_w0 <- lipo_data[which(lipo_data$w_cut == 0),]
lipo_w1 <- lipo_data[which(lipo_data$w_cut == 1),]


dens_both <- density_fun(var = lipo_data$lipocalin2_ngml,
                         weight = lipo_data$sampling_weights)
dens_table_both <- cbind(x = dens_both$x, y = dens_both$y)
write.csv(dens_table_both, "Data/Densities/lipo_RC_density_both.csv", row.names = FALSE)


dens_men <- density_fun(var = lipo_men$lipocalin2_ngml,
                        weight = lipo_men$sampling_weights)
dens_table_men <- cbind(x = dens_men$x, y = dens_men$y)
write.csv(dens_table_men, "Data/Densities/lipo_RC_density_men.csv", row.names = FALSE)


dens_women <- density_fun(var = lipo_women$lipocalin2_ngml,
                          weight = lipo_women$sampling_weights)
dens_table_women <- cbind(x = dens_women$x, y = dens_women$y)
write.csv(dens_table_women, "Data/Densities/lipo_RC_density_women.csv", row.names = FALSE)



dens_w0 <- density_fun(var = lipo_w0$lipocalin2_ngml,
                       weight = lipo_w0$sampling_weights)
plot(dens_w0)
dens_table_w0 <- cbind(x = dens_w0$x, y = dens_w0$y)
write.csv(dens_table_w0, "Data/Densities/lipo_RC_density_w0.csv", row.names = FALSE)

dens_w1 <- density_fun(var = lipo_w1$lipocalin2_ngml,
                       weight = lipo_w1$sampling_weights)
plot(dens_w1)
dens_table_w1 <- cbind(x = dens_w1$x, y = dens_w1$y)
write.csv(dens_table_w1, "Data/Densities/lipo_RC_density_w1.csv", row.names = FALSE)
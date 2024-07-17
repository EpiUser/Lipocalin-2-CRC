
################################################################################
## TITLE: Biomarker attenuation plot.R                                         #
## AUTHOR: Robin Reichmann                                                     #
## R-Version: version 4.1.1                                                    #
#                                                                              #
## PURPOSE: Plotting of the hazard ratios with 95% confidence intervals for    #
#           lipocalin 2 in relation to CRC and subsites, additionally          #
#           adjusting for individual and combinations of other biomarkers      #
#                                                                              #
## OUTPUTS: a panel plot illustrating how the association of lipocalin 2 with  #
#           is attentuated after adjustment for other biomarkers               #
################################################################################


library(ggplot2)
library(ggsci)


input_data <- read.csv("Output/biomarker_attenuation.csv")
input_data$highlight <- as.factor(ifelse(input_data$model == "adj", 1, 0))
input_data$outcome <- factor(input_data$outcome,
                             levels = c("CRC", "CC", "CC_prox", "CC_dist", "RC"),
                             labels = c("Colorectal cancer", "Colon cancer",
                                        "Proximal colon cancer",
                                        "Distal colon cancer",
                                        "Rectal cancer"))
input_data$model <- factor(input_data$model,
                           levels = rev(c("adj", "CRP", "nonHMW_adipo",
                                          "TNFalpha", "HDL_chol", "ROM", 
                                          "Neopt", "allBio")),
                           labels = rev(c("Multivariable adjusted model", 
                                          "+ C-reactive protein",
                                          "+ non-HMW adiponectin",
                                          "+ TNF alpha",
                                          "+ HDL cholesterol",
                                          "+ ROM",
                                          "+ neopterin",
                                          "+ all biomarkers")))
input_data$sex <- factor(input_data$sex,
                         levels = c("0", "1", "2"),
                         labels = c("Both sexes", "Men", "Women"))


ggplot(input_data, aes(x = HR, xmin=LCL, xmax=UCL, y = model, color = highlight)) +
  geom_vline(xintercept = 1) + 
  geom_vline(aes(xintercept = ref), lty = 2, color = pal_jama()(2)[2]) +
  geom_point() +
  geom_errorbar() +
  facet_grid(rows = vars(outcome), cols = vars(sex)) +
  scale_x_continuous(trans = "log10", limits = c(0.7, 3.2), breaks = c(0.7, 1, 1.3, 2, 3)) +
  scale_color_jama() +
  labs(x = "Hazard ratio (95% confidence interval)",
       y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold"))
ggsave("Figs/Biomarker_attenuation.png", dpi = 300, height = 8.1, width = 9)


ggplot(input_data, aes(x = HR, xmin=LCL, xmax=UCL, y = model, color = highlight)) +
  geom_vline(xintercept = 1) + 
  geom_vline(aes(xintercept = ref), lty = 2, color = pal_jama()(2)[2]) +
  geom_point() +
  geom_errorbar() +
  facet_grid(rows = vars(outcome), cols = vars(sex), switch = "y") +
  scale_x_continuous(trans = "log10", limits = c(0.7, 3.2), breaks = c(0.7, 1, 1.3, 2, 3)) +
  scale_color_jama() +
  labs(x = "Hazard ratio (95% confidence interval)",
       y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 10),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, units = "inch"))
ggsave("Figs/Biomarker_attenuation2.png", dpi = 300, height = 7, width = 10)

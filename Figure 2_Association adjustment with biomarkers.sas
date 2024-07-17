
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Figure 2_Association adjustment with biomarkers.sas"         *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: The linear association of lipocalin 2 with CRC and subsites is     *
*           evaluated with additional adjustment for biomarkers.               *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  The linear associations of lipocalin 2 with CRC and subsites is calculated  *
*  using weighted Cox proportional hazard regression. Hazard ratios are        *
*  caculated for the crude model, adjusted only for the matching variable,     *
*  for an adjusted model, and also with additionally added individual and      *
*  combined biomarker information. The resulting hazard ratios and confidence  *
*  intervals were plotted in a panel figure using the separate R script        *
*  "Biomarker attenuation plot.R"                                              *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Analysis macros.sas                  * 
*                                                                              *
** -READ: lipo_data_weighted.csv                                               *
*         lipo_CC_data_weighted.csv                                            *
*         lipo_CC_prox_data_weighted.csv                                       *
*         lipo_CC_dist_data_weighted.csv                                       *
*         lipo_RC_data_weighted.csv                                            *
*******************************************************************************/

* CRC data;
PROC Import datafile = "C:\lipo_data_weighted.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;

DATA Lipo_data; set Lipo_data;
  crp_log = log2(Crp_CLRT_05);
  nonHMW_adipo_log = log2(nonHmw_Adipo_CLRT_04);
  tnfa_log = log2(Tnfa_CLRT_03);
  chol_hdl_log = log2(Chol_Hdl_CLRT_08);
  rom_log = log2(Rom_CLRT_10);
  neopterin_log = log2(Neopterin_CLRT_14);
RUN;
DATA Lipo_data_m; set Lipo_data;
  where sex = 1;
RUN;
DATA Lipo_data_w; set Lipo_data;
  where sex = 2;
RUN;

Title "CRC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "CRP";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "TNFalpha";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "HDL chol";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "ROM";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "neopterin";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;




*=============================================================================;


* Colon cancer data;
PROC Import datafile = "C:\lipo_CC_data_weighted.csv" out = Lipo_CC_data dbms = csv replace;
  guessingrows = max;
RUN;

DATA Lipo_CC_data; set Lipo_CC_data;
  crp_log = log2(Crp_CLRT_05);
  nonHMW_adipo_log = log2(nonHmw_Adipo_CLRT_04);
  tnfa_log = log2(Tnfa_CLRT_03);
  chol_hdl_log = log2(Chol_Hdl_CLRT_08);
  rom_log = log2(Rom_CLRT_10);
  neopterin_log = log2(Neopterin_CLRT_14);
RUN;
DATA Lipo_CC_data_m; set Lipo_CC_data;
  where sex = 1;
RUN;
DATA Lipo_CC_data_w; set Lipo_CC_data;
  where sex = 2;
RUN;



Title "CC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "CRP";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "ROM";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;


*=============================================================================;


* Proximal colon cancer data;
PROC Import datafile = "C:\lipo_CC_prox_data_weighted.csv" out = Lipo_CC_prox_data dbms = csv replace;
  guessingrows = max;
RUN;

DATA Lipo_CC_prox_data; set Lipo_CC_prox_data;
  crp_log = log2(Crp_CLRT_05);
  nonHMW_adipo_log = log2(nonHmw_Adipo_CLRT_04);
  tnfa_log = log2(Tnfa_CLRT_03);
  chol_hdl_log = log2(Chol_Hdl_CLRT_08);
  rom_log = log2(Rom_CLRT_10);
  neopterin_log = log2(Neopterin_CLRT_14);
RUN;
DATA Lipo_CC_prox_data_m; set Lipo_CC_prox_data;
  where sex = 1;
RUN;
DATA Lipo_CC_prox_data_w; set Lipo_CC_prox_data;
  where sex = 2;
RUN;


Title "CC_prox";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "CRP";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "ROM";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;


*=============================================================================;


* Distal colon cancer;
PROC Import datafile = "C:\lipo_CC_dist_data_weighted.csv" out = Lipo_CC_dist_data dbms = csv replace;
  guessingrows = max;
RUN;

DATA Lipo_CC_dist_data; set Lipo_CC_dist_data;
  crp_log = log2(Crp_CLRT_05);
  nonHMW_adipo_log = log2(nonHmw_Adipo_CLRT_04);
  tnfa_log = log2(Tnfa_CLRT_03);
  chol_hdl_log = log2(Chol_Hdl_CLRT_08);
  rom_log = log2(Rom_CLRT_10);
  neopterin_log = log2(Neopterin_CLRT_14);
RUN;
DATA Lipo_CC_dist_data_m; set Lipo_CC_dist_data;
  where sex = 1;
RUN;
DATA Lipo_CC_dist_data_w; set Lipo_CC_dist_data;
  where sex = 2;
RUN;


Title "CC_dist";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "CRP";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "ROM";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title;


*=============================================================================;


* Rectal colon cancer data;
PROC Import datafile = "C:\lipo_RC_data_weighted.csv" out = Lipo_RC_data dbms = csv replace;
  guessingrows = max;
RUN;

DATA Lipo_RC_data; set Lipo_RC_data;
  crp_log = log2(Crp_CLRT_05);
  nonHMW_adipo_log = log2(nonHmw_Adipo_CLRT_04);
  tnfa_log = log2(Tnfa_CLRT_03);
  chol_hdl_log = log2(Chol_Hdl_CLRT_08);
  rom_log = log2(Rom_CLRT_10);
  neopterin_log = log2(Neopterin_CLRT_14);
RUN;
DATA Lipo_RC_data_m; set Lipo_RC_data;
  where sex = 1;
RUN;
DATA Lipo_RC_data_w; set Lipo_RC_data;
  where sex = 2;
RUN;


Title "RC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "CRP";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "TNFalpha";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "HDL chol";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "ROM";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "neopterin";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;



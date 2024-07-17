
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Interactions biomarker.sas"                                  *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Evaluation of potential additive and multiplicative interactions   *
*           between lipocalin 2 and other biomarkers in relation to CRC        *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  Potential interactive effects betwen lipocalin 2 and other biomarkers in    *
*  relation to CRC and subsite cancers are evaluated. Multiplicative           *
*  interaction is checked using the Wald test p value of a multiplicative      *
*  interaction term. For additive interaction, the statistics RERI (relative   *
*  excess risk due to interaction), AP (attributable proportion due to         *
*  interaction), and S (synergy index).                                        *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Analysis macros.sas                  *
*                                                                              * 
** -READ: lipo_data_weighted.csv                                               *
*         lipo_CC_data_weighted.csv                                            *
*         lipo_CC_prox_data_weighted.csv                                       *
*         lipo_CC_dist_data_weighted.csv                                       *
*         lipo_RC_data_weighted.csv                                            *
*******************************************************************************/


* Colorectal cancer;

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

* Covariates;
%let covs = Center Age_Blood Sex T_Bld_Coll Fasting_C Phrt_Bld Menop_Bld;
%let covs_cat = Center Sex Fasting_C Phrt_Bld Menop_Bld;
%let covs2 = Waist_Adj Smoke_Stat Pa_Index Alc_Drinker
             Vegs Fruits Redmeat Procmeat Fish Fibre;
%let covs_cat2 = Smoke_Stat Alc_Drinker Pa_Index;
%let biomarkers = crp_log nonHMW_adipo_log tnfa_log chol_hdl_log rom_log neopterin_log;


Title "CRC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);

Title2 "CRP";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=crp_log, weights=sampling_weights);

Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=nonHMW_adipo_log, weights=sampling_weights);

Title2 "TNFalpha";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=tnfa_log, weights=sampling_weights);

Title2 "HDL chol";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=chol_hdl_log, weights=sampling_weights);

Title2 "ROM";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=rom_log, weights=sampling_weights);

Title2 "neopterin";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_data, caseVar=CACO, timeVar=time_w, var1=lipocalin2_log, var2=neopterin_log, weights=sampling_weights);

Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;




*=============================================================================;


* Colon cancer;

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



Title "CC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);

Title2 "CRP";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=crp_log, weights=sampling_weights);

Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=nonHMW_adipo_log, weights=sampling_weights);

Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=tnfa_log, weights=sampling_weights);

Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=chol_hdl_log, weights=sampling_weights);

Title2 "ROM";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=rom_log, weights=sampling_weights);

Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=neopterin_log, weights=sampling_weights);

Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;


*=============================================================================;


* Proximal colon cancer;

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


Title "CC_prox";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);

Title2 "CRP";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=crp_log, weights=sampling_weights);

Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=nonHMW_adipo_log, weights=sampling_weights);

Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=tnfa_log, weights=sampling_weights);

Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=chol_hdl_log, weights=sampling_weights);

Title2 "ROM";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=rom_log, weights=sampling_weights);

Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var1=lipocalin2_log, var2=neopterin_log, weights=sampling_weights);

Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
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


Title "CC_dist";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);

Title2 "CRP";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=crp_log, weights=sampling_weights);

Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=nonHMW_adipo_log, weights=sampling_weights);

Title2 "TNFalpha";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=tnfa_log, weights=sampling_weights);

Title2 "HDL chol";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=chol_hdl_log, weights=sampling_weights);

Title2 "ROM";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=rom_log, weights=sampling_weights);

Title2 "neopterin";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var1=lipocalin2_log, var2=neopterin_log, weights=sampling_weights);

Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;


*=============================================================================;


* Rectal cancer;

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


Title "RC";

Title2 "Both sexes, crude";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          weights=sampling_weights);
Title2 "Both sexes, model 2";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);

Title2 "CRP";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 crp_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=crp_log, weights=sampling_weights);

Title2 "non-HMW adipo";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 nonHMW_adipo_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=nonHMW_adipo_log, weights=sampling_weights);

Title2 "TNFalpha";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 tnfa_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=tnfa_log, weights=sampling_weights);

Title2 "HDL chol";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 chol_hdl_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=chol_hdl_log, weights=sampling_weights);

Title2 "ROM";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 rom_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=rom_log, weights=sampling_weights);

Title2 "neopterin";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 neopterin_log, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
%addInter_cox(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var1=lipocalin2_log, var2=neopterin_log, weights=sampling_weights);

Title2 "+all biomarkers";
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2 &biomarkers, classVars=&covs_cat &covs_cat2,
          weights=sampling_weights);
title;



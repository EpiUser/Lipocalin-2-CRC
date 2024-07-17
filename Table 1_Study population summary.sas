
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Table 1_Study population summary.sas"                        *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Calculate study population characteristics for cases and controls  *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  The average baseline characteristics are calculated for cases and controls. *
*  Also, the significance of differences between cases and controls are tested.*
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Study population summary macros.sas  * 
*                                                                              *
** -READ: lipo_data_RFimputed.csv                                              *
*******************************************************************************/


PROC Import datafile = "C:\Lipo_data_RFimputed.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;

PROC Freq data = Lipo_data;
 table Country Center CACO;
RUN;

PROC Means data = Lipo_data median mean min max std q1 q3;
  var FollowUp_Time lipocalin2_ngml;
RUN;

PROC Means data = Lipo_data median mean min max std q1 q3;
  class sex;
  var lipocalin2_ngml;
RUN;


%let vars_cat = Sex Smoke_Stat Alc_Drinker L_School Pa_Index Menop_Bld Phrt_Bld;
%let vars_num = Age_blood Bmi_Adj Waist_adj Height_Adj Pa_Score
                Alc_Re Wine Beer Energy Vegs Fruits Dairy Redmeat Procmeat Fish Fibre
                Cheese Yogurt;
%let biomarkers = lipocalin2_ngml Crp_CLRT_05 nonHmw_Adipo_CLRT_04 Tnfa_CLRT_03 Chol_Hdl_CLRT_08 Rom_CLRT_10 Neopterin_CLRT_14;


/*%PopulationSummary(Data = Lipo_data, VarsNum = &vars_num,*/
/*                   VarsCat = &vars_cat);*/
%PopulationSummary(Data = Lipo_data(where=(CACO=0)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);
%PopulationSummary(Data = Lipo_data(where=(CACO=1)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);

%PopulationSummary(Data = Lipo_data(where=(CACO=0)), VarsNum = &biomarkers);
%PopulationSummary(Data = Lipo_data(where=(CACO=1)), VarsNum = &biomarkers);






/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "1a_Data preparation.sas"                                     *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Data preparation and export of CRC case-control data with          *
*           lipocalin 2 information                                            *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  The main dataset for analysis is prepared, i.e. renaming variables,         *
*  definition of the CRC outcome, exclusion of incomplete case-control sets    *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Data preparation macros.sas          * 
*                                                                              *
** -READ: data.sas7bdat                                           *
*                                                                              *
** -WRITE: data.csv                                                       *
*******************************************************************************/

libname library 'xxxx';
/*options fmtsearch = (libnorb);*/

DATA norb_data;
  set data;
RUN;

PROC Freq data = norb_data;
  table country;
RUN;

* Delete entries from Norway;
DATA norb_data; set norb_data;
  if Country = "B" then delete;
RUN;

PROC Freq data = norb_data;
  table country;
RUN;


/*PROC Contents data = norb_data order = varnum;*/
/*RUN;*/

*=============================================================================;


* Rename variables;
DATA norb_data;
  set norb_data(rename = (
    Cncr_Caco_Clrt  = CACO
    QE_ENERGY       = Energy
    QGE02           = Vegs
    QGE0401         = Fruits
    QGE05           = Dairy
    QGE0701         = Redmeat
    QGE0704         = Procmeat
    QGE08           = Fish
    QE_Fibt         = Fibre
    Qge1401         = Wine
    Qge0505         = Cheese
    Qge0503         = Yogurt
    Qge140301       = Beer
  ));
RUN;


* Define new variables;
DATA norb_data; set norb_data;
  FollowUp_Time = Agexit - Age_Blood;

  if CACO = 1 then time_clr = 1; * pseudo time variable for conditional logistic regression using PROC Phreg;
  else if CACO = 0 then time_clr = 2;
  else time_clr = .;

  * Define "missing" as a separate value for matching variables;
  if T_Bld_Coll eq . then T_Bld_Coll = 0;

  if Phrt_Bld eq . AND Sex eq 1 then Phrt_Bld = 9;
  else if Phrt_Bld eq . AND Sex eq 2 then Phrt_Bld = 8;

  if Menop_Bld eq . AND Sex eq 1 then Menop_Bld = 9;
  else if Menop_Bld eq . AND Sex eq 2 then Menop_Bld = 8;

  if Phase_Mnscycle eq . AND Sex eq 1 then Phase_Mnscycle = 9;
  else if Phase_Mnscycle eq . AND Sex eq 2 then Phase_Mnscycle = 8;

  if Fasting_C eq . then Fasting_C = 9;

  * Non-HMW adiponectin;
  nonHmw_Adipo_CLRT_04 = Adipo_CLRT_04 - Hmw_Adipo_CLRT_04;

  * Waist circumference categories;
  if Sex = 1 and Waist_Adj >= 94 then Waist_cut = 1;
  else if Sex = 2 and Waist_adj >= 80 then Waist_cut = 1;
  else Waist_cut = 0;
RUN;


*=============================================================================;


* Exclude participants with missing lipocalin 2;
%CountMissings(Data=norb_data, VarNames=lipocalin2_ngml);
DATA norb_data; 
  set norb_data;
  if lipocalin2_ngml eq . then delete;
RUN;
PROC Freq data = norb_data;
  table CACO;
RUN;

* Handle measurements of lipocalin2 that are below the limit of detection;
DATA norb_data;
  set norb_data;
  if lipocalin2_ngml < 0.02 then lipocalin2_ngml = 0.01; * 1/2 LOD;
  lipocalin2_log = log2(lipocalin2_ngml);
RUN;



* Remove casesets without case;
%CollectMatches(Data=norb_data, OutVar=CACO, MatchVar=Match_Caseset);
PROC Freq data = norb_data;
  table CACO;
RUN;


* Remove remaining unmatched entries;
%CleanIncompleteMatches(Data=norb_data, MatchVar=Match_Caseset);
PROC Freq data = norb_data;
  table CACO;
RUN;


* Remove additional controls;
%RemoveAdditionalControls(Data=norb_data, MatchVar=Match_Caseset, CntlNumVar=Match_Ctrlnum);
PROC Freq data = norb_data;
  table CACO;
RUN;


*=============================================================================;


%let variables_of_interest = CACO FollowUp_Time time_clr
                             /*Country Center*/
                             Fasting_C Phrt_Bld Menop_Bld Phase_Mnscycle 
                             Sex Age_Blood 
                             Height_Adj Weight_Adj Bmi_Adj
                             Waist_Adj L_School Smoke_Stat /*Dur_Smok*/
                             Pa_Score Pa_Index Alc_Drinker Alc_Re
                             Energy Vegs Fruits Dairy Redmeat Procmeat Fish Fibre
                             Wine Cheese Yogurt Beer
/*                             Systol_1 Diastol_1 Systol_2 Diastol_2 */
                             lipocalin2_ngml lipocalin2_log
                             Crp_CLRT_05 Igf1_CLRT_02 Igfbp1_CLRT_02    
                             Igfbp2_CLRT_02 Igfbp3_CLRT_02 Igfbp3i_CLRT_02   
                             Adipo_CLRT_04 Hmw_Adipo_CLRT_04 nonHmw_Adipo_CLRT_04
                             Leptin_CLRT_04 
                             Sleptin_R_CLRT_04 Tnfa_CLRT_03 Chol_CLRT_08
                             Chol_Hdl_CLRT_08 Chol_Ldl_CLRT_08 Tg_CLRT_08
                             Cpeptide_CLRT_02 Hba1c_Ngsp_CLRT_01 Frap_CLRT_10
                             Rom_CLRT_10 Neopterin_CLRT_14 Vitd_CLRT_09;

* Count missing entries according to case status;
%CountMissings(data = Norb_data, varNames = &variables_of_interest);
%CountMissings(data = Norb_data(where=(CACO=0)), varNames = &variables_of_interest);
%CountMissings(data = Norb_data(where=(CACO=1)), varNames = &variables_of_interest);

%CountMissings(data = Norb_data, varNames = Country Center, VarFormat="char");
%CountMissings(data = Norb_data(where=(CACO=0)), varNames = Country Center, VarFormat="char");
%CountMissings(data = Norb_data(where=(CACO=1)), varNames = Country Center, VarFormat="char");


PROC Freq data = Norb_data;
  table Duke_Stat Gradclrt Bdg1clrt Stagclrt;
RUN;


*=============================================================================;


* Export data for missForest imputation in R;
Option missing = "";
PROC Export data = Norb_data outfile = "filepath/data.csv" DBMS = csv replace;
RUN;


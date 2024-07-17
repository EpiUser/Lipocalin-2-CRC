
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "1b_Data preparation for IPW.sas"                             *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Data preparation of full cohort data for                           *
*           inverse probability weighting                                      *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  Data of the full EPIC cohort study is prepared in line with the matched     *
*  case-control data to calculate sampling weights for with                    *
*  inverse-probability weighting                                               *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Data preparation macros.sas          * 
*                                                                              *
** -READ: data.sas7bdat                                   *
*                                                                              *
** -WRITE: data.csv                                                  *
*******************************************************************************/

libname lib_caco 'C:\....';
/*options fmtsearch = (lib_caco);*/

DATA full_cohort;
  set data;
RUN;

PROC Freq data = full_cohort;
  table country;
RUN;

* Delete entries from Norway;
DATA full_cohort; set full_cohort;
  if Country = "B" then delete;
RUN;


PROC Freq data = full_cohort;
  table country;
RUN;

TITLE "Full";
PROC Freq data = full_cohort;
  table Cncr_Mal_clrt;
RUN;
TITLE;


* Define colorectal cancer outcome variables;
DATA full_cohort; set full_cohort;
  if Case_Mal_Colon_Dist = 1 or Case_Mal_Colon_Nos = 1 or Case_Mal_Colon_Prox = 1 then
  CC_outcome = 1;
  else CC_outcome = 0;

  if Case_Mal_Colon_Prox = 1 then
  CC_prox_outcome = 1;
  else CC_prox_outcome = 0;

  if Case_Mal_Colon_Dist = 1 then
  CC_dist_outcome = 1;
  else CC_dist_outcome = 0;

  if Case_Mal_Rectum = 1 then
  RC_outcome = 1;
  else RC_outcome = 0;
RUN;

PROC Freq data = full_cohort;
  table CC_outcome CC_prox_outcome CC_dist_outcome RC_outcome;
RUN;




DATA full_cohort; set full_cohort; 
  if Age_Blood eq . then delete; * Only keep entries with collected blood;
RUN;

TITLE "With Blood";
PROC Freq data = full_cohort;
  table Cncr_Mal_clrt;
RUN;
TITLE;


* Remove entries with missing fundamental data;
%RemoveMissingList(Data = full_cohort,
                   VarList = Cncr_Mal_Clrt Age_Recr Agexit D_Recrui D_Endfup Age_blood D_Bld_Coll);

TITLE "Complete";
PROC Freq data = full_cohort;
  table Cncr_Mal_clrt;
RUN;
TITLE;

* Combined family history variables;
DATA full_cohort; set full_cohort;
  if Bs_Clrt = 1 or Fat_Clrt = 1 or Mot_Clrt = 1 then Fam_Hist = 1;
  else if Bs_Clrt = . or Fat_Clrt = . or Mot_Clrt = . then Fam_Hist = .;
  else Fam_Hist = 0;

  if Fat_Clrt = 1 or Mot_Clrt = 1 then Fam_Hist_p = 1;
  else if Fat_Clrt = . or Mot_Clrt = . then Fam_Hist_p = .;
  else Fam_Hist_p = 0;
RUN;

* Replace missing values within matching variables;
* (since matching include missing matching variables as a matching characteristic);
DATA full_cohort; set full_cohort;
  if Sex eq . then Sex = 3;
  if T_Bld_Coll eq . then T_Bld_Coll = 0;
  if Phrt_Bld eq . AND Sex eq 1 then Phrt_Bld = 9;
  else if Phrt_Bld eq . AND Sex eq 2 then Phrt_Bld = 8;
  if Menop_Bld eq . AND Sex eq 1 then Menop_Bld = 9;
  else if Menop_Bld eq . AND Sex eq 2 then Menop_Bld = 8;
  if Phase_Mnscycle eq . AND Sex eq 1 then Phase_Mnscycle = 9;
  else if Phase_Mnscycle eq . AND Sex eq 2 then Phase_Mnscycle = 8;
  if Fasting_c eq . then Fasting_c = 9;
RUN;

* Export for IPW analysis in R;
DATA full_cohort_export; set full_cohort;
  keep IdEpic_Crypt Cncr_Mal_Clrt CC_outcome CC_prox_outcome CC_dist_outcome RC_outcome 
       Sex Center Country Age_Recr Agexit D_Recrui D_Endfup
       Age_Blood D_Bld_Coll T_Bld_Coll Phrt_Bld Menop_Bld Phase_Mnscycle Fasting_c 
       Fam_Hist Fam_Hist_p Use_Nsaid;
RUN;

Option missing = ""; * no value instead of . for missing values;
PROC Export data = full_cohort_export outfile = "C:/..../data.csv" DBMS = csv replace;
RUN;

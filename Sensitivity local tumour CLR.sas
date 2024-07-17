
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Sensitivity local tumour.sas"                                *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-20-07                                            *
*                                                                              *
** PURPOSE: Sensitivity ananlysis only including cases and matched controls    *
*  with a local colorectal or subsite cancer tumour.                           *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  The linear associations of lipocalin 2 with CRC and subsites is calculated  *
*  using conditional logistic regression after excluding all cases *
*  and matched controls with non-local tumours. Odds ratios are caculated    *
*  for an adjusted model, including all matching factors and additional        *
*  anthropometry and lifestyle factors.                                        *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Analysis macros.sas                  *
*                                                                              *
** -READ: lipo_data_weighted.csv                                               *
*         lipo_CC_data_weighted.csv                                            *
*         lipo_CC_prox_data_weighted.csv                                       *
*         lipo_CC_dist_data_weighted.csv                                       *
*         lipo_RC_data_weighted.csv                                            *
*******************************************************************************/

* Covariates for adjustment;
%let covs = Center Age_Blood Sex T_Bld_Coll Fasting_C Phrt_Bld Menop_Bld;
%let covs_cat = Center Sex Fasting_C Phrt_Bld Menop_Bld;
%let covs2 = Waist_Adj Smoke_Stat Pa_Index Alc_Drinker
             Vegs Fruits Redmeat Procmeat Fish Fibre;
%let covs2_w = Smoke_Stat Pa_Index Alc_Drinker
               Vegs Fruits Redmeat Procmeat Fish Fibre;
%let covs2_cat = Smoke_Stat Alc_Drinker Pa_Index;



* CRC data;
PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\lipo_data_weighted.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_data; set Lipo_data;
  if CACO = 1 and Stagclrt ^= 2 then delete;
RUN;
%CollectMatches(Data=Lipo_data, OutVar=CACO, MatchVar=Match_Caseset);
%CleanIncompleteMatches(Data=Lipo_data, MatchVar=Match_Caseset);
DATA Lipo_data_m; set Lipo_data;
  where sex = 1;
RUN;
DATA Lipo_data_w; set Lipo_data;
  where sex = 2;
RUN;

TITLE "Colorectal cancer";
PROC Freq data = Lipo_data;
  table Waist_cut*CACO / nocol nocum nopercent norow;
RUN; 
TITLE "Both sexes";
%OR_Logistic(data=Lipo_data, caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data(where=(Waist_cut = 0)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data(where=(Waist_cut = 1)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Men";
%OR_Logistic(data=Lipo_data_m, caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data_m(where=(Waist_cut = 0)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data_m(where=(Waist_cut = 1)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Women";
%OR_Logistic(data=Lipo_data_w, caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data_w(where=(Waist_cut = 0)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_data_w(where=(Waist_cut = 1)), caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);



* CC data; 
PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\lipo_CC_data_weighted.csv" out = Lipo_CC_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_data; set Lipo_CC_data;
  if CACO = 1 and Stagclrt ^= 2 then delete;
RUN;
%CollectMatches(Data=Lipo_CC_data, OutVar=CC_outcome, MatchVar=Match_Caseset);
%CleanIncompleteMatches(Data=Lipo_CC_data, MatchVar=Match_Caseset);
DATA Lipo_CC_data_m; set Lipo_CC_data;
  where sex = 1;
RUN;
DATA Lipo_CC_data_w; set Lipo_CC_data;
  where sex = 2;
RUN;

TITLE "Colon cancer";
PROC Freq data = Lipo_CC_data;
  table Waist_cut*CC_outcome / nocol nocum nopercent norow;
RUN; 
TITLE "Both sexes";
%OR_Logistic(data=Lipo_CC_data, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data(where=(Waist_cut = 0)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data(where=(Waist_cut = 1)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Men";
%OR_Logistic(data=Lipo_CC_data_m, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data_m(where=(Waist_cut = 0)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data_m(where=(Waist_cut = 1)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Women";
%OR_Logistic(data=Lipo_CC_data_w, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data_w(where=(Waist_cut = 0)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_data_w(where=(Waist_cut = 1)), caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


* CC proximal data; 
PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\lipo_CC_prox_data_weighted.csv" out = Lipo_CC_prox_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_prox_data; set Lipo_CC_prox_data;
  if CACO = 1 and Stagclrt ^= 2 then delete;
RUN;
%CollectMatches(Data=Lipo_CC_prox_data, OutVar=CC_prox_outcome, MatchVar=Match_Caseset);
%CleanIncompleteMatches(Data=Lipo_CC_prox_data, MatchVar=Match_Caseset);
DATA Lipo_CC_prox_data_m; set Lipo_CC_prox_data;
  where sex = 1;
RUN;
DATA Lipo_CC_prox_data_w; set Lipo_CC_prox_data;
  where sex = 2;
RUN;

TITLE "Proximal colon cancer";
PROC Freq data = Lipo_CC_prox_data;
  table Waist_cut*CC_prox_outcome / nocol nocum nopercent norow;
RUN; 
TITLE "Both sexes";
%OR_Logistic(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Men";
%OR_Logistic(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data_m(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data_m(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Women";
%OR_Logistic(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data_w(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_prox_data_w(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


* CC distal data; 
PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\lipo_CC_dist_data_weighted.csv" out = Lipo_CC_dist_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_dist_data; set Lipo_CC_dist_data;
  if CACO = 1 and Stagclrt ^= 2 then delete;
RUN;
%CollectMatches(Data=Lipo_CC_dist_data, OutVar=CC_dist_outcome, MatchVar=Match_Caseset);
%CleanIncompleteMatches(Data=Lipo_CC_dist_data, MatchVar=Match_Caseset);
DATA Lipo_CC_dist_data_m; set Lipo_CC_dist_data;
  where sex = 1;
RUN;
DATA Lipo_CC_dist_data_w; set Lipo_CC_dist_data;
  where sex = 2;
RUN;

TITLE "Distal colon cancer";
PROC Freq data = Lipo_CC_dist_data;
  table Waist_cut*CC_dist_outcome / nocol nocum nopercent norow;
RUN; 
TITLE "Both sexes";
%OR_Logistic(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Men";
%OR_Logistic(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data_m(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data_m(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Women";
%OR_Logistic(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data_w(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_CC_dist_data_w(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);



* RC data;
PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\lipo_RC_data_weighted.csv" out = Lipo_RC_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_RC_data; set Lipo_RC_data;
  if CACO = 1 and Stagclrt ^= 2 then delete;
RUN;
%CollectMatches(Data=Lipo_RC_data, OutVar=RC_outcome, MatchVar=Match_Caseset);
%CleanIncompleteMatches(Data=Lipo_RC_data, MatchVar=Match_Caseset);
DATA Lipo_RC_data_m; set Lipo_RC_data;
  where sex = 1;
RUN;
DATA Lipo_RC_data_w; set Lipo_RC_data;
  where sex = 2;
RUN;

TITLE "Rectal cancer";
PROC Freq data = Lipo_RC_data;
  table Waist_cut*RC_outcome / nocol nocum nopercent norow;
RUN; 
TITLE "Both sexes";
%OR_Logistic(data=Lipo_RC_data, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data(where=(Waist_cut = 0)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data(where=(Waist_cut = 1)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Men";
%OR_Logistic(data=Lipo_RC_data_m, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data_m(where=(Waist_cut = 0)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data_m(where=(Waist_cut = 1)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
TITLE "Women";
%OR_Logistic(data=Lipo_RC_data_w, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data_w(where=(Waist_cut = 0)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);
%OR_Logistic_subset(data=Lipo_RC_data_w(where=(Waist_cut = 1)), caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);



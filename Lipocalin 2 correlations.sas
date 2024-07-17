
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Lipocalin 2 correlations.sas"                                *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Check correlations between lipocalin 2 and other covariates and    *
*           biomarkers before and after imputation.                            *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  Pairwise Spearman correlations partially adjusted for Age and Sex are       *
*  calculated for the dataset before and the dataset after the imputation. The *
*  resulting Spearman correlation coefficients were plotted using the separate *
*  R script "Biomarker correlation.R"                                          *
*                                                                              *
** -READ: lipo_data.csv                                                        *
*         lipo_data_RFimputed.csv                                              *
*******************************************************************************/


* Calculate partially adjusted pairwise Spearman correlations for the data including missings;
PROC Import datafile = "C:\Lipo_data.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;


DATA Lipo_controls; set Lipo_data;
  where CACO = 0;
RUN;

* Macro to calculate pairwise correlations using pairwise complete data;
%macro SpearmanList(data=,var=,withVarList=,partial=);
  %let NumList = %sysfunc(countw(&withVarList));
  %do i = 1 %to &NumList;
    %let withVar_i = %scan(&withVarList, &i);
    PROC Corr data = &data spearman nosimple;
      var &var;
      with &withVar_i;
      partial &partial;
    RUN;
  %end;
%mend;
%SpearmanList(data=Lipo_controls,
              var=lipocalin2_ngml,
              withVarList=Crp_CLRT_05 Igf1_CLRT_02 Igfbp1_CLRT_02    
                          Igfbp2_CLRT_02 Igfbp3_CLRT_02 Igfbp3i_CLRT_02   
                          Adipo_CLRT_04 Hmw_Adipo_CLRT_04 nonHmw_Adipo_CLRT_04
                          Leptin_CLRT_04 
                          Sleptin_R_CLRT_04 Tnfa_CLRT_03 Chol_CLRT_08
                          Chol_Hdl_CLRT_08 Chol_Ldl_CLRT_08 Tg_CLRT_08
                          Cpeptide_CLRT_02 Hba1c_Ngsp_CLRT_01 Frap_CLRT_10
                          Rom_CLRT_10 Neopterin_CLRT_14 Vitd_CLRT_09,
              partial=Age_blood Sex);


*=============================================================================;


* Calculate partially adjusted pariwise Spearman correlations for the imputed data;
PROC Import datafile = "C:\lipo_data_RFimputed.csv" out = Lipo_data_imp dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_controls_imp; set Lipo_data_imp;
  where CACO = 0;
RUN;

PROC Corr data = Lipo_controls_imp spearman fisher;
  var lipocalin2_ngml;
  with Age_blood;
  partial Sex;
RUN;
PROC Corr data = Lipo_controls_imp spearman fisher;
  var lipocalin2_ngml;
  with Bmi_Adj Waist_adj Crp_CLRT_05 nonHmw_Adipo_CLRT_04 Tnfa_CLRT_03
       Chol_Hdl_CLRT_08 Rom_CLRT_10 Neopterin_CLRT_14;
  partial Age_blood Sex;
RUN;



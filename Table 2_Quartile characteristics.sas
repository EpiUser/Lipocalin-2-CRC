
/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Table 2_Quartile characteristics.sas"                        *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-01-05                                            *
*                                                                              *
** PURPOSE: Calculate study population characteristics of controls by          *
*           quartiles of lipocalin 2                                           *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  The average baseline characteristics are calculated for controls by         *
*  quartiles of lipocalin 2. Also, the significance of trend is tested.        *
*                                                                              *
** -PREREQUISITES: macros from SAS macros/Study population summary macros.sas  * 
*                                                                              *
** -READ: lipo_data_weighted.csv                                               *
*******************************************************************************/


PROC Import datafile = "H:\Projects\RR4_Lipocalin 2 and CRC\Data\Lipo_data_weighted.csv" out = Lipo_data dbms = csv replace;
/*  guessingrows = max;*/
RUN;

/*%controlSexQuartiles(data=Lipo_data, var=lipocalin2_ngml, caseVar=CACO, outVar=lipocalin2_q);*/
%controlQuartiles(data=Lipo_data, var=lipocalin2_ngml, caseVar=CACO, outVar=lipocalin2_q);

PROC Freq data = Lipo_data;
  table lipocalin2_q*CACO / nocol nopercent norow;
RUN;

PROC Means data = Lipo_data median min max;
 where CACO = 0;
 class lipocalin2_q;
 var lipocalin2_ngml;
RUN;

PROC Means data = Lipo_data q1 median q3;
 where CACO = 0;
 class sex;
 var lipocalin2_ngml;
RUN;


*=============================================================================;


DATA Lipo_controls; set Lipo_data;
  where CACO = 0;
RUN;

PROC Freq data = Lipo_controls;
  table lipocalin2_q*CACO / nocol nopercent norow;
RUN;


%let vars_cat = Sex Smoke_Stat /*Alc_Drinker*/ L_School Pa_Index Menop_Bld Phrt_Bld;
%let vars_num = Age_blood Bmi_Adj Waist_adj /*Height_Adj*/ /*Pa_Score*/
                Alc_Re /*Wine Beer*/ Vegs Fruits /*Dairy Cheese Yogurt*/ Redmeat Procmeat Fish Fibre;


%PopulationSummary(Data = Lipo_controls(where=(lipocalin2_q=1)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);
%PopulationSummary(Data = Lipo_controls(where=(lipocalin2_q=2)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);
%PopulationSummary(Data = Lipo_controls(where=(lipocalin2_q=3)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);
%PopulationSummary(Data = Lipo_controls(where=(lipocalin2_q=4)), VarsNum = &vars_num,
                   VarsCat = &vars_cat);


*=============================================================================;


%JT_trendTest(Data=Lipo_controls, numVars=&vars_num, trendVar=lipocalin2_q);


/*%CMH_trendTest(Data=Lipo_controls, CatVars=&vars_cat, trendVar=lipocalin2_q);*/


%DummyVars(data=Lipo_controls, varList=&vars_cat);
%let vars_cat_dummies = Sex1 Sex2 
                        Smoke_Stat1 Smoke_Stat2 Smoke_Stat3 Smoke_Stat4
                        L_School0 L_School1 L_School2 L_School3 L_School4 L_School5
                        Pa_Index1 Pa_Index2 Pa_Index3 Pa_Index4 Pa_Index5
                        Menop_Bld1
                        Phrt_Bld1;
%CA_trendTest(Data=Lipo_controls, CatVars=&vars_cat_dummies, trendVar=lipocalin2_q);
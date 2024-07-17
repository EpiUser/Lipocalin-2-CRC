

%macro controlSexQuintiles(data=, var=, caseVar=, sexVar=sex, outVar=);
  PROC Univariate data = &data(where = (&caseVar=0)) noprint;
    where &sexVar = 1;
    var &var;
    output out = percentiles_m pctlpre=P_ pctlpts=20 to 80 by 20;
  RUN;
  PROC Univariate data = &data(where = (&caseVar=0)) noprint;
    where &sexVar ^= 1;
    var &var;
    output out = percentiles_w pctlpre=P_ pctlpts=20 to 80 by 20;
  RUN;
  DATA _null_;
    set percentiles_m;
    call symputx("p_20_m", p_20);
    call symputx("p_40_m", p_40);
    call symputx("p_60_m", p_60);
    call symputx("p_80_m", p_80);
  RUN;
  %put Quintiles men: &p_20_m &p_40_m &p_60_m &p_80_m;
  DATA _null_;
    set percentiles_w;
    call symputx("p_20_w", p_20);
    call symputx("p_40_w", p_40);
    call symputx("p_60_w", p_60);
    call symputx("p_80_w", p_80);
  RUN;
  %put Quintiles women: &p_20_w &p_40_w &p_60_w &p_80_w;

  DATA &data; set &data;
    if &sexVar = 1 then do;
      if &var < &p_20_m then &outVar = 1;
      else if &var < &p_40_m then &outVar = 2;
      else if &var < &p_60_m then &outVar = 3;
      else if &var < &p_80_m then &outVar = 4;
      else &outVar = 5;
    end;
    else if &sexVar ^= 1 then do;
      if &var < &p_20_w then &outVar = 1;
      else if &var < &p_40_w then &outVar = 2;
      else if &var < &p_60_w then &outVar = 3;
      else if &var < &p_80_w then &outVar = 4;
      else &outVar = 5;
    end;
  RUN;

  PROC Means data = &data median noprint;
    class &outVar &sexVar;
    var &var;
    output out = Out_medians median = median;
  RUN;
  DATA _NULL_; set Out_medians;
    if &outVar = 1 and &sexVar = 1 then call symputx("median_1_m", median);
    else if &outVar = 2 and &sexVar = 1 then call symputx("median_2_m", median);
    else if &outVar = 3 and &sexVar = 1 then call symputx("median_3_m", median);
    else if &outVar = 4 and &sexVar = 1 then call symputx("median_4_m", median);
    else if &outVar = 5 and &sexVar = 1 then call symputx("median_5_m", median);
    else if &outVar = 1 and &sexVar ^= 1 then call symputx("median_1_w", median);
    else if &outVar = 2 and &sexVar ^= 1 then call symputx("median_2_w", median);
    else if &outVar = 3 and &sexVar ^= 1 then call symputx("median_3_w", median);
    else if &outVar = 4 and &sexVar ^= 1 then call symputx("median_4_w", median);
    else if &outVar = 5 and &sexVar ^= 1 then call symputx("median_5_w", median);
  RUN;
  DATA &data; set &data;
    if &outVar = 1 and &sexVar = 1 then &outVar.m = &median_1_m;
    else if &outVar = 2 and &sexVar = 1 then &outVar.m = &median_2_m;
    else if &outVar = 3 and &sexVar = 1 then &outVar.m = &median_3_m;
    else if &outVar = 4 and &sexVar = 1 then &outVar.m = &median_4_m;
    else if &outVar = 5 and &sexVar = 1 then &outVar.m = &median_5_m;
    else if &outVar = 1 and &sexVar ^= 1 then &outVar.m = &median_1_w;
    else if &outVar = 2 and &sexVar ^= 1 then &outVar.m = &median_2_w;
    else if &outVar = 3 and &sexVar ^= 1 then &outVar.m = &median_3_w;
    else if &outVar = 4 and &sexVar ^= 1 then &outVar.m = &median_4_w;
    else if &outVar = 5 and &sexVar ^= 1 then &outVar.m = &median_5_w;
  RUN;

  PROC Delete data = percentiles_m; RUN;
  PROC Delete data = percentiles_w; RUN;
  PROC Delete data = Out_medians; RUN;
%mend;


%macro controlSexQuartiles(data=, var=, caseVar=, sexVar=sex, outVar=);
  PROC Univariate data = &data(where = (&caseVar=0)) noprint;
    where &sexVar = 1;
    var &var;
    output out = percentiles_m pctlpre=P_ pctlpts=25 to 75 by 25;
  RUN;
  PROC Univariate data = &data(where = (&caseVar=0)) noprint;
    where &sexVar ^= 1;
    var &var;
    output out = percentiles_w pctlpre=P_ pctlpts=25 to 75 by 25;
  RUN;
  DATA _null_;
    set percentiles_m;
    call symputx("p_25_m", p_25);
    call symputx("p_50_m", p_50);
    call symputx("p_75_m", p_75);
  RUN;
  %put Quartiles men: &p_25_m &p_50_m &p_75_m;
  DATA _null_;
    set percentiles_w;
    call symputx("p_25_w", p_25);
    call symputx("p_50_w", p_50);
    call symputx("p_75_w", p_75);
  RUN;
  %put Quartiles women: &p_25_w &p_50_w &p_75_w;

  DATA &data; set &data;
    if &sexVar = 1 then do;
      if &var < &p_25_m then &outVar = 1;
      else if &var < &p_50_m then &outVar = 2;
      else if &var < &p_75_m then &outVar = 3;
      else &outVar = 4;
    end;
    else if &sexVar ^= 1 then do;
      if &var < &p_25_w then &outVar = 1;
      else if &var < &p_50_w then &outVar = 2;
      else if &var < &p_75_w then &outVar = 3;
      else &outVar = 4;
    end;
  RUN;

  PROC Means data = &data median noprint;
    class &outVar &sexVar;
    var &var;
    output out = Out_medians median = median;
  RUN;
  DATA _NULL_; set Out_medians;
    if &outVar = 1 and &sexVar = 1 then call symputx("median_1_m", median);
    else if &outVar = 2 and &sexVar = 1 then call symputx("median_2_m", median);
    else if &outVar = 3 and &sexVar = 1 then call symputx("median_3_m", median);
    else if &outVar = 4 and &sexVar = 1 then call symputx("median_4_m", median);
    else if &outVar = 1 and &sexVar ^= 1 then call symputx("median_1_w", median);
    else if &outVar = 2 and &sexVar ^= 1 then call symputx("median_2_w", median);
    else if &outVar = 3 and &sexVar ^= 1 then call symputx("median_3_w", median);
    else if &outVar = 4 and &sexVar ^= 1 then call symputx("median_4_w", median);
  RUN;
  DATA &data; set &data;
    if &outVar = 1 and &sexVar = 1 then &outVar.m = &median_1_m;
    else if &outVar = 2 and &sexVar = 1 then &outVar.m = &median_2_m;
    else if &outVar = 3 and &sexVar = 1 then &outVar.m = &median_3_m;
    else if &outVar = 4 and &sexVar = 1 then &outVar.m = &median_4_m;
    else if &outVar = 1 and &sexVar ^= 1 then &outVar.m = &median_1_w;
    else if &outVar = 2 and &sexVar ^= 1 then &outVar.m = &median_2_w;
    else if &outVar = 3 and &sexVar ^= 1 then &outVar.m = &median_3_w;
    else if &outVar = 4 and &sexVar ^= 1 then &outVar.m = &median_4_w;
  RUN;

  PROC Delete data = percentiles_m; RUN;
  PROC Delete data = percentiles_w; RUN;
  PROC Delete data = Out_medians; RUN;
%mend;



%macro controlQuartiles(data=, var=, caseVar=, outVar=);
  %global p_25 p_50 p_75 median_1 median_2 median_3 median_4;
  PROC Univariate data = &data(where = (&caseVar=0)) noprint;
    var &var;
    output out = percentiles pctlpre=P_ pctlpts=25 to 75 by 25;
  RUN;
  DATA _null_;
    set percentiles;
    call symputx("p_25", p_25);
    call symputx("p_50", p_50);
    call symputx("p_75", p_75);
  RUN;
  %put Quartiles: &p_25 &p_50 &p_75;

  DATA &data; set &data;
    if &var < &p_25 then &outVar = 1;
    else if &var < &p_50 then &outVar = 2;
    else if &var < &p_75 then &outVar = 3;
    else &outVar = 4;
  RUN;

  PROC Means data = &data median noprint;
    class &outVar;
    var &var;
    output out = Out_medians median = median;
  RUN;
  DATA _NULL_; set Out_medians;
    if &outVar = 1 then call symputx("median_1", median);
    else if &outVar = 2 then call symputx("median_2", median);
    else if &outVar = 3 then call symputx("median_3", median);
    else if &outVar = 4 then call symputx("median_4", median);
  RUN;
  DATA &data; set &data;
    if &outVar = 1 then &outVar.m = &median_1;
    else if &outVar = 2 then &outVar.m = &median_2;
    else if &outVar = 3 then &outVar.m = &median_3;
    else if &outVar = 4 then &outVar.m = &median_4;
  RUN;

  PROC Delete data = percentiles; RUN;
  PROC Delete data = Out_medians; RUN;
%mend;



%macro JT_trendTest(Data=, numVars=, trendVar=);
  DATA Out_trend_all;
  RUN;

  ODS Exclude all;
  %let NVars = %sysfunc(countw(&numVars));
  %do i = 1 %to &NVars;
    %let Var_i = %scan(&numVars, &i);

    ODS output JTTest = Out_trend;
    PROC Freq data = &Data;
      table &trendVar * &Var_i / JT noprint;
    RUN;

    DATA Out_trend; set Out_trend;
      where Name1 = "P2_JT";
      p = nValue1;
      var = "&Var_i";
    RUN;

    DATA Out_trend_all;
    set Out_trend_all Out_trend(keep = var p);
    RUN;
  %end;
  ODS Exclude none;
  
  DATA Out_trend_all; set Out_trend_all;
    if _n_ = 1 then delete;
    p = round(p, 0.001);
  RUN;

  PROC Print data = Out_trend_all;
  RUN;

  PROC Datasets nolist;
    delete Out_trend Out_trend_all;
  RUN;
%mend;



%macro CMH_trendTest(Data=, catVars=, trendVar=);
  DATA Out_trend_all;
  RUN;

  ODS Exclude all;
  %let NVars = %sysfunc(countw(&catVars));
  %do i = 1 %to &NVars;
    %let Var_i = %scan(&catVars, &i);

    ODS output CMH = Out_trend;
    PROC Freq data = &Data;
      table &Var_i * &trendVar / CMH noprint;
    RUN;

    DATA Out_trend; set Out_trend;
      length var $10;
      where Statistic = 3;
      p = Prob;
      var = "&Var_i";
    RUN;

    DATA Out_trend_all;
    set Out_trend_all Out_trend(keep = var p);
    RUN;
  %end;
  ODS Exclude none;
  
  DATA Out_trend_all; set Out_trend_all;
    if _n_ = 1 then delete;
    p = round(p, 0.001);
  RUN;

  PROC Print data = Out_trend_all;
  RUN;

  PROC Datasets nolist;
    delete Out_trend Out_trend_all;
  RUN;
%mend;




%macro DummyVars(data=, VarList=);
  DATA &data;
    set &data;
    _Y = 0;
  RUN;
  PROC Glmselect data=&data noprint parmlabelstyle=interlaced(separator='') outdesign(addinputvars)=&data(drop=_Y);
    class &VarList;   
    model _Y = &VarList / noint selection=none;
  RUN;
%mend;


%macro CA_trendTest(Data=, catVars=, trendVar=);
  DATA Out_trend_all;
  RUN;

  ODS Exclude all;
  %let NVars = %sysfunc(countw(&catVars));
  %do i = 1 %to &NVars;
    %let Var_i = %scan(&catVars, &i);

    ODS output TrendTest = Out_trend;
    ODS trace on;
    PROC Freq data = &Data;
      table &Var_i * &trendVar / trend noprint;
    RUN;
    ODS trace off;

    DATA Out_trend; set Out_trend;
      length var $10;
      where Label1 = "Two-sided Pr > |Z|";
      p = nValue1;
      var = "&Var_i";
    RUN;

    DATA Out_trend_all;
    set Out_trend_all Out_trend(keep = var p);
    RUN;
  %end;
  ODS Exclude none;
  
  DATA Out_trend_all; set Out_trend_all;
    if _n_ = 1 then delete;
    p = round(p, 0.001);
  RUN;

  PROC Print data = Out_trend_all;
  RUN;

  PROC Datasets nolist;
    delete Out_trend Out_trend_all;
  RUN;
%mend;


%macro isBlank(param);
  %sysevalf(%superq(param)=,boolean)
%mend;


%macro HR_Phreg(data=, caseVar=, timeVar=, var=, covs=, strata=, classVars=,
                weights=, compact=TRUE, ref=);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;
  ODS Output ParameterEstimates = Out_est;
  PROC Phreg data = &data nosummary;
    %if %isBlank(&classVars) = 0 %then %do;
      %if %isBlank(&ref) = 0 %then %do;
        class &var (ref = "&ref") &classVars;
      %end;
      class &classVars;
    %end;
    model &timeVar * &caseVar (0) = &var &covs / risklimits;
    %if %isBlank(&strata) = 0 %then %do;
      strata &strata;
    %end;
    %if %isBlank(&weights) = 0 %then %do;
      weight &weights;
    %end;
  RUN;
  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;
  DATA Out_est;
    %if %isBlank(&classVars) = 0 %then %do;
      retain Parameter ClassVal0 EST_CI ProbChiSq;
    %end;
    %else %do; 
      retain Parameter EST_CI Prob;
    %end;
    set Out_est;
    where Parameter = "&var";
    EST_CI = cat(strip(put(HazardRatio, comma5.2)), " (",
                 strip(put(HRLowerCL, comma5.2)), " - ",
                 strip(put(HRUpperCL, comma5.2)), ")");
    Prob = put(ProbChiSq, comma6.3);
    %if %isBlank(&classVars) = 0 %then %do;
      keep Parameter ClassVal0 EST_CI Prob;
    %end;
    %else %do; 
      keep Parameter EST_CI Prob;
    %end;
  RUN;
  PROC Print data = Out_est;
  PROC Delete data = Out_est;
  RUN;
  quit;
%mend;


%macro OR_Logistic(data=, caseVar=, var=, covs=, strata=, classVars=,
                compact=TRUE, ref=);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;
  ODS Output ParameterEstimates = Out_est;
  ODS Output OddsRatios = Out_odds;

  PROC Logistic data = &data;
    %if %isBlank(&classVars) = 0 %then %do;
      %if %isBlank(&ref) = 0 %then %do;
        class &var (ref = "&ref") &classVars;
      %end;
      class &classVars;
    %end;
    model &caseVar (event="1") = &var &covs / maxiter = 200;
    %if %isBlank(&strata) = 0 %then %do;
      strata &strata;
    %end;
  RUN;
  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  DATA Out_result;
    merge Out_odds(where = (effect like "&var%")) Out_est(where = (variable="&var"));
    EST_CI = cat(strip(put(OddsRatioEst, comma5.2)), " (",
                 strip(put(LowerCL, comma5.2)), " - ",
                 strip(put(UpperCL, comma5.2)), ")");
    Prob = put(ProbChiSq, comma6.3);
    %if %isBlank(&classVars) = 0 %then %do;
      keep Effect ClassVal0 EST_CI Prob;
    %end;
    %else %do; 
      keep Effect EST_CI Prob;
    %end;
  RUN;
  PROC Print data = Out_result;
  PROC Datasets nodetails nolist;
    delete Out_est Out_odds Out_result;
  RUN;
  quit;
%mend;


%macro OR_Logistic_subset(data=, caseVar=, var=, covs=, strata=, classVars=,
                compact=TRUE, ref=);

  ODS Exclude all;
  %if &compact = FALSE %then %do;
    ODS Exclude none;
  %end;

  DATA newdata;
    set &data;
  RUN;

  %CollectMatches(Data=newdata, OutVar=&caseVar, MatchVar=&strata);
  %CleanIncompleteMatches(Data=newdata, MatchVar=&strata);

  %if &compact = FREQ %then %do;
    ODS Exclude none;
  %end;

  PROC Freq data = newdata;
    table &caseVar / nocol nocum nopercent norow;
  RUN; 

  %if &compact = FREQ %then %do;
    ODS Exclude all;
  %end;

  ODS Output ParameterEstimates = Out_est;
  ODS Output OddsRatios = Out_odds;

  PROC Logistic data = newdata;
    %if %isBlank(&classVars) = 0 %then %do;
      %if %isBlank(&ref) = 0 %then %do;
        class &var (ref = "&ref") &classVars;
      %end;
      class &classVars;
    %end;
    model &caseVar (event="1") = &var &covs / maxiter = 200;
    %if %isBlank(&strata) = 0 %then %do;
      strata &strata;
    %end;
  RUN;

  ODS Exclude none;

  DATA Out_result;
    merge Out_odds(where = (effect like "&var%")) Out_est(where = (variable="&var"));
    EST_CI = cat(strip(put(OddsRatioEst, comma5.2)), " (",
                 strip(put(LowerCL, comma5.2)), " - ",
                 strip(put(UpperCL, comma5.2)), ")");
    Prob = put(ProbChiSq, comma6.3);
    %if %isBlank(&classVars) = 0 %then %do;
      keep Effect ClassVal0 EST_CI Prob;
    %end;
    %else %do; 
      keep Effect EST_CI Prob;
    %end;
  RUN;
  PROC Print data = Out_result;
  PROC Datasets nodetails nolist;
    delete Out_est Out_odds Out_result newdata;
  RUN;
  quit;
%mend;



%macro ICCadjust(value=, icc=1);
  %put %sysfunc(exp(%sysevalf(%sysfunc(log(&value)) / &icc)));
%mend;

/*%ICCadjust(value=1, icc=0.64);*/
/*%ICCadjust(value=2, icc=0.64);*/
/*%ICCadjust(value=4, icc=0.64);*/
/*%ICCadjust(value=0.8, icc=0.64);*/
/*%ICCadjust(value=0.5, icc=0.64);*/


%macro addInter_cox(data=, caseVar=, timeVar=, var1=, var2=, covs=, classVars=, weights=, compact=TRUE);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;
/*  PROC Phreg data = &data nosummary;*/
/*    model &timeVar * &caseVar (0) = &var1 &var2 / risklimits;*/
/*    %if %isBlank(&weights) = 0 %then %do;*/
/*      weight &weights;*/
/*    %end;*/
/*  RUN;*/

  ODS output ParameterEstimates = Out_estimates;
  PROC Phreg data = &data nosummary;
    %if %isBlank(&classVars) = 0 %then %do;
      class &classVars;
    %end;
    model &timeVar * &caseVar (0) = &var1 &var2 &var1*&var2 &covs / risklimits;
    %if %isBlank(&weights) = 0 %then %do;
      weight &weights;
    %end;
  RUN;
  DATA _null_; set Out_estimates;
    if _n_ = 1 then call symputx("Beta1", Estimate);
    if _n_ = 2 then call symputx("Beta2", Estimate);
    if _n_ = 3 then call symputx("Beta3", Estimate);
  RUN;
  %let RERI = %sysevalf(%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - %sysfunc(exp(&Beta1)) - %sysfunc(exp(&Beta2)) + 1);
  %let AP = %sysevalf((%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - %sysfunc(exp(&Beta1)) - %sysfunc(exp(&Beta2)) + 1) / %sysfunc(exp(&Beta1 + &Beta2 + &Beta3)));
  %let S = %sysevalf((%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - 1) / ((%sysfunc(exp(&Beta1)) - 1) + (%sysfunc(exp(&Beta2)) - 1)));
  %put &RERI;
  %put &AP;
  %put &S;
  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  DATA inter_out;
    RERI = &RERI;
    AP = &AP;
    S = &S;
  RUN;
  PROC Print data = inter_out;
  RUN;
%mend;


%macro addInter_CLR(data=, caseVar=, var1=, var2=, covs=, classVars=, strata=, compact=TRUE);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  ODS output ParameterEstimates = Out_estimates;
  PROC Logistic data = &data;
    %if %isBlank(&classVars) = 0 %then %do;
      class &classVars;
    %end;
    model &caseVar (event = "1") = &var1 &var2 &var1*&var2 &covs / maxiter = 200;
    %if %isBlank(&strata) = 0 %then %do;
      strata &strata;
    %end;
  RUN;
  DATA _null_; set Out_estimates;
    if _n_ = 1 then call symputx("Beta1", Estimate);
    if _n_ = 2 then call symputx("Beta2", Estimate);
    if _n_ = 3 then call symputx("Beta3", Estimate);
  RUN;
  %let RERI = %sysevalf(%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - %sysfunc(exp(&Beta1)) - %sysfunc(exp(&Beta2)) + 1);
  %let AP = %sysevalf((%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - %sysfunc(exp(&Beta1)) - %sysfunc(exp(&Beta2)) + 1) / %sysfunc(exp(&Beta1 + &Beta2 + &Beta3)));
  %let S = %sysevalf((%sysfunc(exp(&Beta1 + &Beta2 + &Beta3)) - 1) / ((%sysfunc(exp(&Beta1)) - 1) + (%sysfunc(exp(&Beta2)) - 1)));
  %put &RERI;
  %put &AP;
  %put &S;
  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  DATA inter_out;
    RERI = put(&RERI, comma5.2);
    AP = put(&AP, comma5.2);
    S = put(&S, comma5.2);
  RUN;
  PROC Print data = inter_out;
  RUN;
%mend;



%macro TestNonLinear_RCS5(data=, var=, caseVar=, timeVar=, covs=, classVars=, strata=, weights=, compact=TRUE);

  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  %let knots =5 27_5 50 72_5 95;

  PROC univariate data = &data noprint;
    var &var;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    output out = percentiles pctlpre=P_ pctlpts=1 to 99 by 0.5;
  RUN;

  %do num_i = 1 %to %sysfunc(countw(&knots));
    %let perc_i = %sysfunc(scan(&knots, &num_i));
    DATA _null_; set percentiles;
      call symputx("p_&perc_i", p_&perc_i);
    RUN;
  %end;

  ODS output  TestStmts = test_out;
  PROC Phreg data = &data;
  	%if &classVars ne %then %do;
	    class &classVars;
	  %end;
    model &timeVar * &caseVar (0) = &var spline_1 spline_2 spline_3 &covs /rl ties=efron;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;
    spline_1 = ((&var-&p_5)**3)*(&var>&p_5)
              -((&var-&p_72_5)**3)*(&var>&p_72_5)*(&p_95-&p_5)/(&p_95-&p_72_5)
              +((&var-&p_95)**3)*(&var>&p_95)*(&p_72_5-&p_5)/(&p_95-&p_72_5);
    spline_2 = ((&var-&p_27_5)**3)*(&var>&p_27_5)
              -((&var-&p_72_5)**3)*(&var>&p_72_5)*(&p_95-&p_27_5)/(&p_95-&p_72_5)
              +((&var-&p_95)**3)*(&var>&p_95)*(&p_72_5-&p_27_5)/(&p_95-&p_72_5);
    spline_3 = ((&var-&p_50)**3)*(&var>&p_50)
              -((&var-&p_72_5)**3)*(&var>&p_72_5)*(&p_95-&p_50)/(&p_95-&p_72_5)
              +((&var-&p_95)**3)*(&var>&p_95)*(&p_72_5-&p_50)/(&p_95-&p_72_5);
/*    EFFECT1: TEST &var, spline_1, spline_2, spline_3;*/
    NONLIN1: TEST spline_1, spline_2, spline_3;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
%mend;


%macro TestNonLinear_RCS3(data=, var=, caseVar=, timeVar=, covs=, classVars=, strata=, weights=, compact=TRUE);

  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  %let knots =10 50 90;

  PROC univariate data = &data noprint;
    var &var;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    output out = percentiles pctlpre=P_ pctlpts=1 to 99 by 0.5;
  RUN;

  %do num_i = 1 %to %sysfunc(countw(&knots));
    %let perc_i = %sysfunc(scan(&knots, &num_i));
    DATA _null_; set percentiles;
      call symputx("p_&perc_i", p_&perc_i);
    RUN;
  %end;

  ODS output  TestStmts = test_out;
  PROC Phreg data = &data;
  	%if &classVars ne %then %do;
	    class &classVars;
	  %end;
    model &timeVar * &caseVar (0) = &var spline_1 &covs /rl ties=efron;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;
    spline_1 = ((&var-&p_10)**3)*(&var>&p_10)
              -((&var-&p_50)**3)*(&var>&p_50)*(&p_90-&p_10)/(&p_90-&p_50)
              +((&var-&p_90)**3)*(&var>&p_90)*(&p_50-&p_10)/(&p_90-&p_50);
/*    EFFECT1: TEST &var, spline_1;*/
    NONLIN1: TEST spline_1;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
%mend;



%macro TestNonLinear_CLR_RCS3(data=, var=, caseVar=, covs=, classVars=, strata=, compact=TRUE);

  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  %let knots =10 50 90;

  PROC univariate data = &data noprint;
    var &var;
    output out = percentiles pctlpre=P_ pctlpts=1 to 99 by 0.5;
  RUN;

  %do num_i = 1 %to %sysfunc(countw(&knots));
    %let perc_i = %sysfunc(scan(&knots, &num_i));
    DATA _null_; set percentiles;
      call symputx("p_&perc_i", p_&perc_i);
    RUN;
  %end;

  DATA newdata;
    set &data (keep = &var &caseVar &covs &strata);
    spline_1 = ((&var-&p_10)**3)*(&var>&p_10)
               -((&var-&p_50)**3)*(&var>&p_50)*(&p_90-&p_10)/(&p_90-&p_50)
               +((&var-&p_90)**3)*(&var>&p_90)*(&p_50-&p_10)/(&p_90-&p_50);
  RUN;

  ODS output TestStmts = test_out;
  PROC Logistic data = newdata;
  	%if &classVars ne %then %do;
	    class &classVars;
	  %end;
    model &caseVar (event = "1") = &var spline_1 &covs /rl maxiter = 200;
    %if &strata ne %then %do;
      strata &strata;
    %end;
/*    EFFECT1: TEST &var, spline_1;*/
    NONLIN1: TEST spline_1;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
%mend;



%macro TestNonLinear_squared(data=, var=, caseVar=, timeVar=, covs=, classVars=, strata=, weights=, compact=TRUE);

  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  ODS output  TestStmts = test_out;
  PROC Phreg data = &data;
  	%if &classVars ne %then %do;
	    class &classVars;
	  %end;
    model &timeVar * &caseVar (0) = &var var_sq &covs /rl ties=efron;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;

    var_sq = &var * &var;

    NONLIN1: TEST var_sq;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
%mend;


%macro testInteraction(data=, var1=, var2=, caseVar=, timeVar=, covs=, classVars=, strata=, weights=, compact=TRUE);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  ODS output  TestStmts = test_out;
  PROC Phreg data = &data;
    class &classVars;
    model &timeVar*&caseVar(0) = &var1 &var2 term_interaction &covs;
    weight &weights;
    term_interaction = &var1*&var2;
    interaction: test term_interaction;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
%mend;

%macro testInteraction_CLR(data=, var1=, var2=, caseVar=, covs=, classVars=, strata=, compact=TRUE);
  %if &compact = TRUE %then %do;
    ODS Exclude all;
  %end;

  DATA newdata;
    set &data(keep = &caseVar &var1 &var2 &covs &strata);
    term_interaction = &var1*&var2;
  RUN;

  ODS output  TestStmts = test_out;
  PROC Logistic data = newdata;
    class &classVars;
    model &caseVar(event = "1") = &var1 &var2 term_interaction &covs / maxiter = 200;
    strata &strata;
    interaction: test term_interaction;
  RUN;

  %if &compact = TRUE %then %do;
    ODS Exclude none;
  %end;

  PROC Print data = test_out;
  RUN;
  PROC Datasets nodetails nolist;
    delete newdata test_out;
  RUN;
  quit;
%mend;

/*******************************************************************************
** PROJECT: Lipocalin 2 and CRC                                                *
** PROGRAM NAME: "Figure 1_Splines panel unweighted density.sas"               *
** AUTHOR: Robin Reichmann                                                     *
** SAS-Version: version 9.3_M2                                                 *
** VERSION/DATE: V 1.0 / 2022-20-07                                            *
*                                                                              *
** PURPOSE: The shape of association of lipocalin 2 with CRC and subsites is   *
*           modelled and visualized with restricted cubic splines.             *
*                                                                              *
** PROGRAM SPECIFICATION:                                                      *
*  Lipocalin 2 is non-linearly modelled in conditional logistic    *
*  regression using restricted cubic splines. Associations are modelled for    *
*  CRC and subsites, overall and by sex, also adjusting for matching factors,  *
*  anthropometry, and lifestyle. The modelled splines are visualized ans saved *
*  in a panel figure. Also, p-values for non-linearity are calculated, which   *
*  may be inserted manually into the spline panel, as SAS v9.3_M2 does not yet *
*  allow annotations for PROC SGRender.                                        *
*                                                                              *
** -PREREQUISITES: macro from SAS macros/PlotRCSplines_unlog.sas               * 
*                  macros from SAS macros/Analysis macros.sas                  * 
*                                                                              *
** -READ: lipo_data_weighted.csv                                               *
*         lipo_CC_data_weighted.csv                                            *
*         lipo_CC_prox_data_weighted.csv                                       *
*         lipo_CC_dist_data_weighted.csv                                       *
*         lipo_RC_data_weighted.csv                                            *
*                                                                              *
** -WRITE: Lipocalin2_splinePanel_model2_3k_unweightedDensity.tiff             *
*          Lipocalin2_splinePanel_model2_3k_unweightedDensity.png              *
*******************************************************************************/


* CRC data;
PROC Import datafile = "C:\lipo_data_weighted.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;

* CC data; 
PROC Import datafile = "C:\lipo_CC_data_weighted.csv" out = Lipo_CC_data dbms = csv replace;
  guessingrows = max;
RUN;

* CC_prox data; 
PROC Import datafile = "C:\lipo_CC_prox_data_weighted.csv" out = Lipo_CC_prox_data dbms = csv replace;
  guessingrows = max;
RUN;

* CC_dist data; 
PROC Import datafile = "C:\lipo_CC_dist_data_weighted.csv" out = Lipo_CC_dist_data dbms = csv replace;
  guessingrows = max;
RUN;

* RC data;
PROC Import datafile = "C:\lipo_RC_data_weighted.csv" out = Lipo_RC_data dbms = csv replace;
  guessingrows = max;
RUN;


* Sex specific subsets;
DATA Lipo_data_m; set Lipo_data;
  where sex = 1;
RUN;
DATA Lipo_data_w; set Lipo_data;
  where sex = 2;
RUN;
DATA Lipo_CC_data_m; set Lipo_CC_data;
  where sex = 1;
RUN;
DATA Lipo_CC_data_w; set Lipo_CC_data;
  where sex = 2;
RUN;
DATA Lipo_CC_prox_data_m; set Lipo_CC_prox_data;
  where sex = 1;
RUN;
DATA Lipo_CC_prox_data_w; set Lipo_CC_prox_data;
  where sex = 2;
RUN;
DATA Lipo_CC_dist_data_m; set Lipo_CC_dist_data;
  where sex = 1;
RUN;
DATA Lipo_CC_dist_data_w; set Lipo_CC_dist_data;
  where sex = 2;
RUN;
DATA Lipo_RC_data_m; set Lipo_RC_data;
  where sex = 1;
RUN;
DATA Lipo_RC_data_w; set Lipo_RC_data;
  where sex = 2;
RUN;


* Covariates for adjustment;
%let covs = Center Age_Blood Sex T_Bld_Coll Fasting_C Phrt_Bld Menop_Bld;
%let covs_cat = Center Sex Fasting_C Phrt_Bld Menop_Bld;
%let covs_m = Center Age_Blood T_Bld_Coll Fasting_C;
%let covs_cat_m = Center Fasting_C;
%let covs_w = Center Age_Blood T_Bld_Coll Fasting_C Phrt_Bld Menop_Bld;
%let covs_cat_w = Center Fasting_C Phrt_Bld Menop_Bld;
%let covs2 = Waist_Adj Smoke_Stat Pa_Index Alc_Drinker
             Vegs Fruits Redmeat Procmeat Fish Fibre;
%let covs2_cat = Smoke_Stat Alc_Drinker Pa_Index;

*=============================================================================;


* Model restricted cubic splines and test non-linearity;

%PlotRCSplines_unlog(data = Lipo_data, var = lipocalin2_log, dep_var = CACO,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs &covs2, cov_cat = &covs_cat &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_both2;
  set fig_merge;
  length model $15;
  model = "both2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_data, caseVar=CACO, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_data_m, var = lipocalin2_log, dep_var = CACO,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_m &covs2, cov_cat = &covs_cat_m &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_men2;
  set fig_merge;
  length model $15;
  model = "men2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_data_m, caseVar=CACO, var=lipocalin2_log,
          covs=&covs_m &covs2, classVars=&covs_cat_m &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_data_w, var = lipocalin2_log, dep_var = CACO,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_w &covs2, cov_cat = &covs_cat_w &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_women2;
  set fig_merge;
  length model $15;
  model = "women2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_data_w, caseVar=CACO, var=lipocalin2_log,
          covs=&covs_w &covs2, classVars=&covs_cat_w &covs2_cat,
          strata=Match_caseset);


*=============================================================================;

*CC;
%PlotRCSplines_unlog(data = Lipo_CC_data, var = lipocalin2_log, dep_var = CC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs &covs2, cov_cat = &covs_cat &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_both2;
  set fig_merge;
  length model $15;
  model = "CC_both2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_data, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_data_m, var = lipocalin2_log, dep_var = CC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_m &covs2, cov_cat = &covs_cat_m &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_men2;
  set fig_merge;
  length model $15;
  model = "CC_men2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_data_m, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs_m &covs2, classVars=&covs_cat_m &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_data_w, var = lipocalin2_log, dep_var = CC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_w &covs2, cov_cat = &covs_cat_w &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_women2;
  set fig_merge;
  length model $15;
  model = "CC_women2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_data_w, caseVar=CC_outcome, var=lipocalin2_log,
          covs=&covs_w &covs2, classVars=&covs_cat_w &covs2_cat,
          strata=Match_caseset);


*=============================================================================;

*CC_prox;
%PlotRCSplines_unlog(data = Lipo_CC_prox_data, var = lipocalin2_log, dep_var = CC_prox_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs &covs2, cov_cat = &covs_cat &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_prox_both2;
  set fig_merge;
  length model $15;
  model = "CC_prox_both2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_prox_data_m, var = lipocalin2_log, dep_var = CC_prox_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_m &covs2, cov_cat = &covs_cat_m &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_prox_men2;
  set fig_merge;
  length model $15;
  model = "CC_prox_men2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs_m &covs2, classVars=&covs_cat_m &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_prox_data_w, var = lipocalin2_log, dep_var = CC_prox_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_w &covs2, cov_cat = &covs_cat_w &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_prox_women2;
  set fig_merge;
  length model $15;
  model = "CC_prox_women2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, var=lipocalin2_log,
          covs=&covs_w &covs2, classVars=&covs_cat_w &covs2_cat,
          strata=Match_caseset);


*=============================================================================;

*CC_dist;
%PlotRCSplines_unlog(data = Lipo_CC_dist_data, var = lipocalin2_log, dep_var = CC_dist_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs &covs2, cov_cat = &covs_cat &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_dist_both2;
  set fig_merge;
  length model $15;
  model = "CC_dist_both2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_dist_data_m, var = lipocalin2_log, dep_var = CC_dist_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_m &covs2, cov_cat = &covs_cat_m &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_dist_men2;
  set fig_merge;
  length model $15;
  model = "CC_dist_men2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs_m &covs2, classVars=&covs_cat_m &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_CC_dist_data_w, var = lipocalin2_log, dep_var = CC_dist_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_w &covs2, cov_cat = &covs_cat_w &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_CC_dist_women2;
  set fig_merge;
  length model $15;
  model = "CC_dist_women2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, var=lipocalin2_log,
          covs=&covs_w &covs2, classVars=&covs_cat_w &covs2_cat,
          strata=Match_caseset);


*=============================================================================;

*RC;
%PlotRCSplines_unlog(data = Lipo_RC_data, var = lipocalin2_log, dep_var = RC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs &covs2, cov_cat = &covs_cat &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_RC_both2;
  set fig_merge;
  length model $15;
  model = "RC_both2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_RC_data, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_RC_data_m, var = lipocalin2_log, dep_var = RC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_m &covs2, cov_cat = &covs_cat_m &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_RC_men2;
  set fig_merge;
  length model $15;
  model = "RC_men2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_RC_data_m, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs_m &covs2, classVars=&covs_cat_m &covs2_cat,
          strata=Match_caseset);


%PlotRCSplines_unlog(data = Lipo_RC_data_w, var = lipocalin2_log, dep_var = RC_outcome,
                     strata=Match_caseset, reg_type=logistic,
                     cov = &covs_w &covs2, cov_cat = &covs_cat_w &covs2_cat,
                     knots = 10 50 90,
                     label_variable="LCN2, ng/mL",
                     label_estimate="Adjusted odds ratio (95%-CI)");
DATA Fig_RC_women2;
  set fig_merge;
  length model $15;
  model = "RC_women2";
RUN;
%TestNonLinear_CLR_RCS3(data=Lipo_RC_data_w, caseVar=RC_outcome, var=lipocalin2_log,
          covs=&covs_w &covs2, classVars=&covs_cat_w &covs2_cat,
          strata=Match_caseset);


*=============================================================================;


DATA Fig_merge_all;
  set Fig_both2 Fig_men2 Fig_women2 
      Fig_CC_both2 Fig_CC_men2 Fig_CC_women2 
      Fig_CC_prox_both2 Fig_CC_prox_men2 Fig_CC_prox_women2 
      Fig_CC_dist_both2 Fig_CC_dist_men2 Fig_CC_dist_women2
      Fig_RC_both2 Fig_RC_men2 Fig_RC_women2;
RUN;

%macro SplinePanelTemplate(ModelList=, LabelList=, TemplateName=SplineTemplate, columns=);
  PROC Template;                                     
    define statgraph &TemplateName; 
     begingraph;
      layout lattice / rowgutter=40 columngutter=20 columns=&columns;
        %let NumSplines = %sysfunc(countw(&ModelList, ' '));
          %do i = 1 %to &NumSplines;
            %let model_i = %scan(&ModelList, &i);
            %let label_i = %scan(&LabelList, &i, '|');
            layout lattice / rowweights=(.60 .40) rowgutter=5 columns=1;
              layout overlay
                / xaxisopts = (type=log logopts=(thresholdmin=0 thresholdmax=0) display=none offsetmin=0.01 offsetmax=0.01
                               label=&label_i displaysecondary=(label) labelposition=center labelattrs=(size=10.5 weight=bold))
                  yaxisopts = (type=log logopts=(tickintervalstyle=LOGEXPAND base=2 viewmin=0.35 viewmax=3.1)
                               display=(ticks tickvalues) /*label="OR (95%-CI)" labelattrs=(size=10.5)*/ tickvalueattrs=(size=10));
                bandplot x          = eval(ifn(model eq "&model_i", value_incr,.)) 
                         limitlower = eval(ifn(model eq "&model_i", LowerExp,.))
                         limitupper = eval(ifn(model eq "&model_i", UpperExp,.))
                         / display=(fill outline) outlineattrs=(thickness=0.5 color=grey)
                           fillattrs=(color=lightgrey);
                referenceline y = 1 / lineattrs=(pattern=3 color=black thickness=1);
                seriesplot x = eval(ifn(model eq "&model_i", value_incr,.)) 
                           y = eval(ifn(model eq "&model_i", ExpEstimate,.))
                           / lineattrs=(pattern=1 thickness=2 color=black);
                scatterplot x = eval(ifn(model eq "&model_i", knot_value,.)) 
                            y = eval(ifn(model eq "&model_i", ExpEstimate,.))
                            / markerattrs=(color=black);
              endlayout;
              layout overlay 
                / xaxisopts = (type=log offsetmax=0.01 offsetmin=0.01 label="LCN2, ng/mL" labelattrs=(size=10.5)
                               tickvalueattrs=(size=10) logopts=(thresholdmin=0 thresholdmax=0)) 
                  yaxisopts = (type=linear display=none /*(label) label="Density" labelattrs=(size=10) */
                               offsetmin=0);
                bandplot   x          = eval(ifn(model eq "&model_i", _value_ ,.))
                           limitlower = 0
                           limitupper = eval(ifn(model eq "&model_i", _density_ ,.)) 
                           / datatransparency=0.5 fillattrs=(color=lightgrey);
                seriesplot x = eval(ifn(model eq "&model_i", _value_ ,.)) 
                           y = eval(ifn(model eq "&model_i", _density_ ,.))
                           / lineattrs=(color=grey thickness=1.5);
              endlayout;
              rowheaders;
                layout gridded / columns=1 pad=(right=0);
                  entry "OR (95%-CI)" / textattrs=(size=10.5) rotate=90;
                endlayout;
              endrowheaders;
            endlayout;
          %end;
          rowheaders;
            layout gridded / columns=1 pad=(right=10);
              entry "Colorectal cancer" / textattrs=(size=10.5 weight=bold) rotate=90;
            endlayout;
            layout gridded / columns=1;
              entry "Colon cancer" / textattrs=(size=10.5 weight=bold) rotate=90;
            endlayout;
            layout gridded / columns=1;
              entry "Proximal colon cancer" / textattrs=(size=10.5 weight=bold) rotate=90;
            endlayout;
            layout gridded / columns=1;
              entry "Distal colon cancer" / textattrs=(size=10.5 weight=bold) rotate=90;
            endlayout;
            layout gridded / columns=1;
              entry "Rectal cancer" / textattrs=(size=10.5 weight=bold) rotate=90;
            endlayout;
          endrowheaders;
/*          column2headers;*/
/*            layout gridded / rows=1 pad=(right=10);*/
/*              entry "Both sexes" / textattrs=(size=10.5 weight=bold);*/
/*            endlayout;*/
/*            layout gridded / rows=1;*/
/*              entry "Men" / textattrs=(size=10.5 weight=bold);*/
/*            endlayout;*/
/*            layout gridded / rows=1;*/
/*              entry "Women" / textattrs=(size=10.5 weight=bold);*/
/*            endlayout;*/
/*          endcolumn2headers;*/
        endlayout;
      endgraph;
    end;     
  RUN;
%mend;
%SplinePanelTemplate(ModelList=both2 men2 women2 CC_both2 CC_men2 CC_women2 CC_prox_both2 CC_prox_men2 CC_prox_women2 CC_dist_both2 CC_dist_men2 CC_dist_women2 RC_both2 RC_men2 RC_women2,
                     LabelList="Both sexes"|"Men"|"Women"|"Both sexes"|"Men"|"Women"|"Both sexes"|"Men"|"Women"|"Both sexes"|"Men"|"Women"|"Both sexes"|"Men"|"Women",
                     TemplateName=SplineTemplate, columns=3);


ods _all_ close; 
ods listing gpath = 'H:\Projects\RR4_Lipocalin 2 and CRC\Figs\' image_dpi = 300;
ods graphics / reset = all outputfmt = tiff imagename = "Lipocalin2_splinePanel_model2_3k_unweightedDensity" noborder;
ods graphics / noborder width = 1200 height = 2000 scale = off;
PROC SGRender data = Fig_merge_all template = SplineTemplate; 
RUN;
ods _all_ close; 

ods _all_ close; 
ods listing gpath = 'H:\Projects\RR4_Lipocalin 2 and CRC\Figs\' image_dpi = 300;
ods graphics / reset = all outputfmt = png imagename = "Lipocalin2_splinePanel_model2_3k_unweightedDensity" noborder;
ods graphics / noborder width = 1200 height = 2000 scale = off;
PROC SGRender data = Fig_merge_all template = SplineTemplate; 
RUN;
ods _all_ close; 


*Run SGRender separately in the end to avoid java memory issues;
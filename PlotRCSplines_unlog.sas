%macro PlotRCSplines_unlog(data=, var=, dep_var=,
                     cov=, cov_cat=, strata=, weights=,
                     reg_type=phreg,
                     case_cohort=FALSE, id=,
                     knots=5 27_5 50 72_5 95,
                     spline_res=,
                     ref=50,
                     min_x=, max_x=, min_y=, max_y=,
                     y_axis_type=log,
                     label_variable="Concentration",
                     label_estimate="Hazard ratio (95%-CI)",
                     label_density="Density",
                     color_density=, color_estimate=, color_estimate_CI=,
                     image_name=);

  DATA newdata; set &data;
  RUN;

  * Compute all integer percentiles of the variable of interest;
  PROC univariate data = newdata noprint;
    var &var;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    output out = percentiles pctlpre=P_ pctlpts=1 to 99 by 0.5;
  RUN;

  * ALWAYS extract 3%, 50% and 97% percentiles
    (3% and 97% for x-axis limits, 50% as the spline reference);
  DATA _null_; set percentiles;
    call symputx("p_ref", p_&ref);
    call symputx("p_3", p_3);
/*    call symputx("p_5", p_5);*/
    call symputx("p_50", p_50);
/*    call symputx("p_95", p_95);*/
    call symputx("p_97", p_97);
  RUN;
  %put reference = &p_ref;
  %put p1 = &p_3;
  %put p50 = &p_50;
  %put p99 = &p_97;

  %if &spline_res eq %then %do;
    %let spline_res = %sysevalf(%sysfunc(abs(&p_97)) / 1000);
  %end;

  * Extract knot values of desired percentiles (for the knotmethod=list()
    parameter in the PROC Phreg/Logistic effect statement);
  %global knots_values;
  %let knots_values =;
  %do num_i = 1 %to %sysfunc(countw(&knots));
    %let perc_i = %sysfunc(scan(&knots, &num_i));
    DATA _null_; set percentiles;
      call symputx("p_&perc_i", p_&perc_i);
    RUN;
    %let knots_values = &knots_values &&p_&perc_i;
  %end;
  %let knots_values = %sysfunc(tranwrd(&knots_values, %str( ), %str(,)));

  * Extract density distribution of the specified variable;
  PROC univariate data = newdata(where = (&var >= &p_3 and &var <= &p_97)) noprint;
    var &var;
/*    %if &weights ne %then %do;*/
/*      weight &weights;*/
/*    %end;*/
    histogram &var / kernel outkernel=KerDensity noplot;
  RUN;
  DATA KerDensity; set KerDensity;
    format _density_ 5.3;
    _density_ = _density_ * 10;
  RUN;

  %if &reg_type eq phreg %then %do;
  * Generate a Cox regression model with specified parameters and using splines of the specified variable;
  PROC Phreg data = newdata %if &case_cohort eq TRUE %then %do; covsandwich(aggregate) %end; noprint;
  	%if &cov_cat ne %then %do;
	  class &cov_cat;
	%end;
    effect spl_&var = spline(&var /details naturalcubic basis=tpf(noint) knotmethod=list(&knots_values));
    model &dep_var = spl_&var &cov /rl ties=efron;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;
    %if &id ne %then %do;
      id &id;
    %end;
    store reg_model;
  RUN;
  %end;
  %else %if &reg_type eq logistic %then %do;
  * Generate a logistic regression model with specified parameters and using splines of the specified variable;
  PROC Logistic data = newdata noprint;
  	%if &cov_cat ne %then %do;
	  class &cov_cat;
	%end;
    effect spl_&var = spline(&var /details naturalcubic basis=tpf(noint) knotmethod=list(&knots_values));
    model &dep_var.(event = "1") = spl_&var &cov / rl maxiter = 200;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;
    store reg_model;
  RUN;
  %end;
  %else %if &reg_type eq linear %then %do;
  PROC Glmselect data = newdata noprint;
   	%if &cov_cat ne %then %do;
		class &cov_cat;
	%end;
    effect spl_&var = spline(&var /details naturalcubic basis=tpf(noint) knotmethod=list(&knots_values));
    model &dep_var = spl_&var &cov / selection = none;
    %if &weights ne %then %do;
      weight &weights;
    %end;
    %if &strata ne %then %do;
      strata &strata;
    %end;
    store reg_model;
    output out=SplineOut predicted=Fit;
  RUN;
  %end;

  * Generate spline coordinates for a partial plot (variable concentration VS hazard ratio);
  ods exclude all;
  ODS dataset estimates = estimates;
  PROC Plm restore = reg_model; * from p_3 to p_97 with increments specified in spline_res;
    %let value = &p_3;
    %do %while (%sysevalf(&value <= &p_97));
      estimate "&value." spl_&var [-1, &p_ref] [1, &value] /exp cl; * [-1, &ref] sets reference to [1, &value];
      %let value = %sysevalf(&value + &spline_res);
    %end;
    %do k_i = 1 %to %sysfunc(countw(%quote(&knots_values), %str(,))); * include knot values;
      %let value = %scan(%quote(&knots_values), &k_i, %str(,));
      estimate "&value." spl_&var [-1, &p_ref] [1, &value] /exp cl;
    %end;
  RUN;
  ods exclude none;

  DATA estimates; set estimates;
    value_incr = label * 1; * rename and handle as numeric;
    drop label;
  RUN;
  %do k_i = 1 %to %sysfunc(countw(%quote(&knots_values), %str(,))); * have knot values in separate column;
    %let value = %scan(%quote(&knots_values), &k_i, %str(,));
    DATA estimates; set estimates;
      if value_incr = &value then knot_value = &value;
    RUN;
  %end;
  PROC Sort Data = estimates;
    by value_incr;
  RUN;

  * Define min and max hazard ratios for plotting the splines;
  %if &max_y eq %then %do;
    PROC SQL noprint;
      select UpperExp + (0.05 * UpperExp)
      into :max_y
      from estimates
      having UpperExp = max(UpperExp);
    quit;
  %end;
  %if &min_y eq %then %do;
    PROC SQL noprint;
      select LowerExp - (0.05 * LowerExp)
      into :min_y
      from estimates
      having LowerExp = min(LowerExp);
    quit;
  %end;
  %if %sysevalf(&min_y < 0.01) %then %do;
    %let min_y = 0.01;
  %end;

  * Define min and max concentration for plotting the splines;
  %if &max_x eq %then %do;
/*    %let max_x = %sysfunc(exp(&p_97));*/
      %let max_x = %sysevalf(2**&p_97);
  %end;
  %if &min_x eq %then %do;
/*    %let min_x = %sysfunc(exp(&p_3));*/
      %let min_x = %sysevalf(2**&p_3);
  %end;

  * Define density plot color;
  %if &color_density eq %then %do;
    %let color_density = CX00925B;
  %end;
  %if &color_estimate eq %then %do;
    %let color_estimate = crimson;
  %end;
  %if &color_estimate_CI eq %then %do;
    %let color_estimate_CI = crimson;
  %end;

  * Prepare a dataset for plotting;
  DATA fig_merge;
    set estimates KerDensity(keep = _density_ _VALUE_);
/*    _value_ = exp(_value_);*/
/*    value_incr = exp(value_incr);*/
/*    knot_value = exp(knot_value);*/
    _value_ = 2**_value_;
    value_incr = 2**value_incr;
    knot_value = 2**knot_value;
  RUN;

  * Define partial plot and density plot template;
  PROC template;                                     
    define statgraph temp; 
      begingraph;
        layout lattice / rowweights = (.70 .30) rowgutter = 5 columns = 1;

          layout overlay  / xaxisopts = (type=log logopts=(viewmin=&min_x viewmax=&max_x thresholdmin = 0 thresholdmax = 0)
                                         display=none offsetmin=0 offsetmax=0) 
                            yaxisopts = (%if &y_axis_type eq log %then %do;
                                           type=log logopts= (viewmin=&min_y viewmax=&max_y)
                                         %end; %else %if &y_axis_type eq linear %then %do;
                                           type=linear linearopts= (viewmin=&min_y viewmax=&max_y)
                                         %end; 
                                         display=(tickvalues ticks label) label=&label_estimate labelattrs=(size=10) tickvalueattrs=(size=8)); 
 
            bandplot    x          = value_incr 
                        limitlower = LowerExp
                        limitupper = UpperExp / display = (outline) outlineattrs=(color=&color_estimate_CI pattern=2);

            seriesplot  x = value_incr
                        y = ExpEstimate / lineattrs=(color=&color_estimate pattern=1 thickness=2);

            referenceline y = 1 / lineattrs = (pattern=34);

          endlayout;
          layout overlay / xaxisopts = (type=log logopts=(viewmin=&min_x viewmax=&max_x thresholdmin = 0 thresholdmax = 0)
                                        display=(ticks tickvalues label) label=&label_variable labelattrs=(size=10) 
                                        tickvalueattrs=(size=8) offsetmin=0 offsetmax=0)
                           yaxisopts = (type=linear display=(label) label=&label_density
                                        labelattrs=(size=10) offsetmin=0);
/*                           yaxisopts = (type=linear display=(ticks tickvalues label) label=&label_density*/
/*                                        labelattrs=(size=8) tickvalueattrs=(size=6) offsetmin=0);*/

            seriesplot x = _value_ y = _density_ / lineattrs=(color=&color_density thickness=1.5);

            bandplot   x = _value_ limitlower = 0 limitupper = _density_ / datatransparency=0.5 fillattrs=(color=&color_density);

          endlayout;  
        endlayout;
      endgraph;
    end;     
  RUN;

  * Plot splines;
/*  ods _all_ close; */
  %if &image_name ne %then %do;
    ods listing gpath='D:\R_work\Robin\' image_dpi = 300;
    ods graphics / reset=all outputfmt=tiff imagename=&image_name;
  %end;
  ods graphics / noborder width=600px height=600px scale=off;
  PROC Sgrender data = fig_merge template = temp; 
  RUN; quit;
%mend;

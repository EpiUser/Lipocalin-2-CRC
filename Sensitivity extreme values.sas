

* CRC data;
PROC Import datafile = "C:\lipo_data_weighted.csv" out = Lipo_data dbms = csv replace;
  guessingrows = max;
RUN;

PROC univariate data = Lipo_data noprint;
  class Sex;
  var lipocalin2_ngml;
  output out = percentiles pctlpre=P_ pctlpts=1 99;
RUN;
DATA _null_; set percentiles(where = (Sex = 1));
  call symputx("p_1_m", p_1);
  call symputx("p_99_m", p_99);
RUN;
DATA _null_; set percentiles(where = (Sex = 2));
  call symputx("p_1_w", p_1);
  call symputx("p_99_w", p_99);
RUN;

DATA Lipo_data; set Lipo_data;
  if Sex = 1 and lipocalin2_ngml <= &p_1_m then delete;
  if Sex = 1 and lipocalin2_ngml >= &p_99_m then delete;
  if Sex = 2 and lipocalin2_ngml <= &p_1_w then delete;
  if Sex = 2 and lipocalin2_ngml >= &p_99_w then delete;
RUN;
%put men: p1 = &p_1_m, p99 = &p_99_m;
%put women: p1 = &p_1_w, p99 = &p_99_w;

%CollectMatches(Data=Lipo_data, OutVar=CACO, MatchVar=Match_Caseset); * remove unmatched controls;
%CleanIncompleteMatches(Data=Lipo_data, MatchVar=Match_Caseset); * remove unmatched cases;
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
%HR_Phreg(data=Lipo_data, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data(where=(Waist_cut = 0)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data(where=(Waist_cut = 1)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Men";
%HR_Phreg(data=Lipo_data_m, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m(where=(Waist_cut = 0)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_m(where=(Waist_cut = 1)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Women";
%HR_Phreg(data=Lipo_data_w, caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w(where=(Waist_cut = 0)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_data_w(where=(Waist_cut = 1)), caseVar=CACO, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);



* CC data; 
PROC Import datafile = "C:\lipo_CC_data_weighted.csv" out = Lipo_CC_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_data; set Lipo_CC_data;
  if Sex = 1 and lipocalin2_ngml <= &p_1_m then delete;
  if Sex = 1 and lipocalin2_ngml >= &p_99_m then delete;
  if Sex = 2 and lipocalin2_ngml <= &p_1_w then delete;
  if Sex = 2 and lipocalin2_ngml >= &p_99_w then delete;
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
%HR_Phreg(data=Lipo_CC_data, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data(where=(Waist_cut = 0)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data(where=(Waist_cut = 1)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Men";
%HR_Phreg(data=Lipo_CC_data_m, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m(where=(Waist_cut = 0)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_m(where=(Waist_cut = 1)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Women";
%HR_Phreg(data=Lipo_CC_data_w, caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w(where=(Waist_cut = 0)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_data_w(where=(Waist_cut = 1)), caseVar=CC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);


* CC proximal data; 
PROC Import datafile = "C:\lipo_CC_prox_data_weighted.csv" out = Lipo_CC_prox_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_prox_data; set Lipo_CC_prox_data;
  if Sex = 1 and lipocalin2_ngml <= &p_1_m then delete;
  if Sex = 1 and lipocalin2_ngml >= &p_99_m then delete;
  if Sex = 2 and lipocalin2_ngml <= &p_1_w then delete;
  if Sex = 2 and lipocalin2_ngml >= &p_99_w then delete;
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
%HR_Phreg(data=Lipo_CC_prox_data, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Men";
%HR_Phreg(data=Lipo_CC_prox_data_m, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_m(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Women";
%HR_Phreg(data=Lipo_CC_prox_data_w, caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w(where=(Waist_cut = 0)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_prox_data_w(where=(Waist_cut = 1)), caseVar=CC_prox_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);


* CC distal data; 
PROC Import datafile = "C:\lipo_CC_dist_data_weighted.csv" out = Lipo_CC_dist_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_CC_dist_data; set Lipo_CC_dist_data;
  if Sex = 1 and lipocalin2_ngml <= &p_1_m then delete;
  if Sex = 1 and lipocalin2_ngml >= &p_99_m then delete;
  if Sex = 2 and lipocalin2_ngml <= &p_1_w then delete;
  if Sex = 2 and lipocalin2_ngml >= &p_99_w then delete;
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
%HR_Phreg(data=Lipo_CC_dist_data, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Men";
%HR_Phreg(data=Lipo_CC_dist_data_m, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_m(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Women";
%HR_Phreg(data=Lipo_CC_dist_data_w, caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w(where=(Waist_cut = 0)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_CC_dist_data_w(where=(Waist_cut = 1)), caseVar=CC_dist_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);



* RC data;
PROC Import datafile = "C:\lipo_RC_data_weighted.csv" out = Lipo_RC_data dbms = csv replace;
  guessingrows = max;
RUN;
DATA Lipo_RC_data; set Lipo_RC_data;
  if Sex = 1 and lipocalin2_ngml <= &p_1_m then delete;
  if Sex = 1 and lipocalin2_ngml >= &p_99_m then delete;
  if Sex = 2 and lipocalin2_ngml <= &p_1_w then delete;
  if Sex = 2 and lipocalin2_ngml >= &p_99_w then delete;
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
%HR_Phreg(data=Lipo_RC_data, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data(where=(Waist_cut = 0)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data(where=(Waist_cut = 1)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Men";
%HR_Phreg(data=Lipo_RC_data_m, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m(where=(Waist_cut = 0)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_m(where=(Waist_cut = 1)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
TITLE "Women";
%HR_Phreg(data=Lipo_RC_data_w, caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w(where=(Waist_cut = 0)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);
%HR_Phreg(data=Lipo_RC_data_w(where=(Waist_cut = 1)), caseVar=RC_outcome, timeVar=time_w, var=lipocalin2_log,
          covs=&covs &covs2_w, classVars=&covs_cat &covs2_cat,
          weights=sampling_weights);



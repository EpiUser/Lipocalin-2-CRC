

*----------------------------------------------------------------
MACRO: Generate a summary table with the percentage of categories
       of each categorical characteristic
       -> "Data"_CatSum table with results
       -> save as a local file if desired
----------------------------------------------------------------;
%macro CatSummary(Data=, VarNames=);
  * Calculate percentage of each variables categories;
  ods exclude all;
  ods output OneWayFreqs = CatSum;
  PROC Freq data = &Data;
    table &VarNames;
  RUN;
  ods exclude none;

  * Rearrange results of percentage computation;
  DATA CatSum; set CatSum(rename = (Table = Variable));
    Variable = tranwrd(Variable, "Table ", "");
    format Percent 10.1;
    length Category $50;
  RUN;

  * Gather the values of each F_"Variable" column to one category column
    with category names;
  %let NumVars = %sysfunc(countw(&VarNames));
  %do i = 1 %to &NumVars;
    %let VarName = %scan(&VarNames, &i);
    DATA CatSum; set CatSum;
      Category = catt(Category, F_&VarName);
    RUN;
  %end;

  * Drop unnessecary variables and rearrange the order of columns;
  DATA CatSum;
    retain Variable Category Percent;
    set CatSum;
    keep Variable Category Percent;
  RUN;
%mend;

*---------------------------------------------------------------------
MACRO: Generate a summary table with mean and std of specified numeric 
       characteristics
       -> "Data"_NumSum table with results
       -> save as a local file if desired
---------------------------------------------------------------------;
%macro NumSummary(Data=, VarNames=);
  PROC Univariate data = &data normal outtable = NumSum noprint;
    variables &VarNames;
  RUN;

  * Combine Mean and Std into one string;
  DATA NumSum; set NumSum;
    length Mean_Std $50 Median_IQR $50 Best $50;
    Variable = _Var_;
    Mean_Std = cat(strip(put(_Mean_, comma20.1)), " (", strip(put(_Std_, comma20.1)), ")");
    Median_IQR = cat(strip(put(_Median_, comma20.1)), " (", strip(put(_Q1_, comma20.1)), " - ", strip(put(_Q3_, comma20.1)), ")");
    Normal = strip(put(_NORMAL_, comma20.3));

    if Normal < 0.05 then Best = Mean_Std;
    else Best = Median_IQR;

    keep Variable Mean_Std Median_IQR Normal Best;
  RUN;
%mend;

*---------------------------------------------------------------------
MACRO: Generate a population summary table with summary statistics for 
       categorical and numerical characteristics 
       -> "Data"_FullSum table with results
       -> save as a local file if desired
---------------------------------------------------------------------;
%macro PopulationSummary(Data=, VarsCat=, VarsNum=);
  title "Variable summary: &Data";

  * Generate summary tables for categorical and numerical variables;
  %if &VarsCat ne %then %do; 
    %CatSummary(Data = &Data,
                VarNames = &VarsCat);
  %end;

  %if &VarsNum ne %then %do; 
    %NumSummary(Data = &Data,
                VarNames = &VarsNum);
  %end;

  %if &VarsCat ne and &VarsNum ne %then %do; 
    * Combine categorical and numerical summary tables;
    DATA FullSum; set CatSum;
      length Mean_Std $50 Median_IQR $50 Best $50 Normal $50;
    RUN;
    PROC Append base = FullSum 
                data = NumSum force nowarn;
    RUN;

    * Print summary table for inspection or copy/paste export;
    Proc Print data = FullSum;
    RUN;

    * Remove separate datasets;
    PROC Datasets nolist;
      delete CatSum NumSum;
    RUN;
  %end;

  %else %do;
    %if &VarsCat ne %then %do; 
      Proc Print data = CatSum;
      RUN;
      PROC Datasets nolist;
        delete CatSum;
      RUN;
    %end;
    %if &VarsNum ne %then %do; 
      Proc Print data = NumSum;
      RUN;
      PROC Datasets nolist;
        delete NumSum;
      RUN;
    %end;
  %end;

  title;
  quit;
%mend;





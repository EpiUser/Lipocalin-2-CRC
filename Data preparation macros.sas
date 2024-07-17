
*-----------------------------------------------------------------------
MACRO: Count entries with missing values (VERY FAST for continuous data)
-----------------------------------------------------------------------;
%macro CountMissings(Data=, VarNames=, VarFormat="num");
  * Generate binary format for missing values (for char or num values);
  PROC Format;
    value  $Missing_Format ' '='Missing' other = 'Not Missing';
    value  Missing_Format . = 'Missing' other = 'Not Missing';
  RUN;
  ODS Exclude all;
  ODS output OneWayFreqs = Out_freq;
  * Count missing value frequency (depending on "VarFormat");
  PROC Freq data = &Data;
    %if &VarFormat = "char" %then %do;
      format &VarNames $Missing_Format.;
    %end;
    %if &VarFormat = "num" %then %do;
      format &VarNames Missing_Format.;
    %end;
    tables &VarNames / missing missprint nocum;
  RUN;
  ODS Exclude none;
  DATA Out_freq;                                                                     
    retain table column_value;                                                     
    set Out_freq;                                                                   
    keep table column_value frequency percent;             
    table=scan(table,2,' ');                                                       
    column_value=trim(left(vvaluex(table)));                                                                                
  RUN;
  PROC Print data = out_freq(where = (column_value = "Missing"));
  RUN;
%mend;


*-------------------------------------------------
MACRO: Remove entries with missing "VarName" value
-------------------------------------------------;
%macro RemoveMissing(Data=, VarName=);
  DATA &Data; set &Data;
    where &VarName ne .;
  RUN;
%mend;


*------------------------------------------------------
MACRO: Apply RemoveMissing macro to a list of variables
------------------------------------------------------;
%macro RemoveMissingList(Data=, VarList=);
  * Count the number of variables in "VarList";
  %let NumWords = %sysfunc(countw(&VarList));
  * Use RemoveExtreme macro for every variable;
  %do i = 1 %to &NumWords;
    %let VarName = %scan(&VarList, &i);
    %RemoveMissing(Data = &Data, VarName = &VarName);
  %end;
%mend;


*------------------------------------------------
MACRO: Remove casesets without case
------------------------------------------------;
%macro CollectMatches(Data=, OutVar=, MatchVar=);
  PROC Freq data = &Data noprint;
    table &MatchVar * &OutVar / out = Match_Freqs;
  RUN;
  PROC SQL noprint;
    select &MatchVar
    into :MatchesWithCase separated by ' '
    from Match_Freqs
    where &OutVar = 1;
  RUN;

  DATA data_matched;
  RUN;
  %do i = 1 %to %sysfunc(countw(&MatchesWithCase));
    DATA data_matched;
      set data_matched &Data(where = (Match_Caseset = %scan(&MatchesWithCase, &i)));
    RUN;
  %end;
  DATA &Data; set data_matched;
    if _n_ = 1 then delete;
  RUN;
  PROC Delete data = data_matched; RUN;
%mend;


*------------------------------------------------
MACRO: Remove all incomplete Case/Control matches
------------------------------------------------;
%macro CleanIncompleteMatches(Data=, MatchVar=, IncompleteCount=1);
  * Find unique matching variables;
  PROC Freq data = &Data nlevels noprint;
    tables &MatchVar / out = Unique_MatchVars (where = (count = &IncompleteCount));
  RUN;
  * Gather unique matching variable into a macro variable list;
  %let Incomplete_Matches = 0; * will be overwritten;
  PROC Sql noprint;
    select &MatchVar
      into :Incomplete_Matches separated by ' '
      from Unique_MatchVars
  quit;
  * Remove incomplete matches one by one;
  %do i = 1 %to %sysfunc(countw(&Incomplete_Matches));
    DATA &Data; set &Data;
      where &MatchVar ne %scan(&Incomplete_Matches, &i);
    RUN;
  %end;

  PROC Delete data = Unique_MatchVars; RUN;
%mend;


*------------------------------------------------
MACRO: Remove additional controls
(Run more than once, if more than 3 controls)
------------------------------------------------;
%macro RemoveAdditionalControls(Data=, MatchVar=, CntlNumVar=);
  PROC Freq data = &Data noprint;
    table &MatchVar / out = Casesets_3 (where = (count = 3));
  RUN;
  PROC Sql noprint;
    select &MatchVar
      into :Casesets_3_list separated by ' '
      from Casesets_3
  quit;
  RUN;
  %put Casesets with additional controls: &Casesets_3_list;

  %do i = 1 %to %sysfunc(countw(&Casesets_3_list));
    DATA &Data;
      set &Data;
      if &MatchVar eq %scan(&Casesets_3_list, &i) then do;
        if &CntlNumVar = 2 then delete;
      end;
    RUN;
  %end;
%mend;


*------------------------------------
MACRO: Z tranform specified variables
------------------------------------;
%macro ZTransform(Data=, VarList=);
  %let NumVars = %sysfunc(countw(&VarList));
  %do i = 1 %to &NumVars;
    %let VarName = %scan(&VarList, &i);

    DATA &Data; set &Data;
      &VarName._Z = &VarName;
    RUN;

    PROC Standard data = &Data mean = 0 std = 1 out = &Data;
      var &VarName._Z;
    RUN;
  %end;
%mend;


/***********************************************************************************************
Name: EB_cognition_all_20180320
Purpose: Create empirical bayes (EB) estimates for cognitive domains in SATSA
Study: EWAS of cognition

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB) 
Karolinska Institutet, Stockholm, Sweden

   	Created: 	20180320 by Ida Karlsson (IK), based on code by Malin Ericsson
	Updated: 	20191014 by IK - Add IPT10 data. Split memory into Thurstone and digit span. 
				20191216 by IK - Convert memory (Thurstone and digit span) to T-scores.

  	SAS v9.4

************************************************************************************************/
/**************************************** 1. LIBNAMES ******************************************/
/***********************************************************************************************/

libname resp 'Z:\Data\SATSA_data';
libname raw_cog 'P:\Dementia_IK\Cognition\Cognitive_data\Original_data';
libname eb_cog 'P:\Dementia_IK\Cognition\EB_estimates\Data';

/***********************************************************************************************/
/********************************** 2. XXXXXXXXXXXXXXXXXXX *************************************/
/***********************************************************************************************/

data resp0917; set resp.resp0917;
keep twinnr sex pairid;
run;

/* NOTE! Individuals who were under 50 were omitted when calculating the factor scores. Will do
   the same for thurstone and digit span */
data SATSA_cog_pre; 
	merge raw_cog.cogfac_1to10_10jul19 (in=a)
		  resp0917 (in=b);
	by twinnr;
	if a and b;

	if not50_i1=1 then do;
		thur_z1=.;
		dgsp_z1=.;
	end;
	if not50_i2=1 then do;
		thur_z1=.;
		dgsp_z1=.;
	end;
	if not50_i3=1 then do;
		thur_z1=.;
		dgsp_z1=.;
	end;

	esex=sex-1.5; 
	dsex=sex-1;

		ARRAY tmt (m) thur_z1-thur_z10;
		ARRAY tmd (m) dgsp_z1-dgsp_z10;
		ARRAY tv (m) tverb1-tverb10;
		ARRAY tsd (m) tsped1-tsped10;
		ARRAY tsp (m) tspat1-tspat10;
		Array tpc (m) tpcom1-tpcom10;
		ARRAY dm (m) dem1 dem2 dem3 dem4 dem5 dem6 dem7 dem8 dem9 dem10;
		ARRAY ag (m) iage1 iage2 iage3 iage4 iage5 iage6 iage7 iage8 iage9 iage10;
		do m=1 to 10;
		 tmemot_pre= tmt;
		 tmemod_pre= tmd;
		 tverb= tv;
		 tsped= tsd;
		 tspat= tsp;
		 tpcom=tpc;

		 dem=dm;
		 age = ag;


	
		 /*age variables*/
		  
                 cage75 = ag-75;      
                 cagesq75 = cage75**2;  

				 cage70 = ag-70;      
                 cagesq70 = cage70**2;  
 
                 cage65 = ag-65;      
                 cagesq65 = cage65**2; 


		 /*spline*/
		 cageA75 = cage75; if cage75 gt 0 then cageA75=0; 
 		 cageB75 = cage75; if cage75 lt 0 and cage75 ne . then cageB75=0;

		 cageA70 = cage70; if cage70 gt 0 then cageA70=0; 
 		 cageB70 = cage70; if cage70 lt 0 and cage70 ne . then cageB70=0;

		 cageA65 = cage65; if cage65 gt 0 then cageA65=0; 
 		 cageB65 = cage65; if cage65 lt 0 and cage65 ne . then cageB65=0;


		OUTPUT;
		END;

keep twinnr pairid sex esex dsex dem przygup 
tmemot_pre tmemod_pre tverb tsped tspat tpcom age 
cage75 cagesq75 cage70 cagesq70 cage65 cagesq65 cageA70 cageB70  cageA65 cageB65;
run;

/* Drop those missing across all domains */
/* NOTE! Not including the memory domains since those are less "cleaned" wrt to dementa */
/* Convert Thurstone and digit span to T-scores */
data SATSA_cog;
	set SATSA_cog_pre;
		tmemot=(tmemot_pre*10)+50; 
		tmemod=(tmemod_pre*10)+50; 
		if tverb =. and tsped =. and tspat =. and tpcom =. then delete;
	drop tmemot_pre tmemod_pre;
run;
*n=3836;

/***********************************************************************************************/
/*Verbal*/

/*Verbal - Linear*/
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tverb ne .;
model tverb= cage70 /ddfm=bw solution;
random intercept cage70 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tverb_EB_l;
title 'M1: Tverb linear model Int. at 70';
run;
/* AIC: 21380.6 */

/*Verbal - Quadratic*/
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tverb ne .;
model tverb= cage70 cagesq70 /ddfm=bw solution;
random intercept cage70 cagesq70 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tverb_EB_q;
title 'M1: Tverb quadratic model Int. at 70';
run;
/* AIC: 21147.7 */
/* > Quadratic model has better fit */

/* Intercept*/
proc sort data=tverb_EB_q; by twinnr; run;
data int_tverb;
  set tverb_EB_q (where=(effect="Intercept"))
       ;
  by twinnr;
  qI70_tverb=estimate;
  qI70_tverbse=StdErrPred;
  qI70_tverbse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI70_tverb; 
run;

/* Linear slope */
data A;
  merge tverb_EB_q (where=(effect="cage70"))
        int_tverb;
  by twinnr;
  qL70_tverb=estimate;
  qL70_tverbse=StdErrPred;
  qL70_tverbse_inv=1/StdErrPred;

  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL70_tverb; 
run;

/* Quadratic slope */
data tverb_EB_quad;
  merge tverb_EB_q (where=(effect="cagesq70"))
        A;
  by twinnr;
  qQ70_tverb=estimate;
  qQ70_tverbse=StdErrPred;
  qQ70_tverbse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means data=tverb_EB_quad;
var  qQ70_tverb; 
run;

/***********************************************************************************************/
/*Processing speed*/

/* Processing speed - Linear*/
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tsped ne .;
model tsped= cage65 /ddfm=bw solution;
random intercept cage65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tsped_EB_l;
title 'M1: Tsped linear model Int. at 65';
run;
/* AIC: 23622.3 */

/* Processing speed - Quadratic*/
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tsped ne .;
model tsped= cage65 cagesq65 /ddfm=bw solution;
random intercept cage65 cagesq65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tsped_EB_q;
title 'M1: Tsped quadratic model Int. at 65';
run;
/* AIC: 23481.6 */
/* > Quadratic model has better fit */

/* Intercept*/
proc sort data=tsped_EB_q; by twinnr; run;
data int_tsped;
  set tsped_EB_q (where=(effect="Intercept"))
       ;
  by twinnr;
  qI65_tsped=estimate;
  qI65_tspedse=StdErrPred;
  qI65_tspedse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI65_tsped; 
run;

/* Linear slope */
data B;
 merge tsped_EB_q (where=(effect="cage65"))
        int_tsped;
  by twinnr;
  qL65_tsped=estimate;
  qL65_tspedse=StdErrPred;
  qL65_tspedse_inv=1/StdErrPred;

drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL65_tsped; 
run;

/* Quadratic slope */
data tsped_EB_quad;
  merge tsped_EB_q (where=(effect="cagesq65"))
        B;
  by twinnr;
  qQ65_tsped=estimate;
  qQ65_tspedse=StdErrPred;
  qQ65_tspedse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means data=tsped_EB_quad;
var  qQ65_tsped; 
run;

/***********************************************************************************************/
/*Memory, Thurstone */

/* Memory Linear */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tmemot ne .;
model tmemot= cage65 /ddfm=bw solution;
random intercept cage65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tmemot_EB_l;
title 'M1: Tmemo linear model, Thurstone Int. at 65';
run;
/* AIC: 24599.8 */

/* Memory Quadratic */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tmemot ne .;
model tmemot= cage65 cagesq65 /ddfm=bw solution;
random intercept cage65 cagesq65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tmemot_EB_q;
title 'M1: Tmemo quadratic model, Thurstone Int. at 65';
run;
/* AIC: 24492.8 */
/* > Quadratic model has better fit */

/* Intercept*/
proc sort data=tmemot_EB_q; by twinnr; run;
data int_tmemot;
  set tmemot_EB_q (where=(effect="Intercept"))
       ;
  by twinnr;
  qI65_tmemot=estimate;
  qI65_tmemotse=StdErrPred;
  qI65_tmemotse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI65_tmemot; 
run;

/* Linear slope */
data C;
 merge tmemot_EB_q (where=(effect="cage65"))
        int_tmemot;
  by twinnr;
  qL65_tmemot =estimate;
  qL65_tmemotse=StdErrPred;
  qL65_tmemotse_inv=1/StdErrPred;

drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL65_tmemot; 
run;


/* Quadratic slope */
data tmemot_EB_quad;
  merge tmemot_EB_q (where=(effect="cagesq65"))
        C;
  by twinnr;
  qQ65_tmemot=estimate;
  qQ65_tmemotse=StdErrPred;
  qQ65_tmemotse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means data=tmemot_EB_quad;
var  qQ65_tmemot; 
run;

/***********************************************************************************************/
/*Memory, digit span */

/* Memory Linear */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tmemod ne .;
model tmemod= cage65 /ddfm=bw solution;
random intercept cage65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tmemod_EB_l;
title 'M1: Tmemo linear model, digit span Int. at 65';
run;
/* AIC: 25497.4 */

/* Memory Quadratic */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tmemod ne .;
model tmemod= cage65 cagesq65 /ddfm=bw solution;
random intercept cage65 cagesq65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tmemod_EB_q;
title 'M1: Tmemo quadratic model, digit span Int. at 65';
run;
/* AIC: 25502.6 */
/* > NOTE! Linear model has better fit */

/* Intercept*/
proc sort data=tmemod_EB_l; by twinnr; run;
data int_tmemod;
  set tmemod_EB_l (where=(effect="Intercept"))
       ;
  by twinnr;
  qI65_tmemod=estimate;
  qI65_tmemodse=StdErrPred;
  qI65_tmemodse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI65_tmemod; 
run;

/* Linear slope */
data tmemod_EB_lin;
 merge tmemod_EB_l (where=(effect="cage65"))
        int_tmemod;
  by twinnr;
  qL65_tmemod =estimate;
  qL65_tmemodse=StdErrPred;
  qL65_tmemodse_inv=1/StdErrPred;

drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL65_tmemod; 
run;

/***********************************************************************************************/
/*Spatial*/

/* Spatial Linear */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tspat ne .;
model tspat= cage65 /ddfm=bw solution;
random intercept cage65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tspat_EB_l;
title 'M1: Tspat linear model Int. at 65';
run;
/* AIC: 21426.6 */

/* Spatial Linear */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tspat ne .;
model tspat= cage65 cagesq65 /ddfm=bw solution;
random intercept cage65 cagesq65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tspat_EB_q;
title 'M1: Tspat quadratic model Int. at 65';
run;
/* AIC: 21365.1 */
/* > Quadratic model has better fit */

/* Intercept*/
proc sort data=tspat_EB_q; by twinnr; run;
data int_tspat;
  set tspat_EB_q (where=(effect="Intercept"))
       ;
  by twinnr;
  qI65_tspat=estimate;
  qI65_tspatse=StdErrPred;
  qI65_tspatse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI65_tspat; 
run;

/* Linear slope */
data D;
 merge tspat_EB_q (where=(effect="cage65"))
        int_tspat;
  by twinnr;
  qL65_tspat=estimate;
  qL65_tspatse=StdErrPred;
  qL65_tspatse_inv=1/StdErrPred;

drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL65_tspat; 
run;

/* Quadratic slope */
data tspat_EB_quad;
  merge tspat_EB_q (where=(effect="cagesq65"))
        D;
  by twinnr;
  qQ65_tspat=estimate;
  qQ65_tspatse=StdErrPred;
  qQ65_tspatse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means data=tspat_EB_quad;
var  qQ65_tspat; 
run;

/***********************************************************************************************/
/*General ability*/

/*General ability Linear */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tpcom ne .;
model tpcom= cage65 /ddfm=bw solution;
random intercept cage65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tpcom_EB_l;
title 'M1: Tpcom linear model Int. at 65';
run;
/* AIC: 19141.3 */

/*General ability Quadratic */
proc mixed method=ml covtest data=SATSA_cog IC; class twinnr;
where age ge 50 and dem ne 1 and tpcom ne .;
model tpcom= cage65 cagesq65 /ddfm=bw solution;
random intercept cage65 cagesq65 /subject=twinnr type=un g gcorr solution;
ods output solutionR=tpcom_EB_q;
title 'M1: Tpcom quadratic model Int. at 65';
run;
/* AIC: 18945.0 */
/* > Quadratic model has better fit */

/* Intercept*/
proc sort data=tpcom_EB_q; by twinnr; run;
data int_tpcom;
  set tpcom_EB_q (where=(effect="Intercept"))
       ;
  by twinnr;
  qI65_tpcom=estimate;
  qI65_tpcomse=StdErrPred;
  qI65_tpcomse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means;
var qI65_tpcom; 
run;

/* Linear slope */
data E;
 merge tpcom_EB_q (where=(effect="cage65"))
        int_tpcom;
  by twinnr;
  qL65_tpcom=estimate;
  qL65_tpcomse=StdErrPred;
  qL65_tpcomse_inv=1/StdErrPred;

drop effect estimate StdErrPred--Probt;
run;
proc means;
var qL65_tpcom; 
run;

/* Quadratic slope */
data tpcom_EB_quad;
  merge tpcom_EB_q (where=(effect="cagesq65"))
        E;
  by twinnr;
  qQ65_tpcom=estimate;
  qQ65_tpcomse=StdErrPred;
  qQ65_tpcomse_inv=1/StdErrPred;
  drop effect estimate StdErrPred--Probt;
run;
proc means data=tpcom_EB_quad;
var  qQ65_tpcom; 
run;

/***********************************************************************************************/
/* Merge all together */
data SATSA_Int_L_Q_EB; 
merge tverb_EB_quad tsped_EB_quad tmemot_EB_quad tmemod_EB_lin tspat_EB_quad tpcom_EB_quad;
by twinnr;
run;

/***********************************************************************************************/

/* Saving the data */
data eb_cog.EB_cognition_all_20191216;
	set SATSA_Int_L_Q_EB;
run;

/* As text file */
PROC EXPORT DATA= WORK.SATSA_Int_L_Q_EB 
            OUTFILE= "P:\Dementia_IK\Cognition\EB_estimates\Data\EB_cognition_all_20191216.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=YES;
RUN;

/*************************************************************************************************/
/************************************** END OF FILE **********************************************/
/*************************************************************************************************/

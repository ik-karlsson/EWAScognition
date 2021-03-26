/***********************************************************************************************
# Name: Cognition_EWAS_GrowthCurve_20191126

# Purpose: Run longitudinal analysis of methylation and cognition, selected CpG sites with p<10^-5 in EWAS
# Study: Cognition EWAS
# NOTE! Previously conducted in R, but switch to SAS for convergence issues

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20191126 by Ida Karlsson (IK) 
# Updated:	20200114 by IK - New cognitive data and new CpG list

# Quadratic growth models are applied to Verbal ability, processing speed, spatial ability,
# general cognitive ability, and episodic memory/Thurstone memory test
# Linear growth models are applied to working memory/digit span test

# Input: Files for each cognitive domain in P:\Dementia_IK\Cognition\EWAS\Output, 
# created by Cognition_EWAS_CpG_followup_20181217.R by selecting all CpG sites with p<10^-5 in EWAS 

# Output: For each of the cognitive domains, one table is generated for the main results (beta,
# SE, and p-value for the effect of each CpG on the intercept and slope) and one table for the
# full output (fixed effects, random effects, and model fit info)

************************************************************************************************/
/***********************************************************************************************/
/* Import the cognitive data */
PROC IMPORT OUT= WORK.cog_adata
            DATAFILE= "P:\Dementia_IK\Cognition\Cognitive_data\Cognition_data_long_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=3734;

proc means data=cog_adata;
   var tverb tsped tspat tpcom tmemot tmemod;
run;

proc univariate data=cog_adata;
   var tverb tsped tspat tpcom tmemot tmemod;
   histogram;
run;

/***********************************************************************************************/
									/* 1. Verbal ability */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_70_tverb_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tverb ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-70)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;
*n=2833;

/* Null model */
/* NOTE! Convergence issues (boutndaries) with the quadratic in the random effects for pairid */
/* Hence, dropping the quadratic random effect.. */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tverb = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(8) cg17662953--cg11828983;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tverb = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 CpG_Value*cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG;
		ods output solutionF = fixed_CpG;
        ods output covparms = random_CpG;
run;

/* LRT */
data lrt;
	set fit_cpg;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear CpG_quadratic;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
	if Effect='CpG_Value*cage2_10' then output CpG_quadratic;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			CpG_quadratic (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_quadratic rename=StdErr=SE_CpG_quadratic rename=Probt=P_CpG_quadratic)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
	format P_CpG_quadratic e9.;
run;
	
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tverb.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtverb.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
									/* 2. Processing speed */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_65_tsped_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tsped ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-65)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;

/* Null model */
/* NOTE! Convergence issues (boutndaries) with the quadratic in the random effects for pairid */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tsped = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10/ subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(23) cg00978910--cg05404912;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tsped = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 CpG_Value*cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG;
		ods output solutionF = fixed_CpG;
        ods output covparms = random_CpG;
run;

/* LRT */
data lrt;
	set fit_cpg;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear CpG_quadratic;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
	if Effect='CpG_Value*cage2_10' then output CpG_quadratic;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			CpG_quadratic (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_quadratic rename=StdErr=SE_CpG_quadratic rename=Probt=P_CpG_quadratic)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
	format P_CpG_quadratic e9.;
run;
		
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tsped.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtsped.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
									/* 3. Spatial ability */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_65_tspat_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tspat ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-65)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;

/* Null model */
/* NOTE! Convergence issues (boutndaries) with the quadratic in the random effects for pairid */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tspat = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10/ subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(18) cg13484946--cg18098089;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tspat = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 CpG_Value*cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG_pre;
		ods output solutionF = fixed_CpG_pre;
        ods output covparms = random_CpG_pre;
run;
/* "Estimated G matrix is not positive definite" for cg09605693 and cg18833907. Will be dropped here modelled separately. */

data fit_CpG_1;
	set fit_CpG_pre;
	if CpG ne 'cg09605693' and CpG ne 'cg18833907';
data fixed_CpG_1;
	set fixed_CpG_pre;
	if CpG ne 'cg09605693' and CpG ne 'cg18833907';
data random_CpG_1;
	set random_CpG_pre;
	if CpG ne 'cg09605693' and CpG ne 'cg18833907';
run;

/* LRT */
data lrt_1;
	set fit_cpg_1;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
run;

/* Repeat without random quadratic effet for the two CpGs*/
/* Null model */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tspat = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null_2;
		ods output solutionF = fixed_null_2;
        ods output covparms = random_null_2;
run;

data adata_t_2;
	set adata_t;
	if CpG in ('cg09605693', 'cg18833907');
run;

proc mixed data=adata_t_2 method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tspat = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG_2;
		ods output solutionF = fixed_CpG_2;
        ods output covparms = random_CpG_2;
run;

/* LRT */
data lrt_2;
	set fit_cpg_2;
	if _n_ eq 1 then do;
		set fit_null_2 (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Add the two outputs together */
data lrt;
	set lrt_1 lrt_2;
proc sort;
	by CpG;
data fit_CpG;
	set fit_CpG_1 fit_CpG_2;
proc sort;
	by CpG;
data fixed_CpG;
	set fixed_CpG_1 fixed_CpG_2;
proc sort;
	by CpG;
data random_CpG;
	set random_CpG_1 random_CpG_2;
proc sort;
	by CpG;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear CpG_quadratic;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
	if Effect='CpG_Value*cage2_10' then output CpG_quadratic;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			CpG_quadratic (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_quadratic rename=StdErr=SE_CpG_quadratic rename=Probt=P_CpG_quadratic)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
	format P_CpG_quadratic e9.;
run;
	
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tspat.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtspat.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
									/* 4. General cognition */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_65_tpcom_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tpcom ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-65)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;

/* Null model */
/* NOTE! Convergence issues (boutndaries) with the quadratic in the random effects for pairid */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tpcom = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10/ subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(23) cg00978910--cg13805262;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tpcom = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 CpG_Value*cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 cage2_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG;
		ods output solutionF = fixed_CpG;
        ods output covparms = random_CpG;
run;

/* LRT */
data lrt;
	set fit_cpg;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear CpG_quadratic;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
	if Effect='CpG_Value*cage2_10' then output CpG_quadratic;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			CpG_quadratic (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_quadratic rename=StdErr=SE_CpG_quadratic rename=Probt=P_CpG_quadratic)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
	format P_CpG_quadratic e9.;
run;
	
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tpcom.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtpcom.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
									/* 5. Thurstone memory */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_65_tmemot_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tmemot ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-65)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;

/* Null model */
/* NOTE! Convergence issues (boutndaries) with the quadratic in the random effects, both for twinnr and pairid */
/* Hence, dropping the quadratic random effect.. */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tmemot = cage_10 cage2_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;
/* NOTE! Estimated G matrix is not positive definite" when quadratic age is included as random effect on twinnr, hence dropped. */

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(15) cg12946880--cg22684041;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tmemot = CpG_Value cage_10 cage2_10 CpG_Value*cage_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept cage_10 / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG;
		ods output solutionF = fixed_CpG;
        ods output covparms = random_CpG;
run;

/* LRT */
data lrt;
	set fit_cpg;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
run;
	
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tmemot.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtmemot.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
									/* 5. Digit span memory */
/***********************************************************************************************/
/* Import the data */
PROC IMPORT OUT= WORK.meth_adata
            DATAFILE= "P:\Dementia_IK\Cognition\EWAS\Output\CpG_followup_mval_65_tmemod_20191216.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=535;
proc sort data=meth_adata;
	by twinnr;
run;

proc means data=cog_adata;
	var tmemot;
run;

proc sgplot data=cog_adata;
   scatter y=tmemot x=age;
run;

proc univariate data=cog_adata;
   var tmemod tmemot;
   histogram;
run;
proc sgplot data=cog_adata;
   scatter y=Age x=tmemod;
run;

/* Merge with the cognitive data*/
/* Center age and create binary smoking variable */
data adata;
	merge 	cog_adata (in=a)
			meth_adata (in=b);
	by twinnr;
	if a and b and tmemod ne .;
	/* Center age and divide by 10 to help with convergence*/
	cage_10=(age-65)/10;
	cage2_10=cage_10*cage_10;
	/* Binary smoke, sex, and chip */
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	sex_b=sex-1;
	if chip='EPIC' then chip_b=0;
	else if chip='450K' then chip_b=1;
run;

/* Null model */
/* NOTE! Linear model */
/* NOTE! Convergence issues (boutndaries) with the linear random effects for pairid */
proc mixed data=adata method=ml noclprint noinfo covtest ic;
class twinnr pairid;
model tmemod = cage_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_null;
		ods output solutionF = fixed_null;
        ods output covparms = random_null;
run;

/* Methylation models */
/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set adata;
array meth(44) cg05155472--cg00947744;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    CpG_Value = meth(i);
    output;
end;
proc sort data=adata_t;
	by CpG twinnr age;
run;

proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid;
model tmemod = CpG_Value cage_10 CpG_Value*cage_10 sex_b smoke_b chip_b/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
random intercept cage_10 / subject= pairid (TWINNR) type=UN G Gcorr ;
        ods output InfoCrit = fit_CpG;
		ods output solutionF = fixed_CpG;
        ods output covparms = random_CpG;
run;

/* LRT */
data lrt;
	set fit_cpg;
	if _n_ eq 1 then do;
		set fit_null (keep=Neg2LogLike Parms rename=Neg2LogLike=Null_Neg2LL rename=Parms=Null_Parms);
	end;
   	LLdiff=Null_Neg2LL-Neg2LogLike;
	dgf=Parms-Null_Parms;
	LRT=1-PROBCHI(LLdiff,dgf);   
	keep CpG LRT;
	format LRT e9.;
run;

/* Select out data for output table */
Data CpG_icpt CpG_linear;
	set fixed_cpg;
	if Effect='CpG_Value' then output CpG_icpt;
	if Effect='CpG_Value*cage_10' then output CpG_linear;
run;

/* Merge */
data CpG_results;
	merge 	CpG_icpt (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_icpt rename=StdErr=SE_CpG_icpt rename=Probt=P_CpG_icpt)
			CpG_linear (keep=CpG Estimate StdErr Probt rename=Estimate=Beta_CpG_linear rename=StdErr=SE_CpG_linear rename=Probt=P_CpG_linear)
			LRT;
	by CpG;
	format P_CpG_icpt e9.;
	format P_CpG_linear e9.;
run;
	
/* Export */
PROC EXPORT DATA= CpG_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_tmemod.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Full output table */
data CpG_full;
	retain CpG Effect CovParm Subject;
	set Fixed_CpG Random_cpG Fit_CpG Fixed_null Random_null Fit_null;
	format Probt e9.;
	format Probz e9.;
proc sort;
	by CpG;
run;

/* Export */
PROC EXPORT DATA= CpG_full
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_GrowthCurve_20200114_FULLtmemod.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
/***************************************** END OF FILE *****************************************/
/***********************************************************************************************/


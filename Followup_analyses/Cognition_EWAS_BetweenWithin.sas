/***********************************************************************************************
# Name: Cognition_EWAS_BetweenWithin_20200106

# Purpose: Run between-within models of methylation and cognition, selected CpG sites with 
#		   p<10^-5 in EWAS. CpG methylation as exposure, EB estimates as outcome.
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20200106 by Ida Karlsson (IK) 
# Updated:	

# Input: Files for each cognitive domain in P:\Dementia_IK\Cognition\EWAS\Output, 
# created by Cognition_EWAS_CpG_followup_20181217.R for all CpG sites with p<10^-5 in EWAS 

# Output: For each of the cognitive domains, one table is generated for the main results (beta,
# SE, and p-value for between and within effect of each CpG on the intercept and slope) 

************************************************************************************************/

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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI70_tverb=. and qL70_tverb=. and qQ70_tverb=. then delete;
run;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(8) cg17662953--cg11828983;
array p_meth(8) pairm_cg17662953--pairm_cg11828983;
array dp_meth(8) depm_cg17662953--depm_cg11828983;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg17662953--depm_cg11828983;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI70_tverb = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI70_tverbse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL70_tverb = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL70_tverbse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Quadratic slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qQ70_tverb = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qQ70_tverbse_inv;
        ods output InfoCrit = fit_CpG_Q;
		ods output solutionF = fixed_CpG_Q;
        ods output covparms = random_CpG_Q;
run;

/* Select out data for output table */
Data qQ_pairmean qQ_depmean;
	set fixed_cpg_Q;
	if Effect='pairm_CpG' then output qQ_pairmean;
	if Effect='depm_CpG' then output qQ_depmean;
run;

/* Merge */
data qQ_results;
	merge 	qQ_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qQ_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b) qQ_results (in=c);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
	if c then outcome='Quadratic slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tverb.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI65_tsped=. and qL65_tsped=. and qQ65_tsped=. then delete;
run;
*n=533;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(23) cg00978910--cg05404912;
array p_meth(23) pairm_cg00978910--pairm_cg05404912;
array dp_meth(23) depm_cg00978910--depm_cg05404912;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg00978910--depm_cg05404912;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI65_tsped = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI65_tspedse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL65_tsped = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL65_tspedse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Quadratic slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qQ65_tsped = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qQ65_tspedse_inv;
        ods output InfoCrit = fit_CpG_Q;
		ods output solutionF = fixed_CpG_Q;
        ods output covparms = random_CpG_Q;
run;

/* Select out data for output table */
Data qQ_pairmean qQ_depmean;
	set fixed_cpg_Q;
	if Effect='pairm_CpG' then output qQ_pairmean;
	if Effect='depm_CpG' then output qQ_depmean;
run;

/* Merge */
data qQ_results;
	merge 	qQ_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qQ_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b) qQ_results (in=c);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
	if c then outcome='Quadratic slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tsped.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

/***********************************************************************************************/
									/* 2. Thurstone memory */
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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI65_tmemot=. and qL65_tmemot=. and qQ65_tmemot=. then delete;
run;
*n=533;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(15) cg12946880--cg22684041;
array p_meth(15) pairm_cg12946880--pairm_cg22684041;
array dp_meth(15) depm_cg12946880--depm_cg22684041;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg12946880--depm_cg22684041;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI65_tmemot = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI65_tmemotse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL65_tmemot = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL65_tmemotse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Quadratic slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qQ65_tmemot = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qQ65_tmemotse_inv;
        ods output InfoCrit = fit_CpG_Q;
		ods output solutionF = fixed_CpG_Q;
        ods output covparms = random_CpG_Q;
run;

/* Select out data for output table */
Data qQ_pairmean qQ_depmean;
	set fixed_cpg_Q;
	if Effect='pairm_CpG' then output qQ_pairmean;
	if Effect='depm_CpG' then output qQ_depmean;
run;

/* Merge */
data qQ_results;
	merge 	qQ_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qQ_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b) qQ_results (in=c);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
	if c then outcome='Quadratic slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tmemot.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

/***********************************************************************************************/
									/* 4. Spatial ability */
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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI65_tspat=. and qL65_tspat=. and qQ65_tspat=. then delete;
run;
*n=533;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(18) cg13484946--cg18098089;
array p_meth(18) pairm_cg13484946--pairm_cg18098089;
array dp_meth(18) depm_cg13484946--depm_cg18098089;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg13484946--depm_cg18098089;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI65_tspat = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI65_tspatse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL65_tspat = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL65_tspatse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Quadratic slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qQ65_tspat = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qQ65_tspatse_inv;
        ods output InfoCrit = fit_CpG_Q;
		ods output solutionF = fixed_CpG_Q;
        ods output covparms = random_CpG_Q;
run;

/* Select out data for output table */
Data qQ_pairmean qQ_depmean;
	set fixed_cpg_Q;
	if Effect='pairm_CpG' then output qQ_pairmean;
	if Effect='depm_CpG' then output qQ_depmean;
run;

/* Merge */
data qQ_results;
	merge 	qQ_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qQ_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b) qQ_results (in=c);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
	if c then outcome='Quadratic slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tspat.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

/***********************************************************************************************/
									/* 5. General cognition */
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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI65_tpcom=. and qL65_tpcom=. and qQ65_tpcom=. then delete;
run;
*n=533;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(23) cg00978910--cg13805262;
array p_meth(23) pairm_cg00978910--pairm_cg13805262;
array dp_meth(23) depm_cg00978910--depm_cg13805262;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg00978910--depm_cg13805262;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI65_tpcom = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI65_tpcomse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL65_tpcom = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL65_tpcomse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Quadratic slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qQ65_tpcom = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qQ65_tpcomse_inv;
        ods output InfoCrit = fit_CpG_Q;
		ods output solutionF = fixed_CpG_Q;
        ods output covparms = random_CpG_Q;
run;

/* Select out data for output table */
Data qQ_pairmean qQ_depmean;
	set fixed_cpg_Q;
	if Effect='pairm_CpG' then output qQ_pairmean;
	if Effect='depm_CpG' then output qQ_depmean;
run;

/* Merge */
data qQ_results;
	merge 	qQ_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qQ_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b) qQ_results (in=c);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
	if c then outcome='Quadratic slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tpcom.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

/***********************************************************************************************/
								/* 6. Digit span memory */
								/* NOTE! No quadratic! */
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

/* Modify variables */
data meth_adata2;
	set meth_adata;
	sex_b=sex-1;
	if smoke in(1,2) then smoke_b=0;
	else if smoke=3 then smoke_b=1;
	/* Drop if missing the outcome(s) */
	if qI65_tmemod=. and qL65_tmemod=. then delete;
run;
*n=535;

/* Transpose to just one M-value column to use in the by statement */
data adata_t;
set meth_adata2;
array meth(44) cg05155472--cg00947744;
array p_meth(44) pairm_cg05155472--pairm_cg00947744;
array dp_meth(44) depm_cg05155472--depm_cg00947744;
do i = 1 to dim(meth);
    CpG = vname(meth(i)); 
    pairm_CpG = p_meth(i);
	depm_CpG = dp_meth(i);
    output;
end;
drop cg05155472--depm_cg00947744;
proc sort data=adata_t;
	by CpG twinnr;
run;

/***********************************************************************************************/
/* Intercept */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qI65_tmemod = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qI65_tmemodse_inv;
        ods output InfoCrit = fit_CpG_I;
		ods output solutionF = fixed_CpG_I;
        ods output covparms = random_CpG_I;
run;

/* Select out data for output table */
Data qI_pairmean qI_depmean;
	set fixed_cpg_I;
	if Effect='pairm_CpG' then output qI_pairmean;
	if Effect='depm_CpG' then output qI_depmean;
run;

/* Merge */
data qI_results;
	merge 	qI_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qI_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Linear slope */
proc mixed data=adata_t method=ml noclprint noinfo covtest ic;
by CpG;
class twinnr pairid chip;
model qL65_tmemod = pairm_CpG depm_CpG age sex_b smoke_b chip/ CHISQ S;
random intercept / subject= pairid type=UN G Gcorr ;
weight qL65_tmemodse_inv;
        ods output InfoCrit = fit_CpG_L;
		ods output solutionF = fixed_CpG_L;
        ods output covparms = random_CpG_L;
run;

/* Select out data for output table */
Data qL_pairmean qL_depmean;
	set fixed_cpg_L;
	if Effect='pairm_CpG' then output qL_pairmean;
	if Effect='depm_CpG' then output qL_depmean;
run;

/* Merge */
data qL_results;
	merge 	qL_pairmean (keep=CpG Estimate StdErr Probt rename=Estimate=pairmean_Beta rename=StdErr=pairmean_SE rename=Probt=pairmean_P)
	 		qL_depmean (keep=CpG Estimate StdErr Probt rename=Estimate=depmean_Beta rename=StdErr=depmean_SE rename=Probt=depmean_P);
	by CpG;
	format pairmean_P 10.9;
	format depmean_P 10.9;
run;

/***********************************************************************************************/
/* Stack together */
data all_results;
	format CpG $15.;
	format outcome $15.;
	set qI_results (in=a) qL_results (in=b);
	if a then outcome='Intercept      ';
	if b then outcome='Linear slope';
proc sort;
	by CpG outcome;
run;
	
/* Export */
PROC EXPORT DATA= all_results
	OUTFILE= "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_tmemod.txt" 
    DBMS=TAB REPLACE;
RUN;

/* Remove all work files */
proc datasets library=work kill;
quit;

/***********************************************************************************************/
/***************************************** END OF FILE *****************************************/
/***********************************************************************************************/

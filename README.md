# EWAScognition
Epigenome-wide association study of level and change in cognitive abilities from midlife through late-life

By: Ida Karlsson, 20210326

------------------------------------------------------------------------------------------------------------------
This document provides explanations regarding the codes used in the paper "Epigenome-wide association study of 
level and change in cognitive abilities from midlife through late-life"

Preprint at Research Square: https://doi.org/10.21203/rs.3.rs-167499/v1

P.I: Ida Karlsson

Address: Department of Medical Epidemiology and Biostatistics, Karolinska Institutet, 171 77 Stockholm
        
E-mail: Ida.Karlsson@ki.se

------------------------------------------------------------------------------------------------------------------
ABSTRACT
------------------------------------------------------------------------------------------------------------------
Background: Epigenetic mechanisms are important in aging and may be involved in late-life changes in cognitive 
abilities. We conducted an epigenome-wide association study of leukocyte DNA methylation in relation to level and 
change in cognitive abilities, from midlife through late-life in 535 Swedish twins.

Results: Methylation levels were measured with the Infinium Human Methylation 450K or Infinium MethylationEpic 
array, and all sites passing quality control on both arrays were selected for analysis (n=250 816). Empirical 
Bayes estimates of individual intercept (age 65), linear, and quadratic change were obtained from latent growth 
curve models of cognitive traits and used as outcomes in linear regression models. Significant sites (p<2.4×10-7) 
were followed-up in between-within twin pair models adjusting for familial confounding, and full growth modelling. 
We identified six significant associations between DNA methylation and level of cognitive abilities at age 65: 
cg18064256 (PPP1R13L) with processing speed and spatial ability; cg04549090 (NRXN3) with spatial ability; 
cg09988380 (POGZ), cg25651129 (-), and cg08011941 (ENTPD8) with working memory. The genes are involved in neuro-
inflammation, neuropsychiatric disorders, and ATP metabolism. Within-pair associations were approximately half that 
of between-pair associations across all sites. In full growth curve models, associations between DNA methylation 
and cognitive level at age 65 were of small effect sizes, and associations between DNA methylation and longitudinal 
change in cognitive abilities of very small effect sizes.

Conclusions: Leukocyte DNA methylation was associated with level, but not change in cognitive abilities. The 
associations were substantially attenuated in within pair analyses, indicating they are influenced in part by 
genetic factors.

------------------------------------------------------------------------------------------------------------------
FOLDERS AND CODES
------------------------------------------------------------------------------------------------------------------

---------- Data_preprocessing ----------

Codes for all data preprocessing, described under Methods/DNA methylation measurements.
Performed by Yunzhang Wang.

------------------------------------------------------------------------------------------------------------------
---------- EWAS ----------

Code to compile data for and carry out EWAS analyses, the significant and suggestive findings of which are 
presented in Table 1 and Additional file 1.
Performed by Ida Karlsson.

--- 1: Cognition_EWAS_analysisdata.r 

Prepares the data for analyses

--- 2: Cognition_EWAS_analysis.r 

Conducts the epigenome-wide analyses of DNA methylation and EB-estimates of level and change in cognition.
Results are described in Results/EWAS of empirical Bayes estimates for level and change in cognitive abilities, 
with significant (p<2.4e-7) findings presented in Table 1 and suggestive (p<10e-5) findings in Additional file 2 
(tables created in Followup_analyses step 5 below).

------------------------------------------------------------------------------------------------------------------
---------- Followup_analyses ----------

Codes for all follow-up and supplementary analyses.
Performed by Ida Karlsson.

--- 1: Cognition_EWAS_CpG_followup.r

Selects all significant (p<2.4e-7) and suggestive (p<10e-5) CpGs from the EWAS runs for follow-up analyses.

--- 2: Extract_CpG_Mvalues.bash

Extract unstandardized M-values for significant CpGs, for Table S1 (created in step 3).

--- 3: Cognition_EWAS_descriptives_byCHIP.r

Creates descriptive statistics for full sample and by methylation array.
Results described in Results/ Study population, and presented in Additional file 1, Table S1.
Also checks for differences in associations for cg18064256 by chip, as the methylation levels significantly 
differ between the arrays.
Results described in Results/EWAS of empirical Bayes estimates for level and change in cognitive abilities.

--- 4: Cognition_EWAS_BetweenWithin.sas

Runs between-within (twin pair) models of methylation and cognition.
Results described in Results/ Between-within models of DNA methylation and empirical Bayes estimates for level and 
change in cognitive abilities, presented in Table 1 and Additional file 2.

--- 5: Cognition_EWAS_GrowthCurve.sas

Run latent growth-curve models of methylation and cognition.
Results described in Results/ Latent growth curve models of DNA methylation and level and change in cognitive 
abilities. A null model, without methylation, is presented in Additional file 1, Table S2, and associations with 
DNA methylation at significant and suggested sites are presented in Table 2 and Additional file 3, respectively.

--- 6: Cognition_EWAS_Table1.r

Generate Table 1 for the paper and the basis of Additional file 2 (based on output from the EWAS analyses and 
between-within analyses described above).

--- 7: Plot_trajectories.r

Generate Figure 1 for the paper (based on results from the latent growth-curve models, described above).

--- 8: Grep_mQTLs.bash

Extract data on the two suggested mQTLs (for step 9).

--- 9: Cognition_EWAS_mQTL_fup.r

Follow-up on mQTL effects, analyzing their association with DNA methylation and cognitive ability.
Results described in Results/ Characterization of the CpG sites.

--- 10: Cognition_EWAS_twincorr.r

Look at twin-pair correlation of significant CpG sites, and calculate Falcon heritability.
Results described in Results/ Characterization of the CpG sites.

-- 11: Cognition_EWAS_dementia_by_CpG.r

Run linear regression analyses between DNA methylation and dementia status to make sure associations are not just 
a result of cognitive decline in preclinical dementia.
Results described in Results/ Characterization of the CpG sites.

------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------

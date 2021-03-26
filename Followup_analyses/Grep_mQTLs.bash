# UNIX CODE FOR EXTRACTING 2 mQTLs FROM PSYCHCHIP

# CREATED BY IDA KARLSSON 20200210

# To identify the chunk of Chr in which the SNP of interest sits:
zgrep 'rs144382559' /projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chr10*.gz |cut -d' ' -f1-10
zgrep 'chrX:118976619:I' /projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chrX*.gz |cut -d' ' -f1-10

# Chunk files:
# psych.SE.qc15_chr10_47275278_52177976_impute2.imp.gz:--- rs144382559 50028963 GC G 0 1 0 1 0
# psych.SE.qc15_chrX_females_115667829_120575628_impute2.imp.gz:--- chrX:118976619:I 118976619 C CA 1 0 0 0.974 0.026
# psych.SE.qc15_chrX_males_115667829_120575628_impute2.imp.gz:--- chrX:118976619:I 118976619 C CA 1 0 0 1 0

# Extraction: 
# rs144382559
Adding the header and converting probs to A1 dosages
awk 'BEGIN {OFS="\t"} NR == 1 {printf "GENO_CHR\tSNP\tBP\tA1\tA2"} NR > 2 {printf("\t%s %s", $1, $2)} END {printf("\n")}' \
/projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chr10_47275278_52177976_impute2.imp_samples \
> extracted_rs144382559_A1dosages.txt

zcat /projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chr10_47275278_52177976_impute2.imp.gz \
| awk 'BEGIN {OFS="\t"} $2 == "rs144382559" {printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"); \
for (i = 6; i+2 <= NF; i+=3) {printf("%.3f%s", 2*$i + $(i+1), i+2 == NF ? "\n" : "\t")}}' >> extracted_rs144382559_A1dosages.txt

# chrX:118976619:I
# Females
Adding the header and converting probs to A1 dosages
awk 'BEGIN {OFS="\t"} NR == 1 {printf "GENO_CHR\tSNP\tBP\tA1\tA2"} NR > 2 {printf("\t%s %s", $1, $2)} END {printf("\n")}' \
/projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chrX_females_115667829_120575628_impute2.imp_samples \
> extracted_chrX118976619Ifem_A1dosages.txt

zcat /projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chrX_females_115667829_120575628_impute2.imp.gz \
| awk 'BEGIN {OFS="\t"} $2 == "chrX:118976619:I" {printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"); \
for (i = 6; i+2 <= NF; i+=3) {printf("%.3f%s", 2*$i + $(i+1), i+2 == NF ? "\n" : "\t")}}' >> extracted_chrX118976619Ifem_A1dosages.txt

# Males
Adding the header and converting probs to A1 dosages
awk 'BEGIN {OFS="\t"} NR == 1 {printf "GENO_CHR\tSNP\tBP\tA1\tA2"} NR > 2 {printf("\t%s %s", $1, $2)} END {printf("\n")}' \
/projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chrX_males_115667829_120575628_impute2.imp_samples \
> extracted_chrX118976619Imale_A1dosages.txt

zcat /projects/AGE/twins.psychchip.data/Data/DerivedData/imputation/1kGp1v3/chunks/psych.SE.qc15_chrX_males_115667829_120575628_impute2.imp.gz \
| awk 'BEGIN {OFS="\t"} $2 == "chrX:118976619:I" {printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"); \
for (i = 6; i+2 <= NF; i+=3) {printf("%.3f%s", 2*$i + $(i+1), i+2 == NF ? "\n" : "\t")}}' >> extracted_chrX118976619Imale_A1dosages.txt


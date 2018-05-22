#UKBB_full_QC
#
#pwd:/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC
#
#
#The fam,sample & .dat(relatedness) files for OA are here:
#/lustre/scratch115/projects/ukbiobank_oa/ukb997_cal_chr1_v2_s488374.fam
#/lustre/scratch115/projects/ukbiobank_oa/ukb997_imp_chr1_v2_s487406.sample
#/lustre/scratch115/projects/ukbiobank_oa/ukb997_rel_s488374.dat
#
#The dir typed data are here:
#/lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_cal_chr"1-22,X,Y,XY,MT"_v2.bed
#/lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_snp_chr"1-22,X,Y,XY,MT"_v2.bim
#
#The sample lists for the SR, HD and OP data are here:
#/nfs/humgen01/projects/ukbiobank_oa/for_eleni/July2017extract/SampleLists
#
#The file with all QC metrics is here:
#/lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_sqc_v2.txt
#A description of the genetic data types is here:
#http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531
#I have to paste the IDs from the fam file.
paste <(cut -d " " -f 1,2 /lustre/scratch115/projects/ukbiobank_oa/ukb997_cal_chr1_v2_s488374.fam) /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_sqc_v2.txt  > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/paste_famIDs_to_ukb_sqc_v2.txt
awk 'OFS="\t"{$1=$1} {print $0}' paste_famIDs_to_ukb_sqc_v2.txt > paste_famIDs_to_ukb_sqc_v2_tab.txt
#These are the agreed steps on the QC:
#1.	Remove samples with permission withdrawal. 
#2.	Sex mismatch (compare the self-reported sex and inferred sex data fields)
#3.	Outliers in heterozygosity and missing rates (a list of 968 samples  with unusually high heterozygosity and/or missing rates is provided; for heterozygosity they applied an outlier detection algorithm aberrant and for the missing rate they used a threshold of 5%)
#4.	Ethnicity outliers (a list of 409703 samples who self-reported 'White British' and have very similar genetic ancestry based on a principal components analysis of the genotypes)
#5.	Relatedness (in contrast with the Ethnicity step, here we have some alternative ways of maximizing the number of cases)
#	i.	Exclude 1165 samples (977 samples that don’t have a relatedeness score as they were excluded prior to the analyses due to KING’s strict criteria for missingness and heterozygosity and 188 individuals who have more than 10 putative third-degree relatives in the kinship table).
#	ii.	Use the list with the pairs of individuals related up to the third degree in the data set and:
#		1.	In case-case scenario keep the one with the lowest missingness. We can further keep some more cases by finding the individuals who are part of families,cousins etc and exclude for example the child rather than both parents but it will need even more time.
#		2.	In case-control scenario (always) keep the case.
#		3.	In control-control scenario keep the one with the lowest relatedness and/or missingness.
#		4.	In no-case/control Vs case keep (always) the case.
#		5.	In no-case/control Vs control keep the control.
#		6.	In no-case/control Vs no-case/control  keep the individual with the lowest missingness.
##Create the exclusion list for each of the above steps.
#Sex Mismatch
awk '$12!=$13' paste_famIDs_to_ukb_sqc_v2_tab.txt | head
awk '$12!=$13' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l
#378
awk '$12!=$13' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2,12,13 > Sex-mismatch_378.txt
awk '$12!=$13' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2 > Sex-mismatch_378_plink.txt
#fgrep -c -w -f OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt Sex-mismatch_378_plink.txt
#Het-Missing
#In	total	we	identified	968 samples	with	unusually	high	heterozygosity	and/or	missing	rates.
awk '{print $21}' paste_famIDs_to_ukb_sqc_v2_tab.txt | head
awk '$21==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l
#968
awk '$21==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2,21 > Het-Missing-outliers_968.txt
awk '$21==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2 > Het-Missing-outliers_968_plink.txt
#fgrep -c -w -f OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt  Het-Missing-outliers_968_plink.txt | wc -l
#128
#fgrep -w -f OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt  Het-Missing-outliers_968_plink.txt > HD-Het-Missing.txt
#White British Ancestry
awk '$26==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | head
awk '$26==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l
#409703
awk '$26==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2,26 > British-ancestry_409703.txt
awk '$26==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2 > British-ancestry_409703_plink.txt
#fgrep -c -w -f OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt  British-ancestry_409703_plink.txt
#50621
#fgrep -w -fOA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt  British-ancestry_409703_plink.txt > HD-British-ancestry.txt
#Unrelated
#1a.977 samples that don’t have a relatedeness score as they were excluded prior to the analyses due to KING’s strict criteria for missingness and heterozygosity
awk '$24==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | less -S
awk '$24==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l
#977
awk '$24==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut -f 1,2,24 > Excluded-from-kinship_977.txt
awk '$24==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut -f 1,2 > Excluded-from-kinship_977_plink.txt
#1b. 188 individuals who have more than 10 putative third-degree relatives in the kinship table
awk '$25==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | less -S
awk '$25==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l 
#188
awk '$25==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut -f 1,2,25 > Excess-relatives_188.txt
awk '$25==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut -f 1,2 > Excess-relatives_188_plink.txt
cat Excess-relatives_188_plink.txt Excluded-from-kinship_977_plink.txt | sort -u > Excess-relatives_188_plus_Excluded-from-kinship_977_total-1165_uniqsamples.txt
#We didn't use this list! A list of samples which were used to calculate the PCs, which is a (maximal) subset of unrelated participants after applying some QC filtering.
awk '{print $27}' paste_famIDs_to_ukb_sqc_v2_tab.txt | head
awk '$27==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | wc -l
#407219
awk '$27==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2,27 > Unrelated_407219.txt
awk '$27==1' paste_famIDs_to_ukb_sqc_v2_tab.txt | cut  -f 1,2 > Unrelated_407219_plink.txt
#fgrep -w -f OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt Unrelated_407219_plink.txt > HD-Unrelated.txt
#comm -1 -2 <(sort OA_ICDcodes_IDlist_AllOA_withdrawn_samples.txt) <(sort Unrelated_407219_plink_forgrep.txt) | wc -l
#comm -1 -2 <(sort HD-Unrelated.txt) <(sort HD-British-ancestry.txt) > comm_HD_Unrelated-British.txt
##Start QC
#1.	Remove samples with permission withdrawal (mv to /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Samples_with_permission_withdrawal_01082017)
for i in {1..22} X Y XY MT; do plink --bed /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_cal_chr"$i"_v2.bed --bim /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_snp_chr"$i"_v2.bim --fam /lustre/scratch115/projects/ukbiobank_oa/ukb997_cal_chr1_v2_s488374.fam --remove withdrawn_samples_26072017_plink.txt --make-bed --out UKBB_fr_rm_withdrawn-samples_chr"$i";done
#Missingness (mv to /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Missingness)
for i in {1..22} X Y XY MT; do plink --bfile UKBB_fr_rm_withdrawn-samples_chr"$i" --missing --out UKBB_fr_rm_withdrawn-samples_chr"$i"_missing ; done
#Updating pheno (since it's only sample QC from precalculated values, we will use only chr22)
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do mkdir $i ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do cp $i.txt Controls.txt /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk 'OFS="\t"{print $1,$1,1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/$i.txt > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_pheno.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk 'OFS="\t"{print $1,$1,0}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/Controls.txt > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/Controls_pheno.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do cat /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_pheno.txt /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/Controls_pheno.txt > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"-Controls_pheno.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/UKBB_fr_rm_withdrawn-samples_chr22 --make-pheno /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"-Controls_pheno.txt 1 --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22 ; done
#2.	Sex mismatch
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22 --remove /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Sex-mismatch_378_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex ; done
#3.	Outliers in heterozygosity and missing rates
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex --remove /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Het-Missing-outliers_968_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss ; done
#4.	Ethnicity outliers
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss --keep /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/British-ancestry_409703_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity ; done
#5.	Relatedness
# I checked what's the number of cases if I keep the lagest set of unrelated indi (provided by ukbb)
# for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity --keep /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Unrelated_407219_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_unrelated ; done
# We finally decided to use library(igraph) as UKBB people did in their qc.
# Briefly, these are the steps taken:
# 1.Melt the matrix into long format, excluding pairs below a certain level of relatedness
# 2.Represent the resulting long matrix as a graph
# 3.Extract all connected subgraphs (families)
# 4.For each subgraph, extract all the possible maximal unrelated subsets of individuals
# 5.Among those, select only the subsets that maximise the selection criteria (here we want to maximise the number of cases, then of controls)
# 6.If several optimal independent sets exist, select one at random
# 7.Assemble the list of independent individuals
# We can gain a few thousands of cases by the following way:
#	i.	Exclude 1165 samples not used in PCA.
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity --remove /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Excess-relatives_188_plus_Excluded-from-kinship_977_total-1165_uniqsamples.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples ; done
# Create the lists of cases, controls & others.
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==1) print $1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_CONTROLS.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==2) print $1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_CASES.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==-9) print $1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_OTHERS.txt ; done
# Adding case/control/other status to the "Relatedness" file 
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do perl -e 'open(CASES, "/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/'$i'/'$i'_chr22_sex_het-miss_ethnicity_notPCAsamples_CASES.txt"); my @case=<CASES>; chomp(@case);close(CASES);open(CONTROLS, "/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/'$i'/'$i'_chr22_sex_het-miss_ethnicity_notPCAsamples_CONTROLS.txt"); my @control=<CONTROLS>; chomp(@control);close(CONTROLS);open(REL, "/lustre/scratch115/projects/ukbiobank_oa/ukb997_rel_s488374.dat");my $skip=<REL>;while(<REL>){chomp;my @line=split(/ /, $_);my $one=0; my $two=0;my $three=0; my $four=0;foreach my $cid (@case){if($line[0] eq $cid){$one=1;}if($line[1] eq $cid){$two=1}};foreach my $cid (@control){if($line[0] eq $cid){$three=1;}if($line[1] eq $cid){$four=1}}; print "$one\t$two\t$three\t$four\t", join("\t", @line), "\n";}' > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/${i}_casecontrol;done
#or (it takes 2-3h for each pheno)
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do echo ~/get_cc.sh $i;done | ~ag15/array 2g get_cc
bsub -R"select[mem>2000 && green] rusage[mem=2000]" -M2000 -J "get_cc[1-8]" -q yesterday -o get_cc.%I.o -e get_cc.%I.e './tmp.100870.exe $LSB_JOBINDEX'
#Running the R script which creates the list of unrelated individuals to remove
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Scripts/Relatedness.R /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_casecontrol /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_relatedness_samplestoremove /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples.fam ; done
#Adding one column to the output file for plink
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{print $1,$1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_relatedness_samplestoremove > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_relatedness_samplestoremove_plink ; done
#Excluding the related samples
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples --remove /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_relatedness_samplestoremove_plink --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated; done
# Create the lists of cases, controls & others from the final qced samples.
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==1) print $1,$1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_CONTROLS.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==2) print $1,$1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_CASES.txt ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee; do awk '{if ($6==-9) print $1,$1}' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated.fam > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_OTHERS.txt ; done
##Checking relatedness through plink (in HD only)
#pwd /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/HD/Check-relatedness
#for i in {1..22};do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Samples_with_permission_withdrawal_01082017/UKBB_fr_rm_withdrawn-samples_chr"$i" --keep 266924-controls_47931-cases.txt --make-bed --out chr"$i"_266924-controls_47931-cases_HD; done
#for i in {1..22};do plink --bfile chr"$i"_266924-controls_47931-cases_HD  --make-pheno HD_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_CASES  '*' --make-bed --out chr"$i"_266924-controls_47931-cases_HD_pheno; done
#Extract only the markers that were used in the calculation of kinship
#for i in {1..22};do plink --bfile chr"$i"_266924-controls_47931-cases_HD_pheno  --extract Marker_was_used_in_the_calculation_of_kinship_93511snps.txt --make-bed --out chr"$i"_266924-controls_47931-cases_HD_pheno_kinship-snps; done
#Run genome
#for i in {1..22};do echo plink --bfile chr"$i"_266924-controls_47931-cases_HD_pheno_kinship-snps  --genome gz  --out chr"$i"_266924-controls_47931-cases_HD_pheno_kinship-snps_genome; done | ~ag15/array 300g HD_rel
#bsub -q hugemem  -R"select[mem>300000 && red] rusage[mem=300000]" -M300000 -J "HD_rel[1-22]"  -o HD_rel.%I.o -e HD_rel.%I.e './tmp.954849.exe $LSB_JOBINDEX'
####END OF QC#########

#UKBiobank OA phenotype extraction
#pwd:/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
## Use Wills script to extract a set of phenotypes from the main phenotype file. 
### Use it like this to get the help message displayed:
perl /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Scripts/extract-v1.2.2.pl
#You will need the following files:
#phenotypefile : tab separated main phenotype file with path from current location, eg ukb8868.tab
#dctfile       : the corresponding dct header file with path, eg ukb8868.dct
#inputfilename : a file listing the columns you want to extract, one column per line
#outputfilename: the name of your output file (will be created) 
#samplefile    : only needed if you dont want the whole sample set (eg the initial genotyped ids)
#cp to the pwd the files needed
cp /nfs/humgen01/projects/ukbiobank_oa/Phenotypes/v2_Jan2017/ukb8868.tab /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
cp /nfs/humgen01/projects/ukbiobank_oa/Phenotypes/v2_Jan2017/ukb8868.dct /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
cp /nfs/humgen01/projects/ukbiobank_oa/Phenotypes/977_9377_Mar2017/ukb9377.tab /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
cp /nfs/humgen01/projects/ukbiobank_oa/Phenotypes/977_9377_Mar2017/ukb9377.dct /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
#example files listing phenotypes to extract:
cp /lustre/scratch115/projects/ukbiobank_oa/DataExtracts/for_eleni/columns_to_extract.txt /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/columns_to_extract_ioannas.txt
cp /nfs/team144/it3/UKBiobank_OA/code/columns_to_extract_UKBB_OA_BMI_Waist_Hip.txt /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/columns_to_extract_UKBB_OA_BMI_Waist_Hip.ioannas.txt
cp /nfs/humgen01/projects/ukbiobank_oa/DataExtracts/for_MR_Fernando/columns_to_extract_forMR_Fernando.txt /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction
# to get an overview of what phenotypes are available for each project please have a look at the *.html files 
/nfs/team144/Kostas/UKBB/ukb8868.html
# Creating different formats of the phenotype file
# Each phenotype file provided by UKBB needs to be processed in different steps before it can be used.   
# A detailed description of this process can be found here:
http://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=accessing_data_guide
# The phenotype files provided by UKBB for download have a special format and the ending ".enc". 
# These files need to be unpacked using the UKBB script "ukb_unpack" expecting the file name and 
# a special key file as parameters. The key file is provided by UKBB together with the phenotype 
# file and is specific to each phenotype file.
# The unpacked files have to be converted to different formats. 
# Available formats are: csv, docs, sas, stata, r and bulk; use this command:
# All the above steps have been done for the 1st release so no need to be repeated.
ukb_conv inputfile csv|docs|sas|stata|r|bulk [options]
eg: UKBB_downloads/ukb_conv ukb6152.enc_ukb stata
#UKBiobank OA phenotype extraction (Age, sex, 40PCs, BMI, Waist, Hip, Education)
#40 PCs & genotyping.array UKBL=UK BiLEVE and UKBB=UK Biobank(these are included in the sample qc file)
cut -f 1,2,5,28-67 /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/paste_famIDs_to_ukb_sqc_v2_tab.txt > 40PCs-GenArray.txt
#Age, sex, BMI, Waist, Hip and the traits we used for MR analysis
#perl /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Scripts/extract-v1.2.2.pl -p /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb8868.tab -m /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb8868.dct -i /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/columns_to_extract_UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip.txt -f /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb997_imp_chr1_v2_s487406.sample -o /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip.txt
perl /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Scripts/extract-v1.2.2.pl -p /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb8868.tab -m /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb8868.dct -i /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/columns_to_extract_forMR_Fernando.txt -f /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb997_imp_chr1_v2_s487406.sample -o /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR.txt
#Education attainment
perl /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Scripts/extract-v1.2.2.pl -p /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb9377.tab -m /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb9377.dct -i /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/columns_to_extract_educational-attainment.txt -f /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb997_imp_chr1_v2_s487406.sample -o /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/UKBB_OA_educational-attainment.txt
#Years of schooling (calculated according to Okbay et al, see blw ) 
Educational attainment measure(s)		ISCED mapping		EduYears Mean/SD
Which of the following qualifications do you have? (You can select more than one)?					
1 == College or University degree 
2 == A levels/AS levels or equivalent
3 == O levels/GCSEs or equivalent 
4 == CSEs or equivalent	5 == NVQ or HND or HNC or equivalent
6 == Other prof. qual. eg: nursing, teaching 
7 == None of the above 
8 == Prefer not to answer		
1) ISCED 5 (20 years)
2) ISCED 3
3) ISCED 2
4) ISCED 2	5) ISCED (19 years)
6) ISCED 4
7) ISCED 1
8) Excluded		13.7
(5.1)
			
Supplementary Table 1.2. Mapping from ISCED Level to EduYears and College
ISCED Level	Definition	US years of schooling (EduYears)	College
0	Pre-primary education	1	0
1	Primary education or first stage of basic education	7	0
2	Lower secondary or second stage of basic education	10	0
3	(Upper) secondary education	13	0
4	Post-secondary non-tertiary education	15	0
5	First stage of tertiary education (not leading directly to an advanced research qualification)	19	1
6	Second stage of tertiary education (leading to an advanced research qualification, e.g. a Ph.D.)	22	1
#Merging the two files
diff <(cut -f 1 UKBB_OA_educational-attainment.txt) <(cut -f 1 UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR.txt)
paste UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR.txt <(cut -f 2,3,4,5,6,7,8 UKBB_OA_educational-attainment.txt) > UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR_educational-attainment.txt
#I have to merge the two files containing traits/covariates (40PCs-GenArray.txt, UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR_educational-attainment.txt) with the sample file (ukb997_imp_chr1_v2_s487406.sample) by keeping the order of the samples.
## Remove the 2nd line of the sample file
awk 'NR==1 || NR>2' /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb997_imp_chr1_v2_s487406.sample > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/ukb997_imp_chr1_v2_s487406_no2ndline.sample
# Then, in R
library(data.table)
sample_file=fread("ukb997_imp_chr1_v2_s487406_no2ndline.sample")
to_add_file=fread("UKBB_OA_sex_age-at-attedance_BMI_Waist_Hip_plus-all-traits-for-MR_educational-attainment.txt")
both_merged=merge(sample_file, to_add_file, all.x=T, all.y=F, by.x="ID_1", by.y="Encoded_anonymised_participant_ID")
both_merged=both_merged[match(sample_file$ID_1, both_merged$ID_1),]
write.table(both_merged, "merged.txt", row.names=F, quote=F) 			
both_merged=fread("merged.txt")
pc_array=fread("40PCs-GenArray.txt")
all_merged=merge(both_merged, pc_array, by.x="ID_1", by.y="V1", all.x=T, all.y=F )
all_merged=all_merged[match(sample_file$ID_1, all_merged$ID_1),]
all_merged[,"V2":=NULL]
write.table(all_merged , "all_merged.txt", row.names=F, quote=F)
#After some small modifications, the final sample template file, where the only thing left to be added is the OA pheno, is this:
#/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/OA-template.sample
#
#We decided to run 2 GWAS. One by using all available controls and one with controls selected to be ~4x the number of cases (selecting the older individuals while keeping the number of males and females with the same balance as cases.
#I prepared fam files to pass to Eleni as she is going to run the latter GWAS (select controls)
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee ; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated --keep /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_CONTROLS.txt --pheno /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/Sex-Age_all-imputed-samples_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_Age-Sex_All-QCed-Controls ; done
for i in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee ; do plink --bfile /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated --keep /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_chr22_sex_het-miss_ethnicity_notPCAsamples_unrelated_CASES.txt --pheno /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/Phenotype_extraction/Sex-Age_all-imputed-samples_plink.txt --make-bed --out /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/QC/$i/"$i"_Age-Sex_All-QCed-Cases ; done

###SNPtest###
#pwd /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_all-controls
#pwd /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_4xcontrols
#Stripe the imputed files so that they can be read by SNPtest faster (http://mediawiki.internal.sanger.ac.uk/index.php/LustreStripeSize)
mkdir /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/Stripe
for i in {1..22} ;  do cp /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_imp_chr"$i"_v2.bgen /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/Stripe/ ; done
#Index the imputed files with bgenix
module add hgi/gcc/4.9.1
/software/team144/gavinband-bgen-f2083028f07e/bin/bgenix -g /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_imp_chr22_v2.bgen -list
#Finding variants with MAF>0.001 and info score>=0.4. Then create an inclusion snplist (for checking reason, SNPtest accepts only exclusions list). 
#for i in {1..22} ; do zcat /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_mfi_chr"$i"_v2.txt.gz | awk '{if (($5>0.001) && ($6 >= 0.4)){ print $0}}' > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_all-controls/Scripts/ukb_mfi_chr"$i"_v2_maf-gt0001_AND_infogte04.txt ; done
#Finding variants with MAF<=0.001 and info score<0.4. Then create exclusion snplist to be used by SNPtest.
for i in {1..22} ; do zcat /lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/ukb_mfi_chr"$i"_v2.txt.gz | awk '{if (($5<=0.001) || ($6<0.4)){ print $0}}' > /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_all-controls/Scripts/ukb_mfi_chr"$i"_v2_maf-lte0001_OR_infolt04.txt ; done
#Checking
for i in {1..22}; do comm -1 -2 <(awk '{print '$i'":"$2}' ukb_mfi_chr"$i"_v2_maf-gt0001_AND_infogte04.txt | sort) <(awk '{print $1":"$2}' HRC.r1-1.GRCh37.wgs.mac5.sites.tab | sort) | wc -l; done
for i in {1..22}; do fgrep -w -v -f <(awk '{print '$i'":"$2}' ukb_mfi_chr"$i"_v2_maf-lte0001_OR_infolt04.txt | sort) <(awk '{print $1":"$2}' HRC.r1-1.GRCh37.wgs.mac5.sites.tab | sort) > ukb_mfi_chr"$i"_v2_maf-lte0001_OR_infolt04_rmHRC.txt ; done
#Preparing the structure of the folders folders (in every pheno folder create 1..22 subfolders)
#for i in {1..22}; do for j in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee;do mkdir -p /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_all-controls/$j/chr"$i"; done ; done
for i in {1..22}; do for j in SR HD HD_hip HD_knee HD_hipknee OP_hip OP_knee OP_hipknee;do mkdir -p /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/SNPtest_4xcontrols/$j/chr"$i"; done ; done
#Running SNPtest (the snptest sh and sample file should be in the same folder, usage SNPTEST_chunk_command_with_cov_kh7.sh sample_file.sample name_of_phenotype)
./SNPTEST_chunk_command_with_cov_kh7.sh HD_4xcontrols.sample HD
#Cp all the commands into a single file (eg HD_bsubs_4xcontrols.sh)
cat HD_bsubs_4xcontrols.sh > 
dos2unix HD_bsubs_4xcontrols.sh
chmod 750 HD_bsubs_4xcontrols.sh
./HD_bsubs_4xcontrols.sh
#Checking the jobs:
watch bjobs -A
#Resubmiting failed jobs
for i in {1..1}; do failure=$(grep -L finito chr$i/*.o | sed 's/.o.*//;s/.*\///' | sort -u | tr '\n' ',' | sed 's/.$/\n/'); grep -L finito chr$i/*.o  |  sed 's/:.*//'| while read line; do echo mv $line $line.old; done; head -n $i /lustre/scratch115/projects/ukbiobank_oa/EZ2/eleni_HD/bsub_HD.sh | tail -1| sed 's/\[1-.*\]/\['$failure'\]/;s/long/basement/'; done | grep -v "\[\]"

##BOLT-LMM###

Creating sample/pheno files by excluding the UK cohorts duplicates (Duplicates-100UKHLS-30arco-corex-192arcOGEN_illumina610k.txt)
pwd: /lustre/scratch115/projects/ukbiobank_oa/EZ2/FULL_RELEASE/SNPtest_4xcontrols/SampleFiles
Replace in the existing sample files the disease status ***(1cases0controls-9missing)*** and create files likes this HD_hip_4xcontrols_forBOLT_0cases1controls-9missing.sample.txt

Putting -9 to all duplicates
cat HD_hip_4xcontrols_forBOLT_0cases1controls-9missing.sample.txt | perl -lane 'BEGIN{open(DUPS, "Duplicates-100UKHLS-30arco-corex-192arcOGEN_illumina610k.txt") or die ("no open file 1");our @dups=<DUPS>;chomp(@dups);}{foreach my $dup (@dups){if($F[0]==$dup){$F[$#F]=-9;last;}} print(join("\t", @F));}' > HD_hip_4xcontrols_forBOLT_0cases1controls-9missing.sample.noUKdups.txt
remember to do 2 things
a. Fix the first to columns to read FID IID for the pheno/covar file & ID_1 ID_2 for the sample file
b. find and fix the lines (usually 3 that have been corrupted)

Update the pheno in the fam file
awk '{print $1,$2,$NF}' HD_hip_4xcontrols_forBOLT_0cases1controls-9missing.sample.noUKdups.txt > HD_hip_bolt_noUKdups_fam-pheno.txt
(rm headers)
(remember to replace #N/A with -9)

Create a list of dir typed snps to be excluded by the model (I used a maf of 0.001 & missingness of 0.1 to be excluded as Po-Ru suggests in his example)
for i in {1..22}; do plink --bed /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_cal_chr"$i"_v2.bed --bim /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_snp_chr"$i"_v2.bim --fam /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/HD/ukb997_cal_chr1_v2_s488374_HD_FixCol6.fam --max-maf 0.001 --make-bed --out ukb_snp_chr"$i"_0001; done
cat ukb_snp_chr*_0001.bim | cut -f 2 > autosome_maf_lt_0.001.txt
(rm bed/bim/fam)
for i in {1..22}; do plink --bed /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_cal_chr"$i"_v2.bed --bim /lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_snp_chr"$i"_v2.bim --fam /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/HD/ukb997_cal_chr1_v2_s488374_HD_FixCol6.fam --missing --out ukb_snp_chr"$i"_missing; done
for i in {1..22}; do awk '$5>0.1' ukb_snp_chr"$i"_missing.lmiss; done> ukb_snp_chrALL_missing.lmiss.gt.0.1
awk '{print $2}' ukb_snp_chrALL_missing.lmiss.gt.0.1 > autosome_missing_gt_0.1.txt
(rm bed/bim/fam)

#Create a dataset of the unique cases between SR and HD dataset
#Cases
cat <(awk '$NF==0' /lustre/scratch115/realdata/mdt3/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/HD/HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt) <(awk '$NF==0' /lustre/scratch115/realdata/mdt3/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/SR/SR_4xcontrols_forBOLT_0cases1controls-9missing.sample.noUKdups.txt) | sort -u | cut -f 1 > SRplusHD_cases_68623.txt
#Controls
#I take the all controls files and rm UK duplicates
cat HD_500K_all-controls.sample | perl -lane 'BEGIN{open(DUPS, "Duplicates-100UKHLS-30arco-corex-192arcOGEN_illumina610k.txt") or die ("no open file 1");our @dups=<DUPS>;chomp(@dups);}{foreach my $dup (@dups){if($F[0]==$dup){$F[$#F]=-9;last;}} print(join("\t", @F));}' > HD_500K_all-controls.sample.noUKdups.txt
cat SR_500K_all-controls.sample | perl -lane 'BEGIN{open(DUPS, "Duplicates-100UKHLS-30arco-corex-192arcOGEN_illumina610k.txt") or die ("no open file 1");our @dups=<DUPS>;chomp(@dups);}{foreach my $dup (@dups){if($F[0]==$dup){$F[$#F]=-9;last;}} print(join("\t", @F));}' > SR_500K_all-controls.sample.noUKdups.txt
#Create sample list and keep the in common controls
awk '$NF==0' HD_500K_all-controls.sample.noUKdups.txt | cut -f 1 | cut -d " " -f 1 > HD_all-controls.noUKdups.266792.txt
awk '$NF==0' SR_500K_all-controls.sample.noUKdups.txt | cut -f 1 | cut -d " " -f 1 > SR_all-controls.noUKdups.269365.txt
comm -1 -2 <(sort SR_all-controls.noUKdups.269365.txt) <(sort HD_all-controls.noUKdups.266792.txt) > Incommon_247846controls_of_HD_and_SR.txt

per chr
(Arthur's way)
for i in {1..22}; do echo ./UKBB-full_run_BOLT-LMM_v2.3_perchr.sh; done | ~ag15/array 20g HD | sed 's/red/green/'
bsub -R"select[mem>20000 && green] rusage[mem=20000]" -M20000 -J "HD[1-22]"  -o HD.%I.o -e HD.%I.e './tmp.489390.exe $LSB_JOBINDEX'
bsub -G helic -n8 -q basement -R"span[hosts=1] select[mem>20000 && green] rusage[mem=20000]" -M20000 -J "HD[1-22]"  -o HD.%I.o -e HD.%I.e './tmp.489390.exe $LSB_JOBINDEX'

~ag15/local_programs/gsub 50G -n8 -R"span[hosts=1]" -q basement -G ukbiobank_oa -o SR_bolt_chr14-22.o -e SR_bolt_chr14-22.e ./UKBB-full_run_BOLT-LMM_v2.3_chr14-22.sh
#or Interactively in tmux session
~ag15/local_programs/gsub 50G -n8 -R"span[hosts=1]" -q basement -G ukbiobank_oa -I  ./UKBB-full_run_BOLT-LMM_v2.3_chr14-22.sh

(Ioanna's way; needs to be fixed)
bsub -G ukbiobank_oa -J"HD_UKBB_run_BOLT[1-22]" -q basement -n8 \
-o /lustre/scratch115/projects/ukbiobank_oa/EZ2/FULL_RELEASE/BOLT-LMM_v2.3/HD/HD_UKBB_run_BOLT_LMM.%J-%I.o \
-e /lustre/scratch115/projects/ukbiobank_oa/EZ2/FULL_RELEASE/BOLT-LMM_v2.3/HD/HD_UKBB_run_BOLT_LMM.%J-%I.e \
-M185000 -R'span[hosts=1] select[mem>'185000'] rusage[mem='185000']' /lustre/scratch115/projects/ukbiobank_oa/EZ2/FULL_RELEASE/BOLT-LMM_v2.3/HD/UKBB-full_run_BOLT-LMM_v2.3_perchr.sh

#QQ & Manhattan
#(in gen1a)
~ag15/scripts/man_qq_annotate --chr-col CHR --pos-col BP --auto-label --pval-col P_LINREG --title "P_LINREG-HD-mafgte001.infogt04" --sig-thresh 5e-08 --sig-thresh-line 5e-08 chrALL_bolt_HD-OA.bgen.stats.chr-pos.HRC.mafgte001.infogt04.tab.headers P_LINREG-HD-mafgte001.infogt04


#Summary of the results
#Checking outputs
grep -w Successfully OP_hip.*.o
#Create a column with chr:pos in the BOLT outputs
for i in {1..22}; do  zcat chr"$i"_bolt_OP_hip-OA.bgen.stats.gz | awk 'OFS="\t"{$17=$2":"$3; print $0;}' | gzip > chr"$i"_bolt_OP_hip-OA.bgen.stats.chr-pos.gz; done
#Check that everything went OK
zcat chr*_bolt_OP_hip-OA.bgen.stats.gz | wc -l
zcat chr*_bolt_OP_hip-OA.bgen.stats.chr-pos.gz | wc -l
#Extract the HRC variants
#(in gen 1, one by one the phenos cause crashes to the server)
for i in {1..22}; do LC_ALL=C  zgrep -w -F -f /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/Scripts/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr-pos.snplist  chr"$i"_bolt_HD-OA.bgen.stats.chr-pos.gz | gzip > chr"$i"_bolt_HD-OA.bgen.stats.chr-pos.HRC.gz; done
#(in farm)
~ag15/local_programs/gsub 20G -q normal -G helic -o HRC_grep_OP_hip.o -e HRC_grep_OP_hip.e ./HRC_grep.sh
#cat
zcat chr*_bolt_HD-OA.bgen.stats.chr-pos.HRC.gz | gzip > chrALL_bolt_HD-OA.bgen.stats.chr-pos.HRC.gz
#(add headers)
zcat headers-chr-bp.gz chrALL_bolt_HD_hip-OA.bgen.stats.chr-pos.HRC.gz | gzip  > chrALL_bolt_HD_hip-OA.bgen.stats.chr-pos.HRC.headers.gz
#(filtering for maf; they already contain variants with info > 0.3)
#QQ & Manhattan
#(in gen1a)
~ag15/scripts/man_qq_annotate --chr-col CHR --pos-col BP --auto-label --pval-col P_LINREG --title "P_LINREG-OP_hipknee.infogt03" --sig-thresh 5e-08 --sig-thresh-line 5e-08 chrALL_bolt_OP_hipknee-OA.bgen.stats.chr-pos.HRC.headers.gz P_LINREG-OP_hipknee.HRC.infogt03

#Summary of the results (use the HRC variants only)
MYPATH=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/HD/reformatted
echo "SNPIDwithAlleles SNPID rsID CHR POS EA NEA EAF beta SE pval_bolt-lmm pval_inf INFO MAF betalower betaupper OR ORlower ORupper" > ${MYPATH}/allchr.boltlmm.HD.results
zcat /lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/HD/chr*_bolt_HD-OA.bgen.stats.chr-pos.HRC.gz | grep -v SNP | awk '{if($14!="NA" && $14>0) print $2":"$3"_"$5"_"$6,"chr"$2":"$3,$1,$2,$3,$5,$6,$7,$11,$12,$16,$14,$8}' | awk '{if($8 > 0.5) MAF=(1-$8); else MAF=$8; print $0,MAF}' | awk '{print $0,$9-1.96*$10,$9+1.96*$10}' | awk '{if($9<50 && $9 > -50) {OR=exp($9);ORlower=exp($15);ORupper=exp($16)} else {OR="NA";ORlower="NA";ORupper="NA"}; print $0,OR,ORlower,ORupper}' >> ${MYPATH}/allchr.boltlmm.HD.results
#gzip allchr.boltlmm.HD.results
#Checking for dupliactes
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | awk '{if(NR>1) print $1}' | sort | uniq -d | wc -l
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | awk '{if(NR>1) print $2}' | sort | uniq -d | wc -l
#210265
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | wc -l
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 8 | head -1 # EAF
10:118804703_A_G chr10:118804703 rs182428624 10 118804703 A G 0 -inf -nan 1.0E+00 1.0E+00 1 0 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 8 | tail -1 # EAF
9:99999608_C_T chr9:99999608 rs113309985 9 99999608 C T 1 1.08377e+11 -nan 1.0E+00 1.0E+00 0.959697 0 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 9 | head -1 # beta
10:100597413_C_T chr10:100597413 rs188092515 10 100597413 C T 1 -inf -nan 1.0E+00 1.0E+00 0.973371 0 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 9 | tail -1 # beta
3:21701412_C_T chr3:21701412 rs144595958 3 21701412 C T 1 2.20767e+15 -nan 1.0E+00 1.0E+00 0.532344 0 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 10 | head -1 # se
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 10 | tail -1 # se
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 11 | head -1 # pval
1:176736924_C_T chr1:176736924 rs372068250 1 176736924 C T 0.999974 0.713715 0.11446 4.5E-10 4.6E-10 0.954915 2.6e-05 0.489373 0.938057 2.04156 1.63129 2.55501
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 11 | tail -1 # pval
9:99999811_G_A chr9:99999811 rs529410169 9 99999811 G A 0.999998 2.63559e+09 -nan 1.0E+00 1.0E+00 0.405107 2e-06 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 12 | head -1 # pval_noninf
1:176736924_C_T chr1:176736924 rs372068250 1 176736924 C T 0.999974 0.713715 0.11446 4.5E-10 4.6E-10 0.954915 2.6e-05 0.489373 0.938057 2.04156 1.63129 2.55501
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 12 | tail -1 # pval_noninf
9:99999811_G_A chr9:99999811 rs529410169 9 99999811 G A 0.999998 2.63559e+09 -nan 1.0E+00 1.0E+00 0.405107 2e-06 nan nan NA NA NA
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 13 | head -1 # INFO
11:49506260_T_C chr11:49506260 rs187233193 11 49506260 T C 0.993383 0.00151096 0.0128623 9.1E-01 9.1E-01 0.3 0.006617 -0.0236991 0.0267211 1.00151 0.97658 1.02708
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | sort -g -k 13 | tail -1 # INFO
9:99995514_G_A chr9:99995514 rs10981280 9 99995514 G A 0.548871 -0.000121841 0.00115149 9.2E-01 9.2E-01 1 0.451129 -0.00237876 0.00213508 0.999878 0.997624 1.00214

# check when I get -inf or inf beta & if filtered out thanks to other filters
# check when I get -nan or nan se & if filtered out thanks to other filters
# check if column1 in boltlmm starts with rs or Aff
# must fix it for snptest results at least
# then deal with duplicates snpids
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | awk '{if($9=="-inf" || $9=="inf"){print $0}}' > ${MYPATH}/allchr.boltlmm.HD.results.problematicbeta
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | awk '{if($10=="-nan" || $10=="nan"){print $0}}' > ${MYPATH}/allchr.boltlmm.HD.results.problematicse
wc -l ${MYPATH}/allchr.boltlmm.HD.results.problematicbeta
#19294
wc -l ${MYPATH}/allchr.boltlmm.HD.results.problematicse
#2296990 (all nan)
sort -g -k 10 ${MYPATH}/allchr.boltlmm.HD.results.problematicse | head -1
2:100089703_G_A chr2:100089703 rs572471451 2 100089703 G A 1 -5.25527e+09 -nan 1.0E+00 1.0E+00 0.555176 0 nan nan NA NA NA
sort -g -k 10 ${MYPATH}/allchr.boltlmm.HD.results.problematicse | tail -1
19:9753556_A_T chr19:9753556 rs553339051 19 9753556 A T 0.999999 5.96188e+09 -nan 1.0E+00 1.0E+00 0.615645 1e-06 nan nan NA NA NA
# Check if any rsIDs or Affymetrix IDs present instead of SNPIDs
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | awk '{if($1 ~ /rs/){print $0}}' > ${MYPATH}/allchr.boltlmm.HD.results.rsids
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | awk '{if($1 ~ /Aff/){print $0}}' > ${MYPATH}/allchr.boltlmm.HD.results.affymetricids
zcat ${MYPATH}/allchr.boltlmm.HD.results.gz | grep -v SNPIDwithAlleles | awk '{if($1 ~ /:/){print $0}}' > ${MYPATH}/allchr.boltlmm.HD.results.snpids
### There is a simple “fudge factor” that you can use to covert linear regression coefficients from GEMMA/EMMAX to log-OR.  
### You just have to multiply beta and SE by 1/(pi*(1-pi)) where pi is the proportion of the sampled individuals that are cases.
### So if 50% of the sample are cases you multiply beta and SE from GEMMA by (1/(0.5*0.5)) = 4.
awk '$NF==0' HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt | wc -l
47792
awk '$NF==1' HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt | wc -l
191634
#HD:Therefore Ncases=47792, Ncontrols=191634, N=239426, pcases=0.1996, pcontrols=0.8004, 1/(pcases*(1-pcases))=6.259395
#SR:Therefore Ncases=35041, Ncontrols=140429, N=175470, pcases=0.1997, pcontrols=0.8003, 1/(pcases*(1-pcases))=6.257042
#HD_hip:Therefore Ncases=11286, Ncontrols=45354, N=56640, pcases=0.1992, pcontrols=0.8008, 1/(pcases*(1-pcases))=6.268831
#HD_knee:Therefore Ncases=19502, Ncontrols=78279, N=97781, pcases=0.1994, pcontrols=0.8006, 1/(pcases*(1-pcases))=6.264108
#HD_hipknee:Therefore Ncases=28943, Ncontrols=116223, N=145166, pcases=0.1994, pcontrols=0.8006, 1/(pcases*(1-pcases))=6.264108
#SRplusHD:Therefore Ncases=68623, Ncontrols=247846, N=316469, pcases=0.2168, pcontrols=0.7832, 1/(pcases*(1-pcases))=5.889359
#OP_hip:Therefore Ncases=9592, Ncontrols=38579, N=48171, pcases=0.1991, pcontrols=0.8009, 1/(pcases*(1-pcases))=6.271197
#OP_knee:Therefore Ncases=9062, Ncontrols=36515, N=45577, pcases=0.1988, pcontrols=0.8012, 1/(pcases*(1-pcases))=6.278309
#OP_hipknee:Therefore Ncases=17778, Ncontrols=71536, N=89314, pcases=0.1990, pcontrols=0.801, 1/(pcases*(1-pcases))=6.273565
echo "StrictlySNPID SNPID rsID CHR POS EA NEA EAF beta SE pval_bolt-lmm pval_inf INFO MAF beta_adjusted se_adjusted OR ORlower ORupper" > ${MYPATH}/allchr.boltlmm.HD.results.fudgefactor
awk '{if(NR>1){$15=$9*6.259395;$16=$10*6.259395;if($9<50 && $9 > -50) {$17=exp($9*6.259395);$18=exp($9*6.259395-1.96*$10*6.259395);$19=exp($9*6.259395+1.96*$10*6.259395)} else {$17="NA";$18="NA";$19="NA"}; print $0;}}' allchr.boltlmm.HD.results >> ${MYPATH}/allchr.boltlmm.HD.results.fudgefactor
#QQ and Manhattan
#Unfiltered
~ag15/scripts/man_qq_annotate --chr-col CHR --pos-col POS --auto-label --pval-col pval_bolt-lmm --title "P_BOLT-LMM-HD.infogt03" --sig-thresh 5e-08 --sig-thresh-line 5e-08 allchr.boltlmm.HD.results.fudgefactor P_BOLT-LMM-HD.HRC.infogt03
#maf 0.01
awk '$14>=0.01' allchr.boltlmm.HD.results.fudgefactor > allchr.boltlmm.HD.results.fudgefactor.mafgte001
wc -l allchr.boltlmm.HD.results.fudgefactor.mafgte001
~ag15/scripts/man_qq_annotate --chr-col CHR --pos-col POS --auto-label --pval-col pval_bolt-lmm --title "P_BOLT-LMM-HD.infogt03.mafgte001" --sig-thresh 5e-08 --sig-thresh-line 5e-08 allchr.boltlmm.HD.results.fudgefactor.maf001 P_BOLT-LMM-HD.HRC.infogt03.mafgte001
#7745334
#maf 0.05
awk '$14>=0.05' allchr.boltlmm.HD.results.fudgefactor > allchr.boltlmm.HD.results.fudgefactor.mafgte005
wc -l allchr.boltlmm.HD.results.fudgefactor.mafgte005
#5434184
~ag15/scripts/man_qq_annotate --chr-col CHR --pos-col POS --auto-label --pval-col pval_bolt-lmm --title "P_BOLT-LMM-HD.infogt03.mafgte005" --sig-thresh 5e-08 --sig-thresh-line 5e-08 allchr.boltlmm.HD.results.fudgefactor.mafgte005 P_BOLT-LMM-HD.HRC.infogt03.mafgte005
#p<10^-7 
cat <(grep -w SNPID allchr.boltlmm.HD.results.fudgefactor.mafgte001) <(awk '$11<=0.000001' allchr.boltlmm.HD.results.fudgefactor.mafgte001) > allchr.boltlmm.HD.results.fudgefactor.mafgte001.pval10-7


#Creating snp stats with qctool 2
#Prepare files for cases, controls and missing samples (NA)
awk '{if($NF==0){print $1}}' HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt | wc -l
awk '{if($NF==0){print $1}}' HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt > 47792_HD_cases.txt
#191634_HD_controls.txt, 247983_HD_missing.txt
#Prepare sample file
awk '{print $1,$2,$3,$NF }' HD_4xcontrols_forBOLT_0cases1controls.sample.noUKdups.txt > Samplefile_for_qctool.txt
sed 's/-9/NA/g' Samplefile_for_qctool.txt > Samplefile_for_qctool.final.txt
#(add second line of the sample file if needed)
#run qctool 2
module add hgi/qctool/latest
for i in {1..22}; do echo ./qctool_v2_controls.sh; done | ~ag15/array 20g HDqctool_controls | sed 's/red/green/'
bsub -G helic -q basement -R"select[mem>20000 && green] rusage[mem=20000]" -M20000 -J "HDqctool_controls[1-22]"  -o HDqctool_controls.%I.o -e HDqctool_controls.%I.e './tmp.445788.exe $LSB_JOBINDEX'
#Checking & tidying up
grep -w Successfully HDqctool_casescontrols.*.o
grep -w Successfully HDqctool_cases.*.o
grep -w Successfully HDqctool_controls.*.o
mkdir qctool_snp-stats
mv chr*.cases-controls.snp-stats qctool_snp-stats/
mv chr*.cases.snp-stats qctool_snp-stats/
mv chr*.controls.snp-stats qctool_snp-stats/
mv HDqctool_cases-controls.*.* qctool_snp-stats/
mv HDqctool_cases.*.* qctool_snp-stats/
mv HDqctool_controls.*.* qctool_snp-stats/
mv Samplefile_for_qctool.* qctool_snp-stats/
mv qctool_v2_*.sh qctool_snp-stats/
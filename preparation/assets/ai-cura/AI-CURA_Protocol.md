AI-CURA: 

Procedure: Supporting publications and all of their supplementary materials were submitted together with the prompts listed below to the LLMs. Reformatting Excel files into comma-separated value (CSV) files and Word files into text files (txt) was performed as these document formats are not recognized by both LLMs

Version of DeepSeek-R1: https://poe.com/DeepSeek-R1 (DeepSeek-R1-0528)
Version of o3-mini-high: https://poe.com/o3-mini-high (OpenAI: o3-mini-2025-01-31)

Prompts used:

PVS1(RNA)
Examine the attached paper. Determine if the paper included any RNA study to investigate the splicing effect of the target variant “NM_xxx:c.xxx p.xxx”.
If the paper included an RNA study including RNA-seq,  RT-PCR and mini-gene assay followed by a Sanger sequencing to study the splice effect of the target variant, determine the splice effect of the variant (intron retention, exon skipping, uses of new splice site rather than natural splice site or other effect mentioned by the author) and the effect on the mRNA or protein transcript (whether the variant disrupts the reading frame and causes the mRNA to undergo nonsense-mediated decay, number and proportion of amino acid residues inserted or deleted from the protein product if reading frame is preserved). After determining the splice effect and the mRNA transcript effect, determine whether the proportion of the alternative transcript produced by the target variant is "complete", "near complete", "incomplete" or "not specified" as mentioned by the authors.
Output the results in a table format for the presence or absence of an RNA study, type of RNA study, tissues where the RNA sample were obtained, splice effect of the variant, and variant effect on the mRNA or protein product.

PS2/PM6
The “gene-disease consistency table.csv” was extracted from Supplementary Table 8. The "ClinGen PS2/PM6 recommendation.pdf" was downloaded from https://clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf .
1st prompt:
Examine the publication attached to find out the number of de novo occurrences of the target variant in only the patients and families recruited. The target variant is “NM_xxx:c.xxx p.xxx”.
Data from the main text and the supplementary materials should all be counted. Always check the supplementary materials after checking the main text. Only consider cases recruited by the authors of the publication, do not include those only in citation. 
First, check the number of de novo cases with the target variant. Only count one de novo case from each family. If there are de novo cases carrying the target variant, then check if parental relationships can be confirmed. Parental relationship can be considered as confirmed only if:
-	Trio testing (proband + both parents) with Genome sequencing (GS), exome sequencing (ES), or short tandem repeat (STR) assay
-	Sanger sequencing alone for parents is considered unconfirmed. 
-	If the authors use phrases like 'family-based Sanger confirmation' or 'Sanger validation in parents' without mentioning trio ES/GS or STR, classify parental relationships as unconfirmed. 
-	If the authors only use phrases like "established de novo" without mentioning the exact tests performed, classify parental relationships as unconfirmed.
For cases with parental germline mosaicism, germline mosaicism must be confirmed for the cases to be counted.
For diseases/phenotypes associated with the patient(s) carrying the de novo target variant, list all clinical diagnoses with the specific names of diseases/syndromes associated with them. Tentative diagnoses that were mentioned as inclusion criteria for all patients in the study can be included. If a clinical diagnosis of a specific disease is not available, list the patient’s phenotypes. If a specific DNA methylation signature is found in the patient with the de novo target variant, also include “Specific DNA methylation signature” in the output phenotype.
Output the number of “de novo patients with confirmed parental relationships”, the number of “de novo patients with unconfirmed parental relationships” and the diseases/phenotypes associated with the de novo patients identified in the paper.
2nd prompt
In the uploaded “gene-disease consistency table.csv”, check if the target gene is listed. If it is included in the list, then check the phenotypic consistency of the gene and the diseases/phenotypes reported in the 1st prompt. Only consider the phenotypic consistency listed in “gene-disease consistency table.csv”. Do not judge by existing knowledge not included in the “gene-disease consistency table.csv”. If the target gene is not listed in the table, assign the phenotypic consistency as "Consistent with gene but not highly specific and high genetic heterogeneity". If the target gene is only listed in the row “Specific DNA methylation signature” but the patients reported do not have the specific methylation signature, assign the phenotypic consistency as "Consistent with gene but not highly specific and high genetic heterogeneity".
Gene: xxx (from 1st prompt)
Disease/Phenotype:
The patient was diagnosed with xxx syndrome (from the 1st prompt). The clinical features include: (from the 1st prompt)
Score according to Table 1 in "ClinGen PS2/PM6 recommendation.pdf" and output the sum of the score.

PS3
Examine all uploaded papers and supplementary materials in which the authors confirmed the target variant NM_xxx:c.xxx p.xxx has functional impact by performing experiments. Data from the main text and the supplementary materials should all be counted. Only consider experiments performed by the authors of the publication, do not include those only in citation.    
Check if normal/negative/wild-type control and abnormal/positive/null controls are used in the paper. Abnormal/positive/null controls should not include the target variant. Only consider variants appearing in the main text or the supplementary materials. Do not assume the number of controls if they are not clearly stated in the paper. Do not make any assumptions.  
PS3_Supporting can be applied if  
- the target variant has deleterious functional impact confirmed by experiment AND  
- normal/negative/wild-type control was used AND  
- abnormal/positive/null control was used.  
OR  
- the experiment used multi-cellular model organism such as fly, zebra fish or mice.  
Check if the publication tested a large number of variants using high-throughput assays (e.g. multiplexed assays of variant effect (MAVE) study).    
Output a table listing:  
- deleterious functional impacts of the target variant 
- if PS3_Supporting can be applied  
- whether the evidence is obtained from in vitro experiments, in vivo experiments, or model organisms  
- the normal/negative/wild-type control used  
- the abnormal/positive/null control used  
- if the publication tested a large number of variants (which informs human curators to check possible upgrade of PS3 according to ClinGen recommendation)

PS4
1st prompt (run separately for every publication):
From the attached paper, output a text file of the following format:
"List of authors:
Description of patients included:
Number of patients with the target variant" 
The list of authors should contain both names and their affiliations. Include all authors. The description of patients included should specify where the patients are from and the period when they were recruited.  
Find out the number of patients carrying the target variant in only the patients and families reported. The target variant is “NM_xxx:c.xxx p.xxx”. Data from the main text and the supplementary materials can all be counted. Always check supplementary materials after checking main text. Only consider cases reported or probands directly recruited by the authors of the publication, do not include those only referenced from other literature. For every family with the target variant, count the proband only and do not count other carriers.
2nd prompt (identify and remove publications with overlapping patients):
From the list attached, identify papers with potential overlapping patients carrying the target variant. Consider the list of authors, places where the patients were recruited and when the patients were recruited. Remove potentially overlapping patients and keep the largest sum of the number of patients in each geographic region without overlapping. The attached list has the following format for each paper:
"PMID:
List of authors:
Description of patients included:
Number of patients with the target variant:"
 Output the list of papers following the format below:
"Region:
List of publications (PMID, recruitment period, number of patients with the variant):"
At the end also output the sum of all included patients carrying the target variant

PM3
1st prompt (run separately for main text and every supplementary file):
A summary table of data from supplementary files is included below:
In trans variant/homozygous target variant	Disease	Zygosity	Number of probands with confirmed phasing 	Number of probands with unconfirmed phasing 	Method of phasing
					
The target variant is “NM_xxx:c.xxxx p.xxx”. 
Your task is to update the above-attached summary table with patients carrying another pathogenic variant detected in trans with the target variant or patients reported as homozygotes of the target variant in the uploaded file. For every combination of variants, count only one patient if multiple patients were reported from the same family. Record whether their phasing is confirmed. Only consider data existing in the current uploaded file. Other files will be processed separately. 
Only genotype established from parents or siblings by genetic tests (e.g. parental testing, haplotype analysis, or segregation analysis) can be considered valid confirmation of phasing. If only one parent is tested and found to carry one allele, variants can be counted as in trans. Do not assume the phasing is confirmed just by recessive inheritance. Report probands with unconfirmed phasing as "Number of probands with unconfirmed phasing”. Make sure the proband has been directly recruited by the authors but not quoted from other literature. 
Output the results in table format with the following columns:
In trans variant/homozygous target variant	Disease	Zygosity	Number of probands with confirmed phasing 	Number of probands with unconfirmed phasing 	Method of phasing
					
If no variants can be found or the exact variant cannot be identified in the uploaded file, output an empty table.
Then update the results in the summary table.
2nd prompt (update phasing detail of patients identified in supplementary material):
An overall summary table of data is included below:
In trans variant/homozygous target variant	Disease	Zygosity	Number of probands with confirmed phasing 	Number of probands with unconfirmed phasing 	Method of phasing
variant 1	Disease A	Compound heterozygous	1	0	Sanger
variant 2	Disease A	Compound heterozygous	0	2	NA
variant 3 	Disease A	Homozygous	1	0	NA

The target variant is “NM_xxx:c.xxx p.xxx”.
In the uploaded main text file, update the phasing status of the findings in the summary table if new information is available in the main text. 
Only genotype established from parents or siblings by genetic tests (e.g. parental testing or haplotype analysis) can be considered valid confirmation of phasing. If only one parent is tested and found to carry one allele, variants can be counted as in trans. Do not assume the phasing is confirmed just by recessive inheritance. Report probands with unconfirmed phasing as "Number of probands with unconfirmed phasing”. Make sure the proband is directly recruited by the authors but not quoted from other literature.
 
PP1
"PP1 recommendation.csv" was extracted from Supplementary Table 9.
Using all uploaded papers and supplementary materials, in only the patients and families reported, find out the number of cosegregations of affected individuals and the number of cosegregations of unaffected individuals according to "PP1 recommendation.csv". Always check supplementary materials after checking the main text. The target variant is "NM_xxx:c.xxx p.xxx "
Affected individuals do not need to have all phenotypes in order to be counted as affected. "Affected" refers to patients with phenotypes. Do not count someone as "affected" only because he/she carries the variant. Assume the disease is fully penetrant unless the authors mentioned the disease has incomplete penetrance.
Only count individuals whose genotype is confirmed by genetic tests in the paper. Genotypes must be confirmed by genetic tests. They cannot be deduced from other non-genetic tests. Genetic tests refer to those using DNA or RNA samples as inputs. Do not assume genetic tests were performed and use them to deduce the genotype if the authors did not state it directly. Do not assume the author of the paper made typos.
For a proband carrying a de novo target variant:
- Do not count his/her unaffected parents as "unaffected co-segregation" 
- The offspring of the proband can be counted normally for co-segregation according to the mode of inheritance of the target variant
Do not count the proband (index patient) for affected cosegregation in each family. Each family only has one proband (index patient). If the authors do not specify the proband of the family, then the number of cosegregation equals number of affected patients with the variant minus one. Score according to " PP1 recommendation.csv" and output the sum of the scores. For obligate carriers, they can be counted as having the variant and then further scored based on whether they have the phenotypes. Genetic tests are not required for obligate carriers.

PP4
"PP4 knowledgebase.csv" was extracted from Supplementary Table 10.
Target variant: NM_xxx:c.xxx p.xxx
In the uploaded files, find all individuals carrying the target variant. Then check the features listed in "PP4 knowledgebase.csv". Directly record the phenotypes confirmed by the authors for the individuals analyzed in the publication. Do not record phenotypes in probands only referenced in the literature. Do not assume/deduce whether the individual does or does not have a phenotype. For phenotypes that are specified as inclusion diagnostic criteria for the cohort, consider all included patients as meeting these requirements. For other experiments  confirming a phenotype in a patient, the author must directly state the results at the individual sample level , otherwise do not include them for scoring PP4. Do not deduce the results of experiments of the target patient based on the author's discussions regarding a group of samples. Experimental results from a group of samples are insufficient. Also record the gene panels sequenced for each individual, and assume all genes were checked if WES or WGS were performed. Do not replace sequencing tests with non-genetic tests. 
Calculate the points awarded for every individual according to the "PP4 evidence level to be assigned" column in "PP4.csv". Do not assume any documents contain typos. If the csv assigns PP4 using a point-based system, then in each unrelated family, choose the one with the highest score. Follow the criteria for scoring points strictly. No point should be awarded if the criteria are only partially met. Use the highest score to decide the PP4 evidence to be applied (Supporting, Moderate or Strong). Otherwise, if the csv assigns PP4 using a criteria-based system, choose the proband that meets the criteria for the strongest PP4 and output the strength assigned to that individual.  

Prompt for summarizing graphical pedigree
Summarize the pedigree with the following information for every individual according to the attached figure legend:
-	Carrier state or genotype of the individual
-	Phenotype of the individual
Output as a list and label every individual with the number of generations and the order (from left to right) of the individual. If any information is missing, still label the individual but do not include any information that is not mentioned in the pedigree.


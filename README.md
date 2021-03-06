# Tübingen Danio Sequencing 
This will be the place to find bioinformatic commands, scripts and tools.


# 1. Fins RNA-Seq
# 2. Light Vs Dark Stripes
# 3. Mutant RNA-Seq
# 4. Danio species Skin RNA-Seq
# 5. Danio Hybrid RNA-Seq and DNA-Seq
### Location of FASTQ files:
> `/ebio/ecnv_projects/hybrid_RNA-Seq/data/danio_hybrid/S1906/S1906_1/01_fastq`

### Location of data description:
> `/ebio/ecnv_projects/hybrid_RNA-Seq/doc`

### Samples file for STAR and DeSEQ2 analysis:
> `/ebio/ecnv_projects/hybrid_RNA-Seq/doc/samples_skin_trunk_pair.txt`

### Create working directory
```bash
mkdir your_dir_name
cd your_dir_name
```
### Create soft links to FASTQ files

```bash
ln -s /path/to/your/fatsq/files/file.1.fastq.gz file.1.fastq.gz
```

####	Example

```bash
for i in {1..28};
do
	ln -s /ebio/ecnv_projects/hybrid_RNA-Seq/data/danio_hybrid/S1906/S1906_1/01_fastq/S1906Nr${i}.1.fastq.gz S1906Nr${i}.1.fastq.gz
	ln -s /ebio/ecnv_projects/hybrid_RNA-Seq/data/danio_hybrid/S1906/S1906_1/01_fastq/S1906Nr${i}.2.fastq.gz S1906Nr${i}.2.fastq.gz
done
```


## STAR alignment and DeSeq2 analysis
### RNAseq: alignment to reference genome and differential expression analysis
Modified from:
https://github.com/najasplus/STAR-deseq2


### Input 
Your working directory should contain:
* Paired end reads named SampleName1.1.fastq.gz SampleName1.2.fastq.gz
* sample_description.txt - tab delimited file with at least two columns: sample and condition (primary condition according to which the differential gene expression analysis will be run. See the exemplary sample_description.txt file and follow DESeq2 guidelines for producing sample description file.

sample	condition<br/>
SampleName1	condition1<br/>
SampleName2	condition1<br/>
SampleName3	condition2<br/>
SampleName4	condition2<br/>

The basic executed command is 

```bash
./star_mapping.sh -s sample_description.txt -d working_directory
```
or
 ```bash
path/to/star_mapping.sh -s sample_description.txt -d working_directory
 ```

If you want to pass additional parameters to STAR aligner, put them in the quotes and add with -o flag:

``` bash
./star_mapping_cluster.sh -s sample_description.txt -o "--outFilterScoreMinOverLread 0.3" 
```

The output of STAR alignment will be stored in the STAR_output folder.
The gene counts for individual samples will be stored in gene_counts folder

The consolidated count matrices will be saved in working directory. Depending on the sequencing library preparation method, you should choose one of them:
* col2_df_raw_counts.tsv - for unstranded library
* col3_df_raw_counts.tsv - stranded, forward
* col4_df_raw_counts.tsv - stranded, reverse

#### Analyze on cluster
You can run STAR mapping on compute cluster. For this

* The script star_mapping_cluster.sh and star_mapping.sh should be located in the same directory
* Modify star_mapping_cluster.sh by passing -s -d and (optionally) -o parameters to star_mapping.sh
* Submit your script to the cluster from your current working directory

``` bash
qsub -cwd star_mapping_cluster.sh
```

### Count normalization and differential expression analysis

Depending on the type of sequencing library you should choose one of the raw_counts matrices produced during the previous step. You pass sample description and count_matrix to the R script as positional arguments (order matters!). The normalized counts as well as pairwise condition differential expression analysis and visualizations will be outputted to the whole_matrix_output subdirectory of your working directory

``` bash
Rscript deseq2_analysis.R sample_description.txt count_matrix.tsv
```

## Matrix file prep 
The matrix file header now needs to be modified for downstream analysis.
```bash
#count columns
head -n 1 FILE | awk '{print NF}'
head -n 1 normalized_counts_deseq.tsv | awk '{print NF}'

#there are 31 col

awk 'NR==1; NR > 1 {print $0 | "sort -n -k 1,1"}' col4_df_raw_counts.tsv > col4_df_raw_counts_sorted.tsv

#remove the normilised couts plus annotation but removes the geneID in col 1
cut -d$'\t' -f 2-31 normalized_counts_deseq.tsv > tmp_cut_normalized_counts_deseq_default.tsv
paste col4_df_raw_counts_sorted.tsv tmp_cut_normalized_counts_deseq_default.tsv > raw_normalized_counts_deseq.tsv
rm -r tmp_cut_normalized_counts_deseq_default.tsv
head -n 1 raw_normalized_counts_deseq.tsv > head_raw_normalized_counts_deseq.tsv
#now edit the header!!!

tail -n +2 raw_normalized_counts_deseq.tsv > tmp_tail_raw_normalized_counts_deseq.tsv
cat head_raw_normalized_counts_deseq.tsv tmp_tail_raw_normalized_counts_deseq.tsv > R_raw_normalized_counts_deseq.tsv
```


## Principal Component Analysis (PCA)
```bash
Rscript pca_rnaseq.R your_matrix.tsv your_sample_description.txt
```

##  Creating countplots from RNA-Seq data 


### Coutplot generation matrix generation
Count plots matrix generatioin requires a list of genes to grep out of the complete count matrix

```bash
for i in $(ls /ebio/ecnv/dooley/bio_script/genes/)
do
	grep -f /ebio/ecnv/dooley/bio_script/genes/${i} /ebio/ecnv_projects/fins_rna_seq/data/STAR-RNAseq/rerun/fastq/whole-matrix-output/R_DS_LS_default_counts_normalized_counts_deseq_default.tsv > tmp_R_${i}
	head -n 1 /ebio/ecnv_projects/fins_rna_seq/data/STAR-RNAseq/rerun/fastq/whole-matrix-output/R_DS_LS_default_counts_normalized_counts_deseq_default.tsv | cat - tmp_R_${i} > R_DS_LS_${i}
	rm -f tmp_R_${i}
done
```
### Coutplot script

```bash
for i in $(ls R_*)
do
Rscript /Users/dooley/bio_script/graph_rnaseq_counts_line_mean_sd.R ${i} /Users/dooley/Documents/Tuebingen/hybrid_RNA-Seq/samples_skin_trunk_pair_2.txt counts_line_${i}.pdf
done
```

##  Extracting sequence and calling variants (mpileup) regions from .bam files

```bash
for i in {1..18} 21 22 {25..28};
do
	samtools index S1906Nr${i}Aligned.sortedByCoord.out.bam
	samtools view -b S1906Nr${i}Aligned.sortedByCoord.out.bam 15:40217256-40267485 > kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.bam
	bcftools mpileup -Ou -f /ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.bam | bcftools call -mv -Oz -o kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.vcf.gz
	bcftools index kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.vcf.gz
	bcftools norm -f /ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.vcf.gz -Ob -o kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.norm.bcf
	bcftools filter --IndelGap 5 kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.norm.bcf -Ob -o kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.norm.flt-indels.bcf
	bcftools index kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.norm.flt-indels.bcf
	cat /ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/Danio_rerio.GRCz11.dna.primary_assembly.fa | bcftools consensus kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.norm.flt-indels.bcf > kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus_Danio_rerio.GRCz11.dna.primary_assembly.fa
	bedtools getfasta -fi kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus_Danio_rerio.GRCz11.dna.primary_assembly.fa -bed /ebio/ecnv_projects/danio_species_skin/code/extract_bed/kcnj13_start_157.bed -fo kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.fa
	cat kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.fa | sed -e '1!{/^>.*/d;}' | sed  ':a;N;$!ba;s/\n//2g' > kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.merged.fa
	python3 /ebio/ecnv/dooley/checkouts/bio_script/fa_rev_com.py kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.merged.fa > kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.merged.revcomp.fa
	python3 /ebio/ecnv/dooley/bio_script/dna2proteins.py -i kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.merged.revcomp.fa -o kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus.merged.revcomp_protien.fa
	rm -f kcnj13_S1906Nr${i}Aligned.sortedByCoord.out.consensus_Danio_rerio.GRCz11.dna.primary_assembly.fa
done
```

## Allele specific expression plots

```R
#attempt to make plots of allele frequency

library(ggplot2)
library(reshape2)


#set working directory  
setwd("/Users/dooley/Documents/Tuebingen/hybrid_RNA-Seq/score/allele_specific_exp/kcnj13")
getwd()

allele_table <- read.table("kcnj13_r_vcf_input_skin_p2.tsv", header = TRUE, sep = "\t")
print(allele_table)

#calculate allele frequency of ref allele
allele_table$S1906Nr15_frq <- (allele_table[,4]/(allele_table[,4] + allele_table[,5]))*100

allele_table$S1906Nr17_frq <- (allele_table[,6]/(allele_table[,6] + allele_table[,7]))*100

allele_table$S1906Nr21_frq <- (allele_table[,8]/(allele_table[,8] + allele_table[,9]))*100

allele_table$S1906Nr25_frq <- (allele_table[,10]/(allele_table[,10] + allele_table[,11]))*100

allele_table$S1906Nr27_frq <- (allele_table[,12]/(allele_table[,12] + allele_table[,13]))*100

#allele_table$S1906Nr11_frq <- (allele_table[,14]/(allele_table[,14] + allele_table[,15]))*100

#allele_table$S1906Nr13_frq <- (allele_table[,16]/(allele_table[,16] + allele_table[,17]))*100

#calculate total counts
allele_table$S1906Nr15_count <- (allele_table[,4] + allele_table[,5])
allele_table$S1906Nr17_count <- (allele_table[,6] + allele_table[,7])
allele_table$S1906Nr21_count <- (allele_table[,8] + allele_table[,9])
allele_table$S1906Nr25_count <- (allele_table[,10] + allele_table[,11])
allele_table$S1906Nr27_count <- (allele_table[,12] + allele_table[,13])
#allele_table$S1906Nr11_count <- (allele_table[,14] + allele_table[,15])
#allele_table$S1906Nr13_count <- (allele_table[,16] + allele_table[,17])

#subet out frequencies
frq <- subset(allele_table, select = c(pos, S1906Nr15_frq, S1906Nr17_frq, S1906Nr21_frq, S1906Nr25_frq, S1906Nr27_frq))
print(frq)

#melt for ggplot2
frq_reshape <- melt(frq, id.vars = "pos", variable.name = "Sample", value.name = "Frq", na.rm=TRUE)
frq_reshape$pos <- as.character(frq_reshape$pos)
print(frq_reshape)

#subset out the total counts
count <- subset(allele_table, select = c(pos, S1906Nr15_count, S1906Nr17_count, S1906Nr21_count, S1906Nr25_count, S1906Nr27_count))
print(count)
#melt for ggplot2
count_reshape_2 <- melt(count, id.vars = "pos", variable.name = "Sample", value.name = "Count", na.rm=TRUE)
print(count_reshape_2)
count_reshape_2$pos <- as.character(count_reshape_2$pos)

#produce plot of counts at SNPs
p1 <- ggplot(count_reshape_2, aes(x=pos, y=Count)) + geom_violin(fill="#C0C0C0", adjust=1.0, trim=TRUE)
p1 <- p1 + theme_bw()
p1 <- p1 + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")
p1 <- p1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="blue")
p1 <- p1 + labs(title = ""~italic(kcnj13)~" Pair2 Skin counts at SNPs") + theme_bw() + theme(plot.title = element_text(hjust=0.5))
p1 <- p1 + xlab("SNP Position on Chr")
p1

#produce plot of ref frequency
p2 <- ggplot(frq_reshape, aes(x=pos, y=Frq)) + geom_violin(fill="#C0C0C0", adjust=1.0, trim=TRUE)
p2 <- p2 + theme_bw()
p2 <- p2 + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")
p2 <- p2 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="blue")
p2 <- p2 + labs(title = ""~italic(kcnj13)~" Pair 2 Skin "~italic(D.rerio)~" Allele Frequency") + theme_bw() + theme(plot.title = element_text(hjust=0.5))
p2 <- p2 + ylim(0, 100)
p2 <- p2 + xlab("SNP Position on Chr")
p2 <- p2 + ylab("D.rerio Allele Frequency")
p2
````
	

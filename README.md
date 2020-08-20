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

####Example

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

## Principal Component Analysis (PCA)
```bash
Rscript pca_rnaseq.R your_matrix.tsv your_sample_description.txt
```

##  Creating countplots from RNA-Seq data 


### Coutplot generation

##  Extracting sequence regions from .bam files



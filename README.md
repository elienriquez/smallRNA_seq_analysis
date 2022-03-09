# Small RNA-seq analysis
Here it is detailed the main steps for the analysis of small RNA-Seq data from Trichoderma atroviride using data obtained from the NCBI website under Gene Expression Omnibus (GEO) accession number: GSE190033.

# Objetive

The goal of this repository is to distribute the pipeline used for the analysis of small RNA-seq data from T. atroviride during mycoparasitism.

# Biological significance

To analyze the response of the fungus T. atroviride (wild type and the mutant strain Ddcr2) during three different stages of mycoparasitism against the prey fungus Alternaria
alternata.

# Sequencing technology used and obtained libraries

24 libraries were sequenced with Illumina TruSeq 1x36 single-end.

# Quality analysis, adapter removal and read mapping

The quality of the small RNA-Seq libraries was analyzed using FastQC version 0.11.8.

fastqc &lt;input_file_name.fastq&gt; --outdir /path-to-outputdirectory

Search and removal of adapters was performed using Kraken (Davis et al., 2013). The search for adapters was with the Kraken minion function, using the following command:

minion search-adapter -i &lt;archive-.fastq&gt;

Adapter removal was performed with the following command:

reaper -i &lt;input_file_name. fastq&gt; -geom no-bc -3pa TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATC -tabu GTTCAGAGTTCTACAGTCCGACGATC -mr-tabu 14/2/1 -3p-global 6/1/0/0 -3p-prefix 11/2/1 -3p-head-to-tail 1 -dust-suffix 20/ACTG -format-clean %I%n%C%n%f%n -nnn-check 1/1 -qqq-check 35/10 -basename &lt;output _file_name&gt; --fasta-out

# Mapping
Clean reads were aligned to the newest version of T. atroviride genome IMI206040 (Atriztán-Hernández et al., in preparation) using Bowtie version 1.2.3, using the following command:

bowtie -f -v 3 --best -k 100 --strata --sam --chunkmbs 512 &lt;index_genome&gt;&lt;input_file_name.lane.clean.gz&gt; &lt; output _file_name.sam&gt;

# Annotation, quantification and miRNA prediction

ShortStack version 3.8.5 was used to perform miRNA prediction, annotation of smallRNA loci and quantification of reads. The parameters used were the following:

ShortStack --readfile path_to_fasta_file --dicermin 20 --dicermax 26 --mistmaches 1 --foldsize 600 --show_secundaries --genomefile path_to_genome_file --outdir path_to_output_file

# Differential expression (DE) analysis of small RNA libraries

To perform DE analysis to smallRNA-seq libraries we used R packages: DESeq2 and Rsubread
In order to obtain a matrix of reads the function featurecounts was used to each library. The annotation file and bam files were obtained using ShortStack as pointed above.

library(Rsubread)
library(DESeq2)

counts<-featureCounts("path_to_bam_file",
                         annot.ext = "path_to_annotation_file",
                         isGTFAnnotationFile = TRUE,
                         GTF.featureType = "nc_RNA",
                         GTF.attrType = "ID",
                         countMultiMappingReads = TRUE,
                         fraction = TRUE,
                         allowMultiOverlap = TRUE)
                     
write.csv(counts$counts,file="path_to_destination")

To perfom the DE analysis we used the following script in R
Load reads matrix:
sRNA<-read.csv("path_to_csv_file", header = TRUE,row.names = 1)

Creating the design matrix, in this case for the comparison: WT control vs WT After Contact AC:

colnames(sRNA)<-c("WT_control_R1","WT_control_R2","WT_control_R3",
                  "WT_AC_R1","WT_AC_R2","WT_AC_R3")
colData<- cbind(c(rep("control",3),
                  rep("AC",3)))
row.names(colData)<-colnames(sRNA)
colnames(colData)<-c("treatment")

Because DESeq2 accept rounded numbers we used the function "round":

sRNA    <-round(sRNA)

Creating a DESeq object:

dds <- DESeqDataSetFromMatrix(countData = sRNA,
                              colData = colData,
                              design= ~treatment)
Creating factors levels:
dds$treatment<-factor(dds$treatment,levels=c("control","AC"))

Filtering: 

keep<-rowSums(counts(dds)) >= 10
dds <-dds[keep,]

DE analysis:

dds<-DESeq(dds,fitType="mean")
resultsNames(dds)

Getting results: 

WT_AC<-results(dds,name="treatment_AC_vs_control",alpha = 0.1)
summary(WT_AC)

Recovering differentially expressed clusters, in this particular example we found 784 DE clusters:

WT_AC_DE<-as.data.frame(WT_AC[order(WT_AC$pvalue,decreasing = FALSE),][1:784,])

Saving results for further analysis:

write.csv(WT_AC_DE,file="path_to_destination")

# Prediction of targets 

The search for targets of small RNAs was carried out using TargetFinder and using the TAPIR online platform with the established parameters. The command used in TargetFinder was the following:

targetfinder_threads.pl -f &lt;input_file_name.fasta&gt; -d &lt;Target_Sequence_Database_File.fasta&gt; -c 4 -t 8 -p table -o &lt;output _file_name.txt&gt;


# Pipeline for getting SMPs and SNPs for Arabidopsis thaliana

1. Obtain A.thaliana vcf from https://1001genomes.org/data/GMI-MPI/releases/v3.1/ and filter

```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps
gunzip 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip 1001genomes_snp-short-indel_only_ACGTN.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -v snps --max-alleles 2  -O z -o 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz   1001genomes_snp-short-indel_only_ACGTN.vcf.gz
```
2. Create pseudogenomes
```
#invcf is a file with sample names that are found in the WGBS and genomic vcf file

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr1.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr2.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 3 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr3.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 4 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr4.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 5 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr5.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Mt | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_Mt.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Pt | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_Pt.fa; done

mkdir /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas

/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr1.fa chr1/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr2.fa chr2/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr3.fa chr3/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr4.fa chr4/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr5.fa chr5/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_Pt.fa Pt/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_Mt.fa Mt/

mkdir /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas/merged
cat ../../invcf | while read line ; do cat $line*[0-9].fa >$line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do java -jar /proj/popgen/a.ramesh/software/picard.jar CreateSequenceDictionary R=$file O=${file/.fa/.dict}; done

```

3. Download methyaltion data from NCBI
```
cd  /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/
/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file PRJNA187927_adabidopsis_single_Acc_List.txt  -O data
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done

/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file PRJNA187927_adabidopsis_paired_Acc_List.txt  -O data/paired
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/paired
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done
```
4. Trim data and concatenate or rename files. Renaming in A.thaliana rename scripts folder
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/
for file in *.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 5 $file ${file/.fastq.gz/_trim.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/paired
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

5. Map methylation reads
```
sed 's/^/\/proj\/popgen\/a.ramesh\/projects\/methylomes\/arabidopsis\/pseudogenomes\//' samplenames5 | paste samplenames5 >samplenames6
sed 's/^/\/proj\/popgen\/a.ramesh\/projects\/methylomes\/arabidopsis\/pseudogenomes\//' samplenames2 | paste samplenames2 >samplenames3

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
cat samplenames6 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 $value1.trim.fq.gz  ; done

cat samplenames3 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
```

6. Deduplicate methylation reads
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
for file in *_bismark_hisat2*.bam ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

7. Call methylation variants
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
cat samplenames3 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam  ; done

cat samplenames6 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.trim_bismark_hisat2.deduplicated.bam  ; done

gunzip *.bismark.cov.gz
gunzip *.CpG_report.txt.gz
```

8. Do binomial test
```
setwd("/proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data")
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")

cov <- read.table(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "_")

context <- read.table(file=contextfiles[1,],header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "_")

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub(gsub("\\..*","",covfile$V1)[1],"",cov_context$chromosome)
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov","",covfile[1,])

b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]

colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  print(i)
  
  cov <- read.table(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "_")
  
  context <- read.table(file=contextfiles[i,],header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "_")
  
  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub(gsub("\\..*","",covfile$V1)[i],"",cov_context$chromosome)
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
  if (is.na(noncoversionrate)) {
    noncoversionrate <- 0
  }
  samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov","",covfile[i,])
  
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  

  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
}
cov_context2$chromosome <- gsub("_.*","",cov_context2$ID)
cov_context2$position   <- as.numeric(gsub(".*_","",cov_context2$ID))

write.table(cov_context2,file="cov_context3.txt",row.names = F)

```

9. Split reference genome into genes
```
cd  /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes
sed -e 's/\t/:/' -e  's/\t/-/' gene_pos.bed >gene_pos.list

cd /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/snps
cat /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $line >>genes.fasta; done
mkdir genes_fasta/
/data/proj2/popgen/a.ramesh/software/faSplit byname genes.fasta genes_fasta/
cd genes_fasta/
ls *fa >filenames
```

10. Rscript to count number of cytosines

```
library("methimpute",lib.loc="/data/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")

## Only CG context
files <- read.table(file="filenames")
files <- as.character(files$V1)

cytsosine_count <- ""
for (f in 1:length(files)){
  print(f)
  cytosines <- c()
  try(cytosines <- extractCytosinesFromFASTA(files[f], contexts = 'CG'),silent=T)
  if (length(cytosines) > 0){
    cytsosine_count <- rbind(cytsosine_count,(c(files[f],table(cytosines$context))))
  } else {
    cytsosine_count <- rbind(cytsosine_count,(c(files[f],0)))
  }
  #print(cytsosine_count)
}
cytsosine_count <- cytsosine_count[-c(1),]
write.table(cytsosine_count,file="cytsosine_count.txt",row.names=F, col.names=F,quote=F,sep="\t")
```

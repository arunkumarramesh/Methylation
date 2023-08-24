# Pipeline for getting SMPs and SNPs

Example for Arabidopsis lyrata shown

1. Download data from NCBI
```
## Genomes usually donwloaded from Ensembl, need to ensure it includes the choloplast contig. Index with samtools faidx.

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file  acc_list_bs  -O data
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file  acc_list_rna  -O data_rna

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data
for file in *.sra; do /data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna
for file in *.sra; do /data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done

```

2. Trim data
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data
for file in *_1.fastq.gz; do java -jar /data/proj2/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna
for file in *_1.fastq.gz; do java -jar /data/proj2/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

3. Index reference genome and map RNA-seq data with with STAR 2-pass mode

```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes
/data/proj2/popgen/a.ramesh/software/STAR-2.7.10b/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir STARindex/ --genomeFastaFiles Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa --sjdbGTFfile Arabidopsis_lyrata.v.1.0.55.gtf --sjdbOverhang 201 --runThreadN 20

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna
ulimit -n 10000
for file in *_1.paired.fq.gz; do /data/proj2/popgen/a.ramesh/software/STAR-2.7.10b/source/STAR --runMode alignReads --genomeDir /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/STARindex/ --readFilesIn $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} --readFilesCommand zcat --outFileNamePrefix ${file/_1.paired.fq.gz//}  --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 20  ; done

find . -type f -name "Aligned.sortedByCoord.out.bam" -printf "/%P\n" | while read FILE ; do DIR=$(dirname "$FILE" ); mv ."$FILE" ."$DIR""$DIR".sort.bam;done
find . -name '*.sort.bam' -exec mv {} . \;

```

4. Mark duplicate reads and index
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
for file in *.sort.bam ; do picard MarkDuplicates -I $file -O ${file/.sort.bam/_marked.bam} -M ${file/.sort.bam/_metrics.txt} 2>mark_err; done
for file in *_marked.bam ; do picard BuildBamIndex -I $file 2>index_err; done
```

5. SplitNCigarReads
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
for file in *_marked.bam ; do gatk SplitNCigarReads -R /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa -I $file -O  ${file/_marked.bam/_split.bam}  ; done
```

6. Add read groups and index
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
for file in *_split.bam ; do picard AddOrReplaceReadGroups -I $file -O ${file/_split.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *_readgroup.bam ; do picard BuildBamIndex -I $file; done
ls *_readgroup.bam | sed 's/*_readgroup.bam//' | paste - lc_files >lc_files2 ## for haplotyper caller script (file also called hc_files2)
```

7. Run haplotype caller using perl script to parallelize
```
#!/usr/bin/perl -w


my $dir = "/data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/";
my $input = "hc_files2";
chdir("$dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input") || die "couldn't open the input file $input!";
                    while (my $line = <INPUT_FILE>) {
                        chomp $line;
my @array = split(/\t/,$line);
my $sample = $array[0];
my $file = $array[1];

my $SLURM_header = <<"SLURM";
#\$-cwd
# This option tells gridware to change to the current directory before executing the job
# (default is the home of the user).

#\$-pe serial 2
# Specify this option only for multithreaded jobs that use more than one cpu core.
# IMPORTANT! Don't use more than 4 cores to keep free space for other students!

#\$-l vf=2g
# This option declares the amount of memory the job will use. Specify this value as accurate as possible.
# IMPORTANT! memory request is multiplied by the number of requested CPU cores specified with the -pe.
# Thus, you should divide the overall memory consumption of your job by the number of parallel threads.

#\$-N lc


source ~/.bashrc

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/

SLURM

  my $tmp_file = "lc.$sample";
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "gatk HaplotypeCaller -R /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa  -I $file -O  $sample".".g.lc.vcf.gz -ERC GVCF\n";
  close SLURM;
  system("qsub $tmp_file");

}

            close(INPUT_FILE);
```


8. Combine GVCFs
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
gatk CombineGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa --variant HighSite_GER_2274.g.lc.vcf.gz  --variant HighSite_J1_182.g.lc.vcf.gz   --variant LowSite_GER_2324.g.lc.vcf.gz  --variant LowSite_J3_1601.g.lc.vcf.gz --variant HighSite_GER_2285.g.lc.vcf.gz  --variant HighSite_J1_2023.g.lc.vcf.gz  --variant LowSite_J1_163.g.lc.vcf.gz    --variant LowSite_J3_1651.g.lc.vcf.gz --variant HighSite_GER_2456.g.lc.vcf.gz  --variant LowSite_GER_2187.g.lc.vcf.gz  --variant LowSite_J1_383.g.lc.vcf.gz    --variant LowSite_J3_238.g.lc.vcf.gz --variant HighSite_J1_1490.g.lc.vcf.gz   --variant LowSite_GER_2214.g.lc.vcf.gz  --variant LowSite_J1_761.g.lc.vcf.gz    --variant LowSite_J3_750.g.lc.vcf.gz --variant HighSite_J1_1658.g.lc.vcf.gz   --variant LowSite_GER_2304.g.lc.vcf.gz  --variant LowSite_J1_985.g.lc.vcf.gz -O lyrata.cohort.g.vcf.gz
```

9. Genotype combined GVCFs

```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa -V lyrata.cohort.g.vcf.gz -O lyrata.output.vcf.gz
```

10. Filter vcf, only keep high quality SNPs
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
gatk SelectVariants -V lyrata.output.vcf.gz -select-type SNP -O lyrata.snps.vcf.gz
vcftools --gzvcf lyrata.snps.vcf.gz --out lyrata_snps_filtered --recode --recode-INFO-all --minDP 20 --minGQ 30  --minQ 30
```

11. Get chromosome and sample names, can be done many ways
cut -f 1 ../genomes/PL_genomeandchloroplast_assemblies_S119.fa.fai >chrlist
ls *_marked.bam | sed 's/_marked.bam//' >samplenames


12. Create combined list containing every pair of chromosome and sample names. In R. 
```
file1 <- read.table(file="samplenames")
file2 <- read.table(file="chrlist")

fileall <- expand.grid(file1$V1, file2$V1)  
write.table(fileall, file = "sample_chr", quote = F, sep = "\t", row.names = F, col.names = F)
```

13. Create pseudogenomes by substituting SNPs for each sample from VCF into the reference genome
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
bgzip lyrata_snps_filtered.recode.vcf
tabix lyrata_snps_filtered.recode.vcf.gz
cat sample_chr | while read -r value1 value2 remainder ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/Arabidopsis_lyrata.v.1.0.dna.toplevel_chrloroplast.fa $value2 | bcftools consensus -M N -s $value1 -p ${value1/$/_} -H 1pIu lyrata_snps_filtered.recode.vcf.gz >>${value2/$/_lyrata.fa}; done
cat chrlist | while read line; do /data/proj2/popgen/a.ramesh/software/faSplit byname $line ${line/^/chr_} ; done
mkdir /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/pseudogenomes
mv *fa /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/pseudogenomes
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/pseudogenomes
cat samplenames | while read line ; do cat $line*.fa >$line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do picard CreateSequenceDictionary -R $file -O ${file/.fa/.dict}; done
cat samplenames | while read line ; do mkdir $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa $line  ; done
cat samplenames | while read line ; do mv $line.merged.dict $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa.fai $line  ; done
sed 's/^/\/data\/proj2\/popgen\/a.ramesh\/projects\/methylomes\/lyrata\/pseudogenomes\//' samplenames | paste samplenames - ## not exactly the same, variants exist
```

14. Now map methylation reads with Bismark
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data/
cat samplenames2 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
```

15. Deduplicate methylation reads with Bismark
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data
for file in *_bismark_hisat2_pe.bam; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

16. Call methylation variants
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data
cat samplenames2 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam ; done
```

17. This script is to create methylation vcfs. I need add this after annotating it properly. Lots done here. In R.
```
...
```

18. Get vcf header
```
zcat lyrata_snps_filtered.recode.vcf.gz | grep '##' > vcfheader
cat vcfheader lyrata_meth.vcf >lyrata_meth_all.vcf
```

19. Get summary statistics for methylation and genomic variants. Done on biallelic SNPs. Further NA filtering in R. Only keep gene body variants.
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/
awk '$3 ~ /^gene$/' Arabidopsis_lyrata.v.1.0.55.chr.gff3 | cut -f 1,4,5 >gene_pos.bed
 
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --site-pi --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --het --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data/
vcftools --vcf lyrata_meth_all.vcf  --out lyrata_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --vcf lyrata_meth_all.vcf  --out lyrata_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --site-pi --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed

```

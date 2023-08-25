# Pipeline for getting SMPs and SNPs for maize

1. Get data
```
cd  /proj/popgen/a.ramesh/projects/methylomes/maize/

/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch -X 9999999999999  --option-file PRJNA607675_maize_Acc_List.txt  -O data_wgs
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3 $file; done

## accession file had both WGS and WGBS. WGBS moved to data/ folder.

cd /proj/popgen/a.ramesh/projects/methylomes/maize/genomes
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_902167145.1/download?include_annotation_type=GENOME_GFF,GENOME_GTF&filename=GCF_902167145.1.zip" -H "Accept: application/zip"


grep 'exon' genomic.gtf | cut -f 1,4,5 >gene_pos.bed
```
2. Split reference gneome into genes
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/genomes
sed -e 's/\t/:/' -e  's/\t/-/' gene_pos.bed >gene_pos.list

cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs
cat /proj/popgen/a.ramesh/projects/methylomes/maize/genomes/gene_pos.list | while read -r line ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/maize/genomes/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna $line >>genes.fasta; done
/proj/popgen/a.ramesh/software/faSplit byname genes.fasta genes_fasta/
```

3. Trim data
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 10 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done

cd /proj/popgen/a.ramesh/projects/methylomes/maize/data
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 10 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

4. Index genome
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/genomes
bwa index  Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
picard CreateSequenceDictionary -R Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -O Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.dict
```

5. Map genomic data
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs
for file in *_1.paired.fq.gz; do bwa mem -t 30 /proj/popgen/a.ramesh/projects/methylomes/maize/genomes/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} 2>paired_map_err | samtools view -Sb - > ${file/_1.paired.fq.gz/_mapped.bam} 2>paired_map_err; done
```

6. Add readgroups, filter reads, sort reads, markduplicates and index bams
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs
for file in *_mapped.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar AddOrReplaceReadGroups -I $file -O ${file/_mapped.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *_readgroup.bam; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b $file >${file/_readgroup.bam/_mq20.bam}; done
for file in *_mq20.bam; do samtools sort -@ 5 $file >${file/_mq20.bam/_sort.bam} ; done
for file in *_sort.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I $file -O ${file/_sort.bam/_marked.bam} -M ${file/_sort.bam/_metrics.txt} ; done
for file in *_marked.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar BuildBamIndex -I $file; done
ls *_marked.bam > lc_files
ls *_marked.bam | sed 's/*_marked.bam//' | paste - lc_files >lc_files2 ## for haplotyper caller script (file also called hc_files2)
```

7. Call variants in each sample
```
#!/usr/bin/perl -w


my $dir = "/proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs";
my $input = "lc_files2ah";
chdir("$dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input") || die "couldn't open the input file $input!";
                    while (my $line = <INPUT_FILE>) {
                        chomp $line;


my @array = split(/\t/,$line);
my $sample = $array[0];
my $file = $array[1];

my $SLURM_header = <<"SLURM";
#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=8000
#SBATCH -J lc
#SBATCH -o lc.out
#SBATCH -e lc.err

cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs

SLURM

  my $tmp_file = "hc.$sample";
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk HaplotypeCaller -R /proj/popgen/a.ramesh/projects/methylomes/maize/genomes/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna  -I $file -O  $sample".".g.vcf.gz  -ERC GVCF\n";
  
  close SLURM;
  system("sbatch $tmp_file");

}

            close(INPUT_FILE);


```

8. Combine gvcfs
```
cd /proj/popgen/a.ramesh/projects/methylomes/maize/data_wgs
```


cut -f 1 GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.fai > intervals.list
/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk --java-options "-Xmx18g -Xms1g" GenomicsDBImport -V JRIAL2A.g.vcf.gz -V JRIAL2B.g.vcf.gz  --genomicsdb-workspace-path genomicsdb --tmp-dir /proj/popgen/a.ramesh/projects/methylomes/maize/tmp -L intervals.list
cat samples | while read line; do /proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk --java-options "-Xmx18g" GenomicsDBImport --genomicsdb-update-workspace-path genomicsdb --tmp-dir /proj/popgen/a.ramesh/projects/methylomes/maize/tmp -V $line; done
```


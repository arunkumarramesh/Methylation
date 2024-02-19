## Pipelines to identify polymorphic genetic and methylation variants

Pipeline to process and analyse genome and bisulfite data from Arabidopsis, rice, maize, and soybean. Goals are to identify polymorphic nucleotide and methylation sites for each species and perform evolutionary analyses. 

## List of files
Pipeline_lyrata.md: Pipeline to call genetic and methylation variants in Arabidopsis lyrata

Pipeline_maize.md: Pipeline to call genetic and methylation variants in maize

Pipeline_rice.md: Pipeline to call genetic and methylation variants in rice

Pipeline_soybean.md: Pipeline to call genetic and methylation variants in soybean

Pipeline_thaliana.md: Pipeline to call methylation variants in Arabidopsis thaliana

Calculation_rice.md: Pipeline to for population structure analyses in rice

Calculation_rice.md: Pipeline to for population genetic analyses in rice

Perl_scripts_for_GATK: Perl scripts for calling variants from aligned BAM files in parallel
- Needed as GATK haplotype caller does not have a multithreading option

rna_seq_variant_calling.nf: Nextflow pipeline to call variants from RNA seq data. To use, run
- nextflow rna_seq_variant_calling.nf --reads SAMPLE*_{1,2}.fastq.gz --genome GENOME.fa --gtf GENOME.gtf --dictname GENOME_FILENAME --cpu 8
- reads: paired read files with full path
- genome: genome (not transcriptome!) reference fasta file  with full path
- gtf: gtf file containing genome annotations for reference fasta file  with full path
- dictname: required for picard create seq dictionary. Filename for genome but without the extension (i.e. just GENOME).
- cpu: number of threads



## Overview of pipeline

![Slide1](https://github.com/arunkumarramesh/Methylation/assets/23363383/15cd99fc-5089-497a-b1d5-f63c8d212c20)

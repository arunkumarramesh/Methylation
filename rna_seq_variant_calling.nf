//rna_seq_variant_calling.nf
nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */
params.reads = "/proj/popgen/a.ramesh/projects/methylomes/rice/test_data/SRR*_{1,2}.fastq.gz"
params.genome = "/proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"
params.gtf = "/proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.55.gtf"
params.dictname = "Oryza_sativa.IRGSP-1.0.dna.toplevel"
params.cpus = 8

log.info """\
         RNA SEQ VARIANT CALLING  PIPELINE
         ===================================
         genome: ${params.genome}
         gtf: ${params.gtf}
         reads        : ${params.reads}
         cpus: ${params.cpus}
         """
         .stripIndent()


process GTCALL {

    input:
    tuple val(pair_id), path(reads)
    path genome
    path gtf
    val dictname
    val cpus

    output:
    path(pair_id)

    script:
    """
    samtools faidx $genome
    java -jar /proj/popgen/a.ramesh/software/picard.jar CreateSequenceDictionary -R $genome -O ${dictname}.dict
    /proj/popgen/a.ramesh/software/STAR-2.7.10b/source/STAR --runMode genomeGenerate --genomeDir STARindex/ --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbOverhang 149 --runThreadN $cpus
    java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads $cpus ${reads[0]} ${reads[1]} ${pair_id}_1.paired.fq.gz ${pair_id}_1.unpaired.fq.gz ${pair_id}_2.paired.fq.gz ${pair_id}_2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    ulimit -n 10000
    /proj/popgen/a.ramesh/software/STAR-2.7.10b/source/STAR --runMode alignReads --genomeDir STARindex/ --readFilesIn ${pair_id}_1.paired.fq.gz ${pair_id}_2.paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $pair_id --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN $cpus
    java -jar /proj/popgen/a.ramesh/software/picard.jar AddOrReplaceReadGroups -I ${pair_id}Aligned.sortedByCoord.out.bam -O ${pair_id}_readgroup.bam -LB species -PL illumina -PU 1 -SM $pair_id
    java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I ${pair_id}_readgroup.bam -O ${pair_id}_marked.bam -M ${pair_id}_metrics.txt
    java -jar /proj/popgen/a.ramesh/software/picard.jar BuildBamIndex -I ${pair_id}_marked.bam
    /proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk SplitNCigarReads -R $genome  -I ${pair_id}_marked.bam -O  ${pair_id}_split.bam
    java -jar /proj/popgen/a.ramesh/software/picard.jar BuildBamIndex -I ${pair_id}_split.bam
    /proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk HaplotypeCaller -R $genome  -I ${pair_id}_split.bam -O  ${pair_id}.vcf.gz  -ERC GVCF
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs( params.reads , checkIfExists:true)
    read_pairs_ch.view()
    genome_ch = Channel.fromPath( params.genome , checkIfExists:true)
    gtf_ch = Channel.fromPath( params.gtf , checkIfExists:true)
    dictname_ch = Channel.value(params.dictname)
    cpus_ch = Channel.value(params.cpus)
    gtcall_ch=GTCALL(read_pairs_ch,genome_ch,gtf_ch,dictname_ch,cpus_ch)
}

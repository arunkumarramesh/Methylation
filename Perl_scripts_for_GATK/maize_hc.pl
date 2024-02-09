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
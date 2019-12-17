#!/usr/bin/perl

# This pipeline was built to run on a HPC cluster using the Torque scheduler

use feature qw(say);
use File::Basename;
use Data::Dumper;
use warnings;
use Cwd;
# submit and write $script in $job_path and return $job_id from qsub
sub submit_job
{
    my ($script, $job_path) = @_;
    open my $job, ">", $job_path or die "Cannot write $job_path";
    print $job $script;
    close $job;
    # save the job id to be used as a dependency for alignment job
    chomp(my $job_id= `qsub $job_path`);
    say "Submited $job_path";
    sleep 1;
    return $job_id;
}

# save root directory (my project directory)
my $pwd = getcwd;

$scripts = "$pwd/scripts";      # location of scripts
$varcalls = "$pwd/BluePrint/varcalls/processed"; # location of variation calls

# load WGS metadata. This is a table that associates the individual id to a WGS bam filename
my %id_to_cram = ();
open my $wgs_metadata, "<", "genome_ids_to_filename";
while(<$wgs_metadata>)
{
    ($id, $cram) = split /\s/;
    chomp($cram);
    @candidate_crams = glob "/sf1/project/bws-221-af-004/epiMAP/EGAD00001002663/*$cram";
    $id_to_cram{$id} = $candidate_crams[0];
}

# Command line arguments
# output directory
# read length to crop
# list of sample names to process
# directory containing fastq files
# do broad peak calling?
my ($genome_dir, $crop_size, $scheduled_genomes, $chipseq_path, $broad) = @ARGV;

if(not defined $broad)
{
    say "Doing narrow peak call";
    $broad = "";
} 
elsif ($broad eq "--broad")
{
    say "Doing broad peak call";
}
# get a list of all chipseq bams we have and store in @chipseq_bams
##my $chipseq_path = "$pwd/BluePrint/H3K4me1"; # location of chipseq data
# read available chipseq bams
opendir my $chipseq_dir, $chipseq_path or die "Cannot open BluePrint bam directory";
my @chipseq_bams = readdir $chipseq_dir;
close $chipseq_dir;
mkdir $genome_dir;
open my $genomes, "<", $scheduled_genomes or die "Cannot get scheduled $genome_dir in $scheduled_genomes";
# read available $genome_dir
# for every individual, submit the required jobs
while(my $g = <$genomes>)
{
    next if($g =~ /^#/);
    chomp $g;
    say "\nProcessing $g";
    mkdir "$genome_dir/$g/";
    # create indiviudalized haploids using vcf2diploid
    my $impute_fasta_script = <<"IMPUTE_FASTA_SCRIPT";
    #PBS -A bws-221-ae
    #PBS -d $pwd/$genome_dir/$g
    #PBS -l nodes=1:ppn=1 -l mem=15gb
    #PBS -l walltime=5:00:00
    #PBS -r n
    module load mugqic/java

    VCF2DIPLOID_JAR=/home/cgroza/vcf2diploid/vcf2diploid.jar
    REF_FASTA=/gs/project/mugqic/projects/cgroza/Homo_sapiens.hg19/Homo_sapiens.hg19.fa

    # make a directory for the genome and move to it
    #vcf2diploid will populate the current dir with .chain .map .fa files
    java -Xmx14G -jar \$VCF2DIPLOID_JAR -id $g -chr \$REF_FASTA -vcf $varcalls/*.vcf.gz
    mkdir maps chains maternal paternal
    # merge chromosomes into one fasta (needed to build bwa index)
    cat *_maternal.fa > maternal/maternal.fa
    cat *_paternal.fa > paternal/paternal.fa
    # organize files since they are needed later for .chain file generation
    rm *_maternal.fa *_paternal.fa
    mv *.chain chains
    mv *.map maps
    IMPUTE_FASTA_SCRIPT
    my $impute_fasta_job_path = "$genome_dir/$g/impute_fasta.sh";
    my $impute_fasta_job_id = submit_job $impute_fasta_script, $impute_fasta_job_path;

    # look for chipseq bam file corresponding to this individual;
    my @individual_bams = grep /$g/, @chipseq_bams;

    # if no bam file found, skip
    if(scalar @individual_bams < 1)
    {
		say "Nothing found for $g";
		next;
    }

    # pick the first bam file in list to use
    my $bam_file =  "$pwd/$chipseq_path" . $individual_bams[0];
    # keep sample name for later
    my $sample_name =  $individual_bams[0];
    say "Found $bam_file for $g";

    my ($sample_id) = basename ((split /[_.]/ , $bam_file)[0]);
    say "Sample ID: $sample_id";
    # script to convert bam files to fastq for alignement on individualized $genome_dir
    my $bam2fastq_script = <<"BAM2FASTQ_SCRIPT";
    #PBS -A bws-221-ae
    #PBS -l nodes=1:ppn=1 -l mem=30gb
    #PBS -l walltime=10:00:00
    #PBS -d $pwd/$genome_dir/$g
    #PBS -r n

    module load mugqic/java mugqic/picard
        mkdir temp fastq
        # check if the input is actually a fastq file and not a fastq
        if [[ "$bam_file" =~ "fastq.gz" ]]; then
                                                echo INPUT IS FASTQ: $bam_file
                                                gunzip -c $bam_file > fastq/$sample_name.untrimmed.fastq

                                                else
                                                java -Djava.io.tmpdir=temp  -Xmx24G -jar \$PICARD_JAR SamToFastq\\
                                                VALIDATION_STRINGENCY=LENIENT \\
                                                INPUT=$bam_file \\
                                                FASTQ=fastq/$sample_name.untrimmed.fastq
                                                fi

                                                INPUT_READS=\$($scripts/assign_input_type.sh $sample_id)

                                                $scripts/crop_reads.sh fastq/$sample_name.untrimmed.fastq fastq/$sample_name.1.fastq $crop_size
                                                rm -f fastq/$sample_name.untrimmed.fastq

                                                $scripts/crop_reads.sh \$INPUT_READS  fastq/input.fastq.gz $crop_size
                                                rm -f fastq/input.untrimmed.fastq
                                                BAM2FASTQ_SCRIPT

                                                # write and submit bam to fastq job
                                                my $bam2fastq_job_path = "$genome_dir/$g/bam2fastq.sh";
    # save the job id to be used as a dependency for alignment job
    my $bam2fastq_job_id = submit_job $bam2fastq_script, $bam2fastq_job_path;

    # will store haploid job IDs for dependency purposs
    @haploid_jobs_ids = ();
    # for each haploid, make a reference, align the chipseq data
    for my $haploid ("maternal", "paternal")
    {
    mkdir "$genome_dir/$g/$haploid";
    my $assembly_script= <<"ASSEMBLY_SCRIPT";
    #PBS -A bws-221-ae
    #PBS -l nodes=1:ppn=1 -l mem=15gb
    #PBS -l walltime=3:00:00
    #PBS -d $pwd/$genome_dir/$g/$haploid
    #PBS -W depend=afterok:$impute_fasta_job_id
    #PBS -r n
    module load mugqic/picard mugqic/samtools mugqic/java mugqic/bwa

        java -Xmx14G -jar \$PICARD_HOME/picard.jar CreateSequenceDictionary R=$haploid.fa O=$haploid.dict
        samtools faidx $haploid.fa
        bwa index -a bwtsw $haploid.fa
        ASSEMBLY_SCRIPT
        # submit a job that will build a reference and a bwa_index for use in alignement with bwa
        my $make_ref_job_path = "$genome_dir/$g/$haploid/make_${haploid}_ref.sh";
    my $make_ref_job_id = submit_job $assembly_script, $make_ref_job_path;

    # script that will align to haploid and output a bam file sorted by read name
    my $haploid_alignment_script= <<"HAPLOID_ALIGNMENT_SCRIPT";
    #PBS -A bws-221-ae
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=12:00:00
    #PBS -d $pwd/$genome_dir/$g/
    #PBS -W depend=afterok:$make_ref_job_id:$bam2fastq_job_id
    #PBS -r n

    module load mugqic/bwa mugqic/samtools mugqic/picard

        mkdir ${haploid}_alignment
        bwa mem  \\
        -M -t 12 \\
        $haploid/$haploid.fa \\
        fastq/$sample_name.1.fastq \\
        | samtools view -q 20 -\@12 -h | samtools sort -\@12 -l 5 -o ${haploid}_alignment/$sample_name.allelic.bam -

    java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx14G -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
    REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
    TMP_DIR=/localscratch/ \\
    INPUT=${haploid}_alignment/$sample_name.allelic.bam \\
	OUTPUT=${haploid}_alignment/$sample_name.allelic.dup.bam \\
	METRICS_FILE=${haploid}_alignment/$sample_name.allelic.dup.bam.metrics \\
	MAX_RECORDS_IN_RAM=1000000 && rm ${haploid}_alignment/$sample_name.allelic.bam

bwa mem  \\
  -M -t 12 \\
  $haploid/$haploid.fa \\
  fastq/input.fastq.gz \\
  | samtools view -q 20 -\@12 -h | samtools sort -\@12 -l 5 -o ${haploid}_alignment/input.undup.bam -

java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx14G -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
	REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
	TMP_DIR=/localscratch/ \\
	INPUT=${haploid}_alignment/input.undup.bam \\
	OUTPUT=${haploid}_alignment/input.bam \\
	METRICS_FILE=${haploid}_alignment/input.bam.metrics \\
	MAX_RECORDS_IN_RAM=1000000 && rm ${haploid}_alignment/input.undup.bam
cd ${haploid} && rm *.amb *.ann *.bwt *.pac *.sa
HAPLOID_ALIGNMENT_SCRIPT

	# submit job and save haploid alignments ids
	my $haploid_alignment_job_path = "$genome_dir/$g/align_$haploid.sh";
	my $haploid_alignment_job_id = submit_job $haploid_alignment_script, $haploid_alignment_job_path;
	push @haploid_jobs_ids, $haploid_alignment_job_id;
}

my $reference_alignment_script= <<"REFERENCE_ALIGNMENT_SCRIPT";
#PBS -A bws-221-ae
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -d $pwd/$genome_dir/$g/
#PBS -W depend=afterok:$bam2fastq_job_id
#PBS -r n

module load mugqic/bwa mugqic/samtools mugqic/picard

mkdir reference_alignment
bwa mem  \\
  -M -t 12 \\
  $pwd/hg19_ordered/Homo.sapiens.hg19.good_contigs.fa\\
  fastq/$sample_name.1.fastq \\
  | samtools view -q 20 -\@12 -h | samtools sort -\@12 -l 5 -o reference_alignment/$sample_name.bam -

java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx14G -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
	REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
	TMP_DIR=/localscratch/ \\
	INPUT=reference_alignment/$sample_name.bam \\
	OUTPUT=reference_alignment/$sample_name.dup.bam \\
	METRICS_FILE=reference_alignment/$sample_name.dup.bam.metrics \\
	MAX_RECORDS_IN_RAM=1000000 && rm reference_alignment/$sample_name.bam

bwa mem  \\
  -M -t 12 \\
  $pwd/hg19_ordered/Homo.sapiens.hg19.good_contigs.fa \\
  fastq/input.fastq.gz \\
  | samtools view -q 20 -\@12 -h | samtools sort -\@12 -l 5 -o reference_alignment/input.undup.bam -

java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx14G -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
	REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
	TMP_DIR=/localscratch/ \\
	INPUT=reference_alignment/input.undup.bam \\
	OUTPUT=reference_alignment/input.bam \\
	METRICS_FILE=reference_alignment/input.bam.metrics \\
	MAX_RECORDS_IN_RAM=1000000 && rm reference_alignment/input.undup.bam
REFERENCE_ALIGNMENT_SCRIPT

my $reference_alignment_job_path = "$genome_dir/$g/align_reference.sh";
my $reference_alignment_job_id = submit_job $reference_alignment_script, $reference_alignment_job_path;
# join dependencies with qsub's : separator
my $peak_call_deps = join ":" , @haploid_jobs_ids;

# script to call peaks on bam files with marked duplicates
my $peak_call_script = <<"PEAK_CALL_SCRIPT";
#PBS -A bws-221-ae
#PBS -d $pwd/$genome_dir/$g/
#PBS -l nodes=1:ppn=1 -l mem=10gb
#PBS -W depend=afterok:$peak_call_deps:$reference_alignment_job_id
#PBS -l walltime=10:00:00
#PBS -r n

module load mugqic/python/2.7.13 mugqic/MACS2/2.1.0.20151222 mugqic/mugqic_R_packages
mkdir -p peak_call/paternal peak_call/maternal peak_call/reference
rm -rf fastq
macs2 callpeak --format BAM $broad --nomodel --gsize \$(Rscript $pwd/scripts/effectiveGenomeSize.R paternal/paternal.fa.fai) \\
	--treatment paternal_alignment/$sample_name.allelic.dup.bam --control paternal_alignment/input.bam --nolambda --name peak_call/paternal/paternal && rm paternal_alignment/input.bam 

macs2 callpeak --format BAM $broad --nomodel --gsize \$(Rscript $pwd/scripts/effectiveGenomeSize.R maternal/maternal.fa.fai) \\
	--treatment maternal_alignment/$sample_name.allelic.dup.bam --control maternal_alignment/input.bam --nolambda --name peak_call/maternal/maternal && rm maternal_alignment/input.bam 

macs2 callpeak --format BAM $broad --nomodel --gsize \$(Rscript $pwd/scripts/effectiveGenomeSize.R $pwd/hg19_ordered/Homo.sapiens.hg19.good_contigs.fa.fai) \\
	--treatment reference_alignment/$sample_name.dup.bam --control reference_alignment/input.bam  --nolambda --name peak_call/reference/reference && rm reference_alignment/input.bam 

if [[ "$broad" == "--broad" ]]; then
	echo Called Broad Peaks
else
	# For compatibility with the rest of the scripts, take the first 9 fields of the narrowPeak file and pipe them to a "broadPeak" format
	for peak_call in maternal paternal  reference
	do
		cat peak_call/\$peak_call/\${peak_call}_peaks.narrowPeak | cut -f1-9 > peak_call/\$peak_call/\${peak_call}_peaks.broadPeak
	done
fi
PEAK_CALL_SCRIPT

my $peak_call_job_path = "$genome_dir/$g/call_peaks.sh";
my $peak_call_job_id = submit_job $peak_call_script, $peak_call_job_path;

# produce counts, graphs, csv files
my $analyze_script = <<"ANALYZE_SCRIPT";
#PBS -A bws-221-ae
#PBS -d $pwd/$genome_dir/$g/
#PBS -l nodes=1:ppn=1 -l mem=30gb
#PBS -W depend=afterok:$peak_call_job_id
#PBS -l walltime=10:00:00
#PBS -r n
module load mugqic/bedtools mugqic/vcftools mugqic/ucsc
BASE_DIR=\$(pwd)
$scripts/bam2coverage.sh reference_alignment/$sample_name.dup.bam $pwd/hg19_ordered/hg19.chrom.sizes > reference_coverage.bed
rm -f reference_alignment/$sample_name.dup.bam

# for every haploid
for haploid in maternal paternal
do
echo For \$haploid
mkdir bed_comparison_\$haploid
echo \$(pwd)

$scripts/faiToChromSizes.pl  \$haploid/\$haploid.fa.fai > bed_comparison_\$haploid/Haploid.chrom.sizes
$scripts/bam2coverage.sh \${haploid}_alignment/$sample_name.allelic.dup.bam \$haploid/\$haploid.chrom.sizes > bed_comparison_\$haploid/HaploidCoverage.bed

$scripts/liftAnnotations.sh \$haploid $g $varcalls $chipseq_path $sample_name $pwd

cd bed_comparison_\$haploid/
echo Producing counts
# run split data
$scripts/splitData.sh \$haploid
cd \$BASE_DIR
rm -f \${haploid}_alignment/$sample_name.allelic.dup.bam
rm -rf \$haploid
done
rm -f reference_coverage.bed
ANALYZE_SCRIPT

my $analyze_job_path = "$genome_dir/$g/analyze.sh";
my $analyze_job_id = submit_job $analyze_script, $analyze_job_path;

}

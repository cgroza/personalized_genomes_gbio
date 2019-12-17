#!/usr/bin/perl
# This script follows the template found at: http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt
use feature qw(say);
use File::Basename;
use Cwd;

($orig2bit, $orig_build, $orig_chrom_sizes, $new2bit, $new_build, $new_chrom_sizes, $chain_out) = @ARGV;

# parameter explanation

# $orig2bit: 2bit file relative to $orig_build directory of original genome fasta
# $orig_build: directory holding files of original genome
# $orig_chrom_sizes: .chrom.sizes file relative to $orig_build directory

# $new2bit: 2bit file relative to $new_build directory of new genome fasta
# $new_build: directory holding files of new genome. Must also contain new fasta as contigs in separate files with chr* prefix
# $new_chrom_sizes: .chrom.sizes file relative to $new_build directory

# $chain_out: final chain filename

opendir ($orig_dir, $orig_build) or die "Cannot open $old_build";
opendir ($new_dir, $new_build) or die "Cannot open $new_build";

chdir($orig_build);
$work_dir=cwd();
print $work_dir . "\n";
# create work directories
map {mkdir $_} ("temp", "temp/lift", "temp/split", "temp/chr_chains", "temp/orig_psl", "temp/new_psl", "chain", "net", "over", "jobs");

@job_ids = ();

# create a job for each chromosome
while($new_chrom = readdir($new_dir))
    {
	print $new_chrom;
        $orig_chrom = "orig_$new_chrom";
        # skip non chromosome files
        next if(not $new_chrom =~ /^chr/);
        say  "For $new_chrom and $orig_chrom";
        $job_name = "makeChains_$new_chrom" . basename($orig_build);
        # write job for current chromosme
        open $chrom_script, ">" ,"jobs/$new_chrom.sh" or die "Cannot write jobs/$new_chrom.sh";
        $script = <<"BLAT_JOB";
#!/bin/bash

#PBS -A bws-221-ae
#PBS -d $work_dir
#PBS -l nodes=1:ppn=1 -l mem=5gb
#PBS -l walltime=72:00:00
#PBS -r n
#PBS -N $job_name
#PBS -e jobs/error$new_chrom
#PBS -o jobs/out$new_chrom

module load mugqic/ucsc


# faSplit for efficient blat and .lft file generation
echo "Splitting $new_chrom"
faSplit -lift=temp/lift/$new_chrom.lft size $new_build/$new_chrom -oneFile 3000 temp/split/$new_chrom

# blat each orig chromosome no new chromosome
echo "BLATing chromosomes"
blat $orig2bit temp/split/$new_chrom.fa temp/orig_psl/$orig_chrom.psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap

# lift up coordinates
echo "Lifting up to new coordinates"
liftUp -pslQ temp/new_psl/$new_chrom.psl temp/lift/$new_chrom.lft warn temp/orig_psl/$orig_chrom.psl

# make chromosome chains
echo "axtChaining to chain file"
axtChain -linearGap=medium -psl temp/new_psl/$new_chrom.psl $orig2bit $new2bit temp/chr_chains/$orig_chrom.chain
BLAT_JOB
        print $chrom_script $script;
        close chrom_script;

        # collect job id for dependency list from qsub stdout
        chomp($job_id =  `bqsub -q qfat256 "jobs/$new_chrom.sh"`);
        push @job_ids, $job_id;
    }

#script to merge all the chromosome chains

$jobstring = join ":", @job_ids;
say "Merge dependencies: $jobstring";
$merge_chain = <<"CHAIN_MERGE";
#!/bin/bash
#PBS -A bws-221-ae
#PBS -d $work_dir
#PBS -l nodes=1:ppn=1 -l mem=1gb
#PBS -l walltime=1:00:00
#PBS -r n
#PBS -e jobs/errorChainMerge
#PBS -o jobs/outChainMerge
#PBS -W depend=afterok:$jobstring

module load mugqic/ucsc
echo "Merging chains"
chainMergeSort temp/chr_chains/* | chainSplit chain stdin
echo "Running chainNet"

cd chain
for chrom_chain in *.chain
do
    echo \$chrom_chain
    chainNet \$chrom_chain ../$orig_chrom_sizes $new_chrom_sizes ../net/\$chrom_chain.net /dev/null
done


for chrom_chain in *.chain
do
    netChainSubset ../net/\$chrom_chain.net \$chrom_chain ../over/\$chrom_chain
done

cd ..

cat over/*.chain > $chain_out
CHAIN_MERGE

`bqsub --submit`;

# write merge script
# this file should be submited manually after all the jobs complete
# could be fixed if someome can tell me how job dependency works on mammouth
open $merge_script, ">", "jobs/merge_chains.sh" or die "Cannot write jobs/merge_chains.sh";
print $merge_script $merge_chain;
close $merge_script;

# run this step manually at the end
#qsub jobs/merge_chains.sh`;

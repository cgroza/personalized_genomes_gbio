# graph\_peak\_call.nf

This nextflow pipeline aligns ChIP-seq data to two genome graphs, calls peaks and surjects them back to the linear reference path. It takes the following parameters:
  - --design_file: a tab separated file listing the one treatment fastq and one control fastq per line
  - --ref_graph: directory holding the reference graph. Must contain a subidrectory _graphs_ that contains all vg files for each chromosome.
  - --pop_graph: directory holding the population graph. Must contain a subidrectory _graphs_ that contains all vg files for each chromosome.
  - --ref_name: name prefix of xg, gbwt, gcsa indexes of the reference graph
  - --pop_name: name prefix of xg, gbwt, gcsa indexes of the population graph

  - --genome_size: size in basepairs of the species' genome
  - --chromosomes:  comma separated list of chromosome names
  - --fragment_length: fragment lenght of ChIP-seq library
  - --qvlaue: false discvery rate for peak callling. Ranges between 0 and 1
  - --paired: are fastqs interleaved paired-end reads? Can be true or false
  - --peak_call: perform peak calling step? Can be true of false
  - --altered: substract reference and population graph annoatations. Can be true or false
  - --outDir: output directory
  - --time: maximum walltime for job scheduler during alignment steps
  - --mem: maximum memory allocation for job scheduler during alignment steps

# pop_graph.nf

This is my custom nextflow script for creating population graphs.
You should use the _toil-vg construct_ pipeline instead from the developers of _vg_.

# makeChains.pl

Creates a liftover file between two genomes. Unfortunately, this perl code is not readily portable.

# submit_jobs_MPG and scripts
Perl pipeline for calling personalized peaks in modified reference genomes. It consists of bash and perl code that is not readily portable.

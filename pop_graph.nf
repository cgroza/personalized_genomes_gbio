// path to the reference genome fasta
params.ref = "hg19.fa"
// name of reference genome
params.genome = "POP"
// output directory
params.outdir = workflow.launchDir
// path to variants VCF.gz file
params.vcf = "POP.vcf.gz"

Channel.fromPath(params.vcf).into{vcf_con_ch; vcf_gbwt_ch}

process makeVg {
    time '12h'
    memory '100 GB'
    cpus 6

    publishDir "$params.outdir/graphs", mode: 'copy', pattern: "*.vg"

    input:
    file vcf from vcf_con_ch

    output:
    file "graphs/*.vg" into vgs_ch_gbwt, vgs_ch_gcsa
    file "*.tbi" into vcf_index_ch
    file "*.vcf.gz" into vcf_ch
    file "mapping" into mapping_ch

    script:
    """
module load tabix
tabix -p vcf $vcf
(seq 1 22; echo X; echo Y) | parallel -j 24 "tabix -h $vcf chr{} > chr{}.vcf ; bgzip chr{}.vcf ; tabix chr{}.vcf.gz"

mkdir graphs
(seq 1 22; echo X; echo Y) | parallel -j 6  "vg construct -S -a -p -C -R chr{} -v chr{}.vcf.gz -r $params.ref -t 1 -m 32 > graphs/chr{}.vg"
vg ids -m mapping -j \$(for i in \$(seq 1 22; echo X; echo Y); do echo graphs/chr\$i.vg; done)
"""
}

process indexGBWT_XG {
    cpus 40
    time '2d'
    memory '150 GB'
    publishDir "$params.outdir", mode: 'copy', pattern: "${params.genome}_index.gbwt"
    publishDir "$params.outdir", mode: 'copy', pattern: "${params.genome}_index.xg"

    input:
    file "*" from vgs_ch_gbwt.collect()
    file "*" from vcf_ch.collect()
    file "*" from vcf_index_ch.collect()

    output:
    file "${params.genome}_index.gbwt"
    file "${params.genome}_index.xg" into xg_ch
    file "*.gbwt" into gbwt_ch

    script:
    """
TMPDIR=/home/cgroza/scratch/temp
(seq 1 22; echo X; echo Y) | parallel -j 8 "touch -h chr{}.vcf.gz.tbi ; vg index -G chr{}.gbwt -v chr{}.vcf.gz chr{}.vg"
vg gbwt -m -f -o ${params.genome}_index.gbwt chr*.gbwt
vg index -x ${params.genome}_index.xg *.vg
"""
}


process indexGCSA {
    cpus 40
    time '2d'
    memory '180 GB'
    publishDir "$params.outdir", mode: 'copy'

    input:
    file "*" from vgs_ch_gcsa.collect()
    file "*" from gbwt_ch.collect()
    file mapping from mapping_ch

    output:
    file "${params.genome}_index.gcsa"
    file "${params.genome}_index.gcsa.lcp"

    script:
    """
mkdir graphs
TMPDIR=/home/cgroza/scratch/temp
cp ${mapping} mapping.backup
for i in \$(seq 1 22; echo X; echo Y); do
    vg prune -a -m mapping.backup -u -g chr\${i}.gbwt chr\${i}.vg > graphs/chr\${i}.pruned.vg
done

vg index -g ${params.genome}_index.gcsa -f mapping.backup graphs/*.vg
"""
}

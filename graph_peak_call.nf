params.design_file = "design.tsv"
params.ref_graph = "ref/"
params.pop_graph = "pop/"
params.ref_name = "ref"
params.pop_name = "pop"

params.genome_size = 3100000000
params.fragment_length = 200
params.qvalue = 0.05
params.paired = true
params.time = '60h'
params.mem = '100 GB'
params.sort = false
params.peak_call = true
params.altered = true

params.outDir = workflow.launchDir

chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def design_file = new File(params.design_file)
def design = [:]

def treatments = []
def controls = []

design_file.eachLine {String entry ->
    def (treatment, control) = entry.split()
    treatment_name = (new File(treatment)).getName();
    control_name = (new File(control)).getName();

    controls.add(control)
    treatments.add(treatment)
    // the design mapping outputs the control for each treatment
    design[treatment_name] = control_name
    // the design mapping is identity for control
    design[control_name] = control_name
}

Channel.fromPath("${params.pop_graph}/graphs/*.vg").set{linear_vg_ch}
Channel.fromPath("${params.ref_graph}/graphs/*.vg").set{ref_linear_vg_ch}
Channel.fromPath(treatments).unique().into{fastq_ch; ref_fastq_ch}
Channel.fromPath(controls).unique().into{control_fastq_ch; ref_control_fastq_ch}

println("Treatments:")
println(treatments)
println("Controls:")
println(controls)

println("Treatment/control associations:")
println(design)


Channel.fromPath(
    ["${params.ref_graph}/${params.ref_name}.xg",
     "${params.ref_graph}/${params.ref_name}.gcsa",
     "${params.ref_graph}/${params.ref_name}.gcsa.lcp"]).into{ref_index_treatment_ch; ref_index_control_ch}

Channel.fromPath(
    ["${params.pop_graph}/${params.pop_name}.xg",
     "${params.pop_graph}/${params.pop_name}.gbwt",
     "${params.pop_graph}/${params.pop_name}.gcsa",
     "${params.pop_graph}/${params.pop_name}.gcsa.lcp"]).into{pop_index_treatment_ch; pop_index_control_ch}


vg_flag = ""
if(params.paired) {
    vg_flag = "-i ${vg_flag}"
}

process vgToJsonPop {
    cpus = 40
    memory '120 GB'
    time '24h'

    input:
    file "graphs/*" from linear_vg_ch.collect()

    output:
    file "graphs/*" into vg2json_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -Vj graphs/chr{}.vg > graphs/chr{}.json ; graph_peak_caller create_ob_graph graphs/chr{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_chr{}.txt'
"""
}

process vgToJsonRef {
    cpus = 40
    memory '120 GB'
    time '24h'

    input:
    file "graphs/*" from ref_linear_vg_ch.collect()

    output:
    file "graphs/*" into ref_vg2json_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -Vj graphs/chr{}.vg > graphs/chr{}.json ; graph_peak_caller create_ob_graph graphs/chr{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_chr{}.txt'
"""
}

process linearPathsPop {
    cpus = 40
    memory '120 GB'
    time '24 h'

    input:
    file "graphs/*" from vg2json_ch.collect()

    output:
    file "graphs" into control_linear_ch, treatment_linear_ch, peak_linear_ch

    script:
    """
   (seq 1 22; echo X; echo Y) | parallel -j 3 graph_peak_caller find_linear_path -g graphs/chr{}.nobg graphs/chr{}.json chr{} graphs/chr{}_linear_pathv2.interval
"""
}


process linearPathsRef {
    cpus = 40
    memory '120 GB'
    time '24 h'

    input:
    file "graphs/*" from ref_vg2json_ch.collect()

    output:
    file "graphs" into ref_control_linear_ch, ref_treatment_linear_ch, ref_peak_linear_ch

    script:
    """
     (seq 1 22; echo X; echo Y) | parallel -j 3 graph_peak_caller find_linear_path -g graphs/chr{}.nobg graphs/chr{}.json chr{} graphs/chr{}_linear_pathv2.interval
"""
}

process alignControlRef {
    cpus = 40
    memory "${params.mem}"
    time = "${params.time}"

    input:
    set file(xg), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs") from ref_index_control_ch.collect().combine(ref_control_fastq_ch).combine(ref_control_linear_ch).view()

    output:
    set file(fastq), file("control_json") into ref_control_json_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir control_gam
    vg map $vg_flag -f $fastq -x $xg -g $gcsa -t 40 -u 1 -m 1 > "control_gam/${name}_ref.gam"

    mkdir control_json
    vg view -aj control_gam/${name}_ref.gam > control_json/${name}_ref.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} control_json/${name}_ref.json graphs/
    rm control_json/${name}_ref.json
"""
}


process alignControlPop {
    cpus = 40
    memory "${params.mem}"
    time = "${params.time}"

    input:
    set file(xg), file(gbwt), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs") from pop_index_control_ch.collect().combine(control_fastq_ch).combine(control_linear_ch).view()
    output:
    set file(fastq), file("control_json") into control_json_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir control_gam
    vg map $vg_flag -f $fastq -1 $gbwt -x $xg -g $gcsa -t 40 -u 1 -m 1 > "control_gam/${name}_pop.gam"

    mkdir control_json
    vg view -aj control_gam/${name}_pop.gam > control_json/${name}_pop.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} control_json/${name}_pop.json graphs/
    rm control_json/${name}_pop.json
"""
}

process alignSampleRef {
    cpus = 40
    memory "${params.mem}"
    time = "${params.time}"

    input:
    set file(xg), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs") from ref_index_treatment_ch.collect().combine(ref_fastq_ch).combine(ref_treatment_linear_ch).view()

    output:
    set file(fastq), file("json") into ref_treatment_json_ch
    file("gam/${name}_ref.gam") into ref_treatment_gam_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir gam
    vg map $vg_flag -f $fastq -x $xg -g $gcsa -t 40 -u 1 -m 1 > gam/${name}_ref.gam

    mkdir json
    vg view -aj gam/${name}_ref.gam > json/${name}_ref.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}_ref.json graphs/
    rm json/${name}_ref.json
"""
}


process alignSamplePop {
    cpus = 40
    memory "${params.mem}"
    time = "${params.time}"

    input:
    set file(xg), file(gbwt), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs") from pop_index_treatment_ch.collect().combine(fastq_ch).combine(treatment_linear_ch).view()

    output:
    set file(fastq), file("json") into treatment_json_ch
    file("gam/${name}_pop.gam") into treatment_gam_ch

    script:
    name = fastq.getSimpleName()

    """
    mkdir gam
    vg map $vg_flag -f $fastq -1 $gbwt -x $xg -g $gcsa -t 40 -u 1 -m 1 > gam/${name}_pop.gam

    mkdir json
    vg view -aj gam/${name}_pop.gam > json/${name}_pop.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}_pop.json graphs/
    rm json/${name}_pop.json
"""
}

if(params.sort) {
    process sortSampleRef {
        cpus = 40
        memory '100 GB'
        time = '12h'

        publishDir "$params.outDir", pattern: "ref_${name}.sorted.gam"

        input:
        file(gam) from ref_treatment_gam_ch

        output:
        file "ref_${name}.sorted.gam"
        script:
        name = gam.getSimpleName()
        """
        vg gamsort ${gam} -i ${name}.sorted.gam.gai -t 40 > ref_${name}.sorted.gam
    """
    }

    process sortSamplePop {
        cpus = 40
        memory '100 GB'
        time '12h'

        publishDir "$params.outDir", pattern: "${name}.sorted.gam"

        input:
        file(gam) from treatment_gam_ch

        output:
        file "${name}.sorted.gam"
        script:
        name = gam.getSimpleName()
        """
        vg gamsort ${gam} -i ${name}.sorted.gam.gai -t 40 > ${name}.sorted.gam
    """
    }
}

if(params.peak_call) {
    process callPeaksPop{
        cpus = 40
        memory '120 GB'
        time '24h'

        publishDir "$params.outDir/peaks", pattern: "${name}_peaks.narrowPeak", mode: "copy"

        input:
        set file(fastq), file("json"), file(control_fastq), file("control_json"), file("graphs") from treatment_json_ch.combine(control_json_ch).filter{design[it.get(0).getName()] == design[it.get(2).getName()]}.combine(peak_linear_ch).map{ it.flatten()}.view()

        output:
        set val(name), file("${name}_peaks.narrowPeak") into pop_peaks_ch

        script:
        name = fastq.getSimpleName()
        control_name = control_fastq.getSimpleName()

        """
        (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller count_unique_reads chr{} graphs/ json/${name}_pop_ | tail -n 1 > counted_unique_reads_chr{}.txt'
        read_length=\$(zcat $fastq | head -2 | tail -1 | wc -c)
        unique_reads=\$(awk 'BEGIN{i=0}{i = i + \$1}END{print i}' counted_unique_reads_chr*.txt)
        (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks -q ${params.qvalue} -g graphs/chr{}.nobg -s json/${name}_pop_chr{}.json -c control_json/${control_name}_pop_chr{}.json -G ${params.genome_size} -p True -f ${params.fragment_length} -r \$read_length -u \$unique_reads -n chr{}"
        rename 'touched' '_touched' *touched*
        rename 'background' '_background' *background*
        rename 'direct' '_direct' *direct*
        rename 'fragment' '_fragment' *fragment*
        rename 'pvalues' '_pvalues' *pvalues*

        (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks_whole_genome_from_p_values -q ${params.qvalue} -d graphs/ -n '' -f ${params.fragment_length} -r \$read_length chr{}"
        (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller peaks_to_linear chr{}_max_paths.intervalcollection graphs/chr{}_linear_pathv2.interval chr{} chr{}_linear_peaks.bed'
        cat *_linear_peaks.bed | sort-bed - > ${name}_peaks.narrowPeak
    """
    }
    process callPeaksRef{
        cpus = 40
        memory '120 GB'
        time '24h'

        publishDir "$params.outDir/peaks", pattern: "ref_${name}_peaks.narrowPeak", mode: "copy"

        input:
        set file(fastq), file("json"), file(control_fastq), file("control_json"), file("graphs") from ref_treatment_json_ch.combine(ref_control_json_ch).filter{design[it.get(0).getName()] == design[it.get(2).getName()]}.combine(ref_peak_linear_ch).map{ it.flatten()}.view()

        output:
        set val(name), file("ref_${name}_peaks.narrowPeak") into ref_peaks_ch

        script:
        name = fastq.getSimpleName()
        control_name = control_fastq.getSimpleName()
        """
        (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller count_unique_reads chr{} graphs/ json/${name}_ref_ | tail -n 1 > counted_unique_reads_chr{}.txt'
        read_length=\$(zcat $fastq | head -2 | tail -1 | wc -c)
        unique_reads=\$(awk 'BEGIN{i=0}{i = i + \$1}END{print i}' counted_unique_reads_chr*.txt)
        (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks -q ${params.qvalue} -g graphs/chr{}.nobg -s json/${name}_ref_chr{}.json -c control_json/${control_name}_ref_chr{}.json -G ${params.genome_size} -p True -f ${params.fragment_length} -r \$read_length -u \$unique_reads -n chr{}"
        rename 'touched' '_touched' *touched*
        rename 'background' '_background' *background*
        rename 'direct' '_direct' *direct*
        rename 'fragment' '_fragment' *fragment*
        rename 'pvalues' '_pvalues' *pvalues*

        (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks_whole_genome_from_p_values -q ${params.qvalue} -d graphs/ -n '' -f ${params.fragment_length} -r \$read_length chr{}"
        (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller peaks_to_linear chr{}_max_paths.intervalcollection graphs/chr{}_linear_pathv2.interval chr{} chr{}_linear_peaks.bed'
        cat *_linear_peaks.bed | sort-bed - > ref_${name}_peaks.narrowPeak
    """
    }
}

if(params.altered){
    process alteredPeaks {
        cpus = 1
        publishDir "$params.outDir/peaks", mode: "copy"
        input:
        set val(name), file(pop_peaks), val(ref_name), file(ref_peaks) from pop_peaks_ch.combine(ref_peaks_ch, by: 0).view()
        output:
        file "*.narrowPeak"

        script:
        """
    module load bedtools
    bedtools subtract -A -a ${pop_peaks} -b ${ref_peaks} > ${name}_pers-only.narrowPeak
    bedtools subtract -A -b ${pop_peaks} -a ${ref_peaks} > ${name}_ref-only.narrowPeak
    bedtools intersect -wa -a ${pop_peaks} -b ${ref_peaks} > ${name}_intersected.narrowPeak
    """
    }
}

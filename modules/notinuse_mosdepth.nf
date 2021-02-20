process MOSDEPTH_GENOME {
    tag "$sampleName"

    input:
    tuple val(sampleName), path(bam), path(bai)
    val outdir

    output:
    tuple val(sampleName), path ("${sampleName}.trim.mkD.sorted.bam"), path ("${sampleName}.trim.mkD.sorted.bam.bai"), emit: PICARD_MARKDUP_out
    path "*.{log,sh}"

    script:
    """
    mosdepth \
        --by 200 \
        --fast-mode \
        ${sampleName}.trim.mkD.genome \
        ${sampleName}.trim.mkD.sorted.bam

    plot_mosdepth_regions.r \
        --input_files ${sampleName}.trim.mkD.genome.regions.bed.gz \
        --input_suffix .trim.mkD.genome.regions.bed.gz \
        --output_dir ./ \
        --output_suffix .trim.mkD.genome.regions
    
    cp .command.sh ${sampleName}.mosdepth_genome.sh
    cp .command.log ${sampleName}.mosdepth_genome.log
    """
}

process MOSDEPTH_AMPLICON {
    tag "$sampleName"
    publishDir "$outdir", mode:'link',  overwrite: false

    input:
    tuple val(sampleName), path(bam), path(bai)
    path primer_bed
    val outdir

    output:
    tuple val(sampleName), path ("${sampleName}.trim.mkD.sorted.bam"), path ("${sampleName}.trim.mkD.sorted.bam.bai"), emit: PICARD_MARKDUP_out
    path "*.{metrics.txt,log,sh}"

    script:
    """
    collapse_amplicon_bed.py \
        --left_primer_suffix _LEFT \
        --right_primer_suffix _RIGHT \
        $primer_bed \
        amplicon.collapsed.bed

    mosdepth \
        --by amplicon.collapsed.bed \
        --fast-mode \
        --use-median \
        --thresholds 0,1,10,50,100,500 \
        ${sampleName}.trim.mkD.amplicon \
        ${sampleName}.trim.mkD.sorted.bam

    plot_mosdepth_regions.r \
        --input_files ${sampleName}.trim.mkD.amplicon.regions.bed.gz \
        --input_suffix .trim.mkD.amplicon.regions.bed.gz \
        --output_dir ./ \
        --output_suffix .trim.mkD.amplicon.regions

    cp .command.sh ${sampleName}.mosdepth_amplicon.sh
    cp .command.log ${sampleName}.mosdepth_amplicon.log
    """
}


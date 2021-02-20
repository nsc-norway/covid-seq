process SAMTOOLS_MPILEUP {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(bam), path(bai)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}.trim.mpileup"), emit: SAMTOOLS_MPILEUP_out
    path "*.{log,sh}"

    //publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{mpileup}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    samtools mpileup \
        --count-orphans \
        --no-BAQ \
        --max-depth 0 \
        --fasta-ref $genome \
        --min-BQ 20 \
        --output ${sampleName}.trim.mpileup \
        $bam

    cp .command.sh ${sampleName}.samtools_mpileup.sh
    cp .command.log ${sampleName}.samtools_mpileup.log
    """
}


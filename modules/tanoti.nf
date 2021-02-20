process TANOTI {
    tag "$sampleName"

    label 'medium'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome

    publishDir "${params.outdir}/../for_FHI/2_bam_tanoti", mode: 'link', pattern:'*.bam'
    publishDir "${params.outdir}/../for_FHI/2_bam_tanoti/log", mode: 'link', pattern:'*.{log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.tanoti.bam"), emit: TANOTI_out
    path "*.{log,sh}"

    script:
    """
    gunzip -c ${read1} > ${sampleName}_R1.fastq
    gunzip -c ${read2} > ${sampleName}_R2.fastq
    tanoti -r $genome -i ${sampleName}_R1.fastq ${sampleName}_R2.fastq -o ${sampleName}.sam -p 1 -m 85

    samtools view -O bam -o ${sampleName}.tanoti.bam ${sampleName}.sam

    cp .command.sh ${sampleName}.tanoti.sh
    cp .command.log ${sampleName}.tanoti.log
    """
}
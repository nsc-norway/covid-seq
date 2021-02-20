process PRIMERCLIP {
    tag "$sampleName"

    label 'medium'

    input:
    tuple val(sampleName), path(bam), path(QC)
    path primer_master_file

    output:
    tuple val(sampleName), path ("${sampleName}.trim.sorted.bam"), path ("${sampleName}.trim.sorted.bam.bai"), emit: PRIMERCLIP_out
    path "*.{stats,log,sh}"

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam}*'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'

    script:
    """
    samtools sort -@ $task.cpus -n -O sam $bam > ${sampleName}.coord.sam
    primerclip $primer_master_file ${sampleName}.coord.sam ${sampleName}.trim.sam
    samtools sort -@ ${task.cpus} ${sampleName}.trim.sam -o ${sampleName}.trim.sorted.bam 
    samtools index ${sampleName}.trim.sorted.bam

    samtools stats ${sampleName}.trim.sorted.bam > ${sampleName}.trim.sorted.bam.stats

    mv masterparsefails.log ${sampleName}.primerclip_trim.masterparsefails.log
    cp .command.sh ${sampleName}.primerclip_trim.sh
    cp .command.log ${sampleName}.primerclip_trim.log
    """
}

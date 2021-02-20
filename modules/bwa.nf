process BWA_INDEX {
    tag "$genome"

    label 'small'

    input:
    path genome

    output:
    path "${genome}.*", emit: BWA_INDEX_out
    path "*.{log,sh}"

    script:
    """
    bwa index $genome

    cp .command.log bwa_index.log
    cp .command.sh bwa_index.sh
    """
}

process BWA_ALIGN {
    tag "$sampleName"
    errorStrategy 'ignore'
    
    label 'big'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'
    publishDir "${params.outdir}/../for_FHI/2_bam_bwa", mode: 'link', pattern:'*.bwa.bam'
    publishDir "${params.outdir}/../for_FHI/2_bam_bwa/log", mode: 'link', pattern:'*.{log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.bam"), path ("${sampleName}_QC_PASS"), emit: BWA_ALIGN_out
    tuple val(sampleName), path ("${sampleName}.bwa.bam")
    path "*.{stats,log,sh}"

    script:
    """
    bwa mem \
        -t $task.cpus \
        $genome \
        ${read1} ${read2} \
        2> ${sampleName}.bwa.log \
        | samtools view -O bam -o ${sampleName}.bam  -

    samtools stats ${sampleName}.bam > ${sampleName}.bam.stats

    SEQ=\$(grep 'raw total sequences:' ${sampleName}.bam.stats | cut -f3)
    if [ \$SEQ -gt 10000 ]
    then
        echo 'PASS' > ${sampleName}_QC_PASS
    else
        echo 'FAIL' > ${sampleName}_QC_FAIL
    fi

    cp ${sampleName}.bam ${sampleName}.bwa.bam

    cp .command.sh ${sampleName}.bwa.sh
    cp .command.out ${sampleName}.bwa.log
    """
}


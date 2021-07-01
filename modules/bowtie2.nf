process BOWTIE2_INDEX {
    tag "$genome"

    label 'small'

    input:
    path genome

    output:
    path "${genome}.*", emit: BOWTIE2_INDEX_out
    path "*.{log,sh}"

    script:
    """
    bowtie2-build \
        --seed 1 \
        $genome \
        $genome

    cp .command.log bowtie2_index.log
    cp .command.sh bowtie2_index.sh
    """
}

process BOWTIE2_ALIGN {
    tag "$sampleName"
    errorStrategy 'ignore'
    
    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai"), emit: BOWTIE2_ALIGN_out
    path "*.log", emit: BOWTIE2_log
    path "${sampleName}_QC_PASS"
    path "*.{stats,sh}"

    script:
    """
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        --rg-id NSC-NovaSeq-${sampleName} --rg SM:${sampleName} --rg LB:DUMMYLIB --rg PI:400 --rg PL:ILLUMINA \
        -x $genome \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'raw total sequences:' ${sampleName}.sorted.bam.stats | cut -f3)
    if [ \$SEQ -gt 10000 ]
    then
        echo 'PASS' > ${sampleName}_QC_PASS
    else
        echo 'FAIL' > ${sampleName}_QC_FAIL
    fi

    cp .command.sh ${sampleName}.bowtie2.sh
    """
}


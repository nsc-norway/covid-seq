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
    
    label 'big'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'
    if (params.lab == 'FHI') {
        publishDir "${params.outdir}/../for_FHI/2_bam_bowtie2", mode: 'link', pattern:'*.bowtie2.bam'
        publishDir "${params.outdir}/../for_FHI/2_bam_bowtie2/log", mode: 'link', pattern:'*.{log,sh}'
    }
    output:
    tuple val(sampleName), path ("${sampleName}.bam"), path ("${sampleName}_QC_PASS"), emit: BOWTIE2_ALIGN_out
    tuple val(sampleName), path ("${sampleName}.bowtie2.bam")
    path "*.{stats,log,sh}"

    script:
    """
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x $genome \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.log \
        | samtools view -O bam -o ${sampleName}.bam  -

    samtools stats ${sampleName}.bam > ${sampleName}.bam.stats

    SEQ=\$(grep 'raw total sequences:' ${sampleName}.bam.stats | cut -f3)
    if [ \$SEQ -gt 10000 ]
    then
        echo 'PASS' > ${sampleName}_QC_PASS
    else
        echo 'FAIL' > ${sampleName}_QC_FAIL
    fi

    cp ${sampleName}.bam ${sampleName}.bowtie2.bam

    cp .command.sh ${sampleName}.bowtie2.sh
    """
}


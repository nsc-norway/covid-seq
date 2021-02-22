process BCFTOOLS_VARIANTS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path (bam), path (bai)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_bcftools.vcf.gz"), path ("${sampleName}_bcftools.vcf.gz.tbi"), emit: BCFTOOLS_VARIANTS_out
    path "*.{vcf,txt,log,sh}"

    publishDir "${params.outdir}/3_variants/bcftools", mode: 'link', pattern:'*.{vcf}'
    publishDir "${params.outdir}/3_variants/bcftools/log", mode: 'link', pattern:'*.{txt,log,sh}'

    script:
    """
    echo "${sampleName}_bcftools" > sample_name.list
    bcftools mpileup \
        --count-orphans \
        --no-BAQ \
        --fasta-ref $genome \
        --min-BQ 20 \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        $bam \
        | bcftools call --output-type v --ploidy 1 --keep-masked-ref --multiallelic-caller --variants-only \
        | bcftools reheader --samples sample_name.list \
        | bcftools view --output-file ${sampleName}_bcftools.vcf.gz --output-type z --include 'INFO/DP>=10'
    tabix -p vcf -f ${sampleName}_bcftools.vcf.gz

    bgzip -d -c ${sampleName}_bcftools.vcf.gz > ${sampleName}_bcftools.vcf

    bcftools stats ${sampleName}_bcftools.vcf.gz > ${sampleName}_bcftools.bcftools_stats.txt

    cp .command.sh ${sampleName}.bcftools_variants.sh
    cp .command.log ${sampleName}.bcftools_variants.log
    """
}



process BCFTOOLS_CONSENSUS {
    tag "$sampleName"

    label 'small'
    errorStrategy 'ignore'

    input:
    tuple val(sampleName), path (bam), path (bai), path (vcf), path (tbi)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_bcftools.consensus.masked.fa"), emit: BCFTOOLS_CONSENSUS_out
    tuple val(sampleName), path ("${sampleName}_bcftools.consensus.masked.fa"), emit: FOR_LINEAGE_out
    path "*.{log,sh}"

    publishDir "${params.outdir}/4_consensus/bcftools", mode: 'link', pattern:'*.fa'
    publishDir "${params.outdir}/4_consensus/bcftools/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    bedtools genomecov \\
        -bga \\
        -ibam $bam \\
        -g $genome \\
        | awk '\$4 < 10' | bedtools merge > ${sampleName}.pre_mask.bed

    parse_mask_bed.py $vcf ${sampleName}.pre_mask.bed ${sampleName}.mask.bed
    bedtools maskfasta \\
        -fi $genome \\
        -bed ${sampleName}.mask.bed \\
        -fo ${sampleName}_bcftools.ref.masked.fa
    cat ${sampleName}_bcftools.ref.masked.fa | bcftools consensus $vcf > ${sampleName}_bcftools.consensus.masked.fa

    header=\$(head -n1 ${sampleName}_bcftools.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_bcftools_masked/g" ${sampleName}_bcftools.consensus.masked.fa

    cp .command.sh ${sampleName}.bcftools_consensus.sh
    cp .command.log ${sampleName}.bcftools_consensus.log
    """
}



/*
  bedtools genomecov \\
        -bga \\
        -ibam ${bam[0]} \\
        -g $fasta \\
        | awk '\$4 < $params.min_coverage' | bedtools merge > ${sample}.merged.bed
    parse_mask_bed.py ${vcf[0]} ${sample}.merged.bed ${sample}.mask.bed
    bedtools maskfasta \\
        -fi $fasta \\
        -bed ${sample}.mask.bed \\
        -fo ${index_base}.ref.masked.fa
    cat ${index_base}.ref.masked.fa | bcftools consensus ${vcf[0]} > ${sample}.consensus.masked.fa
    sed -i 's/${index_base}/${sample}/g' ${sample}.consensus.masked.fa
    header=\$(head -n1 ${sample}.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sample}/g" ${sample}.consensus.masked.fa

*/


/*
   cat $genome | bcftools consensus ${sampleName}_bcftools.vcf.gz > ${sampleName}_bcftools.consensus.fa

    bedtools genomecov \
        -bga \
        -ibam $bam \
        -g $genome \
        | awk '\$4 < 10' | bedtools merge > ${sampleName}_bcftools.mask.bed

    bedtools maskfasta \
        -fi ${sampleName}_bcftools.consensus.fa \
        -bed ${sampleName}_bcftools.mask.bed \
        -fo ${sampleName}_bcftools.consensus.masked.fa
    sed -i 's/${params.ref_id}/${sampleName}/g' ${sampleName}_bcftools.consensus.fa
    header=\$(head -n1 ${sampleName}_bcftools.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_bcftools/g" ${sampleName}_bcftools.consensus.masked.fa

*/
process HAPLOTYPECALLER_VARIANTS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path(bam), path(bai)
    path genome
    path genomeIdx

    output:
    tuple val(sampleName), path ("${sampleName}_hapcaller.vcf.gz"), path ("${sampleName}_hapcaller.vcf.gz.tbi"), emit: HAPCALLER_VARIANTS_out
    //path "*.{vcf,txt,log,sh}"
    path "*.{vcf,log,sh}"

    publishDir "${params.outdir}/3_variants/haplotypecaller", mode: 'link', pattern:'*.vcf'
    publishDir "${params.outdir}/3_variants/haplotypecaller/log", mode: 'link', pattern:'*.{log,sh}'
    //publishDir "${params.outdir}/3_variants/haplotypecaller/log", mode: 'link', pattern:'*.{txt,log,sh}'

    script:
    """
    gatk --java-options "-Xmx${task.memory.giga}g" HaplotypeCaller \
        -R $genome \
        -I $bam \
        -O ${sampleName}_hapcaller.vcf.gz \
        --create-output-variant-index \
        --ploidy 1
    gunzip -c ${sampleName}_hapcaller.vcf.gz > ${sampleName}_hapcaller.vcf

    #bcftools stats ${sampleName}_hapcaller.vcf.gz > ${sampleName}_hapcaller.bcftools_stats.txt

    cp .command.sh ${sampleName}.hapcaller_variants.sh
    cp .command.log ${sampleName}.hapcaller_variants.log
    """
}



process HAPLOTYPECALLER_CONSENSUS {
    tag "$sampleName"

    label 'small'

    input:
    tuple val(sampleName), path (bam), path (bai), path (vcf), path (tbi)
    path genome

    output:
    tuple val(sampleName), path ("${sampleName}_hapcaller.consensus.masked.fa"), emit: HAPCALLER_CONSENSUS_out
    tuple val(sampleName), path ("${sampleName}_hapcaller.consensus.masked_Nremoved.fa")
    path "*.{log,sh}"

    publishDir "${params.outdir}/4_consensus/haplotypecaller", mode: 'link', pattern:'*.fa'
    publishDir "${params.outdir}/4_consensus/haplotypecaller/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    bedtools genomecov \
        -bga \
        -ibam $bam \
        -g $genome \
        | awk '\$4 < 10' | bedtools merge > ${sampleName}.pre_mask.bed


    bcftools view -i '2*FORMAT/AD[0:1] >= FORMAT/DP[0]' -O z -o ${sampleName}_variantsForConsensus.vcf.gz $vcf
    tabix -p vcf -f ${sampleName}_variantsForConsensus.vcf.gz

    parse_mask_bed.py ${sampleName}_variantsForConsensus.vcf.gz ${sampleName}.pre_mask.bed ${sampleName}.mask.bed
    
    bedtools maskfasta \
        -fi $genome \
        -bed ${sampleName}.mask.bed \
        -fo ${sampleName}_hapcaller.ref.masked.fa

    cat ${sampleName}_hapcaller.ref.masked.fa \
        | bcftools consensus ${sampleName}_variantsForConsensus.vcf.gz > ${sampleName}_hapcaller.consensus.masked.fa

    header=\$(head -n 1 ${sampleName}_hapcaller.consensus.masked.fa | sed 's/>//g')
    sed -i "s/\${header}/${sampleName}_hapcaller_masked/g" ${sampleName}_hapcaller.consensus.masked.fa

    sed -r '2s/^N{1,}//g' ${sampleName}_hapcaller.consensus.masked.fa \
            | sed -r '\$ s/N{1,}\$//g' \
            > ${sampleName}_hapcaller.consensus.masked_Nremoved.fa

    cp .command.sh ${sampleName}.hapcaller_consensus.sh
    cp .command.log ${sampleName}.hapcaller_consensus.log
    """
}

process NOISE_EXTRACTOR {

    label 'medium'

    input:
    path 'input_bam/*'

    output:
    path '*.{xlsx,pdf}', emit: NOISE_SUMMARY_FILES_out
    path "rawnoise/*.tsv"
    path '*.{log,sh}'

    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'rawnoise/*'
    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'*.{xlsx,pdf}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /home/docker/Scripts/CSAK_NoiseExtractor_docker.R c$task.cpus

    # Move results into work dir, for output
    mv input_bam/*.{xlsx,pdf} .
    mv input_bam/rawnoise .

    cp .command.sh all.noiseextractor.sh
    cp .command.log all.noiseextractor.log
    """
}

process FRAMESHIFT_FINDER {

    label 'large'
    errorStrategy 'ignore'

    input:
    path inputFasta

    output:
    path "*.csv"
    path "*.xlsx"
    path "*.{sh,log}"

    publishDir "${params.outdir}/4_consensus/frameshiftfinder", mode: 'link', pattern: '*.{xlsx,csv}'
    publishDir "${params.outdir}/4_consensus/frameshiftfinder", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /home/docker/Scripts/CSAK_Frameshift_Finder_docker.R c$task.cpus
    
    cp .command.sh all.frameshiftfinder.sh
    cp .command.log all.frameshiftfinder.log
    """
}

process NEXTCLADE_ANALYSIS {
    label 'large'

    input:
    path allSamplesConsensus
    val(caller)

    output:
    path "all_${caller}_Nextclade.results.csv", emit: NEXTCLADE_out
    path "all_${caller}_Nextclade.results2.csv", emit: NEXTCLADE_FOR_FHI_out
    path "*.{sh,log}"

    publishDir "${params.outdir}/5_lineage/nextclade", mode: 'link', pattern:"all_${caller}_Nextclade.results.csv"
    publishDir "${params.outdir}/5_lineage/nextclade/log", mode: 'link', pattern:'*.{sh,log}'

    script:
    """
    nextclade --version

    nextclade -j $task.cpus \
        -i $allSamplesConsensus \
        -c all_${caller}_Nextclade.results.csv

    nextalign -j $task.cpus \
        --sequences=$allSamplesConsensus \
        --reference=/home/docker/CommonFiles/reference_nc.fasta \
        --genemap=/home/docker/CommonFiles/genemap.gff \
        --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
        --output-dir=base \
        --output-basename=run

    Rscript /home/docker/Scripts/InsertionAnalysis.R

    nextclade_output_converter.py all_${caller}_Nextclade.results.csv > all_${caller}_Nextclade.results2.csv

    cp .command.sh all.nextclade.sh
    cp .command.log all.nextclade.log 
    """
}

process NSC4FHI_NOISE_NEXTCLADE {
    label 'tiny'

    input:
    path "results/2_bam/noiseextractor/*"
    path "results/nextclade_for_FHI.tsv"
    path "results/report_v${params.pipeline_version}.tsv"

    output:
    path "nextclade_and_noise_for_FHI.tsv", emit: NSC4FHI_NOISE_NEXTCLADE_out
    path "*.{sh,log}"

    publishDir "${params.outdir}/", mode: 'link', pattern:'nextclade_and_noise_for_FHI.tsv'
    publishDir "${params.outdir}/5_lineage/nextclade/log", mode: 'link', pattern:'*.{sh,log}'

    script:
    """
    Rscript /home/docker/Scripts/NSC4FHI_noise_nextclade.R

    cp results/nextclade_and_noise_for_FHI.tsv .

    cp .command.sh all.NSC4FHI_noise_nextclade.sh
    cp .command.log all.NSC4FHI_noise_nextclade.log 
    """
}

// nextflow.enable.dsl=2

params.OUTDIR = "$PWD"
       
Channel.fromPath(params.samplelist)
    .splitCsv(header:true, sep:",")
    .map{ row -> "${row.sample}" }
    .into{input_base_ch1; input_base_ch2; input_base_ch3}

process Pangolin { 
    tag "${base}"

    input:    
      val base from input_base_ch1
      path viralrecon_result from params.viralrecon_folder

    output:
      val base into report_pangolin_ch1
      file("*.csv") into report_pangolin_ch2
      file("*.log")
    
    publishDir "${params.viralrecon_folder}/results/variants/ivar/pangolin/", mode: 'copy', pattern:'*.csv'
    publishDir "${params.viralrecon_folder}/results/variants/ivar/pangolin/log/", mode: 'copy', pattern:'*.log'

    script:
    """
    ln -s ${viralrecon_result}/results/variants/ivar/consensus/${base}.AF0.75.consensus.fa .
    ln -s ${viralrecon_result}/results/variants/ivar/consensus/${base}.AF0.75.consensus.qual.txt .

    pangolin ${base}.AF0.75.consensus.fa -t 4
    mv lineage_report.csv ${base}.ivar.AF0.75.consensus_pangolin_lineage_report.csv 
    cp .command.log ${base}_pangolin.log
    """
}


process Extract_multiqc_viralrecon {
    tag "${base}"

    input:    
      val base from input_base_ch2
      path viralrecon_result from params.viralrecon_folder + "/results/"

    output:
      file("*mqc.tsv") into report_mqc_ch
//    publishDir "${params.viralrecon_folder}/results/variants/ivar/pangolin/", mode: 'copy', pattern:'*.csv'
//    publishDir "${params.viralrecon_folder}/results/variants/ivar/pangolin/log/", mode: 'copy', pattern:'*.log'

    script:
    """
    multiqc_to_custom_tsv.py -md results/multiqc/multiqc_data/ -s "${base}"
    mv summary_variants_metrics_mqc.tsv ${base}_mqc.tsv
    """
}


process Report_copy {
    tag "${base}"

    input:    
      val base from report_pangolin_ch1
      path viralrecon_result from params.viralrecon_folder
      
    output:
      file("*gz*")
      file("*consensus*")
    
    publishDir "${params.viralrecon_folder}/results_report/${base}", mode: 'copy', pattern:'*gz*'
    publishDir "${params.viralrecon_folder}/results_report/${base}", mode: 'copy', pattern:'*consensus*'
    
    script:
    """
    ln -s ${viralrecon_result}/results/variants/ivar/${base}.AF0.75.vcf.gz .
    ln -s ${viralrecon_result}/results/variants/ivar/${base}.AF0.75.vcf.gz.tbi .

    ln -s ${viralrecon_result}/results/variants/ivar/consensus/${base}.AF0.75.consensus.fa .
    ln -s ${viralrecon_result}/results/variants/ivar/consensus/${base}.AF0.75.consensus.qual.txt .

    ln -s ${viralrecon_result}/results/variants/ivar/pangolin/${base}.ivar.AF0.75.consensus_pangolin_lineage_report.csv .
    """

}


process Report {
//    tag "${base}"

    input:
      file samplelist from Channel.fromPath("$params.samplelist")
      file pangolin_report from report_pangolin_ch2.collect()
      file mqc_report from report_mqc_ch.collect()

    output:
      file("*.xls")
    
    publishDir "${params.viralrecon_folder}/", mode: 'copy', pattern:'*.xls'

    script:
    """
    cat *mqc.tsv > mqc.tsv
    cat *report.csv > pangolin.csv
    extract_report.py
    """
}

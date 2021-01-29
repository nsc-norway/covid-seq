# covid-seq at NSC

## Library prep

Two library preparation methods are being tested now:

1. [Swift Normalase® Amplicon SARS-CoV-2 Panels (SNAP) with Additional Genome Coverage](https://swiftbiosci.com/swif-normalase-amplicon-sars-cov-2-panels/).  
2. [EasySeq™ RC-PCR SARS CoV-2 (novel coronavirus) Whole Genome Sequencing](https://www.nimagen.com/covid19).

## Bioinformatics analysis

Refernces used for the analysis can be found in the folder <util>.
  
Primary analysis is based on [nf-core/viralrecon](https://nf-co.re/viralrecon/1.1.0). Refer to folder <viralrecon> copied from [github repo](https://github.com/nf-core/viralrecon).
  
Secondary analysis is based on [Pangolin](https://cov-lineages.org/) using [docker image](https://hub.docker.com/r/staphb/pangolin) and custom scripts.

Execution can be found in <script.sh>

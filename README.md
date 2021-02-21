![NorSeq logo](http://genomics.no/oslo/uploads/images/NorSeqLogo_Acronym%20Colormix.png)

# Covid-seq at [OUS, NSC node](https://www.sequencing.uio.no/) of [NorSeq](https://www.norseq.org/)

## SARS-CoV-2 whole genome sequencing based on multiplexed amplicon method using short-read Illumina sequencers

## Library prep

Two library preparation methods are being tested now:

1. [Swift Normalase® Amplicon SARS-CoV-2 Panels (SNAP) with Additional Genome Coverage](https://swiftbiosci.com/swif-normalase-amplicon-sars-cov-2-panels/).  
2. [EasySeq™ RC-PCR SARS CoV-2 (novel coronavirus) Whole Genome Sequencing](https://www.nimagen.com/covid19).

150 bp paired-end sequences will be generated primarily using [Illumina](https://www.illumina.com) [NovaSeq 6000](https://www.illumina.com/systems/sequencing-platforms/novaseq.html) with [NextSeq](https://www.illumina.com/systems/sequencing-platforms/nextseq.html), [MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html) and [HiSeq](https://www.illumina.com/systems/sequencing-platforms/hiseq-2500.html) as backup in the said order.

## Bioinformatics analysis

_This is a snapshot of the production code_

References used for the analysis can be found in the folder _util_.
  
Primary analysis is losely based on [nf-core/viralrecon](https://nf-co.re/viralrecon/1.1.0) - [github repo](https://github.com/nf-core/viralrecon).  

Secondary analysis is based on [Pangolin](https://cov-lineages.org/) using [docker image](https://hub.docker.com/r/staphb/pangolin) executed using Nextflow pipeline found in folder _pangolin_ and custom scripts.

Execution:
```bash
nextflow run main.nf --outpath <Output_folder> --samplelist <SampleList.csv> -resume

```

_Nextflow + Singularity (through Docker) + SLURM executed in linux cluster with 500+ cores and 2 TB+ Memory_



The purpose of this repository is to provide the scripts that are used to produce the figures 
for the manuscript by Froehlich & Uyar et al, 2020. 

See the Biorxiv Preprint here:

# crispr-DART pipeline results

The targeted sequencing data of CRISPR-Cas9 treated and control samples were processed using the [crispr-DART pipeline](https://github.com/BIMSBbioinfo/crispr_dart). 

## Pipeline output

### Reports
The reports output of the pipeline from this analysis can be browsed [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/froehlich_uyar_et_al_2020/reports/index.html): 

### Alignment and indel files
The input files to run the pipeline such as `settings.yaml`, `sample_sheet.csv`, `cutsites.bed`, and `comparisons.tsv` along 
with the output files and folders can be downloaded from [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/froehlich_uyar_et_al_2020/crispr_dart_pipeline_output.tgz).

# How to reproduce the figures in the manuscript

## Preparing the inputs

The necessary pipeline outputs needs to be downloaded from [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/froehlich_uyar_et_al_2020/crispr_dart_pipeline_output.tgz).
The downloaded compressed folder needs to be uncompressed and the `settings.yaml` file needs to be modified according to the location of the file paths to 
the various input files. In order to reproduce the figures, the only fields that need to be modified in the `settings.yaml` file are:

- sample_sheet: /path/to/sample_sheet.csv
- cutsites: /path/to/cut_sites.bed
- reference_fasta: /path/to/ce11.fa
- output-dir: /path/to/output
- comparisonsFile: /path/to/comparisons.tsv
- Rscript: /path/to/Rscript 
 
## Scripts in this repository

The scripts in this repository take as input the `settings.yaml` file which contains all the necessary links to other important
files such as the sample sheet, alignment outputs, or the location of the genome sequence file. The pipeline output files are 
parsed and processed to make the figures that were further asthetically processed for publication. However, you can find the 
raw versions of the figures printed by these scripts in this repository. 

### How to run the scripts

#### Summary plots

- To get correlation between sgrna efficiencies and external scores calculated for the designed guides
> cd summary_plots
> /usr/bin/Rscript ../scripts/sgRNA_scores.R */path/to/settings.yaml* ../data/sgRNAscores.txt

- To get various summary plots from the processed pipeline output
> cd summary_plots
> /usr/bin/Rscript ../scripts/summary_plots.R */path/to/settings.yaml* 

#### lin-41 RNA analysis

- To cluster and compare lin-41 pacbio RNA reads:
> cd rna_analysis
> /usr/bin/Rscript ../scripts/cluster_pacbio_reads.R */path/to/settings.yaml* lin41_RNApacbio_L1_all lin41_RNApacbio_L4_all cluster_pacbio_all

- To analyse impact of deletions in L1 vs L4 abundance
> cd rna_analysis
> /usr/bin/Rscript ../scripts/analyse_impact_on_rna_expression.R */path/to/settings.yaml*

#### Fitness/generations analysis

- To make the plots about how deletion or reads with deletions get selected agains over generations from F2 to F5. 
> cd generations_analysis
> /usr/bin/Rscript ../scripts/deletions_impact_on_fitness.R */path/to/settings.yaml*  ../data/analysis_table.tsv



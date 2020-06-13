The purpose of this repository is to provide the scripts that are used to produce the figures 
for the manuscript by Froehlich et al, 2020. 

See the Biorxiv Preprint here:

# CRISPR-DART pipeline results

The targeted sequencing data of CRISPR-Cas9 treated and control samples were processed using the [CRIPSR-DART pipeline](https://github.com/BIMSBbioinfo/pigx_crispr). 

The input files to run the pipeline such as `settings.yaml`, `sample_sheet.csv`, `cutsites.bed`, and `comparisons.tsv` along 
with the output files and folders can be downloaded from [here](link to bimsbstatic). 

# Scripts in this repository

The scripts in this repository take as input the `settings.yaml` file which contains all the necessary links to other important
files such as the sample sheet, alignment outputs, or the location of the genome sequence file. The pipeline output files are 
parsed and processed to make the figures that were further asthetically processed for publication. However, you can find the 
raw versions of the figures printed by these scripts in this repository. 

# How to run the scripts

## Summary plots

- To get correlation between sgrna efficiencies and external scores calculated for the designed guides
/usr/bin/Rscript scripts/sgRNA_scores.R <settings.yaml that was used to run the pipeline> ./data/sgRNAscores.txt

- To get various summary plots from the processed pipeline output
/usr/bin/Rscript scripts/summary_plots.R <settings.yaml that was used to run the pipeline>

## lin-41 RNA analysis

- To cluster and compare lin-41 pacbio RNA reads:
> /usr/bin/Rscript scripts/cluster_pacbio_reads.R <path to pipeline settings.yaml> lin41_RNApacbio_L1_all lin41_RNApacbio_L4_all cluster_pacbio_all

- To analyse impact of deletions in L1 vs L4 abundance
> /usr/bin/Rscript scripts/analyse_impact_on_rna_expression.R <path to pipeline settings.yaml>

## Fitness/generations analysis

- To make the plots about how deletion or reads with deletions get selected agains over generations from F2 to F5. 
> /usr/bin/Rscript scripts/deletions_impact_on_fitness.R <path to pipeline settings.yaml>  data/analysis_table.tsv


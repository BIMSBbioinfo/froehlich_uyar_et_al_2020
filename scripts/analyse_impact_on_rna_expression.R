args = commandArgs(trailingOnly=TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]

#initial setup
library(ggplot2)
library(data.table)
library(rtracklayer)
library(yaml)
library(ggpubr)

source('../utility_functions.R')

settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)
pipelineOutputDir <- settings$`output-dir`

settings <- yaml::read_yaml(settings_file)
pipelineOutputDir <- settings$`output-dir`


# define plotting functions 

# define a function to get relative ratios of indels between combinations of L1, L4, and/or DNA samples
# @param dt: data.table which can be a subset of indels object
# @param samples: list of samples with names 'L1', 'L4', and/or 'dna'. 
# The ratios are calculated in the order provided in the list object: 
# e.g. L1, L4, dna will give log2(ratios) betwen L1/L4, L1/dna, and L4/dna
# if the order is L4, dna, L1, then it gives, L4/dna, L4/L1, and L1/dna
# @param readSupportThreshold: indels that are not supported by at least this many reads are removed
get_relative_ratios <- function(dt, samples, readSupportThreshold = 0) {
  dt <- dt[sample %in% unlist(samples)][ReadSupport >= readSupportThreshold]
  #now summarize indels by sample 
  dt <- data.table::dcast(dt, indelID ~ sample, value.var = 'freq')
  
  # remove the ones don't don't occur in all samples
  dt <- dt[apply(dt, 1 , function(x) sum(is.na(x)) == 0),]
  
  ratios <- do.call(cbind, apply(combn(names(samples), 2), 2, function(x) {
    r <- data.table(log2(dt[,get(samples[[x[1]]])] / dt[,get(samples[[x[2]]])]))
    colnames(r) <- paste0(x[1],'_vs_', x[2])
    return(r)
  }))
  
  return(cbind(dt[,'indelID',drop = F], ratios))
}

get_plots <- function(dt, samples_analysed, sample_order = c('L4', 'L1', 'dna'), 
                      readSupportThreshold = 0) {
  plots <- sapply(simplify = F, names(samples_analysed), function(analysis) {
    message(analysis)
    relative_ratios <- get_relative_ratios(dt = dt, 
                                           samples = samples_analysed[[analysis]][sample_order], 
                                           readSupportThreshold = readSupportThreshold)
    
    site_overlaps <- get_overlaps(as(relative_ratios$indelID, 'GRanges'), sites, 
                                  ignore.strand = T, type = 'any')
    # categorize each deletion by combination of sites that the deletion overlaps
    categories <- apply(site_overlaps, 1, function(x)  
      paste(colnames(site_overlaps)[x], collapse = ':'))
    # this hierarchy of setting categories is important 
    relative_ratios$category <- 'other'
    site_count <- apply(site_overlaps, 1, function(x) sum(x))
    relative_ratios[site_count == 0]$category <- 'no_overlap'
    relative_ratios[site_count == 1]$category <- paste0(categories[site_count == 1], '_only')
    # deletions that overlap both seed1 and seed2
    both_seeds <- apply(site_overlaps[,c('LCS1_seed', 'LCS2_seed')], 1, function(x) sum(x) == 2)
    relative_ratios[both_seeds]$category <- 'both_seeds'
    
    plot_names <- apply(combn(sample_order, 2), 2, function(x) {
      paste0(x[1], '_vs_', x[2])
    })
    
    p <-  sapply(simplify = F, plot_names, function(y) {
      df <- relative_ratios[category != 'other']
      ggboxplot(df, x = 'category', y = y,
                palette = "jco",
                add = "jitter", 
                nrow = 1,  
                add.params = list(size = 1, alpha = 0.5), outlier.shape = NA) +
        labs(x = 'Whether the deletion overlaps the site') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    return(p)
  })
}

#define functionally important sites in the targeted region  
sites <- list('LCS1_seed' = 'chrI:9335255-9335263', 
              'LCS2_seed' = 'chrI:9335208-9335214', 
              'LCS1_3compl' = 'chrI:9335264-9335276',
              'LCS2_3compl' = 'chrI:9335215-9335227')
sites <- as(unlist(sites), 'GRanges')

# define samples to be analysed (trio of dna, l1, and l4 )
samples_analysed = list('lin41_pacbio_sg1516_2627_pool3' = list('dna' = 'lin41_DNA_sg1516_2627_pool3', 
                                               'L1' = 'lin41_RNApacbio_L1_all',
                                               'L4' = 'lin41_RNApacbio_L4_all'))

# get sample sheet, subset by analysed samples 
sampleSheet <- data.table::fread(settings$sample_sheet)
sampleSheet <- sampleSheet[sample_name %in% unlist(samples_analysed)]

genome <- Biostrings::readDNAStringSet(settings$reference_fasta, format = 'fasta')

target <- as(sampleSheet$target_region[1], 'GRanges')

# read indels
indels <- do.call(rbind, getIndels(settings$`output-dir`, sampleSheet$sample_name))
indels$freq <- indels$ReadSupport/indels$coverage
indels$indelID <- paste0(indels$seqname, ':', indels$start, '-', indels$end)

# remove indels that are not on target region
indels <- indels[unique(queryHits(findOverlaps(as(indels, 'GRanges'), target, ignore.strand = T)))]
indels$width <- abs(indels$start - indels$end + 1)
# focus on deletions only 
deletions <- indels[indelType == 'D']

# only consider rna samples, filter for read support > 5
plots <- get_plots(deletions, samples_analysed, sample_order = c('L4', 'L1'), readSupportThreshold = 5)

# save plots
lapply(names(plots), function(analysis) {
  lapply(names(plots[[analysis]]), function(x) {
    outfile <- paste(analysis, x, 'pdf', sep = '.')
    p <- plots[[analysis]][[x]] + labs(title = analysis, subtitle = x)
    ggsave(outfile, p, width = 8, height = 8, units = 'in')
  }) 
})







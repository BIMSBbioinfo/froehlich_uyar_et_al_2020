# script to make figures about the impact of deletions on the fitness of the worms

library(ggplot2)
library(data.table)
library(rtracklayer)
library(yaml)
library(ggpubr)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

settings_file <- args[1] # settings.yaml file that was used to run the pipeline 
analysis_table_file <- args[2] # table of analysed samples (analysis_table.tsv)

source('../utility_functions.R')
analysisTable <- data.table::fread(analysis_table_file, header = T)

# get sample analysis groups: F2 .. F5 generation samples  
samples_analysed <- lapply(split(analysisTable, analysisTable$analysisGroup), 
                           function(x) {
                             lapply(split(x, x$generation), function(y) y$sample)
                           })

#initial setup
settings <- yaml::read_yaml(settings_file)
pipelineOutputDir <- settings$`output-dir`

#define functionally important sites in the targeted region  
sites <- list('LCS1_seed' = 'chrI:9335255-9335263', 
              'LCS2_seed' = 'chrI:9335208-9335214', 
              'LCS1_3compl' = 'chrI:9335264-9335276',
              'LCS2_3compl' = 'chrI:9335215-9335227',
              'polyA' = 'chrI:9334816-9334821',
              'stop_codon' = 'chrI:9335965-9335967')
sites <- as(unlist(sites), 'GRanges')

# get sample sheet, subset by analysed samples 
sampleSheet <- data.table::fread(settings$sample_sheet)
sampleSheet <- sampleSheet[sample_name %in% unlist(samples_analysed)]

genome <- Biostrings::readDNAStringSet(settings$reference_fasta, format = 'fasta')

target <- as(sampleSheet$target_region[1], 'GRanges')

# BEGIN analysis
indels <- do.call(rbind, getIndels(settings$`output-dir`, sampleSheet$sample_name))
indels$freq <- indels$ReadSupport/indels$coverage
indels$indelID <- paste0(indels$seqname, ':', indels$start, '-', indels$end)

# remove indels that are not on target region
indels <- indels[unique(queryHits(findOverlaps(as(indels, 'GRanges'), target, ignore.strand = T)))]
indels$width <- abs(indels$start - indels$end + 1)

# focus on deletions only 
deletions <- indels[indelType == 'D']

# also import reads with deletions 
reads_with_deletions <- do.call(rbind, lapply(unique(sampleSheet$sample_name), function(s) {
  dt <- get_reads_with_indels(settings, s, subsetByRegion = target)[indelType == 'D']
  dt$sample <- s
  return(dt)
}))

message('importing alignments to get lib sizes')
# import alignment files and get library sizes for each sample 
cl <- parallel::makeForkCluster(length(unique(sampleSheet$sample_name)))
lib_sizes <- pbapply::pbsapply(cl = cl, unique(sampleSheet$sample_name), function(s) {
  a <- get_aln(settings, s, subsetByRegion = target)
  length(a)
})
parallel::stopCluster(cl)

# given a reference generation, find out which deletions survive to the offspring generations
# dt: data.table can be a subset of indels/deletions 
# samples: list of samples with generation ids used as names (e.g. list('F2' = a, 'F3' = b))
# reference_sample : which generation is the reference? e.g. 'F2'
# readSupportThreshold: only keep deletions in reference sample with at least this many reads supporting
get_survival_rate <- function(dt, samples, reference_sample = NULL, readSupportThreshold) {
  dt <- dt[sample %in% unlist(samples)][ReadSupport >= readSupportThreshold]
  #now summarize indels by sample 
  dt <- data.table::dcast(dt, indelID ~ sample, value.var = 'freq')
  # remove deletions that don't occur in reference sample 
  if(!is.null(reference_sample)) {
    dt <- dt[!is.na(get(samples[[reference_sample]]))]
  }
  colnames(dt)[2:ncol(dt)] <- names(samples)
  return(dt)
}

# reads_dt: data.table reads withe deletion coordinates
# samples: list where the names of the list represent ordered generations
plot_read_survival <- function(reads_dt, samples, sites, lib_sizes) {
  dt <- reads_dt[sample %in% unlist(samples)]
  site_overlaps <- get_overlaps(as(dt, 'GRanges'), sites)
  
  dt$category <- 'other'
  site_count <- apply(site_overlaps, 1, function(x) sum(x))
  dt[site_count == 0]$category <- 'no_overlap'
  dt[site_count == 1]$category <- paste0(colnames(site_overlaps)[apply(site_overlaps[site_count==1], 1, which)], 
                                         '_only')
  # deletions that overlap both seed1 and seed2
  both_seeds <- apply(site_overlaps[,c('LCS1_seed', 'LCS2_seed')], 1, function(x) sum(x) == 2)
  dt[both_seeds]$category <- 'both_seeds'
  # assign generation to each sample
  dt$generation <- names(samples)[match(dt$sample, unlist(samples))]
  # normalize by library size
  dt$lib_size <- lib_sizes[dt$sample]
  dt <- dt[,length(unique(name)),by = c('category', 'generation', 'lib_size')]
  dt$percent_reads <- dt$V1/dt$lib_size * 100
  p <- ggplot(dt[category != 'other'], aes(x = generation, percent_reads)) + 
    geom_bar(stat = 'identity', aes(fill = generation), position = 'dodge') +
    facet_wrap(~ category, scales = 'free') + 
    labs(y = '% reads with deletions affecting sites normalized by library size') +
    scale_fill_brewer(type = 'qual', palette = 2)
  
  return(p)
}

get_plots <- function(dt, samples_analysed, sample_order = c('F1', 'F2', 'F3', 'F4', 'F5'), 
                      reference_generation = 'F1', 
                      readSupportThreshold = 0, sites) {
  plots <- sapply(simplify = F, names(samples_analysed), function(analysis) {
    message(analysis)
    
    samples <- samples_analysed[[analysis]][sample_order]
    deletion_survival <- get_survival_rate(dt = dt, 
                                           samples = samples, 
                                           reference_sample = reference_generation, 
                                           readSupportThreshold = readSupportThreshold)
    site_overlaps <- get_overlaps(as(deletion_survival$indelID, 'GRanges'),
                                  sites, ignore.strand = T, type = 'any')
    
    # make a heatmap
    M <- as.matrix(data.frame(deletion_survival[,-1], 
                              row.names = deletion_survival$indelID))
    M[is.na(M)] <- 0
    p1 <- pheatmap::pheatmap(M, cluster_cols = F, scale = 'row',
                       show_rownames = F, cellwidth = 50,
                       annotation_row = data.frame(site_overlaps * 1, 
                                                   row.names = rownames(M)), 
                       annotation_legend = FALSE,
                       main = paste('Deletion frequency over generations\n', analysis))
    
    # convert table to boolean (NA means eliminated / not survived)
    deletion_survival <- cbind(deletion_survival[,1], 
                               apply(deletion_survival[,-1], 2, is.na))
    
    # categorize each deletion by combination of sites that the deletion overlaps
    categories <- apply(site_overlaps, 1, function(x)  
      paste(colnames(site_overlaps)[x], collapse = ':'))
    # this hierarchy of setting categories is important 
    # find deletions that overlap at most one sites
    deletion_survival$category <- 'other'
    site_count <- apply(site_overlaps, 1, function(x) sum(x))
    deletion_survival[site_count == 0]$category <- 'no_overlap'
    deletion_survival[site_count == 1]$category <- paste0(categories[site_count == 1], '_only')
    # deletions that overlap both seed1 and seed2
    both_seeds <- apply(site_overlaps[,c('LCS1_seed', 'LCS2_seed')], 1, function(x) sum(x) == 2)
    deletion_survival[both_seeds]$category <- 'both_seeds'
    
    df <- melt(deletion_survival[category != 'other'], id.vars = 'category', measure.vars = names(samples))
    colnames(df) <- c('site_overlap', 'generation', 'eliminated')
    p2 <- ggplot(df, aes(x = generation)) + 
      geom_bar(aes(fill = eliminated), position = 'fill') + 
      facet_grid(~ site_overlap) +
      scale_fill_brewer(type = 'qual', palette = 2)
    return(list('p1' = p1, 'p2' = p2))
  })
}

# Check how many deletions in F2 are eliminated in later generations
# Subset deletions for those that exist in F2. Check the number of deletions that still exist in F3, F4, and F5 - 
# categorized by the deletions' overlap with important sites. 
# 'no_overlap' serves as a control: This is the category of deletions thatdon't overlap any of the other sites. 
plots <- get_plots(dt = deletions, samples_analysed = samples_analysed, reference_generation = 'F2', 
                   sample_order = c('F2', 'F3', 'F4', 'F5'), readSupportThreshold = 0, sites = sites)

lapply(names(plots), function(analysis) {
  ggsave(filename = paste0("deletions_over_generations.start_at_F2.",analysis,'.pdf'), 
         plots[[analysis]][['p2']], width = 12, height = 8, units = 'in')
  pdf(file = paste0("deletions_over_generations.heatmap.start_at_F2.",analysis,'.pdf'))
  print(plots[[analysis]][['p1']])
  dev.off()
})

# same as above, starting at F1
plots <- get_plots(dt = deletions, samples_analysed = samples_analysed, reference_generation = 'F1', 
                   sample_order = c('F1', 'F2', 'F3', 'F4', 'F5'), readSupportThreshold = 0, sites = sites)

lapply(names(plots), function(analysis) {
  ggsave(filename = paste0("deletions_over_generations.start_at_F1.",analysis,'.pdf'), 
         plots[[analysis]][['p2']], width = 12, height = 8, units = 'in')
  pdf(file = paste0("deletions_over_generations.heatmap.start_at_F1.",analysis,'.pdf'))
  print(plots[[analysis]][['p1']])
  dev.off()
})

# Check how many reads with deletions in F1 are eliminated in later generations {.tabset}
# This time, we count reads (rather than grouping reads as deletions/genetypes) that affect different sites (via deletions in reads) and check the percentage/count of those that exist in different generations. 
#'no_overlap' serves as a control: This is the category of reads that have deletions that don't overlap any of the other sites. 
plots <- sapply(simplify = F, names(samples_analysed), function(analysis) {
  message(analysis)
  plot_read_survival(reads_with_deletions, samples_analysed[[analysis]], sites, lib_sizes)
})

lapply(names(plots), function(analysis) {
  ggsave(filename = paste0("reads_with_deletions_over_generations.",analysis,'.pdf'), 
         plots[[analysis]] + labs(title = analysis), width = 12, height = 8, units = 'in')
})






suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(yaml)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggrastr)))

args <- commandArgs(trailingOnly = TRUE)

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]
L1_sample <- args[2] # L1 stage pacbio lin-41 RNA sample
L4_sample <- args[3] # L4 stage pacbio lin-41 RNA sample
prefix <- args[4] # prefix to be used for output files/folders

source('../utility_functions.R')

samples = list('L1' = L1_sample, 
               'L4' = L4_sample)

settings <- yaml::read_yaml(settings_file)
sampleSheet <- data.table::fread(settings$sample_sheet)
pipelineOutputDir <- settings$`output-dir`

sampleSheet <- sampleSheet[sample_name %in% unlist(samples)]

#define let-7 binding site coordinates on the genome  
let7_sites <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame('seqname' = 'chrI', 
                  'start' = c(9335255, 9335208), 
                  'end' = c(9335262,9335214), 'strand' = '-', 
                  row.names = c('let7_site1', 'let7_site2'))
)

genome <- Biostrings::readDNAStringSet(settings$reference_fasta, format = 'fasta')
target <- as(sampleSheet$target_region[1], 'GRanges')

# take the reads that completely cover 
# the highly covered region from around start of target region up to the first intron 
region <- GenomicRanges::GRanges(seqnames = 'chrI', ranges = IRanges(start = 9334840,
                                                                     end = 9336100))
# get total reads aligned 
aln_L1 <- get_aln(settings, samples[['L1']], subsetByRegion = target)
reads_with_indels_L1 <- get_reads_with_indels(settings, samples[['L1']])

aln_L4 <- get_aln(settings, samples[['L4']], subsetByRegion = target)
reads_with_indels_L4 <- get_reads_with_indels(settings, samples[['L4']])

# subset alignments to exclude those that don't cover the region of interest
aln_L1 <- aln_L1[unique(subjectHits(findOverlaps(region, 
                                                 aln_L1, 
                                                 type = 'within', ignore.strand = TRUE)))]
aln_L4 <- aln_L4[unique(subjectHits(findOverlaps(region, 
                                                 aln_L4, 
                                                 type = 'within', ignore.strand = TRUE)))]

reads_with_indels_L1 <- reads_with_indels_L1[name %in% mcols(aln_L1)[['qname']]]
reads_with_indels_L4 <- reads_with_indels_L4[name %in% mcols(aln_L4)[['qname']]]

# Get k-mer profiles of reads and cluster reads, annotate clusters
x1 <- sequenceLayer(x = mcols(aln_L1)$seq, cigar = cigar(aln_L1), 
                    from = 'query', to = 'pairwise')
x4 <- sequenceLayer(x = mcols(aln_L4)$seq, cigar = cigar(aln_L4), 
                    from = 'query', to = 'pairwise')

let7_seqs <- extractAt(genome$chrI, ranges(let7_sites))

# represent read alignments by kmer count matrices and cluster 
kmers <- generateKmers(5)
M <- countPattern(seqs = c(x1, x4), patterns = kmers, maxMismatch = 1, nCores = 50)
rownames(M) <- c(mcols(aln_L1)$qname, mcols(aln_L4)$qname)

if(!file.exists(paste0(prefix, ".read_clusters.seurat.RDS"))) {
  sc <- cluster_reads(t(M))
  saveRDS(sc, file = paste0(prefix, ".read_clusters.seurat.RDS"))
} else {
  sc <- readRDS(paste0(prefix, ".read_clusters.seurat.RDS"))
}

sc$let7_site1_deleted <- Biostrings::vcountPattern(let7_seqs$let7_site1, 
                                                   subject = c(x1, x4), max.mismatch = 0) == 0
sc$let7_site2_deleted <- Biostrings::vcountPattern(let7_seqs$let7_site2, 
                                                   subject = c(x1, x4), max.mismatch = 0) == 0
sc$both_deleted <- sc$let7_site1_deleted & sc$let7_site2_deleted
sc$neither_deleted <- !sc$let7_site1_deleted & !sc$let7_site2_deleted

sc$stage <- c(rep('L1', length(aln_L1)), 
              rep('L4', length(aln_L4)))


# Visualize clustered reads. For each cluster, show deletion profile of reads and compare L1 vs L4 stages. 

get_profile <- function(reads_with_indels, clusters, region) {
  df <- do.call(rbind, sapply(simplify = F, unique(clusters$cl), function(x) {
    reads <- clusters[cl == x]$name
    dt <- reads_with_indels[name %in% reads]
    df <- data.table('pos' = start(region):end(region), 
                     'ins' = as.numeric(coverage(as(dt[indelType == 'I'], 
                                                    'GRanges'))[[1]])[start(region):end(region)],
                     'del' = as.numeric(coverage(as(dt[indelType == 'D'], 
                                                    'GRanges'))[[1]])[start(region):end(region)])
    mdf <- melt(df, id.vars = 'pos')
    mdf$cluster <- x
    return(mdf)  
  }))
  return(df)
}

umapPlot_Raster <- function(seu, group_by = 'seurat_clusters', split_by = NULL) {
  require(ggrastr)
  dt <- cbind(data.table(seu@reductions$umap@cell.embeddings), 
              data.table(seu@meta.data))
  p <- ggplot(dt, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point_rast(aes_string(color = group_by), raster.dpi = 300, alpha = 0.1, 
                    raster.width = 10, raster.height = 8) 
  if(!is.null(split_by)) {
    p <- p + facet_grid(as.formula(paste0('~ ', split_by)))  
  }
  p <- p  + theme_bw()
  return(p)
} 

pdf(file = paste0(prefix, '.umap_plots.pdf'), width = 10, height = 8)
print(DimPlot(sc, split.by = 'stage', label = T, label.size = 5) + 
  labs(title = 'Comparison of read clusters by stage'))
print(DimPlot(sc, group.by = 'let7_site1_deleted' ) + 
  labs(title = 'let7_site1_deleted'))
print(DimPlot(sc, group.by = 'let7_site2_deleted' ) +
  labs(title = 'let7_site2_deleted'))
print(DimPlot(sc, group.by = 'both_deleted' ) +
  labs(title = 'both_sites_deleted'))
dev.off()

# also save rastered dimplots
ggsave(paste0(prefix,'.umap.compare_stages.pdf'), 
       umapPlot_Raster(sc, group_by = 'seurat_clusters', split_by = 'stage') + 
         labs(title = 'Comparison of read clusters by stage'), width = 10, height = 8)
ggsave(paste0(prefix,'.umap.lcs1_deleted.pdf'), 
       umapPlot_Raster(sc, group_by = 'let7_site1_deleted') + 
         labs(title = 'let7_site1_deleted'), width = 10, height = 8)
ggsave(paste0(prefix,'.umap.lcs2_deleted.pdf'), 
       umapPlot_Raster(sc, group_by = 'let7_site2_deleted') + 
         labs(title = 'let7_site2_deleted'), width = 10, height = 8)
ggsave(paste0(prefix,'.umap.both_deleted.pdf'), 
       umapPlot_Raster(sc, group_by = 'both_deleted') + 
         labs(title = 'both_sites_deleted'), width = 10, height = 8)

# Visualize deletion profiles of each cluster
reads_with_indels <- rbind(reads_with_indels_L1, reads_with_indels_L4)
clusters <- data.table('name' = rownames(M), 'cl' = sc$seurat_clusters)
clusters$stage <- c(rep('L1', length(aln_L1)), 
                    rep('L4', length(aln_L4)))

df1 <- get_profile(reads_with_indels, clusters[name %in% mcols(aln_L1)$qname], region)
df2 <- get_profile(reads_with_indels, clusters[name %in% mcols(aln_L4)$qname], region)

df1$value.norm <- df1$value / length(unique(reads_with_indels_L1$name))
df2$value.norm <- df2$value / length(unique(reads_with_indels_L4$name))
df1$stage <- 'L1'
df2$stage <- 'L4'

df <- rbind(df1, df2)

p <- ggplot(df[variable == 'del'], aes(x = pos, y = value)) + 
  geom_point(aes(color = stage)) + 
  geom_vline(xintercept = c(min(start(let7_sites)), max(end(let7_sites)))) + 
  facet_wrap(~ cluster, ncol = 4, scales = 'free') + 
  labs(title = 'Comparison of deletion profiles by cluster and stage',
       y = 'Number of reads', x = 'genomic coordinate')

ggsave(filename = paste0(prefix, '.cluster_deletion_profiles.pdf'), p, 
       width = 12, height = 8, units = 'in')

# print bam files for selected clusters 
bamDir <- paste0(prefix, '_bam')
if(!dir.exists(bamDir)) {
  dir.create(bamDir)
}
dummy <- lapply(unique(clusters$cl), function(x) {
  subset_bam_by_reads(bamFile = get_bamfile_path(settings, samples[['L1']]), 
                      read_names = clusters[stage == 'L1'][cl == x]$name, 
                      outfile = file.path(bamDir, paste0('L1.cl', x, '.bam')))
  subset_bam_by_reads(bamFile = get_bamfile_path(settings, samples[['L4']]), 
                      read_names = clusters[stage == 'L4'][cl == x]$name, 
                      outfile = file.path(bamDir, paste0('L4.cl', x, '.bam')))
})





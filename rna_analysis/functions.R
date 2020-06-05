# some necessary functions
# Define some necessary functions
getIndels <- function(pipeline_output_dir, samples) {
  indels <- sapply(simplify = FALSE, samples, function(s) {
    f <- file.path(pipeline_output_dir, 
                   'indels',
                   s, 
                   paste0(s, ".indels.tsv"))
    if(file.exists(f)) {
      dt <- data.table::fread(f)
      dt$sample <- s
      return(dt)
    } else {
      stop("Can't open indels.tsv file for sample",s,
           "at",f,"\n")
    }})
  return(indels)
}

# given a list of potentially overlapping segments sorted by a certain criteria (Best on topo), 
# pick best ones, remove those that overlap the picked ones 
pick_best_segments <- function(gr) {
  gr.kept <- gr
  kept <- c()
  removed <- c()
  while(length(gr) > 0) {
    #message(length(dt.gr))    
    kept <- c(kept, names(gr)[1])
    removed <- c(removed, names(subsetByOverlaps(gr[-1], gr[1])))
    gr <- gr[!names(gr) %in% c(kept, removed)]
  }
  gr.kept <- gr.kept[names(gr.kept) %in% kept]
  return(gr.kept)
}

# gr: GRanges object with 'score' column
get_segments <- function(gr, window_size = 10) {
  dt <- do.call(rbind, lapply(1:(length(gr)-window_size), function(i) {
    #message(i)
    data.table('seqnames' = as.character(seqnames(gr))[1], 'start' = start(gr)[i], 
               'end' = end(gr)[i] + window_size - 1, 
               'window.median' = median(gr[i:(i+window_size-1),]$score), 
               'window.max' = max(gr[i:(i+window_size-1),]$score),
               'name' = paste0('w_', i))
  }))
  # sort windows by mean score 
  dt <- dt[order(window.median, decreasing = T)]
  # keep windows with high scores, remove those that overlap it
  dt.gr <- as(dt, 'GRanges')
  kept <- c()
  removed <- c()
  while(length(dt.gr) > 0) {
    #message(length(dt.gr))    
    kept <- c(kept, dt.gr[1]$name)
    removed <- c(removed, subsetByOverlaps(dt.gr[-1], dt.gr[1])$name)
    dt.gr <- dt.gr[!dt.gr$name %in% c(kept, removed)]
    #message(length(dt.gr))
  }
  dt <- dt[name %in% kept]
  dt.gr <- as(dt, 'GRanges')
  return(dt.gr)
}

# windows: GRanges object 
# scores: GRanges object with scores per base 
# return: average score per window 
get_window_scores <- function(windows, scores) {
  ov <- findOverlaps(windows, scores)
  ovl <- split(as.data.table(ov), queryHits(ov))
  
  s <- as.numeric(pbapply::pbsapply(ovl, function(x) {
    mean(scores[x$subjectHits]$score, na.rm = T)
  }))
  return(s)
}

get_bamfile_path <- function(settings, sample) {
  file.path(settings$`output-dir`, 
            'aln_merged', 
            paste0(sample, '.bam'))
}

get_aln <- function(settings, sample, subsetByRegion = NULL) {
  aln <- GenomicAlignments::readGAlignments(file.path(settings$`output-dir`, 
                                               'aln_merged', 
                                               paste0(sample, '.bam')), 
                                     param = ScanBamParam(what=c("qname", "seq")))
  if(!is.null(subsetByRegion)) {
    aln <- subsetByOverlaps(aln, subsetByRegion, ignore.strand = T, type = 'any')
  }
  return(aln)
}

get_reads_with_indels <- function(settings, sample, subsetByRegion = NULL) {
  dt <- data.table::fread(file.path(settings$`output-dir`, 
                              'indels', 
                              sample, 
                              paste0(sample, 
                                     '.reads_with_indels.tsv')))
  if(!is.null(subsetByRegion)){
    dt <- dt[overlapsAny(as(dt, 'GRanges'), subsetByRegion, ignore.strand = T, type = 'any')]
  }
  return(dt)
}


countPattern <- function(seqs, patterns, maxMismatch = 0, nCores = 1) {
  
  cl <- parallel::makeForkCluster(nCores)
  M <- do.call(rbind, pbapply::pblapply(cl = cl, X = patterns, 
                                        FUN = function(x) {
                                          Biostrings::vcountPattern(x, seqs, 
                                                                    max.mismatch = maxMismatch)
                                        }))
  parallel::stopCluster(cl)
  rownames(M) <- patterns
  return(t(M))
}

generateKmers <- function(k, letters = c("A", "C", "G", "T")) {
  kmer <- c()
  for(i in 1:k){
    kmer <- unlist(lapply(letters, function(x){paste(kmer, x, sep="")}))
  }
  return(kmer)
}

# given a matrix of k-mer counts for reads,
# apply spectral clustering - use seurat workflow
# rows: kmers
# columns: reads 
cluster_reads <- function(M) {
  require(Seurat)
  sc <- Seurat::CreateSeuratObject(counts = M)
  sc <- FindVariableFeatures(sc)
  all.genes <- rownames(sc)
  sc <- ScaleData(sc, features = all.genes)
  sc <- RunPCA(sc, features = VariableFeatures(object = sc))
  sc <- RunUMAP(sc, dims = 1:10)
  sc <- FindNeighbors(sc, dims = 1:10)
  sc <- FindClusters(sc, resolution = 0.05)  
  return(sc)
}


subset_bam_by_reads <- function(bamFile, read_names, outfile) {
  filter_reads <- FilterRules(list(subset_reads=function(x) x$qname %in% read_names))
  indexBam(filterBam(file = bamFile,
                     destination = BamFile(outfile),
                     filter=filter_reads,
                     param=ScanBamParam(what="qname")))
}




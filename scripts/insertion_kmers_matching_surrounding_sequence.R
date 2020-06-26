# this script takes the "insertedSequences.tsv" from the pipeline output
# and matches kmers from the insertions to the surrounding sequence.
library(data.table)
library(reshape2)
library(ggplot2)
library(Biostrings)
library(BSgenome.Celegans.UCSC.ce11)
ce11 <- BSgenome.Celegans.UCSC.ce11

args = commandArgs(trailingOnly=TRUE)
settings_file <- args[1] # settings.yaml file that was used to run the pipeline 
sample_id <- args[2:length(args)] # sample(s) to be analyzed

# initial setup
settings <- yaml::read_yaml(settings_file)
pipelineOutputDir <- settings$`output-dir`

# define kmer length to be used
kmer.length <- 5
cut.site.window <- 100

# load cut coordinates and input sequences
files <- paste0(pipelineOutputDir, "/indels/", sample_id, "/", sample_id, ".insertedSequences.tsv")

ins.ls <- lapply(files, fread)
names(ins.ls) <- sample_id
ins.ls <- lapply(names(ins.ls), function(x) ins.ls[[x]][, sample := x])

ins.dt = rbindlist(ins.ls)
ins.dt <- ins.dt[insertionWidth >= kmer.length]
ins.dt[, idx := 1:nrow(ins.dt)]


# do the controls (rev comp, shuffle, shuffle rev comp)
library(Biostrings)
ins.dt[, insertedSequenceRevComp := as.character(
                                      reverseComplement(
                                        DNAStringSet(insertedSequence)))]
ins.dt[, insertedSequenceShuffled := sapply(
                                       strsplit(insertedSequence, split=''),
                                         function(x)
                                           paste(sample(x), collapse = ""))]
ins.dt[, insertedSequenceShuffledRevComp := as.character(
                                              reverseComplement(
                                                DNAStringSet(insertedSequenceShuffled)))]


# get windows around cut sites,
# assume everything is on the plus strand
site.gr <- GRanges(seqnames = ins.dt$seqname,
                   ranges   = IRanges(start = ins.dt$start,
                                      end   = ins.dt$start),
                   strand   = "+")
site.gr <- resize(site.gr,  cut.site.window / 2 + 1, "start")
site.gr <- resize(site.gr,  cut.site.window + 1,   "end")

# get site and insert sequences
site.seq <- getSeq(ce11, site.gr)

# split string into kmers of length k
canHazKmers <- function(string, k) {
  substring(string,
            seq(1, nchar(string) - k + 1, 1),
            seq(k, nchar(string)        , 1))
}


# split cut site windows into kmers
site.kmer.ls <- lapply(site.seq,
                       function(x)
                         canHazKmers(as.character(x), k = kmer.length)
)
site.kmer.dt <-
  data.table(idx  = rep(1:length(site.seq), lengths(site.kmer.ls)),
             kmer = unlist(site.kmer.ls),
             pos  = rep(1:length(site.kmer.ls[[1]]), length(site.seq)))


# split inserted sequences into kmers
kmer.ls <- lapply(c("insertedSequence",
                    "insertedSequenceRevComp",
                    "insertedSequenceShuffled",
                    "insertedSequenceShuffledRevComp"),
                  function(x) {
                    tmp.ls <- lapply(ins.dt[[x]],
                                     canHazKmers,
                                     k = kmer.length)
                    data.table(idx = rep(1:length(tmp.ls),
                                         lengths(tmp.ls)),
                               kmer = unlist(tmp.ls))
                    })

# merge cut sites and inserted sequences on indices
# (to get the right cut site => insert xref) and
# kmers
# as a result, position corresponds to the starting 
# position of the matching kmer
kmer.ins.mrg.ls <-  lapply(kmer.ls, function(x) merge(site.kmer.dt, x, by = c("idx", "kmer")))
mrg.summ.ls <- lapply(kmer.ins.mrg.ls, function(x) x[,.(cnt = .N), by = pos])
names(mrg.summ.ls) <- c("insertedSequence",
                        "insertedSequenceRevComp",
                        "insertedSequenceShuffled",
                        "insertedSequenceShuffledRevComp")
lapply(names(mrg.summ.ls), function(x) setnames(mrg.summ.ls[[x]], "cnt", x))

kmer.counts.dt <- Reduce(function(...)
                    merge(...,
                          by = "pos",
                          all = TRUE),
                    mrg.summ.ls
)


kmer.counts.m <- melt(kmer.counts.dt, id.vars = "pos")
kmer.counts.m[, revcomp := grepl("RevComp", variable)]
kmer.counts.m[revcomp == TRUE, value := -value]
kmer.counts.m[, variable := sub("RevComp", "", variable)]


p <- ggplot(NULL , aes(x = pos-(cut.site.window / 2), y = value)) +
      geom_bar(data = kmer.counts.m[variable == "insertedSequence"],
               aes(fill = "insertedSequence"),
              stat = "identity")  +
      geom_bar(data = kmer.counts.m[variable == "insertedSequenceShuffled"],
              aes(fill = "insertedSequenceShuffled"),
              stat = "identity")  +
      scale_fill_manual(values = c("black", "grey50"),
                        labels = c("inserted sequence", "shuffled ctrl")) +
      scale_y_continuous(limits = c(-max(kmer.counts.m$value),
                                    max(kmer.counts.m$value))) +
      xlab("distance form insertion site (bp)") +
      ylab("5-mer matches (count)") +
      geom_hline(yintercept = 0) +
      theme(aspect.ratio = 2,
            legend.title = element_blank())
      

# Print plots to pdf
if (length(sample_id) == 1) {
  out.file <- paste0("insertions_matched_to_surrounding_sequence_5bp_", sample_id, ".pdf")
} else {
  out.file <- paste0("insertions_matched_to_surrounding_sequence_5bp_", length(sample_id), "_samples.pdf")
}
ggsave(out.file, p)
# make plots correlating sgRNA efficiencies with other metrics for sgRNA 

library(ggplot2)
library(ggrepel)
library(data.table)
library(ggpubr)
library(rtracklayer)
library(yaml)
library(ComplexHeatmap)
library(circlize)

args = commandArgs(trailingOnly=TRUE)

ggplot2::theme_set(theme_pubclean())

# settings.yaml file that was used to run the pipeline 
settings_file <- args[1]
sgrna_score_file <- args[2]
settings <- yaml::read_yaml(settings_file)

sampleSheet <- data.table::fread(settings$sample_sheet)

# get a mapping between guides and samples
sampleGuides <- sapply(simplify = F, unique(sampleSheet$sample_name), function(s) {
  sgRNAs <- unlist(strsplit(x = sampleSheet[sampleSheet$sample_name == s,]$sgRNA_ids, 
                            split = ':'))
  return(sgRNAs)
})

# import sgRNA efficiencies from all processed samples
sgrna_stats <- do.call(rbind, lapply(names(sampleGuides), function(s) {
  g <- sampleGuides[[s]]
  if(g[1] != 'none') {
    dt <- data.table::fread(file.path(settings$`output-dir`, 
                                      'indels', 
                                      s, 
                                      paste0(s, '.sgRNA_efficiency.tsv')))
    # only keep scores for the sgRNAs used in the sample
    dt <- dt[sgRNA %in% g]
    return(dt)
  }
}))

# most sgrnas have multiple efficiency scores because they were used in multiple samples.
# we take the average in such cases. 
# which(table(scores$sgRNA) > 1) 
sgrna_stats <- sgrna_stats[,mean(scores),by = sgRNA]
colnames(sgrna_stats)[2] <- 'sgrna_efficiency'


# import sgrna scores collated from various resources 
# these scores are about the sgrna sequence content and design
# keep the relevant columns 
sgrna_scores_boolean <- data.table::fread(sgrna_score_file)[,c(1, 16:27)] 
sgrna_scores_boolean <- cbind(sgrna_scores_boolean[,1], 
                              apply(sgrna_scores_boolean[,-1], 2, as.logical))
colnames(sgrna_scores_boolean)[1] <- 'sgRNA'
colnames(sgrna_scores_boolean) <- gsub("[\\' ]", "_", colnames(sgrna_scores_boolean))

                                

sgrna_scores_numeric <- data.table::fread(sgrna_score_file)[,c(1, 28:38)] 
sgrna_scores_numeric <- cbind(sgrna_scores_numeric[,1], 
                              apply(sgrna_scores_numeric[,-1], 2, as.numeric))
colnames(sgrna_scores_numeric)[1] <- 'sgRNA'
colnames(sgrna_scores_numeric) <- gsub("[\\' ]", "", colnames(sgrna_scores_numeric))

# now merge the sgRNA efficiency scores we obtained from the pipeline and the metrics obtained from external sources

# type: choose 'boolean' or 'numeric'
plot_scores <- function(dt, type) {
  require(ggpubr)
  if(type == 'boolean') {
    mdt <- melt(dt, id.vars = c('sgRNA', 'sgrna_efficiency'), 
                measure.vars = setdiff(colnames(dt), c('sgRNA', 'sgrna_efficiency')))
    plots <- sapply(simplify = F, split(mdt, mdt$variable), function(x){
      if(sum(x$value == TRUE) > -1 & sum(x$value == FALSE) > -1) {
        ggboxplot(x, x = 'value', y = 'sgrna_efficiency',
                  palette = "jco", legend = 'none',
                  add = "jitter", color = 'value',
                  add.params = list(size = 1, alpha = 0.5), outlier.shape = NA) + 
          stat_compare_means() + labs(title = unique(x$variable))
      } else {
        message(paste(c("Not plotting for ",unique(as.character(x$variable)), 
        " the variable doesn't contain two categories with at least 5 members")))
        return(NULL)
      }
    })
    plots <- plots[!sapply(plots, is.null)]
    return(plots)
  } else if(type == 'numeric') {
    mdt <- melt(dt, id.vars = c('sgRNA', 'sgrna_efficiency'), 
                measure.vars = setdiff(colnames(dt), c('sgRNA', 'sgrna_efficiency')))
    plots <- sapply(simplify = F, split(mdt, mdt$variable), function(x) {
      ggplot(x, aes(x = sgrna_efficiency, y = value)) + 
        geom_point() + geom_smooth(method = 'lm') + 
        labs(title = unique(x$variable), 
             subtitle = paste('pearsons_corr:', 
                              round(cor(x$sgrna_efficiency, x$value), 2))) + 
        theme_bw()
    })
    return(plots)
  }
}

dt <- merge(sgrna_stats, sgrna_scores_boolean, by = 'sgRNA')
plots <- plot_scores(dt, 'boolean')
p <- cowplot::plot_grid(plotlist = plots, ncol = 4)
ggsave(filename = 'sgrna_scores.boolean.pdf', plot = p, width = 10, height = 8, units = 'in')

dt <- merge(sgrna_stats, sgrna_scores_numeric, by = 'sgRNA')
dt[is.na(percent_in_injection)]$percent_in_injection <- mean(dt$percent_in_injection, na.rm = T)

plots <- plot_scores(dt[sgrna_efficiency >= 5], 'numeric')
p <- cowplot::plot_grid(plotlist = plots, ncol = 4)
ggsave(filename = 'sgrna_scores.numerical.pdf', plot = p, width = 10, height = 8, units = 'in')


sgRNA_scores <- merge(cbind(sgrna_scores_boolean[,1], 
                            apply(sgrna_scores_boolean[,-1], 2, function(x) {as.factor(as.numeric(x))})), 
                      sgrna_scores_numeric, by = 'sgRNA')

dt <- merge(sgrna_stats, sgRNA_scores, by = 'sgRNA')

df <- data.frame(dt[,-1], row.names = dt$sgRNA)

# remove lethal gene features # discuss with jj # predictions are driven by these features
# df <- df[,grep('lethal', colnames(df), invert = T)]

# one variable has a missing value, just impute that
fit <- ranger::ranger(sgrna_efficiency ~ ., df, num.trees = 1000, importance = 'permutation')

imp <- ranger::importance(fit)

ggplot(data.frame('importance' = imp, 'variable' = names(imp)), aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = 'identity') + coord_flip()

ggpubr::ggscatter(data.frame('pred' = ranger::predictions(fit, df[,-1]), 'actual' = df$sgrna_efficiency), 
                  x = 'actual', y = 'pred', add = 'reg.line', cor.coef = T, cor.coef.coord = c(15, 10), 
                  cor.coef.size = 6)

# pick top variables
top <- names(sort(imp, decreasing = T)[1:5])

fit2 <- ranger::ranger(sgrna_efficiency ~ ., subset(df, select = c('sgrna_efficiency', top)), 
                       num.trees = 1000, importance = 'permutation')

imp2 <- ranger::importance(fit2)
ggplot(data.frame('importance' = imp2, 'variable' = names(imp2)), 
       aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = 'identity') + coord_flip()

ggpubr::ggscatter(data.frame('pred' = ranger::predictions(fit2, df[,-1]), 'actual' = df$sgrna_efficiency), 
                  x = 'actual', y = 'pred', add = 'reg.line', cor.coef = T, cor.coef.coord = c(15, 10), 
                  cor.coef.size = 6)

# TODO: split test/train and report predictions 
# IDEA: if we split sgrna efficiency >5 or <5 then we see some scores are predictive


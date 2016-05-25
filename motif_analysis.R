library('yaml')
library('igraph')
library(digest)
library(data.table)
library(ggplot2)

data_dir = file.path('outputs/stats/31/invertebrate/')
stats_files = dir(data_dir, pattern='*.stats', full.names = TRUE)
stats_names = dir(data_dir, pattern='*.stats')
dir.create('motifs')

data = sapply(stats_files, FUN=function(d) { yaml.load_file(d)}, simplify = FALSE)
names(data) <- sapply(strsplit(stats_names, split='[.]'), "[[", 1)

hash <- function(s) {
  digest(s, algo='xxhash32')
}

getMotifs <- function(name, data) {
  entry <- data[[name]]
  size3motifs = entry$motifs$`3`
  size3graphs = sapply(size3motifs[[1]], hash)
  size3counts = size3motifs[[2]]
  
  size4motifs = entry$motifs$`4`
  size4graphs = sapply(size4motifs[[1]], hash)
  size4counts = size4motifs[[2]]
  
  data.frame(genome=rep(name, length(size3counts) + length(size4counts)),
             motif_size=c(rep(3, length(size3counts)), rep(4, length(size4counts))),
             motif_id=c(size3graphs, size4graphs),
             full_dot=c(names(size3graphs), names(size4graphs)),
             count=c(size3counts, size4counts),
             row.names = NULL)
}

convertDot <- function(dot) {
  # this is dumb and really slow but it works
  dot = as.character(dot)

  gmlfilename = file.path('motifs', paste(hash(dot), '.gml', sep=''))
  dotfilename = file.path('motifs', paste(hash(dot), '.dot', sep=''))

  write(dot, file=dotfilename)
  cmd = paste('gv2gml', dotfilename, '>', gmlfilename, sep = ' ')
  system(cmd)
  read_graph(gmlfilename, format='gml')
}


motif.data = rbindlist(sapply(names(data), getMotifs, data, USE.NAMES = FALSE, simplify = FALSE))
motif.data$graph <- sapply(motif.data$full_dot, convertDot, simplify = FALSE, USE.NAMES = FALSE)
motif.data$iso_class <- as.factor(sapply(motif.data$graph, isomorphism_class, USE.NAMES = FALSE))

p = ggplot(data = motif.data, aes(x=count))
p + facet_grid(iso_class ~ motif_size, labeller = label_both) + 
  scale_x_continuous(trans='log10') +
  geom_histogram(bins=100) +
  xlab('Number of Motifs') +
  ylab('Count')

ggplot(data = motif.data, aes(x=count, fill=iso_class)) +
  facet_grid(motif_size ~ ., labeller = label_both) +
  scale_x_continuous(trans='log10') +
  geom_density(alpha=0.3) +
  theme_bw()


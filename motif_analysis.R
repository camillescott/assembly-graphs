library('yaml')
library(digest)
library(data.table)
library(ggplot2)
library(grid)
library(gridGraphics)
library(gridBase)
library(gridExtra)
library('igraph')
library(Cairo)
library("RColorBrewer")
library(GGally)
library(intergraph)

data_dir = file.path('outputs/stats/31/invertebrate/')
stats_files = dir(data_dir, pattern='*.stats', full.names = TRUE)
stats_names = dir(data_dir, pattern='*.stats')

dir.create('motifs')
data = sapply(stats_files, FUN=function(d) { yaml.load_file(d)}, simplify = FALSE)
names(data) <- sapply(strsplit(stats_names, split='[.]'), "[[", 1)

kmer_files = dir(data_dir, pattern='*.uniquekmers', full.names = TRUE)
kmer_names = sapply(strsplit(dir(data_dir, pattern='*.uniquekmers'), split='[.]'), "[[", 1)

kmer.data = data.frame(genome=kmer_names,
                       unique_kmers=as.integer(sapply(strsplit(sapply(kmer_files, function (f) {readLines(f)[1]}, USE.NAMES = FALSE), 
                                                         split=' '),"[[", 1)))



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
  if (!file.exists(gmlfilename)) {
    write(dot, file=dotfilename)
    cmd = paste('gv2gml', dotfilename, '>', gmlfilename, sep = ' ')
    system(cmd)
  }
  read_graph(gmlfilename, format='gml')
}


motif.data = rbindlist(sapply(names(data), getMotifs, data, USE.NAMES = FALSE, simplify = FALSE))
motif.data$graph <- sapply(motif.data$full_dot, convertDot, simplify = FALSE, USE.NAMES = FALSE)
motif.data$iso_class <- as.factor(sapply(motif.data$graph, isomorphism_class, USE.NAMES = FALSE))

all.data = merge(motif.data, kmer.data, by='genome')

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

motif.graphs <- subset.data.frame(motif.data[ !duplicated(motif.data$iso_class)],
                                  select=c('graph', 'motif_size', 'iso_class'))
                       
mycolors <- brewer.pal(n = 8, name = 'Dark2')
names(mycolors) <- levels(motif.graphs$iso_class)

plotMotifs <- function (data) {
  plot.new()
  nrows = length(data$graph)
  ncols = 1
  par(mai=c(0,0,0,0), mar=c(0,0,0,0), bg=NA)
  layout(matrix(c(1:length(data$graph)), nrows, ncols, byrow=TRUE),
         widths=rep(1, ncols), heights = rep(1, nrows), respect=TRUE)
  layout.show(length(data$graph))
  for (i in 1:length(data$graph)) {
    iso = data$iso_class[[i]]
    g = data$graph[[i]]
    V(g)$label = ''
    E(g)$label = ''
    plot.igraph(g,
                layout=layout_on_grid,
                vertex.size = 30,
                vertex.color=mycolors[[as.character(iso)]],
                vertex.frame.color= "white",
                vertex.label.color = "white",
                vertex.label.family = "sans",
                edge.width=1,
                margin=c(1.0, 0, 1.0, 0),
                edge.color="black",
                asp=.75)
                #main=paste('Iso ', iso))
  }
  grid.echo()
  a <- grid.grab(wrap=TRUE)
  a
}
motifs.grob <- plotMotifs(motif.graphs)

ggplot(data = all.data, aes(x=unique_kmers, y=count)) +
  facet_wrap(~ iso_class, labeller = label_both, ncol=2) +
  geom_point() +
  theme_bw()



k.x.count <- ggplot(data = all.data, aes(x=unique_kmers, y=count, color=iso_class)) +
  geom_smooth(se = FALSE) +
  geom_point(alpha=0.5) +
  theme_minimal() + 
  scale_color_manual(values=mycolors,
                     name='Isoform Class') +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "bottom") +
  xlab('Nodes in Uncompressed Graph (Unique K-mers)') +
  ylab('Motif Count')

grid.arrange(k.x.count, 
             plotMotifs(motif.graphs),
             widths=c(12, 3),
             heights=c(3, 3, 3),
             layout_matrix=rbind(c(1,NA), 
                                 c(1,2),
                                 c(1,NA)))


grid.arrange(grobs=lapply(motif.graphs$graph, 
                          ggnet2, mode='circle',
                          node.size=5,
                          USE.NAMES = FALSE), 
             nrow=4, ncol=2, widths=c(1, 1), heights=rep(1, 4))


library('yaml')
library('digest')
library('data.table')
library('ggplot2')
library('grid')
library('gridGraphics')
library('gridBase')
library('gridExtra')
library('igraph')
library('Cairo')
library('RColorBrewer')
library('GGally')
library('intergraph')
library('psych')
library('reshape2')
library('gplots')
library('grImport')

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

#
# Load data!
#

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

motif.data = rbindlist(sapply(names(data), getMotifs, data, USE.NAMES = FALSE, simplify = FALSE))
motif.data$graph <- sapply(motif.data$full_dot, convertDot, simplify = FALSE, USE.NAMES = FALSE)
motif.data$iso_class <- as.factor(sapply(motif.data$graph, isomorphism_class, USE.NAMES = FALSE))
all.data = merge(motif.data, kmer.data, by='genome')

motif.graphs <- subset.data.frame(motif.data[ !duplicated(motif.data$iso_class)],
                                  select=c('graph', 'motif_size', 'iso_class'))
motif.colors <- brewer.pal(n = 8, name = 'Dark2')
names(motif.colors) <- levels(motif.graphs$iso_class)
motif.graphs$color <- motif.colors
motif.graphs$filename <- sapply(motif.graphs$iso_class, function (item) { paste('Iso', item, '.ps', sep='')})

p = ggplot(data = all.data, aes(x=count))
p + facet_grid(iso_class ~ motif_size, labeller = label_both) + 
  scale_x_continuous(trans='log10') +
  geom_histogram(bins=100) +
  xlab('Number of Motifs') +
  ylab('Count')

ggplot(data = all.data, aes(x=count, fill=iso_class)) +
  facet_grid(motif_size ~ ., labeller = label_both) +
  scale_fill_manual(values=motif.graphs$color, name='Isoform\nClass') +
  scale_x_continuous(trans='log10', name=expression('log'[10]~'Count')) +
  ylab(label='Density') +
  geom_density(alpha=0.3) +
  theme_bw()
ggsave('iso_class_density.png', dpi=300)



plotMotifs <- function (data) {
  plot.new()

  nrows = length(data$graph)
  ncols = 1
  #png(filename ='motifs.png', width=ncols, height=nrows, res=300)
  par(mfrow=c(nrows,1),
      pin=c(nrows, ncols),
      mai=c(0,0,0,0), 
      mar=c(0,0,0,0),
      oma=c(0,0,1,0),
      bg='white',
      xpd=NA)
  #png(filename ='motifs.png', width=ncols, height=nrows, res=300)
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
                edge.color="black",
                asp=1)
    mtext(paste('Iso', iso), side = 1, line = -5, adj = 0.3, cex = .8,col = "black")
  }
  #dev.off()
  #par(bg=NA)
  #plot.new()
  #dev.copy(png, 'motifs.png', width=ncols, height=nrows, res=300)
  #grid.echo()
  #a <- grid.grab(wrap=TRUE)
  #a
}
plotMotifs(motif.graphs)

saveMotifs <- function (data) {

  for (i in 1:length(data$graph)) {
    iso = data$iso_class[[i]]
    g = data$graph[[i]]
    V(g)$label = ''
    E(g)$label = ''
    filename = paste('Iso', iso, '.ps', sep='')
    postscript(file=filename,horiz=TRUE,onefile=TRUE,width=2,height=2)
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
                asp=1)
    #main=paste('Iso ', iso))
    
  }
}

saveMotif <- function(g, filename, node.color, main='') {
  V(g)$label = ''
  E(g)$label = ''
  #postscript(file=filename,horiz=TRUE,onefile=TRUE,width=4,height=4)
  dev.new()
  plot.new()
  par(mar=c(0,0,0,0))
  plot.igraph(g,
              layout=layout_on_grid,
              vertex.size = 30,
              vertex.color=node.color,
              vertex.frame.color= "white",
              vertex.label.color = "white",
              vertex.label.family = "sans",
              edge.width=1,
              margin=c(1.0, 0, 1.0, 0),
              edge.color="black",
              asp=1,
              main=main)

  #dev.copy2eps(file=filename, width=4)
  dev.print(file=filename, width=4)

}

for (i in 1:nrow(motif.graphs)) {
  saveMotif(motif.graphs[[i, 'graph']], 
            motif.graphs[[i, 'filename']], 
            motif.graphs[[i, 'color']])
}

motif.glyphs <- sapply(motif.graphs$filename, function(fn) {
  out_fn <- paste(fn, '.xml', sep='')
  PostScriptTrace(fn, out_fn)
  readPicture(out_fn)
})
names(motif.glyphs) <- sapply(motif.graphs$iso_class, function (n) { paste('Iso', n)})

motif.glyph.axis <- function () {
    function(label, x=0.5, y=0.5, ...) {
      print(label)
      ggplot2:::absoluteGrob(
        do.call("gList", mapply(symbolsGrob, motif.glyphs[label], x, y, SIMPLIFY = FALSE)),
        height = unit(1.5, "cm"))
    }
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

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  cormat
}

all.melted <- melt(all.data, id=c('genome', 'iso_class'), measure.vars=c('count'))
motif.data <- as.data.frame(cast(all.melted, genome ~ iso_class, fill=0)[2:9])
names(motif.data) <- sapply(names(motif.data), function (n) { paste('Iso', n)})
names(motif.data) <- as.factor(names(motif.data))
motif.corr <- reorder_cormat(round(cor(motif.data), 3))

corrplot(motif.corr, 
         method = "circle", 
         tl.cex = 0.93, 
         order = "hclust", 
         ddrect = 3,
         mar=rep(3,4), 
         cl.ratio=0.1, 
         cl.length=5, 
         title='Correlation Between Isoform Classes', 
         tl.col='black', 
         cl.pos='b')

import tempfile

import numpy as np
from pathlib import Path
import yaml

configfile: "config.json"

rule all:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.{suffix}",
             ksize=[31],
             genome=config['genomes'],
             suffix=['info'])

rule stats:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.stats",
             ksize=[31], genome=config['genomes'])

rule kmers:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.uniquekmers",
             ksize=[31], genome=config['genomes'])

rule unique_kmers:
    input: "genomes/refseq/{group}/{genome}.fna.gz"
    output: "outputs/stats/{ksize}/{group}/{genome}.uniquekmers"
    params: ksize="{ksize}"
    threads: 4
    shell: """
        OMP_NUM_THREADS={threads} unique-kmers.py -k {params.ksize} -e 0.01 --report {output} {input}
    """

rule gt_motif_plots:
    input: info="outputs/stats/31/invertebrate/GCF_000004095.1_Hydra_RP_1.0_genomic.info"
    output: "outputs/stats/motifs/{mid}.png".expand(
        mid=["3.0", "3.1",
             "4.0", "4.1", "4.2", "4.3", "4.4", "4.5"]
    run:
        import graph_tool as gt
        import graph_tool.draw

        with open(input.info, 'r') as f:
            data = yaml.load(f)
            motifs = {
                3: [None] * 2,
                4: [None] * 6
            }
            for msize in data['motifs']:
                for i, motif in enumerate(data['motifs'][msize][0]):
                    with tempfile.TemporaryFile() as f:
                        f.write(motif.encode('utf-8'))
                        f.seek(0)
                        mg = gt.load_graph(f, fmt='dot')
                        motifs[msize][i] = mg

        for k in motifs:
            for i, motif in enumerate(motifs[k]):
                pos = gt.draw.arf_layout(motif)
                gt.draw.graph_draw(motif, pos=pos,
                                vertex_size=10,
                                inline=False,
                                output_size=(100,100),
                                output="outputs/stats/motifs/{}.{}.png".format(k, i))


rule gt_describe:
    input: graph="outputs/cdbg/{ksize}/{group}/{genome}.gml.gz"
    output: stats="outputs/stats/{ksize}/{group}/{genome}.desc"
    params: ksize="{ksize}"
    threads: 8
    run:
        import graph_tool as gt
        import graph_tool.stats
        gt.openmp_set_num_threads(threads)

        g = gt.load_graph(input.graph, fmt="gml")
        gt.stats.remove_paralllel_edges(g)
        hist = gt.stats.vertex_hist(g, "total")
        props = {
          'nodes': g.num_vertices(),
          'edges': g.num_edges(),
          'average_degree': gt.stats.vertex_average(g, "total")[0],
          'degree_distribution': hist[0].tolist()
        }

        with open(output.stats, 'w') as stats:
            stats.write(yaml.safe_dump(props))

rule gt_stats:
    input: graph="outputs/cdbg/{ksize}/{group}/{genome}.gml.gz"
    output: stats="outputs/stats/{ksize}/{group}/{genome}.stats"
    params: ksize="{ksize}"
    threads: 4
    run:
        import graph_tool as gt
        import graph_tool.clustering
        gt.openmp_set_num_threads(threads)

        g = gt.load_graph(input.graph, fmt="gml")
        motifs = {}
        for k in (3, 4):
            graphs, counts = gt.clustering.motifs(g, k=k)
            gdots = []
            for sg in graphs:
                with tempfile.TemporaryFile() as f:
                    sg.save(f, fmt='dot')
                    f.seek(0)
                    gdots.append(f.read().decode('utf-8'))

            motifs[k] = (gdots, counts)

        psizes = []
        for v in g.vertices():
            nsize = g.vp['size'][v]
            if nsize > int(params.ksize):
                psizes.append(nsize)

        counts, bins = np.histogram(psizes,
                                    range=(int(params.ksize), 10000),
                                    bins=100)
        props = {
          'motifs': motifs,
          'linear_paths': {
            'counts': counts.tolist(),
            'bins': bins.tolist(),
          }
        }

        with open(output.stats, 'w') as stats:
            stats.write(yaml.safe_dump(props))

rule gt_info:
    input:
        stats="outputs/stats/{ksize}/{group}/{genome}.stats",
        desc="outputs/stats/{ksize}/{group}/{genome}.desc",
        unique="outputs/stats/{ksize}/{group}/{genome}.uniquekmers"
    output: "outputs/stats/{ksize}/{group}/{genome}.info"
    shell: """
        cat {input.stats} {input.desc} \
            <(echo -n "unique_kmers: " & head -1 {input.unique} | cut -d" " -f1) > {output}
    """


rule install_software:
    output:
        "sw/spacegraphcats",
        "sw/khmer"
    shell: """
        mkdir -p sw
        cd sw
        git clone git@github.com:spacegraphcats/spacegraphcats
        git clone git@github.com:dib-lab/khmer -b minhash
    """

rule compress_dbg:
    input: "genomes/refseq/{group}/{genome}.fna.gz"
    output: cdbg="outputs/cdbg/{ksize}/{group}/{genome}.gml.gz"
    params:
        ksize="{ksize}"
    run:
        GML = os.path.splitext(output.cdbg)[0]
        shell("mkdir -p $(dirname {output})")
        shell("zcat {input} | sed -e '/^>/! s/[RYSWKMBDHV]/N/g' > {output}.tmp")
        shell("./sw/khmer/sandbox/extract-compact-dbg.py -o " + GML +
              " -k {params.ksize} -x 2e10 {output}.tmp")
        shell("rm {output}.tmp")
        shell("gzip " + GML)

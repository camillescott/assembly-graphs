import tempfile

import numpy as np
from pathlib import Path
import yaml

configfile: "config.json"

rule all:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.{suffix}",
             ksize=[31], genome=config['genomes'], suffix=['stats', 'uniquekmers'])

rule stats:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.stats",
             ksize=[31], genome=config['genomes'])

rule kmers:
    input:
      expand("outputs/stats/{ksize}/invertebrate/{genome}.uniquekmers",
             ksize=[31], genome=config['genomes'])

rule unique_kmers:
    input: "genomes/refseq/invertebrate/{genome}.fna.gz"
    output: "outputs/stats/{ksize}/invertebrate/{genome}.uniquekmers"
    params: ksize="{ksize}"
    threads: 4
    shell: """
        OMP_NUM_THREADS={threads} unique-kmers.py -k {params.ksize} -e 0.01 --report {output} {input}
    """

rule gt_stats:
    input: graph="outputs/cdbg/{ksize}/invertebrate/{genome}.gml.gz"
    output: stats="outputs/stats/{ksize}/invertebrate/{genome}.stats"
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
    input: "genomes/refseq/invertebrate/{genome}.fna.gz"
    output: cdbg="outputs/cdbg/{ksize}/invertebrate/{genome}.gml.gz"
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

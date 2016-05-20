import networkx as nx
from pathlib import Path

configfile: "config.json"

rule all:
    input:
      expand("outputs/cdbg/{ksize}/invertebrate/{genome}.gml.gz",
             ksize=[31], genome=config['genomes'])

'''
rule unique_kmers:
    output: "outputs/{organism}.{ksize}.unique"
    input: "inputs/{organism}.fa.gz"
    params: ksize="{ksize}"
    threads: 4
    shell: """
        OMP_NUM_THREADS={threads} unique-kmers.py -k {params.ksize} -e 0.01 --report {output} {input}
    """

rule nx_stats:
    input: graph="outputs/cdbg/{ksize}/vertebrate_mammalian/GCF_000001405.33/GCF_000001405.33_GRCh38.p7_genomic.gml.gz"
    output: stats="outputs/{organism}.{ksize}.stats"
    run:
        g = nx.read_graphml(input.graph)
        g.name = "Directed version"
        gu = g.to_undirected()
        gu.name = "Undirected version"
        with open(output.stats, 'w') as stats:
            stats.write(nx.info(g) + "\n\n")
            stats.write(nx.info(gu) + "\n\n")
'''

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
    output: "outputs/cdbg/{ksize}/invertebrate/{genome}.gml.gz"
    params: ksize="{ksize}"
    shell: """
        mkdir -p $(dirname {output})
        zcat {input} | sed -e '/^>/! s/[RYSWKMBDHV]/N/g' > {output}.tmp
        ./sw/khmer/sandbox/extract-compact-dbg.py -o {output} \
            -k {params.ksize} -x 2e10 {output}.tmp
        rm {output}.tmp
    """

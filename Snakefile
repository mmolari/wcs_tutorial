configfile: "config/config.yaml"


with open(config["strains"], "r") as f:
    strains = f.read().splitlines()
strains_subset = config["strains-subset"]
strain_pair = config["marginalize"]


rule download_gbk:
    output:
        "data/ST131_gbk/{acc}.gbk",
    conda:
        "config/conda_env.yml"
    shell:
        """
        ncbi-acc-download {wildcards.acc} -e all -F genbank --out {output}
        """


rule gbk_to_fa:
    input:
        rules.download_gbk.output,
    output:
        "data/ST131_fa/{acc}.fa",
    conda:
        "config/conda_env.yml"
    shell:
        """
        python3 scripts/gbk_to_fa.py --gbk {input} --fa {output}
        """


rule map_bla:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains),
        bla=config["bla-file"],
    output:
        "results/bla15/map.paf",
    conda:
        "config/conda_env.yml"
    shell:
        """
        minimap2 -x asm5 {input.bla} {input.fa} > {output}
        """


rule extract_window:
    input:
        paf=rules.map_bla.output,
        fa=rules.map_bla.input.fa,
    output:
        "results/bla15/extracted_window_{w}.fa",
    conda:
        "config/conda_env.yml"
    params:
        L=int(config["bla-len"] * 0.95),
    shell:
        """
        python3 scripts/extract_matches.py \
            --in_fa {input.fa} \
            --paf {input.paf} \
            --window {wildcards.w} \
            --length {params.L} \
            --out {output}
        """


rule build_window_pangraph:
    input:
        expand(rules.extract_window.output, w=config["window-size"]),
    output:
        "results/bla15/pangraph_window.json",
    conda:
        "config/conda_env.yml"
    shell:
        """
        pangraph build -l 50 -a 10 -b 10 -s 5 {input} > {output}
        """


rule export_window_pangraph:
    input:
        rules.build_window_pangraph.output,
    output:
        directory("results/bla15/export"),
    conda:
        "config/conda_env.yml"
    shell:
        """
        pangraph export \
            -nd \
            --edge-minimum-length 0 \
            -o {output} {input} 
        """


rule extract_alignment:
    input:
        expand(rules.extract_window.output, w=0),
    output:
        "results/bla15/bla_alignment.fa",
    conda:
        "config/conda_env.yml"
    shell:
        """
        mafft --auto --adjustdirection {input} > {output}
        """


rule build_subset_pangraph:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains_subset),
    output:
        "results/pangraph/subset.json",
    conda:
        "config/conda_env.yml"
    shell:
        """
        JULIA_NUM_THREADS=3
        pangraph build --circular -a 20 -b 5 -s 20 {input.fa} > {output}
        """


rule export_subset_pangraph:
    input:
        rules.build_subset_pangraph.output,
    output:
        directory("results/pangraph/export/subset"),
    conda:
        "config/conda_env.yml"
    shell:
        """
        pangraph export -nd -o {output} {input} 
        """


rule marginalize:
    input:
        rules.build_subset_pangraph.output,
    output:
        "results/pangraph/marginalized.json",
    conda:
        "config/conda_env.yml"
    params:
        s1=strain_pair[0],
        s2=strain_pair[1],
    shell:
        """
        pangraph marginalize \
            --strains {params.s1},{params.s2} \
            {input} > {output}
        """


rule export_marginal_graph:
    input:
        rules.marginalize.output,
    output:
        directory("results/pangraph/export/marginalized"),
    conda:
        "config/conda_env.yml"
    shell:
        """
        pangraph export -nd -o {output} {input} 
        """


rule all:
    input:
        rules.export_window_pangraph.output,
        rules.extract_alignment.output,
        rules.export_subset_pangraph.output,
        rules.marginalize.output,
        rules.export_marginal_graph.output,

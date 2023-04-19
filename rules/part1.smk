
rule map_bla:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains),
        bla=config["bla-file"],
    output:
        "results/bla15/map.paf",
    conda:
        "../config/conda_env.yml"
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
        "../config/conda_env.yml"
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
        "../config/conda_env.yml"
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
        "../config/conda_env.yml"
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
        "../config/conda_env.yml"
    shell:
        """
        mafft --auto --adjustdirection {input} > {output}
        """


rule bla_assign_color_to_isolate:
    input:
        paf=rules.map_bla.output,
        tree=config["coregenome-tree"],
    output:
        color="results/bla15/isolate_color.csv",
        fig="figs/bla_coretree.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/assign_leaves_color.py \
            --paf {input.paf} \
            --tree {input.tree} \
            --fig {output.fig} \
            --color_csv {output.color}
        """


rule find_bla_block:
    input:
        exp=rules.export_window_pangraph.output,
        fa=config["bla-file"],
    output:
        "results/bla15/bla_block.txt",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        minimap2 {input.fa} {input.exp}/pangraph.fa \
            | awk '{{print $1}}' \
            > {output}
        """


rule part1_all:
    input:
        rules.export_window_pangraph.output,
        rules.extract_alignment.output,
        rules.bla_assign_color_to_isolate.output,
        rules.find_bla_block.output,


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
            | head -1 | awk '{{print $1}}' \
            > {output}
        """


rule bla_structural_diversity:
    input:
        pan=rules.build_window_pangraph.output,
        block=rules.find_bla_block.output,
        leaves_col=rules.bla_assign_color_to_isolate.output.color,
    output:
        shared_L="results/bla15/shared_bla_length.csv",
        fig_paths="figs/bla_paths_drawing.png",
        fig_matrix="figs/bla_paths_shared_len.png",
        colors="results/bla15/block_colors.csv",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/bla_structural_diversity.py \
            --pangraph {input.pan} \
            --anchor_block $(cat {input.block}) \
            --leaves_colors {input.leaves_col} \
            --fig_paths {output.fig_paths} \
            --fig_matrix {output.fig_matrix} \
            --block_colors {output.colors} \
            --shared_len_df {output.shared_L}
        """


rule bla_structure_vs_coretree:
    input:
        shared_L=rules.bla_structural_diversity.output.shared_L,
        tree=config["coregenome-tree"],
        leaves_col=rules.bla_assign_color_to_isolate.output.color,
    output:
        fig_scatter="figs/bla_shared_len_vs_coretree_scatter.png",
        fig_tree="figs/bla_shared_len_vs_coretree.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/shared_paths_vs_coretree.py \
            --shared_len_df {input.shared_L} \
            --tree {input.tree} \
            --leaves_colors {input.leaves_col} \
            --fig_scatter {output.fig_scatter} \
            --fig_tree {output.fig_tree}
        """


rule part1_all:
    input:
        rules.extract_alignment.output,
        rules.bla_structure_vs_coretree.output,

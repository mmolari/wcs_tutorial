rule marginalize:
    input:
        rules.build_subset_pangraph.output,
    output:
        "results/pangraph/marginalized.json",
    conda:
        "../config/conda_env.yml"
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
        directory("results/pangraph/export_marginalized"),
    conda:
        "../config/conda_env.yml"
    shell:
        """
        pangraph export -nd -ell 0 -o {output} {input} 
        """


rule plot_graph_projection:
    input:
        rules.marginalize.output,
    output:
        "figs/graph_projection.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/plot_graph_projection.py \
            --pangraph {input} \
            --fig {output}
        """


rule plot_bandage_marginal:
    input:
        rules.export_marginal_graph.output,
    output:
        "figs/bandage_marginal.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        bandage image {input}/pangraph.gfa {output} \
            --iter 4 \
            --colour depth \
            --depvallow 1 \
            --depvalhi 2 \
            --width 800
        """


rule plot_selected_pair_in_subtree:
    input:
        tree=config["coregenome-tree"],
    output:
        "figs/selected_pair_in_subtree.png",
    conda:
        "../config/conda_env.yml"
    params:
        clade=" ".join(config["strains-subset"]),
        pair=" ".join(config["marginalize"]),
    shell:
        """
        python3 scripts/plot_selected_pair_in_subtree.py \
            --tree {input.tree} \
            --clade {params.clade} \
            --pair {params.pair} \
            --fig {output}
        """


rule part3_all:
    input:
        rules.plot_selected_pair_in_subtree.output,
        rules.export_marginal_graph.output,
        rules.plot_graph_projection.output,

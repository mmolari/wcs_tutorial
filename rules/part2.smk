
rule build_subset_pangraph:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains_subset),
    output:
        "results/pangraph/subset.json",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        JULIA_NUM_THREADS=3
        pangraph build --circular -a 20 -b 5 -s 20 {input.fa} > {output}
        """


rule export_subset_pangraph:
    input:
        rules.build_subset_pangraph.output,
    output:
        directory("results/pangraph/export_subset"),
    conda:
        "../config/conda_env.yml"
    shell:
        """
        pangraph export -nd -o {output} -ell 0 {input} 
        """


rule plot_block_distr:
    input:
        rules.build_subset_pangraph.output,
    output:
        "figs/block_distr.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/plot_block_distr.py \
            --pangraph {input} \
            --fig {output}
        """


# rule bandage_subset:
#     input:
#         rules.export_subset_pangraph.output,
#     output:
#         "figs/bandage_subset.png",
#     shell:
#         """
#         bandage image {input}/pangraph.gfa {output} \
#             --nodewidth 5 \
#             --iter 4 \
#             --colour depth \
#             --depvallow 1 \
#             --depvalhi 10
#         """


rule private_seq_distance:
    input:
        rules.build_subset_pangraph.output,
    output:
        "results/pangraph/private_seq_distance.csv",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/pairwise_private_seq.py \
            --pangraph {input} \
            --dist_df {output}
        """


rule plot_private_seq:
    input:
        dist_df=rules.private_seq_distance.output,
        tree=config["coregenome-tree"],
    output:
        fig_scatter="figs/private_seq_scatter.png",
        fig_matrix="figs/private_seq_matrix.png",
    conda:
        "../config/conda_env.yml"
    shell:
        """
        python3 scripts/plot_private_seq.py \
            --dist_df {input.dist_df} \
            --tree {input.tree} \
            --fig_scatter {output.fig_scatter} \
            --fig_matrix {output.fig_matrix}
        """


rule part2_all:
    input:
        rules.build_subset_pangraph.output,
        rules.export_subset_pangraph.output,
        rules.plot_block_distr.output,
        rules.private_seq_distance.output,
        rules.plot_private_seq.output,
        # rules.bandage_subset.output,

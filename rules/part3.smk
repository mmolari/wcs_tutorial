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
        directory("results/pangraph/export/marginalized"),
    conda:
        "../config/conda_env.yml"
    shell:
        """
        pangraph export -nd -o {output} {input} 
        """


rule part3_all:
    input:
        rules.export_marginal_graph.output,

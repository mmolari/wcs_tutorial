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


include: "rules/part1.smk"
include: "rules/part2.smk"


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
        rules.part1_all.output,
        rules.part2_all.output,
        # rules.export_marginal_graph.output,

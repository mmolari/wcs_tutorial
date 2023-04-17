configfile: "config/config.yaml"


with open(config["strains"], "r") as f:
    strains = f.read().splitlines()
with open(config["strains-subset"], "r") as f:
    strains_subset = f.read().splitlines()


rule download_gbk:
    output:
        "data/ST131_gbk/{acc}.gbk",
    conda:
        "conda/wcs_tutorial.yml"
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
        "conda/wcs_tutorial.yml"
    shell:
        """
        python3 scripts/gbk_to_fa.py --gbk {input} --fa {output}
        """


rule map_bla:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains),
        bla="data/b_lactam/{bl}.fa",
    output:
        "results/map/{bl}.paf",
    conda:
        "conda/wcs_tutorial.yml"
    shell:
        """
        minimap2 {input.bla} {input.fa} -x asm5 > {output}
        """


rule extract_window:
    input:
        paf=rules.map_bla.output,
        fa=rules.map_bla.input.fa,
    output:
        "results/window/{bl}__{w}.fa",
    conda:
        "conda/wcs_tutorial.yml"
    params:
        L=lambda w: int(config["bla_L"][w.bl] * 0.95),
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
        rules.extract_window.output,
    output:
        "results/pangraph/window__{bl}__{w}.json",
    conda:
        "conda/wcs_tutorial.yml"
    shell:
        """
        pangraph build -l 50 -a 10 -b 10 -s 5 {input} > {output}
        """


rule export_window_pangraph:
    input:
        rules.build_window_pangraph.output,
    output:
        directory("results/pangraph/export/window__{bl}__{w}"),
    conda:
        "conda/wcs_tutorial.yml"
    shell:
        """
        pangraph export \
            -nd \
            --edge-minimum-length 0 \
            -o {output} {input} 
        """


rule extract_alignment:
    input:
        expand(rules.extract_window.output, w="0", allow_missing=True),
    output:
        "results/alignment/{bl}.fa",
    conda:
        "conda/wcs_tutorial.yml"
    shell:
        """
        mafft --auto --adjustdirection {input} > {output}
        """


rule build_alignment_tree:
    input:
        rules.extract_alignment.output,
    output:
        "results/alignment/{bl}.nwk",
    conda:
        "conda/wcs_tutorial.yml"
    shell:
        """
        fasttree -nt {input} > {output}
        """


rule build_subset_pangraph:
    input:
        fa=expand(rules.gbk_to_fa.output, acc=strains_subset),
    output:
        "results/pangraph/subset.json",
    conda:
        "conda/wcs_tutorial.yml"
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
        "conda/wcs_tutorial.yml"
    shell:
        """
        pangraph export -nd -o {output} {input} 
        """


rule all:
    input:
        # expand(rules.download_gbk.output, acc=strains),
        expand(rules.export_window_pangraph.output, bl="bla15", w="5000"),
        expand(rules.build_alignment_tree.output, bl="bla15"),
        # expand(rules.extract_alignment.output, bl="bla15"),
        rules.export_subset_pangraph.output,

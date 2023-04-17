

## setup

- install pangraph and have available in path
- snakemake

## data

bla-ctx-m-15 uniprot: https://www.uniprot.org/uniprotkb/G8FPM5/entry#sequences
downloaded fasta: JN019833.1 Escherichia coli strain CTX_M_925_10 extended spectrum beta-lactamase CTX-M-15 (bla CTX-M-15) gene, partial cds

## workflow

look for matches:
```bash
minimap2 bla_fasta.fa ST131_fa/*.fa > bla_map.paf
```

download genbank:
```bash
mkdir ST131_gbk
while read line
do
    ncbi-acc-download $line -e all -F genbank --out ST131_gbk/$line.gbk
done < 'strains.txt'
```

build a pangraph:
```bash
strains=""
while read line
do
    strains="$strains $line"
done < 'strains.txt'
JULIA_NUM_THREADS=4
pangraph build --circular -a 100 -b 5 -s 10 $strains > pangraph.json
```

export:
```bash
pangraph export -o bla_export pangraph.json
```

cut-bla:
```bash
python3 scripts/extract_matches.py --fld ST131_fa --paf bla_map.paf \
    -w 5000 \
    -l 400 \
    --out bla_window.fa
JULIA_NUM_THREADS=4
pangraph build -a 20 -b 5 -s 20 bla_window.fa > pangraph_bla.json
pangraph export -o bla_export pangraph_bla.json
```
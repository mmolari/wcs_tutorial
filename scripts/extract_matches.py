# script to extract sequences from a fasta file around matches identified by a
# paf file

import argparse
import pathlib
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Extract sequences from a fasta file around matches identified by a paf file"
)
parser.add_argument("--in_fa", type=str, nargs="+", help="input fasta file(s)")
parser.add_argument("--paf", help="paf file")
parser.add_argument("--out", help="output fasta file")
parser.add_argument(
    "-l", "--length", type=int, default=10, help="minimum length of matches"
)
parser.add_argument(
    "-w", "--window", type=int, default=0, help="window size around matches"
)

args = parser.parse_args()

# read paf file
with open(args.paf, "r") as f:
    paf_lines = f.readlines()

# filter paf file
pafs = []
for n_paf, paf in enumerate(paf_lines):
    paf = paf.strip().split("\t")
    seq_id = paf[0]
    match_len = int(paf[10])
    if match_len >= args.length:
        pafs.append(
            {
                "seq_id": seq_id,
                "match_id": n_paf,
                "match_len": match_len,
                "start": int(paf[2]),
                "end": int(paf[3]),
                "L": int(paf[1]),
            }
        )

# read input fasta files
fa_in = {}
for in_file in args.in_fa:
    seq = SeqIO.read(in_file, format="fasta").seq
    fa_in[pathlib.Path(in_file).stem] = seq

# extract relevant part of the fasta file
records = []
for paf in pafs:
    sid = paf["seq_id"]
    nid = paf["match_id"]
    L = paf["L"]
    b, e = paf["start"], paf["end"]
    B, E = b - args.window, e + args.window
    if B < 0:
        B = 0
    if E > L:
        E = L
    seq = fa_in[sid]
    rec = SeqIO.SeqRecord(seq[B:E], id=f"{sid}-[{b}:{e}]", description="")
    records.append(rec)

# write output fasta file
with open(args.out, "w") as f:
    SeqIO.write(records, f, format="fasta")

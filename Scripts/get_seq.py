### This script takes a list of sequence identifiers and a reference fasta file and outputs the fasta sequences of the given identifieres
### odr so...

import argparse, os, sys, re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='input', help="FASTA file", required=True)
parser.add_argument('-id', '--id', dest='id', help="Seqid", required=True)
parser.add_argument('-o', '--output', dest='output', help="name of output file, default = output.txt", default='output.txt')

args = parser.parse_args()


wanted = args.id

fasta_sequences = SeqIO.parse(open(args.input),'fasta')
with open(args.output, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")

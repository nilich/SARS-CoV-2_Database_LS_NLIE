import argparse, os, sys, re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='input', help="FASTA file", required=True)
parser.add_argument('-s', '--sample', dest='sample', help="File with one SequenceID per line", required=True)
parser.add_argument('-o', '--output', dest='output', help="name of output file, default = results.txt", default='results.txt')

args = parser.parse_args()

sample = args.sample
nCount = sum([rec.seq.upper().count("N") for rec in SeqIO.parse(args.input, "fasta")])
l = sum([len(rec) for rec in SeqIO.parse(args.input, "fasta")])

perc = round(nCount/l*100, 2)

output = open(args.output, 'w' )
output.write("{}{}{}{}".format(sample,",",perc,"\n"))
output.close()

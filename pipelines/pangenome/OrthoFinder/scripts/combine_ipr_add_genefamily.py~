#!/usr/bin/env python3

import argparse, re

parser = argparse.ArgumentParser(description='Add some additional cols to Interpro TSV file.')
parser.add_argument('--ipr', help='Interpro TSV file',required=False,default='all.pan.genome.iprout.tsv')
parser.add_argument('--peps',help="Protein sequence file with extra headers to include",required=False,default="rep_seqs.pep")
parser.add_argument('--out',help="Output file for updated IPR results",required=False)


args = parser.parse_args()

if not args.out:
    args.out = args.ipr + ".updated"

pepid = {}
capture = re.compile(r'^>(\S+)\s+(\S+)\s+(\S+)')
with open(args.peps,"r") as inpep:
    for line in inpep:
        if line.startswith(">"):
            m = capture.match(line)
            if m:
                pepid[m.group(1)] = [m.group(2), m.group(3)]

getid = re.compile(r'^(\S+)')
with open(args.ipr,"r") as inipr, open(args.out,"wt") as outipr:
    for line in inipr:
        m = getid.match(line)
        if m:
            genename = m.group(1)
            if genename in pepid:
                outipr.write("\t".join(pepid[genename]) + "\t" + line)
        else:
            print("Not matching a line: {}".format(line))

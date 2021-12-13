#!/usr/bin/env python3

import csv, re, sys, os
import xml.etree.ElementTree as ET
from Bio import Entrez
Entrez.email = 'jason.stajich@ucr.edu'
insamples = "samples.csv"
outsamples="samples_prefix.csv"

if len(sys.argv) > 1:
    insamples = sys.argv[1]

if len(sys.argv) > 2:
    outsamples = sys.argv[2]

seen = {}
if os.path.exists(outsamples):
    with open(outsamples,"rU") as preprocess:
        incsv = csv.reader(preprocess,delimiter=",")
        h = next(incsv)
        for row in incsv:
            seen[row[0]] = row

with open(insamples,"rU") as infh, open(outsamples,"w") as outfh:
    outcsv    = csv.writer(outfh,delimiter=",")
    outcsv.writerow(['SPECIES','STRAIN','JGILIBRARY',
                     'BIOSAMPLE','BIOPROJECT','SRA','LOCUSTAG'])

    samplescsv = csv.reader(infh,delimiter=",")
    for row in samplescsv:
        outrow = [row[5],row[2]]
        strain = row[5]
        if strain in seen:
            outrow = seen[strain]
            outcsv.writerow(outrow)
            continue
        handle = Entrez.esearch(db="biosample",retmax=10,term=strain)
        record = Entrez.read(handle)
        handle.close()
        SRA = ""
        BIOSAMPLE = ""
        BIOPROJECT = ""
        LOCUSTAG = ""
        BIOPROJECTID=""
        for biosampleid in record["IdList"]:
            handle = Entrez.efetch(db="biosample", id=biosampleid)
            tree = ET.parse(handle)
            root = tree.getroot()
            for sample in root:
                BIOSAMPLE = sample.attrib['accession']
                for ids in root.iter('Ids'):
                    for id in ids.iter('Id'):
                        if id.attrib['db'] == "SRA":
                            SRA = id.text
                for links in root.iter('Links'):
                    for link in links:
                        linkdat = link.attrib
                        if linkdat['type'] == 'entrez':
                            BIOPROJECT = linkdat['label']
                            BIOPROJECTID = link.text
        if BIOPROJECTID:
            bioproject_handle = Entrez.efetch(db="bioproject",id = BIOPROJECTID)
            projtree = ET.parse(bioproject_handle)
            projroot = projtree.getroot()

            lt = projroot.iter('LocusTagPrefix')
            for locus in lt:
                LOCUSTAG = locus.text
        outrow.extend([BIOSAMPLE,BIOPROJECT,SRA,LOCUSTAG])
        outcsv.writerow(outrow)

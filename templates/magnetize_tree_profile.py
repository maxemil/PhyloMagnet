#! ${params.python3}
import pandas as pd
import sys
from ete3 import ncbi_taxonomy


def reverse_taxid_dict(name2tax):
    tax2name = {}
    for k,val in name2tax.items():
        for v in val:
            tax2name[v] = k
    return tax2name


def get_lineages(tax_map, rank):
    ncbi = ncbi_taxonomy.NCBITaxa()
    lineages = set()
    for line in open(tax_map):
        lineage_names = [l.strip().strip(';') for l in line.strip().split('\t')[1].split('; ')]
        name2tax = ncbi.get_name_translator(lineage_names)
        tax2name = reverse_taxid_dict(name2tax)
        lineage_taxids = [taxid for l in name2tax.values() for taxid in l]
        for k,v in ncbi.get_rank(lineage_taxids).items():
            if v == rank:
                lineages.add(tax2name[k])
    return lineages


def run(profile, lineages, samplename, genename):
    df = pd.read_csv(profile, sep='\\t')
    df['taxopath'] = df['taxopath'].astype(str)
    with open('decision_%s.txt' % genename, 'w') as out:
        for lineage in lineages:
            lineage_results = df[df['taxopath'].apply(lambda x: lineage in x and not '{};'.format(lineage) in x)]
            if lineage_results['aLWR'].item() >= 1:
                print("{} present in sample {} (tree {})".format(lineage, samplename, genename))
                print("{}\\t{}\\t{}\\tTrue".format(samplename, genename, lineage), file=out)
            else:
                print("{} NOT present in sample {} (tree {})".format(lineage, samplename, genename))
                print("{}\\t{}\\t{}\\tFalse".format(samplename, genename, lineage), file=out)


lineages = get_lineages("$tax_map", "$lineage")
run("$profile", lineages, "${profile.baseName.tokenize('-')[0]}", "${profile.baseName.tokenize('-')[1]}")

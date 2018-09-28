#! ${params.python3}
import pandas as pd
import sys


def is_lineage_in_references(ref_taxonmy, lineage):
    for line in open(ref_taxonmy):
        if lineage in line:
            return True
    return False


def check_references_lineage(ref_taxonmy, lineage, samplename, genename):
    if not is_lineage_in_references(ref_taxonmy, lineage):
        with open('decision_%s.txt' % genename, 'w') as out:
            print("Cannot assess whether {} is present in sample {} (tree {}), as references do not contain this taxonomic label..".format(lineage, samplename, genename))
            print("{}\\t{}\\t{}\\tna".format(samplename, genename, lineage), file=out)
            sys.exit()


def run(profile, lineage, samplename, genename):
    df = pd.read_csv(profile, sep='\\t')
    df['taxopath'] = df['taxopath'].astype(str)
    lineage_results = df[df['taxopath'].apply(lambda x: lineage in x and not '{};'.format(lineage) in x)]
    with open('decision_%s.txt' % genename, 'w') as out, open('decision_%s.log' % genename, 'w') as log:
        try:
            if lineage_results['aLWR'].item() >= ${params.aLWR_threshold}:
                print("{} present in sample {} (tree {})".format(lineage, samplename, genename), file=log)
                print("{}\\t{}\\t{}\\tTrue".format(samplename, genename, lineage), file=out)
            else:
                print("{} NOT present in sample {} (tree {})".format(lineage, samplename, genename), file=log)
                print("{}\\t{}\\t{}\\tFalse".format(samplename, genename, lineage), file=out)
        except ValueError:
            print("{} NOT present in sample {} (tree {})".format(lineage, samplename, genename), file=log)
            print("{}\\t{}\\t{}\\tFalse".format(samplename, genename, lineage), file=out)


check_references_lineage("$tax_map", "$lineage", "${profile.baseName.tokenize('-')[0]}", "${profile.baseName.tokenize('-')[1]}")
run("$profile", "$lineage", "${profile.baseName.tokenize('-')[0]}", "${profile.baseName.tokenize('-')[1]}")

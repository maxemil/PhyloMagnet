#! ${params.python3}
import pandas as pd

def run(profile, lineage, samplename, genename):
    df = pd.read_csv(profile, sep='\\t')
    df['taxopath'] = df['taxopath'].astype(str)
    lineage_results = df[df['taxopath'].apply(lambda x: lineage in x and not '{};'.format(lineage) in x)]
    with open('decision.txt', 'w') as out:
        try:
            if lineage_results['aLWR'].item() >= 1:
                print("{} present in sample {} (tree {})".format(lineage, samplename, genename))
                print("{}\\t{}\\t{}\\tTrue".format(samplename, genename, lineage), file=out)
            else:
                print("{} NOT present in sample {} (tree {})".format(lineage, samplename, genename))
                print("{}\\t{}\\t{}\\tFalse".format(samplename, genename, lineage), file=out)
        except ValueError:
            print("{} NOT present in sample {} (tree {})".format(lineage, samplename, genename))
            print("{}\\t{}\\t{}\\tFalse".format(samplename, genename, lineage), file=out)

run("$profile", "$lineage", "${profile.baseName.tokenize('-')[0]}", "${profile.baseName.tokenize('-')[1]}")

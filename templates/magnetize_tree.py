#! ${params.python3}
import pandas as pd

def run(profile, lineage, basename):
    df = pd.read_csv(profile, sep='\\t')
    df['taxopath'] = df['taxopath'].astype(str)
    lineage_results = df[df['taxopath'].apply(lambda x: lineage in x and not '{};'.format(lineage) in x)]
    with open('decision.txt', 'w') as out:
        try:
            if lineage_results['aLWR'].item() >= 1:
                print("{} present in sample {}".format(lineage, basename))
                print("{}\\t{}\\tTrue".format(basename, lineage), file=out)
            else:
                print("{} NOT present in sample {}".format(lineage, basename))
                print("{}\\t{}\\tFalse".format(basename, lineage), file=out)
        except ValueError:
            print("{} NOT present in sample {}".format(lineage, basename))
            print("{}\\t{}\\tFalse".format(basename, lineage), file=out)

run("$profile", "$lineage", "${profile.baseName}")

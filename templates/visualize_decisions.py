#! ${params.python3}
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
from itertools import product
import matplotlib.pyplot as plt


def get_counts_sample_taxon(df):
    taxa = set(df["Taxon"])
    samples = set(df["Sample"])
    counts = pd.DataFrame(columns=samples, index=taxa)
    counts = counts.fillna(0)

    for t,s in product(taxa, samples):
        counts[s].loc[t] = df[(df["Sample"] == s) & (df["Taxon"] == t)].shape[0]
    return counts

def get_taxa_to_remove(counts, trees):
    low = counts.apply(lambda x: all(x < threshold * len(trees)), axis=1)
    low  = low[low == True].index
    return low

def make_barplot(df):
    ax = sns.factorplot(hue="Taxon", x="Sample", data=df, kind="count", legend_out=True, size=8, aspect=1.1)
    plt.savefig('decision_barplot.pdf')
    plt.clf()
    plt.cla()

def make_heatmap(counts):
    ax = sns.heatmap(counts, annot=True, cmap='Reds', xticklabels=True)
    ax.figure.tight_layout()
    plt.savefig('decision_heatmap.pdf')

def main(threshold, infile, filter_taxa=True):
    df = pd.read_csv(infile, sep='\\t', header=None, names=['Sample', 'Tree', 'Taxon', 'value'], dtype=str)
    df = df[df['value'] == 'True']

    counts = get_counts_sample_taxon(df)
    low = get_taxa_to_remove(counts, set(df['Tree']))

    df = df[df.apply(lambda row: not any([row['Taxon'] == l for l in low]), axis=1)]
    make_barplot(df)

    if filter_taxa:
        counts = counts.loc[[i not in low for i in counts.index]]

    make_heatmap(counts)

if __name__ == '__main__':
    threshold = 0.25
    infile = "tree_decisions.txt"
    main(threshold, infile, False)

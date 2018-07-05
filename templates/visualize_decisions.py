#! ${params.python3}
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
from itertools import product
import matplotlib.pyplot as plt

threshold = 0.25
infile = 'tree_decisions.txt'

df = pd.read_csv(infile, sep='\\t', header=None, names=['Sample', 'Tree', 'Taxon', 'value'])
df = df[df['value'] == True]

taxa = set(df["Taxon"])
samples = set(df["Sample"])
trees = set(df['Tree'])

df2 = pd.DataFrame(columns=samples, index=taxa)
df2 = df2.fillna(0)

for t,s in product(taxa, samples):
    df2[s].loc[t] = df[(df["Sample"] == s) & (df["Taxon"] == t)].shape[0]

low = df2.apply(lambda x: all(x < threshold * len(trees)), axis=1)
low  = low[low == True].index

df = df[df.apply(lambda row: not any([row['Taxon'] == l for l in low]), axis=1)]

ax = sns.factorplot(hue="Taxon", x="Sample", data=df, kind="count", legend_out=True, size=len(samples), aspect=1.1)
plt.savefig('decision_barplot.pdf')
plt.clf()
plt.cla()

df2 = df2.loc[[i not in low for i in df2.index]]
ax = sns.heatmap(df2, annot=True, cmap='Reds', xticklabels=True, center=threshold * len(trees))
ax.figure.tight_layout()
plt.savefig('decision_heatmap.pdf')
